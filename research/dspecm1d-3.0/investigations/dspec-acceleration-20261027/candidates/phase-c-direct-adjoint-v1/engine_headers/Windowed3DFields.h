#ifndef DSPECM1D_WINDOWED_3D_FIELDS_H
#define DSPECM1D_WINDOWED_3D_FIELDS_H

#include "../profile.hpp"
#include "SensitivityKernels.h"
#include "../fast_engine.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <complex>
#include <cstdint>
#include <functional>
#include <exception>
#include <omp.h>
#include <optional>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace SPARSESPEC {

struct Field3DPoint {
  double radiusMetres = 0.0;
  double thetaRadians = 0.0;
  double phiRadians = 0.0;
  // Select a one-sided element at a shared node/source trace; -1 selects the
  // outer element by default. No quadrature weight is implied by a point.
  int radialElement = -1;
};

struct SpinStrainField {
  std::complex<double> rr{}, horizontalTrace{}, radialPlus{}, radialMinus{},
      horizontalPlus{}, horizontalMinus{};
};

struct Field3DSample {
  int radialElement = -1;
  bool fluid = false;
  std::array<std::complex<double>, 3> forwardRawDisplacement{};
  std::array<std::complex<double>, 3>
      observableForwardDisplacementSpectrumMetreSeconds{};
  std::array<std::complex<double>, 3> adjointDualDisplacement{};
  SpinStrainField forwardRawStrain, adjointDualStrain;
  SpinStrainField rScaledForwardStrainMetres, rScaledAdjointDualStrainMetres;
};

/** Compact kernel sample; receiver displacement is present only when requested. */
struct KernelField3DSample {
  int radialElement=-1;
  bool fluid=false;
  std::array<std::complex<double>,3> forwardRawDisplacement{};
  SpinStrainField forwardRawStrain,adjointDualStrain;
};

/** Bounded per-frequency reuse payload for tensor-volume objective batches. */
struct CachedForwardTensorFields3D {
  std::vector<SpinStrainField> volumeStrain;
  std::vector<std::array<std::complex<double>,3>> receiverDisplacement;
  std::size_t payloadBytes() const noexcept {
    return volumeStrain.size()*sizeof(SpinStrainField)+
        receiverDisplacement.size()*sizeof(std::array<std::complex<double>,3>);
  }
};

/** Synchronous view of one solved angular family/degree. References remain
 * valid only during the visitor call; this avoids retaining modal histories. */
struct ModalFrequencyBlock3D {
  SensitivitySolveFamily family;
  int degree = 0;
  int sourceElement = -1;
  const Eigen::MatrixXcd &forward;
  const Eigen::MatrixXcd &adjoint;
  const Eigen::MatrixXcd &sourceLift;
  const Eigen::MatrixXcd &sourceCoefficients;
  const Eigen::MatrixXcd &modeReceiver;
};

struct OwnedModalFrequencyBlock3D {
  SensitivitySolveFamily family;
  int degree = 0;
  int sourceElement = -1;
  Eigen::MatrixXcd forward, adjoint, sourceLift, sourceCoefficients, modeReceiver;
};

namespace windowed3d_detail {

/** Private-workspace modal basis used to synthesize a pole-centred field. */
struct PoleModalFrequencyBlock3D {
  SensitivitySolveFamily family;
  int degree = 0;
  int sourceElement = -1;
  Eigen::MatrixXcd forward, adjointUnit, sourceLift, sourceCoefficients;
  Eigen::MatrixXcd adjointBasisTransform;
  // Direct measurement-direction solutions, one dof-by-m block per receiver.
  // Populated only by the Phase C route; the original basis route is retained
  // as the comparison oracle.
  std::vector<Eigen::MatrixXcd> adjointDirect;
};

/** Geocentric orthonormal frame columns (theta, east, radial). */
struct PoleFrame3D {
  std::array<std::array<double,3>,3> axes{};
};

struct PoleSynthesisBlock3D {
  SensitivitySolveFamily family;
  int degree = 0;
  int sourceElement = -1;
  bool adjoint = false;
  Eigen::MatrixXcd fields, lift;
};

} // namespace windowed3d_detail

struct PoleAngularCoefficientRecord3D {
  windowed3d_detail::PoleFrame3D frame;
  std::uint64_t synthesisSet=0;
  bool adjoint=false;
  int element=-1;
  double radiusModel=0.0;
  std::vector<std::pair<int,int>> labels; // (m,l), sorted by m then l
  std::vector<std::array<std::complex<double>,6>> rawStrain;
};

class Windowed3DFrequencyWorkspace;

/**
 * Stream one frequency of full angular forward/adjoint fields at requested
 * points. `adjointWeights` is one frequency column after
 * WindowedDisplacementMap::adjoint. Forward modal displacement and strain
 * remain raw for later weak-form contractions. The separate observable is
 * the unwindowed damped-Laplace displacement spectrum at physical frequency
 * z/timeNorm, in metre-seconds; it applies L0*T0*(-i/z), and is not the
 * P-windowed output. The adjoint field is a dual weighted solution, not a
 * physical SI displacement. No point or quadrature weights are implied.
 * r-scaled strains are in metres.
 */
class Windowed3DFieldEngine {
  friend class Windowed3DFrequencyWorkspace;
public:
  static constexpr std::size_t poleCoefficientCacheByteBudget=512u*1024u*1024u;
  using Complex = std::complex<double>;
  using Point = Field3DPoint;
  using Sample = Field3DSample;
  using ModalVisitor = std::function<void(const ModalFrequencyBlock3D &)>;
  using CachedAngularWigner=GSHTrans::Wigner<double,GSHTrans::All,GSHTrans::All,
                                               GSHTrans::Multiple>;
  struct AngularCoefficientTable {
    int maxDegree=0,nMax=0;
    std::size_t thetaCount=0;
    std::vector<double> values;
    const double *row(std::size_t theta,int n,int m) const {
      const std::size_t nSlot=static_cast<std::size_t>(n+nMax);
      const std::size_t mSlot=static_cast<std::size_t>(m+2);
      const std::size_t stride=static_cast<std::size_t>(maxDegree+1);
      return values.data()+(((theta*static_cast<std::size_t>(2*nMax+1)+nSlot)*5+mSlot)*stride);
    }
  };

  /** Exact element/radius interpolation weights reusable across frequencies.
   * The explicit byte limit covers this cache only, not field or angular data. */
  struct SamplingPlanCache {
    struct PointAngle {
      double theta=0.0,phi=0.0,a=1.0,b=0.0;
      std::array<Complex,5> orderPhases{};
    };
    struct RadialWeights {
      int element=-1;
      std::vector<double> shape,derivative;
    };
    std::map<std::pair<int,double>,RadialWeights> radialWeights;
    using FrameKey=std::array<std::uint64_t,9>;
    using PointKey=std::pair<std::uint64_t,std::uint64_t>;
    using TensorTraceKey=std::pair<int,std::uint64_t>;
    struct TensorGridLayout {
      std::vector<TensorTraceKey> traces;
      std::vector<PointKey> angleKeys;
      std::vector<std::map<PointKey,std::size_t>> pointForTrace;
      std::size_t payloadBytes=0;
    };
    struct TensorAngularPlan {
      std::vector<double> localTheta,localPhi,rotA,rotB;
      std::shared_ptr<CachedAngularWigner> wigner;
      std::map<std::vector<std::uint64_t>,std::shared_ptr<Eigen::MatrixXcd>> complexBases;
      std::map<std::vector<std::uint64_t>,std::shared_ptr<Eigen::MatrixXd>> realBases;
      std::size_t payloadBytes=0;
    };
    std::map<FrameKey,std::map<PointKey,PointAngle>> pointAngles;
    std::map<std::vector<std::uint64_t>,std::shared_ptr<CachedAngularWigner>> angularTables;
    std::map<std::vector<std::uint64_t>,std::shared_ptr<AngularCoefficientTable>> angularCoefficients;
    // Packed synthesis matrices depend only on the semantic angle grid, spin,
    // and ordered (m,l) columns, so reuse them across radial traces/frequencies.
    std::map<std::vector<std::uint64_t>,std::shared_ptr<Eigen::MatrixXcd>> packedBases;
    // Tensor grid address maps and transformed angular plans are shared across
    // objectives/frequencies. Keys contain exact coordinate/frame bit patterns.
    std::map<std::vector<std::uint64_t>,std::shared_ptr<TensorGridLayout>> tensorGridLayouts;
    std::map<std::vector<std::uint64_t>,std::shared_ptr<TensorAngularPlan>> tensorAngularPlans;
    std::size_t tensorGridHits=0,tensorGridMisses=0;
    std::size_t tensorPlanHits=0,tensorPlanMisses=0;
    std::size_t tensorBasisHits=0,tensorBasisMisses=0;
    std::size_t tensorGridPayloadBytes=0,tensorAngularPayloadBytes=0;
    std::size_t bytes=0;
    static constexpr std::size_t byteBudget=32u*1024u*1024u;
  };

  /** Algebraic r*strain reconstruction, safe at r=0 when modal data are finite. */
  static SpinStrainField rScaledStrainFromModal(
      int degree, Complex U, Complex Up, Complex V, Complex Vp,
      Complex W, Complex Wp, double rModel, Complex x0, Complex xp,
      Complex xm, Complex xpp, Complex xmm) {
    const double lambda = static_cast<double>(degree) * (degree + 1.0);
    const double zeta = std::sqrt(lambda);
    const double eta = degree >= 2 ? std::sqrt((degree - 1.0) * (degree + 2.0)) : 0.0;
    return {rModel*Up*x0, (2.0*U-lambda*V)*x0,
      -zeta/2.0*((rModel*Vp-V+U)+Complex(0,1)*(rModel*Wp-W))*xp,
      zeta/2.0*((rModel*Vp-V+U)-Complex(0,1)*(rModel*Wp-W))*xm,
      zeta*eta*(V+Complex(0,1)*W)*xpp,
      zeta*eta*(V-Complex(0,1)*W)*xmm};
  }

  static std::vector<Sample> evaluateFrequency(
      Full1D::SEM &sem, SourceInfo::EarthquakeCMT &source,
      InputParameters &parameters, Complex z,
      const Eigen::VectorXcd &adjointWeights,
      const std::vector<Point> &points);

  /** Build an owned one-frequency solve workspace for repeated point chunks. */
  static Windowed3DFrequencyWorkspace prepareFrequencyWorkspace(
      Full1D::SEM &sem,SourceInfo::EarthquakeCMT &source,
      InputParameters &parameters,Complex z,const Eigen::VectorXcd &adjointWeights,
      dssolver::BatchResponseCache *batchCache=nullptr,bool poleOnly=false,
      bool forwardOnly=false,double *basisReconstructionRelative=nullptr,int solveThreads=1,std::optional<Complex> attenuation=std::nullopt);

  static std::vector<Sample> sampleModalBlocks(
      Full1D::SEM &sem,InputParameters &parameters,Complex z,
      const std::vector<OwnedModalFrequencyBlock3D> &blocks,
      const std::vector<Point> &points,bool allowCenter=false,
      bool kernelRoute=false,std::size_t displacementBegin=0) {
    if(points.empty())throw std::invalid_argument("3D field chunk must be nonempty.");
    const double scale=sem.lengthScaleMetres();
    if(!(std::isfinite(scale)&&scale>0.0&&std::isfinite(sem.timeScaleSeconds())&&
         sem.timeScaleSeconds()>0.0))throw std::invalid_argument("3D SI scales must be positive and finite.");
    std::vector<Sample> result(points.size());
    std::vector<double> radii(points.size());std::vector<int> elements(points.size());
    std::vector<std::vector<std::size_t>> radialGroups;
    std::map<std::pair<int,double>,std::size_t> groupIndex;
    const auto tol=[](double r,double b){return 8.0*std::numeric_limits<double>::epsilon()*
        std::max({1.0,std::abs(r),std::abs(b)});};
    for(std::size_t p=0;p<points.size();++p){
      const auto &point=points[p];
      if(!(std::isfinite(point.radiusMetres)&&(allowCenter?point.radiusMetres>=0.0:point.radiusMetres>0.0))||
         !std::isfinite(point.thetaRadians)||point.thetaRadians<0.0||point.thetaRadians>EIGEN_PI||
         !std::isfinite(point.phiRadians))throw std::invalid_argument("Invalid 3D field chunk point.");
      double r=point.radiusMetres/scale;int e=point.radialElement;
      if(e< -1||e>=static_cast<int>(sem.mesh().NE()))throw std::invalid_argument("Invalid one-sided field element.");
      if(e<0)for(int candidate=0;candidate<sem.mesh().NE();++candidate)
        if(r>=sem.mesh().ELR(candidate)-tol(r,sem.mesh().ELR(candidate))&&
           r<=sem.mesh().EUR(candidate)+tol(r,sem.mesh().EUR(candidate)))e=candidate;
      if(e<0)throw std::invalid_argument("3D field chunk point lies outside the SEM mesh.");
      const double bounds[2]={sem.mesh().ELR(e),sem.mesh().EUR(e)};
      for(double b:bounds)if(b>0.0&&std::abs(r-b)<=tol(r,b))r=b;
      if(r<sem.mesh().ELR(e)||r>sem.mesh().EUR(e)||
         r>sem.mesh().EUR(sem.mesh().NE()-1))throw std::invalid_argument("One-sided field element does not contain the point.");
      radii[p]=r;elements[p]=e;result[p].radialElement=e;result[p].fluid=sem.elementIsFluid(e);
      const auto [group,added]=groupIndex.try_emplace({e,r},radialGroups.size());
      if(added)radialGroups.emplace_back();
      radialGroups[group->second].push_back(p);
    }
    int maxAngularDegree=-1;
    for(const auto &block:blocks)if(block.family!=SensitivitySolveFamily::Radial)
      maxAngularDegree=std::max(maxAngularDegree,block.degree);
    if(maxAngularDegree<0){
      for(const auto &owned:blocks){
        const ModalFrequencyBlock3D block{owned.family,owned.degree,owned.sourceElement,
          owned.forward,owned.adjoint,owned.sourceLift,owned.sourceCoefficients,owned.modeReceiver};
        accumulateModalFields(sem,block,points,radii,elements,radialGroups,scale,result,nullptr,0,{},kernelRoute,displacementBegin);
      }
    }else{
      std::map<std::uint64_t,std::size_t> angleMap;
      std::vector<double> uniqueTheta;
      std::vector<std::size_t> angleForPoint(points.size());
      for(std::size_t p=0;p<points.size();++p){
        const auto [it,added]=angleMap.try_emplace(
            std::bit_cast<std::uint64_t>(points[p].thetaRadians),uniqueTheta.size());
        if(added)uniqueTheta.push_back(points[p].thetaRadians);
        angleForPoint[p]=it->second;
      }
      constexpr std::size_t tableBudget=32u*1024u*1024u;
      const std::size_t perAngle=angularWignerBytesPerAngle(maxAngularDegree,maxAngularDegree);
      const std::size_t angleChunkSize=std::max<std::size_t>(1,
          std::min<std::size_t>(64,tableBudget/std::max<std::size_t>(1,perAngle)));
      const std::size_t chunkCount=(uniqueTheta.size()+angleChunkSize-1)/angleChunkSize;
      std::vector<AngularPointChunk> chunks(chunkCount);
      for(std::size_t c=0;c<chunkCount;++c){
        const std::size_t begin=c*angleChunkSize;
        const std::size_t end=std::min(uniqueTheta.size(),begin+angleChunkSize);
        chunks[c].theta.assign(uniqueTheta.begin()+begin,uniqueTheta.begin()+end);
      }
      for(std::size_t g=0;g<radialGroups.size();++g)for(const auto p:radialGroups[g]){
        const std::size_t c=angleForPoint[p]/angleChunkSize;
        auto &chunk=chunks[c];
        if(chunk.lastRadialGroup!=g){chunk.radialGroups.emplace_back();chunk.lastRadialGroup=g;}
        chunk.radialGroups.back().push_back(p);
      }
      for(std::size_t c=0;c<chunks.size();++c){
        const std::size_t angleBegin=c*angleChunkSize;
        AngularWigner w(maxAngularDegree,maxAngularDegree,std::min(2,maxAngularDegree),chunks[c].theta);
        for(const auto &owned:blocks){
          const ModalFrequencyBlock3D block{owned.family,owned.degree,owned.sourceElement,
            owned.forward,owned.adjoint,owned.sourceLift,owned.sourceCoefficients,owned.modeReceiver};
          accumulateModalFields(sem,block,points,radii,elements,chunks[c].radialGroups,
              scale,result,&w,angleBegin,angleForPoint,kernelRoute,displacementBegin);
        }
      }
    }
    const Complex fz=-Complex(0.0,1.0)/z;
    for(auto &sample:result)for(int c=0;c<3;++c)
      sample.observableForwardDisplacementSpectrumMetreSeconds[c]=
          sem.timeScaleSeconds()*fz*scale*sample.forwardRawDisplacement[c];
    return result;
  }

  /** Center-safe reconstruction for an explicitly weighted volume rule.
   * Raw point strains are intentionally unavailable at r=0; algebraic
   * r-scaled strains remain defined for the caller's r^2-weighted integral.
   */
  static std::vector<Sample> evaluateFrequencyForWeightedVolume(
      Full1D::SEM &sem, SourceInfo::EarthquakeCMT &source,
      InputParameters &parameters, Complex z,
      const Eigen::VectorXcd &adjointWeights,
      const std::vector<Point> &points);

private:
  using PoleModalBlock=windowed3d_detail::PoleModalFrequencyBlock3D;
  using PoleFrame=windowed3d_detail::PoleFrame3D;
  using PoleSynthesisBlock=windowed3d_detail::PoleSynthesisBlock3D;
  struct RadialMode { Complex U{}, Up{}, V{}, Vp{}, W{}, Wp{}; };
  struct PoleStrainCoefficients {
    std::array<Complex,6> raw{},rScaled{};
  };
  struct PoleDegreeRow {
    std::vector<RadialMode> values;
    std::vector<PoleStrainCoefficients> strain;
    std::vector<unsigned char> present;
  };
  using PoleCoefficientRows=std::map<int,PoleDegreeRow>;
  struct PoleCoefficientCache {
    using Key=std::tuple<SamplingPlanCache::FrameKey,std::uint64_t,int,std::uint64_t>;
    std::map<Key,PoleCoefficientRows> rows;
    std::size_t bytes=0;
    static constexpr std::size_t byteBudget=Windowed3DFieldEngine::poleCoefficientCacheByteBudget;
  };
  using AngularWigner = GSHTrans::Wigner<double,GSHTrans::All,GSHTrans::All,
                                           GSHTrans::Multiple>;

  struct AngularPointChunk {
    std::vector<double> theta;
    std::vector<std::vector<std::size_t>> radialGroups;
    std::size_t lastRadialGroup=std::numeric_limits<std::size_t>::max();
  };

  using RadialSamplingPlan=SamplingPlanCache::RadialWeights;

  static RadialSamplingPlan makeRadialSamplingPlan(
      const Full1D::SEM &sem,double radius,int element) {
    const auto &mesh=sem.mesh();
    std::vector<double> nodes(static_cast<std::size_t>(mesh.NN()));
    for(int q=0;q<mesh.NN();++q)nodes[static_cast<std::size_t>(q)]=mesh.NodeRadius(element,q);
    const auto basis=Interpolation::LagrangeBasis(nodes);
    RadialSamplingPlan plan;plan.element=element;
    plan.shape.resize(nodes.size());plan.derivative.resize(nodes.size());
    const char *nativeDerivativeEnv=std::getenv("DSFAST_NATIVE_DERIVATIVE");
    const bool physicalNodeDerivative=nativeDerivativeEnv&&
        std::string(nativeDerivativeEnv)=="physical";
    for(int i=0;i<mesh.NN();++i){
      const double tol=16.0*std::numeric_limits<double>::epsilon()*
          std::max({1.0,std::abs(radius),std::abs(nodes[static_cast<std::size_t>(i)])});
      if(std::abs(radius-nodes[static_cast<std::size_t>(i)])<=tol){
        if(physicalNodeDerivative){
          // Compatibility alternative: preserve the physical-coordinate
          // Lagrange derivative evaluated at the nodal radius. These weights
          // are cached by (element,radius) just like the Dref/J branch below.
          std::fill(plan.shape.begin(),plan.shape.end(),0.0);
          plan.shape[static_cast<std::size_t>(i)]=1.0;
          for(int j=0;j<mesh.NN();++j)
            plan.derivative[static_cast<std::size_t>(j)]=basis.template Evaluate<1>(j,radius);
          return plan;
        }
        const auto dref=sensitivity_detail::derivativeMatrix(sem);
        const double jacobian=0.5*(mesh.EUR(element)-mesh.ELR(element));
        if(!(std::isfinite(jacobian)&&jacobian>0.0))
          throw std::runtime_error("Native SEM element has an invalid radial Jacobian.");
        std::fill(plan.shape.begin(),plan.shape.end(),0.0);
        std::fill(plan.derivative.begin(),plan.derivative.end(),0.0);
        plan.shape[static_cast<std::size_t>(i)]=1.0;
        for(int j=0;j<mesh.NN();++j)
          plan.derivative[static_cast<std::size_t>(j)]=dref(i,j)/jacobian;
        return plan;
      }
    }
    for(int q=0;q<mesh.NN();++q){
      plan.shape[static_cast<std::size_t>(q)]=basis.Evaluate(q,radius);
      plan.derivative[static_cast<std::size_t>(q)]=basis.template Evaluate<1>(q,radius);
    }
    return plan;
  }

  static Eigen::MatrixXcd interpolatePlanned(const Full1D::SEM &sem,
      SensitivitySolveFamily family,const Eigen::MatrixXcd &fields,
      const Eigen::MatrixXcd &lift,const RadialSamplingPlan &plan,
      int column,int liftElement) {
    const int element=plan.element;
    const bool useLift=element==liftElement&&lift.rows()==fields.rows()&&
                       lift.cols()==fields.cols();
    const auto value=[&](int dof){
      return fields(dof,column)+(useLift?lift(dof,column):Complex{});
    };
    Eigen::MatrixXcd out=Eigen::MatrixXcd::Zero(6,1);
    for(int q=0;q<sem.mesh().NN();++q){
      const auto index=static_cast<std::size_t>(q);
      const double shape=plan.shape[index],derivative=plan.derivative[index];
      if(family==SensitivitySolveFamily::Spheroidal){
        const auto iu=sem.ltgS(0,element,q),iv=sem.ltgS(1,element,q);
        out(0,0)+=shape*value(iu);out(1,0)+=derivative*value(iu);
        out(2,0)+=shape*value(iv);out(3,0)+=derivative*value(iv);
      }else if(family==SensitivitySolveFamily::Toroidal){
        const auto iw=sem.ltgT(element,q);
        out(4,0)+=shape*value(iw);out(5,0)+=derivative*value(iw);
      }else{
        const auto ir=sem.ltgR(0,element,q);
        out(0,0)+=shape*value(ir);out(1,0)+=derivative*value(ir);
      }
    }
    return out;
  }

  static std::size_t angularWignerBytesPerAngle(int maxDegree,int mMax) {
    const int nMax=std::min(2,maxDegree);
    mMax=std::min(mMax,maxDegree);
    std::size_t values=0;
    for(int n=-nMax;n<=nMax;++n)
      values+=GSHTrans::GSHIndices<GSHTrans::All>(maxDegree,mMax,n).Size();
    return values*sizeof(double)+static_cast<std::size_t>(2*nMax+1)*sizeof(std::size_t);
  }

  static std::size_t angularCoefficientBytesPerAngle(int maxDegree) {
    const int nMax=std::min(2,maxDegree);
    return static_cast<std::size_t>(2*nMax+1)*5u*
        static_cast<std::size_t>(maxDegree+1)*sizeof(double);
  }

  static std::size_t angularTableBytesPerAngle(int maxDegree,int mMax) {
    return angularWignerBytesPerAngle(maxDegree,mMax)+angularCoefficientBytesPerAngle(maxDegree);
  }

  static Complex harmonic(int degree, int m, int n, double theta, double phi) {
    GSHTrans::Wigner<double, GSHTrans::All, GSHTrans::All,
                     GSHTrans::Single> w(degree, degree, std::min(2, degree), theta);
    return w[n][degree, m] * std::exp(Complex(0.0, m * phi));
  }

  template <class Wigner>
  static Complex harmonic(const Wigner &w, int degree, int m, int n,
                          std::size_t thetaIndex, Complex phase) {
    return w[n, static_cast<int>(thetaIndex)][degree, m] * phase;
  }

  static Eigen::MatrixXcd embedReceiver(const Full1D::SEM &sem,
      const Eigen::MatrixXcd &local, std::size_t low, int dofs) {
    if (low + local.cols() > dofs)
      throw std::invalid_argument("3D field receiver slice exceeds the SEM operator.");
    Eigen::MatrixXcd full = Eigen::MatrixXcd::Zero(local.rows(), dofs);
    full.block(0, low, local.rows(), local.cols()) = local;
    return full;
  }

  static PoleFrame poleFrame(double latitudeDegrees,double longitudeDegrees) {
    const double theta=(90.0-latitudeDegrees)*EIGEN_PI/180.0;
    const double phi=longitudeDegrees*EIGEN_PI/180.0;
    const double st=std::sin(theta),ct=std::cos(theta),sp=std::sin(phi),cp=std::cos(phi);
    return {{{{ct*cp,ct*sp,-st},{-sp,cp,0.0},{st*cp,st*sp,ct}}}};
  }

  static std::array<double,3> frameCoordinates(const PoleFrame &frame,
      const std::array<double,3> &global) {
    std::array<double,3> local{};
    for(int j=0;j<3;++j)for(int i=0;i<3;++i)local[j]+=frame.axes[j][i]*global[i];
    return local;
  }

  static std::array<double,3> frameVector(const PoleFrame &frame,
      const std::array<double,3> &local) {
    std::array<double,3> global{};
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)global[i]+=frame.axes[j][i]*local[j];
    return global;
  }

  static std::array<double,3> tangentRotation(const PoleFrame &frame,
      double theta,double phi,double &a,double &b) {
    const auto x=frameCoordinates(frame,unitPoint(theta,phi));
    const double horizontal=std::hypot(x[0],x[1]);
    const double localTheta=std::atan2(horizontal,x[2]);
    const double localPhi=horizontal==0.0?0.0:std::atan2(x[1],x[0]);
    const auto localEt=frameVector(frame,{std::cos(localTheta)*std::cos(localPhi),
        std::cos(localTheta)*std::sin(localPhi),-std::sin(localTheta)});
    const std::array<double,3> globalEt={std::cos(theta)*std::cos(phi),
        std::cos(theta)*std::sin(phi),-std::sin(theta)};
    const std::array<double,3> globalEp={-std::sin(phi),std::cos(phi),0.0};
    a=0.0;b=0.0;
    for(int i=0;i<3;++i){a+=globalEt[i]*localEt[i];b+=globalEp[i]*localEt[i];}
    return {localTheta,localPhi,0.0};
  }

  // Match the validated sampler's boundary snap before forming bitwise keys.
  // Keep the element in the key so opposite sides of an interface stay distinct.
  static double canonicalSamplingRadius(const Full1D::SEM &sem,int element,double r) {
    if(element<0||element>=sem.mesh().NE())
      throw std::invalid_argument("Tensor sampling requires an explicit valid radial element.");
    for(double b:{sem.mesh().ELR(element),sem.mesh().EUR(element)}) {
      const double tol=8.0*std::numeric_limits<double>::epsilon()*
          std::max({1.0,std::abs(r),std::abs(b)});
      if(b>0.0&&std::abs(r-b)<=tol)r=b;
    }
    return r;
  }

  static std::array<Complex,3> rotateDisplacement(const std::array<Complex,3> &v,
      double a,double b) {
    return {v[0],a*v[1]+b*v[2],-b*v[1]+a*v[2]};
  }

  static SpinStrainField rotateStrain(const SpinStrainField &v,double a,double b) {
    const Complex c(a,-b),c2=c*c;
    return {v.rr,v.horizontalTrace,c*v.radialPlus,std::conj(c)*v.radialMinus,
        c2*v.horizontalPlus,std::conj(c2)*v.horizontalMinus};
  }

  static void rotateContribution(Sample &sample,double a,double b,bool adjoint) {
    if(adjoint){
      sample.adjointDualDisplacement=rotateDisplacement(sample.adjointDualDisplacement,a,b);
      sample.adjointDualStrain=rotateStrain(sample.adjointDualStrain,a,b);
      sample.rScaledAdjointDualStrainMetres=rotateStrain(sample.rScaledAdjointDualStrainMetres,a,b);
    }else{
      sample.forwardRawDisplacement=rotateDisplacement(sample.forwardRawDisplacement,a,b);
      sample.forwardRawStrain=rotateStrain(sample.forwardRawStrain,a,b);
      sample.rScaledForwardStrainMetres=rotateStrain(sample.rScaledForwardStrainMetres,a,b);
    }
  }

  static std::array<double,3> unitPoint(double theta,double phi) {
    return {std::sin(theta)*std::cos(phi),std::sin(theta)*std::sin(phi),std::cos(theta)};
  }

  static Eigen::MatrixXcd receiverAtPole(SensitivitySolveFamily family,int degree) {
    const int mmax=std::min(1,degree),orders=2*mmax+1;
    Eigen::MatrixXcd result=Eigen::MatrixXcd::Zero(3,orders);
    const double norm=std::sqrt((2.0*degree+1.0)/(4.0*EIGEN_PI));
    for(int m=-mmax;m<=mmax;++m){
      const int col=m+mmax;
      const Complex y0=m==0?Complex(norm,0.0):Complex{};
      const Complex ym=m==-1?Complex(norm,0.0):Complex{};
      const Complex yp=m==1?Complex(norm,0.0):Complex{};
      if(family==SensitivitySolveFamily::Spheroidal)
        result(0,col)=y0;
      if(family==SensitivitySolveFamily::Spheroidal){
        result(1,col)=yp-ym;result(2,col)=-Complex(0,1)*(ym+yp);
      }else{
        result(1,col)=Complex(0,1)*(ym+yp);result(2,col)=yp-ym;
      }
    }
    return result;
  }

  static Eigen::MatrixXcd interpolate(const Full1D::SEM &sem,
      SensitivitySolveFamily family, const Eigen::MatrixXcd &fields,
      const Eigen::MatrixXcd &lift, double radius, int element, int column = 0,
      int liftElement = -1) {
    const bool useLift = element == liftElement && lift.rows() == fields.rows() &&
                         lift.cols() == fields.cols();
    const auto value = [&](int dof) {
      return fields(dof, column) + (useLift ? lift(dof, column) : Complex{});
    };
    return interpolateOne(sem, family, radius, element, value);
  }

  static Eigen::MatrixXcd interpolateCombination(const Full1D::SEM &sem,
      SensitivitySolveFamily family, const Eigen::MatrixXcd &fields,
      const Eigen::MatrixXcd &lift, const Eigen::VectorXcd &coefficients,
      double radius, int element, int liftElement) {
    if (coefficients.size() != fields.cols())
      throw std::invalid_argument("Modal source coefficients do not match SEM columns.");
    const bool useLift = element == liftElement && lift.rows() == fields.rows() &&
                         lift.cols() == fields.cols();
    const auto value = [&](int dof) {
      Complex result = fields.row(dof).cwiseProduct(coefficients.transpose()).sum();
      if (useLift)
        result += lift.row(dof).cwiseProduct(coefficients.transpose()).sum();
      return result;
    };
    return interpolateOne(sem, family, radius, element, value);
  }

  template <class Value>
  static Eigen::MatrixXcd interpolateOne(const Full1D::SEM &sem,
      SensitivitySolveFamily family, double radius, int element, Value value) {
    const auto &mesh = sem.mesh();
    const int nq = mesh.NN();
    std::vector<double> nodes(nq);
    for (int q = 0; q < nq; ++q) nodes[q] = mesh.NodeRadius(element, q);
    const auto basis = Interpolation::LagrangeBasis(nodes);
    Eigen::MatrixXcd out = Eigen::MatrixXcd::Zero(6, 1);
    for (int q = 0; q < nq; ++q) {
      const double shape = basis.Evaluate(q, radius);
      const double derivative = basis.template Evaluate<1>(q, radius);
      if (family == SensitivitySolveFamily::Radial) {
        const auto i = sem.ltgR(0, element, q);
        out(0, 0) += shape * value(i);
        out(1, 0) += derivative * value(i);
      } else if (family == SensitivitySolveFamily::Spheroidal) {
        const auto iu = sem.ltgS(0, element, q), iv = sem.ltgS(1, element, q);
        out(0, 0) += shape * value(iu); out(1, 0) += derivative * value(iu);
        out(2, 0) += shape * value(iv); out(3, 0) += derivative * value(iv);
      } else {
        const auto iw = sem.ltgT(element, q);
        out(4, 0) += shape * value(iw); out(5, 0) += derivative * value(iw);
      }
    }
    return out;
  }

  static RadialMode mode(const Eigen::MatrixXcd &v) {
    return {v(0,0), v(1,0), v(2,0), v(3,0), v(4,0), v(5,0)};
  }

  static std::array<Complex,3> radialToDisplacement(int l,
      const Eigen::MatrixXcd &values, double r, Complex x0, Complex xp,
      Complex xm, Complex, Complex) {
    return radialToDisplacement(l,mode(values),x0,xp,xm);
  }

  static std::array<Complex,3> radialToDisplacement(int l,const RadialMode &a,
      Complex x0,Complex xp,Complex xm) {
    const double zeta = std::sqrt(l * (l + 1.0));
    return {a.U*x0, zeta/2.0*a.V*(xp-xm), -Complex(0,1)*zeta/2.0*a.V*(xp+xm)};
  }

  static std::array<Complex,3> toroidalDisplacement(int l,
      const Eigen::MatrixXcd &values, Complex xp, Complex xm) {
    const auto a = mode(values); const double zeta = std::sqrt(l * (l + 1.0));
    return {Complex{}, Complex(0,1)*zeta/2.0*a.W*(xp+xm), zeta/2.0*a.W*(xp-xm)};
  }

  static SpinStrainField radialToStrain(int l, const Eigen::MatrixXcd &values,
      double r, Complex x0, Complex xp, Complex xm, Complex xpp, Complex xmm) {
    return radialToStrain(l,mode(values),r,x0,xp,xm,xpp,xmm);
  }

  static SpinStrainField radialToStrain(int l,const RadialMode &a,double r,
      Complex x0,Complex xp,Complex xm,Complex xpp,Complex xmm) {
    const double lambda = l * (l + 1.0);
    const double zeta = std::sqrt(lambda);
    const double eta = l >= 2 ? std::sqrt((l - 1.0) * (l + 2.0)) : 0.0;
    return {a.Up*x0, (2.0*a.U-lambda*a.V)/r*x0,
      -zeta/2.0*(a.Vp-(a.V-a.U)/r+Complex(0,1)*(a.Wp-a.W/r))*xp,
      zeta/2.0*(a.Vp-(a.V-a.U)/r-Complex(0,1)*(a.Wp-a.W/r))*xm,
      zeta*eta/r*(a.V+Complex(0,1)*a.W)*xpp,
      zeta*eta/r*(a.V-Complex(0,1)*a.W)*xmm};
  }

  static PoleStrainCoefficients poleStrainCoefficients(
      int degree,const RadialMode &mode,double radius,double radiusScale) {
    const auto raw0=radius>0.0?radialToStrain(degree,mode,radius,1.0,{},{},{},{}):SpinStrainField{};
    const auto rawP=radius>0.0?radialToStrain(degree,mode,radius,{},1.0,{},{},{}):SpinStrainField{};
    const auto rawM=radius>0.0?radialToStrain(degree,mode,radius,{},{},1.0,{},{}):SpinStrainField{};
    const auto rawPP=radius>0.0?radialToStrain(degree,mode,radius,{},{},{},1.0,{}):SpinStrainField{};
    const auto rawMM=radius>0.0?radialToStrain(degree,mode,radius,{},{},{},{},1.0):SpinStrainField{};
    const auto scaled0=scale(rScaledStrainFromModal(degree,mode.U,mode.Up,mode.V,mode.Vp,
        mode.W,mode.Wp,radius,1.0,{},{},{},{}),radiusScale);
    const auto scaledP=scale(rScaledStrainFromModal(degree,mode.U,mode.Up,mode.V,mode.Vp,
        mode.W,mode.Wp,radius,{},1.0,{},{},{}),radiusScale);
    const auto scaledM=scale(rScaledStrainFromModal(degree,mode.U,mode.Up,mode.V,mode.Vp,
        mode.W,mode.Wp,radius,{},{},1.0,{},{}),radiusScale);
    const auto scaledPP=scale(rScaledStrainFromModal(degree,mode.U,mode.Up,mode.V,mode.Vp,
        mode.W,mode.Wp,radius,{},{},{},1.0,{}),radiusScale);
    const auto scaledMM=scale(rScaledStrainFromModal(degree,mode.U,mode.Up,mode.V,mode.Vp,
        mode.W,mode.Wp,radius,{},{},{},{},1.0),radiusScale);
    return {{raw0.rr,raw0.horizontalTrace,rawP.radialPlus,rawM.radialMinus,
             rawPP.horizontalPlus,rawMM.horizontalMinus},
            {scaled0.rr,scaled0.horizontalTrace,scaledP.radialPlus,scaledM.radialMinus,
             scaledPP.horizontalPlus,scaledMM.horizontalMinus}};
  }

public:
  struct BatchResponseFootprint {
    std::uint64_t entries=0, bytes=0, maximumDofs=0;
  };

  /** Exact dense-response retention estimate for the batch cache. This
   * mirrors the forward RHS widths and exact receiver-row basis widths used
   * by solveFrequencyModes; no sparse factors are retained by that cache. */
  static BatchResponseFootprint estimateBatchResponseFootprint(
      Full1D::SEM &sem,SourceInfo::EarthquakeCMT &source,InputParameters &parameters,
      Complex z,std::uint64_t frequencyCount) {
    if(frequencyCount==0)throw std::invalid_argument("response-cache frequency count must be positive");
    const int lmax=parameters.lmax();
    const auto flags=resolveModeFlags(parameters.type(),parameters.lmin(),lmax);
    BatchResponseFootprint result;
    const auto add=[&](const Eigen::MatrixXcd &forwardRhs,const Eigen::MatrixXcd &basisRhs) {
      if(forwardRhs.rows()!=basisRhs.rows()||forwardRhs.cols()==0||basisRhs.cols()==0)
        throw std::runtime_error("cannot estimate invalid response-cache RHS dimensions");
      const std::uint64_t n=static_cast<std::uint64_t>(forwardRhs.rows());
      const std::uint64_t columns=static_cast<std::uint64_t>(forwardRhs.cols()+basisRhs.cols());
      if(n>std::numeric_limits<std::uint64_t>::max()/columns)
        throw std::overflow_error("response-cache dense matrix size overflows");
      const std::uint64_t values=n*columns;
      if(values>(std::numeric_limits<std::uint64_t>::max()-32)/48)
        throw std::overflow_error("response-cache byte estimate overflows");
      const std::uint64_t bytes=48*values+32;
      if(result.bytes>std::numeric_limits<std::uint64_t>::max()-bytes)
        throw std::overflow_error("response-cache total byte estimate overflows");
      result.bytes+=bytes;++result.entries;
      result.maximumDofs=std::max(result.maximumDofs,n);
    };
    for(int degree=std::max(1,parameters.lmin());degree<=lmax;++degree) {
      if(flags.inc_sph) {
        const auto treatment=sem.sourceTreatedSpheroidal(source,degree,z);
        const auto receivers=sem.receiverElements(parameters,false);
        const auto receiverBase=sem.rvBaseFull(parameters,degree,false);
        const auto receiver=embedReceiver(sem,receiverBase,sem.ltgS(0,receivers.front(),0),
            static_cast<int>(sem.hS(degree).rows()));
        const auto basis=dsfast::exactReceiverRowBasis(receiver);
        add(treatment.force,basis.rightHandSides);
      }
      if(flags.inc_tor&&!sem.elementIsFluid(sem.sourceElement(source))) {
        const auto treatment=sem.sourceTreatedToroidal(source,degree,z);
        const auto receivers=sem.receiverElements(parameters,true);
        const auto receiverBase=sem.rvBaseFullT(parameters,degree);
        const auto receiver=embedReceiver(sem,receiverBase,sem.ltgT(receivers.front(),0),
            static_cast<int>(sem.hTk(degree).rows()));
        const auto basis=dsfast::exactReceiverRowBasis(receiver);
        add(treatment.force,basis.rightHandSides);
      }
    }
    if(result.entries>std::numeric_limits<std::uint64_t>::max()/frequencyCount||
       result.bytes>std::numeric_limits<std::uint64_t>::max()/frequencyCount)
      throw std::overflow_error("response-cache frequency-scaled estimate overflows");
    result.entries*=frequencyCount;result.bytes*=frequencyCount;
    return result;
  }

private:
  static void solveFrequencyModes(Full1D::SEM &sem,
      SourceInfo::EarthquakeCMT &source,InputParameters &parameters,Complex z,
      const Eigen::VectorXcd &adjointWeights,const ModalVisitor &visitor,
      std::vector<PoleModalBlock> *poleBlocks=nullptr,
      PoleFrame *sourceFrame=nullptr,
      dssolver::BatchResponseCache *batchCache=nullptr,bool poleOnly=false,
      bool forwardOnly=false,double *basisReconstructionRelative=nullptr,int solveThreads=1,std::optional<Complex> attenuation=std::nullopt) {
    mathprofile::Scope total("modal_projection_source_lift_and_frequency_solve");
    const char *restrictionEnv=std::getenv("DSPEC_PHASE_B_RESTRICTED");
    const bool restrictDomain=restrictionEnv&&std::string(restrictionEnv)=="1";
    if(restrictDomain&&solveThreads!=1)
      throw std::invalid_argument("Phase B restriction requires one inner solve thread per frequency worker.");
    const auto sourceElement=sem.sourceElement(source);
    const auto receiverElements=sem.receiverElements(parameters,false);
    if(receiverElements.empty())throw std::invalid_argument("Phase B restriction requires a receiver element.");
    FixedMeshAdjointSolver::ScopedRestrictionRequest restrictionScope(
        restrictDomain,sourceElement,receiverElements.front());
    const Complex dispersion=attenuation.value_or(Complex{});
    if((parameters.attenuation()&&!attenuation)||!std::isfinite(dispersion.real())||!std::isfinite(dispersion.imag()))throw std::invalid_argument("explicit finite dispersion required for attenuated kernels");
    if(solveThreads<1||solveThreads>48)
      throw std::invalid_argument("solveThreads must be between 1 and 48");
    if(solveThreads>1&&(!poleOnly&&!forwardOnly))
      throw std::invalid_argument("parallel modal solves require the pole-frame or forward-only route");
    if(solveThreads>1&&omp_in_parallel())
      throw std::invalid_argument("nested modal solve parallelism is not supported");

    mathprofile::count("modal_lmax",static_cast<std::uint64_t>(parameters.lmax()));
    mathprofile::count("modal_receiver_rhs_components",static_cast<std::uint64_t>(adjointWeights.size()));
    mathprofile::count("full_weighted_measurement_adjoint_components",static_cast<std::uint64_t>(adjointWeights.size()));
    mathprofile::count("pole_only_nonradial_weighted_adjoint_route",poleOnly?1:0);
    mathprofile::count("true_forward_only_route",forwardOnly?1:0);
    if(poleOnly||forwardOnly)mathprofile::count("global_paraminfo_setup_skipped");
    if(!sem.sourceMatches(source)||
       parameters.cowlingFrequencyMhz()!=0.0||parameters.gravitation()!=2||
       parameters.output_type()!=0||!std::isfinite(z.real())||
       !std::isfinite(z.imag())||z==Complex{}||
       adjointWeights.size()!=3*parameters.num_receivers()||!adjointWeights.allFinite())
      throw std::invalid_argument("Invalid one-frequency 3D modal workspace input.");
    const int lmax=parameters.lmax();
    if(basisReconstructionRelative)*basisReconstructionRelative=0.0;
    if(parameters.lmin()<0||lmax<parameters.lmin()||lmax>sem.maxDegree())
      throw std::invalid_argument("3D modal degree range is invalid.");
    const auto flags=resolveModeFlags(parameters.type(),parameters.lmin(),lmax);
    if(!flags.inc_rad&&!flags.inc_sph&&!flags.inc_tor)
      throw std::invalid_argument("3D field mode selection is empty.");
    const Complex fz=-Complex(0.0,1.0)/z;
    const char *directAdjointEnv=std::getenv("DSPEC_PHASE_C_DIRECT_ADJOINT");
    const bool directAdjointRoute=directAdjointEnv&&std::string(directAdjointEnv)=="1";
    const char *directValidationEnv=std::getenv("DSPEC_PHASE_C_VALIDATE_DIRECT_ADJOINT");
    const bool validateDirectAdjoint=directValidationEnv&&std::string(directValidationEnv)=="1";
    if(directAdjointRoute&&!poleOnly)
      throw std::invalid_argument("Phase C direct measurement adjoint requires pole-frame workspace mode.");
    if(sourceFrame)*sourceFrame=poleFrame(source.Latitude(),source.Longitude());
    std::unique_ptr<ParamInfo> angular;
    if(!poleOnly&&!forwardOnly){
      mathprofile::Scope t("global_receiver_paraminfo_harmonic_setup");
      angular=std::make_unique<ParamInfo>(parameters,lmax);
    }
    if(flags.inc_rad){
      mathprofile::Scope radial("radial_l0_source_and_adjoint_branch");
      auto treatment=sem.sourceTreatedRadial(source,z,dispersion);
      Eigen::MatrixXcd forward,adjoint,unusedBasis;
      if(forwardOnly){
        Eigen::MatrixXcd empty;
        FixedMeshAdjointSolver::solveForwardAdjointFrequency(
            sem,SensitivitySolveFamily::Radial,0,z,treatment.force,empty,empty,
            forward,adjoint,unusedBasis,dispersion);
      }else{
        Eigen::MatrixXcd receiver=Eigen::MatrixXcd::Zero(
            3*parameters.num_receivers(),treatment.force.rows());
        for(int rec=0;rec<parameters.num_receivers();++rec)
          receiver.row(3*rec)=sem.rvZR(parameters,rec).transpose();
        const Eigen::MatrixXcd rhs=receiver.adjoint()*(std::conj(fz)*adjointWeights);
        mathprofile::count("full_weighted_modal_adjoint_rhs_columns",static_cast<std::uint64_t>(rhs.cols()));
        if(batchCache)mathprofile::count("response_cache_radial_l0_bypass_calls");
        FixedMeshAdjointSolver::solveForwardAdjointFrequency(
            sem,SensitivitySolveFamily::Radial,0,z,treatment.force,rhs,{},
            forward,adjoint,unusedBasis,dispersion);
      }
      const Eigen::MatrixXcd noData;
      visitor({SensitivitySolveFamily::Radial,0,treatment.sourceElement,
          forward,adjoint,treatment.lift,noData,noData});
    }
    const auto solveDegree=[&](int degree,std::vector<PoleModalBlock> *poleBlocks,
                                double *basisReconstructionRelative){
      mathprofile::count("modal_degree_iterations");
      const int orders=2*degree+1;
      if(flags.inc_sph&&poleOnly)mathprofile::count("pole_only_spheroidal_weighted_rhs_columns_avoided",static_cast<std::uint64_t>(orders));
      if(flags.inc_sph&&forwardOnly)mathprofile::count("forward_only_spheroidal_unit_adjoint_modes_skipped");
      if(flags.inc_sph&&poleOnly)mathprofile::count("pole_only_spheroidal_global_harmonic_setup_skips");
      if(flags.inc_sph&&forwardOnly)mathprofile::count("forward_only_spheroidal_global_harmonic_setup_skips");
      if(flags.inc_sph){
        mathprofile::Scope sph("spheroidal_source_lift_projection_and_solver");
        mathprofile::Scope sourceTreatment("spheroidal_source_treatment_and_lift");
        auto treatment=sem.sourceTreatedSpheroidal(source,degree,z,dispersion);
        sourceTreatment.stop();
        Eigen::MatrixXcd adjoint,adjointUnit;
        Eigen::MatrixXcd forward;
        Eigen::MatrixXcd receiverBasisRhs,reconstruction,rhs,modeReceiver,sourceCoefficients;
        std::vector<Eigen::MatrixXcd> directAdjointByReceiver;
        if(!forwardOnly){
          mathprofile::Scope receiverBasisBuild(directAdjointRoute&&!validateDirectAdjoint?
              "direct_measurement_receiver_map_assembly":"spheroidal_receiver_basis_assembly_and_qr");
          const auto receivers=sem.receiverElements(parameters,false);
          const auto receiverBase=sem.rvBaseFull(parameters,degree,false);
          const auto receiver=embedReceiver(sem,receiverBase,
              sem.ltgS(0,receivers.front(),0),static_cast<int>(sem.hS(degree).rows()));
          dsfast::ReceiverRowBasis receiverBasis;
          if(!directAdjointRoute||validateDirectAdjoint)
            receiverBasis=dsfast::exactReceiverRowBasis(receiver);
          receiverBasisBuild.stop();
          if(!directAdjointRoute){
            if(basisReconstructionRelative)*basisReconstructionRelative=std::max(
                *basisReconstructionRelative,receiverBasis.reconstructionRelative);
            receiverBasisRhs=receiverBasis.rightHandSides;
            reconstruction=receiverBasis.reconstruction;
          }
          if(!poleOnly){
            {mathprofile::Scope t("global_receiver_harmonic_coefficients");modeReceiver=angular->rvFullSph(degree);}
            {mathprofile::Scope t("global_source_harmonic_coefficients");sourceCoefficients=sem.calculateForceCoefficients(source,degree);}
            rhs=Eigen::MatrixXcd::Zero(receiver.cols(),orders);
            mathprofile::count("full_weighted_modal_adjoint_rhs_columns",static_cast<std::uint64_t>(rhs.cols()));
            for(int m=-degree;m<=degree;++m){
              Eigen::MatrixXcd weighted=receiver;
              for(int row=0;row<weighted.rows();++row)weighted.row(row)*=modeReceiver(row,m+degree);
              rhs.col(m+degree)=weighted.adjoint()*(std::conj(fz)*adjointWeights);
            }
          }
          if(directAdjointRoute){
            const auto poleReceiver=receiverAtPole(SensitivitySolveFamily::Spheroidal,degree);
            const Eigen::Index receiverCount=parameters.num_receivers();
            const Eigen::Index orderCount=poleReceiver.cols();
            Eigen::MatrixXcd localWeights=Eigen::MatrixXcd::Zero(3*receiverCount,
                receiverCount*orderCount);
            for(Eigen::Index rec=0;rec<receiverCount;++rec)
              for(int comp=0;comp<3;++comp)
                for(Eigen::Index m=0;m<orderCount;++m)
                  localWeights(3*rec+comp,rec*orderCount+m)=
                      std::conj(poleReceiver(comp,m))*(std::conj(fz)*adjointWeights[3*rec+comp]);
            const Eigen::MatrixXcd directRhs=receiver.adjoint()*localWeights;
            Eigen::MatrixXcd basisTransform,basisRhs;
            if(validateDirectAdjoint){
              basisRhs=receiverBasis.rightHandSides;
              basisTransform=receiverBasis.reconstruction.adjoint()*localWeights;
              const Eigen::MatrixXcd reconstructedRhs=basisRhs*basisTransform;
              const double rhsError=(directRhs-reconstructedRhs).norm()/
                  std::max(directRhs.norm(),1e-300);
              mathprofile::count("phase_c_direct_rhs_equivalence_scaled_1e12",
                  static_cast<std::uint64_t>(std::min(1e12,rhsError*1e12)));
              if(!std::isfinite(rhsError)||rhsError>1e-10)
                throw std::runtime_error("Phase C direct RHS differs from the accepted receiver-basis map.");
            }
            Eigen::MatrixXcd directSolution,basisSolution;
            FixedMeshAdjointSolver::solveForwardAdjointFrequency(
                sem,SensitivitySolveFamily::Spheroidal,degree,z,treatment.force,directRhs,
                basisRhs,forward,directSolution,basisSolution,dispersion);
            if(validateDirectAdjoint){
              const double solutionError=(directSolution-basisSolution*basisTransform).norm()/
                  std::max(directSolution.norm(),1e-300);
              mathprofile::count("phase_c_direct_solution_equivalence_scaled_1e12",
                  static_cast<std::uint64_t>(std::min(1e12,solutionError*1e12)));
              if(!std::isfinite(solutionError)||solutionError>1e-9)
                throw std::runtime_error("Phase C direct adjoint differs from the accepted receiver-basis projection.");
            }
            directAdjointByReceiver.reserve(static_cast<std::size_t>(receiverCount));
            for(Eigen::Index rec=0;rec<receiverCount;++rec)
              directAdjointByReceiver.push_back(directSolution.middleCols(rec*orderCount,orderCount));
            mathprofile::count("phase_c_direct_adjoint_rhs_columns",static_cast<std::uint64_t>(directRhs.cols()));
          }
        }
        if(forwardOnly){
          Eigen::MatrixXcd empty;
          FixedMeshAdjointSolver::solveForwardAdjointFrequency(
              sem,SensitivitySolveFamily::Spheroidal,degree,z,treatment.force,empty,empty,
              forward,adjoint,adjointUnit,dispersion);
        }else if(directAdjointRoute){
          // Direct H^H solutions were formed above; do not consult or populate
          // the unit-basis response cache on this path.
        }else if(batchCache){
          const auto key=dssolver::BatchResponseKey::make(&sem,
              static_cast<int>(SensitivitySolveFamily::Spheroidal),degree,z,dispersion);
          const auto response=batchCache->getOrCompute(key,treatment.force,
              receiverBasisRhs,[&](const auto &forwardRhs,const auto &basisRhs){
                Eigen::MatrixXcd f,a,b,empty(0,0);
                FixedMeshAdjointSolver::solveForwardAdjointFrequency(
                    sem,SensitivitySolveFamily::Spheroidal,degree,z,forwardRhs,empty,basisRhs,f,a,b,dispersion);
                return std::pair{std::move(f),std::move(b)};
              });
          if(poleOnly){
            adjointUnit=response.unitAdjoint;
          }else{
            const auto weighted=batchCache->deriveWeightedAdjoint(key,receiverBasisRhs,rhs);
            adjoint=weighted.weightedAdjoint;adjointUnit=response.unitAdjoint;
          }
          forward=response.forward;
        }else FixedMeshAdjointSolver::solveForwardAdjointFrequency(
            sem,SensitivitySolveFamily::Spheroidal,degree,z,treatment.force,rhs,
            receiverBasisRhs,forward,adjoint,adjointUnit,dispersion);
        if(!poleOnly&&!forwardOnly)visitor({SensitivitySolveFamily::Spheroidal,degree,treatment.sourceElement,
          forward,adjoint,treatment.lift,sourceCoefficients,modeReceiver});
        if(poleBlocks){
          const auto reduced=sem.calculateForceRedCoefficients(source,degree,0.0);
          poleBlocks->push_back({SensitivitySolveFamily::Spheroidal,degree,
              treatment.sourceElement,forward,adjointUnit,treatment.lift,reduced,reconstruction});
          if(directAdjointRoute)
            poleBlocks->back().adjointDirect=std::move(directAdjointByReceiver);
        }
      }
      if(flags.inc_tor&&poleOnly)mathprofile::count("pole_only_toroidal_weighted_rhs_columns_avoided",static_cast<std::uint64_t>(orders));
      if(flags.inc_tor&&forwardOnly)mathprofile::count("forward_only_toroidal_unit_adjoint_modes_skipped");
      if(flags.inc_tor&&poleOnly)mathprofile::count("pole_only_toroidal_global_harmonic_setup_skips");
      if(flags.inc_tor&&forwardOnly)mathprofile::count("forward_only_toroidal_global_harmonic_setup_skips");
      if(flags.inc_tor&&!sem.elementIsFluid(sem.sourceElement(source))){
        mathprofile::Scope tor("toroidal_source_lift_projection_and_solver");
        mathprofile::Scope sourceTreatment("toroidal_source_treatment_and_lift");
        auto treatment=sem.sourceTreatedToroidal(source,degree,z,dispersion);
        sourceTreatment.stop();
        Eigen::MatrixXcd adjoint,adjointUnit;
        Eigen::MatrixXcd forward;
        Eigen::MatrixXcd receiverBasisRhs,reconstruction,rhs,modeReceiver,sourceCoefficients;
        std::vector<Eigen::MatrixXcd> directAdjointByReceiver;
        if(!forwardOnly){
          mathprofile::Scope receiverBasisBuild(directAdjointRoute&&!validateDirectAdjoint?
              "direct_measurement_receiver_map_assembly":"toroidal_receiver_basis_assembly_and_qr");
          const auto receivers=sem.receiverElements(parameters,true);
          const auto receiverBase=sem.rvBaseFullT(parameters,degree);
          const auto receiver=embedReceiver(sem,receiverBase,
              sem.ltgT(receivers.front(),0),static_cast<int>(sem.hTk(degree).rows()));
          dsfast::ReceiverRowBasis receiverBasis;
          if(!directAdjointRoute||validateDirectAdjoint)
            receiverBasis=dsfast::exactReceiverRowBasis(receiver);
          receiverBasisBuild.stop();
          if(!directAdjointRoute){
            if(basisReconstructionRelative)*basisReconstructionRelative=std::max(
                *basisReconstructionRelative,receiverBasis.reconstructionRelative);
            receiverBasisRhs=receiverBasis.rightHandSides;
            reconstruction=receiverBasis.reconstruction;
          }
          if(!poleOnly){
            {mathprofile::Scope t("global_receiver_harmonic_coefficients");modeReceiver=angular->rvFullTor(degree);}
            {mathprofile::Scope t("global_source_harmonic_coefficients");sourceCoefficients=sem.calculateForceCoefficientsT(source,degree);}
            rhs=Eigen::MatrixXcd::Zero(receiver.cols(),orders);
            mathprofile::count("full_weighted_modal_adjoint_rhs_columns",static_cast<std::uint64_t>(rhs.cols()));
            for(int m=-degree;m<=degree;++m){
              Eigen::MatrixXcd weighted=receiver;
              for(int row=0;row<weighted.rows();++row)weighted.row(row)*=modeReceiver(row,m+degree);
              rhs.col(m+degree)=weighted.adjoint()*(std::conj(fz)*adjointWeights);
            }
          }
          if(directAdjointRoute){
            const auto poleReceiver=receiverAtPole(SensitivitySolveFamily::Toroidal,degree);
            const Eigen::Index receiverCount=parameters.num_receivers();
            const Eigen::Index orderCount=poleReceiver.cols();
            Eigen::MatrixXcd localWeights=Eigen::MatrixXcd::Zero(3*receiverCount,
                receiverCount*orderCount);
            for(Eigen::Index rec=0;rec<receiverCount;++rec)
              for(int comp=0;comp<3;++comp)
                for(Eigen::Index m=0;m<orderCount;++m)
                  localWeights(3*rec+comp,rec*orderCount+m)=
                      std::conj(poleReceiver(comp,m))*(std::conj(fz)*adjointWeights[3*rec+comp]);
            const Eigen::MatrixXcd directRhs=receiver.adjoint()*localWeights;
            Eigen::MatrixXcd basisTransform,basisRhs;
            if(validateDirectAdjoint){
              basisRhs=receiverBasis.rightHandSides;
              basisTransform=receiverBasis.reconstruction.adjoint()*localWeights;
              const Eigen::MatrixXcd reconstructedRhs=basisRhs*basisTransform;
              const double rhsError=(directRhs-reconstructedRhs).norm()/
                  std::max(directRhs.norm(),1e-300);
              mathprofile::count("phase_c_direct_rhs_equivalence_scaled_1e12",
                  static_cast<std::uint64_t>(std::min(1e12,rhsError*1e12)));
              if(!std::isfinite(rhsError)||rhsError>1e-10)
                throw std::runtime_error("Phase C direct RHS differs from the accepted receiver-basis map.");
            }
            Eigen::MatrixXcd directSolution,basisSolution;
            FixedMeshAdjointSolver::solveForwardAdjointFrequency(
                sem,SensitivitySolveFamily::Toroidal,degree,z,treatment.force,directRhs,
                basisRhs,forward,directSolution,basisSolution,dispersion);
            if(validateDirectAdjoint){
              const double solutionError=(directSolution-basisSolution*basisTransform).norm()/
                  std::max(directSolution.norm(),1e-300);
              mathprofile::count("phase_c_direct_solution_equivalence_scaled_1e12",
                  static_cast<std::uint64_t>(std::min(1e12,solutionError*1e12)));
              if(!std::isfinite(solutionError)||solutionError>1e-9)
                throw std::runtime_error("Phase C direct adjoint differs from the accepted receiver-basis projection.");
            }
            directAdjointByReceiver.reserve(static_cast<std::size_t>(receiverCount));
            for(Eigen::Index rec=0;rec<receiverCount;++rec)
              directAdjointByReceiver.push_back(directSolution.middleCols(rec*orderCount,orderCount));
            mathprofile::count("phase_c_direct_adjoint_rhs_columns",static_cast<std::uint64_t>(directRhs.cols()));
          }
        }
        if(forwardOnly){
          Eigen::MatrixXcd empty;
          FixedMeshAdjointSolver::solveForwardAdjointFrequency(
              sem,SensitivitySolveFamily::Toroidal,degree,z,treatment.force,empty,empty,
              forward,adjoint,adjointUnit,dispersion);
        }else if(directAdjointRoute){
          // Direct H^H solutions were formed above; do not consult or populate
          // the unit-basis response cache on this path.
        }else if(batchCache){
          const auto key=dssolver::BatchResponseKey::make(&sem,
              static_cast<int>(SensitivitySolveFamily::Toroidal),degree,z,dispersion);
          const auto response=batchCache->getOrCompute(key,treatment.force,
              receiverBasisRhs,[&](const auto &forwardRhs,const auto &basisRhs){
                Eigen::MatrixXcd f,a,b,empty(0,0);
                FixedMeshAdjointSolver::solveForwardAdjointFrequency(
                    sem,SensitivitySolveFamily::Toroidal,degree,z,forwardRhs,empty,basisRhs,f,a,b,dispersion);
                return std::pair{std::move(f),std::move(b)};
              });
          if(poleOnly){
            adjointUnit=response.unitAdjoint;
          }else{
            const auto weighted=batchCache->deriveWeightedAdjoint(key,receiverBasisRhs,rhs);
            adjoint=weighted.weightedAdjoint;adjointUnit=response.unitAdjoint;
          }
          forward=response.forward;
        }else FixedMeshAdjointSolver::solveForwardAdjointFrequency(
            sem,SensitivitySolveFamily::Toroidal,degree,z,treatment.force,rhs,
            receiverBasisRhs,forward,adjoint,adjointUnit,dispersion);
        if(!poleOnly&&!forwardOnly)visitor({SensitivitySolveFamily::Toroidal,degree,treatment.sourceElement,
          forward,adjoint,treatment.lift,sourceCoefficients,modeReceiver});
        if(poleBlocks){
          const auto reduced=sem.calculateForceRedCoefficientsT(source,degree,0.0);
          poleBlocks->push_back({SensitivitySolveFamily::Toroidal,degree,
              treatment.sourceElement,forward,adjointUnit,treatment.lift,reduced,reconstruction});
          if(directAdjointRoute)
            poleBlocks->back().adjointDirect=std::move(directAdjointByReceiver);
        }
      }
    };
    const int firstDegree=std::max(1,parameters.lmin());
    if(solveThreads==1){
      for(int degree=firstDegree;degree<=lmax;++degree)
        solveDegree(degree,poleBlocks,basisReconstructionRelative);
    }else{
      // SEM matrices, mesh, source and input getters are read-only on this
      // route. Sparse factors and profiling state are thread-local. Each
      // degree owns one result slot; no worker invokes a shared visitor.
      struct DegreeResult {
        std::vector<PoleModalBlock> blocks;
        double residual=0.0;
        std::exception_ptr error;
      };
      const int workers=std::min(solveThreads,std::max(1,lmax-firstDegree+1));
      std::vector<DegreeResult> results(static_cast<std::size_t>(lmax+1));
      std::vector<mathprofile::WorkerData> profiles(static_cast<std::size_t>(workers));
      const bool profileWorkers=mathprofile::active;
      const int backend=FixedMeshAdjointSolver::diagnosticBackend();
      mathprofile::Scope parallelRegion("modal_degree_parallel_region_wall");
      #pragma omp parallel num_threads(workers)
      {
        mathprofile::WorkerCapture profile(profileWorkers,profiles.at(omp_get_thread_num()));
        const int previousBackend=FixedMeshAdjointSolver::diagnosticBackend();
        FixedMeshAdjointSolver::setDiagnosticBackend(backend);
        #pragma omp for schedule(static)
        for(int degree=firstDegree;degree<=lmax;++degree){
          auto &result=results[static_cast<std::size_t>(degree)];
          try {
            solveDegree(degree,poleBlocks?&result.blocks:nullptr,
                        basisReconstructionRelative?&result.residual:nullptr);
          }catch(...){result.error=std::current_exception();}
        }
        FixedMeshAdjointSolver::setDiagnosticBackend(previousBackend);
      }
      parallelRegion.stop();
      mathprofile::Scope profileMerge("modal_degree_worker_profile_merge");
      for(const auto &profile:profiles)mathprofile::mergeWorker(profile);
      profileMerge.stop();
      // Ordered collection preserves the serial harmonic accumulation and
      // spool order. Exceptions cross the OpenMP boundary only after joining.
      mathprofile::Scope orderedCollection("modal_degree_ordered_result_collection");
      for(int degree=firstDegree;degree<=lmax;++degree){
        auto &result=results[static_cast<std::size_t>(degree)];
        if(result.error)std::rethrow_exception(result.error);
        if(basisReconstructionRelative)*basisReconstructionRelative=std::max(
            *basisReconstructionRelative,result.residual);
        if(poleBlocks)for(auto &block:result.blocks)poleBlocks->push_back(std::move(block));
      }
      orderedCollection.stop();
    }
  }

  static void accumulateModalFields(const Full1D::SEM &sem,
      const ModalFrequencyBlock3D &block,const std::vector<Point> &points,
      const std::vector<double> &radii,const std::vector<int> &elements,
      const std::vector<std::vector<std::size_t>> &radialGroups,
      double radiusScale,std::vector<Sample> &result,const AngularWigner *w,
      std::size_t angleBegin,const std::vector<std::size_t> &angleForPoint,
      bool kernelRoute=false,std::size_t displacementBegin=0) {
    const int l=block.degree;
    if(block.family==SensitivitySolveFamily::Radial){
      const Complex y00=harmonic(0,0,0,0.0,0.0);
      for(const auto &group:radialGroups){
        const auto p=group.front();
        const auto fRadial=interpolate(sem,block.family,block.forward,
            block.sourceLift,radii[p],elements[p],0,block.sourceElement);
        Eigen::MatrixXcd aRadial=Eigen::MatrixXcd::Zero(2,1);
        if(block.adjoint.cols()!=0)aRadial=interpolate(sem,block.family,block.adjoint,
            Eigen::MatrixXcd{},radii[p],elements[p]);
        Eigen::MatrixXcd f=Eigen::MatrixXcd::Zero(6,1),a=f;
        f(0,0)=fRadial(0,0);f(1,0)=fRadial(1,0);
        a(0,0)=aRadial(0,0);a(1,0)=aRadial(1,0);
        for(const auto q:group){
          if(!kernelRoute||q>=displacementBegin)
            add(result[q].forwardRawDisplacement,radialToDisplacement(0,f,1.0,y00,{},{},{},{}));
          if(!kernelRoute)
            add(result[q].adjointDualDisplacement,radialToDisplacement(0,a,1.0,y00,{},{},{},{}));
          if(radii[q]>0.0){
            add(result[q].forwardRawStrain,radialToStrain(0,f,radii[q],y00,{},{},{},{}));
            add(result[q].adjointDualStrain,radialToStrain(0,a,radii[q],y00,{},{},{},{}));}
          if(!kernelRoute){
            add(result[q].rScaledForwardStrainMetres,scale(rScaledStrainFromModal(
                0,f(0,0),f(1,0),{},{},{},{},radii[q],y00,{},{},{},{}),radiusScale));
            add(result[q].rScaledAdjointDualStrainMetres,scale(rScaledStrainFromModal(
                0,a(0,0),a(1,0),{},{},{},{},radii[q],y00,{},{},{},{}),radiusScale));
          }
        }
      }
      return;
    }
    if(w==nullptr||angleForPoint.size()!=points.size())
      throw std::logic_error("Missing angular table for 3D field block.");
    for(const auto &group:radialGroups){
      const auto representative=group.front();
      if(block.family==SensitivitySolveFamily::Toroidal&&
         (sem.elementIsFluid(elements[representative])||elements[representative]<sem.el()||elements[representative]>=sem.eu()))continue;
      for(int m=-l;m<=l;++m){
        const auto fmode=interpolateCombination(sem,block.family,block.forward,
            block.sourceLift,block.sourceCoefficients.row(m+l).transpose(),
            radii[representative],elements[representative],block.sourceElement);
        const auto amode=interpolate(sem,block.family,block.adjoint,
            Eigen::MatrixXcd{},radii[representative],elements[representative],m+l);
        for(const auto p:group){
          const Complex phase=std::exp(Complex(0.0,m*points[p].phiRadians));
          const std::size_t thetaIndex=angleForPoint[p]-angleBegin;
          const Complex x0=block.family==SensitivitySolveFamily::Spheroidal?
              harmonic(*w,l,m,0,thetaIndex,phase):Complex{};
          const Complex xp=harmonic(*w,l,m,1,thetaIndex,phase);
          const Complex xm=harmonic(*w,l,m,-1,thetaIndex,phase);
          const Complex xpp=l>=2?harmonic(*w,l,m,2,thetaIndex,phase):Complex{};
          const Complex xmm=l>=2?harmonic(*w,l,m,-2,thetaIndex,phase):Complex{};
          const auto fd=block.family==SensitivitySolveFamily::Spheroidal?
              radialToDisplacement(l,fmode,radii[p],x0,xp,xm,xpp,xmm):
              toroidalDisplacement(l,fmode,xp,xm);
          add(result[p].forwardRawDisplacement,fd);
          if(radii[p]>0.0)add(result[p].forwardRawStrain,radialToStrain(
              l,fmode,radii[p],x0,xp,xm,xpp,xmm));
          const auto fm=mode(fmode);
          add(result[p].rScaledForwardStrainMetres,scale(rScaledStrainFromModal(
              l,fm.U,fm.Up,fm.V,fm.Vp,fm.W,fm.Wp,radii[p],x0,xp,xm,xpp,xmm),radiusScale));
          const auto ad=block.family==SensitivitySolveFamily::Spheroidal?
              radialToDisplacement(l,amode,radii[p],x0,xp,xm,xpp,xmm):
              toroidalDisplacement(l,amode,xp,xm);
          add(result[p].adjointDualDisplacement,ad);
          if(radii[p]>0.0)add(result[p].adjointDualStrain,radialToStrain(
              l,amode,radii[p],x0,xp,xm,xpp,xmm));
          const auto am=mode(amode);
          add(result[p].rScaledAdjointDualStrainMetres,scale(rScaledStrainFromModal(
              l,am.U,am.Up,am.V,am.Vp,am.W,am.Wp,radii[p],x0,xp,xm,xpp,xmm),radiusScale));
        }
      }
    }
  }

  static void synthesizePoleFrame(const Full1D::SEM &sem,
      const std::vector<PoleSynthesisBlock> &blocks,const PoleFrame &frame,
      const std::vector<Point> &points,const std::vector<double> &radii,
      const std::vector<int> &elements,double radiusScale,
      std::vector<Sample> &result,bool kernelRoute=false,
      std::size_t displacementBegin=0,SamplingPlanCache *samplingPlan=nullptr,
      PoleCoefficientCache *frequencyCoefficientCache=nullptr,
      std::uint64_t synthesisSet=0,int contractionBackend=0,
      const std::function<void(const PoleAngularCoefficientRecord3D&)> &coefficientVisitor={}) {
    auditprofile::Scope whole("pole_synthesis");
#ifndef DSFAST_DISABLE_PROFILE
    auto geometryStart=std::chrono::steady_clock::now();
#endif
    if(blocks.empty())return;
    using PointAngle=SamplingPlanCache::PointAngle;
    SamplingPlanCache::FrameKey frameKey{};
    for(int axis=0;axis<3;++axis)for(int coordinate=0;coordinate<3;++coordinate)
      frameKey[static_cast<std::size_t>(3*axis+coordinate)]=
          std::bit_cast<std::uint64_t>(frame.axes[axis][coordinate]);
    std::vector<PointAngle> angles(points.size());
    std::map<std::uint64_t,std::size_t> angleMap;
    std::vector<double> uniqueTheta;
    std::vector<std::size_t> angleForPoint(points.size());
    for(std::size_t p=0;p<points.size();++p){
      const SamplingPlanCache::PointKey pointKey{
          std::bit_cast<std::uint64_t>(points[p].thetaRadians),
          std::bit_cast<std::uint64_t>(points[p].phiRadians)};
      SamplingPlanCache::PointAngle angle;
      bool found=false;
      if(samplingPlan){
        const auto frameIt=samplingPlan->pointAngles.find(frameKey);
        if(frameIt!=samplingPlan->pointAngles.end()){
          const auto pointIt=frameIt->second.find(pointKey);
          if(pointIt!=frameIt->second.end()){angle=pointIt->second;found=true;++auditprofile::counters["geometry_hit"];}
        }
      }
      if(!found){++auditprofile::counters["geometry_miss"];
        double a=1.0,b=0.0;
        const auto local=tangentRotation(frame,points[p].thetaRadians,points[p].phiRadians,a,b);
        angle.theta=local[0];angle.phi=local[1];angle.a=a;angle.b=b;
        for(int m=-2;m<=2;++m)
          angle.orderPhases[static_cast<std::size_t>(m+2)]=std::exp(Complex(0.0,m*local[1]));
        constexpr std::size_t estimatedAngleEntryBytes=256;
        if(samplingPlan&&samplingPlan->bytes+estimatedAngleEntryBytes<=SamplingPlanCache::byteBudget){
          samplingPlan->bytes+=estimatedAngleEntryBytes;
          samplingPlan->pointAngles[frameKey].emplace(pointKey,angle);
        }
      }
      angles[p]=angle;
      const auto [it,added]=angleMap.try_emplace(std::bit_cast<std::uint64_t>(angle.theta),uniqueTheta.size());
      if(added)uniqueTheta.push_back(angle.theta);
      angleForPoint[p]=it->second;
    }
    std::vector<std::vector<std::size_t>> radialGroups;
    std::map<std::pair<int,double>,std::size_t> groupIndex;
    for(std::size_t p=0;p<points.size();++p){
      const auto [it,added]=groupIndex.try_emplace({elements[p],radii[p]},radialGroups.size());
      if(added)radialGroups.emplace_back();
      radialGroups[it->second].push_back(p);
    }
    int maxDegree=0;
    for(const auto &block:blocks)maxDegree=std::max(maxDegree,block.degree);
    constexpr std::size_t tableBudget=32u*1024u*1024u;
    const std::size_t perAngle=angularTableBytesPerAngle(maxDegree,std::min(2,maxDegree));
    const std::size_t angleChunkSize=std::max<std::size_t>(1,
        std::min<std::size_t>(64,tableBudget/std::max<std::size_t>(1,perAngle)));
    const std::size_t chunkCount=(uniqueTheta.size()+angleChunkSize-1)/angleChunkSize;
    std::vector<AngularPointChunk> chunks(chunkCount);
    for(std::size_t c=0;c<chunkCount;++c){
      const std::size_t begin=c*angleChunkSize;
      const std::size_t end=std::min(uniqueTheta.size(),begin+angleChunkSize);
      chunks[c].theta.assign(uniqueTheta.begin()+begin,uniqueTheta.begin()+end);
    }
    for(std::size_t g=0;g<radialGroups.size();++g)for(const auto p:radialGroups[g]){
      auto &chunk=chunks[angleForPoint[p]/angleChunkSize];
      if(chunk.lastRadialGroup!=g){chunk.radialGroups.emplace_back();chunk.lastRadialGroup=g;}
      chunk.radialGroups.back().push_back(p);
    }
    std::vector<std::size_t> pointForScratch(points.size(),
        std::numeric_limits<std::size_t>::max());
#ifndef DSFAST_DISABLE_PROFILE
    auditprofile::data["pole_geometry_grouping"].seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-geometryStart).count();
    ++auditprofile::data["pole_geometry_grouping"].calls;
#endif
    for(std::size_t chunk=0;chunk<chunkCount;++chunk){
      const std::size_t begin=chunk*angleChunkSize;
      std::vector<std::uint64_t> angularKey;
      angularKey.reserve(chunks[chunk].theta.size()+1);
      angularKey.push_back(static_cast<std::uint64_t>(maxDegree));
      for(const double theta:chunks[chunk].theta)
        angularKey.push_back(std::bit_cast<std::uint64_t>(theta));
      std::shared_ptr<CachedAngularWigner> angularTable;
      if(samplingPlan){
        const auto tableIt=samplingPlan->angularTables.find(angularKey);
        if(tableIt!=samplingPlan->angularTables.end()){angularTable=tableIt->second;++auditprofile::counters["wigner_hit"];mathprofile::count("packed_wigner_cache_hits");}
      }
      if(!angularTable){auditprofile::Scope t("wigner_build");mathprofile::Scope mt("packed_Wigner_d_table_build");mathprofile::count("packed_wigner_cache_misses");++auditprofile::counters["wigner_miss"];
        angularTable=std::make_shared<CachedAngularWigner>(
            maxDegree,std::min(2,maxDegree),std::min(2,maxDegree),chunks[chunk].theta);
        const std::size_t bytes=angularWignerBytesPerAngle(maxDegree,std::min(2,maxDegree))*chunks[chunk].theta.size();
        if(samplingPlan&&samplingPlan->bytes+bytes<=SamplingPlanCache::byteBudget){
          samplingPlan->bytes+=bytes;
          samplingPlan->angularTables.emplace(angularKey,angularTable);
        }
      }
      const auto &w=*angularTable;
      std::shared_ptr<AngularCoefficientTable> angularCoefficients;
      if(samplingPlan){
        const auto coefficientsIt=samplingPlan->angularCoefficients.find(angularKey);
        if(coefficientsIt!=samplingPlan->angularCoefficients.end())
          {angularCoefficients=coefficientsIt->second;++auditprofile::counters["angular_coefficient_hit"];mathprofile::count("packed_angular_coeff_cache_hits");}
      }
      if(!angularCoefficients){auditprofile::Scope t("angular_coefficient_build");mathprofile::Scope mt("packed_angular_coefficients_and_basis_build");mathprofile::count("packed_angular_coeff_cache_misses");++auditprofile::counters["angular_coefficient_miss"];
        auto prepared=std::make_shared<AngularCoefficientTable>();
        prepared->maxDegree=maxDegree;
        prepared->nMax=std::min(2,maxDegree);
        prepared->thetaCount=chunks[chunk].theta.size();
        const std::size_t stride=static_cast<std::size_t>(maxDegree+1);
        prepared->values.assign(prepared->thetaCount*
            static_cast<std::size_t>(2*prepared->nMax+1)*5*stride,0.0);
        for(std::size_t theta=0;theta<prepared->thetaCount;++theta)
          for(int n=-prepared->nMax;n<=prepared->nMax;++n)
            for(int m=-2;m<=2;++m)
              for(int degree=std::max(std::abs(m),std::abs(n));degree<=maxDegree;++degree)
                prepared->values[(((theta*static_cast<std::size_t>(2*prepared->nMax+1)+
                    static_cast<std::size_t>(n+prepared->nMax))*5+
                    static_cast<std::size_t>(m+2))*stride)+static_cast<std::size_t>(degree)]=
                    w[n,static_cast<int>(theta)][degree,m];
        const std::size_t bytes=prepared->values.size()*sizeof(double);
        if(samplingPlan&&samplingPlan->bytes+bytes<=SamplingPlanCache::byteBudget){
          samplingPlan->bytes+=bytes;
          angularCoefficients=samplingPlan->angularCoefficients.emplace(
              angularKey,std::move(prepared)).first->second;
        }else angularCoefficients=std::move(prepared);
      }
      std::vector<std::size_t> chunkPoints;
      for(const auto &group:chunks[chunk].radialGroups)
        chunkPoints.insert(chunkPoints.end(),group.begin(),group.end());
      std::vector<Sample> accumulated(chunkPoints.size());
      for(std::size_t i=0;i<chunkPoints.size();++i){
        const auto p=chunkPoints[i];pointForScratch[p]=i;
        accumulated[i].radialElement=elements[p];
        accumulated[i].fluid=sem.elementIsFluid(elements[p]);
      }
      std::map<std::pair<int,double>,RadialSamplingPlan> radialPlans;
      std::size_t radialPlanBytes=0;
      for(const auto &group:chunks[chunk].radialGroups){
        const auto representative=group.front();
        const auto radialKey=std::pair{elements[representative],radii[representative]};
        RadialSamplingPlan fallbackPlan;
        const RadialSamplingPlan *radialPlanPtr=nullptr;
        auto radialPlanIt=radialPlans.find(radialKey);
        if(samplingPlan){
          const auto sharedIt=samplingPlan->radialWeights.find(radialKey);
          if(sharedIt!=samplingPlan->radialWeights.end())radialPlanPtr=&sharedIt->second;
        }
        if(!radialPlanPtr){
          if(radialPlanIt==radialPlans.end()){
            fallbackPlan=makeRadialSamplingPlan(sem,radii[representative],elements[representative]);
            const std::size_t bytes=(fallbackPlan.shape.size()+fallbackPlan.derivative.size())*sizeof(double);
            if(samplingPlan&&samplingPlan->bytes+bytes<=SamplingPlanCache::byteBudget){
              samplingPlan->bytes+=bytes;
              radialPlanPtr=&samplingPlan->radialWeights.emplace(radialKey,std::move(fallbackPlan)).first->second;
            }else if(radialPlanBytes+bytes<=SamplingPlanCache::byteBudget){
              radialPlanBytes+=bytes;
              radialPlanIt=radialPlans.emplace(radialKey,std::move(fallbackPlan)).first;
            }
          }
          if(!radialPlanPtr&&radialPlanIt!=radialPlans.end())radialPlanPtr=&radialPlanIt->second;
        }
        if(!radialPlanPtr)radialPlanPtr=&fallbackPlan;
        const RadialSamplingPlan &radialPlan=*radialPlanPtr;
        const SamplingPlanCache::FrameKey coefficientFrameKey=frameKey;
        const auto coefficientKey=PoleCoefficientCache::Key{
            coefficientFrameKey,synthesisSet,elements[representative],
            std::bit_cast<std::uint64_t>(radii[representative])};
        PoleCoefficientRows fallbackCoefficients;
        const PoleCoefficientRows *coefficientsPtr=nullptr;
        if(frequencyCoefficientCache){
          const auto cached=frequencyCoefficientCache->rows.find(coefficientKey);
          if(cached!=frequencyCoefficientCache->rows.end()){coefficientsPtr=&cached->second;++auditprofile::counters["radial_coefficients_hit"];mathprofile::count("radial_coefficient_cache_hits");}
        }
        const auto buildCoefficients=[&](PoleCoefficientRows &coefficients){auditprofile::Scope t("radial_coefficients_build");mathprofile::Scope mt("radial_interpolation_and_strain_numerator_build");mathprofile::count("radial_coefficient_cache_misses");++auditprofile::counters["radial_coefficients_miss"];
          for(const auto &block:blocks){
            if(block.family==SensitivitySolveFamily::Toroidal&&
               (sem.elementIsFluid(elements[representative])||elements[representative]<sem.el()||elements[representative]>=sem.eu()))continue;
            const int mmax=(static_cast<int>(block.fields.cols())-1)/2;
            for(int m=-mmax;m<=mmax;++m){
              const auto values=interpolatePlanned(sem,block.family,block.fields,
                  block.lift,radialPlan,m+mmax,block.sourceElement);
              const auto modal=mode(values);
              auto &row=coefficients[m];
              if(row.values.size()<=static_cast<std::size_t>(block.degree)){
                row.values.resize(static_cast<std::size_t>(block.degree)+1);
                row.present.resize(static_cast<std::size_t>(block.degree)+1,0);
              }
              auto &combined=row.values[static_cast<std::size_t>(block.degree)];
              row.present[static_cast<std::size_t>(block.degree)]=1;
              if(block.family==SensitivitySolveFamily::Spheroidal){
                combined.U=modal.U;combined.Up=modal.Up;
                combined.V=modal.V;combined.Vp=modal.Vp;
              }else if(block.family==SensitivitySolveFamily::Toroidal){
                combined.W=modal.W;combined.Wp=modal.Wp;
              }else{
                combined.U=modal.U;combined.Up=modal.Up;
              }
            }
          }
          for(auto &[m,row]:coefficients){
            (void)m;
            row.strain.resize(row.values.size());
            for(std::size_t l=0;l<row.values.size();++l)if(row.present[l]){
              const int degree=static_cast<int>(l);
              row.strain[l]=poleStrainCoefficients(
                  degree,row.values[l],radii[representative],radiusScale);
              const char *bridgeEnv=std::getenv("DSFAST_RSTRAIN_BRIDGE");
              if(kernelRoute&&bridgeEnv&&std::string(bridgeEnv)=="1"&&
                 radii[representative]>0.0){
                // Exact nodal alternative: form r*strain from the SEM nodal
                // displacement/derivative values, interpolate that polynomial
                // with shape functions only, then divide by the output radius.
                std::array<Complex,6> numerator{};
                const int element=elements[representative];
                for(int q=0;q<sem.mesh().NN();++q){
                  const double nodeRadius=sem.mesh().NodeRadius(element,q);
                  const auto nodePlan=makeRadialSamplingPlan(sem,nodeRadius,element);
                  RadialMode nodalMode{};
                  for(const auto &block:blocks){
                    if(block.degree!=degree)continue;
                    if(block.family==SensitivitySolveFamily::Toroidal&&
                       (sem.elementIsFluid(element)||element<sem.el()||element>=sem.eu()))continue;
                    const int mmax=(static_cast<int>(block.fields.cols())-1)/2;
                    if(m < -mmax || m > mmax)continue;
                    const auto nodal=mode(interpolatePlanned(sem,block.family,block.fields,
                        block.lift,nodePlan,m+mmax,block.sourceElement));
                    if(block.family==SensitivitySolveFamily::Spheroidal){
                      nodalMode.U=nodal.U;nodalMode.Up=nodal.Up;
                      nodalMode.V=nodal.V;nodalMode.Vp=nodal.Vp;
                    }else if(block.family==SensitivitySolveFamily::Toroidal){
                      nodalMode.W=nodal.W;nodalMode.Wp=nodal.Wp;
                    }else{nodalMode.U=nodal.U;nodalMode.Up=nodal.Up;}
                  }
                  const auto rstrain=rScaledStrainFromModal(degree,nodalMode.U,nodalMode.Up,
                      nodalMode.V,nodalMode.Vp,nodalMode.W,nodalMode.Wp,nodeRadius,
                      1.0,1.0,1.0,1.0,1.0);
                  const std::array<Complex,6> values{rstrain.rr,rstrain.horizontalTrace,
                      rstrain.radialPlus,rstrain.radialMinus,rstrain.horizontalPlus,
                      rstrain.horizontalMinus};
                  for(std::size_t c=0;c<6;++c)
                    numerator[c]+=radialPlan.shape[static_cast<std::size_t>(q)]*values[c];
                }
                const double invRadius=1.0/radii[representative];
                for(std::size_t c=0;c<6;++c)row.strain[l].raw[c]=numerator[c]*invRadius;
              }
            }
          }
        };
        if(!coefficientsPtr){
          buildCoefficients(fallbackCoefficients);
          std::size_t bytes=0;
          for(const auto &[m,row]:fallbackCoefficients){
            (void)m;
            bytes+=sizeof(PoleDegreeRow)+row.values.capacity()*sizeof(RadialMode)+
                row.strain.capacity()*sizeof(PoleStrainCoefficients)+
                row.present.capacity()*sizeof(unsigned char);
          }
          if(frequencyCoefficientCache&&
             frequencyCoefficientCache->bytes+bytes<=PoleCoefficientCache::byteBudget){
            frequencyCoefficientCache->bytes+=bytes;
            coefficientsPtr=&frequencyCoefficientCache->rows.emplace(
                coefficientKey,std::move(fallbackCoefficients)).first->second;
          }else coefficientsPtr=&fallbackCoefficients;
        }
        const auto &coefficients=*coefficientsPtr;
        const bool adjoint=blocks.front().adjoint;
        if(coefficientVisitor){
          PoleAngularCoefficientRecord3D record;
          record.frame=frame;record.synthesisSet=synthesisSet;record.adjoint=adjoint;
          record.element=elements[representative];record.radiusModel=radii[representative];
          for(const auto &[m,row]:coefficients)for(std::size_t l=0;l<row.present.size();++l)
            if(row.present[l]){record.labels.emplace_back(m,static_cast<int>(l));
              record.rawStrain.push_back(row.strain[l].raw);}
          coefficientVisitor(record);
        }
        std::vector<Eigen::MatrixXcd> packedStrains;
        std::vector<std::pair<int,std::size_t>> packedKeys;
        if(kernelRoute&&contractionBackend>0){
          const bool fluidTrace=sem.elementIsFluid(elements[representative]);
          for(const auto &[m,row]:coefficients)
            for(std::size_t l=0;l<row.present.size();++l)
              if(row.present[l])packedKeys.emplace_back(m,l);
          packedStrains.resize(fluidTrace?1:6);
          const std::array<int,6> spins{0,0,1,-1,2,-2};
          for(std::size_t c=0;c<packedStrains.size();++c){
            const int spin=spins[c];
            std::vector<std::pair<int,std::size_t>> componentKeys;
            for(const auto &key:packedKeys)if(key.second>=static_cast<std::size_t>(std::abs(spin)))
              componentKeys.push_back(key);
            std::vector<std::uint64_t> basisKey;
            basisKey.reserve(3+group.size()*2+componentKeys.size()*2);
            basisKey.push_back(static_cast<std::uint64_t>(maxDegree));
            basisKey.push_back(static_cast<std::uint64_t>(spin+2));
            basisKey.push_back(static_cast<std::uint64_t>(group.size()));
            for(const auto p:group){
              basisKey.push_back(std::bit_cast<std::uint64_t>(angles[p].theta));
              basisKey.push_back(std::bit_cast<std::uint64_t>(angles[p].phi));
            }
            for(const auto &[m,l]:componentKeys){
              basisKey.push_back(static_cast<std::uint64_t>(static_cast<std::int64_t>(m)));
              basisKey.push_back(static_cast<std::uint64_t>(l));
            }
            std::shared_ptr<Eigen::MatrixXcd> cachedBasis;
            if(samplingPlan){
              const auto it=samplingPlan->packedBases.find(basisKey);
              if(it!=samplingPlan->packedBases.end()){
                cachedBasis=it->second;
                ++auditprofile::counters["packed_basis_hit"];
              }
            }
            if(!cachedBasis){
              ++auditprofile::counters["packed_basis_miss"];
              cachedBasis=std::make_shared<Eigen::MatrixXcd>(
                  static_cast<Eigen::Index>(group.size()),
                  static_cast<Eigen::Index>(componentKeys.size()));
              auto &basis=*cachedBasis;
              for(std::size_t q=0;q<componentKeys.size();++q){
                const auto [m,l]=componentKeys[q];
                for(std::size_t j=0;j<group.size();++j){
                  const auto p=group[j];
                  const std::size_t thetaIndex=angleForPoint[p]-begin;
                  const double angular=angularCoefficients->row(thetaIndex,spin,m)[l];
                  const Complex phase=std::abs(m)<=2?
                      angles[p].orderPhases[static_cast<std::size_t>(m+2)]:
                      std::exp(Complex(0.0,m*angles[p].phi));
                  basis(static_cast<Eigen::Index>(j),static_cast<Eigen::Index>(q))=angular*phase;
                }
              }
              const std::size_t bytes=static_cast<std::size_t>(cachedBasis->size())*sizeof(Complex);
              if(samplingPlan&&samplingPlan->bytes+bytes<=SamplingPlanCache::byteBudget){
                samplingPlan->bytes+=bytes;
                samplingPlan->packedBases.emplace(std::move(basisKey),cachedBasis);
              }
            }
            Eigen::MatrixXcd values(static_cast<Eigen::Index>(componentKeys.size()),1);
            for(std::size_t q=0;q<componentKeys.size();++q){
              const auto [m,l]=componentKeys[q];
              const auto &row=coefficients.at(m);
              values(static_cast<Eigen::Index>(q),0)=fluidTrace?
                  row.strain[l].raw[0]+row.strain[l].raw[1]:row.strain[l].raw[c];
            }
            packedStrains[c]=dsfast::synthesizeComplex(*cachedBasis,values,
                static_cast<dsfast::DotBackend>(contractionBackend));
          }
        }
        {auditprofile::Scope t("degree_sums_phases");
        for(const auto& [m,row]:coefficients)for(std::size_t l=0;l<row.present.size();++l)if(row.present[l])auditprofile::counters["strain_component_macs"]+=6*group.size();
        for(std::size_t j=0;j<group.size();++j){
          const auto p=group[j];
          const std::size_t thetaIndex=angleForPoint[p]-begin;
          auto &contribution=accumulated[pointForScratch[p]];
          const bool needDisplacement=kernelRoute?!adjoint&&p>=displacementBegin:true;
          bool packedStrainAdded=false;
          for(const auto &[m,row]:coefficients){
            std::array<Complex,6> strainSum{},rScaledSum{};
            std::array<Complex,3> displacementSum{};
            for(std::size_t l=0;l<row.values.size();++l){
              if(!row.present[l])continue;
              const auto &modal=row.values[l];
              const int degree=static_cast<int>(l);
              const double x0=angularCoefficients->row(thetaIndex,0,m)[l];
              const double xp=angularCoefficients->row(thetaIndex,1,m)[l];
              const double xm=angularCoefficients->row(thetaIndex,-1,m)[l];
              const double xpp=degree>=2?angularCoefficients->row(thetaIndex,2,m)[l]:0.0;
              const double xmm=degree>=2?angularCoefficients->row(thetaIndex,-2,m)[l]:0.0;
              const auto &strain=row.strain[l];
              const std::array<Complex,6> angular{x0,x0,xp,xm,xpp,xmm};
              for(std::size_t c=0;c<6;++c){
                if(!(kernelRoute&&contractionBackend>0))strainSum[c]+=strain.raw[c]*angular[c];
                if(!kernelRoute)rScaledSum[c]+=strain.rScaled[c]*angular[c];
              }
              if(needDisplacement){
                const auto spheroidal=radialToDisplacement(degree,modal,x0,xp,xm);
                const double zeta=std::sqrt(degree*(degree+1.0));
                displacementSum[0]+=spheroidal[0];
                displacementSum[1]+=spheroidal[1]+Complex(0,1)*zeta/2.0*modal.W*(xp+xm);
                displacementSum[2]+=-Complex(0,1)*zeta/2.0*modal.V*(xp+xm)+
                    zeta/2.0*modal.W*(xp-xm);
              }
            }
            const Complex phase=std::abs(m)<=2?
                angles[p].orderPhases[static_cast<std::size_t>(m+2)]:
                std::exp(Complex(0.0,m*angles[p].phi));
            const auto displacement=needDisplacement?std::array<Complex,3>{
                phase*displacementSum[0],phase*displacementSum[1],phase*displacementSum[2]}:
                std::array<Complex,3>{};
            SpinStrainField strain;
            if(kernelRoute&&contractionBackend>0){
              if(!packedStrainAdded){
                if(packedStrains.size()==1){
                  strainSum[0]=packedStrains[0](static_cast<Eigen::Index>(j),0);
                  strainSum[1]=0.0;
                }else for(std::size_t c=0;c<6;++c)
                  strainSum[c]=packedStrains[c](static_cast<Eigen::Index>(j),0);
                strain={strainSum[0],strainSum[1],strainSum[2],strainSum[3],strainSum[4],strainSum[5]};
                packedStrainAdded=true;
              }
            }else strain={phase*strainSum[0],phase*strainSum[1],phase*strainSum[2],
                phase*strainSum[3],phase*strainSum[4],phase*strainSum[5]};
            SpinStrainField rScaled{phase*rScaledSum[0],phase*rScaledSum[1],
                phase*rScaledSum[2],phase*rScaledSum[3],phase*rScaledSum[4],phase*rScaledSum[5]};
            if(adjoint){
              add(contribution.adjointDualDisplacement,displacement);
              add(contribution.adjointDualStrain,strain);
              add(contribution.rScaledAdjointDualStrainMetres,rScaled);
            }else{
              add(contribution.forwardRawDisplacement,displacement);
              add(contribution.forwardRawStrain,strain);
              add(contribution.rScaledForwardStrainMetres,rScaled);
            }
          }
        }
      }
      }
      auditprofile::Scope rotationTimer("rotate_add_scratch");
      for(std::size_t i=0;i<chunkPoints.size();++i){
        const auto p=chunkPoints[i];
        const bool adjoint=blocks.front().adjoint;
        if(adjoint){
          rotateContribution(accumulated[i],angles[p].a,angles[p].b,true);
          add(result[p].adjointDualDisplacement,accumulated[i].adjointDualDisplacement);
          add(result[p].adjointDualStrain,accumulated[i].adjointDualStrain);
          add(result[p].rScaledAdjointDualStrainMetres,accumulated[i].rScaledAdjointDualStrainMetres);
        }else{
          rotateContribution(accumulated[i],angles[p].a,angles[p].b,false);
          add(result[p].forwardRawDisplacement,accumulated[i].forwardRawDisplacement);
          add(result[p].forwardRawStrain,accumulated[i].forwardRawStrain);
          add(result[p].rScaledForwardStrainMetres,accumulated[i].rScaledForwardStrainMetres);
        }
        pointForScratch[p]=std::numeric_limits<std::size_t>::max();
      }
    }
  }

  static void add(std::array<Complex,3> &a, const std::array<Complex,3> &b) {
    for (int i=0;i<3;++i) a[i] += b[i];
  }
  static void add(SpinStrainField &a, const SpinStrainField &b) {
    a.rr += b.rr; a.horizontalTrace += b.horizontalTrace;
    a.radialPlus += b.radialPlus; a.radialMinus += b.radialMinus;
    a.horizontalPlus += b.horizontalPlus; a.horizontalMinus += b.horizontalMinus;
  }

  static SpinStrainField scale(const SpinStrainField &a, double x) {
    return {x*a.rr, x*a.horizontalTrace, x*a.radialPlus, x*a.radialMinus,
            x*a.horizontalPlus, x*a.horizontalMinus};
  }
};

/** Owned modal solutions for one frequency; callers may sample arbitrarily
 * sized chunks without repeating a factorization or retaining volume fields. */
class Windowed3DFrequencyWorkspace {
  friend class Windowed3DFieldEngine;
public:
  using Complex=std::complex<double>;
  /** Synchronous view of one post-projection radial mode at a native SEM
   * node. `values` is [U,U',V,V',W,W'] in normalized-radius units; it is valid
   * only during the callback. The visitor exposes already weighted adjoint
   * duals and applies the source lift only on the recorded source element. */
  using ProjectedRadialModeVisitor = std::function<void(
      SensitivitySolveFamily,int,bool,int,int,int,int,double,bool,bool,
      const std::array<std::array<double,3>,3> &,
      const Eigen::MatrixXcd &)>;
  struct AdjointSynthesis {
    windowed3d_detail::PoleFrame3D frame;
    std::vector<windowed3d_detail::PoleSynthesisBlock3D> blocks;
    std::uint64_t synthesisSet=0;
  };
  static Windowed3DFrequencyWorkspace prepare(Full1D::SEM &sem,
      SourceInfo::EarthquakeCMT &source,InputParameters &parameters,Complex z,
      const Eigen::VectorXcd &adjointWeights,
      dssolver::BatchResponseCache *batchCache=nullptr,bool poleOnly=false,
      bool forwardOnly=false,double *basisReconstructionRelative=nullptr,int solveThreads=1,std::optional<Complex> attenuation=std::nullopt) {
    mathprofile::Scope profileScope("frequency_workspace_prepare_and_modal_capture");
    sensitivity_detail::requireConstructedSEM(sem);
    auditprofile::Scope t("prepare_workspace");
    std::vector<OwnedModalFrequencyBlock3D> blocks;
    std::vector<windowed3d_detail::PoleModalFrequencyBlock3D> poleBlocks;
    windowed3d_detail::PoleFrame3D sourceFrame;
    const auto capture=[&](const ModalFrequencyBlock3D &b){
      if((poleOnly||forwardOnly)&&b.family!=SensitivitySolveFamily::Radial)return;
      mathprofile::Scope captureCopy("modal_response_copy_and_capture");
#ifndef MATHPROFILE_DISABLED
      mathprofile::count("modal_response_capture_copy_bytes",static_cast<std::uint64_t>(
          (b.forward.size()+b.adjoint.size()+b.sourceLift.size()+
           b.sourceCoefficients.size()+b.modeReceiver.size())*sizeof(Complex)));
#endif
      blocks.push_back({b.family,b.degree,b.sourceElement,b.forward,b.adjoint,
                        b.sourceLift,b.sourceCoefficients,b.modeReceiver});
    };
    Windowed3DFieldEngine::solveFrequencyModes(
        sem,source,parameters,z,adjointWeights,capture,&poleBlocks,&sourceFrame,batchCache,poleOnly,forwardOnly,
        basisReconstructionRelative,solveThreads,attenuation);
    return Windowed3DFrequencyWorkspace(sem,parameters,z,adjointWeights,
        std::move(blocks),std::move(poleBlocks),sourceFrame,forwardOnly);
  }
  std::vector<Field3DSample> sample(const std::vector<Field3DPoint> &points) const {
    return samplePoleFrames(points,false);
  }
  std::vector<Field3DSample> sampleWeightedChunk(
      const std::vector<Field3DPoint> &points) const {
    return samplePoleFrames(points,true);
  }
  /** Read-only Phase A diagnostic hook. It streams every native element/node
   * side, including the separate radial l=0 block, without retaining a modal
   * history or changing field sampling, cache, or solver state. Arguments are
   * (family, degree, adjointRole, receiver, m, element, localNode, r_nd, fluid,
   * sourceLiftApplied, frameAxes, radialValues). */
  void visitProjectedRadialModes(const ProjectedRadialModeVisitor &visitor,
                                double minimumRadiusMetres=0.0) const {
    if (!visitor) throw std::invalid_argument("Projected radial mode visitor is empty.");
    if (!(std::isfinite(minimumRadiusMetres) && minimumRadiusMetres>=0.0))
      throw std::invalid_argument("Projected radial mode visitor radius bound is invalid.");
    const auto &mesh=sem_->mesh();
    Eigen::MatrixXcd noLift;
    for(int element=0;element<mesh.NE();++element){
      const bool fluid=sem_->elementIsFluid(element);
      for(int q=0;q<mesh.NN();++q){
        const double radius=mesh.NodeRadius(element,q);
        if (radius*sem_->lengthScaleMetres()<minimumRadiusMetres) continue;
        const auto plan=Windowed3DFieldEngine::makeRadialSamplingPlan(*sem_,radius,element);
        const auto emit=[&](SensitivitySolveFamily family,int degree,bool adjoint,
                            int receiver,int m,const Eigen::MatrixXcd &fields,
                            const Eigen::MatrixXcd &lift,int liftElement){
          if(family==SensitivitySolveFamily::Toroidal &&
             (fluid || element<sem_->el() || element>=sem_->eu())) return;
          const int mmax=(static_cast<int>(fields.cols())-1)/2;
          if(m < -mmax || m > mmax) return;
          const bool liftApplied=!adjoint && element==liftElement &&
              lift.rows()==fields.rows() && lift.cols()==fields.cols();
          const auto values=Windowed3DFieldEngine::interpolatePlanned(
              *sem_,family,fields,lift,plan,m+mmax,liftElement);
          const auto &axes=adjoint&&receiver>=0
              ?adjointSynthesis_.at(static_cast<std::size_t>(receiver)).frame.axes
              :sourceFrame_.axes;
          visitor(family,degree,adjoint,receiver,m,element,q,radius,fluid,liftApplied,axes,values);
        };
        // l=0 has no pole-frame synthesis vector. Export its forward and
        // weighted-adjoint radial coefficients directly from captured blocks.
        for(const auto &block:blocks_)if(block.family==SensitivitySolveFamily::Radial){
          emit(block.family,block.degree,false,-1,0,block.forward,block.sourceLift,block.sourceElement);
          if(block.adjoint.cols()!=0)
            emit(block.family,block.degree,true,-1,0,block.adjoint,noLift,-1);
        }
        for(const auto &block:sourceSynthesis_){
          const int mmax=(static_cast<int>(block.fields.cols())-1)/2;
          for(int m=-mmax;m<=mmax;++m)
            emit(block.family,block.degree,false,-1,m,block.fields,block.lift,block.sourceElement);
        }
        for(const auto &receiver:adjointSynthesis_)
          for(const auto &block:receiver.blocks){
            const int mmax=(static_cast<int>(block.fields.cols())-1)/2;
            for(int m=-mmax;m<=mmax;++m)
              emit(block.family,block.degree,true,static_cast<int>(receiver.synthesisSet)-1,
                   m,block.fields,noLift,-1);
          }
      }
    }
  }
  /** Return local (theta,phi) and accepted frame-rotation (a,b) values for
   * the exact field frame used by the modal visitor. */
  std::array<double,4> phaseAngularCoordinates(bool adjoint,int receiver,
                                               double theta,double phi) const {
    const windowed3d_detail::PoleFrame3D *frame=&sourceFrame_;
    if(adjoint){
      if(receiver<0||static_cast<std::size_t>(receiver)>=adjointSynthesis_.size())
        throw std::out_of_range("Phase A adjoint receiver frame index");
      frame=&adjointSynthesis_[static_cast<std::size_t>(receiver)].frame;
    }
    double a=0.0,b=0.0;
    const auto local=Windowed3DFieldEngine::tangentRotation(*frame,theta,phi,a,b);
    return {local[0],local[1],a,b};
  }
  /** Kernel route materializes strain everywhere and forward displacement only
   * at the appended receiver suffix. Adjoint displacement is not needed by
   * the local material contraction. */
  std::vector<KernelField3DSample> sampleForKernels(
      const std::vector<Field3DPoint> &points,std::size_t receiverBegin,
      Windowed3DFieldEngine::SamplingPlanCache *samplingPlan=nullptr) const {
    auditprofile::Scope t("sample_compact_total");
    if(receiverBegin>points.size())
      throw std::invalid_argument("Kernel receiver suffix starts outside the point list.");
    const auto full=samplePoleFrames(points,false,true,receiverBegin,samplingPlan);
    auditprofile::Scope copyTimer("compact_copy");
    std::vector<KernelField3DSample> compact(full.size());
    for(std::size_t p=0;p<full.size();++p){
      compact[p].radialElement=full[p].radialElement;compact[p].fluid=full[p].fluid;
      compact[p].forwardRawDisplacement=full[p].forwardRawDisplacement;
      compact[p].forwardRawStrain=full[p].forwardRawStrain;
      compact[p].adjointDualStrain=full[p].adjointDualStrain;
    }
    return compact;
  }
  std::vector<KernelField3DSample> sampleForKernelsPacked(
      const std::vector<Field3DPoint> &points,std::size_t receiverBegin,
      Windowed3DFieldEngine::SamplingPlanCache *samplingPlan,int contractionBackend,
      const std::function<void(const PoleAngularCoefficientRecord3D&)> &coefficientVisitor={}) const {
    if(receiverBegin>points.size()||contractionBackend<1||contractionBackend>3)
      throw std::invalid_argument("Packed kernel sample arguments are invalid.");
    auditprofile::Scope t("sample_packed_total");
    const auto full=samplePoleFrames(points,false,true,receiverBegin,samplingPlan,contractionBackend,coefficientVisitor);
    mathprofile::Scope compactCopy("packed_to_compact_kernel_sample_copy");
#ifndef MATHPROFILE_DISABLED
    mathprofile::count("packed_to_compact_kernel_sample_copy_bytes",
        static_cast<std::uint64_t>(full.size()*sizeof(KernelField3DSample)));
#endif
    std::vector<KernelField3DSample> compact(full.size());
    for(std::size_t p=0;p<full.size();++p){
      compact[p].radialElement=full[p].radialElement;compact[p].fluid=full[p].fluid;
      compact[p].forwardRawDisplacement=full[p].forwardRawDisplacement;
      compact[p].forwardRawStrain=full[p].forwardRawStrain;
      compact[p].adjointDualStrain=full[p].adjointDualStrain;
    }
    compactCopy.stop();
    return compact;
  }

  /** Tensor contraction prototype: gather solver-derived radial coefficients
   * once per (element,radius), then contract angular bases against bounded
   * tiles of radial columns using genuine matrix-matrix products. The sampled
   * volume must be a complete Cartesian product of one shared angle set and
   * one-sided radial traces; use sampleForKernelsPacked for arbitrary clouds. */
  std::vector<KernelField3DSample> sampleForKernelsTensor(
      const std::vector<Field3DPoint> &points,std::size_t receiverBegin,
      Windowed3DFieldEngine::SamplingPlanCache *samplingPlan,int contractionBackend,
      const std::function<void(const PoleAngularCoefficientRecord3D&)> &coefficientVisitor={},
      bool realDegreeLayout=false,
      const CachedForwardTensorFields3D *cachedForward=nullptr) const {
    mathprofile::Scope all("tensor_sampling_total");
    mathprofile::count("tensor_volume_point_count",receiverBegin);
    mathprofile::count("tensor_receiver_suffix_count",points.size()-receiverBegin);
    if(receiverBegin>points.size()||contractionBackend<1||contractionBackend>3)
      throw std::invalid_argument("Tensor kernel sample arguments are invalid.");
    if(receiverBegin==0)throw std::invalid_argument("Tensor sample requires a receiver suffix.");
    if(cachedForward&&(cachedForward->volumeStrain.size()!=receiverBegin||
        cachedForward->receiverDisplacement.size()!=points.size()-receiverBegin))
      throw std::invalid_argument("Cached forward tensor fields do not match the volume/receiver grid.");
    using TraceKey=Windowed3DFieldEngine::SamplingPlanCache::TensorTraceKey;
    using AngleKey=Windowed3DFieldEngine::SamplingPlanCache::PointKey;
    const auto traceKey=[&](const Field3DPoint &p){return TraceKey{p.radialElement,
        std::bit_cast<std::uint64_t>(Windowed3DFieldEngine::canonicalSamplingRadius(
            *sem_,p.radialElement,p.radiusMetres/sem_->lengthScaleMetres()))};};
    const auto angleKey=[](const Field3DPoint &p){return AngleKey{
        std::bit_cast<std::uint64_t>(p.thetaRadians),std::bit_cast<std::uint64_t>(p.phiRadians)};};
    using TensorGridLayout=Windowed3DFieldEngine::SamplingPlanCache::TensorGridLayout;
    std::vector<std::uint64_t> gridKey;
    TensorGridLayout localGrid;
    std::shared_ptr<TensorGridLayout> gridOwner;
    const TensorGridLayout *gridPtr=nullptr;
    {
      mathprofile::Scope gridKeyLookup("tensor_grid_semantic_key_build_and_lookup");
      gridKey.reserve(2+4*receiverBegin);gridKey.push_back(0x4453464153544752ULL);
      gridKey.push_back(static_cast<std::uint64_t>(receiverBegin));
      for(std::size_t p=0;p<receiverBegin;++p){
        const auto tk=traceKey(points[p]);const auto ak=angleKey(points[p]);
        gridKey.push_back(static_cast<std::uint64_t>(static_cast<std::int64_t>(tk.first)));
        gridKey.push_back(tk.second);gridKey.push_back(ak.first);gridKey.push_back(ak.second);
      }
      if(samplingPlan){
        const auto cached=samplingPlan->tensorGridLayouts.find(gridKey);
        if(cached!=samplingPlan->tensorGridLayouts.end()){
          gridOwner=cached->second;gridPtr=gridOwner.get();++samplingPlan->tensorGridHits;mathprofile::count("tensor_grid_cache_hits");
        }
      }
    }
    if(!gridPtr){
      mathprofile::Scope gridCold("tensor_grid_cold_layout_map_build");
      if(samplingPlan)++samplingPlan->tensorGridMisses;mathprofile::count("tensor_grid_cache_misses");
      auto &grid=localGrid;
      std::map<TraceKey,std::size_t> traceIndex;
      std::map<AngleKey,std::size_t> angleIndex;
      for(std::size_t p=0;p<receiverBegin;++p){
        const auto tk=traceKey(points[p]);
        auto [ti,added]=traceIndex.try_emplace(tk,grid.traces.size());
        if(added){grid.traces.push_back(tk);grid.pointForTrace.emplace_back();}
        const auto ak=angleKey(points[p]);
        auto [ai,newAngle]=angleIndex.try_emplace(ak,grid.angleKeys.size());
        if(newAngle)grid.angleKeys.push_back(ak);
        auto [pi,pointAdded]=grid.pointForTrace[ti->second].try_emplace(ak,p);
        (void)pi;
        if(!pointAdded)throw std::invalid_argument("Tensor grid has duplicate angle samples on a radial trace.");
      }
      if(grid.traces.empty()||grid.angleKeys.empty())throw std::invalid_argument("Tensor grid is empty.");
      for(const auto &trace:grid.pointForTrace)if(trace.size()!=grid.angleKeys.size())
        throw std::invalid_argument("Tensor grid is not a shared angle-set by radial-trace product.");
      std::size_t estimate=sizeof(TensorGridLayout)+gridKey.capacity()*sizeof(std::uint64_t)+
          grid.traces.capacity()*sizeof(TraceKey)+grid.angleKeys.capacity()*sizeof(AngleKey)+
          grid.pointForTrace.capacity()*sizeof(std::map<AngleKey,std::size_t>);
      for(const auto &trace:grid.pointForTrace)
        estimate+=trace.size()*(sizeof(std::pair<const AngleKey,std::size_t>)+4*sizeof(void*));
      if(samplingPlan&&samplingPlan->bytes+estimate<=
          Windowed3DFieldEngine::SamplingPlanCache::byteBudget){
        grid.payloadBytes=estimate;samplingPlan->bytes+=estimate;
        samplingPlan->tensorGridPayloadBytes+=estimate;
        gridOwner=std::make_shared<TensorGridLayout>(std::move(localGrid));
        gridPtr=gridOwner.get();
        samplingPlan->tensorGridLayouts.emplace(std::move(gridKey),gridOwner);
      }else gridPtr=&localGrid;
    }
    const auto &traces=gridPtr->traces;
    const auto &angleKeys=gridPtr->angleKeys;
    const auto &pointForTrace=gridPtr->pointForTrace;

    std::vector<Field3DPoint> representatives;
    std::map<std::tuple<std::uint64_t,bool,int,std::uint64_t>,PoleAngularCoefficientRecord3D> records;
    {
      mathprofile::Scope capture("tensor_representative_and_coefficient_record_capture");
      representatives.reserve(traces.size());
      for(std::size_t t=0;t<traces.size();++t)
        representatives.push_back(points[pointForTrace[t].begin()->second]);
      const auto visitor=[&](const PoleAngularCoefficientRecord3D &record){
        const auto key=std::tuple{record.synthesisSet,record.adjoint,record.element,
            std::bit_cast<std::uint64_t>(record.radiusModel)};
        records.try_emplace(key,record);
        if(coefficientVisitor)coefficientVisitor(record);
      };
      (void)samplePoleFrames(representatives,false,true,representatives.size(),
          samplingPlan,contractionBackend,visitor);
    }

    auto result=[&]{
      mathprofile::Scope allocation("tensor_result_allocation_and_volume_tagging");
      std::vector<KernelField3DSample> samples(points.size());
      for(std::size_t p=0;p<receiverBegin;++p){
        samples[p].radialElement=points[p].radialElement;
        samples[p].fluid=sem_->elementIsFluid(points[p].radialElement);
      }
      return samples;
    }();
    const int maxDegree=parameters_->lmax();
    constexpr std::size_t radialTile=16;
    const auto synthesisOne=[&](std::uint64_t set,bool adjoint){
      std::vector<const PoleAngularCoefficientRecord3D*> byTrace(traces.size(),nullptr);
      for(std::size_t t=0;t<traces.size();++t){
        const auto it=records.find({set,adjoint,traces[t].first,traces[t].second});
        if(it!=records.end())byTrace[t]=&it->second;
      }
      if(std::none_of(byTrace.begin(),byTrace.end(),[](auto p){return p!=nullptr;}))return;
      const auto *frameRecord=*std::find_if(byTrace.begin(),byTrace.end(),[](auto p){return p!=nullptr;});
      Windowed3DFieldEngine::PoleFrame frame{frameRecord->frame.axes};
      using TensorAngularPlan=Windowed3DFieldEngine::SamplingPlanCache::TensorAngularPlan;
      std::vector<std::uint64_t> tensorPlanKey;
      TensorAngularPlan localPlan;
      std::shared_ptr<TensorAngularPlan> planOwner;
      TensorAngularPlan *planPtr=nullptr;
      {
        mathprofile::Scope angularKeyLookup("tensor_angular_plan_key_build_and_lookup");
        tensorPlanKey.reserve(11+2*angleKeys.size());
        tensorPlanKey.push_back(0x445346415354414eULL);
        tensorPlanKey.push_back(static_cast<std::uint64_t>(maxDegree));
        tensorPlanKey.push_back(static_cast<std::uint64_t>(angleKeys.size()));
        for(const auto &axis:frame.axes)for(const double x:axis)
          tensorPlanKey.push_back(std::bit_cast<std::uint64_t>(x));
        for(const auto &[thetaBits,phiBits]:angleKeys){tensorPlanKey.push_back(thetaBits);tensorPlanKey.push_back(phiBits);}
        if(samplingPlan){
          const auto cached=samplingPlan->tensorAngularPlans.find(tensorPlanKey);
          if(cached!=samplingPlan->tensorAngularPlans.end()){
            planOwner=cached->second;planPtr=planOwner.get();++samplingPlan->tensorPlanHits;mathprofile::count("tensor_angular_plan_cache_hits");
          }
        }
      }
      if(!planPtr){
        mathprofile::Scope angularCold("tensor_cold_frame_geometry_and_wigner_build");
        if(samplingPlan)++samplingPlan->tensorPlanMisses;mathprofile::count("tensor_angular_plan_cache_misses");
        localPlan.localTheta.resize(angleKeys.size());localPlan.localPhi.resize(angleKeys.size());
        localPlan.rotA.resize(angleKeys.size());localPlan.rotB.resize(angleKeys.size());
        for(std::size_t a=0;a<angleKeys.size();++a){
          const auto &point=points[pointForTrace.front().at(angleKeys[a])];
          double aa=1.0,bb=0.0;
          const auto local=Windowed3DFieldEngine::tangentRotation(frame,
              point.thetaRadians,point.phiRadians,aa,bb);
          localPlan.localTheta[a]=local[0];localPlan.localPhi[a]=local[1];
          localPlan.rotA[a]=aa;localPlan.rotB[a]=bb;
        }
        localPlan.wigner=std::make_shared<Windowed3DFieldEngine::CachedAngularWigner>(
            maxDegree,std::min(2,maxDegree),std::min(2,maxDegree),localPlan.localTheta);
        const std::size_t estimate=tensorPlanKey.capacity()*sizeof(std::uint64_t)+
            4*angleKeys.size()*sizeof(double)+
            angleKeys.size()*Windowed3DFieldEngine::angularWignerBytesPerAngle(
                maxDegree,std::min(2,maxDegree));
        if(samplingPlan&&samplingPlan->bytes+estimate<=
            Windowed3DFieldEngine::SamplingPlanCache::byteBudget){
          localPlan.payloadBytes=estimate;samplingPlan->bytes+=estimate;
          samplingPlan->tensorAngularPayloadBytes+=estimate;
          planOwner=std::make_shared<TensorAngularPlan>(std::move(localPlan));
          planPtr=planOwner.get();
          samplingPlan->tensorAngularPlans.emplace(std::move(tensorPlanKey),planOwner);
        }else planPtr=&localPlan;
      }
      const auto &localTheta=planPtr->localTheta;
      const auto &localPhi=planPtr->localPhi;
      const auto &rotA=planPtr->rotA;
      const auto &rotB=planPtr->rotB;
      const auto &wigner=*planPtr->wigner;
      const auto complexBasisFor=[&](int spin,const std::vector<std::pair<int,int>> &keys){
        auto basisKey=[&]{
          mathprofile::Scope keyBuild("tensor_complex_basis_key_construction");
          std::vector<std::uint64_t> key{0x4453464153544342ULL,
              static_cast<std::uint64_t>(spin),static_cast<std::uint64_t>(keys.size())};
          for(const auto &[m,l]:keys){key.push_back(static_cast<std::uint64_t>(static_cast<std::int64_t>(m)));key.push_back(static_cast<std::uint64_t>(l));}
          return key;
        }();
        {
          mathprofile::Scope basisKeyLookup("tensor_complex_basis_cache_lookup");
          const auto found=planPtr->complexBases.find(basisKey);
          if(found!=planPtr->complexBases.end()){
            if(samplingPlan)++samplingPlan->tensorBasisHits;
            mathprofile::count("tensor_basis_cache_hits");
            return found->second;
          }
        }
        mathprofile::Scope basisFill("tensor_complex_basis_cold_fill_and_cache");
        if(samplingPlan)++samplingPlan->tensorBasisMisses;
        mathprofile::count("tensor_basis_cache_misses");
        auto basis=std::make_shared<Eigen::MatrixXcd>(
            static_cast<Eigen::Index>(angleKeys.size()),static_cast<Eigen::Index>(keys.size()));
        for(std::size_t a=0;a<angleKeys.size();++a)for(std::size_t q=0;q<keys.size();++q){
          const auto [m,l]=keys[q];
          (*basis)(static_cast<Eigen::Index>(a),static_cast<Eigen::Index>(q))=
              wigner[spin,static_cast<int>(a)][l,m]*std::exp(Complex(0.0,m*localPhi[a]));
        }
        const std::size_t estimate=basisKey.capacity()*sizeof(std::uint64_t)+
            static_cast<std::size_t>(basis->size())*sizeof(Complex);
        if(planOwner&&samplingPlan&&samplingPlan->bytes+estimate<=
            Windowed3DFieldEngine::SamplingPlanCache::byteBudget){
          samplingPlan->bytes+=estimate;planOwner->payloadBytes+=estimate;
          samplingPlan->tensorAngularPayloadBytes+=estimate;
          planOwner->complexBases.emplace(std::move(basisKey),basis);
        }
        return basis;
      };
      const auto realBasisFor=[&](int spin,int m,const std::vector<int> &degrees){
        auto basisKey=[&]{
          mathprofile::Scope keyBuild("tensor_real_basis_key_construction");
          std::vector<std::uint64_t> key{0x4453464153545242ULL,
              static_cast<std::uint64_t>(spin),static_cast<std::uint64_t>(static_cast<std::int64_t>(m)),
              static_cast<std::uint64_t>(degrees.size())};
          for(const int degree:degrees)key.push_back(static_cast<std::uint64_t>(degree));
          return key;
        }();
        {
          mathprofile::Scope basisKeyLookup("tensor_real_basis_cache_lookup");
          const auto found=planPtr->realBases.find(basisKey);
          if(found!=planPtr->realBases.end()){
            if(samplingPlan)++samplingPlan->tensorBasisHits;
            mathprofile::count("tensor_basis_cache_hits");
            return found->second;
          }
        }
        mathprofile::Scope basisFill("tensor_real_basis_cold_fill_and_cache");
        if(samplingPlan)++samplingPlan->tensorBasisMisses;
        mathprofile::count("tensor_basis_cache_misses");
        auto basis=std::make_shared<Eigen::MatrixXd>(
            static_cast<Eigen::Index>(angleKeys.size()),static_cast<Eigen::Index>(degrees.size()));
        for(std::size_t a=0;a<angleKeys.size();++a)for(std::size_t q=0;q<degrees.size();++q)
          (*basis)(static_cast<Eigen::Index>(a),static_cast<Eigen::Index>(q))=
              wigner[spin,static_cast<int>(a)][degrees[q],m];
        const std::size_t estimate=basisKey.capacity()*sizeof(std::uint64_t)+
            static_cast<std::size_t>(basis->size())*sizeof(double);
        if(planOwner&&samplingPlan&&samplingPlan->bytes+estimate<=
            Windowed3DFieldEngine::SamplingPlanCache::byteBudget){
          samplingPlan->bytes+=estimate;planOwner->payloadBytes+=estimate;
          samplingPlan->tensorAngularPayloadBytes+=estimate;
          planOwner->realBases.emplace(std::move(basisKey),basis);
        }
        return basis;
      };
      std::vector<std::array<Complex,6>> localFields(angleKeys.size()*traces.size());
      const std::array<int,6> spins{0,0,1,-1,2,-2};
      using CoefficientLabel=std::pair<int,int>;
      constexpr std::size_t missing=std::numeric_limits<std::size_t>::max();
      std::vector<std::vector<CoefficientLabel>> uniqueLabelLayouts;
      std::vector<std::size_t> labelLayoutForTrace(traces.size(),missing);
      std::array<std::vector<CoefficientLabel>,3> keysByAbsSpin;
      std::array<std::vector<std::vector<std::size_t>>,3> permutationsByAbsSpin;
      {
      mathprofile::Scope layoutBuild("tensor_coefficient_layout_dedup_and_permutation_build");
      std::unordered_map<std::uint64_t,std::vector<std::size_t>> layoutHashBuckets;
      const auto labelHash=[](const std::vector<CoefficientLabel> &labels){
        std::uint64_t hash=1469598103934665603ULL;
        for(const auto &[m,l]:labels){
          hash^=static_cast<std::uint64_t>(static_cast<std::int64_t>(m));hash*=1099511628211ULL;
          hash^=static_cast<std::uint64_t>(static_cast<std::int64_t>(l));hash*=1099511628211ULL;
        }
        hash^=static_cast<std::uint64_t>(labels.size());hash*=1099511628211ULL;
        return hash;
      };
      for(std::size_t t=0;t<traces.size();++t){
        const auto *record=byTrace[t];if(!record)continue;
        const auto hash=labelHash(record->labels);
        auto &bucket=layoutHashBuckets[hash];
        auto found=std::find_if(bucket.begin(),bucket.end(),[&](std::size_t id){
          return uniqueLabelLayouts[id]==record->labels;
        });
        if(found==bucket.end()){
          labelLayoutForTrace[t]=uniqueLabelLayouts.size();
          uniqueLabelLayouts.push_back(record->labels);bucket.push_back(labelLayoutForTrace[t]);
        }else labelLayoutForTrace[t]=*found;
      }
      mathprofile::count("tensor_packing_unique_label_layouts",
          static_cast<std::uint64_t>(uniqueLabelLayouts.size()));
      for(std::size_t absSpin=0;absSpin<keysByAbsSpin.size();++absSpin){
        auto &keys=keysByAbsSpin[absSpin];
        for(const auto &layout:uniqueLabelLayouts)
          for(const auto &[m,l]:layout)if(l>=static_cast<int>(absSpin))keys.emplace_back(m,l);
        mathprofile::count("tensor_packing_union_label_inputs",
            static_cast<std::uint64_t>(keys.size()));
        std::sort(keys.begin(),keys.end());keys.erase(std::unique(keys.begin(),keys.end()),keys.end());
        mathprofile::count("tensor_packing_unique_keys",static_cast<std::uint64_t>(keys.size()));
        auto &permutations=permutationsByAbsSpin[absSpin];
        permutations.resize(uniqueLabelLayouts.size());
        std::size_t entries=0;
        for(std::size_t layoutId=0;layoutId<uniqueLabelLayouts.size();++layoutId){
          const auto &layout=uniqueLabelLayouts[layoutId];
          auto &permutation=permutations[layoutId];
          permutation.assign(layout.size(),missing);
          std::size_t key=0;
          for(std::size_t q=0;q<layout.size();++q){
            const auto &label=layout[q];if(label.second<static_cast<int>(absSpin))continue;
            while(key<keys.size()&&keys[key]<label)++key;
            if(key<keys.size()&&keys[key]==label)permutation[q]=key;
          }
          entries+=layout.size();
        }
        mathprofile::count("tensor_packing_record_permutation_entries",
            static_cast<std::uint64_t>(entries));
        mathprofile::count("tensor_packing_tree_record_lookups_avoided",
            static_cast<std::uint64_t>(entries));
      }
      }
      for(std::size_t component=0;component<6;++component){
        mathprofile::Scope packing("tensor_coefficient_packing_and_key_reorder");
        const std::size_t absSpin=static_cast<std::size_t>(std::abs(spins[component]));
        const auto &keys=keysByAbsSpin[absSpin];
        const auto &recordKeyIndex=permutationsByAbsSpin[absSpin];
        if(keys.empty())continue;
        if(!realDegreeLayout){
          const auto basis=complexBasisFor(spins[component],keys);
          for(std::size_t first=0;first<traces.size();first+=radialTile){
            const std::size_t cols=std::min(radialTile,traces.size()-first);
            Eigen::MatrixXcd coefficients=Eigen::MatrixXcd::Zero(
                static_cast<Eigen::Index>(keys.size()),static_cast<Eigen::Index>(cols));
            for(std::size_t r=0;r<cols;++r){
              const auto *record=byTrace[first+r];if(!record)continue;
              const auto &permutation=recordKeyIndex[labelLayoutForTrace[first+r]];
              const bool fluid=sem_->elementIsFluid(traces[first+r].first);
              for(std::size_t q=0;q<permutation.size();++q){
                const std::size_t key=permutation[q];if(key==missing)continue;
                const auto &values=record->rawStrain[q];
                coefficients(static_cast<Eigen::Index>(key),static_cast<Eigen::Index>(r))=
                    fluid?(component==0?values[0]+values[1]:component==1?Complex{}:Complex{}):
                          values[component];
              }
            }
            const auto synthesized=dsfast::synthesizeComplex(*basis,coefficients,
                static_cast<dsfast::DotBackend>(contractionBackend));
            for(std::size_t a=0;a<angleKeys.size();++a)for(std::size_t r=0;r<cols;++r)
              localFields[a*traces.size()+first+r][component]=
                  synthesized(static_cast<Eigen::Index>(a),static_cast<Eigen::Index>(r));
          }
        }else{
          // Reorder by m so each DGEMM uses only the real Wigner-d matrix
          // X[angle,l]. Apply exp(+i*m*phi) after the degree contraction.
          std::vector<int> orders;orders.reserve(keys.size());
          for(const auto &[m,l]:keys){(void)l;if(orders.empty()||orders.back()!=m)orders.push_back(m);}
          for(const int m:orders){
            const auto keyBegin=std::lower_bound(keys.begin(),keys.end(),std::pair{m,std::numeric_limits<int>::min()});
            const auto keyEnd=std::upper_bound(keys.begin(),keys.end(),std::pair{m,std::numeric_limits<int>::max()});
            std::vector<int> degrees;degrees.reserve(static_cast<std::size_t>(keyEnd-keyBegin));
            for(auto key=keyBegin;key!=keyEnd;++key)degrees.push_back(key->second);
            const auto basis=realBasisFor(spins[component],m,degrees);
            for(std::size_t first=0;first<traces.size();first+=radialTile){
              const std::size_t cols=std::min(radialTile,traces.size()-first);
              Eigen::MatrixXcd coefficients=Eigen::MatrixXcd::Zero(
                  static_cast<Eigen::Index>(degrees.size()),static_cast<Eigen::Index>(cols));
              for(std::size_t r=0;r<cols;++r){
                const auto *record=byTrace[first+r];if(!record)continue;
                const auto &permutation=recordKeyIndex[labelLayoutForTrace[first+r]];
                const bool fluid=sem_->elementIsFluid(traces[first+r].first);
                for(std::size_t q=0;q<permutation.size();++q){
                  const std::size_t key=permutation[q];
                  if(key==missing||key<static_cast<std::size_t>(keyBegin-keys.begin())||
                     key>=static_cast<std::size_t>(keyEnd-keys.begin()))continue;
                  const auto &values=record->rawStrain[q];
                  coefficients(static_cast<Eigen::Index>(key-static_cast<std::size_t>(keyBegin-keys.begin())),static_cast<Eigen::Index>(r))=
                      fluid?(component==0?values[0]+values[1]:Complex{}):values[component];
                }
              }
              const auto degreeSums=dsfast::synthesize(*basis,coefficients,
                  static_cast<dsfast::DotBackend>(contractionBackend));
              for(std::size_t a=0;a<angleKeys.size();++a){
                const Complex phase=std::exp(Complex(0.0,m*localPhi[a]));
                for(std::size_t r=0;r<cols;++r)
                  localFields[a*traces.size()+first+r][component]+=
                      phase*degreeSums(static_cast<Eigen::Index>(a),static_cast<Eigen::Index>(r));
              }
            }
          }
        }
      }
      {mathprofile::Scope scatter("tensor_phase_spin_rotation_and_scatter");for(std::size_t t=0;t<traces.size();++t)for(std::size_t a=0;a<angleKeys.size();++a){
        const auto p=pointForTrace[t].at(angleKeys[a]);
        const auto &v=localFields[a*traces.size()+t];
        const SpinStrainField local{v[0],v[1],v[2],v[3],v[4],v[5]};
        const auto global=Windowed3DFieldEngine::rotateStrain(local,rotA[a],rotB[a]);
        auto &target=adjoint?result[p].adjointDualStrain:result[p].forwardRawStrain;
        target.rr+=global.rr;target.horizontalTrace+=global.horizontalTrace;
        target.radialPlus+=global.radialPlus;target.radialMinus+=global.radialMinus;
        target.horizontalPlus+=global.horizontalPlus;target.horizontalMinus+=global.horizontalMinus;
      }
      }
    };
    if(!cachedForward)synthesisOne(0,false);
    for(int rec=0;rec<parameters_->num_receivers();++rec)
      synthesisOne(static_cast<std::uint64_t>(rec)+1,true);

    // Radial l=0 is already synthesized by sampleModalBlocks, and is kept as
    // a separate reserved record so it can be added without refactoring it.
    {
    mathprofile::Scope l0Scatter("tensor_radial_l0_scatter");
    for(std::size_t t=0;t<traces.size();++t){
      const auto *fwd=records.contains({std::numeric_limits<std::uint64_t>::max(),false,
          traces[t].first,traces[t].second})?&records.at({std::numeric_limits<std::uint64_t>::max(),false,
          traces[t].first,traces[t].second}):nullptr;
      const auto *adj=records.contains({std::numeric_limits<std::uint64_t>::max()-1,true,
          traces[t].first,traces[t].second})?&records.at({std::numeric_limits<std::uint64_t>::max()-1,true,
          traces[t].first,traces[t].second}):nullptr;
      for(const auto &[ak,p]:pointForTrace[t]){
        (void)ak;
        if(!cachedForward&&fwd&&!fwd->rawStrain.empty()){
          result[p].forwardRawStrain.rr+=fwd->rawStrain[0][0];
          result[p].forwardRawStrain.horizontalTrace+=fwd->rawStrain[0][1];
        }
        if(adj&&!adj->rawStrain.empty()){
          result[p].adjointDualStrain.rr+=adj->rawStrain[0][0];
          result[p].adjointDualStrain.horizontalTrace+=adj->rawStrain[0][1];
        }
      }
    }
    }
    if(cachedForward){
      mathprofile::Scope cacheCopy("tensor_forward_cache_copy_and_receiver_suffix");
      for(std::size_t p=0;p<receiverBegin;++p)
        result[p].forwardRawStrain=cachedForward->volumeStrain[p];
      for(std::size_t i=0;i<cachedForward->receiverDisplacement.size();++i){
        const std::size_t p=receiverBegin+i;
        result[p].radialElement=points[p].radialElement;
        result[p].fluid=sem_->elementIsFluid(points[p].radialElement);
        result[p].forwardRawDisplacement=cachedForward->receiverDisplacement[i];
      }
    }else if(receiverBegin<points.size()){
      mathprofile::Scope receiverSuffix("tensor_receiver_suffix_sampling_and_scatter");
      std::vector<Field3DPoint> receivers(points.begin()+receiverBegin,points.end());
      const auto sampled=sampleForKernelsPacked(receivers,0,samplingPlan,contractionBackend);
      for(std::size_t i=0;i<receivers.size();++i){
        const auto p=receiverBegin+i;
        result[p].radialElement=sampled[i].radialElement;result[p].fluid=sampled[i].fluid;
        result[p].forwardRawDisplacement=sampled[i].forwardRawDisplacement;
      }
    }
    return result;
  }
  const std::vector<OwnedModalFrequencyBlock3D> &modalBlocks() const { return blocks_; }
  Full1D::SEM &sem() const { return *sem_; }
  InputParameters &parameters() const { return *parameters_; }
  const Eigen::VectorXcd &adjointWeights() const { return adjointWeights_; }
  Complex z() const { return z_; }
  std::size_t radialCoefficientCacheBytes() const { return frequencyCoefficientCache_.bytes; }
private:
  Windowed3DFrequencyWorkspace(Full1D::SEM &sem,InputParameters &parameters,
      Complex z,const Eigen::VectorXcd &adjointWeights,
      std::vector<OwnedModalFrequencyBlock3D> blocks,
      std::vector<windowed3d_detail::PoleModalFrequencyBlock3D> poleBlocks,
      windowed3d_detail::PoleFrame3D sourceFrame,bool forwardOnly)
      :sem_(&sem),parameters_(&parameters),z_(z),adjointWeights_(adjointWeights),
       blocks_(std::move(blocks)),poleBlocks_(std::move(poleBlocks)),sourceFrame_(sourceFrame){
    mathprofile::Scope synthesisVectorConstruction("workspace_synthesis_vector_construction");
    const Complex fz=-Complex(0.0,1.0)/z_;
    sourceSynthesis_.reserve(poleBlocks_.size());
    for(const auto &block:poleBlocks_){
      const auto coeff=block.sourceCoefficients.transpose();
      sourceSynthesis_.push_back({block.family,block.degree,block.sourceElement,false,
          block.forward*coeff,block.sourceLift*coeff});
    }
    if(forwardOnly)return;
    adjointSynthesis_.reserve(static_cast<std::size_t>(parameters_->num_receivers()));
    for(int rec=0;rec<parameters_->num_receivers();++rec){
      AdjointSynthesis synthesis;
      synthesis.synthesisSet=static_cast<std::uint64_t>(rec)+1;
      const auto station=parameters_->receivers()[rec];
      synthesis.frame=Windowed3DFieldEngine::poleFrame(station.first,station.second);
      synthesis.blocks.reserve(poleBlocks_.size());
      for(const auto &block:poleBlocks_){
        Eigen::MatrixXcd poleReceiver;
        {mathprofile::Scope t("receiver_at_pole_closed_form");
          poleReceiver=Windowed3DFieldEngine::receiverAtPole(block.family,block.degree);}
        const int mmax=(static_cast<int>(poleReceiver.cols())-1)/2;
        Eigen::MatrixXcd localWeights;
        {mathprofile::Scope t("receiver_local_weight_materialization");
          localWeights=Eigen::MatrixXcd::Zero(3,2*mmax+1);
          for(int comp=0;comp<3;++comp)for(int m=-mmax;m<=mmax;++m)
            localWeights(comp,m+mmax)=std::conj(poleReceiver(comp,m+mmax))*
                (std::conj(fz)*adjointWeights_[3*rec+comp]);}
        Eigen::MatrixXcd localAdjoint;
        if(!block.adjointDirect.empty()){
          if(static_cast<std::size_t>(rec)>=block.adjointDirect.size())
            throw std::runtime_error("Phase C direct adjoint receiver inventory mismatch.");
          mathprofile::Scope t("direct_measurement_adjoint_solution_use");
          localAdjoint=block.adjointDirect[static_cast<std::size_t>(rec)];
        }else{
          mathprofile::Scope t("receiver_basis_mix_and_unit_adjoint_projection");
          const auto receiverBasisMix=block.adjointBasisTransform.middleRows(3*rec,3).adjoint()*localWeights;
          localAdjoint=block.adjointUnit*receiverBasisMix;
        }
        synthesis.blocks.push_back({block.family,block.degree,block.sourceElement,true,
            localAdjoint,Eigen::MatrixXcd{}});
      }
      adjointSynthesis_.push_back(std::move(synthesis));
    }
    synthesisVectorConstruction.stop();
  }

  std::vector<Field3DSample> samplePoleFrames(
      const std::vector<Field3DPoint> &points,bool allowCenter,
      bool kernelRoute=false,std::size_t displacementBegin=0,
      Windowed3DFieldEngine::SamplingPlanCache *samplingPlan=nullptr,
      int contractionBackend=0,
      const std::function<void(const PoleAngularCoefficientRecord3D&)> &coefficientVisitor={}) const {
    mathprofile::Scope all("pole_frame_sampling_total");
    auditprofile::Scope t("sample_full_total");
    std::vector<OwnedModalFrequencyBlock3D> radialBlocks;
    for(const auto &block:blocks_)if(block.family==SensitivitySolveFamily::Radial)
      radialBlocks.push_back(block);
#ifndef DSFAST_DISABLE_PROFILE
    auto radialStart=std::chrono::steady_clock::now();
#endif
    auto result=[&]{mathprofile::Scope t("radial_l0_sampling_and_interpolation");return Windowed3DFieldEngine::sampleModalBlocks(
        *sem_,*parameters_,z_,radialBlocks,points,allowCenter,
        kernelRoute,displacementBegin);}();
    if(kernelRoute&&coefficientVisitor){
      std::set<std::pair<int,std::uint64_t>> radialTraces;
      for(std::size_t p=0;p<points.size();++p){
        const double radius=Windowed3DFieldEngine::canonicalSamplingRadius(
            *sem_,result[p].radialElement,points[p].radiusMetres/sem_->lengthScaleMetres());
        const auto key=std::pair{result[p].radialElement,std::bit_cast<std::uint64_t>(radius)};
        if(!radialTraces.insert(key).second)continue;
        for(int adj=0;adj<2;++adj){
          PoleAngularCoefficientRecord3D record;
          record.frame.axes={{{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}}};
          record.synthesisSet=adj?std::numeric_limits<std::uint64_t>::max()-1:
                                    std::numeric_limits<std::uint64_t>::max();
          record.adjoint=adj!=0;
          record.element=result[p].radialElement;
          record.radiusModel=radius;
          record.labels.emplace_back(0,0);
          const auto &strain=adj?result[p].adjointDualStrain:result[p].forwardRawStrain;
          record.rawStrain.push_back({strain.rr,strain.horizontalTrace,Complex{},Complex{},Complex{},Complex{}});
          coefficientVisitor(record);
        }
      }
    }
#ifndef DSFAST_DISABLE_PROFILE
    auditprofile::data["radial_sampling_and_point_validation"].seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-radialStart).count();
    ++auditprofile::data["radial_sampling_and_point_validation"].calls;
#endif
    if(points.empty())return result;
    const double scale=sem_->lengthScaleMetres();
    const auto tol=[](double r,double b){return 8.0*std::numeric_limits<double>::epsilon()*
        std::max({1.0,std::abs(r),std::abs(b)});};
    std::vector<double> radii(points.size());std::vector<int> elements(points.size());
    for(std::size_t p=0;p<points.size();++p){
      double r=points[p].radiusMetres/scale;const int e=result[p].radialElement;
      for(double b:{sem_->mesh().ELR(e),sem_->mesh().EUR(e)})
        if(b>0.0&&std::abs(r-b)<=tol(r,b))r=b;
      radii[p]=r;elements[p]=e;
    }
    const Complex fz=-Complex(0.0,1.0)/z_;
    for(const auto &adjoint:adjointSynthesis_){
      {mathprofile::Scope t("receiver_frame_geometry_phase_spin_rotation_and_scatter");Windowed3DFieldEngine::synthesizePoleFrame(*sem_,adjoint.blocks,
          adjoint.frame,points,radii,elements,scale,result,kernelRoute,displacementBegin,
          samplingPlan,&frequencyCoefficientCache_,adjoint.synthesisSet,contractionBackend,coefficientVisitor);}
    }
    {mathprofile::Scope t("source_frame_geometry_phase_spin_rotation_and_scatter");Windowed3DFieldEngine::synthesizePoleFrame(*sem_,sourceSynthesis_,
        sourceFrame_,points,radii,elements,scale,result,kernelRoute,displacementBegin,
        samplingPlan,&frequencyCoefficientCache_,0,contractionBackend,coefficientVisitor);}
    for(auto &sample:result)for(int c=0;c<3;++c)
      sample.observableForwardDisplacementSpectrumMetreSeconds[c]=
          sem_->timeScaleSeconds()*fz*scale*sample.forwardRawDisplacement[c];
    return result;
  }

  Full1D::SEM *sem_;
  InputParameters *parameters_;
  Complex z_;
  Eigen::VectorXcd adjointWeights_;
  std::vector<OwnedModalFrequencyBlock3D> blocks_;
  std::vector<windowed3d_detail::PoleModalFrequencyBlock3D> poleBlocks_;
  windowed3d_detail::PoleFrame3D sourceFrame_;
  std::vector<windowed3d_detail::PoleSynthesisBlock3D> sourceSynthesis_;
  std::vector<AdjointSynthesis> adjointSynthesis_;
  mutable Windowed3DFieldEngine::PoleCoefficientCache frequencyCoefficientCache_;
};

inline Windowed3DFrequencyWorkspace Windowed3DFieldEngine::prepareFrequencyWorkspace(
    Full1D::SEM &sem,SourceInfo::EarthquakeCMT &source,InputParameters &parameters,
    Complex z,const Eigen::VectorXcd &adjointWeights,
    dssolver::BatchResponseCache *batchCache,bool poleOnly,bool forwardOnly,
    double *basisReconstructionRelative,int solveThreads,std::optional<Complex> attenuation) {
  return Windowed3DFrequencyWorkspace::prepare(sem,source,parameters,z,adjointWeights,batchCache,poleOnly,forwardOnly,
      basisReconstructionRelative,solveThreads,attenuation);
}

inline std::vector<Windowed3DFieldEngine::Sample>
Windowed3DFieldEngine::evaluateFrequency(Full1D::SEM &sem,
    SourceInfo::EarthquakeCMT &source,InputParameters &parameters,Complex z,
    const Eigen::VectorXcd &adjointWeights,const std::vector<Point> &points) {
  auto workspace=prepareFrequencyWorkspace(sem,source,parameters,z,adjointWeights);
  return workspace.sample(points);
}

inline std::vector<Windowed3DFieldEngine::Sample>
Windowed3DFieldEngine::evaluateFrequencyForWeightedVolume(Full1D::SEM &sem,
    SourceInfo::EarthquakeCMT &source,InputParameters &parameters,Complex z,
    const Eigen::VectorXcd &adjointWeights,const std::vector<Point> &points) {
  auto workspace=prepareFrequencyWorkspace(sem,source,parameters,z,adjointWeights);
  return workspace.sampleWeightedChunk(points);
}

} // namespace SPARSESPEC

#endif
