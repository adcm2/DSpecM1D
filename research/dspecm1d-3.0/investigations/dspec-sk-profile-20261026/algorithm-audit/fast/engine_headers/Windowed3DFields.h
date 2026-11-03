#ifndef DSPECM1D_WINDOWED_3D_FIELDS_H
#define DSPECM1D_WINDOWED_3D_FIELDS_H

#include "SensitivityKernels.h"
#include "../fast_engine.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <complex>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <tuple>
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
    std::map<FrameKey,std::map<PointKey,PointAngle>> pointAngles;
    std::map<std::vector<std::uint64_t>,std::shared_ptr<CachedAngularWigner>> angularTables;
    std::map<std::vector<std::uint64_t>,std::shared_ptr<AngularCoefficientTable>> angularCoefficients;
    // Packed synthesis matrices depend only on the semantic angle grid, spin,
    // and ordered (m,l) columns, so reuse them across radial traces/frequencies.
    std::map<std::vector<std::uint64_t>,std::shared_ptr<Eigen::MatrixXcd>> packedBases;
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
      dssolver::BatchResponseCache *batchCache=nullptr);

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
      const std::size_t perAngle=angularTableBytesPerAngle(maxAngularDegree);
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
    static constexpr std::size_t byteBudget=32u*1024u*1024u;
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

  static std::size_t angularTableBytesPerAngle(int maxDegree) {
    const int nMax=std::min(2,maxDegree);
    std::size_t values=0;
    for(int n=-nMax;n<=nMax;++n)
      for(int degree=std::abs(n);degree<=maxDegree;++degree)
        values+=static_cast<std::size_t>(2*degree+1);
    return values*sizeof(double)+static_cast<std::size_t>(2*nMax+1)*sizeof(std::size_t);
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
    for(int m=-mmax;m<=mmax;++m){
      const int col=m+mmax;
      if(family==SensitivitySolveFamily::Spheroidal)
        result(0,col)=harmonic(degree,m,0,0.0,0.0);
      const Complex ym=harmonic(degree,m,-1,0.0,0.0);
      const Complex yp=harmonic(degree,m,1,0.0,0.0);
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

  static void solveFrequencyModes(Full1D::SEM &sem,
      SourceInfo::EarthquakeCMT &source,InputParameters &parameters,Complex z,
      const Eigen::VectorXcd &adjointWeights,const ModalVisitor &visitor,
      std::vector<PoleModalBlock> *poleBlocks=nullptr,
      PoleFrame *sourceFrame=nullptr,
      dssolver::BatchResponseCache *batchCache=nullptr) {
    if(!sem.sourceMatches(source)||parameters.attenuation()||
       parameters.cowlingFrequencyMhz()!=0.0||parameters.gravitation()!=2||
       parameters.output_type()!=0||!std::isfinite(z.real())||
       !std::isfinite(z.imag())||z==Complex{}||
       adjointWeights.size()!=3*parameters.num_receivers()||!adjointWeights.allFinite())
      throw std::invalid_argument("Invalid one-frequency 3D modal workspace input.");
    const int lmax=parameters.lmax();
    if(parameters.lmin()<0||lmax<parameters.lmin()||lmax>sem.maxDegree())
      throw std::invalid_argument("3D modal degree range is invalid.");
    const auto flags=resolveModeFlags(parameters.type(),parameters.lmin(),lmax);
    if(!flags.inc_rad&&!flags.inc_sph&&!flags.inc_tor)
      throw std::invalid_argument("3D field mode selection is empty.");
    const Complex fz=-Complex(0.0,1.0)/z;
    if(sourceFrame)*sourceFrame=poleFrame(source.Latitude(),source.Longitude());
    ParamInfo angular(parameters,lmax);
    if(flags.inc_rad){
      auto treatment=sem.sourceTreatedRadial(source,z);
      Eigen::MatrixXcd receiver=Eigen::MatrixXcd::Zero(
          3*parameters.num_receivers(),treatment.force.rows());
      for(int rec=0;rec<parameters.num_receivers();++rec)
        receiver.row(3*rec)=sem.rvZR(parameters,rec).transpose();
      const Eigen::MatrixXcd rhs=receiver.adjoint()*(std::conj(fz)*adjointWeights);
      const auto receiverBasis=dsfast::exactReceiverRowBasis(receiver);
      Eigen::MatrixXcd forward,adjoint,unusedBasis;
      if(batchCache){
        const auto key=dssolver::BatchResponseKey::make(&sem,
            static_cast<int>(SensitivitySolveFamily::Radial),0,z);
        const auto response=batchCache->getOrCompute(key,treatment.force,
            receiverBasis.rightHandSides,[&](const auto &forwardRhs,const auto &basisRhs){
              Eigen::MatrixXcd f,a,b,empty(0,0);
              FixedMeshAdjointSolver::solveForwardAdjointFrequency(
                  sem,SensitivitySolveFamily::Radial,0,z,forwardRhs,empty,basisRhs,f,a,b);
              return std::pair{std::move(f),std::move(b)};
            });
        const auto weighted=batchCache->deriveWeightedAdjoint(
            key,receiverBasis.rightHandSides,rhs);
        forward=response.forward;adjoint=weighted.weightedAdjoint;unusedBasis=response.unitAdjoint;
      }else FixedMeshAdjointSolver::solveForwardAdjointFrequency(
          sem,SensitivitySolveFamily::Radial,0,z,treatment.force,rhs,{},
          forward,adjoint,unusedBasis);
      const Eigen::MatrixXcd noData;
      visitor({SensitivitySolveFamily::Radial,0,treatment.sourceElement,
          forward,adjoint,treatment.lift,noData,noData});
    }
    for(int degree=std::max(1,parameters.lmin());degree<=lmax;++degree){
      const int orders=2*degree+1;
      if(flags.inc_sph){
        auto treatment=sem.sourceTreatedSpheroidal(source,degree,z);
        const Eigen::MatrixXcd sourceCoefficients=sem.calculateForceCoefficients(source,degree);
        const Eigen::MatrixXcd modeReceiver=angular.rvFullSph(degree);
        const auto receivers=sem.receiverElements(parameters,false);
        const auto receiverBase=sem.rvBaseFull(parameters,degree,false);
        const auto receiver=embedReceiver(sem,receiverBase,
            sem.ltgS(0,receivers.front(),0),static_cast<int>(sem.hS(degree).rows()));
        Eigen::MatrixXcd rhs=Eigen::MatrixXcd::Zero(receiver.cols(),orders);
        for(int m=-degree;m<=degree;++m){
          Eigen::MatrixXcd weighted=receiver;
          for(int row=0;row<weighted.rows();++row)weighted.row(row)*=modeReceiver(row,m+degree);
          rhs.col(m+degree)=weighted.adjoint()*(std::conj(fz)*adjointWeights);
        }
        const auto receiverBasis=dsfast::exactReceiverRowBasis(receiver);
        Eigen::MatrixXcd adjoint,adjointUnit;
        Eigen::MatrixXcd forward;
        if(batchCache){
          const auto key=dssolver::BatchResponseKey::make(&sem,
              static_cast<int>(SensitivitySolveFamily::Spheroidal),degree,z);
          const auto response=batchCache->getOrCompute(key,treatment.force,
              receiverBasis.rightHandSides,[&](const auto &forwardRhs,const auto &basisRhs){
                Eigen::MatrixXcd f,a,b,empty(0,0);
                FixedMeshAdjointSolver::solveForwardAdjointFrequency(
                    sem,SensitivitySolveFamily::Spheroidal,degree,z,forwardRhs,empty,basisRhs,f,a,b);
                return std::pair{std::move(f),std::move(b)};
              });
          const auto weighted=batchCache->deriveWeightedAdjoint(
              key,receiverBasis.rightHandSides,rhs);
          forward=response.forward;adjoint=weighted.weightedAdjoint;adjointUnit=response.unitAdjoint;
        }else FixedMeshAdjointSolver::solveForwardAdjointFrequency(
            sem,SensitivitySolveFamily::Spheroidal,degree,z,treatment.force,rhs,
            receiverBasis.rightHandSides,forward,adjoint,adjointUnit);
        visitor({SensitivitySolveFamily::Spheroidal,degree,treatment.sourceElement,
          forward,adjoint,treatment.lift,sourceCoefficients,modeReceiver});
        if(poleBlocks){
          const auto reduced=sem.calculateForceRedCoefficients(source,degree,0.0);
          poleBlocks->push_back({SensitivitySolveFamily::Spheroidal,degree,
              treatment.sourceElement,forward,adjointUnit,treatment.lift,reduced,
              receiverBasis.reconstruction});
        }
      }
      if(flags.inc_tor&&!sem.elementIsFluid(sem.sourceElement(source))){
        auto treatment=sem.sourceTreatedToroidal(source,degree,z);
        const Eigen::MatrixXcd sourceCoefficients=sem.calculateForceCoefficientsT(source,degree);
        const Eigen::MatrixXcd modeReceiver=angular.rvFullTor(degree);
        const auto receivers=sem.receiverElements(parameters,true);
        const auto receiverBase=sem.rvBaseFullT(parameters,degree);
        const auto receiver=embedReceiver(sem,receiverBase,
            sem.ltgT(receivers.front(),0),static_cast<int>(sem.hTk(degree).rows()));
        Eigen::MatrixXcd rhs=Eigen::MatrixXcd::Zero(receiver.cols(),orders);
        for(int m=-degree;m<=degree;++m){
          Eigen::MatrixXcd weighted=receiver;
          for(int row=0;row<weighted.rows();++row)weighted.row(row)*=modeReceiver(row,m+degree);
          rhs.col(m+degree)=weighted.adjoint()*(std::conj(fz)*adjointWeights);
        }
        const auto receiverBasis=dsfast::exactReceiverRowBasis(receiver);
        Eigen::MatrixXcd adjoint,adjointUnit;
        Eigen::MatrixXcd forward;
        if(batchCache){
          const auto key=dssolver::BatchResponseKey::make(&sem,
              static_cast<int>(SensitivitySolveFamily::Toroidal),degree,z);
          const auto response=batchCache->getOrCompute(key,treatment.force,
              receiverBasis.rightHandSides,[&](const auto &forwardRhs,const auto &basisRhs){
                Eigen::MatrixXcd f,a,b,empty(0,0);
                FixedMeshAdjointSolver::solveForwardAdjointFrequency(
                    sem,SensitivitySolveFamily::Toroidal,degree,z,forwardRhs,empty,basisRhs,f,a,b);
                return std::pair{std::move(f),std::move(b)};
              });
          const auto weighted=batchCache->deriveWeightedAdjoint(
              key,receiverBasis.rightHandSides,rhs);
          forward=response.forward;adjoint=weighted.weightedAdjoint;adjointUnit=response.unitAdjoint;
        }else FixedMeshAdjointSolver::solveForwardAdjointFrequency(
            sem,SensitivitySolveFamily::Toroidal,degree,z,treatment.force,rhs,
            receiverBasis.rightHandSides,forward,adjoint,adjointUnit);
        visitor({SensitivitySolveFamily::Toroidal,degree,treatment.sourceElement,
          forward,adjoint,treatment.lift,sourceCoefficients,modeReceiver});
        if(poleBlocks){
          const auto reduced=sem.calculateForceRedCoefficientsT(source,degree,0.0);
          poleBlocks->push_back({SensitivitySolveFamily::Toroidal,degree,
              treatment.sourceElement,forward,adjointUnit,treatment.lift,reduced,
              receiverBasis.reconstruction});
        }
      }
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
        const auto aRadial=interpolate(sem,block.family,block.adjoint,
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
    const std::size_t perAngle=angularTableBytesPerAngle(maxDegree);
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
        if(tableIt!=samplingPlan->angularTables.end()){angularTable=tableIt->second;++auditprofile::counters["wigner_hit"];}
      }
      if(!angularTable){auditprofile::Scope t("wigner_build");++auditprofile::counters["wigner_miss"];
        angularTable=std::make_shared<CachedAngularWigner>(
            maxDegree,std::min(2,maxDegree),std::min(2,maxDegree),chunks[chunk].theta);
        const std::size_t bytes=angularTableBytesPerAngle(maxDegree)*chunks[chunk].theta.size();
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
          {angularCoefficients=coefficientsIt->second;++auditprofile::counters["angular_coefficient_hit"];}
      }
      if(!angularCoefficients){auditprofile::Scope t("angular_coefficient_build");++auditprofile::counters["angular_coefficient_miss"];
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
          if(cached!=frequencyCoefficientCache->rows.end()){coefficientsPtr=&cached->second;++auditprofile::counters["radial_coefficients_hit"];}
        }
        const auto buildCoefficients=[&](PoleCoefficientRows &coefficients){auditprofile::Scope t("radial_coefficients_build");++auditprofile::counters["radial_coefficients_miss"];
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
  struct AdjointSynthesis {
    windowed3d_detail::PoleFrame3D frame;
    std::vector<windowed3d_detail::PoleSynthesisBlock3D> blocks;
    std::uint64_t synthesisSet=0;
  };
  static Windowed3DFrequencyWorkspace prepare(Full1D::SEM &sem,
      SourceInfo::EarthquakeCMT &source,InputParameters &parameters,Complex z,
      const Eigen::VectorXcd &adjointWeights,
      dssolver::BatchResponseCache *batchCache=nullptr) {
    sensitivity_detail::requireConstructedSEM(sem);
    auditprofile::Scope t("prepare_workspace");
    std::vector<OwnedModalFrequencyBlock3D> blocks;
    std::vector<windowed3d_detail::PoleModalFrequencyBlock3D> poleBlocks;
    windowed3d_detail::PoleFrame3D sourceFrame;
    const auto capture=[&](const ModalFrequencyBlock3D &b){
      blocks.push_back({b.family,b.degree,b.sourceElement,b.forward,b.adjoint,
                        b.sourceLift,b.sourceCoefficients,b.modeReceiver});
    };
    Windowed3DFieldEngine::solveFrequencyModes(
        sem,source,parameters,z,adjointWeights,capture,&poleBlocks,&sourceFrame,batchCache);
    return Windowed3DFrequencyWorkspace(sem,parameters,z,adjointWeights,
        std::move(blocks),std::move(poleBlocks),sourceFrame);
  }
  std::vector<Field3DSample> sample(const std::vector<Field3DPoint> &points) const {
    return samplePoleFrames(points,false);
  }
  std::vector<Field3DSample> sampleWeightedChunk(
      const std::vector<Field3DPoint> &points) const {
    return samplePoleFrames(points,true);
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
    std::vector<KernelField3DSample> compact(full.size());
    for(std::size_t p=0;p<full.size();++p){
      compact[p].radialElement=full[p].radialElement;compact[p].fluid=full[p].fluid;
      compact[p].forwardRawDisplacement=full[p].forwardRawDisplacement;
      compact[p].forwardRawStrain=full[p].forwardRawStrain;
      compact[p].adjointDualStrain=full[p].adjointDualStrain;
    }
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
      bool realDegreeLayout=false) const {
    if(receiverBegin>points.size()||contractionBackend<1||contractionBackend>3)
      throw std::invalid_argument("Tensor kernel sample arguments are invalid.");
    if(receiverBegin==0)throw std::invalid_argument("Tensor sample requires a receiver suffix.");
    using TraceKey=std::pair<int,std::uint64_t>;
    using AngleKey=std::pair<std::uint64_t,std::uint64_t>;
    const auto traceKey=[&](const Field3DPoint &p){return TraceKey{p.radialElement,
        std::bit_cast<std::uint64_t>(p.radiusMetres/sem_->lengthScaleMetres())};};
    const auto angleKey=[](const Field3DPoint &p){return AngleKey{
        std::bit_cast<std::uint64_t>(p.thetaRadians),std::bit_cast<std::uint64_t>(p.phiRadians)};};
    std::vector<TraceKey> traces;
    std::map<TraceKey,std::size_t> traceIndex;
    std::vector<AngleKey> angleKeys;
    std::map<AngleKey,std::size_t> angleIndex;
    std::vector<std::map<AngleKey,std::size_t>> pointForTrace;
    for(std::size_t p=0;p<receiverBegin;++p){
      const auto tk=traceKey(points[p]);
      auto [ti,added]=traceIndex.try_emplace(tk,traces.size());
      if(added){traces.push_back(tk);pointForTrace.emplace_back();}
      const auto ak=angleKey(points[p]);
      auto [ai,newAngle]=angleIndex.try_emplace(ak,angleKeys.size());
      if(newAngle)angleKeys.push_back(ak);
      auto [pi,pointAdded]=pointForTrace[ti->second].try_emplace(ak,p);
      (void)pi;
      if(!pointAdded)throw std::invalid_argument("Tensor grid has duplicate angle samples on a radial trace.");
    }
    if(traces.empty()||angleKeys.empty())throw std::invalid_argument("Tensor grid is empty.");
    for(const auto &grid:pointForTrace)if(grid.size()!=angleKeys.size())
      throw std::invalid_argument("Tensor grid is not a shared angle-set by radial-trace product.");

    std::vector<Field3DPoint> representatives;
    representatives.reserve(traces.size());
    for(std::size_t t=0;t<traces.size();++t)
      representatives.push_back(points[pointForTrace[t].begin()->second]);
    std::map<std::tuple<std::uint64_t,bool,int,std::uint64_t>,PoleAngularCoefficientRecord3D> records;
    const auto visitor=[&](const PoleAngularCoefficientRecord3D &record){
      const auto key=std::tuple{record.synthesisSet,record.adjoint,record.element,
          std::bit_cast<std::uint64_t>(record.radiusModel)};
      records.try_emplace(key,record);
      if(coefficientVisitor)coefficientVisitor(record);
    };
    (void)samplePoleFrames(representatives,false,true,representatives.size(),
        samplingPlan,contractionBackend,visitor);

    std::vector<KernelField3DSample> result(points.size());
    for(std::size_t p=0;p<receiverBegin;++p){
      result[p].radialElement=points[p].radialElement;
      result[p].fluid=sem_->elementIsFluid(points[p].radialElement);
    }
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
      std::vector<double> localTheta(angleKeys.size()),localPhi(angleKeys.size()),rotA(angleKeys.size()),rotB(angleKeys.size());
      for(std::size_t a=0;a<angleKeys.size();++a){
        const auto &point=points[pointForTrace.front().at(angleKeys[a])];
        double aa=1.0,bb=0.0;
        const auto local=Windowed3DFieldEngine::tangentRotation(frame,
            point.thetaRadians,point.phiRadians,aa,bb);
        localTheta[a]=local[0];localPhi[a]=local[1];rotA[a]=aa;rotB[a]=bb;
      }
      Windowed3DFieldEngine::CachedAngularWigner wigner(maxDegree,
          std::min(2,maxDegree),std::min(2,maxDegree),localTheta);
      std::vector<std::array<Complex,6>> localFields(angleKeys.size()*traces.size());
      const std::array<int,6> spins{0,0,1,-1,2,-2};
      for(std::size_t component=0;component<6;++component){
        std::set<std::pair<int,int>> keySet;
        for(const auto *record:byTrace)if(record)
          for(const auto &[m,l]:record->labels)
            if(l>=std::abs(spins[component]))keySet.emplace(m,l);
        if(keySet.empty())continue;
        std::vector<std::pair<int,int>> keys(keySet.begin(),keySet.end());
        if(!realDegreeLayout){
          Eigen::MatrixXcd basis(static_cast<Eigen::Index>(angleKeys.size()),
                                 static_cast<Eigen::Index>(keys.size()));
          for(std::size_t a=0;a<angleKeys.size();++a)for(std::size_t q=0;q<keys.size();++q){
            const auto [m,l]=keys[q];
            const double d=wigner[spins[component],static_cast<int>(a)][l,m];
            basis(static_cast<Eigen::Index>(a),static_cast<Eigen::Index>(q))=
                d*std::exp(Complex(0.0,m*localPhi[a]));
          }
          for(std::size_t first=0;first<traces.size();first+=radialTile){
            const std::size_t cols=std::min(radialTile,traces.size()-first);
            Eigen::MatrixXcd coefficients=Eigen::MatrixXcd::Zero(
                static_cast<Eigen::Index>(keys.size()),static_cast<Eigen::Index>(cols));
            for(std::size_t r=0;r<cols;++r){
              const auto *record=byTrace[first+r];if(!record)continue;
              std::map<std::pair<int,int>,std::size_t> recordIndex;
              for(std::size_t q=0;q<record->labels.size();++q)recordIndex.emplace(record->labels[q],q);
              const bool fluid=sem_->elementIsFluid(traces[first+r].first);
              for(std::size_t q=0;q<keys.size();++q){
                const auto it=recordIndex.find(keys[q]);if(it==recordIndex.end())continue;
                const auto &values=record->rawStrain[it->second];
                coefficients(static_cast<Eigen::Index>(q),static_cast<Eigen::Index>(r))=
                    fluid?(component==0?values[0]+values[1]:component==1?Complex{}:Complex{}):
                          values[component];
              }
            }
            const auto synthesized=dsfast::synthesizeComplex(basis,coefficients,
                static_cast<dsfast::DotBackend>(contractionBackend));
            for(std::size_t a=0;a<angleKeys.size();++a)for(std::size_t r=0;r<cols;++r)
              localFields[a*traces.size()+first+r][component]=
                  synthesized(static_cast<Eigen::Index>(a),static_cast<Eigen::Index>(r));
          }
        }else{
          // Reorder by m so each DGEMM uses only the real Wigner-d matrix
          // X[angle,l]. Apply exp(+i*m*phi) after the degree contraction.
          std::set<int> orders;for(const auto &[m,l]:keys){(void)l;orders.insert(m);}
          for(const int m:orders){
            std::vector<int> degrees;
            for(const auto &[km,l]:keys)if(km==m)degrees.push_back(l);
            Eigen::MatrixXd basis(static_cast<Eigen::Index>(angleKeys.size()),
                                  static_cast<Eigen::Index>(degrees.size()));
            for(std::size_t a=0;a<angleKeys.size();++a)
              for(std::size_t q=0;q<degrees.size();++q)
                basis(static_cast<Eigen::Index>(a),static_cast<Eigen::Index>(q))=
                    wigner[spins[component],static_cast<int>(a)][degrees[q],m];
            for(std::size_t first=0;first<traces.size();first+=radialTile){
              const std::size_t cols=std::min(radialTile,traces.size()-first);
              Eigen::MatrixXcd coefficients=Eigen::MatrixXcd::Zero(
                  static_cast<Eigen::Index>(degrees.size()),static_cast<Eigen::Index>(cols));
              for(std::size_t r=0;r<cols;++r){
                const auto *record=byTrace[first+r];if(!record)continue;
                std::map<std::pair<int,int>,std::size_t> recordIndex;
                for(std::size_t q=0;q<record->labels.size();++q)recordIndex.emplace(record->labels[q],q);
                const bool fluid=sem_->elementIsFluid(traces[first+r].first);
                for(std::size_t q=0;q<degrees.size();++q){
                  const auto it=recordIndex.find({m,degrees[q]});if(it==recordIndex.end())continue;
                  const auto &values=record->rawStrain[it->second];
                  coefficients(static_cast<Eigen::Index>(q),static_cast<Eigen::Index>(r))=
                      fluid?(component==0?values[0]+values[1]:Complex{}):values[component];
                }
              }
              const auto degreeSums=dsfast::synthesize(basis,coefficients,
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
      for(std::size_t t=0;t<traces.size();++t)for(std::size_t a=0;a<angleKeys.size();++a){
        const auto p=pointForTrace[t].at(angleKeys[a]);
        const auto &v=localFields[a*traces.size()+t];
        const SpinStrainField local{v[0],v[1],v[2],v[3],v[4],v[5]};
        const auto global=Windowed3DFieldEngine::rotateStrain(local,rotA[a],rotB[a]);
        auto &target=adjoint?result[p].adjointDualStrain:result[p].forwardRawStrain;
        target.rr+=global.rr;target.horizontalTrace+=global.horizontalTrace;
        target.radialPlus+=global.radialPlus;target.radialMinus+=global.radialMinus;
        target.horizontalPlus+=global.horizontalPlus;target.horizontalMinus+=global.horizontalMinus;
      }
    };
    synthesisOne(0,false);
    for(int rec=0;rec<parameters_->num_receivers();++rec)
      synthesisOne(static_cast<std::uint64_t>(rec)+1,true);

    // Radial l=0 is already synthesized by sampleModalBlocks, and is kept as
    // a separate reserved record so it can be added without refactoring it.
    for(std::size_t t=0;t<traces.size();++t){
      const auto *fwd=records.contains({std::numeric_limits<std::uint64_t>::max(),false,
          traces[t].first,traces[t].second})?&records.at({std::numeric_limits<std::uint64_t>::max(),false,
          traces[t].first,traces[t].second}):nullptr;
      const auto *adj=records.contains({std::numeric_limits<std::uint64_t>::max()-1,true,
          traces[t].first,traces[t].second})?&records.at({std::numeric_limits<std::uint64_t>::max()-1,true,
          traces[t].first,traces[t].second}):nullptr;
      for(const auto &[ak,p]:pointForTrace[t]){
        (void)ak;
        if(fwd&&!fwd->rawStrain.empty()){
          result[p].forwardRawStrain.rr+=fwd->rawStrain[0][0];
          result[p].forwardRawStrain.horizontalTrace+=fwd->rawStrain[0][1];
        }
        if(adj&&!adj->rawStrain.empty()){
          result[p].adjointDualStrain.rr+=adj->rawStrain[0][0];
          result[p].adjointDualStrain.horizontalTrace+=adj->rawStrain[0][1];
        }
      }
    }
    if(receiverBegin<points.size()){
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
private:
  Windowed3DFrequencyWorkspace(Full1D::SEM &sem,InputParameters &parameters,
      Complex z,const Eigen::VectorXcd &adjointWeights,
      std::vector<OwnedModalFrequencyBlock3D> blocks,
      std::vector<windowed3d_detail::PoleModalFrequencyBlock3D> poleBlocks,
      windowed3d_detail::PoleFrame3D sourceFrame)
      :sem_(&sem),parameters_(&parameters),z_(z),adjointWeights_(adjointWeights),
       blocks_(std::move(blocks)),poleBlocks_(std::move(poleBlocks)),sourceFrame_(sourceFrame){
    const Complex fz=-Complex(0.0,1.0)/z_;
    sourceSynthesis_.reserve(poleBlocks_.size());
    for(const auto &block:poleBlocks_){
      const auto coeff=block.sourceCoefficients.transpose();
      sourceSynthesis_.push_back({block.family,block.degree,block.sourceElement,false,
          block.forward*coeff,block.sourceLift*coeff});
    }
    adjointSynthesis_.reserve(static_cast<std::size_t>(parameters_->num_receivers()));
    for(int rec=0;rec<parameters_->num_receivers();++rec){
      AdjointSynthesis synthesis;
      synthesis.synthesisSet=static_cast<std::uint64_t>(rec)+1;
      const auto station=parameters_->receivers()[rec];
      synthesis.frame=Windowed3DFieldEngine::poleFrame(station.first,station.second);
      synthesis.blocks.reserve(poleBlocks_.size());
      for(const auto &block:poleBlocks_){
        const auto poleReceiver=Windowed3DFieldEngine::receiverAtPole(block.family,block.degree);
        const int mmax=(static_cast<int>(poleReceiver.cols())-1)/2;
        Eigen::MatrixXcd localWeights=Eigen::MatrixXcd::Zero(3,2*mmax+1);
        for(int comp=0;comp<3;++comp)for(int m=-mmax;m<=mmax;++m)
          localWeights(comp,m+mmax)=std::conj(poleReceiver(comp,m+mmax))*
              (std::conj(fz)*adjointWeights_[3*rec+comp]);
        const auto receiverBasisMix=block.adjointBasisTransform.middleRows(3*rec,3).adjoint()*localWeights;
        const auto localAdjoint=block.adjointUnit*receiverBasisMix;
        synthesis.blocks.push_back({block.family,block.degree,block.sourceElement,true,
            localAdjoint,Eigen::MatrixXcd{}});
      }
      adjointSynthesis_.push_back(std::move(synthesis));
    }
  }

  std::vector<Field3DSample> samplePoleFrames(
      const std::vector<Field3DPoint> &points,bool allowCenter,
      bool kernelRoute=false,std::size_t displacementBegin=0,
      Windowed3DFieldEngine::SamplingPlanCache *samplingPlan=nullptr,
      int contractionBackend=0,
      const std::function<void(const PoleAngularCoefficientRecord3D&)> &coefficientVisitor={}) const {
    auditprofile::Scope t("sample_full_total");
    std::vector<OwnedModalFrequencyBlock3D> radialBlocks;
    for(const auto &block:blocks_)if(block.family==SensitivitySolveFamily::Radial)
      radialBlocks.push_back(block);
#ifndef DSFAST_DISABLE_PROFILE
    auto radialStart=std::chrono::steady_clock::now();
#endif
    auto result=Windowed3DFieldEngine::sampleModalBlocks(
        *sem_,*parameters_,z_,radialBlocks,points,allowCenter,
        kernelRoute,displacementBegin);
    if(kernelRoute&&coefficientVisitor){
      std::set<std::pair<int,std::uint64_t>> radialTraces;
      for(std::size_t p=0;p<points.size();++p){
        const auto key=std::pair{points[p].radialElement,
            std::bit_cast<std::uint64_t>(points[p].radiusMetres/sem_->lengthScaleMetres())};
        if(!radialTraces.insert(key).second)continue;
        for(int adj=0;adj<2;++adj){
          PoleAngularCoefficientRecord3D record;
          record.frame.axes={{{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}}};
          record.synthesisSet=adj?std::numeric_limits<std::uint64_t>::max()-1:
                                    std::numeric_limits<std::uint64_t>::max();
          record.adjoint=adj!=0;
          record.element=points[p].radialElement;
          record.radiusModel=points[p].radiusMetres/sem_->lengthScaleMetres();
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
      Windowed3DFieldEngine::synthesizePoleFrame(*sem_,adjoint.blocks,
          adjoint.frame,points,radii,elements,scale,result,kernelRoute,displacementBegin,
          samplingPlan,&frequencyCoefficientCache_,adjoint.synthesisSet,contractionBackend,coefficientVisitor);
    }
    Windowed3DFieldEngine::synthesizePoleFrame(*sem_,sourceSynthesis_,
        sourceFrame_,points,radii,elements,scale,result,kernelRoute,displacementBegin,
        samplingPlan,&frequencyCoefficientCache_,0,contractionBackend,coefficientVisitor);
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
    dssolver::BatchResponseCache *batchCache) {
  return Windowed3DFrequencyWorkspace::prepare(sem,source,parameters,z,adjointWeights,batchCache);
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
