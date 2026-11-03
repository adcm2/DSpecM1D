#ifndef DSPECM1D_WINDOWED_3D_FIELDS_H
#define DSPECM1D_WINDOWED_3D_FIELDS_H

#include "SensitivityKernels.h"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <complex>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <stdexcept>
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
      InputParameters &parameters,Complex z,const Eigen::VectorXcd &adjointWeights);

  static std::vector<Sample> sampleModalBlocks(
      Full1D::SEM &sem,InputParameters &parameters,Complex z,
      const std::vector<OwnedModalFrequencyBlock3D> &blocks,
      const std::vector<Point> &points,bool allowCenter=false) {
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
        accumulateModalFields(sem,block,points,radii,elements,radialGroups,scale,result,nullptr,0,{});
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
              scale,result,&w,angleBegin,angleForPoint);
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
  using AngularWigner = GSHTrans::Wigner<double,GSHTrans::All,GSHTrans::All,
                                           GSHTrans::Multiple>;

  struct AngularPointChunk {
    std::vector<double> theta;
    std::vector<std::vector<std::size_t>> radialGroups;
    std::size_t lastRadialGroup=std::numeric_limits<std::size_t>::max();
  };

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
    const auto a = mode(values); const double zeta = std::sqrt(l * (l + 1.0));
    return {a.U*x0, zeta/2.0*a.V*(xp-xm), -Complex(0,1)*zeta/2.0*a.V*(xp+xm)};
  }

  static std::array<Complex,3> toroidalDisplacement(int l,
      const Eigen::MatrixXcd &values, Complex xp, Complex xm) {
    const auto a = mode(values); const double zeta = std::sqrt(l * (l + 1.0));
    return {Complex{}, Complex(0,1)*zeta/2.0*a.W*(xp+xm), zeta/2.0*a.W*(xp-xm)};
  }

  static SpinStrainField radialToStrain(int l, const Eigen::MatrixXcd &values,
      double r, Complex x0, Complex xp, Complex xm, Complex xpp, Complex xmm) {
    const auto a = mode(values); const double lambda = l * (l + 1.0);
    const double zeta = std::sqrt(lambda);
    const double eta = l >= 2 ? std::sqrt((l - 1.0) * (l + 2.0)) : 0.0;
    return {a.Up*x0, (2.0*a.U-lambda*a.V)/r*x0,
      -zeta/2.0*(a.Vp-(a.V-a.U)/r+Complex(0,1)*(a.Wp-a.W/r))*xp,
      zeta/2.0*(a.Vp-(a.V-a.U)/r-Complex(0,1)*(a.Wp-a.W/r))*xm,
      zeta*eta/r*(a.V+Complex(0,1)*a.W)*xpp,
      zeta*eta/r*(a.V-Complex(0,1)*a.W)*xmm};
  }

  static void solveFrequencyModes(Full1D::SEM &sem,
      SourceInfo::EarthquakeCMT &source,InputParameters &parameters,Complex z,
      const Eigen::VectorXcd &adjointWeights,const ModalVisitor &visitor,
      std::vector<PoleModalBlock> *poleBlocks=nullptr,
      PoleFrame *sourceFrame=nullptr) {
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
      const auto forward=FixedMeshAdjointSolver::solveForwardFrequency(
          sem,SensitivitySolveFamily::Radial,0,z,treatment.force);
      Eigen::MatrixXcd receiver=Eigen::MatrixXcd::Zero(
          3*parameters.num_receivers(),treatment.force.rows());
      for(int rec=0;rec<parameters.num_receivers();++rec)
        receiver.row(3*rec)=sem.rvZR(parameters,rec).transpose();
      const Eigen::MatrixXcd rhs=receiver.adjoint()*(std::conj(fz)*adjointWeights);
      const auto adjoint=FixedMeshAdjointSolver::solveAdjointFrequency(
          sem,SensitivitySolveFamily::Radial,0,z,rhs);
      const Eigen::MatrixXcd noData;
      visitor({SensitivitySolveFamily::Radial,0,treatment.sourceElement,
          forward,adjoint,treatment.lift,noData,noData});
    }
    for(int degree=std::max(1,parameters.lmin());degree<=lmax;++degree){
      const int orders=2*degree+1;
      if(flags.inc_sph){
        auto treatment=sem.sourceTreatedSpheroidal(source,degree,z);
        const auto forward=FixedMeshAdjointSolver::solveForwardFrequency(
            sem,SensitivitySolveFamily::Spheroidal,degree,z,treatment.force);
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
        Eigen::MatrixXcd adjoint,adjointUnit;
        FixedMeshAdjointSolver::solveAdjointFrequencyWithSupplementalBasis(
            sem,SensitivitySolveFamily::Spheroidal,degree,z,rhs,
            receiver.topRows(3).adjoint(),adjoint,adjointUnit);
        visitor({SensitivitySolveFamily::Spheroidal,degree,treatment.sourceElement,
          forward,adjoint,treatment.lift,sourceCoefficients,modeReceiver});
        if(poleBlocks){
          const auto reduced=sem.calculateForceRedCoefficients(source,degree,0.0);
          poleBlocks->push_back({SensitivitySolveFamily::Spheroidal,degree,
              treatment.sourceElement,forward,adjointUnit,treatment.lift,reduced});
        }
      }
      if(flags.inc_tor&&!sem.elementIsFluid(sem.sourceElement(source))){
        auto treatment=sem.sourceTreatedToroidal(source,degree,z);
        const auto forward=FixedMeshAdjointSolver::solveForwardFrequency(
            sem,SensitivitySolveFamily::Toroidal,degree,z,treatment.force);
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
        Eigen::MatrixXcd adjoint,adjointUnit;
        FixedMeshAdjointSolver::solveAdjointFrequencyWithSupplementalBasis(
            sem,SensitivitySolveFamily::Toroidal,degree,z,rhs,
            receiver.topRows(3).adjoint(),adjoint,adjointUnit);
        visitor({SensitivitySolveFamily::Toroidal,degree,treatment.sourceElement,
          forward,adjoint,treatment.lift,sourceCoefficients,modeReceiver});
        if(poleBlocks){
          const auto reduced=sem.calculateForceRedCoefficientsT(source,degree,0.0);
          poleBlocks->push_back({SensitivitySolveFamily::Toroidal,degree,
              treatment.sourceElement,forward,adjointUnit,treatment.lift,reduced});
        }
      }
    }
  }

  static void accumulateModalFields(const Full1D::SEM &sem,
      const ModalFrequencyBlock3D &block,const std::vector<Point> &points,
      const std::vector<double> &radii,const std::vector<int> &elements,
      const std::vector<std::vector<std::size_t>> &radialGroups,
      double radiusScale,std::vector<Sample> &result,const AngularWigner *w,
      std::size_t angleBegin,const std::vector<std::size_t> &angleForPoint) {
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
          add(result[q].forwardRawDisplacement,radialToDisplacement(0,f,1.0,y00,{},{},{},{}));
          add(result[q].adjointDualDisplacement,radialToDisplacement(0,a,1.0,y00,{},{},{},{}));
          if(radii[q]>0.0){
            add(result[q].forwardRawStrain,radialToStrain(0,f,radii[q],y00,{},{},{},{}));
            add(result[q].adjointDualStrain,radialToStrain(0,a,radii[q],y00,{},{},{},{}));}
          add(result[q].rScaledForwardStrainMetres,scale(rScaledStrainFromModal(
              0,f(0,0),f(1,0),{},{},{},{},radii[q],y00,{},{},{},{}),radiusScale));
          add(result[q].rScaledAdjointDualStrainMetres,scale(rScaledStrainFromModal(
              0,a(0,0),a(1,0),{},{},{},{},radii[q],y00,{},{},{},{}),radiusScale));
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
      std::vector<Sample> &result) {
    if(blocks.empty())return;
    struct PointAngle { double theta=0.0,phi=0.0,a=1.0,b=0.0; };
    std::vector<PointAngle> angles(points.size());
    std::map<std::uint64_t,std::size_t> angleMap;
    std::vector<double> uniqueTheta;
    std::vector<std::size_t> angleForPoint(points.size());
    for(std::size_t p=0;p<points.size();++p){
      double a=1.0,b=0.0;
      const auto local=tangentRotation(frame,points[p].thetaRadians,points[p].phiRadians,a,b);
      angles[p]={local[0],local[1],a,b};
      const auto [it,added]=angleMap.try_emplace(std::bit_cast<std::uint64_t>(local[0]),uniqueTheta.size());
      if(added)uniqueTheta.push_back(local[0]);
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
    for(std::size_t chunk=0;chunk<chunkCount;++chunk){
      const std::size_t begin=chunk*angleChunkSize;
      AngularWigner w(maxDegree,std::min(2,maxDegree),std::min(2,maxDegree),chunks[chunk].theta);
      for(const auto &block:blocks){
        const int l=block.degree,mmax=(static_cast<int>(block.fields.cols())-1)/2;
        for(const auto &group:chunks[chunk].radialGroups){
          const auto representative=group.front();
          if(block.family==SensitivitySolveFamily::Toroidal&&
             (sem.elementIsFluid(elements[representative])||elements[representative]<sem.el()||elements[representative]>=sem.eu()))continue;
          for(int m=-mmax;m<=mmax;++m){
            const int column=m+mmax;
            const auto values=interpolate(sem,block.family,block.fields,block.lift,
                radii[representative],elements[representative],column,block.sourceElement);
            for(const auto p:group){
              const std::size_t thetaIndex=angleForPoint[p]-begin;
              const Complex phase=std::exp(Complex(0.0,m*angles[p].phi));
              const Complex x0=block.family==SensitivitySolveFamily::Spheroidal?
                  harmonic(w,l,m,0,thetaIndex,phase):Complex{};
              const Complex xp=harmonic(w,l,m,1,thetaIndex,phase);
              const Complex xm=harmonic(w,l,m,-1,thetaIndex,phase);
              const Complex xpp=l>=2?harmonic(w,l,m,2,thetaIndex,phase):Complex{};
              const Complex xmm=l>=2?harmonic(w,l,m,-2,thetaIndex,phase):Complex{};
              Sample contribution;contribution.radialElement=elements[p];
              contribution.fluid=sem.elementIsFluid(elements[p]);
              const auto displacement=block.family==SensitivitySolveFamily::Spheroidal?
                  radialToDisplacement(l,values,radii[p],x0,xp,xm,xpp,xmm):
                  toroidalDisplacement(l,values,xp,xm);
              const auto strain=radii[p]>0.0?
                  radialToStrain(l,values,radii[p],x0,xp,xm,xpp,xmm):SpinStrainField{};
              const auto modal=mode(values);
              const auto rScaled=scale(rScaledStrainFromModal(
                  l,modal.U,modal.Up,modal.V,modal.Vp,modal.W,modal.Wp,radii[p],
                  x0,xp,xm,xpp,xmm),radiusScale);
              if(block.adjoint){
                contribution.adjointDualDisplacement=displacement;
                contribution.adjointDualStrain=strain;
                contribution.rScaledAdjointDualStrainMetres=rScaled;
              }else{
                contribution.forwardRawDisplacement=displacement;
                contribution.forwardRawStrain=strain;
                contribution.rScaledForwardStrainMetres=rScaled;
              }
              rotateContribution(contribution,angles[p].a,angles[p].b,block.adjoint);
              if(block.adjoint){
                add(result[p].adjointDualDisplacement,contribution.adjointDualDisplacement);
                add(result[p].adjointDualStrain,contribution.adjointDualStrain);
                add(result[p].rScaledAdjointDualStrainMetres,contribution.rScaledAdjointDualStrainMetres);
              }else{
                add(result[p].forwardRawDisplacement,contribution.forwardRawDisplacement);
                add(result[p].forwardRawStrain,contribution.forwardRawStrain);
                add(result[p].rScaledForwardStrainMetres,contribution.rScaledForwardStrainMetres);
              }
            }
          }
        }
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
  static Windowed3DFrequencyWorkspace prepare(Full1D::SEM &sem,
      SourceInfo::EarthquakeCMT &source,InputParameters &parameters,Complex z,
      const Eigen::VectorXcd &adjointWeights) {
    sensitivity_detail::requireConstructedSEM(sem);
    std::vector<OwnedModalFrequencyBlock3D> blocks;
    std::vector<windowed3d_detail::PoleModalFrequencyBlock3D> poleBlocks;
    windowed3d_detail::PoleFrame3D sourceFrame;
    const auto capture=[&](const ModalFrequencyBlock3D &b){
      blocks.push_back({b.family,b.degree,b.sourceElement,b.forward,b.adjoint,
                        b.sourceLift,b.sourceCoefficients,b.modeReceiver});
    };
    Windowed3DFieldEngine::solveFrequencyModes(
        sem,source,parameters,z,adjointWeights,capture,&poleBlocks,&sourceFrame);
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
       blocks_(std::move(blocks)),poleBlocks_(std::move(poleBlocks)),sourceFrame_(sourceFrame){}

  std::vector<Field3DSample> samplePoleFrames(
      const std::vector<Field3DPoint> &points,bool allowCenter) const {
    std::vector<OwnedModalFrequencyBlock3D> radialBlocks;
    for(const auto &block:blocks_)if(block.family==SensitivitySolveFamily::Radial)
      radialBlocks.push_back(block);
    auto result=Windowed3DFieldEngine::sampleModalBlocks(
        *sem_,*parameters_,z_,radialBlocks,points,allowCenter);
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
    std::vector<windowed3d_detail::PoleSynthesisBlock3D> sourceSynthesis;
    const Complex fz=-Complex(0.0,1.0)/z_;
    for(const auto &block:poleBlocks_){
      const auto coeff=block.sourceCoefficients.transpose();
      sourceSynthesis.push_back({block.family,block.degree,block.sourceElement,false,
          block.forward*coeff,block.sourceLift*coeff});
    }
    for(int rec=0;rec<parameters_->num_receivers();++rec){
      std::vector<windowed3d_detail::PoleSynthesisBlock3D> adjointSynthesis;
      const auto station=parameters_->receivers()[rec];
      const auto receiverFrame=Windowed3DFieldEngine::poleFrame(station.first,station.second);
      for(const auto &block:poleBlocks_){
        const auto poleReceiver=Windowed3DFieldEngine::receiverAtPole(block.family,block.degree);
        const int mmax=(static_cast<int>(poleReceiver.cols())-1)/2;
        Eigen::MatrixXcd localWeights=Eigen::MatrixXcd::Zero(3,2*mmax+1);
        for(int comp=0;comp<3;++comp)for(int m=-mmax;m<=mmax;++m)
          localWeights(comp,m+mmax)=std::conj(poleReceiver(comp,m+mmax))*
              (std::conj(fz)*adjointWeights_[3*rec+comp]);
        const auto localAdjoint=block.adjointUnit*localWeights;
        adjointSynthesis.push_back({block.family,block.degree,block.sourceElement,true,
            localAdjoint,Eigen::MatrixXcd{}});
      }
      Windowed3DFieldEngine::synthesizePoleFrame(*sem_,adjointSynthesis,
          receiverFrame,points,radii,elements,scale,result);
    }
    Windowed3DFieldEngine::synthesizePoleFrame(*sem_,sourceSynthesis,
        sourceFrame_,points,radii,elements,scale,result);
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
};

inline Windowed3DFrequencyWorkspace Windowed3DFieldEngine::prepareFrequencyWorkspace(
    Full1D::SEM &sem,SourceInfo::EarthquakeCMT &source,InputParameters &parameters,
    Complex z,const Eigen::VectorXcd &adjointWeights) {
  return Windowed3DFrequencyWorkspace::prepare(sem,source,parameters,z,adjointWeights);
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
