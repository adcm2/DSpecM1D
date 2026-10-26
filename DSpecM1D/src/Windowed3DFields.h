#ifndef DSPECM1D_WINDOWED_3D_FIELDS_H
#define DSPECM1D_WINDOWED_3D_FIELDS_H

#include "SensitivityKernels.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
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
    }
    for(const auto &owned:blocks){
      const ModalFrequencyBlock3D block{owned.family,owned.degree,owned.sourceElement,
        owned.forward,owned.adjoint,owned.sourceLift,owned.sourceCoefficients,owned.modeReceiver};
      accumulateModalFields(sem,block,points,radii,elements,scale,result);
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
  struct RadialMode { Complex U{}, Up{}, V{}, Vp{}, W{}, Wp{}; };

  static Complex harmonic(int degree, int m, int n, double theta, double phi) {
    GSHTrans::Wigner<double, GSHTrans::All, GSHTrans::All,
                     GSHTrans::Single> w(degree, degree, std::min(2, degree), theta);
    return w[n][degree, m] * std::exp(Complex(0.0, m * phi));
  }

  template <class Wigner>
  static Complex harmonic(const Wigner &w, int degree, int m, int n,
                          const Point &point, std::size_t index) {
    return w[n, static_cast<int>(index)][degree, m] *
           std::exp(Complex(0.0, m * point.phiRadians));
  }

  static Eigen::MatrixXcd embedReceiver(const Full1D::SEM &sem,
      const Eigen::MatrixXcd &local, std::size_t low, int dofs) {
    if (low + local.cols() > dofs)
      throw std::invalid_argument("3D field receiver slice exceeds the SEM operator.");
    Eigen::MatrixXcd full = Eigen::MatrixXcd::Zero(local.rows(), dofs);
    full.block(0, low, local.rows(), local.cols()) = local;
    return full;
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
      const Eigen::VectorXcd &adjointWeights,const ModalVisitor &visitor) {
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
        const auto receiver=embedReceiver(sem,sem.rvBaseFull(parameters,degree,false),
            sem.ltgS(0,receivers.front(),0),static_cast<int>(sem.hS(degree).rows()));
        Eigen::MatrixXcd rhs=Eigen::MatrixXcd::Zero(receiver.cols(),orders);
        for(int m=-degree;m<=degree;++m){
          Eigen::MatrixXcd weighted=receiver;
          for(int row=0;row<weighted.rows();++row)weighted.row(row)*=modeReceiver(row,m+degree);
          rhs.col(m+degree)=weighted.adjoint()*(std::conj(fz)*adjointWeights);
        }
        const auto adjoint=FixedMeshAdjointSolver::solveAdjointFrequency(
            sem,SensitivitySolveFamily::Spheroidal,degree,z,rhs);
        visitor({SensitivitySolveFamily::Spheroidal,degree,treatment.sourceElement,
          forward,adjoint,treatment.lift,sourceCoefficients,modeReceiver});
      }
      if(flags.inc_tor&&!sem.elementIsFluid(sem.sourceElement(source))){
        auto treatment=sem.sourceTreatedToroidal(source,degree,z);
        const auto forward=FixedMeshAdjointSolver::solveForwardFrequency(
            sem,SensitivitySolveFamily::Toroidal,degree,z,treatment.force);
        const Eigen::MatrixXcd sourceCoefficients=sem.calculateForceCoefficientsT(source,degree);
        const Eigen::MatrixXcd modeReceiver=angular.rvFullTor(degree);
        const auto receivers=sem.receiverElements(parameters,true);
        const auto receiver=embedReceiver(sem,sem.rvBaseFullT(parameters,degree),
            sem.ltgT(receivers.front(),0),static_cast<int>(sem.hTk(degree).rows()));
        Eigen::MatrixXcd rhs=Eigen::MatrixXcd::Zero(receiver.cols(),orders);
        for(int m=-degree;m<=degree;++m){
          Eigen::MatrixXcd weighted=receiver;
          for(int row=0;row<weighted.rows();++row)weighted.row(row)*=modeReceiver(row,m+degree);
          rhs.col(m+degree)=weighted.adjoint()*(std::conj(fz)*adjointWeights);
        }
        const auto adjoint=FixedMeshAdjointSolver::solveAdjointFrequency(
            sem,SensitivitySolveFamily::Toroidal,degree,z,rhs);
        visitor({SensitivitySolveFamily::Toroidal,degree,treatment.sourceElement,
          forward,adjoint,treatment.lift,sourceCoefficients,modeReceiver});
      }
    }
  }

  static void accumulateModalFields(const Full1D::SEM &sem,
      const ModalFrequencyBlock3D &block,const std::vector<Point> &points,
      const std::vector<double> &radii,const std::vector<int> &elements,
      double radiusScale,std::vector<Sample> &result) {
    const int l=block.degree;
    if(block.family==SensitivitySolveFamily::Radial){
      const Complex y00=harmonic(0,0,0,0.0,0.0);
      for(std::size_t p=0;p<points.size();++p){
        const auto fRadial=interpolate(sem,block.family,block.forward,
            block.sourceLift,radii[p],elements[p],0,block.sourceElement);
        const auto aRadial=interpolate(sem,block.family,block.adjoint,
            Eigen::MatrixXcd{},radii[p],elements[p]);
        Eigen::MatrixXcd f=Eigen::MatrixXcd::Zero(6,1),a=f;
        f(0,0)=fRadial(0,0);f(1,0)=fRadial(1,0);
        a(0,0)=aRadial(0,0);a(1,0)=aRadial(1,0);
        add(result[p].forwardRawDisplacement,radialToDisplacement(0,f,1.0,y00,{},{},{},{}));
        add(result[p].adjointDualDisplacement,radialToDisplacement(0,a,1.0,y00,{},{},{},{}));
        if(radii[p]>0.0){
          add(result[p].forwardRawStrain,radialToStrain(0,f,radii[p],y00,{},{},{},{}));
          add(result[p].adjointDualStrain,radialToStrain(0,a,radii[p],y00,{},{},{},{}));}
        add(result[p].rScaledForwardStrainMetres,scale(rScaledStrainFromModal(
            0,f(0,0),f(1,0),{},{},{},{},radii[p],y00,{},{},{},{}),radiusScale));
        add(result[p].rScaledAdjointDualStrainMetres,scale(rScaledStrainFromModal(
            0,a(0,0),a(1,0),{},{},{},{},radii[p],y00,{},{},{},{}),radiusScale));
      }
      return;
    }
    std::vector<double> theta;theta.reserve(points.size());
    for(const auto &p:points)theta.push_back(p.thetaRadians);
    GSHTrans::Wigner<double,GSHTrans::All,GSHTrans::All,
                     GSHTrans::Multiple> w(l,l,std::min(2,l),theta);
    for(std::size_t p=0;p<points.size();++p){
      if(block.family==SensitivitySolveFamily::Toroidal&&
         (sem.elementIsFluid(elements[p])||elements[p]<sem.el()||elements[p]>=sem.eu()))continue;
      for(int m=-l;m<=l;++m){
        const auto fmode=interpolateCombination(sem,block.family,block.forward,
            block.sourceLift,block.sourceCoefficients.row(m+l).transpose(),
            radii[p],elements[p],block.sourceElement);
        const auto amode=interpolate(sem,block.family,block.adjoint,
            Eigen::MatrixXcd{},radii[p],elements[p],m+l);
        const Complex x0=block.family==SensitivitySolveFamily::Spheroidal?
            harmonic(w,l,m,0,points[p],p):Complex{};
        const Complex xp=harmonic(w,l,m,1,points[p],p);
        const Complex xm=harmonic(w,l,m,-1,points[p],p);
        const Complex xpp=l>=2?harmonic(w,l,m,2,points[p],p):Complex{};
        const Complex xmm=l>=2?harmonic(w,l,m,-2,points[p],p):Complex{};
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
    const auto capture=[&](const ModalFrequencyBlock3D &b){
      blocks.push_back({b.family,b.degree,b.sourceElement,b.forward,b.adjoint,
                        b.sourceLift,b.sourceCoefficients,b.modeReceiver});
    };
    Windowed3DFieldEngine::solveFrequencyModes(
        sem,source,parameters,z,adjointWeights,capture);
    return Windowed3DFrequencyWorkspace(sem,parameters,z,adjointWeights,std::move(blocks));
  }
  std::vector<Field3DSample> sample(const std::vector<Field3DPoint> &points) const {
    return Windowed3DFieldEngine::sampleModalBlocks(*sem_,*parameters_,z_,blocks_,points,false);
  }
  std::vector<Field3DSample> sampleWeightedChunk(
      const std::vector<Field3DPoint> &points) const {
    return Windowed3DFieldEngine::sampleModalBlocks(*sem_,*parameters_,z_,blocks_,points,true);
  }
  const std::vector<OwnedModalFrequencyBlock3D> &modalBlocks() const { return blocks_; }
  Full1D::SEM &sem() const { return *sem_; }
  InputParameters &parameters() const { return *parameters_; }
  const Eigen::VectorXcd &adjointWeights() const { return adjointWeights_; }
  Complex z() const { return z_; }
private:
  Windowed3DFrequencyWorkspace(Full1D::SEM &sem,InputParameters &parameters,
      Complex z,const Eigen::VectorXcd &adjointWeights,
      std::vector<OwnedModalFrequencyBlock3D> blocks)
      :sem_(&sem),parameters_(&parameters),z_(z),adjointWeights_(adjointWeights),
       blocks_(std::move(blocks)){}
  Full1D::SEM *sem_;
  InputParameters *parameters_;
  Complex z_;
  Eigen::VectorXcd adjointWeights_;
  std::vector<OwnedModalFrequencyBlock3D> blocks_;
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
