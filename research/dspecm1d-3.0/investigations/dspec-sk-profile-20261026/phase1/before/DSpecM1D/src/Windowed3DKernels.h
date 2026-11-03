#ifndef DSPECM1D_WINDOWED_3D_KERNELS_H
#define DSPECM1D_WINDOWED_3D_KERNELS_H

#include "Windowed3DFields.h"

#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

namespace SPARSESPEC {

struct ModulusKernel3D {
  double A = 0.0, C = 0.0, F = 0.0, L = 0.0, N = 0.0;
  double kappa = 0.0;
  bool fluid = false;
};

struct WeightedField3DPoint {
  Field3DPoint point;
  // Explicit SI quadrature weights. The radial weight is dr; radius squared
  // and the angular weight dOmega are applied by the integrator.
  double radialWeightMetres = 0.0;
  double angularWeightSteradians = 0.0;
};

struct IntegratedModulusKernels3D {
  double A = 0.0, C = 0.0, F = 0.0, L = 0.0, N = 0.0;
  double kappa = 0.0;
};

struct SourceShellContact3D {
  double A = 0.0, C = 0.0, F = 0.0, L = 0.0, N = 0.0;
};

struct Windowed3DObjectiveResult {
  Eigen::MatrixXcd rawDisplacementSpectrum;
  Eigen::MatrixXcd processedSISpectrum;
  double objective = 0.0;
};

struct Windowed3DFrequencyResult {
  std::vector<Field3DSample> fields;
  std::vector<SourceShellContact3D> sourceContact;
};

struct Windowed3DEvaluationResult {
  Eigen::MatrixXcd rawDisplacementSpectrum;
  Eigen::MatrixXcd processedSISpectrum;
  std::vector<ModulusKernel3D> pointKernels;
  std::vector<SourceShellContact3D> sourceContact;
  double objective = 0.0;
};

/** Opt-in local material kernels built from fully reconstructed angular fields. */
class Windowed3DKernelEngine {
public:
  using Complex = std::complex<double>;
  using Sample = Field3DSample;

  /** High-level P-map objective and frequency-summed point/contact kernels.
   * The measurement adjoint is formed once, then bins are solved sequentially. */
  static Windowed3DEvaluationResult evaluate(
      Full1D::SEM &sem, SourceInfo::EarthquakeCMT &source,
      InputParameters &parameters, const SpectraSolver::FreqFull &frequency,
      const Eigen::VectorXd &timeWindow, const Eigen::MatrixXcd &weightsSI,
      int firstBin, int endBin, const std::vector<Field3DPoint> &points,
      const std::vector<std::array<double,2>> &angles = {}) {
    const int nfreq=frequency.nt()/2+1;
    if(weightsSI.rows()!=3*parameters.num_receivers()||weightsSI.cols()!=nfreq||
       !weightsSI.allFinite()||firstBin<0||firstBin>=endBin||endBin>nfreq||points.empty())
      throw std::invalid_argument("3D evaluation dimensions or bin interval are invalid.");
    sensitivity_detail::requireConstructedSEM(sem);
    if(!sem.sourceMatches(source))
      throw std::invalid_argument("3D evaluation requires a source-aligned SEM.");
    sensitivity_detail::requireMatchingTimeScale(sem,frequency);
    std::vector<Field3DPoint> allPoints=points;
    const auto receiverElements=sem.receiverElements(parameters,false);
    if(receiverElements.empty())
      throw std::invalid_argument("3D evaluation found no radial receiver trace.");
    const int receiverElement=receiverElements.front();
    const double radians=static_cast<double>(EIGEN_PI)/180.0;
    for(int rec=0;rec<parameters.num_receivers();++rec){
      const auto coordinate=parameters.receivers().at(rec);
      const double radius=(sem.mesh().PR()-1000.0*parameters.receiver_depth()/
          sem.lengthScaleMetres())*sem.lengthScaleMetres();
      allPoints.push_back({radius,(90.0-coordinate.first)*radians,
          coordinate.second*radians,receiverElement});
    }
    WindowedDisplacementMap measurement(frequency,sem.lengthScaleMetres(),timeWindow);
    const Eigen::MatrixXcd qhat=measurement.adjoint(weightsSI);
    Windowed3DEvaluationResult result;
    result.rawDisplacementSpectrum=Eigen::MatrixXcd::Zero(3*parameters.num_receivers(),nfreq);
    result.pointKernels.resize(points.size());
    result.sourceContact.resize(angles.size());
    const auto w=frequency.w();
    for(int bin=firstBin;bin<endBin;++bin){
      const Complex z=w.at(bin)+Complex(0.0,-frequency.ep());
      const auto fields=angles.empty()
          ? Windowed3DFrequencyResult{Windowed3DFieldEngine::evaluateFrequency(
                sem,source,parameters,z,qhat.col(bin),allPoints),{}}
          : evaluateFrequencyWithContact(sem,source,parameters,z,qhat.col(bin),allPoints,angles);
      for(std::size_t p=0;p<points.size();++p){
        const auto kernel=pointKernel(sem,fields.fields[p]);
        addPoint(result.pointKernels[p],kernel);
      }
      for(std::size_t p=0;p<fields.sourceContact.size();++p){
        result.sourceContact[p].A+=fields.sourceContact[p].A;
        result.sourceContact[p].C+=fields.sourceContact[p].C;
        result.sourceContact[p].F+=fields.sourceContact[p].F;
        result.sourceContact[p].L+=fields.sourceContact[p].L;
        result.sourceContact[p].N+=fields.sourceContact[p].N;
      }
      const Complex fz=-Complex(0.0,1.0)/z;
      for(int rec=0;rec<parameters.num_receivers();++rec)
        for(int c=0;c<3;++c)
          result.rawDisplacementSpectrum(3*rec+c,bin)=
              fz*fields.fields[points.size()+static_cast<std::size_t>(rec)].forwardRawDisplacement[c];
    }
    result.processedSISpectrum=measurement.forward(result.rawDisplacementSpectrum);
    result.objective=weightsSI.conjugate().cwiseProduct(result.processedSISpectrum).sum().real();
    return result;
  }

  /** Reconstruct requested samples and contact densities from one set of
   * forward/adjoint modal solves for this frequency. */
  static Windowed3DFrequencyResult evaluateFrequencyWithContact(
      Full1D::SEM &sem, SourceInfo::EarthquakeCMT &source,
      InputParameters &parameters, Complex z,
      const Eigen::VectorXcd &adjointWeights,
      const std::vector<Field3DPoint> &points,
      const std::vector<std::array<double,2>> &angles) {
    sensitivity_detail::requireConstructedSEM(sem);
    const int e=sem.sourceElement(source);
    if(!sem.sourceMatches(source)||angles.empty()||sem.elementIsFluid(e))
      throw std::invalid_argument("3D contact requires a source-aligned SEM, solid source, and angles.");
    for(const auto &angle:angles)
      if(!std::isfinite(angle[0])||angle[0]<0.0||angle[0]>EIGEN_PI||
         !std::isfinite(angle[1]))
        throw std::invalid_argument("Source-shell angles must be finite geocentric coordinates.");
    auto workspace=Windowed3DFieldEngine::prepareFrequencyWorkspace(
        sem,source,parameters,z,adjointWeights);
    Windowed3DFrequencyResult result;
    result.fields=workspace.sample(points);
    result.sourceContact=sourceShellContact(workspace,source,angles);
    return result;
  }

  /** Contact density from an existing modal workspace (no additional solve). */
  static std::vector<SourceShellContact3D> sourceShellContact(
      const Windowed3DFrequencyWorkspace &workspace,
      SourceInfo::EarthquakeCMT &source,
      const std::vector<std::array<double,2>> &angles) {
    auto &sem=workspace.sem();auto &parameters=workspace.parameters();
    sensitivity_detail::requireConstructedSEM(sem);
    const int e=sem.sourceElement(source);
    if(!sem.sourceMatches(source)||angles.empty()||sem.elementIsFluid(e))
      throw std::invalid_argument("3D contact requires a source-aligned SEM, solid source, and angles.");
    const double C=sem.meshModel().C(e,0),L=sem.meshModel().L(e,0),P0=sem.pressureScalePascals();
    if(!(sem.sourceRadius()>0.0&&C>0.0&&L>0.0&&P0>0.0))
      throw std::invalid_argument("Source-shell radius and material/SI scales must be positive.");
    for(const auto &angle:angles)
      if(!std::isfinite(angle[0])||angle[0]<0.0||angle[0]>EIGEN_PI||!std::isfinite(angle[1]))
        throw std::invalid_argument("Source-shell angles must be finite geocentric coordinates.");
    std::vector<AngularDual> accum(angles.size());
    for(const auto &owned:workspace.modalBlocks()){
      const ModalFrequencyBlock3D block{owned.family,owned.degree,owned.sourceElement,
        owned.forward,owned.adjoint,owned.sourceLift,owned.sourceCoefficients,owned.modeReceiver};
      accumulateContactBlock(sem,source,parameters,workspace.z(),workspace.adjointWeights(),
                             angles,block,accum);
    }
    return finishContact(sem,source,accum);
  }

  /** Cached-dual angular material contact density at explicit angles. */
  static std::vector<SourceShellContact3D> sourceShellContactFrequency(
      Full1D::SEM &sem, SourceInfo::EarthquakeCMT &source,
      InputParameters &parameters, Complex z,
      const Eigen::VectorXcd &adjointWeights,
      const std::vector<std::array<double,2>> &angles) {
    auto workspace=Windowed3DFieldEngine::prepareFrequencyWorkspace(
        sem,source,parameters,z,adjointWeights);
    return sourceShellContact(workspace,source,angles);
  }

  /** Forward-only objective path for validation and streamed C evaluations. */
  static Windowed3DObjectiveResult evaluateForwardObjective(
      Full1D::SEM &sem, SourceInfo::EarthquakeCMT &source,
      InputParameters &parameters, const SpectraSolver::FreqFull &frequency,
      const Eigen::VectorXd &timeWindow, const Eigen::MatrixXcd &weightsSI,
      int firstBin, int endBin) {
    if (weightsSI.rows() != 3 * parameters.num_receivers() ||
        weightsSI.cols() != frequency.nt()/2+1 || !weightsSI.allFinite() ||
        firstBin < 0 || firstBin >= endBin || endBin > weightsSI.cols())
      throw std::invalid_argument("3D forward objective dimensions or bin interval are invalid.");
    SPARSESPEC::sensitivity_detail::requireConstructedSEM(sem);
    if (!sem.sourceMatches(source))
      throw std::invalid_argument("3D forward objective requires a source-aligned SEM.");
    SPARSESPEC::sensitivity_detail::requireMatchingTimeScale(sem,frequency);
    SPARSESPEC::WindowedDisplacementMap measurement(
        frequency, sem.lengthScaleMetres(), timeWindow);
    std::vector<Field3DPoint> points;
    points.reserve(parameters.num_receivers());
    const auto receiverElements = sem.receiverElements(parameters, false);
    if (receiverElements.empty())
      throw std::invalid_argument("3D forward objective found no radial receiver trace.");
    const int receiverElement = receiverElements.front();
    const double radians = static_cast<double>(EIGEN_PI) / 180.0;
    for (int rec = 0; rec < parameters.num_receivers(); ++rec) {
      const auto coordinate = parameters.receivers().at(rec);
      const double radius = (sem.mesh().PR() -
          1000.0 * parameters.receiver_depth() / sem.lengthScaleMetres()) *
          sem.lengthScaleMetres();
      points.push_back({radius, (90.0-coordinate.first)*radians,
                        coordinate.second*radians,
                        receiverElement});
    }
    const int nfreq = frequency.nt()/2+1;
    Windowed3DObjectiveResult result;
    result.rawDisplacementSpectrum = Eigen::MatrixXcd::Zero(3*parameters.num_receivers(), nfreq);
    Eigen::VectorXcd zeroWeights = Eigen::VectorXcd::Zero(3*parameters.num_receivers());
    const auto w = frequency.w();
    for (int bin = firstBin; bin < endBin; ++bin) {
      const Complex z = w.at(bin) + Complex(0.0, -frequency.ep());
      const Complex displacementFactor = -Complex(0.0, 1.0) / z;
      const auto fields = Windowed3DFieldEngine::evaluateFrequency(
          sem, source, parameters, z, zeroWeights, points);
      for (int rec = 0; rec < parameters.num_receivers(); ++rec)
        for (int component = 0; component < 3; ++component)
          result.rawDisplacementSpectrum(3*rec+component, bin) =
              displacementFactor * fields[rec].forwardRawDisplacement[component];
    }
    result.processedSISpectrum = measurement.forward(result.rawDisplacementSpectrum);
    result.objective = weightsSI.conjugate().cwiseProduct(
        result.processedSISpectrum).sum().real();
    return result;
  }

  static ModulusKernel3D pointKernel(const Full1D::SEM &sem,
                                     const Sample &sample) {
    const double scale = objectiveScale(sem);
    ModulusKernel3D result;
    result.fluid = sample.fluid;
    if (sample.fluid) {
      const Complex divForward = sample.forwardRawStrain.rr +
                                 sample.forwardRawStrain.horizontalTrace;
      const Complex divAdjoint = sample.adjointDualStrain.rr +
                                 sample.adjointDualStrain.horizontalTrace;
      result.kappa = -std::real(std::conj(divAdjoint) * divForward) / scale;
      return result;
    }
    contract(sample.forwardRawStrain, sample.adjointDualStrain,
             scale, result);
    return result;
  }

  /** Evaluate unweighted point samples at one frequency; no weights are inferred. */
  static std::vector<ModulusKernel3D> evaluateFrequency(
      Full1D::SEM &sem, SourceInfo::EarthquakeCMT &source,
      InputParameters &parameters, Complex z,
      const Eigen::VectorXcd &adjointWeights,
      const std::vector<Field3DPoint> &points) {
    const auto fields = Windowed3DFieldEngine::evaluateFrequency(
        sem, source, parameters, z, adjointWeights, points);
    std::vector<ModulusKernel3D> result;
    result.reserve(fields.size());
    for (const auto &field : fields)
      result.push_back(pointKernel(sem, field));
    return result;
  }

  static std::vector<ModulusKernel3D> evaluateFrequency(
      const Windowed3DFrequencyWorkspace &workspace,
      const std::vector<Field3DPoint> &points) {
    const auto fields=workspace.sample(points);
    std::vector<ModulusKernel3D> result;result.reserve(fields.size());
    for(const auto &field:fields)result.push_back(pointKernel(workspace.sem(),field));
    return result;
  }

  /**
   * Explicit weighted volume rule. Each sample contributes r^2 K dr dOmega.
   * The center is accepted only here and uses algebraic r*strain bilinears;
   * pointKernel is never called at r=0.
   */
  static IntegratedModulusKernels3D integrateWeightedVolumeFrequency(
      Full1D::SEM &sem, SourceInfo::EarthquakeCMT &source,
      InputParameters &parameters, Complex z,
      const Eigen::VectorXcd &adjointWeights,
      const std::vector<WeightedField3DPoint> &quadrature) {
    if (quadrature.empty())
      throw std::invalid_argument("A weighted 3D volume rule must be nonempty.");
    for (const auto &q : quadrature) {
      if (!(std::isfinite(q.radialWeightMetres) &&
            q.radialWeightMetres > 0.0 &&
            std::isfinite(q.angularWeightSteradians) &&
            q.angularWeightSteradians > 0.0))
        throw std::invalid_argument("3D radial and angular quadrature weights must be positive and finite.");
    }
    auto workspace=Windowed3DFieldEngine::prepareFrequencyWorkspace(
        sem,source,parameters,z,adjointWeights);
    return integrateWeightedVolumeFrequency(workspace,quadrature);
  }

  /** Stream an explicit quadrature through an existing frequency workspace. */
  static IntegratedModulusKernels3D integrateWeightedVolumeFrequency(
      const Windowed3DFrequencyWorkspace &workspace,
      const std::vector<WeightedField3DPoint> &quadrature,
      std::size_t chunkSize=256) {
    if(quadrature.empty()||chunkSize==0)
      throw std::invalid_argument("A weighted volume rule and positive chunk size are required.");
    auto &sem=workspace.sem();
    for(const auto &q:quadrature)
      if(!(std::isfinite(q.radialWeightMetres)&&q.radialWeightMetres>0.0&&
           std::isfinite(q.angularWeightSteradians)&&q.angularWeightSteradians>0.0))
        throw std::invalid_argument("3D radial and angular quadrature weights must be positive and finite.");
    const double scale = objectiveScale(sem);
    IntegratedModulusKernels3D result;
    for(std::size_t begin=0;begin<quadrature.size();begin+=chunkSize){
      const std::size_t end=std::min(quadrature.size(),begin+chunkSize);
      std::vector<Field3DPoint> points;points.reserve(end-begin);
      for(std::size_t i=begin;i<end;++i)points.push_back(quadrature[i].point);
      const auto fields=workspace.sampleWeightedChunk(points);
      for(std::size_t local=0;local<fields.size();++local){
        const std::size_t i=begin+local;
        const double measure=quadrature[i].radialWeightMetres*quadrature[i].angularWeightSteradians;
        const auto &field=fields[local];
        if(quadrature[i].point.radiusMetres==0.0){
          ModulusKernel3D scaled;scaled.fluid=field.fluid;
          if(field.fluid){
            const Complex f=field.rScaledForwardStrainMetres.rr+field.rScaledForwardStrainMetres.horizontalTrace;
            const Complex a=field.rScaledAdjointDualStrainMetres.rr+field.rScaledAdjointDualStrainMetres.horizontalTrace;
            scaled.kappa=-std::real(std::conj(a)*f)/scale;
          }else contract(field.rScaledForwardStrainMetres,
                         field.rScaledAdjointDualStrainMetres,scale,scaled);
          add(result,scaled,measure);
        }else{
          const auto kernel=pointKernel(sem,field);
          add(result,kernel,measure*quadrature[i].point.radiusMetres*quadrature[i].point.radiusMetres);
        }
      }
    }
    return result;
  }

private:
  struct AngularDual {
    Complex dC{}, dF{}, dLn{}, dLe{}, jr{}, jn{}, je{};
  };

  struct TangentialValues {
    Complex spheroidalNorth{}, spheroidalEast{};
    Complex toroidalNorth{}, toroidalEast{};
  };

  static Eigen::MatrixXcd embedReceiver(const Full1D::SEM &sem,
      const Eigen::MatrixXcd &local, std::size_t low, int dofs) {
    if (low + static_cast<std::size_t>(local.cols()) >
        static_cast<std::size_t>(dofs))
      throw std::invalid_argument("3D contact receiver slice exceeds the SEM operator.");
    Eigen::MatrixXcd full=Eigen::MatrixXcd::Zero(local.rows(),dofs);
    full.block(0,low,local.rows(),local.cols())=local;
    return full;
  }

  static Complex scalarHarmonic(int degree,int m,double theta,double phi) {
    GSHTrans::Wigner<double,GSHTrans::All,GSHTrans::All,
                     GSHTrans::Single> w(degree,degree,std::min(2,degree),theta);
    return w[0][degree,m]*std::exp(Complex(0.0,m*phi));
  }

  static TangentialValues tangentialBasis(int degree,int m,
                                           double theta,double phi) {
    const double zeta=std::sqrt(static_cast<double>(degree)*(degree+1.0));
    GSHTrans::Wigner<double,GSHTrans::All,GSHTrans::All,
                     GSHTrans::Single> w(degree,degree,std::min(2,degree),theta);
    const Complex phase=std::exp(Complex(0.0,m*phi));
    const Complex xp=w[1][degree,m]*phase, xm=w[-1][degree,m]*phase;
    const Complex i(0.0,1.0);
    return {zeta/2.0*(xp-xm), -i*zeta/2.0*(xp+xm),
            i*zeta/2.0*(xp+xm), zeta/2.0*(xp-xm)};
  }

  static void accumulateContactBlock(
      Full1D::SEM &sem, SourceInfo::EarthquakeCMT &source,
      InputParameters &parameters, Complex z,
      const Eigen::VectorXcd &adjointWeights,
      const std::vector<std::array<double,2>> &angles,
      const ModalFrequencyBlock3D &block, std::vector<AngularDual> &accum) {
    const int e=block.sourceElement, nq=sem.mesh().NN(), l=block.degree;
    const double rs=sem.sourceRadius(), F=sem.meshModel().F(e,0);
    const Complex fz=-Complex(0.0,1.0)/z;
    const int modes=2*l+1;
    std::vector<std::vector<Complex>> Y(angles.size(),std::vector<Complex>(modes));
    std::vector<std::vector<TangentialValues>> B(angles.size(),
        std::vector<TangentialValues>(modes));
    const double zeta=std::sqrt(static_cast<double>(l)*(l+1.0));
    const Complex imag(0.0,1.0);
    for(std::size_t p=0;p<angles.size();++p){
      GSHTrans::Wigner<double,GSHTrans::All,GSHTrans::All,
                       GSHTrans::Single> w(l,l,std::min(2,l),angles[p][0]);
      for(int m=-l;m<=l;++m){
        const int column=m+l;
        const Complex phase=std::exp(Complex(0.0,m*angles[p][1]));
        Y[p][column]=w[0][l,m]*phase;
        if(l>0){
          const Complex xp=w[1][l,m]*phase,xm=w[-1][l,m]*phase;
          B[p][column]={zeta/2.0*(xp-xm),-imag*zeta/2.0*(xp+xm),
                        imag*zeta/2.0*(xp+xm),zeta/2.0*(xp-xm)};
        }
      }
    }
    const auto addRadial=[&](int degree,int m,Complex D,Complex T,Complex jump){
      for(std::size_t p=0;p<angles.size();++p){
        const Complex y=Y[p][m+l];
        accum[p].dC+=(D+F*T)*std::conj(y);
        accum[p].dF+=T*std::conj(y);
        accum[p].jr+=jump*y;
      }
    };
    const auto addTangential=[&](int degree,int m,Complex D,Complex jump,bool tor){
      const double lambda=static_cast<double>(degree)*(degree+1.0);
      for(std::size_t p=0;p<angles.size();++p){
        const auto &basis=B[p][m+l];
        const Complex bn=tor?basis.toroidalNorth:basis.spheroidalNorth;
        const Complex be=tor?basis.toroidalEast:basis.spheroidalEast;
        accum[p].dLn+=(D/lambda)*std::conj(bn);
        accum[p].dLe+=(D/lambda)*std::conj(be);
        accum[p].jn+=jump*bn; accum[p].je+=jump*be;
      }
    };
    if(block.family==SensitivitySolveFamily::Radial){
      const int u=sem.ltgR(0,e,0);
      Eigen::VectorXcd unit=Eigen::VectorXcd::Zero(block.sourceLift.rows());
      unit[u]=1.0;
      Eigen::VectorXcd local=Eigen::VectorXcd::Zero(2*nq);local[0]=1.0;
      const auto op=Full1D::source_treatment_detail::radialElementOperator(sem,e,z);
      const auto action=op*local;
      Eigen::VectorXcd adjLocal(2*nq);
      for(int q=0;q<nq;++q){adjLocal[q]=block.adjoint(sem.ltgR(0,e,q),0);
        adjLocal[nq+q]=block.adjoint(sem.ltgR(1,e,q),0);}
      const auto direct=sem.sourceLiftReceiverRadial(parameters,unit,e);
      const Complex D=-adjLocal.dot(action)+adjointWeights.dot(fz*direct.col(0));
      const Complex T=-2.0*rs*std::conj(block.adjoint(u,0));
      addRadial(0,0,D,T,block.sourceLift(u,0));
    }else if(block.family==SensitivitySolveFamily::Spheroidal){
      const auto op=Full1D::source_treatment_detail::
          spheroidalElementOperator(sem,e,l,z,{},false);
      const int u=sem.ltgS(0,e,0),v=sem.ltgS(1,e,0);
      Eigen::VectorXcd unitU=Eigen::VectorXcd::Zero(block.sourceLift.rows()),unitV=unitU;
      unitU[u]=1.0;unitV[v]=1.0;
      Eigen::VectorXcd localU=Eigen::VectorXcd::Zero(3*nq),localV=localU;
      localU[0]=1.0;localV[nq]=1.0;
      const auto actionU=op*localU,actionV=op*localV;
      const auto directU=sem.sourceLiftReceiver(parameters,l,unitU,e,false);
      const auto directV=sem.sourceLiftReceiver(parameters,l,unitV,e,false);
      const double lambda=static_cast<double>(l)*(l+1.0);
      for(int m=-l;m<=l;++m){
        Eigen::VectorXcd localAdj(3*nq);
        for(int q=0;q<nq;++q){localAdj[q]=block.adjoint(sem.ltgS(0,e,q),m+l);
          localAdj[nq+q]=block.adjoint(sem.ltgS(1,e,q),m+l);
          localAdj[2*nq+q]=block.adjoint(sem.ltgS(2,e,q),m+l);}
        Eigen::VectorXcd receiverU=directU.col(0),receiverV=directV.col(0);
        for(int row=0;row<receiverU.size();++row){
          receiverU[row]*=block.modeReceiver(row,m+l);
          receiverV[row]*=block.modeReceiver(row,m+l);}
        const Complex Dr=-localAdj.dot(actionU)+adjointWeights.dot(fz*receiverU);
        const Complex Dh=-localAdj.dot(actionV)+adjointWeights.dot(fz*receiverV);
        const Complex T=-2.0*rs*std::conj(block.adjoint(u,m+l))+
            lambda*rs*std::conj(block.adjoint(v,m+l));
        addRadial(l,m,Dr,T,block.sourceLift(u,0)*block.sourceCoefficients(m+l,0));
        addTangential(l,m,Dh,block.sourceLift(v,2)*block.sourceCoefficients(m+l,2),false);
      }
    }else{
      if(e<sem.el()||e>=sem.eu())return;
      const auto op=Full1D::source_treatment_detail::toroidalElementOperator(sem,e,l,z);
      Eigen::VectorXcd unit=Eigen::VectorXcd::Zero(block.sourceLift.rows());
      unit[sem.ltgT(e,0)]=1.0;
      Eigen::VectorXcd localUnit=Eigen::VectorXcd::Zero(nq);localUnit[0]=1.0;
      const auto action=op*localUnit;
      const auto direct=sem.sourceLiftReceiverToroidal(parameters,l,unit,e);
      for(int m=-l;m<=l;++m){
        Eigen::VectorXcd localAdj(nq);
        for(int q=0;q<nq;++q)localAdj[q]=block.adjoint(sem.ltgT(e,q),m+l);
        Eigen::VectorXcd receiverMode=direct.col(0);
        for(int row=0;row<receiverMode.size();++row)
          receiverMode[row]*=block.modeReceiver(row,m+l);
        const Complex D=-localAdj.dot(action)+adjointWeights.dot(fz*receiverMode);
        addTangential(l,m,D,block.sourceLift(sem.ltgT(e,0),0)*
            block.sourceCoefficients(m+l,0),true);
      }
    }
  }

  static std::vector<SourceShellContact3D> finishContact(
      const Full1D::SEM &sem,SourceInfo::EarthquakeCMT &source,
      const std::vector<AngularDual> &accum) {
    const int e=sem.sourceElement(source);
    const double C=sem.meshModel().C(e,0),L=sem.meshModel().L(e,0),
                 P0=sem.pressureScalePascals();
    std::vector<SourceShellContact3D> result(accum.size());
    for(std::size_t p=0;p<accum.size();++p){
      result[p].C=-std::real(accum[p].dC*accum[p].jr)/(C*P0);
      result[p].F= std::real(accum[p].dF*accum[p].jr)/P0;
      result[p].L=-std::real(accum[p].dLn*accum[p].jn+
                             accum[p].dLe*accum[p].je)/(L*P0);
    }
    return result;
  }

  static double objectiveScale(const Full1D::SEM &sem) {
    const double L0 = sem.lengthScaleMetres();
    const double P0 = sem.pressureScalePascals();
    if (!(std::isfinite(L0) && L0 > 0.0 &&
          std::isfinite(P0) && P0 > 0.0))
      throw std::invalid_argument("3D kernel SI scales must be positive and finite.");
    return P0 * L0 * L0 * L0;
  }

  static void contract(const SpinStrainField &forward,
                       const SpinStrainField &adjoint,
                       double scale, ModulusKernel3D &out) {
    const Complex hF = forward.horizontalTrace, hA = adjoint.horizontalTrace;
    const Complex rF = forward.rr, rA = adjoint.rr;
    out.A = -std::real(std::conj(hA) * hF) / scale;
    out.C = -std::real(std::conj(rA) * rF) / scale;
    out.F = -std::real(std::conj(rA) * hF + std::conj(hA) * rF) / scale;
    out.L = -2.0 * std::real(
        std::conj(adjoint.radialPlus) * forward.radialPlus +
        std::conj(adjoint.radialMinus) * forward.radialMinus) / scale;
    out.N = -std::real(
        0.5 * (std::conj(adjoint.horizontalPlus) * forward.horizontalPlus +
               std::conj(adjoint.horizontalMinus) * forward.horizontalMinus) -
        std::conj(hA) * hF) / scale;
  }

  static void add(IntegratedModulusKernels3D &sum,
                  const ModulusKernel3D &term, double weight) {
    if (term.fluid) {
      sum.kappa += weight * term.kappa;
    } else {
      sum.A += weight * term.A; sum.C += weight * term.C;
      sum.F += weight * term.F; sum.L += weight * term.L;
      sum.N += weight * term.N;
    }
  }

  static void addPoint(ModulusKernel3D &sum,const ModulusKernel3D &term) {
    sum.fluid=term.fluid;
    if(term.fluid)sum.kappa+=term.kappa;
    else{sum.A+=term.A;sum.C+=term.C;sum.F+=term.F;sum.L+=term.L;sum.N+=term.N;}
  }
};

} // namespace SPARSESPEC

#endif
