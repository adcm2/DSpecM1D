#include "sem_cases.hpp"
#include "profile_headers/SensitivityKernels.h"
#include <DSpecM1D/All>
#include <DSpecM1D/ModelInput>
#include <DSpecM1D/src/InputParser.h>
#include <DSpecM1D/src/SourceInfo.h>
#include <DSpecM1D/src/ParamInfo.h>

#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <sys/resource.h>

namespace fs=std::filesystem;
using namespace dssolver;
using Clock=std::chrono::steady_clock;
constexpr int nt=128, nfreq=nt/2+1;

std::string paramsText(const fs::path &model,const fs::path &output,int lmax,int rawBins) {
  const double f2mHz=std::max(0.8,rawBins*0.390625+0.4);
  return "\""+(output/"unused").string()+"\"\n\""+model.string()+
      "\"\n4\n0\n2\n0\n0\n1e-5\n0\n"+std::to_string(lmax)+
      "\n0.1\n"+std::to_string(f2mHz)+"\n30\n20\n0\n0\n0\n0\n50\n10\n20\n1e19\n2e18\n-3e18\n-5e18\n4e18\n-5e18\n3\n1\n45 90\n";
}
Eigen::VectorXd windowHann() {
  Eigen::VectorXd w=Eigen::VectorXd::Zero(nt);
  for(int n=0;n<=90;++n)w[n]=0.5-0.5*std::cos(2.0*std::numbers::pi*n/90.0);
  return w;
}
long peakRssBytes() {
  rusage usage{};if(getrusage(RUSAGE_SELF,&usage)!=0)throw std::runtime_error("getrusage failed");
#if defined(__APPLE__)
  return usage.ru_maxrss;
#else
  return usage.ru_maxrss*1024L;
#endif
}
const char *familyName(Family f) {
  switch(f){case Family::Radial:return "radial";case Family::Spheroidal:return "spheroidal";case Family::Toroidal:return "toroidal";}
  return "unknown";
}

Dense makeActualAdjointRhs(Full1D::SEM &sem,InputParameters &parameters,
    ParamInfo &angular,Family family,int degree,Complex z,const Eigen::VectorXcd &qhat) {
  const Complex fz=-Complex(0.0,1.0)/z;
  if(family==Family::Radial){
    Dense receiver=Dense::Zero(3*parameters.num_receivers(),sem.hR().rows());
    for(int rec=0;rec<parameters.num_receivers();++rec)
      receiver.row(3*rec)=sem.rvZR(parameters,rec).col(0).transpose();
    return receiver.adjoint()*(std::conj(fz)*qhat);
  }
  if(family==Family::Spheroidal){
    const auto elems=sem.receiverElements(parameters,false);
    const auto base=sem.rvBaseFull(parameters,degree,false);
    const int low=static_cast<int>(sem.ltgS(0,elems.front(),0));
    Dense receiver=Dense::Zero(base.rows(),sem.hS(degree).rows());
    receiver.block(0,low,base.rows(),base.cols())=base;
    const Dense modes=angular.rvFullSph(degree);
    Dense rhs=Dense::Zero(receiver.cols(),2*degree+1);
    for(int m=-degree;m<=degree;++m){
      Dense weighted=receiver;
      for(Eigen::Index row=0;row<weighted.rows();++row)weighted.row(row)*=modes(row,m+degree);
      rhs.col(m+degree)=weighted.adjoint()*(std::conj(fz)*qhat);
    }
    return rhs;
  }
  const auto elems=sem.receiverElements(parameters,true);
  const auto base=sem.rvBaseFullT(parameters,degree);
  const int low=static_cast<int>(sem.ltgT(elems.front(),0));
  Dense receiver=Dense::Zero(base.rows(),sem.hTk(degree).rows());
  receiver.block(0,low,base.rows(),base.cols())=base;
  const Dense modes=angular.rvFullTor(degree);
  Dense rhs=Dense::Zero(receiver.cols(),2*degree+1);
  for(int m=-degree;m<=degree;++m){
    Dense weighted=receiver;
    for(Eigen::Index row=0;row<weighted.rows();++row)weighted.row(row)*=modes(row,m+degree);
    rhs.col(m+degree)=weighted.adjoint()*(std::conj(fz)*qhat);
  }
  return rhs;
}

void quote(const std::string &s) {
  std::cout<<'"';for(char c:s){if(c=='"'||c=='\\')std::cout<<'\\';if(c=='\n'||c=='\r')std::cout<<' ';else std::cout<<c;}std::cout<<'"';
}
void emitFailure(const char *backend,const ActualOperatorCase &c,const std::string &message) {
  std::cout<<R"({"backend":)";quote(backend);std::cout<<R"(,"ok":false,"family":)";quote(familyName(c.family));
  std::cout<<R"(,"degree":)"<<c.degree<<R"(,"raw_bin":)"<<c.rawBin<<R"(,"z_re":)"<<c.z.real()<<R"(,"z_im":)"<<c.z.imag()<<R"(,"error":)";quote(message);std::cout<<"}\n";
}

void emitComparison(const char *name,const ActualOperatorCase &c,const Dense &adjRhs,
    const Dense &forward,const Dense &adjoint,const Dense &refForward,const Dense &refAdjoint,
    double factorSeconds,double forwardSeconds,double adjointSeconds,
    double scalingSeconds,double packSeconds,double analysisSeconds,double numericSeconds,
    double substitutionSeconds,double recoverySeconds,double schurAssemblySeconds,
    int traceSize=-1,int interiorSize=-1,bool symbolicReused=false) {
  const auto bw=measureBandwidth(c.matrix);
  std::cout<<std::setprecision(17)<<R"({"backend":)";quote(name);
  std::cout<<R"(,"ok":true,"family":)";quote(familyName(c.family));
  std::cout<<R"(,"degree":)"<<c.degree<<R"(,"raw_bin":)"<<c.rawBin<<R"(,"z_re":)"<<c.z.real()<<R"(,"z_im":)"<<c.z.imag()
           <<R"(,"n":)"<<c.matrix.rows()<<R"(,"nnz":)"<<c.matrix.nonZeros()
           <<R"(,"kl":)"<<bw.lower<<R"(,"ku":)"<<bw.upper
           <<R"(,"trace":)"<<traceSize<<R"(,"interior":)"<<interiorSize
           <<R"(,"forward_residual":)"<<relativeResidual(c.matrix,forward,c.sourceRhs)
           <<R"(,"adjoint_residual":)"<<relativeResidual(c.matrix.adjoint(),adjoint,adjRhs)
           <<R"(,"forward_field_relerr":)"<<relativeFieldDifference(forward,refForward)
           <<R"(,"adjoint_field_relerr":)"<<relativeFieldDifference(adjoint,refAdjoint)
           <<R"(,"factor_seconds":)"<<factorSeconds<<R"(,"forward_solve_seconds":)"<<forwardSeconds
           <<R"(,"adjoint_solve_seconds":)"<<adjointSeconds<<R"(,"scaling_seconds":)"<<scalingSeconds
           <<R"(,"pack_seconds":)"<<packSeconds<<R"(,"analysis_seconds":)"<<analysisSeconds
           <<R"(,"numeric_factor_seconds":)"<<numericSeconds
           <<R"(,"substitution_kernel_seconds":)"<<substitutionSeconds
           <<R"(,"recovery_seconds":)"<<recoverySeconds
           <<R"(,"schur_assembly_seconds":)"<<schurAssemblySeconds
           <<R"(,"symbolic_reused":)"<<(symbolicReused?"true":"false")
           <<R"(,"peak_rss_bytes":)"<<peakRssBytes()<<"}\n";
}

template<class Backend>
bool runBackend(const char *name,Backend &solver,const ActualOperatorCase &c,
    const Dense &adjRhs,const Dense &refForward,const Dense &refAdjoint,
    double scale,double pack,double analysis,double numeric,int trace=-1,int interior=-1,
    bool useElementPartition=false) {
  try{
    const auto t0=Clock::now();
    if constexpr(requires { solver.factorize(c.matrix,c.elements); }) {
      if(useElementPartition)solver.factorize(c.matrix,c.elements);
      else solver.compute(c.matrix);
    } else if constexpr(requires { solver.factorize(c.matrix); }) solver.factorize(c.matrix);
    else solver.compute(c.matrix);
    const double factor=std::chrono::duration<double>(Clock::now()-t0).count();
    auto t1=Clock::now();const Dense f=solver.solve(c.sourceRhs);
    const double ft=std::chrono::duration<double>(Clock::now()-t1).count();
    t1=Clock::now();const Dense a=solver.solveAdjoint(adjRhs);
    const double at=std::chrono::duration<double>(Clock::now()-t1).count();
    const auto stage=solver.stats();
    double scaled=scale,packed=pack,analyzed=analysis,numericFactor=numeric;
    if constexpr(requires { stage.scalingSeconds; })scaled=stage.scalingSeconds;
    if constexpr(requires { stage.packSeconds; })packed=stage.packSeconds;
    if constexpr(requires { stage.analysisSeconds; })analyzed=stage.analysisSeconds;
    if constexpr(requires { stage.factorSeconds; })numericFactor=stage.factorSeconds;
    double subs=0.0,recovery=0.0,schurAssembly=0.0;bool symbolicReused=false;
    if constexpr(requires { stage.substitutionSeconds; })subs=stage.substitutionSeconds;
    if constexpr(requires { stage.recoverySeconds; })recovery=stage.recoverySeconds;
    if constexpr(requires { stage.schurAssemblySeconds; })schurAssembly=stage.schurAssemblySeconds;
    if constexpr(requires { stage.symbolicReused; })symbolicReused=stage.symbolicReused;
    emitComparison(name,c,adjRhs,f,a,refForward,refAdjoint,factor,ft,at,
        scaled,packed,analyzed,numericFactor,subs,recovery,schurAssembly,trace,interior,symbolicReused);
    return true;
  }catch(const std::exception &e){emitFailure(name,c,e.what());return false;}
}

int main(int argc,char **argv) try {
  if(argc!=5)throw std::runtime_error("usage: sem_solver_benchmark MODEL OUTPUT LMAX RAW_BINS");
  const fs::path model=fs::absolute(argv[1]),output=fs::absolute(argv[2]);
  const int lmax=std::stoi(argv[3]),rawBins=std::stoi(argv[4]);
  if((lmax!=3&&lmax!=6&&lmax!=12)||(rawBins!=1&&rawBins!=3&&rawBins!=12&&rawBins!=48))
    throw std::runtime_error("LMAX must be 3/6/12 and RAW_BINS 1/3/12/48");
  if(!fs::is_directory(output)||fs::directory_iterator(output)!=fs::directory_iterator{})
    throw std::runtime_error("OUTPUT must be an existing empty directory");
  const auto parameterFile=output/"parameters.txt";
  {std::ofstream f(parameterFile);if(!f)throw std::runtime_error("cannot write temporary parameters");f<<paramsText(model,output,lmax,rawBins);}
  InputParametersNew input(parameterFile.string(),8,1,0.05,1.0,0.1,0.0,1);
  auto &parameters=input.inputParameters();auto &frequency=input.freqFull();auto &source=input.cmt();auto &base=input.earthModel();
  Full1D::SEM sem(base,0.05,8,lmax,source);
  if(frequency.nt()!=nt||frequency.i1()>1||frequency.i2()<rawBins+1)
    throw std::runtime_error("resolved frequency grid does not contain requested raw bins");
  Eigen::MatrixXcd weights=Eigen::MatrixXcd::Zero(3*parameters.num_receivers(),nfreq);
  weights(0,1)=Complex(0.7,0.2);weights(0,2)=Complex(-0.25,0.45);
  SPARSESPEC::WindowedDisplacementMap measurement(frequency,sem.lengthScaleMetres(),windowHann());
  const Eigen::MatrixXcd qhat=measurement.adjoint(weights);
  ParamInfo angular(parameters,lmax);
  const auto omega=frequency.w();
  const double ieps=-frequency.ep();
  std::unordered_map<int,std::unique_ptr<EquilibratedNaturalSparseLU>> naturalByFamilyDegree;
  std::unordered_map<int,std::unique_ptr<ExactStaticCondensation<>>> condensedByFamilyDegree;
  const auto run=[&](Family family,int degree,int bin){
    const Complex z=omega[bin]+Complex(0.0,ieps);
    const auto caseStart=Clock::now();
    ActualOperatorCase c=makeActualOperatorCase(sem,source,family,degree,z);c.rawBin=bin;
    const Dense adjRhs=makeActualAdjointRhs(sem,parameters,angular,family,degree,z,qhat.col(bin));
    const double setup=std::chrono::duration<double>(Clock::now()-caseStart).count();
    if(c.sourceRhs.rows()!=c.matrix.rows()||adjRhs.rows()!=c.matrix.rows()||adjRhs.cols()==0)
      throw std::runtime_error("actual source/receiver RHS dimensions differ from H");
    SPARSESPEC::EquilibratedSparseLU<Complex> reference;
    const auto t0=Clock::now();reference.compute(c.matrix);
    const double refFactor=std::chrono::duration<double>(Clock::now()-t0).count();
    auto t1=Clock::now();const Dense refF=reference.solve(c.sourceRhs);
    const double refFs=std::chrono::duration<double>(Clock::now()-t1).count();
    t1=Clock::now();const Dense refA=reference.solveAdjoint(adjRhs);
    const double refAs=std::chrono::duration<double>(Clock::now()-t1).count();
    emitComparison("colamd_reference",c,adjRhs,refF,refA,refF,refA,refFactor,refFs,refAs,0,0,0,0,refFs+refAs,0,0);

    const int backendKey=static_cast<int>(family)*100+degree;
    auto &natural=naturalByFamilyDegree[backendKey];
    if(!natural)natural=std::make_unique<EquilibratedNaturalSparseLU>();
    runBackend("natural_ordering",*natural,c,adjRhs,refF,refA,0,0,0,0);
    EquilibratedGeneralBandLU band;runBackend("equilibrated_band_lu",band,c,adjRhs,refF,refA,0,0,0,0);
    auto &condensed=condensedByFamilyDegree[backendKey];
    if(!condensed)condensed=std::make_unique<ExactStaticCondensation<>>();
    try{
      const auto t=Clock::now();condensed->factorize(c.matrix,c.elements);
      const double factor=std::chrono::duration<double>(Clock::now()-t).count();
      auto ts=Clock::now();const Dense cf=condensed->solve(c.sourceRhs);
      const double cfs=std::chrono::duration<double>(Clock::now()-ts).count();
      ts=Clock::now();const Dense ca=condensed->solveAdjoint(adjRhs);
      const double cas=std::chrono::duration<double>(Clock::now()-ts).count();
      emitComparison("element_static_condensation",c,adjRhs,cf,ca,refF,refA,factor,cfs,cas,
          condensed->stats().scalingSeconds,0.0,
          condensed->stats().analysisSeconds,condensed->stats().schurFactorSeconds,
          condensed->stats().substitutionSeconds,condensed->stats().recoverySeconds,
          condensed->stats().schurAssemblySeconds,
          condensed->stats().traceSize,condensed->stats().interiorSize,
          condensed->stats().symbolicReused);
    }catch(const std::exception &e){emitFailure("element_static_condensation",c,e.what());}
    std::cerr<<"case_setup_seconds="<<std::setprecision(9)<<setup<<" family="<<familyName(family)
             <<" degree="<<degree<<" bin="<<bin<<" ref_factor_seconds="<<refFactor<<'\n';
  };
  for(int bin=1;bin<=rawBins;++bin){
    run(Family::Radial,0,bin);
    if(!sem.elementIsFluid(sem.sourceElement(source)))
      for(int degree=1;degree<=lmax;++degree){
        run(Family::Spheroidal,degree,bin);
        run(Family::Toroidal,degree,bin);
      }
    else for(int degree=1;degree<=lmax;++degree)run(Family::Spheroidal,degree,bin);
  }
  fs::remove(parameterFile);
  return 0;
}catch(const std::exception &e){std::cerr<<"SEM solver benchmark: "<<e.what()<<'\n';return 1;}
