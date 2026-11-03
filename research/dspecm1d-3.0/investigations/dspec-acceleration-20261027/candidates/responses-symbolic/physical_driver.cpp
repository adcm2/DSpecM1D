#include "engine_headers/SensitivityKernels.h"
#include "engine_headers/Windowed3DKernels.h"
#include "profile.hpp"
#include "global/global_kernel_consumer.hpp"
#include <DSpecM1D/All>
#include <DSpecM1D/ModelInput>
#include <DSpecM1D/src/InputParser.h>
#include <DSpecM1D/src/NormClass.h>
#include <DSpecM1D/src/SourceInfo.h>
#include <GSHTrans/GaussLegendreGrid.hpp>
#include <GSHTrans/Utility.hpp>
#include <GaussQuad/All>

#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <cstdint>
#include <limits>
#include <memory>
#include <sstream>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using Complex = std::complex<double>;
using Clock = std::chrono::steady_clock;
using namespace std::chrono_literals;

namespace {
constexpr double kPi = std::numbers::pi;
constexpr double kLengthScaleMetres = 6371000.0;
constexpr double kSourceLatitude = 10.0, kSourceLongitude = 20.0;
constexpr double kReceiverLatitude = 45.0, kReceiverLongitude = 75.0;
constexpr double kSourceDepthKm = 50.0, kReceiverDepthKm = 3.0;

struct Options {
  fs::path model, output;
  std::string mode = "forward";
  int lmax = 64, semDegree = 8, qex = 1, radialTraceLimit = 1;
  std::size_t probeAngleCount = 0;
  double dt = 5.0, duration = 3600.0, radialRefinement = 0.05, maxOutputGiB = 8.0;
  double responseCacheGiB = 32.0;
  double responseSpoolGiB = 128.0;
  std::size_t globalWorkers = 1, globalAngleTile = 128;
  bool includeDc = true, includeContacts = false, poleOnly = true, forwardOnly = true;
  bool responseReuse = true;
  bool responseRecordSpool = false;
};

Options parse(int argc, char **argv) {
  if (argc < 4) throw std::runtime_error(
      "usage: physical_driver MODEL EMPTY_OUTPUT forward|kernel|kernel-probe [--lmax N] [--nq N] [--radial-refinement X] [--dt-s S] [--duration-s S] [--qex N] [--radial-traces N|0-for-all] [--probe-angle-count N] [--dc include|exclude] [--contacts on|off] [--pole-only on|off] [--forward-only on|off] [--response-reuse on|off] [--response-cache-gib G] [--response-record-spool on|off] [--response-spool-gib G] [--global-workers N] [--global-angle-tile N] [--max-output-gib G]");
  Options o; o.model = fs::absolute(argv[1]); o.output = fs::absolute(argv[2]); o.mode = argv[3];
  if (o.mode != "forward" && o.mode != "kernel" && o.mode != "kernel-probe")
    throw std::runtime_error("mode must be forward, kernel, or kernel-probe");
  for (int i=4; i<argc; ++i) {
    if (i+1>=argc) throw std::runtime_error("option is missing its value");
    const std::string key=argv[i++], value=argv[i];
    if (key=="--lmax") o.lmax=std::stoi(value);
    else if (key=="--nq" || key=="--sem-degree") o.semDegree=std::stoi(value);
    else if (key=="--qex") o.qex=std::stoi(value);
    else if (key=="--radial-traces") o.radialTraceLimit=std::stoi(value);
    else if (key=="--probe-angle-count") o.probeAngleCount=static_cast<std::size_t>(std::stoull(value));
    else if (key=="--radial-refinement") o.radialRefinement=std::stod(value);
    else if (key=="--dt-s") o.dt=std::stod(value);
    else if (key=="--duration-s") o.duration=std::stod(value);
    else if (key=="--dc") { if(value!="include"&&value!="exclude") throw std::runtime_error("--dc must be include or exclude"); o.includeDc=value=="include"; }
    else if (key=="--contacts") { if(value!="on"&&value!="off") throw std::runtime_error("--contacts must be on or off"); o.includeContacts=value=="on"; }
    else if (key=="--pole-only") { if(value!="on"&&value!="off") throw std::runtime_error("--pole-only must be on or off"); o.poleOnly=value=="on"; }
    else if (key=="--forward-only") { if(value!="on"&&value!="off") throw std::runtime_error("--forward-only must be on or off"); o.forwardOnly=value=="on"; }
    else if (key=="--response-reuse") { if(value!="on"&&value!="off") throw std::runtime_error("--response-reuse must be on or off"); o.responseReuse=value=="on"; }
    else if (key=="--response-cache-gib") o.responseCacheGiB=std::stod(value);
    else if (key=="--response-record-spool") { if(value!="on"&&value!="off") throw std::runtime_error("--response-record-spool must be on or off"); o.responseRecordSpool=value=="on"; }
    else if (key=="--response-spool-gib") o.responseSpoolGiB=std::stod(value);
    else if (key=="--global-workers") o.globalWorkers=static_cast<std::size_t>(std::stoull(value));
    else if (key=="--global-angle-tile") o.globalAngleTile=static_cast<std::size_t>(std::stoull(value));
    else if (key=="--max-output-gib") o.maxOutputGiB=std::stod(value);
    else throw std::runtime_error("unknown option: "+key);
  }
  if (!fs::exists(o.model)) throw std::runtime_error("MODEL path does not exist");
  fs::create_directories(o.output);
  if (fs::directory_iterator(o.output)!=fs::directory_iterator{})
    throw std::runtime_error("EMPTY_OUTPUT must be empty");
  if (o.lmax<1 || o.semDegree<2 || o.qex<1 || o.radialTraceLimit<0 || !(o.dt>0.0) || !(o.duration>0.0) || !(o.radialRefinement>0.0) || !(std::isfinite(o.maxOutputGiB)&&o.maxOutputGiB>0.0) || !(std::isfinite(o.responseCacheGiB)&&o.responseCacheGiB>0.0) || !(std::isfinite(o.responseSpoolGiB)&&o.responseSpoolGiB>0.0) || o.globalWorkers==0 || o.globalWorkers>48 || o.globalAngleTile==0)
    throw std::runtime_error("degree, time, duration, and radial refinement must be positive");
  if(o.mode=="kernel-probe"&&(o.probeAngleCount==0||o.radialTraceLimit!=0))
    throw std::runtime_error("kernel-probe requires --probe-angle-count N and --radial-traces 0 (all native radial sides)");
  if(o.mode!="kernel-probe"&&o.probeAngleCount!=0)
    throw std::runtime_error("--probe-angle-count is valid only with kernel-probe mode");
  return o;
}

std::string paramsText(const Options &o) {
  // Legacy parser units: t_out is minutes and frequency cutoffs are mHz.
  std::ostringstream s;
  s << '"' << (o.output/"unused").string() << "\"\n\"" << o.model.string() << "\"\n"
    << "4\n0\n2\n0\n0\n1e-5\n0\n" << o.lmax << "\n0\n10\n"
    << o.duration/60.0 << '\n' << o.dt << "\n0\n0\n0\n0\n"
    << kSourceDepthKm << '\n' << kSourceLatitude << '\n' << kSourceLongitude << '\n'
    << "1e19\n2e18\n-3e18\n-5e18\n4e18\n-5e18\n"
    << kReceiverDepthKm << "\n1\n" << kReceiverLatitude << ' ' << kReceiverLongitude << "\n";
  return s.str();
}

Eigen::VectorXd sqrtHann(const SpectraSolver::FreqFull &freq) {
  Eigen::VectorXd w(freq.nt());
  const double dt=freq.dt()*freq.timeNorm();
  for(int n=0;n<freq.nt();++n) {
    const double t=n*dt;
    if(t<1000.0 || t>1200.0) w[n]=0.0;
    else w[n]=std::sqrt(0.5*(1.0-std::cos(2.0*kPi*(t-1000.0)/200.0)));
  }
  return w;
}

double bandWeight(double f) {
  if(f<0.008) return 1.0;
  if(f<0.010) return 0.5*(1.0+std::cos(kPi*(f-0.008)/0.002));
  return 0.0;
}

void projectRealFftEndpoints(Eigen::MatrixXcd &s) {
  if(s.cols()==0) return;
  s.col(0)=s.col(0).real().cast<Complex>();
  if((s.cols()-1)>0) s.col(s.cols()-1)=s.col(s.cols()-1).real().cast<Complex>();
}

struct PhysicalContext {
  std::unique_ptr<InputParametersNew> input;
  std::unique_ptr<Full1D::SEM> sem;
};

PhysicalContext prepare(const Options &o) {
  PhysicalContext c;
  mathprofile::Scope setup("dspec_setup");
  { mathprofile::Scope t("parameter_file_generation"); std::ofstream f(o.output/"parameters.txt"); f<<paramsText(o); }
  { mathprofile::Scope t("input_parameters_model_and_frequency_initialization");
    c.input=std::make_unique<InputParametersNew>((o.output/"parameters.txt").string(),o.semDegree,1,o.radialRefinement,1.0,0.1,0.0,o.qex); }
  auto &base=c.input->earthModel(); auto &source=c.input->cmt();
  { mathprofile::Scope t("SEM_mesh_and_matrix_setup"); c.sem=std::make_unique<Full1D::SEM>(base,o.radialRefinement,o.semDegree,o.lmax,source); }
  const auto &freq=c.input->freqFull();
  const int expected=static_cast<int>(std::ceil(o.duration/o.dt));
  if(freq.nt()<expected || std::abs(freq.dt()*freq.timeNorm()-o.dt)>1e-12)
    throw std::runtime_error("resolved DSpec FFT grid does not cover requested dt/duration");
  return c;
}

struct RadialTrace { int element, node; double radius, width; bool fluid; };
std::vector<RadialTrace> nativeTraces(const Full1D::SEM &sem) {
  std::vector<RadialTrace> out;
  for(int e=0;e<sem.mesh().NE();++e) for(int q=0;q<sem.mesh().NN();++q) {
    const double r=sem.mesh().NodeRadius(e,q)*sem.lengthScaleMetres();
    if(r>0.0) out.push_back({e,q,r,(sem.mesh().EUR(e)-sem.mesh().ELR(e))*sem.lengthScaleMetres(),sem.elementIsFluid(e)});
  }
  return out;
}

Eigen::MatrixXcd baseline(PhysicalContext &c, const Options &o, int endBin, double &dcImag,
                          dssolver::BatchResponseCache *batchCache=nullptr,
                          const std::function<void(int,const SPARSESPEC::PoleAngularCoefficientRecord3D&)> &recordSink={},
                          std::vector<double> *basisResiduals=nullptr) {
  auto &sem=*c.sem; auto &source=c.input->cmt(); auto &parameters=c.input->inputParameters();
  const auto &freq=c.input->freqFull();
  const int nf=freq.nt()/2+1;
  Eigen::VectorXd zeroWindow=Eigen::VectorXd::Zero(freq.nt());
  SPARSESPEC::WindowedDisplacementMap measurement(freq,sem.lengthScaleMetres(),zeroWindow);
  Eigen::MatrixXcd q=Eigen::MatrixXcd::Zero(3*parameters.num_receivers(),nf);
  Eigen::MatrixXcd out=Eigen::MatrixXcd::Zero(3*parameters.num_receivers(),nf);
  const int e=sem.receiverElements(parameters,false).at(0);
  const auto station=parameters.receivers().at(0);
  const double rad=kPi/180.0;
  const double r=(sem.mesh().PR()-1000.0*parameters.receiver_depth()/sem.lengthScaleMetres())*sem.lengthScaleMetres();
  const std::vector<SPARSESPEC::Field3DPoint> receiver{{r,(90.0-station.first)*rad,station.second*rad,e}};
  const auto w=freq.w();
  const auto traces=recordSink?nativeTraces(sem):std::vector<RadialTrace>{};
  const auto stationTheta=(90.0-station.first)*rad;
  const auto stationPhi=station.second*rad;
  const std::vector<SPARSESPEC::Field3DPoint> receiverPoints{{r,stationTheta,stationPhi,e}};
  for(int k=o.includeDc?0:1;k<endBin;++k) {
    mathprofile::Scope bin("forward_frequency_bin");
    const Complex z=w.at(k)+Complex(0.0,-freq.ep());
    Eigen::VectorXcd unitWeights=Eigen::VectorXcd::Zero(3*parameters.num_receivers());
    if(recordSink)unitWeights[0]=Complex(1.0,0.0);
    const Eigen::VectorXcd weights=recordSink?unitWeights:Eigen::VectorXcd(q.col(k));
    double basisResidual=0.0;
    auto workspace=SPARSESPEC::Windowed3DFieldEngine::prepareFrequencyWorkspace(
        sem,source,parameters,z,weights,batchCache,o.poleOnly,
        (batchCache||recordSink)?false:o.forwardOnly,recordSink?&basisResidual:nullptr);
    if(recordSink&&basisResiduals)basisResiduals->at(static_cast<std::size_t>(k))=basisResidual;
    const auto field=workspace.sampleForKernelsPacked(receiver,0,nullptr,2).at(0);
    const Complex fz=Complex(0.0,-1.0)/z;
    for(int component=0;component<3;++component)
      out(component,k)=fz*field.forwardRawDisplacement[component];
    if(recordSink) {
      constexpr double captureTheta=1.23456789012345,capturePhi=0.45678901234567;
      std::vector<SPARSESPEC::Field3DPoint> points;
      points.reserve(traces.size()+receiverPoints.size());
      for(const auto &trace:traces)
        points.push_back({trace.radius,captureTheta,capturePhi,trace.element});
      points.insert(points.end(),receiverPoints.begin(),receiverPoints.end());
      auto capture=[&](const SPARSESPEC::PoleAngularCoefficientRecord3D &record){recordSink(k,record);};
      (void)workspace.sampleForKernelsTensor(points,traces.size(),nullptr,2,capture,false);
    }
  }
  if(!out.allFinite()) throw std::runtime_error("baseline displacement spectrum contains nonfinite values");
  dcImag=o.includeDc?std::abs(out(0,0).imag()):std::numeric_limits<double>::quiet_NaN();
  return out;
}

Eigen::VectorXd filteredVerticalTime(const Eigen::MatrixXcd &raw,const SpectraSolver::FreqFull &freq,
                                    double scale,int usedBins) {
  Eigen::MatrixXcd band=Eigen::MatrixXcd::Zero(1,freq.nt()/2+1);
  band.block(0,0,1,usedBins)=raw.block(0,0,1,usedBins);
  for(int k=0;k<usedBins;++k) band(0,k)*=bandWeight(freq.f(k)/freq.timeNorm());
  projectRealFftEndpoints(band);
  Eigen::MatrixXd time=processfunctions::freq2time(band,freq);
  return time.row(0).transpose()*scale;
}

double energy(const Eigen::VectorXd &u,const Eigen::VectorXd &sqrtWindow,double dt) {
  return 0.5*dt*(u.array().square()*sqrtWindow.array().square()).sum();
}

void writeTrace(const fs::path &path,const Eigen::VectorXd &u,const Eigen::VectorXd &window,double dt) {
  mathprofile::Scope io("trace_csv_write"); std::ofstream f(path); f<<"time_s,vertical_m,window_weight\n"<<std::setprecision(17);
  for(Eigen::Index i=0;i<u.size();++i) f<<i*dt<<','<<u[i]<<','<<window[i]*window[i]<<'\n';
}

void writeNativeTraces(const fs::path &path,const std::vector<RadialTrace> &traces) {
  mathprofile::Scope io("native_radial_trace_metadata_write");
  std::ofstream f(path); f<<"element,node,radius_m,element_width_m,fluid\n"<<std::setprecision(17);
  for(const auto &r:traces) f<<r.element<<','<<r.node<<','<<r.radius<<','<<r.width<<','<<(r.fluid?1:0)<<'\n';
}

struct EnergyAdjointAudit { Eigen::MatrixXcd weights; double qdotp=0.0, relativeError=0.0; };

EnergyAdjointAudit physicalEnergyAdjoint(PhysicalContext &c,const Eigen::MatrixXcd &raw,
                                         int usedBins,const Eigen::VectorXd &win,double energyValue) {
  auto &sem=*c.sem; auto &freq=c.input->freqFull();
  SPARSESPEC::WindowedDisplacementMap measurement(freq,sem.lengthScaleMetres(),win);
  Eigen::MatrixXcd band=Eigen::MatrixXcd::Zero(raw.rows(),raw.cols());
  for(int k=0;k<usedBins;++k) band.col(k)=raw.col(k)*bandWeight(freq.f(k)/freq.timeNorm());
  const Eigen::MatrixXcd processed=measurement.forward(band);
  Eigen::MatrixXcd q=Eigen::MatrixXcd::Zero(raw.rows(),raw.cols());
  const double dfSI=freq.df()/freq.timeNorm();
  for(int k=0;k<q.cols();++k) q(0,k)=dfSI*((k==0 || ((freq.nt()%2)==0 && k==raw.cols()-1))?1.0:2.0)*processed(0,k);
  projectRealFftEndpoints(q);
  const double qdotp=(q.conjugate().cwiseProduct(processed)).sum().real();
  const double relativeError=std::abs(qdotp-2.0*energyValue)/std::abs(2.0*energyValue);
  if(!std::isfinite(qdotp)||!std::isfinite(relativeError)||!(energyValue>0.0)||relativeError>=1e-10)
    throw std::runtime_error("Energy derivative identity q dot P(B raw) = 2E failed its 1e-10 relative gate.");
  Eigen::MatrixXcd qhat=measurement.adjoint(q);
  for(int k=0;k<qhat.cols();++k) qhat.col(k)*=bandWeight(freq.f(k)/freq.timeNorm());
  return {std::move(qhat),qdotp,relativeError};
}

EnergyAdjointAudit runKernel(PhysicalContext &c,const Options &o,const Eigen::MatrixXcd &raw,
                             int usedBins,const Eigen::VectorXd &win,double energyValue,
                             dssolver::BatchResponseCache *batchCache=nullptr) {
  auto &sem=*c.sem; auto &source=c.input->cmt(); auto &parameters=c.input->inputParameters();
  auto &freq=c.input->freqFull();
  auto traces=nativeTraces(sem);
  const std::size_t take=o.radialTraceLimit==0?traces.size():std::min<std::size_t>(o.radialTraceLimit,traces.size());
  std::vector<RadialTrace> selected;
  for(std::size_t i=0;i<take;++i) selected.push_back(traces[take==1?traces.size()-1:(i*(traces.size()-1))/std::max<std::size_t>(1,take-1)]);
  if(selected.empty()) throw std::runtime_error("no native radial element-node sides are available");
  auto energyAdjoint=[&]{mathprofile::Scope t("physical_energy_adjoint_all_output_bins");
    return physicalEnergyAdjoint(c,raw,usedBins,win,energyValue);}();
  const int nTheta=2*o.lmax+1, nPhi=GSHTrans::FastFFTSize(4*o.lmax+1);
  auto quadrature=GaussQuad::Quadrature1D<double>(GaussQuad::LegendrePolynomial<double>{}.GaussQuadrature(nTheta));
  quadrature.Transform([](double x){return std::acos(-x);},[](double){return 1.0;});
  const std::size_t totalAngles=static_cast<std::size_t>(nTheta)*static_cast<std::size_t>(nPhi);
  const bool probe=o.mode=="kernel-probe";
  if(probe&&o.probeAngleCount>totalAngles)
    throw std::runtime_error("kernel-probe angle count exceeds the complete native angular grid");
  const std::size_t sampledAngles=probe?o.probeAngleCount:totalAngles;
  if(totalAngles>std::numeric_limits<std::size_t>::max()/selected.size()/6/sizeof(double))
    throw std::runtime_error("kernel result size overflows size_t");
  if(sampledAngles>std::numeric_limits<std::size_t>::max()/selected.size()/6/sizeof(double))
    throw std::runtime_error("sampled kernel result size overflows size_t");
  const std::size_t targetPointCount=totalAngles*selected.size();
  const std::size_t fullTargetResultBytes=targetPointCount*6*sizeof(double);
  const std::size_t pointCount=sampledAngles*selected.size(), resultBytes=pointCount*6*sizeof(double);
  if(totalAngles>std::numeric_limits<std::size_t>::max()/(5*sizeof(double)))
    throw std::runtime_error("source-contact result size overflows size_t");
  const std::size_t fullTargetContactBytes=o.includeContacts?totalAngles*5*sizeof(double):0;
  const std::size_t contactBytes=o.includeContacts?sampledAngles*5*sizeof(double):0;
  if(contactBytes>std::numeric_limits<std::size_t>::max()-resultBytes)
    throw std::runtime_error("combined kernel/contact result size overflows size_t");
  const std::size_t requiredBytes=resultBytes+contactBytes;
  if(fullTargetContactBytes>std::numeric_limits<std::size_t>::max()-fullTargetResultBytes)
    throw std::runtime_error("full-target kernel/contact result size overflows size_t");
  const std::size_t fullTargetCombinedBytes=fullTargetResultBytes+fullTargetContactBytes;
  const double bytesPerGiB=1024.0*1024.0*1024.0;
  if(o.maxOutputGiB>static_cast<double>(std::numeric_limits<std::size_t>::max())/bytesPerGiB)
    throw std::runtime_error("--max-output-gib exceeds addressable file size");
  const std::size_t outputCapBytes=static_cast<std::size_t>(o.maxOutputGiB*bytesPerGiB);
  const auto space=fs::space(o.output);
  const std::uintmax_t availableBytes=space.available;
  if(requiredBytes>outputCapBytes || requiredBytes>availableBytes)
    throw std::runtime_error("kernel output preflight rejected combined binary size "+std::to_string(requiredBytes)+
        " bytes (configured cap "+std::to_string(outputCapBytes)+", available disk "+
        std::to_string(availableBytes)+"); select fewer radial sides or raise --max-output-gib");
  constexpr std::size_t kAngleTile=64;
  const fs::path kernelPath=o.output/(probe?"kernel_probe.bin":"kernel.bin");
  const fs::path contactPath=o.output/(probe?"source_contact_probe.bin":"source_contact.bin");
  { mathprofile::Scope io("kernel_result_file_open_and_zero_initialization");
    std::ofstream f(kernelPath,std::ios::binary|std::ios::trunc);
    std::vector<char> zeros(1u<<20,0); std::size_t remaining=resultBytes;
    while(remaining){const std::size_t n=std::min(remaining,zeros.size());f.write(zeros.data(),n);remaining-=n;}
    if(!f) throw std::runtime_error("failed to initialize streamed kernel result");
  }
  if(o.includeContacts) {
    mathprofile::Scope io("source_contact_result_file_zero_initialization");
    std::ofstream f(contactPath,std::ios::binary|std::ios::trunc);
    std::vector<char> zeros(1u<<20,0); std::size_t remaining=contactBytes;
    while(remaining){const std::size_t n=std::min(remaining,zeros.size());f.write(zeros.data(),n);remaining-=n;}
    if(!f) throw std::runtime_error("failed to initialize source-contact result");
  }
  { mathprofile::Scope io("kernel_grid_and_radial_metadata_write");
    std::ofstream a(o.output/"angular_grid.csv"); a<<"theta_index,theta_rad,theta_weight\n"<<std::setprecision(17);
    for(int i=0;i<nTheta;++i) a<<i<<','<<quadrature.Points()[i]<<','<<quadrature.Weights()[i]<<'\n';
    if(probe) {
      std::ofstream g(o.output/"kernel_probe_angles.csv");
      g<<"sample_index,theta_index,phi_index,theta_rad,phi_rad,theta_weight,phi_weight\n"<<std::setprecision(17);
      for(std::size_t ai=0;ai<sampledAngles;++ai) {
        const std::size_t itheta=ai/static_cast<std::size_t>(nPhi), iphi=ai%static_cast<std::size_t>(nPhi);
        g<<ai<<','<<itheta<<','<<iphi<<','<<quadrature.Points()[itheta]<<','<<2.0*kPi*iphi/nPhi<<','<<quadrature.Weights()[itheta]<<','<<2.0*kPi/nPhi<<'\n';
      }
    }
    std::ofstream r(o.output/"kernel_radial_nodes.csv"); r<<"radial_index,element,node,radius_m,element_width_m,fluid\n"<<std::setprecision(17);
    for(std::size_t i=0;i<selected.size();++i){const auto &x=selected[i];r<<i<<','<<x.element<<','<<x.node<<','<<x.radius<<','<<x.width<<','<<(x.fluid?1:0)<<'\n';}
  }
  std::fstream out(kernelPath,std::ios::binary|std::ios::in|std::ios::out);
  if(!out) throw std::runtime_error("cannot reopen streamed kernel result");
  std::fstream contactOut;
  if(o.includeContacts){contactOut.open(contactPath,std::ios::binary|std::ios::in|std::ios::out);
    if(!contactOut) throw std::runtime_error("cannot reopen source-contact result");}
  SPARSESPEC::Windowed3DFieldEngine::SamplingPlanCache plan;
  std::size_t maxRadialCoefficientCacheBytes=0;
  for(int k=0;k<usedBins;++k) {
    if((k==0&&!o.includeDc)||bandWeight(freq.f(k)/freq.timeNorm())==0.0) continue;
    const Complex z=freq.w().at(k)+Complex(0.0,-freq.ep());
    auto workspace=[&]{ mathprofile::Scope t("frequency_workspace_modal_solve_and_adjoint_construction");
      return SPARSESPEC::Windowed3DFieldEngine::prepareFrequencyWorkspace(sem,source,parameters,z,energyAdjoint.weights.col(k),batchCache,o.poleOnly); }();
    for(std::size_t angleBegin=0;angleBegin<sampledAngles;angleBegin+=kAngleTile) {
      const std::size_t angleEnd=std::min(sampledAngles,angleBegin+kAngleTile), angleN=angleEnd-angleBegin;
      std::vector<SPARSESPEC::Field3DPoint> points; points.reserve(angleN*selected.size());
      std::vector<std::array<double,2>> angles; angles.reserve(angleN);
      for(std::size_t ai=angleBegin;ai<angleEnd;++ai) {
        const std::size_t itheta=ai/static_cast<std::size_t>(nPhi), iphi=ai%static_cast<std::size_t>(nPhi);
        const double theta=quadrature.Points()[itheta], phi=2.0*kPi*iphi/nPhi;
        angles.push_back({theta,phi});
        for(const auto &r:selected) points.push_back({r.radius,theta,phi,r.element});
      }
      auto fields=[&]{ mathprofile::Scope t("harmonic_synthesis_native_node_sampling");
        return workspace.sampleForKernelsTensor(points,points.size(),&plan,2,{},false); }();
      std::vector<double> sums(fields.size()*6);
      { mathprofile::Scope t("pointwise_strain_rotation_kernel_contraction");
        for(std::size_t i=0;i<fields.size();++i){const auto v=SPARSESPEC::Windowed3DKernelEngine::pointKernel(sem,fields[i]);
          const std::array<double,6> x{v.A,v.C,v.F,v.L,v.N,v.kappa};
          std::copy(x.begin(),x.end(),sums.begin()+i*6);}
      }
      if(o.includeContacts) {
        const auto contacts=[&]{mathprofile::Scope t("source_contact_separate_case");
          return SPARSESPEC::Windowed3DKernelEngine::sourceShellContact(workspace,source,angles);}();
        std::vector<double> values(contacts.size()*5);
        for(std::size_t i=0;i<contacts.size();++i){const auto &v=contacts[i];
          values[i*5+0]=v.A;values[i*5+1]=v.C;values[i*5+2]=v.F;values[i*5+3]=v.L;values[i*5+4]=v.N;}
        const std::size_t nBytes=values.size()*sizeof(double), offset=angleBegin*5*sizeof(double);
        std::vector<double> previous(values.size());
        {mathprofile::Scope t("source_contact_frequency_sum_read");contactOut.seekg(static_cast<std::streamoff>(offset));
          contactOut.read(reinterpret_cast<char*>(previous.data()),static_cast<std::streamsize>(nBytes));
          if(!contactOut)throw std::runtime_error("failed reading source-contact accumulation tile");}
        for(std::size_t i=0;i<values.size();++i)values[i]+=previous[i];
        {mathprofile::Scope t("source_contact_frequency_sum_write");contactOut.seekp(static_cast<std::streamoff>(offset));
          contactOut.write(reinterpret_cast<const char*>(values.data()),static_cast<std::streamsize>(nBytes));
          if(!contactOut)throw std::runtime_error("failed writing source-contact accumulation tile");}
      }
      const std::size_t firstPoint=angleBegin*selected.size(), nBytes=sums.size()*sizeof(double), offset=firstPoint*6*sizeof(double);
      std::vector<double> previous(sums.size());
      { mathprofile::Scope t("kernel_frequency_sum_read"); out.seekg(static_cast<std::streamoff>(offset));out.read(reinterpret_cast<char*>(previous.data()),static_cast<std::streamsize>(nBytes));
        if(!out) throw std::runtime_error("failed reading kernel accumulation tile"); }
      for(std::size_t i=0;i<sums.size();++i) sums[i]+=previous[i];
      { mathprofile::Scope t("kernel_frequency_sum_write"); out.seekp(static_cast<std::streamoff>(offset));out.write(reinterpret_cast<const char*>(sums.data()),static_cast<std::streamsize>(nBytes));
        if(!out) throw std::runtime_error("failed writing kernel accumulation tile"); }
    }
    maxRadialCoefficientCacheBytes=std::max(maxRadialCoefficientCacheBytes,workspace.radialCoefficientCacheBytes());
  }
  out.close(); contactOut.close();
  if(fs::file_size(kernelPath)!=resultBytes || (o.includeContacts&&fs::file_size(contactPath)!=contactBytes))
    throw std::runtime_error("streamed kernel output size differs from the preflighted sampled subset size");
  mathprofile::count("sampling_plan_cache_bytes_final",plan.bytes);
  mathprofile::count("tensor_radial_grid_cache_bytes",plan.tensorGridPayloadBytes);
  mathprofile::count("tensor_angular_plan_cache_bytes",plan.tensorAngularPayloadBytes);
  mathprofile::count("radial_source_coefficient_cache_bytes_max",maxRadialCoefficientCacheBytes);
  { std::ofstream m(o.output/(probe?"kernel_probe.json":"kernel_domain.json")); m<<"{\"domain\":\"selected_native_SEM_element_node_sides\",\"selected_radial_trace_count\":"<<selected.size()
    <<",\"available_radial_trace_sides\":"<<traces.size()<<",\"angular_rule\":\"GSHTrans ForBand(L,2,2.0) equivalent nodes, streamed without a global transform grid\",\"n_theta\":"<<nTheta<<",\"n_phi\":"<<nPhi<<",\"full_target_angular_points_per_trace\":"<<totalAngles
    <<",\"sampled_angular_points_per_trace\":"<<sampledAngles<<",\"sampled_angle_subset_is_complete\":"<<(sampledAngles==totalAngles?"true":"false")
    <<",\"full_target_point_count\":"<<targetPointCount<<",\"actual_sampled_point_count\":"<<pointCount<<",\"full_target_kernel_binary_bytes\":"<<fullTargetResultBytes<<",\"actual_kernel_binary_bytes\":"<<resultBytes
    <<",\"full_target_contact_binary_bytes\":"<<fullTargetContactBytes<<",\"actual_contact_binary_bytes\":"<<contactBytes<<",\"full_target_combined_output_bytes\":"<<fullTargetCombinedBytes<<",\"actual_combined_output_bytes\":"<<requiredBytes
    <<",\"configured_output_cap_bytes\":"<<outputCapBytes<<",\"available_disk_bytes_at_preflight\":"<<availableBytes<<",\"binary_file\":\""<<(probe?"kernel_probe.bin":"kernel.bin")<<"\",\"binary_layout\":\"point-major; actual sampled angle_index outer, selected radial side inner, doubles[A,C,F,L,N,kappa] summed over frequency\",\"kernel_is_full_earth_volume_integral\":false,\"diagnostic_angle_subset_only\":"<<(probe?"true":"false")
    <<",\"source_contact_enabled\":"<<(o.includeContacts?"true":"false")<<",\"source_contact_file\":\""<<(probe?"source_contact_probe.bin":"source_contact.bin")<<"; actual sampled angles only; angle-major doubles[A,C,F,L,N] summed over frequency; separate source-shell contact, not volume\",\"energy_adjoint_qdotp\":"<<energyAdjoint.qdotp<<",\"energy_adjoint_relative_error_vs_2E\":"<<energyAdjoint.relativeError
    <<",\"sampling_plan_cache_budget_bytes\":"<<SPARSESPEC::Windowed3DFieldEngine::SamplingPlanCache::byteBudget<<",\"sampling_plan_cache_bytes_final\":"<<plan.bytes<<",\"tensor_grid_cache_bytes\":"<<plan.tensorGridPayloadBytes<<",\"tensor_angular_cache_bytes\":"<<plan.tensorAngularPayloadBytes
    <<",\"pole_coefficient_cache_budget_bytes\":"<<SPARSESPEC::Windowed3DFieldEngine::poleCoefficientCacheByteBudget<<",\"radial_source_coefficient_cache_bytes_max\":"<<maxRadialCoefficientCacheBytes
    <<",\"max_angle_tile\":"<<kAngleTile<<",\"known_driver_tile_point_count_max\":"<<std::min(kAngleTile,sampledAngles)*selected.size()<<",\"known_driver_tile_kernel_sample_payload_bytes_max\":"<<std::min(kAngleTile,sampledAngles)*selected.size()*sizeof(SPARSESPEC::KernelField3DSample)<<",\"known_driver_tile_kernel_sum_payload_bytes_max\":"<<std::min(kAngleTile,sampledAngles)*selected.size()*6*sizeof(double)<<",\"memory_note\":\"cache budget/retention and known driver tile payloads are reported separately; these are not a process peak-memory estimate and exclude transform workspaces and transient field internals\"}\n"; }
  return energyAdjoint;
}

} // namespace

int main(int argc,char **argv) try {
  const auto processStart=Clock::now(); const Options o=parse(argc,argv);
  mathprofile::Root topProfile("dspec",-1);
  mathprofile::Scope totalProfile("physical_case_total");
  auto c=prepare(o); auto &freq=c.input->freqFull(); auto &sem=*c.sem;
  const int nf=freq.nt()/2+1;
  int last=0; while(last+1<nf && freq.f(last+1)/freq.timeNorm()<0.010) ++last;
  const int endBin=std::min(nf,last+2); // retain first out-of-band guard bin
  if(!o.includeDc && (o.mode=="kernel"||o.mode=="kernel-probe")) throw std::runtime_error("kernel modes require --dc include for the physical energy derivative");
  if(o.includeContacts && o.poleOnly) throw std::runtime_error("--contacts on requires --pole-only off because source contacts use global modal blocks");
  double dcImag=0.0;
  Eigen::MatrixXcd raw;
  const bool kernelMode=o.mode=="kernel"||o.mode=="kernel-probe";
  if(o.responseRecordSpool&&(!kernelMode||o.includeContacts||o.radialTraceLimit!=0))
    throw std::runtime_error("--response-record-spool on requires kernel mode, contacts off, and --radial-traces 0");
  const bool responseSpoolUsed=kernelMode&&o.responseRecordSpool;
  bool responseReuseEnabled=kernelMode&&o.responseReuse&&!responseSpoolUsed;
  bool responseReuseCapacityFallback=false;
  bool responseReusePreflightFallback=false;
  dssolver::BatchResponseCache::Limits responseLimits;
  const double bytesPerGiB=1024.0*1024.0*1024.0;
  if(o.responseCacheGiB>static_cast<double>(std::numeric_limits<std::size_t>::max())/bytesPerGiB)
    throw std::runtime_error("--response-cache-gib exceeds addressable size");
  responseLimits.maximumBytes=static_cast<std::size_t>(o.responseCacheGiB*bytesPerGiB);
  std::optional<SPARSESPEC::Windowed3DFieldEngine::BatchResponseFootprint> responseFootprint;
  std::unique_ptr<dssolver::BatchResponseCache> responseCache;
  std::unique_ptr<dspec_global_gsh::PoleResponseSpoolWriter> responseSpool;
  std::vector<double> responseBasisResiduals(static_cast<std::size_t>(endBin),0.0);
  if(responseSpoolUsed) {
    if(o.responseSpoolGiB>static_cast<double>(std::numeric_limits<std::uint64_t>::max())/bytesPerGiB)
      throw std::runtime_error("--response-spool-gib exceeds addressable size");
    const auto traces=nativeTraces(sem);
    std::vector<dspec_global_gsh::NativeSide> sides;
    sides.reserve(traces.size());
    for(const auto &trace:traces)sides.push_back({trace.element,trace.node,trace.radius,
        trace.radius/sem.lengthScaleMetres(),trace.width,trace.fluid});
    const auto cap=static_cast<std::uint64_t>(o.responseSpoolGiB*bytesPerGiB);
    responseSpool=std::make_unique<dspec_global_gsh::PoleResponseSpoolWriter>(
        (o.output/"pole_response_spool.bin").string(),o.lmax,static_cast<std::size_t>(endBin),
        std::move(sides),cap);
  }
  if(responseSpoolUsed) {
    const auto sink=[&](int bin,const SPARSESPEC::PoleAngularCoefficientRecord3D &record){
      responseSpool->append(static_cast<std::size_t>(bin),record);
    };
    mathprofile::Scope t("baseline_forward_and_compact_response_capture");
    raw=baseline(c,o,endBin,dcImag,nullptr,sink,&responseBasisResiduals);
  } else if(responseReuseEnabled) {
    const Complex firstZ=freq.w().at(o.includeDc?0:1)+Complex(0.0,-freq.ep());
    responseFootprint=SPARSESPEC::Windowed3DFieldEngine::estimateBatchResponseFootprint(
        sem,c.input->cmt(),c.input->inputParameters(),firstZ,
        static_cast<std::uint64_t>(endBin-(o.includeDc?0:1)));
    if(responseFootprint->entries>std::numeric_limits<std::size_t>::max()||
       responseFootprint->bytes>responseLimits.maximumBytes) {
      responseReuseEnabled=false;
      responseReuseCapacityFallback=true;
      responseReusePreflightFallback=true;
      mathprofile::count("response_reuse_preflight_ineligible");
    } else {
      responseLimits.maximumEntries=static_cast<std::size_t>(responseFootprint->entries);
      responseCache=std::make_unique<dssolver::BatchResponseCache>(responseLimits);
      mathprofile::count("response_cache_preflight_estimated_bytes",responseFootprint->bytes);
      mathprofile::count("response_cache_preflight_entries",responseFootprint->entries);
      mathprofile::count("response_cache_preflight_max_dofs",responseFootprint->maximumDofs);
    }
  }
  if(responseReuseEnabled) {
    try {
      mathprofile::Scope t("baseline_forward_with_bounded_response_capture");
      raw=baseline(c,o,endBin,dcImag,responseCache.get());
    } catch(const dssolver::BatchResponseCache::CapacityExceeded &) {
      // Discard a partial cache and rerun the baseline through the reference
      // forward-only path; no partial response set is used by the kernel.
      responseCache->clear();
      responseCache.reset();
      responseReuseEnabled=false;
      responseReuseCapacityFallback=true;
      mathprofile::count("response_reuse_capacity_fallbacks");
      mathprofile::Scope t("baseline_forward_capacity_fallback");
      raw=baseline(c,o,endBin,dcImag,nullptr);
    }
  } else {
    mathprofile::Scope t("baseline_forward_receiver_only");
    raw=baseline(c,o,endBin,dcImag,nullptr);
  }
  const auto win=sqrtHann(freq); const double dt=freq.dt()*freq.timeNorm();
  writeNativeTraces(o.output/"radial_nodes.csv",nativeTraces(sem));
  Eigen::VectorXd u=[&]{mathprofile::Scope t("measurement_filter_and_inverse_fft");
    return filteredVerticalTime(raw,freq,sem.lengthScaleMetres(),endBin);}();
  Eigen::VectorXd unwindowed=u;
  const double e=[&]{mathprofile::Scope t("windowed_energy_evaluation");return energy(u,win,dt);}();
  const double dcScale=raw.row(0).head(endBin).cwiseAbs().maxCoeff();
  const double dcImagRelative=o.includeDc?dcImag/std::max(dcScale,std::numeric_limits<double>::min()):
      std::numeric_limits<double>::quiet_NaN();
  if(o.includeDc&&(!std::isfinite(dcImagRelative)||dcImagRelative>=1e-10))
    throw std::runtime_error("normalized raw DC imaginary residual failed its 1e-10 scale-aware gate");
  std::optional<EnergyAdjointAudit> energyAudit;
  std::optional<dspec_global_gsh::ConsumerResult> spoolKernelResult;
  if(responseSpoolUsed) {
    energyAudit=physicalEnergyAdjoint(c,raw,endBin,win,e);
    for(int k=0;k<endBin;++k) {
      const Complex z=freq.w().at(k)+Complex(0.0,-freq.ep());
      responseSpool->setAudit(static_cast<std::size_t>(k),{k,z,energyAudit->weights(0,k),
          responseBasisResiduals.at(static_cast<std::size_t>(k))});
    }
    responseSpool->finalize();
    mathprofile::count("response_spool_bytes",responseSpool->bytesWritten());
    dspec_global_gsh::ConsumerOptions consumerOptions;
    consumerOptions.sampleAngleCount=o.mode=="kernel-probe"?o.probeAngleCount:0;
    consumerOptions.workers=o.globalWorkers;
    consumerOptions.longitudeTile=o.globalAngleTile;
    consumerOptions.maxOutputGiB=o.maxOutputGiB;
    spoolKernelResult=dspec_global_gsh::runPoleResponseSpoolKernel(
        (o.output/"pole_response_spool.bin").string(),sem,(o.output/"kernel.bin").string(),consumerOptions);
  } else if(kernelMode) energyAudit=runKernel(c,o,raw,endBin,win,e,
      responseReuseEnabled?responseCache.get():nullptr);
  { mathprofile::Scope io("physical_result_and_spectrum_csv");
    std::ofstream spectrum(o.output/"spectrum.csv"); spectrum<<"bin,frequency_hz,real,imag\n"<<std::setprecision(17);
    for(int k=o.includeDc?0:1;k<=last;++k) spectrum<<k<<','<<freq.f(k)/freq.timeNorm()<<','<<raw(0,k).real()*sem.lengthScaleMetres()*freq.timeNorm()<<','<<raw(0,k).imag()*sem.lengthScaleMetres()*freq.timeNorm()<<'\n';
  }
  writeTrace(o.output/"trace.csv",unwindowed,win,dt);
  { std::ofstream m(o.output/"result.json"); m<<std::setprecision(17)
    <<"{\"mode\":\""<<o.mode<<"\",\"degree\":"<<o.lmax<<",\"sem_degree\":"<<o.semDegree
    <<",\"nq_gll_points_per_element\":"<<o.semDegree<<",\"polynomial_degree\":"<<o.semDegree-1
    <<",\"radial_refinement_maxstep\":"<<o.radialRefinement<<",\"qex\":"<<o.qex<<",\"dt_s\":"<<dt<<",\"duration_s\":"<<o.duration
    <<",\"fft_samples\":"<<freq.nt()<<",\"fft_period_s\":"<<freq.nt()*dt<<",\"frequency_bins\":"<<last+1
    <<",\"epsilon_si_s_inverse\":"<<freq.ep()/freq.timeNorm()<<",\"df_si_hz\":"<<freq.df()/freq.timeNorm()
    <<",\"model_length_scale_m\":"<<sem.lengthScaleMetres()<<",\"model_time_scale_s\":"<<sem.timeScaleSeconds()
    <<",\"sem_elements\":"<<sem.mesh().NE()<<",\"gll_points_per_element\":"<<sem.mesh().NN()
    <<",\"native_element_node_sides\":"<<sem.mesh().NE()*sem.mesh().NN()
    <<",\"radial_dofs\":"<<sem.hR().rows()<<",\"spheroidal_dofs_max_l\":"<<sem.hS(o.lmax).rows()<<",\"toroidal_dofs_max_l\":"<<sem.hTk(o.lmax).rows()
    <<",\"source_element\":"<<sem.sourceElement(c.input->cmt())<<",\"source_element_fluid\":"<<(sem.elementIsFluid(sem.sourceElement(c.input->cmt()))?"true":"false")
    <<",\"receiver_element\":"<<sem.receiverElements(c.input->inputParameters(),false).at(0)<<",\"receiver_element_fluid\":"<<(sem.elementIsFluid(sem.receiverElements(c.input->inputParameters(),false).at(0))?"true":"false")
    <<",\"dc_included\":"<<(o.includeDc?"true":"false")<<",\"dc_imaginary_residual_raw\":";
    if(o.includeDc)m<<dcImag;else m<<"null";
    m<<",\"dc_imaginary_residual_relative_to_spectrum_scale\":";
    if(o.includeDc)m<<dcImagRelative;else m<<"null";
    m<<",\"nyquist_raw_solved\":false,\"nyquist_raw_imaginary_residual\":null"
    <<",\"pole_only_nonradial_adjoint\":"<<(o.poleOnly?"true":"false")<<",\"forward_only_baseline\":"<<((o.forwardOnly&&!responseReuseEnabled&&!responseSpoolUsed)?"true":"false")<<",\"source_contact_enabled\":"<<(o.includeContacts?"true":"false")
    <<",\"response_reuse_requested\":"<<(o.responseReuse?"true":"false")<<",\"response_reuse_used\":"<<(responseReuseEnabled?"true":"false")<<",\"response_reuse_capacity_fallback\":"<<(responseReuseCapacityFallback?"true":"false")<<",\"response_reuse_preflight_fallback\":"<<(responseReusePreflightFallback?"true":"false")<<",\"response_cache_budget_bytes\":"<<responseLimits.maximumBytes<<",\"response_cache_preflight_entries\":"<<(responseFootprint?responseFootprint->entries:0)<<",\"response_cache_preflight_estimated_bytes\":"<<(responseFootprint?responseFootprint->bytes:0)<<",\"response_cache_preflight_max_dofs\":"<<(responseFootprint?responseFootprint->maximumDofs:0)<<",\"response_cache_entries\":"<<(responseCache?responseCache->entryCount():0)<<",\"response_cache_bytes\":"<<(responseCache?responseCache->retainedBytes():0)<<",\"response_cache_hits\":"<<(responseCache?responseCache->cacheHits():0)<<",\"response_cache_misses\":"<<(responseCache?responseCache->cacheMisses():0)<<",\"response_basis_reconstruction_max_relative\":"<<(responseCache?responseCache->maximumBasisReconstructionRelative():0.0)
    <<",\"symbolic_reuse_enabled\":"<<(SPARSESPEC::FixedMeshAdjointSolver::symbolicReuseEnabled()?"true":"false")<<",\"sparse_symbolic_analysis_calls\":"<<mathprofile::counters["sparse_symbolic_analysis_calls"]<<",\"sparse_symbolic_analysis_reuses\":"<<mathprofile::counters["sparse_symbolic_analysis_reuses"]
    <<",\"response_record_spool_used\":"<<(responseSpoolUsed?"true":"false")<<",\"response_spool_bytes\":"<<(responseSpool?responseSpool->bytesWritten():0)<<",\"response_spool_kernel_elapsed_s\":"<<(spoolKernelResult?spoolKernelResult->elapsedSeconds:0.0)<<",\"response_spool_kernel_output_bytes\":"<<(spoolKernelResult?spoolKernelResult->outputBytes:0)<<",\"response_spool_kernel_sampled_angles\":"<<(spoolKernelResult?spoolKernelResult->angularPointCount:0)<<",\"global_consumer_workers\":"<<o.globalWorkers<<",\"global_consumer_angle_tile\":"<<o.globalAngleTile
    <<",\"spectrum_units\":\"m s; vertical radial displacement, SI damping convention\",\"source\":[10,20,50000],\"receiver\":[45,75,3000]"
    <<",\"vertical_component\":\"radial displacement, positive outward/up, component 0\",\"energy_m2_s\":"<<e<<",\"energy_definition\":\"0.5*dt*sum(Hann1000-1200*vertical_m^2)\",\"energy_adjoint_qdotp\":";
    if(energyAudit)m<<energyAudit->qdotp;else m<<"null";
    m<<",\"energy_adjoint_relative_error_vs_2E\":";
    if(energyAudit)m<<energyAudit->relativeError;else m<<"null";
    m<<",\"elapsed_s\":"<<std::chrono::duration<double>(Clock::now()-processStart).count()<<"}\n"; }
  fs::remove(o.output/"parameters.txt");
  return 0;
} catch(const std::exception &e) { std::cerr<<"FAIL: "<<e.what()<<'\n'; return 1; }
