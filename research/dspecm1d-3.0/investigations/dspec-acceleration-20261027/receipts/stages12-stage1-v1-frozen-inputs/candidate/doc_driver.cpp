// Matched documentation-case adapter: common model/measurement, independent solvers.
#include <YSpec/YSpec.hpp>
#include "Spectra.hpp"
#include "KernelParts.hpp"
#include "Sphere.hpp"
#include "matched_mesh_model.hpp"
#include "engine_headers/Windowed3DKernels.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <sstream>
#include <omp.h>
namespace fs=std::filesystem;
using Complex=std::complex<double>;
using Clock=std::chrono::steady_clock;
constexpr double pi=std::numbers::pi;
double seconds(Clock::time_point t){return std::chrono::duration<double>(Clock::now()-t).count();}
#ifndef MATHPROFILE_DISABLED
void writeProfileAccounting(const fs::path &path,
    std::vector<std::pair<std::string,mathprofile::Snapshot>> phases,
    const mathprofile::Snapshot &diagnostics,const mathprofile::Snapshot &final,
    double snapshotOverheadSeconds) {
 mathprofile::Scope reportScope("profile_snapshot_report_serialization_outside_pipeline");
 for(auto &[phase,delta]:phases)if(phase=="kernel_production_including_separately_isolated_diagnostics"){
  for(const auto &[name,value]:diagnostics.entries){
   auto found=delta.entries.find(name);if(found==delta.entries.end())std::terminate();
   if(found->second.inclusive<value.inclusive||found->second.exclusive<value.exclusive||found->second.calls<value.calls)std::terminate();
   found->second.inclusive-=value.inclusive;found->second.exclusive-=value.exclusive;found->second.calls-=value.calls;
   if(found->second.inclusive==0&&found->second.exclusive==0&&found->second.calls==0)delta.entries.erase(found);
  }
  for(const auto &[name,value]:diagnostics.counters){auto found=delta.counters.find(name);if(found==delta.counters.end()||found->second<value)std::terminate();found->second-=value;if(found->second==0)delta.counters.erase(found);}
  phase="kernel_production_excluding_isolated_diagnostics";
 }
 const auto pipeline=final.entries.find("profile_pipeline_wall");if(pipeline==final.entries.end())throw std::runtime_error("missing profile pipeline scope");
 const double pipelineWall=pipeline->second.inclusive*1e-9;
 long double partition=static_cast<long double>(pipeline->second.exclusive);
 for(const auto &[phase,delta]:phases){(void)phase;for(const auto &[name,value]:delta.entries)if(!name.starts_with("parallel_worker."))partition+=value.exclusive;}
 for(const auto &[name,value]:diagnostics.entries)if(!name.starts_with("parallel_worker."))partition+=value.exclusive;
 const double partitionSeconds=static_cast<double>(partition)*1e-9;
 const double residual=partitionSeconds-pipelineWall;
 const double tolerance=std::max(.01,.001*pipelineWall);
 if(!std::isfinite(pipelineWall)||!std::isfinite(partitionSeconds)||std::abs(residual)>tolerance)
  throw std::runtime_error("profile exclusive wall partition failed closure");
 auto component=[](const std::string &name){
  if(name.find("diagnostic")!=std::string::npos)return "diagnostic";
  if(name.find("kernel_fused")!=std::string::npos)return "material_evaluation_weak_form_contraction_accumulation";
  if(name.starts_with("batch_response_cache")||name.starts_with("batch_cache_"))return "response_cache";
  if(name.find("radial")!=std::string::npos)return "radial_sampling_reconstruction";
  if(name.find("receiver_basis")!=std::string::npos||name.find("receiver_basis_mix")!=std::string::npos||name.find("unit_adjoint_projection")!=std::string::npos)return "receiver_basis_projection";
  if(name.find("sparse_")!=std::string::npos||name.find("substitution")!=std::string::npos||name.find("HH_solve")!=std::string::npos)return "modal_solver";
  if(name.find("cache")!=std::string::npos||name.find("lookup")!=std::string::npos)return "sampling_plan_or_basis_cache";
  if(name.find("GEMM")!=std::string::npos||name.find("gemm")!=std::string::npos||name.find("synthesis")!=std::string::npos||name.find("tensor_")!=std::string::npos||name.find("basis")!=std::string::npos||name.find("Wigner")!=std::string::npos)return "tensor_angular_synthesis";
  if(name.find("io")!=std::string::npos||name.find("output")!=std::string::npos||name.find("trace")!=std::string::npos||name.find("spectrum")!=std::string::npos)return "output_io";
  if(name.find("driver_")!=std::string::npos||name.find("model")!=std::string::npos||name.find("sem_constructor")!=std::string::npos)return "driver_setup";
  return "unclassified";
 };
 std::ofstream out(path);if(!out)throw std::runtime_error("cannot create profile accounting output");
 out<<std::setprecision(17)<<"{\"schema\":1,\"snapshot_overhead_s\":"<<snapshotOverheadSeconds
    <<",\"report_serialization_scope\":\"profile_snapshot_report_serialization_outside_pipeline\""
    <<",\"wall_partition\":{\"profile_pipeline_wall_s\":"<<pipelineWall<<",\"master_exclusive_partition_s\":"<<partitionSeconds
    <<",\"residual_s\":"<<residual<<",\"tolerance_s\":"<<tolerance<<",\"worker_elapsed_excluded\":true},\"rows\":[";
 bool first=true;
 auto writeRows=[&](const std::string &phase,const mathprofile::Snapshot &snapshot,bool diagnostic){
  for(const auto &[name,value]:snapshot.entries){if(!first)out<<',';first=false;
   const bool worker=name.starts_with("parallel_worker.");
   out<<"{\"phase\":\""<<phase<<"\",\"scope_name\":\""<<name<<"\",\"execution_domain\":\""<<(worker?"worker_elapsed_sum":"master_wall")
      <<"\",\"inclusive_s\":"<<value.inclusive*1e-9<<",\"exclusive_s\":"<<value.exclusive*1e-9<<",\"calls\":"<<value.calls
      <<",\"logical_component\":\""<<component(name)<<"\",\"diagnostic\":"<<(diagnostic?"true":"false")
      <<",\"parent_overlap_description\":\""<<(worker?"worker elapsed sum overlaps the degree-region wall bucket and is excluded from master wall partition":name=="driver_kernel_pipeline_including_diagnostics"?"raw inclusive umbrella includes diagnostic parent; exclusive value excludes child diagnostics and is the wall partition term":"exclusive time is partitioned; inclusive time overlaps children and is not additive")<<"\"}";
  }
 };
 for(const auto &[phase,snapshot]:phases)writeRows(phase,snapshot,false);
 writeRows("diagnostic_tensor_validation",diagnostics,true);
 mathprofile::Snapshot pipelineOnly;pipelineOnly.entries.emplace("profile_pipeline_wall",pipeline->second);
 writeRows("profile_pipeline_total",pipelineOnly,false);
 out<<"],\"counters\":[";first=true;
 auto writeCounters=[&](const std::string &phase,const mathprofile::Snapshot &snapshot,bool diagnostic){for(const auto &[name,value]:snapshot.counters){if(!first)out<<',';first=false;out<<"{\"phase\":\""<<phase<<"\",\"name\":\""<<name<<"\",\"value\":"<<value<<",\"diagnostic\":"<<(diagnostic?"true":"false")<<"}";}};
 for(const auto &[phase,snapshot]:phases)writeCounters(phase,snapshot,false);
 writeCounters("diagnostic_tensor_validation",diagnostics,true);
 out<<"]}\n";out.close();reportScope.stop();
}
#endif
double profile(double r){return r>5.3e6&&r<5.8e6?std::pow(std::sin(pi*(r-5.3e6)/.5e6),2):0.;}
struct CommonModel {
 const YSpec::Model &model;double epsilon=0;
 auto values(int l,double r) const {
  auto m=model(l,r);double d=m.mu*std::expm1(2*epsilon*profile(r*LengthNorm()));
  m.F-=2*d;m.L+=d;m.N+=d;m.kappa-=4*d/3;m.mu+=d;return m;
 }
 double LengthNorm() const{return model.GetScales().Length();}
 double TimeNorm() const{return model.GetScales().Time();}
 double DensityNorm() const{return model.GetScales().Density();}
 double MassNorm() const{return DensityNorm()*std::pow(LengthNorm(),3);}
 double GravitationalConstant() const{return 1/(DensityNorm()*TimeNorm()*TimeNorm());}
 int NumberOfLayers() const{return model.NumberOfLayers();}
 double OuterRadius() const{return model.SurfaceRadius();}
 double LowerRadius(int l) const{return model.GetLayer(l).Lower();}
 double UpperRadius(int l) const{return model.GetLayer(l).Upper();}
 bool IsFluid(int l) const{return model.GetLayer(l).IsFluid();}
 bool IsSolid(int l) const{return !IsFluid(l);}
 auto LayerRadii(int l) const {auto r=model.GetLayer(l).Knots();return std::vector<double>(r.begin(),r.end());}
 #define PROPERTY(name,expression) auto name(int l) const{return [this,l](double r){auto m=values(l,r);return (expression);};}
 PROPERTY(Density,m.rho) PROPERTY(A,m.A) PROPERTY(C,m.C) PROPERTY(F,m.F) PROPERTY(L,m.L) PROPERTY(N,m.N)
 PROPERTY(Kappa,m.kappa) PROPERTY(Mu,m.mu)
 PROPERTY(VPV,std::sqrt(m.C/m.rho)) PROPERTY(VPH,std::sqrt(m.A/m.rho))
 PROPERTY(VSV,std::sqrt(std::max(0.,m.L/m.rho))) PROPERTY(VSH,std::sqrt(std::max(0.,m.N/m.rho)))
 PROPERTY(VP,std::sqrt((m.kappa+4*m.mu/3)/m.rho)) PROPERTY(VS,std::sqrt(std::max(0.,m.mu/m.rho)))
 PROPERTY(QKappa,m.qKappa>0?1/m.qKappa:0.) PROPERTY(QMu,m.qMu>0?1/m.qMu:0.)
 PROPERTY(Eta,m.F/(m.A-2*m.L)) PROPERTY(AAtten,m.A*m.qA) PROPERTY(CAtten,m.C*m.qA)
 PROPERTY(FAtten,m.F*m.qF) PROPERTY(LAtten,m.L*m.qMu) PROPERTY(NAtten,m.N*m.qMu)
 #undef PROPERTY
};
struct Node {int element,node,layer;double radius,dr;bool solid;};
std::string params(const YSpec::KernelInput &in,int degree){
 const auto &s=in.source;const auto &r=in.receiver;const auto &m=s.moment;
 auto lat=[](double x){return 90-YSpec::Sphere::Colatitude(x)*180/pi;};
 std::ostringstream o;o<<std::setprecision(17)<<"\"unused\"\n\""<<in.model.string()<<"\"\n4\n"<<(in.calculation.attenuation?1:0)<<"\n2\n0\n0\n1e-5\n0\n"<<degree<<"\n2\n10\n40\n2\n2\n4\n8\n10\n"<<s.depth/1000<<'\n'<<lat(s.location.latitude)<<'\n'<<s.location.longitude<<'\n'<<m.rr<<'\n'<<m.rt<<'\n'<<m.rp<<'\n'<<m.tt<<'\n'<<m.tp<<'\n'<<m.pp<<'\n'<<r.depth/1000<<"\n1\n"<<lat(r.location.latitude)<<' '<<r.location.longitude<<'\n';return o.str();
}
struct Measured {YSpec::Measurement measurement;std::vector<double> velocity,rate,window;double norm=0,delayAudit=0;};
Measured measure(const YSpec::KernelInput &in,const YSpec::Model &model,const YSpec::FrequencyGrid &grid,const std::vector<Complex> &raw){
 Measured out;auto rate=raw;for(std::size_t k=0;k<rate.size();++k)rate[k]*=Complex(0,1)*grid.Omega(k);
 out.velocity=YSpec::InverseTransform(grid,YSpec::MakeTaper(in.calculation,model.GetScales()),raw,1);
 out.rate=YSpec::InverseTransform(grid,YSpec::MakeTaper(in.calculation,model.GetScales()),rate,1);
 const auto &s=model.GetScales();double dt=in.calculation.sampling.timeStep;
 std::size_t samples=YSpec::SampleCount(grid,in.calculation.sampling);
 out.velocity.resize(samples);out.rate.resize(samples);out.measurement={in.receiver,in.measurements[0].direction,{}};
 for(std::size_t n=0;n<samples;++n){out.velocity[n]*=s.Velocity();out.rate[n]*=s.Velocity()/s.Time();double w=in.measurements[0].window(n*dt);out.window.push_back(w);out.norm+=w*out.rate[n]*out.rate[n]*dt;out.measurement.adjointSource.push_back(-w*out.rate[n]);}
 if(!(out.norm>0))throw std::runtime_error("zero traveltime normalisation");
 for(auto &j:out.measurement.adjointSource)j/=out.norm;
 auto adjoint=YSpec::MakeAdjointSpectrum(model,in.calculation,out.measurement);
 for(std::size_t k=0;k<raw.size();++k){Complex z=grid.Omega(k),q=std::conj(adjoint.weight[k]*adjoint.force[k]);out.delayAudit+=std::real(std::conj(q)*Complex(0,-1)/z*(-Complex(0,1)*z*raw[k]/s.Time()));}
 if(std::abs(out.delayAudit-1)>1e-10)throw std::runtime_error("discrete time-shift adjoint identity failed");
 return out;
}
// Complex weak-form contractions, with fixed effective Q, before real projection.
double kernelVs(const Full1D::SEM &sem,const SPARSESPEC::KernelField3DSample &sample,const YSpec::Material &m,Complex dispersion){
 if(sample.fluid)return 0;
 const auto &f=sample.forwardRawStrain;const auto &a=sample.adjointDualStrain;
 Complex F=-(std::conj(a.rr)*f.horizontalTrace+std::conj(a.horizontalTrace)*f.rr);
 Complex L=-2.*(std::conj(a.radialPlus)*f.radialPlus+std::conj(a.radialMinus)*f.radialMinus);
 Complex N=-(.5*(std::conj(a.horizontalPlus)*f.horizontalPlus+std::conj(a.horizontalMinus)*f.horizontalMinus)-std::conj(a.horizontalTrace)*f.horizontalTrace);
 double length=sem.lengthScaleMetres();
 return 2*m.mu*std::real(-2.*F*(1.+m.qF*dispersion)+(L+N)*(1.+m.qMu*dispersion))/std::pow(length,3);
}
int main(int argc,char **argv)try{
 if(argc!=10)throw std::runtime_error("usage: doc-driver dspec|yspec|dspec-forward INPUT EMPTY_OUTPUT L NQ REFINEMENT WORKERS sections|audit EPSILON");
 auto start=Clock::now();std::string backend=argv[1],gridKind=argv[8];fs::path out=fs::absolute(argv[3]);
 int degree=std::stoi(argv[4]),nq=std::stoi(argv[5]),workers=std::stoi(argv[7]);double refinement=std::stod(argv[6]),epsilon=std::stod(argv[9]);
 if(workers<1||workers>8||degree<1||degree>200||nq<6||nq>16||!(refinement>0)||std::abs(epsilon)>.001)throw std::runtime_error("bounded documentation scope exceeded");
 fs::create_directories(out);if(fs::directory_iterator(out)!=fs::directory_iterator{})throw std::runtime_error("output must be empty");
 Eigen::setNbThreads(1);omp_set_dynamic(0);omp_set_num_threads(1);
 mathprofile::Root root("documentation-wedge",-1);
 mathprofile::Scope pipelineScope("profile_pipeline_wall");
 mathprofile::Scope inputAndGridSetup("driver_input_and_grid_setup");
#ifndef MATHPROFILE_DISABLED
 std::vector<std::pair<std::string,mathprofile::Snapshot>> phaseSnapshots;
 mathprofile::Snapshot previousSnapshot,diagnosticTotals;
 double profileSnapshotOverheadSeconds=0.0;
 auto capturePhase=[&](const char *name){const auto before=Clock::now();auto current=mathprofile::snapshot();auto delta=mathprofile::difference(current,previousSnapshot);previousSnapshot=std::move(current);phaseSnapshots.emplace_back(name,std::move(delta));profileSnapshotOverheadSeconds+=seconds(before);};
#endif
 auto in=YSpec::ReadKernelInput(argv[2]);in.calculation.degrees.max=degree;
 if(in.calculation.quantity!=YSpec::Quantity::Velocity||in.measurements.size()!=1||in.measurements[0].component!="T")throw std::runtime_error("expected documentation transverse-velocity case");
 mathprofile::Scope modelSetup("driver_model_deck_and_common_material_setup");
 auto model=YSpec::Model::ReadDeck(in.model);CommonModel common{model,epsilon};
 modelSetup.stop();
 const auto &scales=model.GetScales();mathprofile::Scope gridSetup("driver_frequency_grid_setup");auto grid=YSpec::MakeGrid(in.calculation,scales);gridSetup.stop();
 auto geoLat=[](double lat){return 90-YSpec::Sphere::Colatitude(lat)*180/pi;};const auto &moment=in.source.moment;
 SourceInfo::EarthquakeCMT source(geoLat(in.source.location.latitude),in.source.location.longitude,in.source.depth/1000,moment.rr,moment.rt,moment.rp,moment.tt,moment.tp,moment.pp);
 mathprofile::Scope semSetup("driver_sem_constructor_mesh_static_matrices_and_model_setup");Full1D::SEM sem(common,refinement,nq,degree,source);semSetup.stop();
 mathprofile::Scope parameterIO("driver_parameters_file_io_and_input_parameters");std::ofstream pf(out/"parameters.txt");pf<<params(in,degree);pf.close();InputParameters parameters((out/"parameters.txt").string());parameterIO.stop();
 mathprofile::Scope nativeSideAudit("driver_native_side_sampling_and_gravity_audit");
 std::vector<Node> nodes;auto gll=GaussQuad::GaussLobattoLegendreQuadrature1D<double>(nq);double gravityError=0;
 for(int e=0;e<sem.mesh().NE();++e)for(int q=0;q<nq;++q){double r=sem.mesh().NodeRadius(e,q);int layer=sem.mesh().LayerNumber(e);gravityError=std::max(gravityError,std::abs(sem.meshModel().Gravity(e,q)-model(layer,r).g));if(r*scales.Length()>=3e6)nodes.push_back({e,q,layer,r*scales.Length(),.5*sem.mesh().EW(e)*gll.W(q)*scales.Length(),!sem.elementIsFluid(e)});}
 if(gravityError>1e-10)throw std::runtime_error("common-model gravity mismatch");
 nativeSideAudit.stop();
 mathprofile::Scope sectionAndReceiverSetup("driver_section_sampling_and_receiver_construction");
 std::vector<YSpec::KernelSampling> samplings(2);std::vector<double> angularWeights;
 const double spacing=360./(4*degree+1);
 if(gridKind=="sections"){
  YSpec::SectionOptions options{0,4,spacing,500e3};samplings[0]=YSpec::PathSection(model,in.source.location,in.receiver.location,options);samplings[1]=YSpec::CrossSection(model,in.source.location,in.receiver.location,options);
 }else if(gridKind=="audit"){
  auto quad=GaussQuad::Quadrature1D<double>(GaussQuad::LegendrePolynomial<double>{}.GaussQuadrature(2*degree+1));int np=4*degree+1;
  for(int t=0;t<quad.N();++t)for(int p=0;p<np;++p){samplings[0].places.push_back({std::acos(quad.X(t)),2*pi*p/np,0,0});angularWeights.push_back(quad.W(t)*2*pi/np);}samplings.resize(1);
 }else throw std::runtime_error("unknown sampling kind");
 for(auto &s:samplings){s.nodes.clear();for(auto n:nodes)s.nodes.push_back({n.layer,n.radius});}
 const auto dir=in.measurements[0].direction;Eigen::Vector3cd direction;direction<<dir.up,dir.north,dir.east;
 int receiverElement=sem.receiverElements(parameters,true).front();double receiverRadius=scales.Length()*model.SurfaceRadius()-in.receiver.depth;
 std::vector<SPARSESPEC::Field3DPoint> receiver{{receiverRadius,YSpec::Sphere::Colatitude(in.receiver.location.latitude),in.receiver.location.longitude*pi/180,receiverElement}};
 sectionAndReceiverSetup.stop();
 mathprofile::Scope cacheAndResultSetup("driver_result_cache_and_sampling_plan_initialisation");
 std::vector<Complex> raw(grid.Count());std::vector<std::array<Complex,3>> components(grid.Count());double forwardSeconds=0,kernelSeconds=0,measurementSeconds=0,tensorReferenceError=0,tensorValidationSeconds=0;
 dssolver::BatchResponseCache cache({static_cast<std::size_t>(2*degree)*grid.Count()+1,16ull<<30,1e-12});
 auto dispersion=[&](Complex z){return in.calculation.attenuation?2/pi*std::log(Complex(0,1)*z/model.ReferenceFrequency()):Complex{};};
 cacheAndResultSetup.stop();
 auto timer=Clock::now();
 inputAndGridSetup.stop();
#ifndef MATHPROFILE_DISABLED
 capturePhase("setup");
#endif
 mathprofile::Scope forwardPipeline("driver_initial_frequency_workspace_and_forward_capture");
 if(backend=="yspec"){
  if(epsilon!=0)throw std::runtime_error("YSpec perturbation not requested");omp_set_num_threads(workers);
  auto spectra=YSpec::SeismogramSpectra(model,in.calculation,std::span{&in.source,1},std::span{&in.receiver,1});
  for(std::size_t k=0;k<grid.Count();++k){components[k]={spectra(0,0,YSpec::Vertical,k),spectra(0,0,YSpec::North,k),spectra(0,0,YSpec::East,k)};raw[k]=dir.up*spectra(0,0,YSpec::Vertical,k)+dir.north*spectra(0,0,YSpec::North,k)+dir.east*spectra(0,0,YSpec::East,k);}
 }else if(backend=="dspec"||backend=="dspec-forward"){
  for(std::size_t k=0;k<grid.Count();++k){auto z=grid.Omega(k);auto workspace=SPARSESPEC::Windowed3DFieldEngine::prepareFrequencyWorkspace(sem,source,parameters,z,Eigen::Vector3cd::Zero(),&cache,true,false,nullptr,workers,dispersion(z));auto field=workspace.sampleForKernelsPacked(receiver,0,nullptr,2).front();components[k]={field.forwardRawDisplacement[0],field.forwardRawDisplacement[1],field.forwardRawDisplacement[2]};for(int c=0;c<3;++c)raw[k]+=direction[c]*field.forwardRawDisplacement[c];}
 }else throw std::runtime_error("unknown backend");
 forwardPipeline.stop();
#ifndef MATHPROFILE_DISABLED
 capturePhase("initial_frequency_workspace_and_capture");
#endif
 forwardSeconds=seconds(timer);timer=Clock::now();mathprofile::Scope measurementPipeline("driver_measurement_fft_window_adjoint_and_diagnostic_construction");auto measured=measure(in,model,grid,raw);measurementSeconds=seconds(timer);measurementPipeline.stop();
 mathprofile::Scope outputAndMeasurementIO("driver_initial_trace_spectrum_components_and_nodes_io");
 {std::ofstream f(out/"trace.csv");f<<"time_s,velocity_m_s,rate_m_s2,window,adjoint_source\n"<<std::setprecision(17);for(std::size_t n=0;n<measured.velocity.size();++n)f<<n*in.calculation.sampling.timeStep<<','<<measured.velocity[n]<<','<<measured.rate[n]<<','<<measured.window[n]<<','<<measured.measurement.adjointSource[n]<<'\n';}
 {std::ofstream f(out/"spectrum.csv");f<<"bin,frequency_hz,real_m,imag_m\n"<<std::setprecision(17);for(std::size_t k=0;k<raw.size();++k)f<<grid.First()+k<<','<<grid.Frequency(k)/scales.Time()<<','<<raw[k].real()*scales.Length()<<','<<raw[k].imag()*scales.Length()<<'\n';}
 {std::ofstream f(out/"components.csv");f<<"bin,up_re,up_im,north_re,north_im,east_re,east_im\n"<<std::setprecision(17);for(std::size_t k=0;k<raw.size();++k){f<<grid.First()+k;for(auto v:components[k])f<<','<<v.real()*scales.Length()<<','<<v.imag()*scales.Length();f<<'\n';}}
 {std::ofstream f(out/"nodes.csv");f<<"element,node,layer,radius_m,dr_m,solid\n"<<std::setprecision(17);for(auto n:nodes)f<<n.element<<','<<n.node<<','<<n.layer<<','<<n.radius<<','<<n.dr<<','<<n.solid<<'\n';}
 outputAndMeasurementIO.stop();
#ifndef MATHPROFILE_DISABLED
 capturePhase("measurement_and_initial_output");
#endif
 std::vector<std::vector<double>> values;for(auto &s:samplings)values.emplace_back(s.nodes.size()*s.places.size());
 timer=Clock::now();
 mathprofile::Scope kernelPipeline("driver_kernel_pipeline_including_diagnostics");
 if(backend=="yspec"){
  auto kernels=YSpec::ComputeSampledKernels(model,in.calculation,in.source,std::span{&measured.measurement,1},YSpec::ParameterSet::Velocities,YSpec::Scaling::Relative,samplings,[](std::size_t n,std::size_t total){std::cout<<"YSpec frequencies "<<n<<'/'<<total<<std::endl;});
  for(std::size_t s=0;s<samplings.size();++s){auto &k=kernels[0][s];auto names=k.Names();auto it=std::find(names.begin(),names.end(),"vs");if(it==names.end())throw std::runtime_error("missing vs");auto parameter=it-names.begin();for(std::size_t n=0;n<nodes.size();++n){auto v=k.Values(parameter,n);std::copy(v.begin(),v.end(),values[s].begin()+n*v.size());}}
 }else if(backend=="dspec"){
  mathprofile::Scope adjointAndPlanSetup("driver_kernel_adjoint_spectrum_and_plan_setup");auto adjoint=YSpec::MakeAdjointSpectrum(model,in.calculation,measured.measurement);
  std::vector<SPARSESPEC::Windowed3DFieldEngine::SamplingPlanCache> plans(samplings.size());adjointAndPlanSetup.stop();
  for(std::size_t k=0;k<grid.Count();++k){
   mathprofile::Scope frequencySetup("driver_kernel_frequency_weights_and_workspace");
   auto z=grid.Omega(k),d=dispersion(z);auto weights=direction*std::conj(adjoint.weight[k]*adjoint.force[k]);auto workspace=SPARSESPEC::Windowed3DFieldEngine::prepareFrequencyWorkspace(sem,source,parameters,z,weights,&cache,true,false,nullptr,workers,d);
   frequencySetup.stop();
   for(std::size_t s=0;s<samplings.size();++s){const auto &sampling=samplings[s];mathprofile::Scope pointsConstruction("driver_section_points_construction");std::vector<SPARSESPEC::Field3DPoint> points;for(auto p:sampling.places)for(auto n:nodes)points.push_back({n.radius,p.colatitude,p.longitude,n.element});pointsConstruction.stop();
   mathprofile::Scope tensorCall("driver_tensor_sampling_call");auto fields=workspace.sampleForKernelsTensor(points,points.size(),&plans[s],2,{},false);tensorCall.stop();
   if(k<2){
#ifndef MATHPROFILE_DISABLED
    const auto beforeDiagnosticSnapshot=Clock::now();auto diagnosticBefore=mathprofile::snapshot();profileSnapshotOverheadSeconds+=seconds(beforeDiagnosticSnapshot);
#endif
    auto validationStart=Clock::now();
    mathprofile::Scope tensorValidation("driver_tensor_vs_packed_diagnostic_validation");
    std::vector<SPARSESPEC::Field3DPoint> probe;std::vector<std::size_t> ids;
    for(std::size_t a=0;a<sampling.places.size();a+=std::max<std::size_t>(1,sampling.places.size()/7))
      for(std::size_t n=0;n<nodes.size();++n){ids.push_back(a*nodes.size()+n);probe.push_back(points[ids.back()]);}
    auto reference=workspace.sampleForKernelsPacked(probe,probe.size(),nullptr,2);
    auto array=[](auto v){return std::array<Complex,6>{v.rr,v.horizontalTrace,v.radialPlus,v.radialMinus,v.horizontalPlus,v.horizontalMinus};};
    double error=0,norm=0;
    for(std::size_t i=0;i<ids.size();++i)for(int adj=0;adj<2;++adj){
      auto x=array(adj?fields[ids[i]].adjointDualStrain:fields[ids[i]].forwardRawStrain);
      auto y=array(adj?reference[i].adjointDualStrain:reference[i].forwardRawStrain);
      for(int c=0;c<6;++c){error+=std::norm(x[c]-y[c]);norm+=std::norm(y[c]);}
    }
    tensorReferenceError=std::max(tensorReferenceError,std::sqrt(error/std::max(norm,1e-300)));
    if(tensorReferenceError>1e-11)throw std::runtime_error("tensor versus packed native-side regression failed");
    tensorValidation.stop();
    tensorValidationSeconds+=seconds(validationStart);
#ifndef MATHPROFILE_DISABLED
    const auto afterDiagnosticSnapshot=Clock::now();auto diagnosticAfter=mathprofile::snapshot();auto diagnosticDelta=mathprofile::difference(diagnosticAfter,diagnosticBefore);mathprofile::accumulate(diagnosticTotals,diagnosticDelta);profileSnapshotOverheadSeconds+=seconds(afterDiagnosticSnapshot);
#endif
   }
   mathprofile::Scope fusedMaterialContraction("kernel_fused_material_evaluation_weak_form_contraction_frequency_accumulation");
#ifndef MATHPROFILE_DISABLED
   const auto fluidNodes=static_cast<std::uint64_t>(std::count_if(nodes.begin(),nodes.end(),[](const Node &node){return !node.solid;}));
   mathprofile::count("kernel_point_samples",static_cast<std::uint64_t>(points.size()));
   mathprofile::count("kernel_fluid_point_samples",fluidNodes*static_cast<std::uint64_t>(sampling.places.size()));
#endif
   for(std::size_t a=0;a<sampling.places.size();++a)for(std::size_t n=0;n<nodes.size();++n)values[s][n*sampling.places.size()+a]+=kernelVs(sem,fields[a*nodes.size()+n],model(nodes[n].layer,nodes[n].radius/scales.Length()),d);
   fusedMaterialContraction.stop();}
   std::cout<<"DSpec frequencies "<<k+1<<'/'<<grid.Count()<<std::endl;
  }
 }
 kernelSeconds=seconds(timer)-tensorValidationSeconds;
 kernelPipeline.stop();
#ifndef MATHPROFILE_DISABLED
 capturePhase("kernel_production_including_separately_isolated_diagnostics");
#endif
 mathprofile::Scope kernelOutput("driver_kernel_output_and_result_io");
 for(std::size_t s=0;s<samplings.size();++s){std::string name=s==0?"along":"across";std::ofstream f(out/(name+".bin"),std::ios::binary);f.write(reinterpret_cast<char*>(values[s].data()),values[s].size()*sizeof(double));std::ofstream a(out/(name+"-angles.csv"));a<<"colatitude_rad,longitude_rad,along_deg,across_deg,angular_weight_sr\n"<<std::setprecision(17);for(std::size_t p=0;p<samplings[s].places.size();++p){auto x=samplings[s].places[p];a<<x.colatitude<<','<<x.longitude<<','<<x.along<<','<<x.across<<','<<(gridKind=="audit"?angularWeights[p]:0)<<'\n';}}
 double integrated=0;if(gridKind=="audit")for(std::size_t n=0;n<nodes.size();++n)for(std::size_t a=0;a<angularWeights.size();++a)integrated+=values[0][n*angularWeights.size()+a]*profile(nodes[n].radius)*nodes[n].radius*nodes[n].radius*nodes[n].dr*angularWeights[a];
 {std::ofstream f(out/"result.json");f<<std::setprecision(17)<<"{\"backend\":\""<<backend<<"\",\"degree\":"<<degree<<",\"nq\":"<<nq<<",\"refinement\":"<<refinement<<",\"workers\":"<<workers<<",\"attenuation\":"<<(in.calculation.attenuation?"true":"false")<<",\"epsilon_perturbation\":"<<epsilon<<",\"frequency_count\":"<<grid.Count()<<",\"first_bin\":"<<grid.First()<<",\"last_bin\":"<<grid.Last()<<",\"fft_samples\":"<<grid.Samples()<<",\"fft_period_s\":"<<grid.Samples()*in.calculation.sampling.timeStep<<",\"damping_si\":"<<grid.Damping()/scales.Time()<<",\"native_side_count\":"<<nodes.size()<<",\"along_angles\":"<<samplings[0].places.size()<<",\"across_angles\":"<<(samplings.size()>1?samplings[1].places.size():0)<<",\"gravity_max_abs_nd_error\":"<<gravityError<<",\"time_shift_adjoint_identity\":"<<measured.delayAudit<<",\"measurement_rate_norm\":"<<measured.norm<<",\"integrated_profile_kernel\":"<<integrated<<",\"forward_s\":"<<forwardSeconds<<",\"measurement_s\":"<<measurementSeconds<<",\"kernel_s\":"<<kernelSeconds<<",\"tensor_packed_relative_l2\":"<<tensorReferenceError<<",\"tensor_validation_s\":"<<tensorValidationSeconds<<",\"cache_bytes\":"<<cache.retainedBytes()<<",\"elapsed_s\":"<<seconds(start)<<",\"kernel_parameter\":\"ln_vs_at_fixed_vp_density_and_effective_Q\",\"kernel_units\":\"s m^-3\",\"layout\":\"native-side-major then angle, float64\",\"source_shell_contacts\":false}\n";}
 kernelOutput.stop();
#ifndef MATHPROFILE_DISABLED
 capturePhase("kernel_output_and_result_io");
 pipelineScope.stop();
 const auto finalSnapshotStart=Clock::now();auto finalProfileSnapshot=mathprofile::snapshot();profileSnapshotOverheadSeconds+=seconds(finalSnapshotStart);
 writeProfileAccounting(out/"profile.json",std::move(phaseSnapshots),diagnosticTotals,finalProfileSnapshot,profileSnapshotOverheadSeconds);
#endif
 return 0;
}catch(const std::exception &e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}
