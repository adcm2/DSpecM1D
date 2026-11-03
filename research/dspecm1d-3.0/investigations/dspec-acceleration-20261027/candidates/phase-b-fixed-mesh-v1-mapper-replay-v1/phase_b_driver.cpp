// Matched documentation-case adapter: common model/measurement, independent solvers.
#include <YSpec/YSpec.hpp>
#include "Spectra.hpp"
#include "KernelParts.hpp"
#include "Sphere.hpp"
#include "matched_mesh_model.hpp"
#include "phase_a_contract.hpp"
#include "phase_b_controls.hpp"
#include "engine_headers/Windowed3DKernels.h"
#include "frequency_scheduler.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
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
    double snapshotOverheadSeconds,bool frequencyRoute) {
 mathprofile::Scope reportScope("profile_snapshot_report_serialization_outside_pipeline");
 for(auto &[phase,delta]:phases)if(frequencyRoute&&phase=="kernel_frequency_parallel_pass_and_ordered_reduction_including_diagnostics"){
  for(const auto &[name,value]:diagnostics.entries){if(!name.starts_with("frequency_worker."))continue;auto found=delta.entries.find(name);if(found==delta.entries.end())std::terminate();
   if(found->second.inclusive<value.inclusive||found->second.exclusive<value.exclusive||found->second.calls<value.calls)std::terminate();
   found->second.inclusive-=value.inclusive;found->second.exclusive-=value.exclusive;found->second.calls-=value.calls;
   if(found->second.inclusive==0&&found->second.exclusive==0&&found->second.calls==0)delta.entries.erase(found);
  }
  for(const auto &[name,value]:diagnostics.counters){auto found=delta.counters.find(name);if(found==delta.counters.end()||found->second<value)std::terminate();found->second-=value;if(found->second==0)delta.counters.erase(found);}
 }
 for(auto &[phase,delta]:phases)if(!frequencyRoute&&phase=="kernel_production_including_separately_isolated_diagnostics"){
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
 for(const auto &[phase,delta]:phases){(void)phase;for(const auto &[name,value]:delta.entries)if(!name.starts_with("parallel_worker.")&&!name.starts_with("frequency_worker."))partition+=value.exclusive;}
 for(const auto &[name,value]:diagnostics.entries)if(!name.starts_with("parallel_worker.")&&!name.starts_with("frequency_worker."))partition+=value.exclusive;
 const double partitionSeconds=static_cast<double>(partition)*1e-9;
 const double residual=partitionSeconds-pipelineWall;
 const double tolerance=std::max(.01,.001*pipelineWall);
 if(!std::isfinite(pipelineWall)||!std::isfinite(partitionSeconds)||std::abs(residual)>tolerance)
  throw std::runtime_error("profile exclusive wall partition failed closure");
 auto component=[](const std::string &rawName){
  std::string name=rawName;
  if(name.starts_with("frequency_worker."))name.erase(0,std::string("frequency_worker.").size());
  else if(name.starts_with("parallel_worker."))name.erase(0,std::string("parallel_worker.").size());
  static const std::map<std::string,std::string> scopeComponents={
  {"batch_cache_locked_hit_response_return_copies","batch_response_cache"},
  {"batch_cache_locked_hit_signature_validation","batch_response_cache"},
  {"batch_cache_locked_initial_lookup","batch_response_cache"},
  {"batch_response_cache_get_or_compute","batch_response_cache"},
  {"batch_response_cache_miss_admission_signatures_and_copies","batch_response_cache"},
  {"driver_frequency_grid_setup","driver_setup"},
  {"driver_frequency_ordered_kernel_contribution_reduction","ordered_kernel_reduction"},
  {"driver_frequency_parallel_kernel_pass_including_diagnostics","kernel_pipeline_including_diagnostics"},
  {"driver_initial_frequency_workspace_and_forward_capture","driver_setup"},
  {"driver_initial_trace_spectrum_components_and_nodes_io","driver_setup"},
  {"driver_input_and_grid_setup","driver_setup"},
  {"driver_kernel_adjoint_spectrum_and_plan_setup","driver_setup"},
  {"driver_kernel_frequency_weights_and_workspace","driver_setup"},
  {"driver_kernel_output_and_result_io","driver_setup"},
  {"driver_kernel_pipeline_including_diagnostics","kernel_pipeline_including_diagnostics"},
  {"driver_measurement_fft_window_adjoint_and_diagnostic_construction","measurement_and_adjoint_construction"},
  {"driver_model_deck_and_common_material_setup","driver_setup"},
  {"driver_native_side_sampling_and_gravity_audit","driver_setup"},
  {"driver_parameters_file_io_and_input_parameters","driver_setup"},
  {"driver_result_cache_and_sampling_plan_initialisation","driver_setup"},
  {"driver_section_points_construction","driver_setup"},
  {"driver_section_sampling_and_receiver_construction","driver_setup"},
  {"driver_sem_constructor_mesh_static_matrices_and_model_setup","driver_setup"},
  {"driver_tensor_sampling_call","tensor_angular_synthesis"},
  {"driver_tensor_vs_packed_diagnostic_validation","diagnostic_validation"},
  {"dynamic_H_assembly","modal_solver"},
  {"forward_substitution","modal_solver"},
  {"frequency_workspace_prepare_and_modal_capture","unclassified_by_scope_name"},
  {"initial_frequency_workspace_forward_and_unit_adjoint_capture","frequency_initial_task"},
  {"kernel_frequency_workspace_sampling_validation_and_contribution","frequency_kernel_task"},
  {"kernel_fused_material_evaluation_weak_form_contraction_frequency_accumulation","material_evaluation_weak_form_contraction_accumulation"},
  {"minimal_unit_basis_HH_solve","modal_solver"},
  {"modal_degree_ordered_result_collection","parallel_ordered_result_collection"},
  {"modal_degree_parallel_region_wall","parallel_degree_region_wall"},
  {"modal_degree_worker_profile_merge","parallel_worker_profile_merge"},
  {"modal_projection_source_lift_and_frequency_solve","unclassified_by_scope_name"},
  {"modal_response_copy_and_capture","unclassified_by_scope_name"},
  {"packed_Wigner_d_table_build","tensor_angular_synthesis"},
  {"packed_angular_coefficients_and_basis_build","tensor_angular_synthesis"},
  {"packed_to_compact_kernel_sample_copy","tensor_angular_synthesis"},
  {"pole_frame_sampling_total","tensor_angular_synthesis"},
  {"profile_pipeline_wall","pipeline_total"},
  {"radial_interpolation_and_strain_numerator_build","radial_sampling_reconstruction"},
  {"radial_l0_sampling_and_interpolation","radial_sampling_reconstruction"},
  {"radial_l0_source_and_adjoint_branch","radial_sampling_reconstruction"},
  {"receiver_at_pole_closed_form","unclassified_by_scope_name"},
  {"receiver_basis_and_QR_reconstruction","receiver_basis_projection"},
  {"receiver_basis_mix_and_unit_adjoint_projection","receiver_basis_projection"},
  {"receiver_frame_geometry_phase_spin_rotation_and_scatter","tensor_angular_synthesis"},
  {"receiver_local_weight_materialization","unclassified_by_scope_name"},
  {"shared_forward_and_adjoint_solver","unclassified_by_scope_name"},
  {"solver_factorization_dispatch","modal_solver"},
  {"source_frame_geometry_phase_spin_rotation_and_scatter","tensor_angular_synthesis"},
  {"sparse_adjoint_eigen_solve","modal_solver"},
  {"sparse_adjoint_result_unscaling_and_validation","modal_solver"},
  {"sparse_adjoint_result_validation","modal_solver"},
  {"sparse_adjoint_rhs_scaling","modal_solver"},
  {"sparse_adjoint_rhs_validation","modal_solver"},
  {"sparse_adjoint_solve_including_rhs_scaling","modal_solver"},
  {"sparse_equilibration_and_scaling","modal_solver"},
  {"sparse_factorization_total","modal_solver"},
  {"sparse_forward_eigen_solve","modal_solver"},
  {"sparse_forward_result_unscaling_and_validation","modal_solver"},
  {"sparse_forward_result_validation","modal_solver"},
  {"sparse_forward_rhs_scaling","modal_solver"},
  {"sparse_forward_rhs_validation","modal_solver"},
  {"sparse_forward_solve_including_rhs_scaling","modal_solver"},
  {"sparse_numerical_factorization","modal_solver"},
  {"sparse_symbolic_analysis","modal_solver"},
  {"spheroidal_receiver_basis_assembly_and_qr","receiver_basis_projection"},
  {"spheroidal_source_lift_projection_and_solver","source_projection_solver_wrapper_remainder"},
  {"spheroidal_source_treatment_and_lift","source_treatment_and_lift"},
  {"tensor_angular_gemm_complex_and_packing","tensor_angular_synthesis"},
  {"tensor_angular_plan_key_build_and_lookup","tensor_angular_synthesis"},
  {"tensor_coefficient_layout_dedup_and_permutation_build","tensor_angular_synthesis"},
  {"tensor_coefficient_packing_and_key_reorder","tensor_angular_synthesis"},
  {"tensor_cold_frame_geometry_and_wigner_build","tensor_angular_synthesis"},
  {"tensor_complex_basis_cache_lookup","tensor_angular_synthesis"},
  {"tensor_complex_basis_cold_fill_and_cache","tensor_angular_synthesis"},
  {"tensor_complex_basis_key_construction","tensor_angular_synthesis"},
  {"tensor_grid_cold_layout_map_build","tensor_angular_synthesis"},
  {"tensor_grid_semantic_key_build_and_lookup","tensor_angular_synthesis"},
  {"tensor_phase_spin_rotation_and_scatter","tensor_angular_synthesis"},
  {"tensor_radial_l0_scatter","radial_sampling_reconstruction"},
  {"tensor_representative_and_coefficient_record_capture","tensor_angular_synthesis"},
  {"tensor_result_allocation_and_volume_tagging","tensor_angular_synthesis"},
  {"tensor_sampling_total","tensor_angular_synthesis"},
  {"toroidal_receiver_basis_assembly_and_qr","receiver_basis_projection"},
  {"toroidal_source_lift_projection_and_solver","source_projection_solver_wrapper_remainder"},
  {"toroidal_source_treatment_and_lift","source_treatment_and_lift"},
  {"workspace_synthesis_vector_construction","tensor_angular_synthesis"},
  };
  auto found=scopeComponents.find(name);
  return found==scopeComponents.end()?std::string("unclassified"):found->second;
 };
 std::ofstream out(path);if(!out)throw std::runtime_error("cannot create profile accounting output");
 out<<std::setprecision(17)<<"{\"schema\":1,\"snapshot_overhead_s\":"<<snapshotOverheadSeconds
    <<",\"report_serialization_scope\":\"profile_snapshot_report_serialization_outside_pipeline\""
    <<",\"frequency_route\":"<<(frequencyRoute?"true":"false")<<",\"diagnostic_worker_elapsed_is_not_subtracted_from_wall\":"<<(frequencyRoute?"true":"false")<<",\"worker_diagnostic_partition\":{\"production_rows_exclude_diagnostic_overlay\":"<<(frequencyRoute?"true":"false")<<",\"raw_inclusive_umbrellas_preserve_diagnostics\":"<<(frequencyRoute?"true":"false")<<",\"diagnostic_overlay_phase\":\"diagnostic_tensor_validation\"},\"wall_partition\":{\"profile_pipeline_wall_s\":"<<pipelineWall<<",\"master_exclusive_partition_s\":"<<partitionSeconds
    <<",\"residual_s\":"<<residual<<",\"tolerance_s\":"<<tolerance<<",\"worker_elapsed_excluded\":true},\"rows\":[";
 bool first=true;
 auto writeRows=[&](const std::string &phase,const mathprofile::Snapshot &snapshot,bool diagnostic){
  for(const auto &[name,value]:snapshot.entries){if(!first)out<<',';first=false;
   const bool worker=name.starts_with("parallel_worker.")||name.starts_with("frequency_worker.");
   std::string localName=name;if(localName.starts_with("frequency_worker."))localName.erase(0,std::string("frequency_worker.").size());else if(localName.starts_with("parallel_worker."))localName.erase(0,std::string("parallel_worker.").size());
   const bool pureDiagnosticScope=localName=="driver_tensor_vs_packed_diagnostic_validation";
   out<<"{\"phase\":\""<<phase<<"\",\"scope_name\":\""<<name<<"\",\"execution_domain\":\""<<(worker?"worker_elapsed_sum":"master_wall")
      <<"\",\"inclusive_s\":"<<value.inclusive*1e-9<<",\"exclusive_s\":"<<value.exclusive*1e-9<<",\"calls\":"<<value.calls
      <<",\"logical_component\":\""<<component(name)<<"\",\"diagnostic\":"<<((diagnostic||pureDiagnosticScope)?"true":"false")
      <<",\"parent_overlap_description\":\""<<(name.starts_with("frequency_worker.")?"frequency worker elapsed sum overlaps frequency-pass wall and is excluded from master wall partition":worker?"degree worker elapsed sum overlaps the degree-region wall bucket and is excluded from master wall partition":(name=="driver_kernel_pipeline_including_diagnostics"||name=="driver_frequency_parallel_kernel_pass_including_diagnostics")?"raw inclusive umbrella includes diagnostic work; exclusive value is retained as its master wall partition term":"exclusive time is partitioned; inclusive time overlaps children and is not additive")<<"\"}";
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
struct CommonComparisonPoint {
 int section=-1,layer=-1,referenceElement=-1,referenceNode=-1,sectionPoint=-1;
 double radiusMetres=0;bool fluid=false,hasForward=false,hasAdjoint=false;
 SPARSESPEC::Field3DPoint point{};
 std::array<Complex,6> forward{},adjoint{};
};
int commonSideElement(const Full1D::SEM &sem,const CommonModel &model,
                      const phasea::ReferenceFieldRecord &r,double tolerance) {
 const auto &mesh=sem.mesh();int chosen=-1;
 const double radius=r.radiusNormalized;
 const bool physicalBreak=std::abs(radius-sem.sourceRadius())<=tolerance||
     std::abs(radius-model.LowerRadius(r.layer))<=tolerance||
     std::abs(radius-model.UpperRadius(r.layer))<=tolerance;
 for(int e=0;e<mesh.NE();++e){
  if(mesh.LayerNumber(e)!=r.layer)continue;
  const double lo=mesh.ELR(e),hi=mesh.EUR(e);
  // Preserve a reference's explicit one-sided element at any coincident
  // boundary, including artificial boundaries in the old reference mesh.
  if(r.localNode==0&&std::abs(lo-radius)<=tolerance){chosen=e;break;}
  if(r.localNode==9&&std::abs(hi-radius)<=tolerance){chosen=e;break;}
 }
 // A noncoincident artificial reference boundary may lie inside a candidate
 // element. Fall back only within the same material; never cross a physical
 // source or material break.
 if(chosen<0&&!physicalBreak){
  for(int e=0;e<mesh.NE();++e){
   if(mesh.LayerNumber(e)!=r.layer)continue;
   if(radius>=mesh.ELR(e)-tolerance&&radius<=mesh.EUR(e)+tolerance){chosen=e;break;}
  }
 }
 if(chosen<0||sem.elementIsFluid(chosen)!=(r.fluid!=0))
  throw std::runtime_error("Phase B common-radius side has no same-material candidate element.");
 return chosen;
}
struct FrequencyWorkerContext {
 mathprofile::WorkerCapture profile;
 int previousBackend;
 FrequencyWorkerContext(bool profiling,mathprofile::WorkerData &data,int backend) noexcept
  :profile(profiling,data),previousBackend(SPARSESPEC::FixedMeshAdjointSolver::diagnosticBackend()) {
   try{SPARSESPEC::FixedMeshAdjointSolver::setDiagnosticBackend(backend);}catch(...){std::terminate();}
 }
 ~FrequencyWorkerContext() noexcept {try{SPARSESPEC::FixedMeshAdjointSolver::setDiagnosticBackend(previousBackend);}catch(...){std::terminate();}}
 FrequencyWorkerContext(const FrequencyWorkerContext&)=delete;
};
void mergeWorkerData(mathprofile::WorkerData &target,const mathprofile::WorkerData &source){
 for(const auto &[name,value]:source.entries){auto &entry=target.entries[name];entry.inclusive+=value.inclusive;entry.exclusive+=value.exclusive;entry.calls+=value.calls;}
 for(const auto &[name,value]:source.counters)target.counters[name]+=value;
}
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
 if(argc<10||argc>12)throw std::runtime_error("usage: phase-b-driver dspec INPUT EMPTY_OUTPUT L NQ REFINEMENT WORKERS sections|audit EPSILON [FROZEN_TRACE_OR_- [REFERENCE_PACKET_DIR]]");
 auto start=Clock::now();std::string backend=argv[1],gridKind=argv[8];fs::path out=fs::absolute(argv[3]);
 int degree=std::stoi(argv[4]),nq=std::stoi(argv[5]),workers=std::stoi(argv[7]);double refinement=std::stod(argv[6]),epsilon=std::stod(argv[9]);
 const char *frequencyEnv=std::getenv("DSPEC_FREQUENCY_WORKERS");const bool frequencyRoute=frequencyEnv!=nullptr;int frequencyWorkers=0;if(frequencyRoute){std::string setting(frequencyEnv);std::size_t consumed=0;frequencyWorkers=std::stoi(setting,&consumed);if(consumed!=setting.size())throw std::runtime_error("DSPEC_FREQUENCY_WORKERS must be a plain integer");}
 if(workers<1||workers>8||degree<1||degree>200||nq<5||nq>16||!(refinement>0)||std::abs(epsilon)>.001)throw std::runtime_error("bounded documentation scope exceeded");
 if(frequencyRoute&&(frequencyWorkers<1||frequencyWorkers>8||workers!=1||(backend!="dspec"&&backend!="dspec-forward")))throw std::runtime_error("frequency route requires dspec backend, CLI solve workers=1, and frequency workers in [1,8]");
 const int solveWorkers=frequencyRoute?1:workers;
#ifndef MATHPROFILE_DISABLED
 const bool workerProfiling=true;
#else
 const bool workerProfiling=false;
#endif
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
 {std::ofstream f(out/"frequency-grid.csv");if(!f)throw std::runtime_error("Cannot create Phase B exact frequency grid receipt.");f<<"bin,z_real,z_imag\n"<<std::setprecision(17);for(std::size_t k=0;k<grid.Count();++k){const auto z=grid.Omega(k);f<<grid.First()+static_cast<int>(k)<<','<<z.real()<<','<<z.imag()<<'\n';}}
 auto geoLat=[](double lat){return 90-YSpec::Sphere::Colatitude(lat)*180/pi;};const auto &moment=in.source.moment;
 SourceInfo::EarthquakeCMT source(geoLat(in.source.location.latitude),in.source.location.longitude,in.source.depth/1000,moment.rr,moment.rt,moment.rp,moment.tt,moment.tp,moment.pp);
 mathprofile::Scope parameterIO("driver_parameters_file_io_and_input_parameters");std::ofstream pf(out/"parameters.txt");pf<<params(in,degree);pf.close();InputParameters parameters((out/"parameters.txt").string());parameterIO.stop();
 const char *meshSelectionEnv=std::getenv("DSPEC_PHASE_B_MESH_SELECTION");
 phaseb::MeshSelection meshSelection;
 if(meshSelectionEnv&&std::string(meshSelectionEnv)=="1"){
  const char *fineEnv=std::getenv("DSPEC_PHASE_B_FINER");
  meshSelection=phaseb::selectMesh(nq,parameters.relative_error(),scales.Time(),
      grid.Omega(grid.Count()-1),fineEnv&&std::string(fineEnv)=="1");
  refinement=meshSelection.selectedStep;
 }
 mathprofile::Scope semSetup("driver_sem_constructor_mesh_static_matrices_and_model_setup");Full1D::SEM sem(common,refinement,nq,degree,source);semSetup.stop();
 mathprofile::Scope nativeSideAudit("driver_native_side_sampling_and_gravity_audit");
 {std::ofstream f(out/"mesh-elements.csv");if(!f)throw std::runtime_error("Cannot create Phase B mesh element-boundary receipt.");f<<"element,layer,lower_radius_m,upper_radius_m,fluid\n"<<std::setprecision(17);for(int e=0;e<sem.mesh().NE();++e)f<<e<<','<<sem.mesh().LayerNumber(e)<<','<<sem.mesh().ELR(e)*scales.Length()<<','<<sem.mesh().EUR(e)*scales.Length()<<','<<(sem.elementIsFluid(e)?1:0)<<'\n';}
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
 std::map<std::tuple<int,int,int,int,int>,CommonComparisonPoint> commonPointMap;
 std::map<int,std::vector<CommonComparisonPoint>> commonPointsByBin;
 if(argc>=12){
  for(const int bin:{9,25,40}){
   phasea::PacketHeader header;
   const auto records=phasea::readPacket<phasea::ReferenceFieldRecord>(
       fs::path(argv[11])/("reference-bin"+std::to_string(bin)+".bin"),header);
   if(header.bin!=bin||header.recordCount!=records.size()||records.size()!=2u*293u*(135u+19u)||
      header.lengthScaleMetres!=scales.Length())
    throw std::runtime_error("Phase B frozen common-field packet does not match the accepted NQ10 reference contract.");
   for(const auto &record:records){
    if(record.bin!=bin||record.section<0||record.section>1||record.sectionPoint<0||record.layer<0||
       !std::isfinite(record.radiusNormalized)||!std::isfinite(record.thetaRadians)||!std::isfinite(record.phiRadians))
      throw std::runtime_error("Phase B common-field reference row is invalid.");
    const auto key=std::tuple{bin,record.section,record.element,record.localNode,record.sectionPoint};
    auto [it,added]=commonPointMap.try_emplace(key);
    auto &entry=it->second;
    if(added){
     entry.section=record.section;entry.layer=record.layer;entry.referenceElement=record.element;
     entry.referenceNode=record.localNode;entry.sectionPoint=record.sectionPoint;
     entry.radiusMetres=record.radiusNormalized*scales.Length();entry.fluid=record.fluid!=0;
     const double tol=32*std::numeric_limits<double>::epsilon()*
         std::max({1.0,std::abs(record.radiusNormalized),std::abs(sem.mesh().PR())});
     const int element=commonSideElement(sem,common,record,tol);
     entry.point={entry.radiusMetres,record.thetaRadians,record.phiRadians,element};
    }else if(entry.layer!=record.layer||entry.fluid!=(record.fluid!=0)||
             entry.radiusMetres!=record.radiusNormalized*scales.Length())
     throw std::runtime_error("Phase B common-side metadata changed between packet roles.");
    auto &target=record.role==phasea::FieldRole::ForwardSource?entry.forward:entry.adjoint;
    auto &present=record.role==phasea::FieldRole::ForwardSource?entry.hasForward:entry.hasAdjoint;
    if(present)throw std::runtime_error("Phase B common reference has duplicate role rows.");
    present=true;
    for(std::size_t c=0;c<6;++c)target[c]=record.strain[c].value();
   }
  }
  for(const auto &[key,entry]:commonPointMap){(void)key;if(!entry.hasForward||!entry.hasAdjoint)throw std::runtime_error("Phase B common reference lacks one field role.");commonPointsByBin[std::get<0>(key)].push_back(entry);}
  if(commonPointsByBin.size()!=3||commonPointsByBin[9].size()!=293u*(135u+19u)||
     commonPointsByBin[25].size()!=293u*(135u+19u)||commonPointsByBin[40].size()!=293u*(135u+19u))
   throw std::runtime_error("Phase B common one-sided location count differs from the frozen 293-side grid.");
 }
 sectionAndReceiverSetup.stop();
 mathprofile::Scope cacheAndResultSetup("driver_result_cache_and_sampling_plan_initialisation");
 std::vector<Complex> raw(grid.Count());std::vector<std::array<Complex,3>> components(grid.Count());double forwardSeconds=0,kernelSeconds=0,measurementSeconds=0,tensorReferenceError=0,tensorValidationSeconds=0;
 double kernelParallelPassSeconds=0,reductionSeconds=0,diagnosticWorkerElapsedSeconds=0;std::size_t contributionBufferBytes=0;
 dssolver::BatchResponseCache cache({static_cast<std::size_t>(2*degree)*grid.Count()+1,16ull<<30,1e-12});
 auto dispersion=[&](Complex z){return in.calculation.attenuation?2/pi*std::log(Complex(0,1)*z/model.ReferenceFrequency()):Complex{};};
 cacheAndResultSetup.stop();
 auto timer=Clock::now();
 inputAndGridSetup.stop();
#ifndef MATHPROFILE_DISABLED
 capturePhase("setup");
#endif
 const char *phaseBRestrictionEnv=std::getenv("DSPEC_PHASE_B_RESTRICTED");
 const bool phaseBRestrictedRequested=phaseBRestrictionEnv&&std::string(phaseBRestrictionEnv)=="1";
 const char *replayEnv=std::getenv("DSPEC_PHASE_B_DIAGNOSTIC_REPLAY");
 const bool diagnosticReplay=replayEnv&&std::string(replayEnv)=="1";
 if(diagnosticReplay){
  const char *directEnv=std::getenv("DSPEC_PHASE_C_DIRECT_ADJOINT");
  if(argc<12||frequencyRoute||workers!=1||(directEnv&&std::string(directEnv)=="1"))
   throw std::runtime_error("B frozen replay requires one worker, packet+trace inputs, and the original rank-basis adjoint route.");
  std::ifstream frozen(argv[10]);std::string line;
  if(!frozen||!std::getline(frozen,line)||line!="time_s,velocity_m_s,rate_m_s2,window,adjoint_source")
   throw std::runtime_error("B replay cannot read the hash-bound frozen waveform adjoint source.");
  YSpec::Measurement replayMeasurement{in.receiver,dir,{}};
  while(std::getline(frozen,line)){
   std::stringstream row(line);std::string field;std::array<std::string,5> columns;
   for(auto &column:columns)if(!std::getline(row,column,','))throw std::runtime_error("B replay frozen trace has an invalid row.");
   if(std::getline(row,field,','))throw std::runtime_error("B replay frozen trace has trailing columns.");
   replayMeasurement.adjointSource.push_back(std::stod(columns[4]));
  }
  if(replayMeasurement.adjointSource.size()!=YSpec::SampleCount(grid,in.calculation.sampling))
   throw std::runtime_error("B replay frozen adjoint-source sample count differs from the case grid.");
  const auto adjoint=YSpec::MakeAdjointSpectrum(model,in.calculation,replayMeasurement);
  mathprofile::Scope replayScope("phase_b_frozen_three_bin_rank_basis_diagnostic_replay");
  for(const int bin:{9,25,40}){
   const auto index=static_cast<std::size_t>(bin-grid.First());
   if(index>=grid.Count())throw std::runtime_error("B replay target bin is outside the retained full mesh-selection grid.");
   const auto z=grid.Omega(index),d=dispersion(z);
   const auto weights=direction*std::conj(adjoint.weight[index]*adjoint.force[index]);
   auto workspace=SPARSESPEC::Windowed3DFieldEngine::prepareFrequencyWorkspace(
       sem,source,parameters,z,weights,nullptr,true,false,nullptr,solveWorkers,d);
   const auto receiverControl=workspace.sampleForKernelsPacked(receiver,0,nullptr,2).front();
   {std::ofstream receiverCsv(out/("kernel-workspace-receiver-bin"+std::to_string(bin)+".csv"));
    if(!receiverCsv)throw std::runtime_error("Cannot create B replay receiver spectrum.");
    receiverCsv<<"bin,up_re,up_im,north_re,north_im,east_re,east_im\n"<<std::setprecision(17)<<bin;
    for(int c=0;c<3;++c)receiverCsv<<','<<receiverControl.forwardRawDisplacement[c].real()*scales.Length()<<','<<receiverControl.forwardRawDisplacement[c].imag()*scales.Length();receiverCsv<<'\n';}
   const auto &entries=commonPointsByBin.at(bin);std::vector<SPARSESPEC::Field3DPoint> commonPoints;commonPoints.reserve(entries.size());
   for(const auto &entry:entries)commonPoints.push_back(entry.point);
   const auto commonFields=workspace.sampleForKernelsPacked(commonPoints,commonPoints.size(),nullptr,2);
   std::ofstream commonCsv(out/("common-fields-bin"+std::to_string(bin)+".csv"));
   if(!commonCsv)throw std::runtime_error("Cannot create B replay common-side diagnostics.");
   commonCsv<<"bin,section,reference_element,reference_node,section_point,layer,radius_m,candidate_element,fluid,location_class,kernel_vs";
   for(int c=0;c<6;++c)commonCsv<<",forward_ref"<<c<<"_re,forward_ref"<<c<<"_im,forward_candidate"<<c<<"_re,forward_candidate"<<c<<"_im,adjoint_ref"<<c<<"_re,adjoint_ref"<<c<<"_im,adjoint_candidate"<<c<<"_re,adjoint_candidate"<<c<<"_im";
   commonCsv<<'\n'<<std::setprecision(17);
   const double sourceRadius=scales.Length()*model.SurfaceRadius()-in.source.depth;
   const double receiverRadiusMetres=scales.Length()*model.SurfaceRadius()-in.receiver.depth;
   for(std::size_t i=0;i<entries.size();++i){
    const auto &entry=entries[i];const auto &field=commonFields[i];
    const std::array<Complex,6> f{field.forwardRawStrain.rr,field.forwardRawStrain.horizontalTrace,field.forwardRawStrain.radialPlus,field.forwardRawStrain.radialMinus,field.forwardRawStrain.horizontalPlus,field.forwardRawStrain.horizontalMinus};
    const std::array<Complex,6> a{field.adjointDualStrain.rr,field.adjointDualStrain.horizontalTrace,field.adjointDualStrain.radialPlus,field.adjointDualStrain.radialMinus,field.adjointDualStrain.horizontalPlus,field.adjointDualStrain.horizontalMinus};
    const auto material=model(entry.layer,entry.radiusMetres/scales.Length());const double product=kernelVs(sem,field,material,d);
    const auto &layer=model.GetLayer(entry.layer);const bool interface=std::abs(entry.radiusMetres/scales.Length()-layer.Lower())<1e-11||std::abs(entry.radiusMetres/scales.Length()-layer.Upper())<1e-11;
    const char *kind=std::abs(entry.radiusMetres-sourceRadius)<1e-3?"source":std::abs(entry.radiusMetres-receiverRadiusMetres)<1e-3?"receiver":interface?"material_interface":"bulk";
    commonCsv<<bin<<','<<entry.section<<','<<entry.referenceElement<<','<<entry.referenceNode<<','<<entry.sectionPoint<<','<<entry.layer<<','<<entry.radiusMetres<<','<<entry.point.radialElement<<','<<(field.fluid?1:0)<<','<<kind<<','<<product;
    for(int c=0;c<6;++c)commonCsv<<','<<entry.forward[c].real()<<','<<entry.forward[c].imag()<<','<<f[c].real()<<','<<f[c].imag()<<','<<entry.adjoint[c].real()<<','<<entry.adjoint[c].imag()<<','<<a[c].real()<<','<<a[c].imag();commonCsv<<'\n';
   }
  }
  replayScope.stop();
  {std::ofstream result(out/"phase-b-diagnostic-replay.json");if(!result)throw std::runtime_error("Cannot write B replay metadata.");
   result<<std::setprecision(17)<<"{\"status\":\"replay_only\",\"waveform_solve_performed\":false,\"weights\":\"hash-bound frozen accepted adjoint source\",\"degree\":"<<degree<<",\"nq\":"<<nq<<",\"refinement\":"<<refinement<<",\"relative_error\":"<<parameters.relative_error()<<",\"mesh_selection_zmax_real\":"<<meshSelection.highestComplexOmega.real()<<",\"mesh_selection_zmax_imag\":"<<meshSelection.highestComplexOmega.imag()<<",\"element_count\":"<<sem.mesh().NE()<<",\"source_element\":"<<sem.sourceElement(source)<<",\"receiver_element\":"<<receiverElement<<",\"restricted\":"<<(phaseBRestrictedRequested?"true":"false")<<",\"bins\":[9,25,40],\"radial_dofs_R\":"<<sem.hR().rows()<<",\"radial_dofs_S\":"<<sem.hS(degree).rows()<<",\"radial_dofs_T\":"<<sem.hTk(degree).rows()<<"}\n";}
  pipelineScope.stop();
  const auto snapshot=mathprofile::snapshot();
  {std::ofstream profileFile(out/"profile.json");if(!profileFile)throw std::runtime_error("Cannot write B replay profile receipt.");
   profileFile<<"{\"schema\":1,\"rows\":[";bool first=true;
   for(const auto &[name,value]:snapshot.entries){if(!first)profileFile<<',';first=false;profileFile<<"{\"phase\":\"diagnostic_replay\",\"scope_name\":\""<<name<<"\",\"calls\":"<<value.calls<<",\"inclusive_s\":"<<value.inclusive*1e-9<<",\"exclusive_s\":"<<value.exclusive*1e-9<<"}";}
   profileFile<<"],\"counters\":[";first=true;
   for(const auto &[name,value]:snapshot.counters){if(!first)profileFile<<',';first=false;profileFile<<"{\"phase\":\"diagnostic_replay\",\"name\":\""<<name<<"\",\"value\":"<<value<<"}";}
   profileFile<<"]}\n";}
  return 0;
 }
 // The paired control keeps the complete 35-bin forward waveform untruncated;
 // enable the scoped operator restriction only after that solve and measurement.
 if(phaseBRestrictedRequested&&::setenv("DSPEC_PHASE_B_RESTRICTED","0",1)!=0)
  throw std::runtime_error("Cannot disable Phase B restriction for the full-domain waveform control.");
 mathprofile::Scope forwardPipeline("driver_initial_frequency_workspace_and_forward_capture");
 if(backend=="yspec"){
  if(epsilon!=0)throw std::runtime_error("YSpec perturbation not requested");omp_set_num_threads(workers);
  auto spectra=YSpec::SeismogramSpectra(model,in.calculation,std::span{&in.source,1},std::span{&in.receiver,1});
  for(std::size_t k=0;k<grid.Count();++k){components[k]={spectra(0,0,YSpec::Vertical,k),spectra(0,0,YSpec::North,k),spectra(0,0,YSpec::East,k)};raw[k]=dir.up*spectra(0,0,YSpec::Vertical,k)+dir.north*spectra(0,0,YSpec::North,k)+dir.east*spectra(0,0,YSpec::East,k);}
 }else if(backend=="dspec"||backend=="dspec-forward"){
  if(frequencyRoute){
   std::vector<mathprofile::WorkerData> workerData(static_cast<std::size_t>(frequencyWorkers));
   const int diagnosticBackend=SPARSESPEC::FixedMeshAdjointSolver::diagnosticBackend();
   dspec::frequency::runIndexedWithContext(grid.Count(),frequencyWorkers,grid.First(),"initial",
    [&](int worker) noexcept {return FrequencyWorkerContext(workerProfiling,workerData[static_cast<std::size_t>(worker)],diagnosticBackend);},
    [&](std::size_t k,auto &){mathprofile::Scope taskScope("initial_frequency_workspace_forward_and_unit_adjoint_capture");
      auto z=grid.Omega(k);auto workspace=SPARSESPEC::Windowed3DFieldEngine::prepareFrequencyWorkspace(sem,source,parameters,z,Eigen::Vector3cd::Zero(),&cache,true,false,nullptr,solveWorkers,dispersion(z));
      auto field=workspace.sampleForKernelsPacked(receiver,0,nullptr,2).front();components[k]={field.forwardRawDisplacement[0],field.forwardRawDisplacement[1],field.forwardRawDisplacement[2]};for(int c=0;c<3;++c)raw[k]+=direction[c]*field.forwardRawDisplacement[c];
    });
   for(const auto &data:workerData)mathprofile::mergeWorker(data,"frequency_worker.");
  }else for(std::size_t k=0;k<grid.Count();++k){auto z=grid.Omega(k);auto workspace=SPARSESPEC::Windowed3DFieldEngine::prepareFrequencyWorkspace(sem,source,parameters,z,Eigen::Vector3cd::Zero(),&cache,true,false,nullptr,solveWorkers,dispersion(z));auto field=workspace.sampleForKernelsPacked(receiver,0,nullptr,2).front();components[k]={field.forwardRawDisplacement[0],field.forwardRawDisplacement[1],field.forwardRawDisplacement[2]};for(int c=0;c<3;++c)raw[k]+=direction[c]*field.forwardRawDisplacement[c];}
 }else throw std::runtime_error("unknown backend");
 forwardPipeline.stop();
 if(phaseBRestrictedRequested&&::setenv("DSPEC_PHASE_B_RESTRICTED","1",1)!=0)
  throw std::runtime_error("Cannot enable Phase B restriction for kernel solves.");
#ifndef MATHPROFILE_DISABLED
 capturePhase("initial_frequency_workspace_and_capture");
#endif
 forwardSeconds=seconds(timer);timer=Clock::now();mathprofile::Scope measurementPipeline("driver_measurement_fft_window_adjoint_and_diagnostic_construction");auto measured=measure(in,model,grid,raw);measurementSeconds=seconds(timer);measurementPipeline.stop();
 YSpec::Measurement kernelMeasurement=measured.measurement;
 bool frozenKernelWeights=false;
 if(argc>=11&&std::string(argv[10])!="-"){
  std::ifstream frozen(argv[10]);std::string line;
  if(!frozen||!std::getline(frozen,line)||line!="time_s,velocity_m_s,rate_m_s2,window,adjoint_source")
    throw std::runtime_error("Phase B cannot read frozen accepted adjoint-source trace");
  kernelMeasurement.adjointSource.clear();
  while(std::getline(frozen,line)){
    std::stringstream row(line);std::string field;std::array<std::string,5> columns;
    for(auto &column:columns)if(!std::getline(row,column,','))throw std::runtime_error("Phase B frozen trace has an invalid row");
    if(std::getline(row,field,','))throw std::runtime_error("Phase B frozen trace has trailing columns");
    kernelMeasurement.adjointSource.push_back(std::stod(columns[4]));
  }
  if(kernelMeasurement.adjointSource.size()!=YSpec::SampleCount(grid,in.calculation.sampling))
    throw std::runtime_error("Phase B frozen adjoint-source sample count mismatch");
  frozenKernelWeights=true;
 }
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
  auto kernels=YSpec::ComputeSampledKernels(model,in.calculation,in.source,std::span{&kernelMeasurement,1},YSpec::ParameterSet::Velocities,YSpec::Scaling::Relative,samplings,[](std::size_t n,std::size_t total){std::cout<<"YSpec frequencies "<<n<<'/'<<total<<std::endl;});
  for(std::size_t s=0;s<samplings.size();++s){auto &k=kernels[0][s];auto names=k.Names();auto it=std::find(names.begin(),names.end(),"vs");if(it==names.end())throw std::runtime_error("missing vs");auto parameter=it-names.begin();for(std::size_t n=0;n<nodes.size();++n){auto v=k.Values(parameter,n);std::copy(v.begin(),v.end(),values[s].begin()+n*v.size());}}
 }else if(backend=="dspec"){
  mathprofile::Scope adjointAndPlanSetup("driver_kernel_adjoint_spectrum_and_plan_setup");auto adjoint=YSpec::MakeAdjointSpectrum(model,in.calculation,kernelMeasurement);
  if(frequencyRoute){
   using Plan=SPARSESPEC::Windowed3DFieldEngine::SamplingPlanCache;
   std::vector<std::vector<Plan>> plans;plans.reserve(static_cast<std::size_t>(frequencyWorkers));
   for(int w=0;w<frequencyWorkers;++w){std::vector<Plan> sectionPlans;sectionPlans.reserve(samplings.size());for(std::size_t q=0;q<samplings.size();++q)sectionPlans.emplace_back();plans.push_back(std::move(sectionPlans));}
   adjointAndPlanSetup.stop();
   std::size_t cellsPerFrequency=0;
   for(const auto &sampling:samplings){if(nodes.size()&&sampling.places.size()>std::numeric_limits<std::size_t>::max()/nodes.size())throw std::overflow_error("frequency contribution dimensions overflow");const auto cells=nodes.size()*sampling.places.size();if(cells>std::numeric_limits<std::size_t>::max()-cellsPerFrequency)throw std::overflow_error("frequency contribution dimensions overflow");cellsPerFrequency+=cells;}
   if(cellsPerFrequency&&grid.Count()>std::numeric_limits<std::size_t>::max()/cellsPerFrequency)throw std::overflow_error("frequency contribution count overflow");
   const std::size_t contributionValues=grid.Count()*cellsPerFrequency;
   if(contributionValues>std::numeric_limits<std::size_t>::max()/sizeof(double))throw std::overflow_error("frequency contribution byte count overflow");
   contributionBufferBytes=contributionValues*sizeof(double);
   constexpr std::size_t contributionBufferAllowanceBytes=512ull<<20;
   if(contributionBufferBytes>contributionBufferAllowanceBytes)
    throw std::runtime_error("frequency contribution buffer exceeds the Phase B 512 MiB reservation within the existing 16 GiB solve cap");
   std::vector<std::vector<std::vector<double>>> contributions(grid.Count());
   for(std::size_t k=0;k<grid.Count();++k){contributions[k].reserve(samplings.size());for(const auto &sampling:samplings)contributions[k].emplace_back(nodes.size()*sampling.places.size(),0.0);}
   std::vector<double> validationErrorByBin(grid.Count(),0.0),validationSecondsByBin(grid.Count(),0.0);
#ifndef MATHPROFILE_DISABLED
   std::vector<mathprofile::Snapshot> diagnosticByBin(grid.Count());std::vector<double> snapshotOverheadByBin(grid.Count(),0.0);
#endif
   std::vector<mathprofile::WorkerData> workerData(static_cast<std::size_t>(frequencyWorkers));
   const int diagnosticBackend=SPARSESPEC::FixedMeshAdjointSolver::diagnosticBackend();
   const auto passStart=Clock::now();mathprofile::Scope parallelPassScope("driver_frequency_parallel_kernel_pass_including_diagnostics");
   dspec::frequency::runIndexedWithContext(grid.Count(),frequencyWorkers,grid.First(),"kernel",
    [&](int worker) noexcept {return FrequencyWorkerContext(workerProfiling,workerData[static_cast<std::size_t>(worker)],diagnosticBackend);},
    [&](std::size_t k,auto &){
     mathprofile::Scope workerKernel("kernel_frequency_workspace_sampling_validation_and_contribution");
     mathprofile::Scope frequencySetup("driver_kernel_frequency_weights_and_workspace");
     auto z=grid.Omega(k),d=dispersion(z);auto weights=direction*std::conj(adjoint.weight[k]*adjoint.force[k]);
     auto workspace=SPARSESPEC::Windowed3DFieldEngine::prepareFrequencyWorkspace(sem,source,parameters,z,weights,phaseBRestrictedRequested?nullptr:&cache,true,false,nullptr,solveWorkers,d);
     const int bin=grid.First()+static_cast<int>(k);
     const auto commonIt=commonPointsByBin.find(bin);
     if(commonIt!=commonPointsByBin.end()){
      mathprofile::Scope commonSideEvaluation("phase_b_direct_common_side_complex_fields_and_products");
      const auto receiverControl=workspace.sampleForKernelsPacked(receiver,0,nullptr,2).front();
      std::ofstream receiverCsv(out/("kernel-workspace-receiver-bin"+std::to_string(bin)+".csv"));
      if(!receiverCsv)throw std::runtime_error("Cannot create Phase B kernel-workspace receiver spectrum output.");
      receiverCsv<<"bin,up_re,up_im,north_re,north_im,east_re,east_im\n"<<std::setprecision(17)<<bin;
      for(int c=0;c<3;++c)receiverCsv<<','<<receiverControl.forwardRawDisplacement[c].real()*scales.Length()<<','<<receiverControl.forwardRawDisplacement[c].imag()*scales.Length();
      receiverCsv<<'\n';
      std::vector<SPARSESPEC::Field3DPoint> commonPoints;commonPoints.reserve(commonIt->second.size());
      for(const auto &entry:commonIt->second)commonPoints.push_back(entry.point);
      // Reference packets deliberately retain duplicate angular coordinates on
      // opposite material/source sides. Packed evaluation accepts these locations
      // independently and still evaluates the candidate radial polynomial and
      // its derivative at each explicit radius.
      const auto commonFields=workspace.sampleForKernelsPacked(
          commonPoints,commonPoints.size(),nullptr,2);
      std::ofstream commonCsv(out/("common-fields-bin"+std::to_string(bin)+".csv"));
      if(!commonCsv)throw std::runtime_error("Cannot create Phase B common-side comparison output.");
      commonCsv<<"bin,section,reference_element,reference_node,section_point,layer,radius_m,candidate_element,fluid,location_class,kernel_vs";
      for(int c=0;c<6;++c)commonCsv<<",forward_ref"<<c<<"_re,forward_ref"<<c<<"_im,forward_candidate"<<c<<"_re,forward_candidate"<<c<<"_im,adjoint_ref"<<c<<"_re,adjoint_ref"<<c<<"_im,adjoint_candidate"<<c<<"_re,adjoint_candidate"<<c<<"_im";
      commonCsv<<'\n'<<std::setprecision(17);
      const double sourceRadius=scales.Length()*model.SurfaceRadius()-in.source.depth;
      const double receiverRadiusMetres=scales.Length()*model.SurfaceRadius()-in.receiver.depth;
      for(std::size_t i=0;i<commonIt->second.size();++i){
       const auto &entry=commonIt->second[i];const auto &field=commonFields[i];
       const std::array<Complex,6> f{field.forwardRawStrain.rr,field.forwardRawStrain.horizontalTrace,
        field.forwardRawStrain.radialPlus,field.forwardRawStrain.radialMinus,
        field.forwardRawStrain.horizontalPlus,field.forwardRawStrain.horizontalMinus};
       const std::array<Complex,6> a{field.adjointDualStrain.rr,field.adjointDualStrain.horizontalTrace,
        field.adjointDualStrain.radialPlus,field.adjointDualStrain.radialMinus,
        field.adjointDualStrain.horizontalPlus,field.adjointDualStrain.horizontalMinus};
       const auto material=model(entry.layer,entry.radiusMetres/scales.Length());
       const double product=kernelVs(sem,field,material,d);
       const auto &layer=model.GetLayer(entry.layer);
       const bool interface=std::abs(entry.radiusMetres/scales.Length()-layer.Lower())<1e-11||
                            std::abs(entry.radiusMetres/scales.Length()-layer.Upper())<1e-11;
       const char *kind=std::abs(entry.radiusMetres-sourceRadius)<1e-3?"source":
                        std::abs(entry.radiusMetres-receiverRadiusMetres)<1e-3?"receiver":
                        interface?"material_interface":"bulk";
       commonCsv<<bin<<','<<entry.section<<','<<entry.referenceElement<<','<<entry.referenceNode<<','
        <<entry.sectionPoint<<','<<entry.layer<<','<<entry.radiusMetres<<','<<entry.point.radialElement<<','
        <<(field.fluid?1:0)<<','<<kind<<','<<product;
       for(int c=0;c<6;++c)commonCsv<<','<<entry.forward[c].real()<<','<<entry.forward[c].imag()<<','
        <<f[c].real()<<','<<f[c].imag()<<','<<entry.adjoint[c].real()<<','<<entry.adjoint[c].imag()<<','
        <<a[c].real()<<','<<a[c].imag();
       commonCsv<<'\n';
      }
     }
     frequencySetup.stop();mathprofile::Snapshot diagnosticForBin;
     for(std::size_t section=0;section<samplings.size();++section){const auto &sampling=samplings[section];
      mathprofile::Scope pointsConstruction("driver_section_points_construction");std::vector<SPARSESPEC::Field3DPoint> points;for(auto point:sampling.places)for(auto node:nodes)points.push_back({node.radius,point.colatitude,point.longitude,node.element});pointsConstruction.stop();
      mathprofile::Scope tensorCall("driver_tensor_sampling_call");auto fields=workspace.sampleForKernelsTensor(points,points.size(),&plans.at(static_cast<std::size_t>(omp_get_thread_num())).at(section),2,{},false);tensorCall.stop();
      if(k<2){
#ifndef MATHPROFILE_DISABLED
       const auto beforeDiagnosticSnapshot=Clock::now();auto diagnosticBefore=mathprofile::snapshot();snapshotOverheadByBin[k]+=seconds(beforeDiagnosticSnapshot);
#endif
       const auto validationStart=Clock::now();mathprofile::Scope tensorValidation("driver_tensor_vs_packed_diagnostic_validation");
       std::vector<SPARSESPEC::Field3DPoint> probe;std::vector<std::size_t> ids;
       for(std::size_t a=0;a<sampling.places.size();a+=std::max<std::size_t>(1,sampling.places.size()/7))for(std::size_t n=0;n<nodes.size();++n){ids.push_back(a*nodes.size()+n);probe.push_back(points[ids.back()]);}
       auto reference=workspace.sampleForKernelsPacked(probe,probe.size(),nullptr,2);
       auto array=[](auto v){return std::array<Complex,6>{v.rr,v.horizontalTrace,v.radialPlus,v.radialMinus,v.horizontalPlus,v.horizontalMinus};};double error=0,norm=0;
       for(std::size_t i=0;i<ids.size();++i)for(int a=0;a<2;++a){auto x=array(a?fields[ids[i]].adjointDualStrain:fields[ids[i]].forwardRawStrain);auto y=array(a?reference[i].adjointDualStrain:reference[i].forwardRawStrain);for(int c=0;c<6;++c){error+=std::norm(x[c]-y[c]);norm+=std::norm(y[c]);}}
       validationErrorByBin[k]=std::max(validationErrorByBin[k],std::sqrt(error/std::max(norm,1e-300)));
       if(validationErrorByBin[k]>1e-11)throw std::runtime_error("tensor versus packed native-side regression failed");
       tensorValidation.stop();validationSecondsByBin[k]+=seconds(validationStart);
#ifndef MATHPROFILE_DISABLED
       const auto afterDiagnosticSnapshot=Clock::now();auto diagnosticAfter=mathprofile::snapshot();auto diagnosticDelta=mathprofile::difference(diagnosticAfter,diagnosticBefore);mathprofile::accumulate(diagnosticForBin,diagnosticDelta);snapshotOverheadByBin[k]+=seconds(afterDiagnosticSnapshot);
#endif
      }
      mathprofile::Scope fusedMaterialContraction("kernel_fused_material_evaluation_weak_form_contraction_frequency_accumulation");
#ifndef MATHPROFILE_DISABLED
      const auto fluidNodes=static_cast<std::uint64_t>(std::count_if(nodes.begin(),nodes.end(),[](const Node &node){return !node.solid;}));
      mathprofile::count("kernel_point_samples",static_cast<std::uint64_t>(points.size()));mathprofile::count("kernel_fluid_point_samples",fluidNodes*static_cast<std::uint64_t>(sampling.places.size()));
#endif
      auto &contribution=contributions[k][section];for(std::size_t a=0;a<sampling.places.size();++a)for(std::size_t n=0;n<nodes.size();++n)contribution[n*sampling.places.size()+a]=kernelVs(sem,fields[a*nodes.size()+n],model(nodes[n].layer,nodes[n].radius/scales.Length()),d);
      fusedMaterialContraction.stop();
     }
#ifndef MATHPROFILE_DISABLED
     diagnosticByBin[k]=std::move(diagnosticForBin);
#endif
    });
   parallelPassScope.stop();kernelParallelPassSeconds=seconds(passStart);
   for(const auto &data:workerData)mathprofile::mergeWorker(data,"frequency_worker.");
   for(std::size_t k=0;k<grid.Count();++k){tensorReferenceError=std::max(tensorReferenceError,validationErrorByBin[k]);tensorValidationSeconds+=validationSecondsByBin[k];diagnosticWorkerElapsedSeconds+=validationSecondsByBin[k];}
#ifndef MATHPROFILE_DISABLED
   for(std::size_t k=0;k<grid.Count();++k){mathprofile::Snapshot prefixed;for(const auto &[name,value]:diagnosticByBin[k].entries){auto &entry=prefixed.entries["frequency_worker."+name];entry.inclusive+=value.inclusive;entry.exclusive+=value.exclusive;entry.calls+=value.calls;}for(const auto &[name,value]:diagnosticByBin[k].counters)prefixed.counters[name]+=value;mathprofile::accumulate(diagnosticTotals,prefixed);profileSnapshotOverheadSeconds+=snapshotOverheadByBin[k];}
#endif
   const auto reductionStart=Clock::now();mathprofile::Scope reductionScope("driver_frequency_ordered_kernel_contribution_reduction");
   for(std::size_t k=0;k<grid.Count();++k)for(std::size_t section=0;section<samplings.size();++section){auto &destination=values[section];const auto &contribution=contributions[k][section];for(std::size_t a=0;a<samplings[section].places.size();++a)for(std::size_t n=0;n<nodes.size();++n)destination[n*samplings[section].places.size()+a]+=contribution[n*samplings[section].places.size()+a];}
   reductionScope.stop();reductionSeconds=seconds(reductionStart);
  }else{
   std::vector<SPARSESPEC::Windowed3DFieldEngine::SamplingPlanCache> plans(samplings.size());adjointAndPlanSetup.stop();
   for(std::size_t k=0;k<grid.Count();++k){mathprofile::Scope frequencySetup("driver_kernel_frequency_weights_and_workspace");auto z=grid.Omega(k),d=dispersion(z);auto weights=direction*std::conj(adjoint.weight[k]*adjoint.force[k]);auto workspace=SPARSESPEC::Windowed3DFieldEngine::prepareFrequencyWorkspace(sem,source,parameters,z,weights,phaseBRestrictedRequested?nullptr:&cache,true,false,nullptr,solveWorkers,d);frequencySetup.stop();
    for(std::size_t section=0;section<samplings.size();++section){const auto &sampling=samplings[section];mathprofile::Scope pointsConstruction("driver_section_points_construction");std::vector<SPARSESPEC::Field3DPoint> points;for(auto point:sampling.places)for(auto node:nodes)points.push_back({node.radius,point.colatitude,point.longitude,node.element});pointsConstruction.stop();mathprofile::Scope tensorCall("driver_tensor_sampling_call");auto fields=workspace.sampleForKernelsTensor(points,points.size(),&plans[section],2,{},false);tensorCall.stop();
     if(k<2){
#ifndef MATHPROFILE_DISABLED
      const auto beforeDiagnosticSnapshot=Clock::now();auto diagnosticBefore=mathprofile::snapshot();profileSnapshotOverheadSeconds+=seconds(beforeDiagnosticSnapshot);
#endif
      auto validationStart=Clock::now();mathprofile::Scope tensorValidation("driver_tensor_vs_packed_diagnostic_validation");std::vector<SPARSESPEC::Field3DPoint> probe;std::vector<std::size_t> ids;for(std::size_t a=0;a<sampling.places.size();a+=std::max<std::size_t>(1,sampling.places.size()/7))for(std::size_t n=0;n<nodes.size();++n){ids.push_back(a*nodes.size()+n);probe.push_back(points[ids.back()]);}auto reference=workspace.sampleForKernelsPacked(probe,probe.size(),nullptr,2);auto array=[](auto v){return std::array<Complex,6>{v.rr,v.horizontalTrace,v.radialPlus,v.radialMinus,v.horizontalPlus,v.horizontalMinus};};double error=0,norm=0;for(std::size_t i=0;i<ids.size();++i)for(int a=0;a<2;++a){auto x=array(a?fields[ids[i]].adjointDualStrain:fields[ids[i]].forwardRawStrain);auto y=array(a?reference[i].adjointDualStrain:reference[i].forwardRawStrain);for(int c=0;c<6;++c){error+=std::norm(x[c]-y[c]);norm+=std::norm(y[c]);}}tensorReferenceError=std::max(tensorReferenceError,std::sqrt(error/std::max(norm,1e-300)));if(tensorReferenceError>1e-11)throw std::runtime_error("tensor versus packed native-side regression failed");tensorValidation.stop();tensorValidationSeconds+=seconds(validationStart);
#ifndef MATHPROFILE_DISABLED
      const auto afterDiagnosticSnapshot=Clock::now();auto diagnosticAfter=mathprofile::snapshot();auto diagnosticDelta=mathprofile::difference(diagnosticAfter,diagnosticBefore);mathprofile::accumulate(diagnosticTotals,diagnosticDelta);profileSnapshotOverheadSeconds+=seconds(afterDiagnosticSnapshot);
#endif
     }
     mathprofile::Scope fusedMaterialContraction("kernel_fused_material_evaluation_weak_form_contraction_frequency_accumulation");
#ifndef MATHPROFILE_DISABLED
     const auto fluidNodes=static_cast<std::uint64_t>(std::count_if(nodes.begin(),nodes.end(),[](const Node &node){return !node.solid;}));mathprofile::count("kernel_point_samples",static_cast<std::uint64_t>(points.size()));mathprofile::count("kernel_fluid_point_samples",fluidNodes*static_cast<std::uint64_t>(sampling.places.size()));
#endif
     for(std::size_t a=0;a<sampling.places.size();++a)for(std::size_t n=0;n<nodes.size();++n)values[section][n*sampling.places.size()+a]+=kernelVs(sem,fields[a*nodes.size()+n],model(nodes[n].layer,nodes[n].radius/scales.Length()),d);fusedMaterialContraction.stop();
    }
    std::cout<<"DSpec frequencies "<<k+1<<'/'<<grid.Count()<<std::endl;
   }
  }
 }
 kernelSeconds=frequencyRoute?seconds(timer):seconds(timer)-tensorValidationSeconds;
 kernelPipeline.stop();
#ifndef MATHPROFILE_DISABLED
 capturePhase(frequencyRoute?"kernel_frequency_parallel_pass_and_ordered_reduction_including_diagnostics":"kernel_production_including_separately_isolated_diagnostics");
#endif
 mathprofile::Scope kernelOutput("driver_kernel_output_and_result_io");
 for(std::size_t s=0;s<samplings.size();++s){std::string name=s==0?"along":"across";std::ofstream f(out/(name+".bin"),std::ios::binary);f.write(reinterpret_cast<char*>(values[s].data()),values[s].size()*sizeof(double));std::ofstream a(out/(name+"-angles.csv"));a<<"colatitude_rad,longitude_rad,along_deg,across_deg,angular_weight_sr\n"<<std::setprecision(17);for(std::size_t p=0;p<samplings[s].places.size();++p){auto x=samplings[s].places[p];a<<x.colatitude<<','<<x.longitude<<','<<x.along<<','<<x.across<<','<<(gridKind=="audit"?angularWeights[p]:0)<<'\n';}}
 double integrated=0;if(gridKind=="audit")for(std::size_t n=0;n<nodes.size();++n)for(std::size_t a=0;a<angularWeights.size();++a)integrated+=values[0][n*angularWeights.size()+a]*profile(nodes[n].radius)*nodes[n].radius*nodes[n].radius*nodes[n].dr*angularWeights[a];
 {std::ofstream f(out/"result.json");f<<std::setprecision(17)<<"{\"backend\":\""<<backend<<"\",\"degree\":"<<degree<<",\"nq\":"<<nq<<",\"refinement\":"<<refinement<<",\"relative_error\":"<<parameters.relative_error()<<",\"mesh_selection\":\""<<(meshSelectionEnv?"standard FullSpecMultiSem highest-complex-frequency h, optional h/2":"caller-specified maxstep")<<"\",\"zmax_real\":"<<(meshSelectionEnv?meshSelection.highestComplexOmega.real():grid.Omega(grid.Count()-1).real())<<",\"zmax_imag\":"<<(meshSelectionEnv?meshSelection.highestComplexOmega.imag():grid.Omega(grid.Count()-1).imag())<<",\"omega_reference\":"<<(meshSelectionEnv?meshSelection.omegaReference:std::abs(grid.Omega(grid.Count()-1)))<<",\"nominal_h\":"<<(meshSelectionEnv?meshSelection.nominalStep:refinement)<<",\"finer_mesh\":"<<(meshSelectionEnv&&meshSelection.finer?"true":"false")<<",\"kernel_adjoint_weight_source\":\""<<(frozenKernelWeights?"frozen accepted Phase A waveform adjoint source":"this mesh's own waveform-derived adjoint source")<<"\",\"element_count\":"<<sem.mesh().NE()<<",\"radial_dofs_R\":"<<sem.hR().rows()<<",\"radial_dofs_S\":"<<sem.hS(degree).rows()<<",\"radial_dofs_T\":"<<sem.hTk(degree).rows()<<",\"source_element\":"<<sem.sourceElement(source)<<",\"receiver_element\":"<<receiverElement<<",\"workers\":"<<workers<<",\"schedule_mode\":\""<<(frequencyRoute?"frequency":"degree")<<"\",\"frequency_workers\":"<<(frequencyRoute?frequencyWorkers:1)<<",\"solve_workers\":"<<solveWorkers<<",\"kernel_timing_semantics\":\""<<(frequencyRoute?"kernel_s includes setup, parallel pass with diagnostics, and parent ordered reduction":"legacy degree kernel_s excludes tensor validation")<<"\",\"kernel_parallel_pass_wall_s\":"<<kernelParallelPassSeconds<<",\"reduction_wall_s\":"<<reductionSeconds<<",\"diagnostic_worker_elapsed_s\":"<<diagnosticWorkerElapsedSeconds<<",\"frequency_contribution_buffer_bytes\":"<<contributionBufferBytes<<",\"attenuation\":"<<(in.calculation.attenuation?"true":"false")<<",\"epsilon_perturbation\":"<<epsilon<<",\"frequency_count\":"<<grid.Count()<<",\"first_bin\":"<<grid.First()<<",\"last_bin\":"<<grid.Last()<<",\"fft_samples\":"<<grid.Samples()<<",\"fft_period_s\":"<<grid.Samples()*in.calculation.sampling.timeStep<<",\"damping_si\":"<<grid.Damping()/scales.Time()<<",\"native_side_count\":"<<nodes.size()<<",\"along_angles\":"<<samplings[0].places.size()<<",\"across_angles\":"<<(samplings.size()>1?samplings[1].places.size():0)<<",\"gravity_max_abs_nd_error\":"<<gravityError<<",\"time_shift_adjoint_identity\":"<<measured.delayAudit<<",\"measurement_rate_norm\":"<<measured.norm<<",\"integrated_profile_kernel\":"<<integrated<<",\"forward_s\":"<<forwardSeconds<<",\"measurement_s\":"<<measurementSeconds<<",\"kernel_s\":"<<kernelSeconds<<",\"tensor_packed_relative_l2\":"<<tensorReferenceError<<",\"tensor_validation_s\":"<<tensorValidationSeconds<<",\"cache_bytes\":"<<cache.retainedBytes()<<",\"elapsed_s\":"<<seconds(start)<<",\"kernel_parameter\":\"ln_vs_at_fixed_vp_density_and_effective_Q\",\"kernel_units\":\"s m^-3\",\"layout\":\"native-side-major then angle, float64\",\"source_shell_contacts\":false}\n";}
 kernelOutput.stop();
#ifndef MATHPROFILE_DISABLED
 capturePhase("kernel_output_and_result_io");
 pipelineScope.stop();
 const auto finalSnapshotStart=Clock::now();auto finalProfileSnapshot=mathprofile::snapshot();profileSnapshotOverheadSeconds+=seconds(finalSnapshotStart);
 writeProfileAccounting(out/"profile.json",std::move(phaseSnapshots),diagnosticTotals,finalProfileSnapshot,profileSnapshotOverheadSeconds,frequencyRoute);
#endif
 return 0;
}catch(const std::exception &e){std::cerr<<"FAIL "<<e.what()<<'\n';return 1;}
