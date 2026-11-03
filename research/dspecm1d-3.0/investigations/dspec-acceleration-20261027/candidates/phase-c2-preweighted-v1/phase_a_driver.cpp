// Matched documentation-case adapter: common model/measurement, independent solvers.
#include <YSpec/YSpec.hpp>
#include "Spectra.hpp"
#include "KernelParts.hpp"
#include "Sphere.hpp"
#include "matched_mesh_model.hpp"
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

#include "phase_a_reconstruction.hpp"
#include "phase_a_standard_forward.hpp"

namespace {
using RadialKey=std::tuple<int,int,int,int>; // bin, degree, element, local node

std::uint8_t familyBit(SPARSESPEC::SensitivitySolveFamily family) {
  return static_cast<std::uint8_t>(1u<<static_cast<unsigned>(family));
}

void setFrame(std::array<double,9> &target,
              const std::array<std::array<double,3>,3> &axes) {
  for(std::size_t r=0;r<3;++r)for(std::size_t c=0;c<3;++c)
    target[3*r+c]=axes[r][c];
}

std::vector<phasea::RadialSideDegreeRecord> collectRadial(
    SPARSESPEC::Windowed3DFrequencyWorkspace &workspace,Full1D::SEM &sem,int bin) {
  std::map<RadialKey,phasea::RadialSideDegreeRecord> bySideDegree;
  workspace.visitProjectedRadialModes(
    [&](SPARSESPEC::SensitivitySolveFamily family,int degree,bool adjoint,int receiver,
        int m,int element,int localNode,double radius,bool fluid,bool liftApplied,
        const std::array<std::array<double,3>,3> &axes,const Eigen::MatrixXcd &values) {
      if(values.rows()!=6||values.cols()!=1)throw std::runtime_error("visitor radial shape changed");
      const RadialKey key{bin,degree,element,localNode};
      auto [it,inserted]=bySideDegree.try_emplace(key);
      auto &record=it->second;
      if(inserted){
        record.bin=bin;record.receiver=0;record.degree=degree;
        record.layer=sem.mesh().LayerNumber(element);record.element=element;
        record.localNode=localNode;record.radiusNormalized=radius;
        record.radiusMetres=radius*sem.lengthScaleMetres();record.fluid=fluid?1:0;
        record.adjointFrameIsIsotropic=0;
        record.elementLowerNormalized=sem.mesh().ELR(element);
        record.elementUpperNormalized=sem.mesh().EUR(element);
        record.elementJacobianNormalized=.5*(record.elementUpperNormalized-record.elementLowerNormalized);
      } else if(record.fluid!=(fluid?1:0)||record.radiusNormalized!=radius||
                record.elementLowerNormalized!=sem.mesh().ELR(element)||
                record.elementUpperNormalized!=sem.mesh().EUR(element)||
                record.elementJacobianNormalized!=.5*(sem.mesh().EUR(element)-sem.mesh().ELR(element))||
                record.layer!=sem.mesh().LayerNumber(element)) {
        throw std::runtime_error("inconsistent metadata while merging radial family modes");
      }
      if(adjoint&&receiver<0){
        record.adjointFrameIsIsotropic=1;
        record.adjointFrameAxes={1,0,0,0,1,0,0,0,1};
        record.adjointFrameSet=1;
      }else{
        auto &stored=adjoint?record.adjointFrameAxes:record.sourceFrameAxes;
        auto &set=adjoint?record.adjointFrameSet:record.sourceFrameSet;
        if(set){
          for(std::size_t r=0;r<3;++r)for(std::size_t c=0;c<3;++c)
            if(stored[3*r+c]!=axes[r][c])
              throw std::runtime_error("inconsistent pole-frame axes while merging radial modes");
        }else {setFrame(stored,axes);set=1;}
      }
      const std::size_t slot=adjoint?phasea::SourceOrderSlots+static_cast<std::size_t>(m+1):static_cast<std::size_t>(m+2);
      if(slot>=phasea::ModalOrderSlots)throw std::runtime_error("radial mode order outside padded packet");
      const auto bit=familyBit(family);
      if(record.familyMask[slot]&bit)throw std::runtime_error("duplicate radial family mode slot");
      if(!adjoint&&record.familyMask[slot]!=0&&
         record.sourceLiftApplied[slot]!=(liftApplied?1:0))
        throw std::runtime_error("merged source-family modes disagree on source-lift application");
      if(!adjoint)record.sourceLiftApplied[slot]=liftApplied?1:0;
      record.familyMask[slot]|=bit;
      for(std::size_t c=0;c<6;++c){
        auto value=record.radial[slot*6+c].value()+values(static_cast<Eigen::Index>(c),0);
        record.radial[slot*6+c]=phasea::ComplexPair::from(value);
      }
    },3.0e6);
  std::vector<phasea::RadialSideDegreeRecord> result;result.reserve(bySideDegree.size());
  for(auto &[key,record]:bySideDegree){(void)key;result.push_back(record);}
  return result;
}

std::vector<phasea::ReferenceFieldRecord> collectReference(
    SPARSESPEC::Windowed3DFrequencyWorkspace &workspace,
    const std::vector<YSpec::KernelSampling> &samplings,
    const std::vector<Node> &nodes,int bin) {
  std::vector<phasea::ReferenceFieldRecord> records;
  for(std::size_t section=0;section<samplings.size();++section){
    const auto &sampling=samplings[section];
    std::vector<SPARSESPEC::Field3DPoint> points;
    points.reserve(sampling.places.size()*nodes.size());
    for(const auto &point:sampling.places)for(const auto &node:nodes)
      points.push_back({node.radius,point.colatitude,point.longitude,node.element});
    auto fields=workspace.sampleForKernelsPacked(points,points.size(),nullptr,2);
    for(std::size_t a=0;a<sampling.places.size();++a)for(std::size_t n=0;n<nodes.size();++n){
      const std::size_t p=a*nodes.size()+n;
      for(int role=0;role<2;++role){
        const auto &strain=role?fields[p].adjointDualStrain:fields[p].forwardRawStrain;
        phasea::ReferenceFieldRecord record;record.bin=bin;
        record.role=role?phasea::FieldRole::AdjointDual:phasea::FieldRole::ForwardSource;
        record.receiver=role?0:-1;record.section=static_cast<int>(section);
        record.layer=nodes[n].layer;
        record.element=nodes[n].element;record.localNode=nodes[n].node;
        record.sectionPoint=static_cast<int>(a);record.fluid=fields[p].fluid?1:0;
        record.radiusNormalized=nodes[n].radius/workspace.sem().lengthScaleMetres();
        record.thetaRadians=sampling.places[a].colatitude;
        record.phiRadians=sampling.places[a].longitude;
        const std::array<Complex,6> values{strain.rr,strain.horizontalTrace,strain.radialPlus,
            strain.radialMinus,strain.horizontalPlus,strain.horizontalMinus};
        for(std::size_t c=0;c<6;++c)record.strain[c]=phasea::ComplexPair::from(values[c]);
        records.push_back(record);
      }
    }
  }
  return records;
}

int runA1(const fs::path &input,const fs::path &out,int degree,int nq,
          double refinement,int solveWorkers,const fs::path &tracePath,
          const fs::path &resultPath) {
  if(fs::exists(out)&&fs::directory_iterator(out)!=fs::directory_iterator{})
    throw std::runtime_error("A1 output directory must be empty");
  fs::create_directories(out);
  Eigen::setNbThreads(1);omp_set_dynamic(0);omp_set_num_threads(1);
  mathprofile::Root root("phase-a-radial-packet",-1);
  mathprofile::Scope pipeline("phase_a_a1_pipeline");
  auto in=YSpec::ReadKernelInput(input);in.calculation.degrees.max=degree;
  if(in.calculation.quantity!=YSpec::Quantity::Velocity||in.measurements.size()!=1||
     in.measurements[0].component!="T")throw std::runtime_error("expected fixed transverse-velocity case");
  auto model=YSpec::Model::ReadDeck(in.model);CommonModel common{model,0.0};
  const auto &scales=model.GetScales();auto grid=YSpec::MakeGrid(in.calculation,scales);
  auto geoLat=[](double lat){return 90-YSpec::Sphere::Colatitude(lat)*180/pi;};
  const auto &moment=in.source.moment;
  SourceInfo::EarthquakeCMT source(geoLat(in.source.location.latitude),in.source.location.longitude,
      in.source.depth/1000,moment.rr,moment.rt,moment.rp,moment.tt,moment.tp,moment.pp);
  Full1D::SEM sem(common,refinement,nq,degree,source);
  {std::ofstream pf(out/"parameters.txt");pf<<params(in,degree);}
  InputParameters parameters((out/"parameters.txt").string());
  std::vector<Node> nodes;auto gll=GaussQuad::GaussLobattoLegendreQuadrature1D<double>(nq);
  double gravityError=0.0;
  for(int e=0;e<sem.mesh().NE();++e)for(int q=0;q<nq;++q){
    const double r=sem.mesh().NodeRadius(e,q);const int layer=sem.mesh().LayerNumber(e);
    gravityError=std::max(gravityError,std::abs(sem.meshModel().Gravity(e,q)-model(layer,r).g));
    if(r*scales.Length()>=3.0e6)nodes.push_back({e,q,layer,r*scales.Length(),
       .5*sem.mesh().EW(e)*gll.W(q)*scales.Length(),!sem.elementIsFluid(e)});
  }
  if(gravityError>1e-10||nodes.size()!=293)throw std::runtime_error("fixed mesh/gravity/side inventory mismatch");
  const double spacing=360.0/(4*degree+1);
  YSpec::SectionOptions options{0,4,spacing,500e3};
  std::vector<YSpec::KernelSampling> samplings(2);
  samplings[0]=YSpec::PathSection(model,in.source.location,in.receiver.location,options);
  samplings[1]=YSpec::CrossSection(model,in.source.location,in.receiver.location,options);
  for(auto &s:samplings){s.nodes.clear();for(const auto &n:nodes)s.nodes.push_back({n.layer,n.radius});}
  const auto dir=in.measurements[0].direction;Eigen::Vector3cd direction;direction<<dir.up,dir.north,dir.east;
  const int receiverElement=sem.receiverElements(parameters,true).front();
  if(parameters.num_receivers()!=1)
    throw std::runtime_error("Phase A packet schema expects the frozen single receiver");
  const double receiverRadius=scales.Length()*model.SurfaceRadius()-in.receiver.depth;
  std::vector<SPARSESPEC::Field3DPoint> receiver{{receiverRadius,
      YSpec::Sphere::Colatitude(in.receiver.location.latitude),
      in.receiver.location.longitude*pi/180,receiverElement}};
  auto dispersion=[&](Complex z){return in.calculation.attenuation?
      2/pi*std::log(Complex(0,1)*z/model.ReferenceFrequency()):Complex{};};
  YSpec::Measurement measurement{in.receiver,in.measurements[0].direction,{}};
  std::ifstream trace(tracePath);std::string line;
  if(!trace||!std::getline(trace,line)||line!="time_s,velocity_m_s,rate_m_s2,window,adjoint_source")
    throw std::runtime_error("A1 cannot read the frozen accepted adjoint-source trace");
  while(std::getline(trace,line)){
    std::stringstream row(line);std::string field;std::array<std::string,5> columns;
    for(std::size_t i=0;i<columns.size();++i)
      if(!std::getline(row,columns[i],','))throw std::runtime_error("A1 frozen trace row has an invalid column count");
    if(std::getline(row,field,','))throw std::runtime_error("A1 frozen trace row has trailing columns");
    measurement.adjointSource.push_back(std::stod(columns[4]));
  }
  if(measurement.adjointSource.size()!=YSpec::SampleCount(grid,in.calculation.sampling))
    throw std::runtime_error("A1 frozen adjoint-source sample count differs from the case");
  auto adjoint=YSpec::MakeAdjointSpectrum(model,in.calculation,measurement);
  const auto readJsonNumber=[&](const std::string &name){
    std::ifstream file(resultPath);if(!file)throw std::runtime_error("A1 missing frozen result metadata");
    const std::string content((std::istreambuf_iterator<char>(file)),{});
    const auto key=content.find("\""+name+"\":");if(key==std::string::npos)
      throw std::runtime_error("A1 missing frozen result number: "+name);
    const char *begin=content.c_str()+key+name.size()+3;char *end=nullptr;
    const double value=std::strtod(begin,&end);
    if(end==begin||!std::isfinite(value))throw std::runtime_error("A1 invalid frozen result number: "+name);
    return value;
  };
  const double delayIdentity=readJsonNumber("time_shift_adjoint_identity");
  const double measurementNorm=readJsonNumber("measurement_rate_norm");
  if(std::abs(delayIdentity-1.0)>1e-10||!(measurementNorm>0.0)||
     readJsonNumber("degree")!=degree||readJsonNumber("nq")!=nq||
     readJsonNumber("refinement")!=refinement||readJsonNumber("frequency_count")!=35||
     readJsonNumber("first_bin")!=7||readJsonNumber("last_bin")!=41||
     readJsonNumber("fft_samples")!=2048||readJsonNumber("fft_period_s")!=4096||
     readJsonNumber("native_side_count")!=293||readJsonNumber("along_angles")!=135||
     readJsonNumber("across_angles")!=19)
    throw std::runtime_error("A1 frozen Stage 2 trace/result case or time identity gate failed");
  {std::ifstream file(resultPath);const std::string content((std::istreambuf_iterator<char>(file)),{});
    if(content.find("\"attenuation\":true")==std::string::npos)
      throw std::runtime_error("A1 frozen trace result does not match attenuated case");}
  std::ofstream active(out/"active_bins.csv");
  active<<"bin,z_re,z_im,weight_norm\n"<<std::setprecision(17);
  for(const int bin:{9,25,40}){
    const std::size_t k=static_cast<std::size_t>(bin-grid.First());
    if(k>=grid.Count())throw std::runtime_error("active bin outside exact grid");
    const Complex z=grid.Omega(k);
    const Eigen::Vector3cd weights=direction*std::conj(adjoint.weight[k]*adjoint.force[k]);
    if(!(weights.norm()>0.0))throw std::runtime_error("active adjoint weight unexpectedly zero");
    dssolver::BatchResponseCache cache({static_cast<std::size_t>(2*degree)*grid.Count()+1,16ull<<30,1e-12});
    auto workspace=SPARSESPEC::Windowed3DFieldEngine::prepareFrequencyWorkspace(
        sem,source,parameters,z,weights,&cache,true,false,nullptr,solveWorkers,dispersion(z));
    auto radial=collectRadial(workspace,sem,bin);
    auto reference=collectReference(workspace,samplings,nodes,bin);
    phasea::PacketHeader header;header.bin=bin;header.zReal=z.real();header.zImag=z.imag();
    header.lengthScaleMetres=scales.Length();header.timeScaleSeconds=scales.Time();
    phasea::writePacket(out/("radial-bin"+std::to_string(bin)+".bin"),header,radial);
    phasea::writePacket(out/("reference-bin"+std::to_string(bin)+".bin"),header,reference);
    active<<bin<<','<<z.real()<<','<<z.imag()<<','<<weights.norm()<<'\n';
  }
  std::ofstream meta(out/"packet-metadata.json");
  meta<<std::setprecision(17)<<"{\"case\":\"L200_NQ10_refinement0.025\",\"degree\":"<<degree
      <<",\"native_sides\":"<<nodes.size()<<",\"along_points\":"<<samplings[0].places.size()
      <<",\"across_points\":"<<samplings[1].places.size()<<",\"gravity_max_abs_nd_error\":"<<gravityError
      <<",\"delay_identity\":"<<delayIdentity<<",\"rate_norm\":"<<measurementNorm
      <<",\"receiver_count\":"<<parameters.num_receivers()
      <<",\"weight_source\":\"accepted Stage 2 trace adjoint_source plus unchanged MakeAdjointSpectrum\",\"frequency_solves\":3"
      <<",\"derivative_convention\":\"existing Dref/J and canonical native node\",\"adjoint_field\":\"Hermitian weighted dual; existing projection\"}\n";
  pipeline.stop();
  return 0;
}

using SideKey=std::tuple<int,int>;
using TableKey=std::tuple<int,int,int,int>; // element, node, role, receiver
using AngularKey=std::tuple<int,int,int,int>; // section, point, role, receiver
struct AngularEval {
  std::unique_ptr<phasea::LegendreTable> legendre;
  std::vector<Eigen::MatrixXcd> eigenBasis;
  int orders=0;
  double theta=0.0,phi=0.0,a=0.0,b=0.0;
};
using Axes=std::array<double,9>;

std::array<double,4> localAngularCoordinates(const Axes &axes,double theta,double phi) {
  const std::array<double,3> global{std::sin(theta)*std::cos(phi),
      std::sin(theta)*std::sin(phi),std::cos(theta)};
  std::array<double,3> local{};
  for(std::size_t j=0;j<3;++j)for(std::size_t i=0;i<3;++i)
    local[j]+=axes[3*j+i]*global[i];
  const double horizontal=std::hypot(local[0],local[1]);
  const double localTheta=std::atan2(horizontal,local[2]);
  const double localPhi=horizontal==0.0?0.0:std::atan2(local[1],local[0]);
  const std::array<double,3> localEt{std::cos(localTheta)*std::cos(localPhi),
      std::cos(localTheta)*std::sin(localPhi),-std::sin(localTheta)};
  std::array<double,3> globalEt{},globalEp{-std::sin(phi),std::cos(phi),0.0};
  for(std::size_t i=0;i<3;++i)for(std::size_t j=0;j<3;++j)
    globalEt[i]+=axes[3*j+i]*localEt[j];
  const std::array<double,3> referenceEt{std::cos(theta)*std::cos(phi),
      std::cos(theta)*std::sin(phi),-std::sin(theta)};
  double a=0.0,b=0.0;
  for(std::size_t i=0;i<3;++i){a+=referenceEt[i]*globalEt[i];b+=globalEp[i]*globalEt[i];}
  return {localTheta,localPhi,a,b};
}

phasea::Strain evaluateEigen(const phasea::StrainTable &table,std::size_t node,
                             const AngularEval &angular) {
  phasea::Strain out{};
  for(int m=-table.orders();m<=table.orders();++m)
    for(std::size_t component=0;component<phasea::Components;++component){
      const std::size_t basisIndex=component*static_cast<std::size_t>(2*angular.orders+1)+
          static_cast<std::size_t>(m+angular.orders);
      const auto &basis=angular.eigenBasis.at(basisIndex);
      Eigen::Map<const Eigen::VectorXcd> coefficients(table.row(node,component,m),
          static_cast<Eigen::Index>(table.degrees()));
      out[component]+=(basis*coefficients)(0,0);
    }
  return out;
}

double kernelVsAtLength(double length,const SPARSESPEC::KernelField3DSample &sample,
                        const YSpec::Material &m,Complex dispersion) {
  if(sample.fluid)return 0.0;
  const auto &f=sample.forwardRawStrain;const auto &a=sample.adjointDualStrain;
  const Complex F=-(std::conj(a.rr)*f.horizontalTrace+
                    std::conj(a.horizontalTrace)*f.rr);
  const Complex L=-2.0*(std::conj(a.radialPlus)*f.radialPlus+
                        std::conj(a.radialMinus)*f.radialMinus);
  const Complex N=-(.5*(std::conj(a.horizontalPlus)*f.horizontalPlus+
                        std::conj(a.horizontalMinus)*f.horizontalMinus)-
                    std::conj(a.horizontalTrace)*f.horizontalTrace);
  return 2*m.mu*std::real(-2.0*F*(1.0+m.qF*dispersion)+
      (L+N)*(1.0+m.qMu*dispersion))/std::pow(length,3);
}

int runA2(const fs::path &input,const fs::path &out,
          const std::array<fs::path,3> &radialPaths,
          const std::array<fs::path,3> &referencePaths) {
  if(fs::exists(out)&&fs::directory_iterator(out)!=fs::directory_iterator{})
    throw std::runtime_error("A2 output directory must be empty");
  fs::create_directories(out);
  Eigen::setNbThreads(1);omp_set_dynamic(0);omp_set_num_threads(1);
  mathprofile::Root root("phase-a-reconstruction-replay",-1);
  mathprofile::Scope pipeline("phase_a_a2_pipeline");
  const auto in=YSpec::ReadKernelInput(input);auto model=YSpec::Model::ReadDeck(in.model);
  const auto &scales=model.GetScales();
  std::ofstream summary(out/"a2-summary.json");
  std::ofstream csv(out/"a2-field-comparison.csv");
  std::ofstream kernelCsv(out/"a2-kernel-products.csv");
  kernelCsv<<"bin,section,element,node,section_point,fluid,packed_reference_kernel_vs,old_ds_kernel_vs,new_common_kernel_vs,ds_vs_common_abs_error,packed_vs_common_abs_error\n"<<std::setprecision(17);
  csv<<"bin,section,role,receiver,element,node,section_point,fluid,component,reference_projected_re,reference_projected_im,old_ds_replay_re,old_ds_replay_im,new_projected_re,new_projected_im,new_raw_re,new_raw_im,new_vs_reference_abs_error\n"<<std::setprecision(17);
  struct Metric {long double error2=0,norm2=0;double maxAbs=0;std::size_t count=0;
    int maxElement=-1,maxNode=-1,maxPoint=-1,maxComponent=-1;};
  using MetricKey=std::tuple<int,int,int,int>; // bin, section, role, fluid
  std::map<MetricKey,Metric> metrics;
  using ApprovedMetricKey=std::tuple<int,int,int>; // bin, section, role; approved gate spans solid+fluid
  std::map<ApprovedMetricKey,Metric> approvedMetrics,approvedDSMetrics;
  using PairKey=std::tuple<int,int,int,int,int>; // bin, section, element, node, point
  struct FieldPair {
    std::array<std::array<Complex,6>,2> reference{},dsReplay{},fresh{};
    std::array<bool,2> seen{};
    int layer=-1;double radius=0.0;bool fluid=false;Complex z{};
  };
  std::map<PairKey,FieldPair> pairedFields;
  std::map<MetricKey,Metric> dsMetrics;
  struct KernelMetric {long double error2=0,norm2=0,referenceError2=0,referenceNorm2=0;
    double maxAbs=0,peak=0,referenceMaxAbs=0,referencePeak=0;std::size_t count=0;};
  using KernelMetricKey=std::tuple<int,int>; // bin, section
  std::map<KernelMetricKey,KernelMetric> kernelMetrics;
  double packetReadSeconds=0.0,radialPreparationSeconds=0.0,angularPreparationSeconds=0.0;
  double oldReplaySeconds=0.0,commonReplaySeconds=0.0,referenceComparisonSeconds=0.0;
  double kernelProductSeconds=0.0;
  double maximumRelative=0.0,maximumDSCommonRelative=0.0;
  std::uint64_t totalReferences=0;
  std::map<AngularKey,AngularEval> angularCache;
  for(std::size_t task=0;task<3;++task){
    phasea::PacketHeader radialHeader,referenceHeader;
    const auto packetReadStart=Clock::now();
    auto radial=phasea::readPacket<phasea::RadialSideDegreeRecord>(radialPaths[task],radialHeader);
    auto references=phasea::readPacket<phasea::ReferenceFieldRecord>(referencePaths[task],referenceHeader);
    packetReadSeconds+=seconds(packetReadStart);
    if(radial.size()!=201u*293u)
      throw std::runtime_error("A2 radial packet does not contain all 201 degrees and 293 native sides");
    std::set<std::tuple<int,int,int>> radialKeys;
    for(const auto &r:radial)if(!radialKeys.emplace(r.degree,r.element,r.localNode).second)
      throw std::runtime_error("A2 radial packet has duplicate side-degree records");
    if(references.size()!=2u*293u*(135u+19u))
      throw std::runtime_error("A2 packed-reference packet does not cover both roles and both sections");
    if(radialHeader.bin!=referenceHeader.bin||radialHeader.bin!=std::array<int,3>{9,25,40}[task]||
       radialHeader.zReal!=referenceHeader.zReal||radialHeader.zImag!=referenceHeader.zImag||
       radialHeader.lengthScaleMetres!=referenceHeader.lengthScaleMetres||
       radialHeader.timeScaleSeconds!=referenceHeader.timeScaleSeconds)
      throw std::runtime_error("A2 packet headers do not bind to the same active-bin response");
    const int bin=radialHeader.bin;
    std::map<TableKey,std::unique_ptr<phasea::StrainTable>> tables;
    std::map<TableKey,std::unique_ptr<phasea::StrainTable>> dsTables;
    std::map<TableKey,bool> isotropic;
    std::map<TableKey,Axes> frames;
    const auto radialPreparationStart=Clock::now();
    for(const auto &record:radial){
      for(int role=0;role<2;++role){
        const std::size_t roleBegin=role?phasea::SourceOrderSlots:0;
        const std::size_t roleEnd=role?phasea::ModalOrderSlots:phasea::SourceOrderSlots;
        bool hasRoleData=false;
        for(std::size_t slot=roleBegin;slot<roleEnd;++slot)
          hasRoleData=hasRoleData||record.familyMask[slot]!=0;
        if(!hasRoleData)continue;
        const int receiver=role?(record.degree==0&&record.adjointFrameIsIsotropic?-1:0):-1;
        const TableKey key{record.element,record.localNode,role,receiver};
        const int orders=role?1:2;
        auto [it,added]=tables.try_emplace(key);
        auto [dsIt,dsAdded]=dsTables.try_emplace(key);
        if(added)it->second=std::make_unique<phasea::StrainTable>(1,orders,200);
        if(dsAdded)dsIt->second=std::make_unique<phasea::StrainTable>(1,orders,200);
        if(role?(!record.adjointFrameSet):(!record.sourceFrameSet))
          throw std::runtime_error("A2 populated packet role has no frame metadata");
        isotropic[key]=role&&record.adjointFrameIsIsotropic!=0;
        frames[key]=role?record.adjointFrameAxes:record.sourceFrameAxes;
        for(std::size_t slot=roleBegin;slot<roleEnd;++slot)
          if(record.familyMask[slot]!=0){
            phasea::assignRadialMode(*it->second,0,record,slot,
                role?phasea::FieldRole::AdjointDual:phasea::FieldRole::ForwardSource);
            phasea::assignRadialModeDS(*dsIt->second,0,record,slot,
                role?phasea::FieldRole::AdjointDual:phasea::FieldRole::ForwardSource);
          }
      }
    }
    radialPreparationSeconds+=seconds(radialPreparationStart);
    const auto angularPreparationStart=Clock::now();
    for(const auto &reference:references){
      const int role=reference.role==phasea::FieldRole::AdjointDual?1:0;
      const TableKey sourceKey{reference.element,reference.localNode,role,-1};
      const TableKey receiverKey{reference.element,reference.localNode,role,0};
      for(const auto &key:{sourceKey,receiverKey}){
        if(!tables.contains(key))continue;
        const bool iso=isotropic.at(key);
        const AngularKey ak{reference.section,reference.sectionPoint,role,std::get<3>(key)};
        double localTheta=reference.thetaRadians,localPhi=reference.phiRadians,a=0.0,b=0.0;
        if(!iso){
          const auto transformed=localAngularCoordinates(frames.at(key),reference.thetaRadians,reference.phiRadians);
          localTheta=transformed[0];localPhi=transformed[1];a=transformed[2];b=transformed[3];
        }
        const auto cached=angularCache.find(ak);
        if(cached!=angularCache.end()){
          const auto &old=cached->second;
          if(old.theta!=localTheta||old.phi!=localPhi||old.a!=a||old.b!=b)
            throw std::runtime_error("A2 immutable angular key changed between active bins");
          continue;
        }
        AngularEval eval;eval.theta=localTheta;eval.phi=localPhi;eval.a=a;eval.b=b;
        eval.orders=role?1:2;
        eval.legendre=std::make_unique<phasea::LegendreTable>(200,eval.orders,
            std::vector<double>{localTheta});
        eval.eigenBasis.reserve(phasea::Components*static_cast<std::size_t>(2*eval.orders+1));
        for(std::size_t component=0;component<phasea::Components;++component)
          for(int m=-eval.orders;m<=eval.orders;++m){
            const Complex phase=std::exp(Complex(0.0,m*localPhi));
            Eigen::MatrixXcd basis(1,201);
            const auto *values=eval.legendre->row(0,phasea::Spin[component],m);
            for(int l=0;l<=200;++l)basis(0,l)=values[l]*phase;
            eval.eigenBasis.push_back(std::move(basis));
          }
        angularCache.emplace(ak,std::move(eval));
      }
    }
    angularPreparationSeconds+=seconds(angularPreparationStart);

    using ReferenceReplayKey=std::tuple<int,int,int,int,int>; // section,point,role,element,node
    std::map<ReferenceReplayKey,phasea::Strain> commonFields,oldDSFields;
    auto replay=[&](const phasea::ReferenceFieldRecord &reference,bool oldDS){
      const int role=reference.role==phasea::FieldRole::AdjointDual?1:0;
      const TableKey sourceKey{reference.element,reference.localNode,role,-1};
      const TableKey receiverKey{reference.element,reference.localNode,role,0};
      phasea::Strain result{};
      for(const auto &key:{sourceKey,receiverKey}){
        const auto found=tables.find(key);if(found==tables.end())continue;
        const AngularKey ak{reference.section,reference.sectionPoint,role,std::get<3>(key)};
        const auto &eval=angularCache.at(ak);
        phasea::Strain local{};
        if(oldDS){
          const auto dsFound=dsTables.find(key);
          if(dsFound==dsTables.end())throw std::runtime_error("A2 DS replay table missing");
          local=evaluateEigen(*dsFound->second,0,eval);
        } else {
          local=phasea::evaluate(*found->second,0,*eval.legendre,0,eval.phi);
        }
        if(!isotropic.at(key))local=phasea::rotateToGlobal(local,eval.a,eval.b);
        for(std::size_t c=0;c<phasea::Components;++c)result[c]+=local[c];
      }
      return result;
    };
    for(int section=0;section<2;++section)for(int role=0;role<2;++role){
      const auto commonStart=Clock::now();
      for(const auto &reference:references)
        if(reference.section==section&&
           static_cast<int>(reference.role==phasea::FieldRole::AdjointDual)==role){
          const ReferenceReplayKey key{section,reference.sectionPoint,role,
                                       reference.element,reference.localNode};
          commonFields.emplace(key,replay(reference,false));
        }
      commonReplaySeconds+=seconds(commonStart);
    }
    for(int section=0;section<2;++section)for(int role=0;role<2;++role){
      const auto oldStart=Clock::now();
      for(const auto &reference:references)
        if(reference.section==section&&
           static_cast<int>(reference.role==phasea::FieldRole::AdjointDual)==role){
          const ReferenceReplayKey key{section,reference.sectionPoint,role,
                                       reference.element,reference.localNode};
          oldDSFields.emplace(key,replay(reference,true));
        }
      oldReplaySeconds+=seconds(oldStart);
    }

    const auto comparisonStart=Clock::now();
    for(const auto &reference:references){
      const int role=reference.role==phasea::FieldRole::AdjointDual?1:0;
      const ReferenceReplayKey replayKey{reference.section,reference.sectionPoint,role,
                                          reference.element,reference.localNode};
      const auto rawReconstructed=commonFields.at(replayKey);
      const auto rawDSReplay=oldDSFields.at(replayKey);
      auto reconstructed=rawReconstructed;
      auto oldDSReplay=rawDSReplay;
      std::array<Complex,6> old{};
      for(std::size_t c=0;c<6;++c)old[c]=reference.strain[c].value();
      if(reference.fluid){
        reconstructed=phasea::inheritedFluidProjection(reconstructed);
        oldDSReplay=phasea::inheritedFluidProjection(oldDSReplay);
        old=phasea::inheritedFluidProjection(old);
      }
      const auto key=MetricKey{bin,reference.section,role,reference.fluid};
      auto &metric=metrics[key];++metric.count;++totalReferences;
      auto &dsMetric=dsMetrics[key];
      auto &approvedMetric=approvedMetrics[{bin,reference.section,role}];
      auto &approvedDSMetric=approvedDSMetrics[{bin,reference.section,role}];
      ++approvedMetric.count;++approvedDSMetric.count;
      auto &pair=pairedFields[{bin,reference.section,reference.element,
                               reference.localNode,reference.sectionPoint}];
      pair.reference[role]=old;pair.dsReplay[role]=oldDSReplay;
      pair.fresh[role]=reconstructed;pair.seen[role]=true;
      pair.layer=reference.layer;pair.radius=reference.radiusNormalized;
      pair.fluid=reference.fluid!=0;pair.z={radialHeader.zReal,radialHeader.zImag};
      for(std::size_t c=0;c<6;++c){
        const double error=std::abs(reconstructed[c]-old[c]);
        metric.error2+=std::norm(reconstructed[c]-old[c]);metric.norm2+=std::norm(old[c]);
        approvedMetric.error2+=std::norm(reconstructed[c]-old[c]);
        approvedMetric.norm2+=std::norm(old[c]);
        if(error>metric.maxAbs){metric.maxAbs=error;metric.maxElement=reference.element;
          metric.maxNode=reference.localNode;metric.maxPoint=reference.sectionPoint;
          metric.maxComponent=static_cast<int>(c);}
        const double dsError=std::abs(reconstructed[c]-oldDSReplay[c]);
        dsMetric.maxAbs=std::max(dsMetric.maxAbs,dsError);
        dsMetric.error2+=std::norm(reconstructed[c]-oldDSReplay[c]);
        dsMetric.norm2+=std::norm(oldDSReplay[c]);++dsMetric.count;
        approvedDSMetric.error2+=std::norm(reconstructed[c]-oldDSReplay[c]);
        approvedDSMetric.norm2+=std::norm(oldDSReplay[c]);
        csv<<bin<<','<<reference.section<<','<<role<<','<<reference.receiver<<','
           <<reference.element<<','<<reference.localNode<<','<<reference.sectionPoint<<','
           <<static_cast<int>(reference.fluid)<<','<<c<<','<<old[c].real()<<','<<old[c].imag()
           <<','<<oldDSReplay[c].real()<<','<<oldDSReplay[c].imag()<<','
           <<reconstructed[c].real()<<','<<reconstructed[c].imag()<<','
           <<rawReconstructed[c].real()<<','<<rawReconstructed[c].imag()<<','<<error<<'\n';
      }
    }
    referenceComparisonSeconds+=seconds(comparisonStart);
  }
  const auto kernelProductStart=Clock::now();
  for(const auto &[key,pair]:pairedFields){
    if(!pair.seen[0]||!pair.seen[1])throw std::runtime_error("A2 lacks paired forward/adjoint fields");
    SPARSESPEC::KernelField3DSample referenceField,oldDSField,newField;
    referenceField.fluid=oldDSField.fluid=newField.fluid=pair.fluid;
    referenceField.forwardRawStrain={pair.reference[0][0],pair.reference[0][1],pair.reference[0][2],pair.reference[0][3],pair.reference[0][4],pair.reference[0][5]};
    referenceField.adjointDualStrain={pair.reference[1][0],pair.reference[1][1],pair.reference[1][2],pair.reference[1][3],pair.reference[1][4],pair.reference[1][5]};
    oldDSField.forwardRawStrain={pair.dsReplay[0][0],pair.dsReplay[0][1],pair.dsReplay[0][2],pair.dsReplay[0][3],pair.dsReplay[0][4],pair.dsReplay[0][5]};
    oldDSField.adjointDualStrain={pair.dsReplay[1][0],pair.dsReplay[1][1],pair.dsReplay[1][2],pair.dsReplay[1][3],pair.dsReplay[1][4],pair.dsReplay[1][5]};
    newField.forwardRawStrain={pair.fresh[0][0],pair.fresh[0][1],pair.fresh[0][2],pair.fresh[0][3],pair.fresh[0][4],pair.fresh[0][5]};
    newField.adjointDualStrain={pair.fresh[1][0],pair.fresh[1][1],pair.fresh[1][2],pair.fresh[1][3],pair.fresh[1][4],pair.fresh[1][5]};
    const Complex dispersion=in.calculation.attenuation?
        2/pi*std::log(Complex(0,1)*pair.z/model.ReferenceFrequency()):Complex{};
    const auto material=model(pair.layer,pair.radius);
    const double referenceKernel=kernelVsAtLength(scales.Length(),referenceField,material,dispersion);
    const double oldKernel=kernelVsAtLength(scales.Length(),oldDSField,material,dispersion);
    const double newKernel=kernelVsAtLength(scales.Length(),newField,material,dispersion);
    const double absError=std::abs(oldKernel-newKernel);
    const double referenceError=std::abs(referenceKernel-newKernel);
    auto &km=kernelMetrics[{std::get<0>(key),std::get<1>(key)}];
    km.error2+=static_cast<long double>(absError)*absError;
    km.norm2+=static_cast<long double>(oldKernel)*oldKernel;
    km.referenceError2+=static_cast<long double>(referenceError)*referenceError;
    km.referenceNorm2+=static_cast<long double>(referenceKernel)*referenceKernel;
    km.maxAbs=std::max(km.maxAbs,absError);km.peak=std::max(km.peak,std::abs(oldKernel));
    km.referenceMaxAbs=std::max(km.referenceMaxAbs,referenceError);
    km.referencePeak=std::max(km.referencePeak,std::abs(referenceKernel));++km.count;
    const auto [bin,section,element,node,point]=key;
    kernelCsv<<bin<<','<<section<<','<<element<<','<<node<<','<<point<<','
        <<static_cast<int>(pair.fluid)<<','<<referenceKernel<<','<<oldKernel<<','<<newKernel<<','
        <<absError<<','<<referenceError<<'\n';
  }
  kernelProductSeconds=seconds(kernelProductStart);
  for(const auto &[key,metric]:approvedMetrics){
    const auto [bin,section,role]=key;
    const std::size_t expected=static_cast<std::size_t>(section?19:135)*293;
    if(metric.count!=expected||!(metric.norm2>0.0L))
      throw std::runtime_error("A2 approved packed-reference group has incomplete coverage or zero norm");
    const double relative=std::sqrt(static_cast<double>(metric.error2/metric.norm2));
    maximumRelative=std::max(maximumRelative,relative);
    if(relative>1e-11)throw std::runtime_error("A2 packed-reference aggregate complex field gate failed");
  }
  for(const auto &[key,metric]:approvedDSMetrics){
    const auto [bin,section,role]=key;
    const std::size_t expected=static_cast<std::size_t>(section?19:135)*293;
    if(metric.count!=expected||!(metric.norm2>0.0L))
      throw std::runtime_error("A2 approved DS/common group has incomplete coverage or zero norm");
    const double relative=std::sqrt(static_cast<double>(metric.error2/metric.norm2));
    maximumDSCommonRelative=std::max(maximumDSCommonRelative,relative);
    if(relative>1e-11)throw std::runtime_error("A2 DS/common aggregate complex field gate failed");
  }
  for(int bin:{9,25,40})for(int section=0;section<2;++section){
    const std::size_t expected=static_cast<std::size_t>(section?19:135)*293;
    std::size_t forwardCount=0,adjointCount=0;long double refAdjNorm=0,dsAdjNorm=0;
    for(const auto &[key,metric]:metrics){
      const auto [b,s,role,fluid]=key;(void)fluid;
      if(b==bin&&s==section){if(role==0)forwardCount+=metric.count;
        else {adjointCount+=metric.count;refAdjNorm+=metric.norm2;}}
    }
    for(const auto &[key,metric]:dsMetrics){
      const auto [b,s,role,fluid]=key;(void)fluid;
      if(b==bin&&s==section&&role==1)dsAdjNorm+=metric.norm2;
    }
    if(forwardCount!=expected||adjointCount!=expected||!(refAdjNorm>0.0L)||!(dsAdjNorm>0.0L))
      throw std::runtime_error("A2 section/role coverage or nonzero adjoint-field gate failed");
  }
  if(totalReferences==0||maximumRelative>1e-11)
    throw std::runtime_error("A2 reconstruction did not cover the bound reference packets");
  summary<<std::setprecision(17)<<"{\"bins\":[9,25,40],\"reference_records\":"<<totalReferences
      <<",\"max_reference_vs_common_field_l2\":"<<maximumRelative
      <<",\"max_same_packet_ds_vs_common_field_l2\":"<<maximumDSCommonRelative
      <<",\"field_gate\":1e-11,\"field_groups\":[";
  bool firstField=true;
  for(const auto &[key,metric]:metrics){
    if(!firstField)summary<<',';firstField=false;
    const auto [bin,section,role,fluid]=key;
    summary<<"{\"bin\":"<<bin<<",\"section\":"<<section<<",\"role\":"<<role
       <<",\"fluid\":"<<fluid<<",\"records\":"<<metric.count<<",\"relative_l2\":";
    if(metric.norm2>0)summary<<std::sqrt(static_cast<double>(metric.error2/metric.norm2));
    else summary<<"null";
    summary<<",\"max_abs\":"<<metric.maxAbs<<",\"max_location\":["
       <<metric.maxElement<<','<<metric.maxNode<<','<<metric.maxPoint<<','<<metric.maxComponent<<"]}";
  }
  summary<<"],\"approved_field_groups\":[";
  bool firstApproved=true;
  for(const auto &[key,metric]:approvedMetrics){
    if(!firstApproved)summary<<',';firstApproved=false;
    const auto [bin,section,role]=key;
    summary<<"{\"bin\":"<<bin<<",\"section\":"<<section<<",\"role\":"<<role
       <<",\"records\":"<<metric.count<<",\"relative_l2\":"
       <<std::sqrt(static_cast<double>(metric.error2/metric.norm2))<<"}";
  }
  summary<<"],\"ds_common_field_groups\":[";
  firstField=true;
  for(const auto &[key,metric]:dsMetrics){
    if(!firstField)summary<<',';firstField=false;
    const auto [bin,section,role,fluid]=key;
    summary<<"{\"bin\":"<<bin<<",\"section\":"<<section<<",\"role\":"<<role
       <<",\"fluid\":"<<fluid<<",\"records\":"<<metric.count<<",\"relative_l2\":";
    if(metric.norm2>0)summary<<std::sqrt(static_cast<double>(metric.error2/metric.norm2));
    else summary<<"null";
    summary<<",\"max_abs\":"<<metric.maxAbs<<"}";
  }
  summary<<"],\"approved_ds_common_field_groups\":[";
  firstApproved=true;
  for(const auto &[key,metric]:approvedDSMetrics){
    if(!firstApproved)summary<<',';firstApproved=false;
    const auto [bin,section,role]=key;
    summary<<"{\"bin\":"<<bin<<",\"section\":"<<section<<",\"role\":"<<role
       <<",\"records\":"<<metric.count<<",\"relative_l2\":"
       <<std::sqrt(static_cast<double>(metric.error2/metric.norm2))<<"}";
  }
  summary<<"],\"phase_timing_s\":{\"packet_read_io\":"<<packetReadSeconds
      <<",\"radial_strain_table_preparation\":"<<radialPreparationSeconds
      <<",\"immutable_angular_table_geometry_build\":"<<angularPreparationSeconds
      <<",\"old_ds_complex_eigen_replay_harness_elapsed_s\":"<<oldReplaySeconds
      <<",\"common_degree_inner_replay_harness_elapsed_s\":"<<commonReplaySeconds
      <<",\"packed_reference_comparison_and_projection\":"<<referenceComparisonSeconds
      <<",\"kernel_product_check\":"<<kernelProductSeconds<<"},\"kernel_products\":[";
  bool firstKernel=true;
  for(const auto &[key,km]:kernelMetrics){
    if(!firstKernel)summary<<',';firstKernel=false;
    const double rel=km.norm2>0?std::sqrt(static_cast<double>(km.error2/km.norm2)):0.0;
    const double refRel=km.referenceNorm2>0?std::sqrt(static_cast<double>(km.referenceError2/km.referenceNorm2)):0.0;
    const auto [bin,section]=key;
    summary<<"{\"bin\":"<<bin<<",\"section\":"<<section<<",\"samples\":"<<km.count
       <<",\"old_ds_common_relative_l2\":";
    if(km.norm2>0)summary<<rel;else summary<<"null";
    summary<<",\"old_ds_common_section_peak_scaled_max_abs\":";
    if(km.peak>0)summary<<km.maxAbs/km.peak;else summary<<"null";
    summary<<",\"packed_reference_common_relative_l2\":";
    if(km.referenceNorm2>0)summary<<refRel;else summary<<"null";
    summary<<",\"packed_reference_common_section_peak_scaled_max_abs\":";
    if(km.referencePeak>0)summary<<km.referenceMaxAbs/km.referencePeak;else summary<<"null";
    summary<<",\"old_ds_common_max_abs\":"<<km.maxAbs
       <<",\"packed_reference_common_max_abs\":"<<km.referenceMaxAbs<<"}";
  }
  summary<<"],\"kernel_product_relative_floor_gate\":null,\"fluid_projection\":\"old and new rr+h to rr, other channels zero\",\"source_lift\":\"packet values include accepted lift once on source element\",\"radial_solves\":0,\"old_reference_route\":\"same-workspace PACKED oracle, not Stage 2 TENSOR timing\",\"same_packet_old_route\":\"separate DS radial formula plus complex Eigen basis contraction\",\"replay_scope\":\"post-radial only; excludes solve/capture/full-pipeline costs\"}\n";
  pipeline.stop();

  return 0;
}

struct AuditTrace {
  std::vector<double> time,velocity,adjointSource;
};
AuditTrace readAuditTrace(const fs::path &path) {
  std::ifstream in(path);std::string line;
  if(!in||!std::getline(in,line)||line!="time_s,velocity_m_s,rate_m_s2,window,adjoint_source")
    throw std::runtime_error("A3 invalid frozen trace header");
  AuditTrace out;
  while(std::getline(in,line)){
    std::stringstream row(line);std::array<std::string,5> c;
    for(auto &v:c)if(!std::getline(row,v,','))throw std::runtime_error("A3 malformed frozen trace row");
    std::string extra;if(std::getline(row,extra,','))throw std::runtime_error("A3 frozen trace has trailing columns");
    out.time.push_back(std::stod(c[0]));out.velocity.push_back(std::stod(c[1]));
    out.adjointSource.push_back(std::stod(c[4]));
  }
  if(out.time.empty()||out.velocity.size()!=out.time.size()||
     out.adjointSource.size()!=out.time.size())throw std::runtime_error("A3 frozen trace sample count mismatch");
  return out;
}
std::map<std::string,double> readAuditResult(const fs::path &path) {
  std::ifstream in(path);if(!in)throw std::runtime_error("A3 missing frozen result.json");
  const std::string text((std::istreambuf_iterator<char>(in)),{});
  const std::array<std::string,19> names{"degree","nq","refinement","frequency_count","first_bin",
    "last_bin","fft_samples","fft_period_s","native_side_count","along_angles","across_angles",
    "gravity_max_abs_nd_error","time_shift_adjoint_identity","measurement_rate_norm",
    "integrated_profile_kernel","epsilon_perturbation","frequency_workers",
    "solve_workers","damping_si"};
  std::map<std::string,double> out;
  for(const auto &name:names){
    const auto key=text.find("\""+name+"\":");if(key==std::string::npos)
      throw std::runtime_error("A3 result metadata missing "+name);
    const char *begin=text.c_str()+key+name.size()+3;char *end=nullptr;
    const double value=std::strtod(begin,&end);
    if(end==begin||!std::isfinite(value))throw std::runtime_error("A3 invalid result metadata "+name);
    out[name]=value;
  }
  if(text.find("\"attenuation\":true")==std::string::npos)
    throw std::runtime_error("A3 frozen L3 fixture is not attenuated");
  return out;
}

int runA3(const fs::path &input,const fs::path &out,
          const fs::path &baseTracePath,const fs::path &baseResultPath,
          const fs::path &plusTracePath,const fs::path &plusResultPath,
          const fs::path &minusTracePath,const fs::path &minusResultPath,
          int degree,int nq,double refinement,int solveWorkers) {
  if(fs::exists(out)&&fs::directory_iterator(out)!=fs::directory_iterator{})
    throw std::runtime_error("A3 output directory must be empty");
  if(degree!=3||nq!=8||refinement!=0.05||solveWorkers!=1)
    throw std::runtime_error("A3 is restricted to the frozen L3/NQ8/refinement .05 audit");
  fs::create_directories(out);Eigen::setNbThreads(1);omp_set_dynamic(0);omp_set_num_threads(1);
  mathprofile::Root root("phase-a-l3-fd-kernel",-1);
  mathprofile::Scope pipeline("phase_a_a3_l3_audit_kernel");
  auto in=YSpec::ReadKernelInput(input);
  in.calculation.degrees.max=degree;
  if(in.calculation.quantity!=YSpec::Quantity::Velocity||in.measurements.size()!=1||
     in.measurements[0].component!="T")throw std::runtime_error("A3 case is not the frozen transverse-velocity audit");
  auto model=YSpec::Model::ReadDeck(in.model);CommonModel common{model,0.0};
  const auto &scales=model.GetScales();auto grid=YSpec::MakeGrid(in.calculation,scales);
  if(grid.Count()!=35||grid.First()!=7||grid.Last()!=41)
    throw std::runtime_error("A3 L3 case frequency grid changed");
  auto geoLat=[](double lat){return 90-YSpec::Sphere::Colatitude(lat)*180/pi;};
  const auto &moment=in.source.moment;
  SourceInfo::EarthquakeCMT source(geoLat(in.source.location.latitude),in.source.location.longitude,
      in.source.depth/1000,moment.rr,moment.rt,moment.rp,moment.tt,moment.tp,moment.pp);
  Full1D::SEM sem(common,refinement,nq,degree,source);
  {std::ofstream pf(out/"parameters.txt");pf<<params(in,degree);}
  InputParameters parameters((out/"parameters.txt").string());
  const auto baseMeta=readAuditResult(baseResultPath),plusMeta=readAuditResult(plusResultPath),
             minusMeta=readAuditResult(minusResultPath);
  for(const auto *meta:{&baseMeta,&plusMeta,&minusMeta}){
    if(meta->at("degree")!=degree||meta->at("nq")!=nq||meta->at("refinement")!=refinement||
       meta->at("frequency_count")!=35||meta->at("first_bin")!=7||meta->at("last_bin")!=41||
       meta->at("fft_samples")!=2048||meta->at("fft_period_s")!=4096||
       meta->at("native_side_count")!=157||meta->at("along_angles")!=91||
       meta->at("across_angles")!=0||std::abs(meta->at("time_shift_adjoint_identity")-1)>1e-10||
       meta->at("gravity_max_abs_nd_error")>1e-10||!(meta->at("measurement_rate_norm")>0))
      throw std::runtime_error("A3 frozen baseline/FD fixture metadata differs from this exact L3 case");
  }
  const double epsilon=plusMeta.at("epsilon_perturbation");
  if(!(epsilon>0)||minusMeta.at("epsilon_perturbation")!=-epsilon||
     baseMeta.at("epsilon_perturbation")!=0)
    throw std::runtime_error("A3 frozen +/- perturbations are not a matched central pair");
  const auto baseTrace=readAuditTrace(baseTracePath),plusTrace=readAuditTrace(plusTracePath),
             minusTrace=readAuditTrace(minusTracePath);
  const std::size_t expectedSamples=YSpec::SampleCount(grid,in.calculation.sampling);
  if(baseTrace.time.size()!=expectedSamples||baseTrace.time!=plusTrace.time||
     baseTrace.time!=minusTrace.time)
    throw std::runtime_error("A3 frozen FD traces do not share the exact case sample grid");
  const double dt=in.calculation.sampling.timeStep;
  for(std::size_t i=0;i<baseTrace.time.size();++i)
    if(baseTrace.time[i]!=static_cast<double>(i)*dt)
      throw std::runtime_error("A3 frozen FD trace sample time differs from the common case grid");
  long double finiteDifference=0.0L;
  for(std::size_t i=0;i<baseTrace.time.size();++i)
    finiteDifference+=dt*baseTrace.adjointSource[i]*
        (plusTrace.velocity[i]-minusTrace.velocity[i])/(2*epsilon);
  const auto direction=in.measurements[0].direction;
  Eigen::Vector3cd directionVector;directionVector<<direction.up,direction.north,direction.east;
  YSpec::Measurement measurement{in.receiver,in.measurements[0].direction,baseTrace.adjointSource};
  const auto adjoint=YSpec::MakeAdjointSpectrum(model,in.calculation,measurement);
  std::vector<Node> nodes;auto gll=GaussQuad::GaussLobattoLegendreQuadrature1D<double>(nq);
  double gravityError=0.0;
  for(int e=0;e<sem.mesh().NE();++e)for(int q=0;q<nq;++q){
    const double r=sem.mesh().NodeRadius(e,q);const int layer=sem.mesh().LayerNumber(e);
    gravityError=std::max(gravityError,std::abs(sem.meshModel().Gravity(e,q)-model(layer,r).g));
    if(r*scales.Length()>=3e6)nodes.push_back({e,q,layer,r*scales.Length(),
       .5*sem.mesh().EW(e)*gll.W(q)*scales.Length(),!sem.elementIsFluid(e)});
  }
  if(gravityError>1e-10||nodes.size()!=157)
    throw std::runtime_error("A3 native radial-side/gravity inventory differs from frozen L3 audit");
  auto quad=GaussQuad::Quadrature1D<double>(
      GaussQuad::LegendrePolynomial<double>{}.GaussQuadrature(2*degree+1));
  const int nPhi=4*degree+1;
  std::vector<std::array<double,3>> angles;
  for(int t=0;t<quad.N();++t)for(int p=0;p<nPhi;++p)
    angles.push_back({std::acos(quad.X(t)),2*pi*p/nPhi,quad.W(t)*2*pi/nPhi});
  const auto dispersion=[&](Complex z){return in.calculation.attenuation?
      2/pi*std::log(Complex(0,1)*z/model.ReferenceFrequency()):Complex{};};
  const auto bulkProfile=[](double r){return r>5.3e6&&r<5.8e6?
      std::pow(std::sin(pi*(r-5.3e6)/.5e6),2):0.0;};
  std::vector<double> commonKernelValues(nodes.size()*angles.size(),0.0);
  using AuditFrameKey=std::tuple<Axes,bool,int>;
  std::map<AuditFrameKey,std::vector<AngularEval>> auditAngularCache;
  std::map<TableKey,Axes> frozenFrames;
  std::map<TableKey,bool> frozenIsotropic;
  long double weightedAdjointNorm2=0.0L;
  double solveSeconds=0.0,reconstructionSeconds=0.0;
  for(std::size_t k=0;k<grid.Count();++k){
    const int bin=grid.First()+static_cast<int>(k);const Complex z=grid.Omega(k);
    const Eigen::Vector3cd weights=directionVector*std::conj(adjoint.weight[k]*adjoint.force[k]);
    dssolver::BatchResponseCache cache({static_cast<std::size_t>(2*degree)*grid.Count()+1,4ull<<30,1e-12});
    const auto solveStart=Clock::now();
    auto workspace=SPARSESPEC::Windowed3DFieldEngine::prepareFrequencyWorkspace(
        sem,source,parameters,z,weights,&cache,true,false,nullptr,solveWorkers,dispersion(z));
    solveSeconds+=seconds(solveStart);
    const auto radial=collectRadial(workspace,sem,bin);
    if(radial.size()!=static_cast<std::size_t>(degree+1)*nodes.size())
      throw std::runtime_error("A3 radial export omitted an L3 degree/native side");
    std::map<TableKey,std::unique_ptr<phasea::StrainTable>> tables;
    std::map<TableKey,Axes> frames;std::map<TableKey,bool> isotropic;
    const auto reconstructionStart=Clock::now();
    for(const auto &record:radial)for(int role=0;role<2;++role){
      const std::size_t begin=role?phasea::SourceOrderSlots:0;
      const std::size_t end=role?phasea::ModalOrderSlots:phasea::SourceOrderSlots;
      bool present=false;for(std::size_t slot=begin;slot<end;++slot)present|=record.familyMask[slot]!=0;
      if(!present)continue;
      const int receiver=role?(record.degree==0&&record.adjointFrameIsIsotropic?-1:0):-1;
      const TableKey key{record.element,record.localNode,role,receiver};
      auto [it,added]=tables.try_emplace(key);
      if(added)it->second=std::make_unique<phasea::StrainTable>(1,role?1:2,degree);
      frames[key]=role?record.adjointFrameAxes:record.sourceFrameAxes;
      isotropic[key]=role&&record.adjointFrameIsIsotropic!=0;
      const auto [frozenFrameIt,frameInserted]=frozenFrames.emplace(key,frames[key]);
      if(!frameInserted&&frozenFrameIt->second!=frames[key])
        throw std::runtime_error("A3 native frame axes changed across frequencies");
      const auto [frozenIsoIt,isoInserted]=frozenIsotropic.emplace(key,isotropic[key]);
      if(!isoInserted&&frozenIsoIt->second!=isotropic[key])
        throw std::runtime_error("A3 isotropic-frame tag changed across frequencies");
      for(std::size_t slot=begin;slot<end;++slot)if(record.familyMask[slot])
        phasea::assignRadialMode(*it->second,0,record,slot,
            role?phasea::FieldRole::AdjointDual:phasea::FieldRole::ForwardSource);
    }
    std::set<AuditFrameKey> activeAngularKeys;
    for(const auto &[key,table]:tables){
      const AuditFrameKey angularKey{frames.at(key),isotropic.at(key),table->orders()};
      activeAngularKeys.insert(angularKey);
      if(k==0){
        auto [cacheIt,inserted]=auditAngularCache.try_emplace(angularKey);
        if(!inserted)continue;
        auto &cache=cacheIt->second;cache.reserve(angles.size());
        const bool iso=isotropic.at(key);
        for(const auto &angle:angles){
          const double theta=angle[0],phi=angle[1];
          double localTheta=theta,localPhi=phi,a=0.0,b=0.0;
          if(!iso){const auto transformed=localAngularCoordinates(frames.at(key),theta,phi);
            localTheta=transformed[0];localPhi=transformed[1];a=transformed[2];b=transformed[3];}
          AngularEval eval;eval.theta=localTheta;eval.phi=localPhi;eval.a=a;eval.b=b;
          eval.orders=table->orders();
          eval.legendre=std::make_unique<phasea::LegendreTable>(
              degree,eval.orders,std::vector<double>{localTheta});
          cache.push_back(std::move(eval));
        }
      }
    }
    if(k==0){
      if(activeAngularKeys.size()!=auditAngularCache.size())
        throw std::runtime_error("A3 shared immutable angular cache does not cover all unique frames");
    }else if(activeAngularKeys.size()!=auditAngularCache.size()||
             !std::all_of(activeAngularKeys.begin(),activeAngularKeys.end(),
                 [&](const auto &key){return auditAngularCache.contains(key);} )){
      throw std::runtime_error("A3 shared immutable angular cache changed across frequencies");
    }
    for(const auto &node:nodes){
      const TableKey fKey{node.element,node.node,0,-1};
      const TableKey aSourceKey{node.element,node.node,1,-1};
      const TableKey aReceiverKey{node.element,node.node,1,0};
      const auto fIt=tables.find(fKey);
      const auto aSourceIt=tables.find(aSourceKey);
      const auto aReceiverIt=tables.find(aReceiverKey);
      if(fIt==tables.end()||(aSourceIt==tables.end()&&aReceiverIt==tables.end()))
        throw std::runtime_error("A3 field pairing lacks a forward or weighted-adjoint table");
      const auto material=model(node.layer,node.radius/scales.Length());
      const auto synthesize=[&](const phasea::StrainTable &table,const TableKey &key,
                                std::size_t angleIndex){
        const AuditFrameKey angularKey{frames.at(key),isotropic.at(key),table.orders()};
        const auto &eval=auditAngularCache.at(angularKey).at(angleIndex);
        auto value=phasea::evaluate(table,0,*eval.legendre,0,eval.phi);
        if(!isotropic.at(key))value=phasea::rotateToGlobal(value,eval.a,eval.b);
        return value;
      };
      for(std::size_t ai=0;ai<angles.size();++ai){
        SPARSESPEC::KernelField3DSample sample;sample.fluid=!node.solid;
        const auto forward=synthesize(*fIt->second,fKey,ai);
        phasea::Strain adjointField{};
        for(const auto &adjointKey:{aSourceKey,aReceiverKey}){
          const auto adjointIt=tables.find(adjointKey);
          if(adjointIt==tables.end())continue;
          const auto local=synthesize(*adjointIt->second,adjointKey,ai);
          for(std::size_t c=0;c<phasea::Components;++c)adjointField[c]+=local[c];
        }
        for(const auto value:adjointField)weightedAdjointNorm2+=std::norm(value);
        sample.forwardRawStrain={forward[0],forward[1],forward[2],forward[3],forward[4],forward[5]};
        sample.adjointDualStrain={adjointField[0],adjointField[1],adjointField[2],
            adjointField[3],adjointField[4],adjointField[5]};
        commonKernelValues[static_cast<std::size_t>(&node-nodes.data())*angles.size()+ai]+=
            kernelVsAtLength(scales.Length(),sample,material,dispersion(z));
      }
    }
    reconstructionSeconds+=seconds(reconstructionStart);
  }
  if(!(weightedAdjointNorm2>0.0L))
    throw std::runtime_error("A3 reconstructed weighted-adjoint field has zero aggregate norm");
  double integrated=0.0;
  for(std::size_t n=0;n<nodes.size();++n){
    const double radialWeight=bulkProfile(nodes[n].radius)*nodes[n].radius*nodes[n].radius*nodes[n].dr;
    for(std::size_t ai=0;ai<angles.size();++ai)
      integrated+=commonKernelValues[n*angles.size()+ai]*radialWeight*angles[ai][2];
  }
  const double relativeFD=std::abs(static_cast<double>(finiteDifference)-integrated)/
      std::max(std::abs(static_cast<double>(finiteDifference)),std::abs(integrated));
  const double oldBase=baseMeta.at("integrated_profile_kernel");
  const double commonToFrozen=std::abs(integrated-oldBase)/std::max(std::abs(integrated),std::abs(oldBase));
  if(relativeFD>2e-4)
    throw std::runtime_error("A3 common reconstruction fails the existing 2e-4 central finite-difference gate");
  const double timeIdentity=baseMeta.at("time_shift_adjoint_identity");
  pipeline.stop();
  std::ofstream result(out/"a3-summary.json");
  result<<std::setprecision(17)<<"{\"degree\":3,\"nq\":8,\"refinement\":0.05,\"frequency_count\":35,"
      <<"\"radial_sides\":"<<nodes.size()<<",\"angular_points\":"<<angles.size()
      <<",\"shared_angular_frame_caches\":"<<auditAngularCache.size()
      <<",\"gravity_max_abs_nd_error\":"<<gravityError<<",\"time_shift_adjoint_identity\":"<<timeIdentity
      <<",\"time_shift_adjoint_identity_source\":\"frozen L3 baseline result.json\""
      <<",\"fd_seconds\":"<<static_cast<double>(finiteDifference)
      <<",\"common_integrated_kernel_seconds\":"<<integrated
      <<",\"frozen_ds_integrated_kernel_seconds\":"<<oldBase
      <<",\"common_vs_frozen_ds_relative\":"<<commonToFrozen
      <<",\"fd_relative_error\":"<<relativeFD<<",\"fd_limit\":0.0002"
      <<",\"solve_elapsed_s\":"<<solveSeconds<<",\"reconstruction_elapsed_s\":"<<reconstructionSeconds
      <<",\"source_integrity\":\"same frozen +/- traces and same L3 model/mesh/weights\","
      <<"\"quadrature\":\"Gauss-Legendre(2L+1) colatitude x (4L+1) uniform longitude; native GLL radial weights\","
      <<"\"profile\":\"sin^2 bump on 5300-5800 km; source/receiver support excluded\"}\n";
  return 0;
}

std::vector<std::array<Complex,3>> readComponentSpectrum(const fs::path &path,
                                                         int firstBin,int count) {
  std::ifstream in(path);std::string line;
  if(!in||!std::getline(in,line)||line!="bin,up_re,up_im,north_re,north_im,east_re,east_im")
    throw std::runtime_error("A4 invalid frozen component-spectrum header");
  std::vector<std::array<Complex,3>> out(static_cast<std::size_t>(count));
  std::vector<unsigned char> seen(static_cast<std::size_t>(count));
  while(std::getline(in,line)){
    std::stringstream row(line);std::array<std::string,7> c;
    for(auto &v:c)if(!std::getline(row,v,','))throw std::runtime_error("A4 malformed spectrum row");
    const int bin=std::stoi(c[0]);const int k=bin-firstBin;
    if(k<0||k>=count||seen[static_cast<std::size_t>(k)])throw std::runtime_error("A4 duplicate/outside-bin spectrum row");
    for(int comp=0;comp<3;++comp)out[static_cast<std::size_t>(k)][comp]={
        std::stod(c[1+2*comp]),std::stod(c[2+2*comp])};
    seen[static_cast<std::size_t>(k)]=1;
  }
  if(std::find(seen.begin(),seen.end(),0)!=seen.end())throw std::runtime_error("A4 incomplete frozen component spectrum");
  return out;
}

void writeProfileSnapshot(const fs::path &path,const mathprofile::Snapshot &snapshot,
                          const std::string &route) {
  std::ofstream out(path);if(!out)throw std::runtime_error("cannot create Phase A profile evidence");
  out<<std::setprecision(17)<<"{\"route\":\""<<route<<"\",\"scope_rows\":[";
  bool first=true;
  for(const auto &[name,value]:snapshot.entries){if(!first)out<<',';first=false;
    const bool worker=name.starts_with("standard_worker.")||name.starts_with("frequency_worker.")||name.starts_with("parallel_worker.");
    out<<"{\"name\":\""<<name<<"\",\"execution_domain\":\""<<(worker?"worker_elapsed_sum":"master_wall")<<"\",\"inclusive_s\":"<<value.inclusive*1e-9
       <<",\"exclusive_s\":"<<value.exclusive*1e-9<<",\"calls\":"<<value.calls<<"}";}
  out<<"],\"counters\":[";first=true;
  for(const auto &[name,value]:snapshot.counters){if(!first)out<<',';first=false;
    out<<"{\"name\":\""<<name<<"\",\"value\":"<<value<<"}";}
  out<<"]}\n";
}

int runA4(const fs::path &input,const fs::path &out,const fs::path &referenceComponents,
          int degree,int nq,double refinement,int frequencyWorkers) {
  if(fs::exists(out)&&fs::directory_iterator(out)!=fs::directory_iterator{})
    throw std::runtime_error("A4 output directory must be empty");
  fs::create_directories(out);Eigen::setNbThreads(1);omp_set_dynamic(0);omp_set_num_threads(1);
  mathprofile::Root root("phase-a-forward-only-control",-1);
  mathprofile::Scope pipeline("phase_a_a4_forward_only_pipeline");
  const auto in=YSpec::ReadKernelInput(input);auto model=YSpec::Model::ReadDeck(in.model);
  CommonModel common{model,0.0};const auto &scales=model.GetScales();
  auto grid=YSpec::MakeGrid(in.calculation,scales);
  auto geoLat=[](double lat){return 90-YSpec::Sphere::Colatitude(lat)*180/pi;};
  const auto &moment=in.source.moment;
  SourceInfo::EarthquakeCMT source(geoLat(in.source.location.latitude),in.source.location.longitude,
      in.source.depth/1000,moment.rr,moment.rt,moment.rp,moment.tt,moment.tp,moment.pp);
  Full1D::SEM sem(common,refinement,nq,degree,source);
  {std::ofstream pf(out/"parameters.txt");pf<<params(in,degree);}
  InputParameters parameters((out/"parameters.txt").string());
  if(parameters.num_receivers()!=1||grid.Count()!=35||grid.First()!=7||grid.Last()!=41)
    throw std::runtime_error("A4 case/receiver/frequency grid differs from the frozen contract");
  const auto reference=readComponentSpectrum(referenceComponents,grid.First(),static_cast<int>(grid.Count()));
  const int receiverElement=sem.receiverElements(parameters,true).front();
  const double receiverRadius=scales.Length()*model.SurfaceRadius()-in.receiver.depth;
  std::vector<SPARSESPEC::Field3DPoint> receiver{{receiverRadius,
      YSpec::Sphere::Colatitude(in.receiver.location.latitude),
      in.receiver.location.longitude*pi/180,receiverElement}};
  auto dispersion=[&](Complex z){return in.calculation.attenuation?
      2/pi*std::log(Complex(0,1)*z/model.ReferenceFrequency()):Complex{};};
  dssolver::BatchResponseCache cache({static_cast<std::size_t>(2*degree)*grid.Count()+1,16ull<<30,1e-12});
  std::vector<std::array<Complex,3>> components(grid.Count());
  const int workerCount=std::min<int>(frequencyWorkers,static_cast<int>(grid.Count()));
  if(workerCount<1||workerCount>8)throw std::runtime_error("A4 frequency worker bound exceeded");
  std::vector<mathprofile::WorkerData> workerData(static_cast<std::size_t>(workerCount));
  const int diagnosticBackend=SPARSESPEC::FixedMeshAdjointSolver::diagnosticBackend();
  const auto solveStart=Clock::now();
  dspec::frequency::runIndexedWithContext(grid.Count(),workerCount,grid.First(),"phase_a_a4",
    [&](int worker) noexcept {return FrequencyWorkerContext(true,workerData[static_cast<std::size_t>(worker)],diagnosticBackend);},
    [&](std::size_t k,auto &){
      auto z=grid.Omega(k);
      auto workspace=SPARSESPEC::Windowed3DFieldEngine::prepareFrequencyWorkspace(
          sem,source,parameters,z,Eigen::Vector3cd::Zero(),&cache,true,true,nullptr,1,dispersion(z));
      const auto field=workspace.sampleForKernelsPacked(receiver,0,nullptr,2).front();
      for(int c=0;c<3;++c)components[k][static_cast<std::size_t>(c)]=
          scales.Length()*field.forwardRawDisplacement[static_cast<std::size_t>(c)];
    });
  for(const auto &data:workerData)mathprofile::mergeWorker(data,"frequency_worker.");
  const double forwardElapsed=seconds(solveStart);
  double error2=0.0,norm2=0.0,maxAbs=0.0;
  std::ofstream csv(out/"a4-components.csv");
  csv<<"bin,up_re,up_im,north_re,north_im,east_re,east_im\n"<<std::setprecision(17);
  for(std::size_t k=0;k<grid.Count();++k){
    csv<<grid.First()+k;
    for(std::size_t c=0;c<3;++c){
      const Complex diff=components[k][c]-reference[k][c];error2+=std::norm(diff);
      norm2+=std::norm(reference[k][c]);maxAbs=std::max(maxAbs,std::abs(diff));
      csv<<','<<components[k][c].real()<<','<<components[k][c].imag();
    }
    csv<<'\n';
  }
  const double relative=std::sqrt(error2/std::max(norm2,1e-300));
  if(relative>1e-11)throw std::runtime_error("A4 forward-only receiver spectrum parity gate failed");
  pipeline.stop();
  const auto profile=mathprofile::snapshot();
  const auto count=[&](const std::string &name){auto it=profile.counters.find(name);return it==profile.counters.end()?std::uint64_t{0}:it->second;};
  if(count("true_forward_only_route")!=grid.Count()||count("solver_forward_rhs_columns")==0||
     count("solver_unit_basis_rhs_columns")!=0||count("solver_weighted_adjoint_rhs_columns")!=0||
     count("full_weighted_modal_adjoint_rhs_columns")!=0)
    throw std::runtime_error("A4 forward-only profile counters do not match the required route");
  for(const auto &[name,entry]:profile.entries){(void)entry;
    if(name.find("minimal_unit_basis_HH_solve")!=std::string::npos||
       name.find("full_weighted_HH_solve")!=std::string::npos||
       name.find("receiver_basis_and_QR_reconstruction")!=std::string::npos||
       name.find("receiver_basis_assembly_and_qr")!=std::string::npos)
      throw std::runtime_error("A4 forward-only path entered an adjoint/receiver-QR phase");
  }
  writeProfileSnapshot(out/"a4-profile.json",profile,"forwardOnly=true,poleOnly=true");
  std::ofstream result(out/"a4-summary.json");
  result<<std::setprecision(17)<<"{\"bins\":35,\"workers\":"<<workerCount
      <<",\"forward_only_route_calls\":"<<count("true_forward_only_route")
      <<",\"forward_elapsed_s_instrumented\":"<<forwardElapsed
      <<",\"receiver_relative_l2\":"<<relative<<",\"receiver_max_abs\":"<<maxAbs
      <<",\"unit_hh_rhs_columns\":"<<count("solver_unit_basis_rhs_columns")
      <<",\"weighted_hh_rhs_columns\":"<<count("solver_weighted_adjoint_rhs_columns")
      <<",\"receiver_qr_scopes\":0}\n";

  return 0;
}

int runA5(const fs::path &input,const fs::path &out,const fs::path &referenceComponents,
          int degree,int nq,double refinement,int workerCount) {
  if(fs::exists(out)&&fs::directory_iterator(out)!=fs::directory_iterator{})
    throw std::runtime_error("A5 output directory must be empty");
  if(workerCount<1||workerCount>8)throw std::runtime_error("A5 worker bound exceeded");
  fs::create_directories(out);Eigen::setNbThreads(1);omp_set_dynamic(0);
  omp_set_num_threads(workerCount);
  mathprofile::Root root("phase-a-standard-forward-batch",-1);
  mathprofile::Scope pipeline("phase_a_a5_standard_forward_batch");
  const auto in=YSpec::ReadKernelInput(input);auto model=YSpec::Model::ReadDeck(in.model);
  CommonModel common{model,0.0};const auto &scales=model.GetScales();
  const auto grid=YSpec::MakeGrid(in.calculation,scales);
  auto geoLat=[](double lat){return 90-YSpec::Sphere::Colatitude(lat)*180/pi;};
  const auto &moment=in.source.moment;
  SourceInfo::EarthquakeCMT source(geoLat(in.source.location.latitude),in.source.location.longitude,
      in.source.depth/1000,moment.rr,moment.rt,moment.rp,moment.tt,moment.tp,moment.pp);
  Full1D::SEM sem(common,refinement,nq,degree,source);
  {std::ofstream pf(out/"parameters.txt");pf<<params(in,degree);}
  InputParameters legacy((out/"parameters.txt").string());
  SpectraSolver::FreqFull freq(legacy.f1(),legacy.f2(),legacy.f11(),legacy.f12(),
      legacy.f21(),legacy.f22(),legacy.time_step_sec(),legacy.t_out()/60.0,
      1.0,0.1,0.0,legacy.t_out()/60.0,1,scales.Time());
  if(legacy.num_receivers()!=1||freq.i1()!=7||freq.i2()!=42||
     degree!=200||nq!=10||refinement!=0.025)
    throw std::runtime_error("A5 case, receiver, mesh, or frequency grid differs from the frozen contract");
  const auto reference=readComponentSpectrum(referenceComponents,7,35);
  if(!sem.sourceMatches(source))throw std::runtime_error("A5 standard request source differs from the common case");
  if(std::abs(legacy.receiver_depth()*1000.0-in.receiver.depth)>1e-10)
    throw std::runtime_error("A5 standard request receiver depth differs from the common case");
  const auto receiverCoordinates=legacy.receivers().at(0);
  if(receiverCoordinates.first!=geoLat(in.receiver.location.latitude)||
     receiverCoordinates.second!=in.receiver.location.longitude)
    throw std::runtime_error("A5 standard request receiver coordinates differ from the common case");
  std::vector<Complex> exactZ;exactZ.reserve(grid.Count());
  for(std::size_t k=0;k<grid.Count();++k)exactZ.push_back(grid.Omega(k));
  const double referencePeriodND=2.0*pi/model.ReferenceFrequency();
  constexpr int standardNskip=10;
  SPARSESPEC::SpectraRunContext request(freq,source,legacy,referencePeriodND,standardNskip);
  SPARSESPEC::SpecConstants standardConstants(freq.ep(),referencePeriodND);
  if(standardConstants.ieps.imag()!=exactZ.front().imag())
    throw std::runtime_error("A5 standard/common damping differs in the shared time scale");
  auto writeSpectrum=[&](const std::string &stem,const Eigen::MatrixXcd &matrix){
    if(matrix.rows()!=3||matrix.cols()<42)throw std::runtime_error("A5 standard spectrum shape mismatch");
    std::ofstream csv(out/(stem+".csv"));
    csv<<"bin,up_re,up_im,north_re,north_im,east_re,east_im\n"<<std::setprecision(17);
    double error2=0.0,norm2=0.0,maxAbs=0.0;
    for(int bin=7;bin<=41;++bin){
      csv<<bin;
      for(int c=0;c<3;++c){
        const Complex z=grid.Omega(static_cast<std::size_t>(bin-7));
        const Complex outputMultiplier=SPARSESPEC::outputFactor(legacy.output_type(),z,standardConstants.myi);
        if(outputMultiplier==Complex{})throw std::runtime_error("A5 standard output conversion multiplier is zero");
        const Complex value=scales.Length()*matrix(c,bin)/outputMultiplier;
        const Complex difference=value-reference[static_cast<std::size_t>(bin-7)][static_cast<std::size_t>(c)];
        error2+=std::norm(difference);norm2+=std::norm(reference[static_cast<std::size_t>(bin-7)][static_cast<std::size_t>(c)]);
        maxAbs=std::max(maxAbs,std::abs(difference));
        csv<<','<<value.real()<<','<<value.imag();
      }
      csv<<'\n';
    }
    return std::array<double,2>{std::sqrt(error2/std::max(norm2,1e-300)),maxAbs};
  };
  auto previousProfile=mathprofile::snapshot();
  std::vector<SPARSESPEC::PhaseACutIndexRecord> untruncatedCutIndices;
  const auto untruncatedStart=Clock::now();
  Eigen::MatrixXcd untruncated;
  {mathprofile::Scope submode("phase_a_a5_untruncated_standard_forward");
    untruncated=SPARSESPEC::phaseAStandardForward(request,sem,false,true,exactZ,&untruncatedCutIndices);}
  const double untruncatedSeconds=seconds(untruncatedStart);
  const auto untruncatedSnapshot=mathprofile::snapshot();
  const auto untruncatedProfile=mathprofile::difference(untruncatedSnapshot,previousProfile);
  previousProfile=untruncatedSnapshot;
  const auto untruncatedError=writeSpectrum("a5-untruncated-components",untruncated);
  if(untruncatedError[0]>1e-11)
    throw std::runtime_error("A5 untruncated standard forward disagrees with A4 receiver spectrum");
  const auto truncatedStart=Clock::now();
  Eigen::MatrixXcd truncated;
  std::vector<SPARSESPEC::PhaseACutIndexRecord> truncatedCutIndices;
  {mathprofile::Scope submode("phase_a_a5_existing_truncated_standard_forward");
    truncated=SPARSESPEC::phaseAStandardForward(request,sem,false,false,exactZ,&truncatedCutIndices);}
  const double truncatedSeconds=seconds(truncatedStart);
  const auto truncatedSnapshot=mathprofile::snapshot();
  const auto truncatedProfile=mathprofile::difference(truncatedSnapshot,previousProfile);
  const auto truncatedError=writeSpectrum("a5-truncated-components",truncated);
  const auto profileCounter=[](const mathprofile::Snapshot &snapshot,const std::string &name){
    const auto found=snapshot.counters.find(name);return found==snapshot.counters.end()?std::uint64_t{0}:found->second;
  };
  const auto untruncatedFactors=profileCounter(untruncatedProfile,"sparse_factorization_calls");
  const auto untruncatedForwardRhs=profileCounter(untruncatedProfile,"forward_rhs_columns");
  const auto untruncatedAdjointRhs=profileCounter(untruncatedProfile,"adjoint_rhs_columns");
  const auto truncatedFactors=profileCounter(truncatedProfile,"sparse_factorization_calls");
  const auto truncatedForwardRhs=profileCounter(truncatedProfile,"forward_rhs_columns");
  const auto truncatedAdjointRhs=profileCounter(truncatedProfile,"adjoint_rhs_columns");
  if(untruncatedFactors!=14035||untruncatedForwardRhs!=42035||untruncatedAdjointRhs!=0||
     truncatedFactors!=14035||truncatedForwardRhs!=42035||truncatedAdjointRhs!=0)
    throw std::runtime_error("A5 standard forward counters do not cover all factors/RHS or include an adjoint solve");
  const auto writeCutIndices=[&](const fs::path &path,
      const std::vector<SPARSESPEC::PhaseACutIndexRecord> &records,bool isUntruncated){
    std::ofstream csv(path);if(!csv)throw std::runtime_error("cannot create A5 cut-index audit");
    csv<<"family,degree,bin,formulation,start_dof,retained_rows,full_rows,nskip,policy\n";
    for(const auto &record:records){
      if(record.family<0)continue;
      const char *family=record.family==0?"radial":record.family==1?"toroidal":"spheroidal";
      const char *formulation=record.formulation==1?"cowling":record.formulation==2?"full":"standard";
      csv<<family<<','<<record.degree<<','<<record.bin<<','<<formulation<<','
         <<record.startDof<<','<<record.retainedRows<<','<<record.fullRows<<','
         <<record.nskip<<','<<(isUntruncated?"untruncated":"source_aware")<<'\n';
    }
  };
  writeCutIndices(out/"a5-untruncated-cut-indices.csv",untruncatedCutIndices,true);
  writeCutIndices(out/"a5-truncated-cut-indices.csv",truncatedCutIndices,false);
  const auto activeCutRows=[](const auto &records){return std::count_if(records.begin(),records.end(),
      [](const auto &record){return record.family>=0;});};
  pipeline.stop();
  writeProfileSnapshot(out/"a5-untruncated-profile.json",untruncatedProfile,
                       "pinned_standard_forward_source_start_zero; OpenMP degree-parallel worker domain");
  writeProfileSnapshot(out/"a5-truncated-profile.json",truncatedProfile,
                       "pinned_standard_forward_existing_source_aware_truncation; OpenMP degree-parallel worker domain");
  std::ofstream summary(out/"a5-summary.json");
  summary<<std::setprecision(17)<<"{\"workers\":"<<workerCount
      <<",\"worker_domain\":\"standard_degree_parallel\",\"bins\":35,\"untruncated_instrumented_elapsed_s\":"<<untruncatedSeconds
      <<",\"nskip\":"<<standardNskip<<",\"prepared_single_sem_policy\":true"
      <<",\"untruncated_receiver_relative_l2\":"<<untruncatedError[0]
      <<",\"untruncated_receiver_max_abs\":"<<untruncatedError[1]
      <<",\"truncated_instrumented_elapsed_s\":"<<truncatedSeconds
      <<",\"truncated_receiver_relative_l2\":"<<truncatedError[0]
      <<",\"truncated_receiver_max_abs\":"<<truncatedError[1]
      <<",\"untruncated_factorizations\":"<<untruncatedFactors
      <<",\"untruncated_forward_rhs_columns\":"<<untruncatedForwardRhs
      <<",\"untruncated_adjoint_rhs_columns\":"<<untruncatedAdjointRhs
      <<",\"truncated_factorizations\":"<<truncatedFactors
      <<",\"truncated_forward_rhs_columns\":"<<truncatedForwardRhs
      <<",\"truncated_adjoint_rhs_columns\":"<<truncatedAdjointRhs
      <<",\"untruncated_cut_index_rows\":"<<activeCutRows(untruncatedCutIndices)
      <<",\"truncated_cut_index_rows\":"<<activeCutRows(truncatedCutIndices)
      <<",\"instrumented_elapsed_includes_cut_index_capture_and_allocation\":true"
      <<",\"instrumented_elapsed_includes_cut_index_csv_serialization\":false"
      <<",\"same_mesh_and_forcing\":true,\"full_pipeline_speed_claim\":false}\n";
  return 0;
}
}

int main(int argc,char **argv) try {
  if(argc>=2&&std::string(argv[1])=="a1"){
    if(argc!=10)throw std::runtime_error("usage: phase-a-driver a1 INPUT OUTPUT L NQ REFINEMENT SOLVE_WORKERS FROZEN_TRACE FROZEN_RESULT");
    return runA1(argv[2],argv[3],std::stoi(argv[4]),std::stoi(argv[5]),
        std::stod(argv[6]),std::stoi(argv[7]),argv[8],argv[9]);
  }
  if(argc>=2&&std::string(argv[1])=="a2"){
    if(argc!=10)throw std::runtime_error("usage: phase-a-driver a2 INPUT OUTPUT R9 REF9 R25 REF25 R40 REF40");
    return runA2(argv[2],argv[3],{argv[4],argv[6],argv[8]},{argv[5],argv[7],argv[9]});
  }
  if(argc>=2&&std::string(argv[1])=="a3"){
    if(argc!=14)throw std::runtime_error("usage: phase-a-driver a3 INPUT OUTPUT BASE_TRACE BASE_RESULT PLUS_TRACE PLUS_RESULT MINUS_TRACE MINUS_RESULT L NQ REFINEMENT SOLVE_WORKERS");
    return runA3(argv[2],argv[3],argv[4],argv[5],argv[6],argv[7],argv[8],argv[9],
        std::stoi(argv[10]),std::stoi(argv[11]),std::stod(argv[12]),std::stoi(argv[13]));
  }
  if(argc>=2&&std::string(argv[1])=="a4"){
    if(argc!=9)throw std::runtime_error("usage: phase-a-driver a4 INPUT OUTPUT REFERENCE_COMPONENTS L NQ REFINEMENT FREQUENCY_WORKERS");
    return runA4(argv[2],argv[3],argv[4],std::stoi(argv[5]),std::stoi(argv[6]),
        std::stod(argv[7]),std::stoi(argv[8]));
  }
  if(argc>=2&&std::string(argv[1])=="a5"){
    if(argc!=9)throw std::runtime_error("usage: phase-a-driver a5 INPUT OUTPUT REFERENCE_COMPONENTS L NQ REFINEMENT WORKERS");
    return runA5(argv[2],argv[3],argv[4],std::stoi(argv[5]),std::stoi(argv[6]),
        std::stod(argv[7]),std::stoi(argv[8]));
  }
  throw std::runtime_error("unknown Phase A mode");
} catch(const std::exception &e) {
  std::cerr<<"FAIL "<<e.what()<<'\n';return 1;
}
