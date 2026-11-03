"""Exclusive mathematical-stage accounting; preserve native and instrumented boundaries."""
from pathlib import Path
import csv,json,statistics as st
R=Path(__file__).resolve().parent
D={
'measurement_adjoint':'measurement',
'dynamic_H_assembly':'operator',
'sparse_equilibration_and_scaling':'equilibration',
'sparse_symbolic_analysis':'symbolic',
'sparse_numerical_factorization':'factorization',
'sparse_forward_solve_including_rhs_scaling':'forward_substitution',
'forward_substitution':'forward_substitution',
'sparse_adjoint_solve_including_rhs_scaling':'adjoint_substitution',
'full_weighted_HH_solve':'adjoint_substitution',
'minimal_unit_basis_HH_solve':'adjoint_substitution',
'receiver_basis_and_QR_reconstruction':'source_receiver_recovery',
'batch_weighted_adjoint_reconstruction':'source_receiver_recovery',
'spheroidal_source_lift_projection_and_solver':'source_receiver_recovery',
'toroidal_source_lift_projection_and_solver':'source_receiver_recovery',
'radial_l0_source_and_adjoint_branch':'source_receiver_recovery',
'shared_forward_and_adjoint_solver':'source_receiver_recovery',
'modal_projection_source_lift_and_frequency_solve':'source_receiver_recovery',
'frequency_workspace_prepare_and_modal_capture':'source_receiver_recovery',
'frequency_workspace_and_modal_solve':'source_receiver_recovery',
'batch_response_cache_get_or_compute':'response_cache',
'solver_factorization_dispatch':'factor_overhead',
'sparse_factorization_total':'factor_overhead',
'radial_interpolation_and_strain_numerator_build':'radial_strain',
'radial_l0_sampling_and_interpolation':'radial_strain',
'packed_Wigner_d_table_build':'angular_plan_basis',
'packed_angular_coefficients_and_basis_build':'angular_plan_basis',
'source_frame_geometry_phase_spin_rotation_and_scatter':'representative_receiver_sampling',
'receiver_frame_geometry_phase_spin_rotation_and_scatter':'representative_receiver_sampling',
'pole_frame_sampling_total':'representative_receiver_sampling',
'actual_ZGEMM_complex':'angular_gemm',
'actual_DGEMM_real_part':'angular_gemm',
'actual_DGEMM_imag_part':'angular_gemm',
'tensor_angular_gemm_complex_and_packing':'gemm_buffers_copies',
'tensor_angular_gemm_real_and_packing':'gemm_buffers_copies',
'tensor_coefficient_packing_and_key_reorder':'coefficient_packing',
'tensor_phase_spin_rotation_and_scatter':'rotation_scatter',
'tensor_grid_semantic_key_build_and_lookup':'grid_key_lookup',
'tensor_grid_cold_layout_map_build':'grid_cold_layout',
'tensor_representative_and_coefficient_record_capture':'representative_receiver_sampling',
'tensor_result_allocation_and_volume_tagging':'result_allocation',
'tensor_angular_plan_key_build_and_lookup':'angular_plan_basis',
'tensor_cold_frame_geometry_and_wigner_build':'angular_plan_basis',
'tensor_complex_basis_key_build_and_lookup':'angular_plan_basis',
 'tensor_complex_basis_key_construction':'angular_plan_basis',
 'tensor_complex_basis_cache_lookup':'angular_plan_basis',
 'tensor_real_basis_key_construction':'angular_plan_basis',
 'tensor_real_basis_cache_lookup':'angular_plan_basis',
'tensor_complex_basis_cold_fill_and_cache':'angular_plan_basis',
'tensor_real_basis_key_build_and_lookup':'angular_plan_basis',
'tensor_real_basis_cold_fill_and_cache':'angular_plan_basis',
'tensor_radial_l0_scatter':'l0_scatter',
'tensor_forward_cache_copy_and_receiver_suffix':'forward_cache',
'tensor_receiver_suffix_sampling_and_scatter':'representative_receiver_sampling',
'forward_angular_cache_lookup':'forward_cache',
'forward_angular_cache_store':'forward_cache',
'receiver_suffix_and_frequency_scatter':'receiver_scatter',
'point_kernel_contraction_and_scatter':'kernel_products',
'source_contact_separate_case':'source_contacts',
'Pforward_fft_and_window':'Pforward',
'final_frequency_objective_accumulation':'objective',
'tensor_sampling_total':'tensor_misc_cleanup',
'tensor_grid_semantic_plan_and_sampling':'tensor_misc_cleanup',
'candidate_body':'api_misc_cleanup',
'candidate_api_including_local_cleanup':'api_misc_cleanup',
}
Y={
'radial.initial_conditions_and_upward_integration':'ode_up_integration',
'radial.initial_conditions_and_downward_integration':'ode_down_integration',
'radial.spheroidal_b_vector_transport':'b_vector_transport',
'radial.determinant_jumps_matching_and_response_transport':'matching_transport',
'radial.radial_physical_variable_recovery':'radial_recovery',
'radial.toroidal_physical_variable_recovery':'radial_recovery',
'radial.spheroidal_physical_variable_recovery':'radial_recovery',
'radial.coefficient_combine_and_strain_recovery':'radial_strain',
'kernel.angular_harmonic_synthesis.source':'angular_source',
'kernel.angular_harmonic_synthesis.adjoint':'angular_adjoint',
'kernel.relative_spin_rotation':'rotation',
'kernel.elastic_energy_pair_contraction':'energy_pairs',
'kernel.objective_weight_and_parameter_products':'parameter_products',
'kernel.frequency_block_accumulation_copy':'frequency_accumulation',
'kernel.angular_block_batch':'angular_block_misc_cleanup',
'kernel.measurement_frequency_basis_weights':'measurement',
'kernel.modulus_dispersion_frequency_weights':'modulus_weights',
'output.parameter_transform_scaling_and_construction':'output_scaling',
'api.ComputeSampledKernels':'api_misc_cleanup',
}
def bucket(solver,name):
 if solver=='dspec':return D[name]
 if name in Y:return Y[name]
 if name.startswith('setup.'):return 'setup_api'
 if name.startswith('angular.'):return 'angular_tables'
 if name.startswith('source_weights.'):return 'source_receiver_weights'
 if name.startswith('radial.'):return 'radial_driver_overhead'
 raise KeyError(name)
def api_profiles(r):
 return [p for p in r['profiles'] if p['accounting_root'] in ['candidate_api_including_local_cleanup','api.ComputeSampledKernels']]
def api_native(r):
 m=r['metrics'];return sum(m['candidate_call_wall_s']) if r['solver']=='dspec' else sum(m['api_s'])
def aggregate_profiles(ps,solver):
 stages={};counters={};groups={}
 for p in ps:
  for n,v in p['stages'].items():
   dest=stages.setdefault(n,{'inclusive_s':0.,'exclusive_s':0.,'calls':0})
   for k in dest:dest[k]+=v[k]
   g=bucket(solver,n);groups[g]=groups.get(g,0)+v['exclusive_s']
  for n,v in p['counters'].items():counters[n]=counters.get(n,0)+v
 total=sum(p['stages'][p['accounting_root']]['inclusive_s'] for p in ps)
 assert abs(sum(groups.values())-total)<1e-8
 return {'api_s':total,'groups':groups,'stages':stages,'counters':counters}
all_records=json.loads((R/'runs-01/records.json').read_text())
if (R/'contacts-01/records.json').exists():all_records+=json.loads((R/'contacts-01/records.json').read_text())
summaries=[];stage_rows=[]
keys=sorted({(r['degree'],r['bins'],r['objectives'],r['contacts'],r['solver']) for r in all_records})
for l,f,m,c,solver in keys:
 runs=[r for r in all_records if (r['degree'],r['bins'],r['objectives'],r['contacts'],r['solver'])==(l,f,m,c,solver)]
 native=[r for r in runs if r['mode']=='native'];profile=[r for r in runs if r['mode']=='profile'];original=[r for r in runs if r['mode']=='original']
 assert len(native)==len(profile)==3 and len(original)==1
 native_total=[api_native(r) for r in native]
 prof_total=[aggregate_profiles(api_profiles(r),solver)['api_s'] for r in profile]
 representative=sorted(zip(prof_total,profile),key=lambda x:x[0])[1][1]
 ag=aggregate_profiles(api_profiles(representative),solver)
 def med(fn,rs=native):return st.median(fn(r) for r in rs)
 s={'degree':l,'bins':f,'objectives':m,'contacts':c,'solver':solver,'processes_native':3,'processes_profiled':3,
 'native_api_median_s':st.median(native_total),'native_api_range_s':[min(native_total),max(native_total)],
 'profile_api_median_s':st.median(prof_total),'profile_api_range_s':[min(prof_total),max(prof_total)],
 'profile_over_native_percent':100*(st.median(prof_total)/st.median(native_total)-1),
 'original_api_single_s':sum(original[0]['metrics'].get('candidate_s',original[0]['metrics'].get('api_s'))),
 'native_peak_hwm_bytes':max(r['peak_current_executable_hwm_bytes'] for r in native),
 'representative_repeat':representative['repeat'],'representative_data':representative['data'],
 'representative_profile':ag,'representative_per_objective':[aggregate_profiles([p],solver) for p in api_profiles(representative)],
 'setup_native_median_s':med(lambda r:r['metrics']['main_setup_s' if solver=='dspec' else 'setup_s']),
 'model_to_arrays_native_median_s':med(lambda r:r['metrics']['main_setup_s' if solver=='dspec' else 'setup_s']+api_native(r)),
 'all_original_equivalence_pass':all(r['original_equivalence']['passed'] for r in native+profile),
 'all_clone_output_byte_identical':all(r['output_sha256']==original[0]['output_sha256'] for r in native+profile),
 }
 s['setup_profile']=[p for p in representative['profiles'] if p not in api_profiles(representative)]
 if solver=='dspec':
  s['binary_export_native_median_s']=med(lambda r:sum(r['metrics']['binary_array_export_s']))
  s['reference_benchmark_only_native_median_s']=med(lambda r:sum(r['metrics']['reference_s']))
  s['validation_benchmark_only_native_median_s']=med(lambda r:sum(r['metrics']['reference_validation_compare_s']))
  s['native_per_objective_median_s']=[st.median(r['metrics']['candidate_call_wall_s'][i] for r in native) for i in range(m)]
 else:
  s['CSV_export_native_median_s']=med(lambda r:r['metrics']['end_to_end_s']-r['metrics']['first_model_to_arrays_s'])
 summaries.append(s)
 for r in profile:
  for p in r['profiles']:
    for name,v in p['stages'].items():stage_rows.append({'degree':l,'bins':f,'objectives':m,'contacts':c,'solver':solver,'instrumentation':'fine','repeat':r['repeat'],'profile_root':p['accounting_root'],'objective':p['objective'],'stage':name,**v})
coarse=json.loads((R/'coarse-runs-01/records.json').read_text())
for s in summaries:
 if s['solver']!='yspec':continue
 cr=[r for r in coarse if (r['degree'],r['bins'],r['objectives'])==(s['degree'],s['bins'],s['objectives'])]
 assert len(cr)==3 and all(r['original_equivalence']['passed'] for r in cr)
 ct=[aggregate_profiles(api_profiles(r),'yspec')['api_s'] for r in cr]
 representative=sorted(zip(ct,cr),key=lambda x:x[0])[1][1]
 s['coarse_api_median_s']=st.median(ct);s['coarse_api_range_s']=[min(ct),max(ct)]
 s['coarse_over_native_percent']=100*(st.median(ct)/s['native_api_median_s']-1)
 s['coarse_representative_repeat']=representative['repeat']
 s['coarse_representative_profile']=aggregate_profiles(api_profiles(representative),'yspec')
 s['coarse_all_byte_identical']=all(r['output_sha256']==next(x['output_sha256'] for x in all_records if x['solver']=='yspec' and x['mode']=='original' and (x['degree'],x['bins'],x['objectives'])==(s['degree'],s['bins'],s['objectives'])) for r in cr)
 for r in cr:
  for p in r['profiles']:
   for name,v in p['stages'].items():stage_rows.append({'degree':s['degree'],'bins':s['bins'],'objectives':s['objectives'],'contacts':0,'solver':'yspec','instrumentation':'coarse','repeat':r['repeat'],'profile_root':p['accounting_root'],'objective':p['objective'],'stage':name,**v})
(R/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
with (R/'all-stage-times.csv').open('w') as h:
 w=csv.DictWriter(h,fieldnames=list(stage_rows[0]));w.writeheader();w.writerows(stage_rows)
print('rows',len(summaries),'stages',len(stage_rows),'all outputs pass',all(x['all_original_equivalence_pass'] for x in summaries))
for s in summaries:print(s['solver'],s['degree'],s['bins'],s['objectives'],s['contacts'],round(s['native_api_median_s'],6),round(s['profile_api_median_s'],6),round(s['profile_over_native_percent'],2))
