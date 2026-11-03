#!/usr/bin/env python3
"""Summarise a completed, hash-verified parallel-solve campaign snapshot."""
from pathlib import Path
import argparse,json,hashlib
p=argparse.ArgumentParser();p.add_argument('snapshot',type=Path);p.add_argument('controls',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
s=json.loads((a.snapshot/'state.json').read_text());assert s['status']=='completed'
rows=[];checks=[]
for job in s['jobs']:
 assert job['status']=='completed' and job['returncode']==0 and job['inputs_unchanged'] and job['cleanup']['verified_gone']
 for rel,digest in job['outputs_sha256'].items():
  f=a.snapshot/job['id']/rel
  if f.exists():assert hashlib.sha256(f.read_bytes()).hexdigest()==digest,str(f)
 f=a.snapshot/job['id']/'checks.json'
 if f.exists():
  c=json.loads(f.read_text());assert c['status']=='PASS';checks.append({'job':job['id'],'status':c['status'],'byte_exact_files':c.get('byte_exact_files')});continue
 f=a.snapshot/job['id']/'data/result.json'
 if not f.exists():continue
 result=json.loads(f.read_text());log=(a.snapshot/job['id']/job['log']).read_text();profiles=[json.loads(x.split('MATHPROFILE ',1)[1]) for x in log.splitlines() if x.startswith('MATHPROFILE ')]
 row={'job':job['id'],'degree':result['degree'],'solve_threads_requested':result['solve_threads_requested'],'mode':result['mode'],'response_reuse':result['response_reuse_used'],'response_spool':result['response_record_spool_used'],'elapsed_s':result['elapsed_s'],'energy_m2_s':result['energy_m2_s'],'peak_observed_threads_including_supervisor':job['peak_observed_threads'],'peak_observed_rss_gib_userspace':job['peak_observed_rss_gib_userspace'],'factorizations':None}
 if profiles:
  assert len(profiles)==1
  profile=profiles[0];stages=profile['stages'];counters=profile['counters'];n=result['degree'];bins=result['frequency_bins'];total_modes=2*n+1
  expected=total_modes*(bins+1) if row['response_spool'] else 2*n*(bins+1)+2*bins+1 if row['response_reuse'] else total_modes*(2*bins+1)
  assert counters['sparse_factorization_calls']==expected,(job['id'],counters['sparse_factorization_calls'],expected)
  assert counters['sparse_symbolic_analysis_calls']+counters['sparse_symbolic_analysis_reuses']==expected
  row['factorizations']=expected;row['symbolic_analyses']=counters['sparse_symbolic_analysis_calls']
  for name,key in [('modal_solve_wall_s','modal_projection_source_lift_and_frequency_solve'),('baseline_capture_wall_s','baseline_forward_with_bounded_response_capture'),('spool_capture_wall_s','baseline_forward_and_compact_response_capture'),('kernel_workspace_wall_s','frequency_workspace_modal_solve_and_adjoint_construction')]:row[name]=stages.get(key,{}).get('inclusive_s')
  row['worker_stage_durations_are_summed_elapsed_not_wall']=True
  row['worker_components']={key:{'inclusive_s':sum(stages.get(prefix+key,{}).get('inclusive_s',0) for prefix in ('','parallel_worker.')),'calls':sum(stages.get(prefix+key,{}).get('calls',0) for prefix in ('','parallel_worker.'))} for key in ('dynamic_H_assembly','solver_factorization_dispatch','forward_substitution','minimal_unit_basis_HH_solve','full_weighted_HH_solve')}
 rows.append(row)
controls=json.loads((a.controls/'runtime-controls/audit/controls.json').read_text());cc={x['name']:x for x in controls['cases']};assert cc['serial-all']['returncode']==0 and cc['minimal-openmp']['race_reported'] and cc['serial-solver-parallel-setup']['race_reported']
report={'status':'PASS_scientific_equivalence_and_cache_stress','full_process_tsan_status':'INCONCLUSIVE_OPENMP_RUNTIME_SYNCHRONIZATION','full_process_tsan_clean_not_claimed':True,'minimal_openmp_control_race_reported':True,'serial_solver_parallel_setup_reproduces_warning':True,'serial_all_tsan_passes':True,'thread_cap':50,'campaign_watchdog_thread_cap':s.get('thread_cap',40),'resource_enforcement':'USER_SPACE_SESSION_WATCHDOG','rows':rows,'checks':checks,'snapshot_verification':json.loads((a.snapshot/'snapshot-verification.json').read_text())}
with a.output.open('x') as f:json.dump(report,f,indent=2,allow_nan=False)
print(json.dumps({'rows':len(rows),'checks':len(checks),'max_observed_threads':max(r['peak_observed_threads_including_supervisor'] for r in rows)}))
