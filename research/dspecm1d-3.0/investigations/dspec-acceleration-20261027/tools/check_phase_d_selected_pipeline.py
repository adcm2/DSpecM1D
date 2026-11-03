#!/usr/bin/env python3
"""Check matched selected-mesh scalar/batched DSpec and one YSpec timing run."""
from __future__ import annotations
import argparse,csv,json,math,pathlib,struct

def rows(p):
    with pathlib.Path(p).open(newline='') as f:return list(csv.DictReader(f))
def obj(p):return json.loads(pathlib.Path(p).read_text())
def need(x,msg):
    if not x:raise ValueError(msg)
def finite(x,msg):
    if not math.isfinite(x):raise ValueError(f'nonfinite {msg}')
def rel_l2(x,y):
    need(len(x)==len(y),'numeric vector lengths differ');e=n=0.0
    for a,b in zip(x,y):
        a=float(a);b=float(b);finite(a,'candidate channel');finite(b,'reference channel');e+=(a-b)**2;n+=b*b
    return math.sqrt(e/n) if n else (0.0 if e==0 else math.inf)
def compare_csv(a,b,exact,values):
    x,y=rows(a),rows(b);need(len(x)==len(y),f'row count differs: {a} {b}')
    metrics={}
    for k in exact:need([r[k] for r in x]==[r[k] for r in y],f'exact grid/key mismatch {k}')
    for k in values:metrics[k]=rel_l2([r[k] for r in x],[r[k] for r in y])
    return metrics
def compare_bin(a,b):
    x=pathlib.Path(a).read_bytes();y=pathlib.Path(b).read_bytes()
    need(len(x)==len(y) and len(x)%8==0,f'binary output shape mismatch {a} {b}')
    xv=struct.unpack('<'+str(len(x)//8)+'d',x);yv=struct.unpack('<'+str(len(y)//8)+'d',y)
    return rel_l2(xv,yv)
def compare_science(a,b):
    out={}
    out['components']=compare_csv(a/'components.csv',b/'components.csv',['bin'],['up_re','up_im','north_re','north_im','east_re','east_im'])
    out['spectrum']=compare_csv(a/'spectrum.csv',b/'spectrum.csv',['bin','frequency_hz'],['real_m','imag_m'])
    out['trace']=compare_csv(a/'trace.csv',b/'trace.csv',['time_s','window'],['velocity_m_s','rate_m_s2','adjoint_source'])
    out['kernels']={name:compare_bin(a/name,b/name) for name in ('along.bin','across.bin')}
    for fn in ('nodes.csv','frequency-grid.csv','along-angles.csv','across-angles.csv'):
        x,y=rows(a/fn),rows(b/fn);need(x==y,f'exact selected grid geometry differs in {fn}')
    return out
def gate(metrics,limit,label):
    for group,channels in metrics.items():
        for channel,value in channels.items():
            finite(value,f'{label}.{group}.{channel}')
            need(value<=limit,f'{label}.{group}.{channel} relative L2 {value:.17g} exceeds {limit}')
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--scalar',type=pathlib.Path,required=True);ap.add_argument('--batched',type=pathlib.Path,required=True)
    ap.add_argument('--profile',type=pathlib.Path,required=True);ap.add_argument('--yspec',type=pathlib.Path,required=True)
    ap.add_argument('--accuracy',type=pathlib.Path,required=True);ap.add_argument('--out',type=pathlib.Path,required=True);a=ap.parse_args()
    scalar,batch,profile,yspec,accepted=[obj(p/'result.json') for p in (a.scalar,a.batched,a.profile,a.yspec,a.accuracy)]
    selected_h=0.02035382429229197
    for label,d in (('scalar',scalar),('batched',batch),('profile',profile)):
        need(d.get('backend')=='dspec' and d.get('degree')==200 and d.get('nq')==6 and d.get('frequency_count')==35 and d.get('first_bin')==7 and d.get('last_bin')==41,f'{label} DSpec workload differs')
        need(abs(float(d.get('refinement',-1))-selected_h)<1e-14 and abs(float(d.get('phase_d_mesh_scale',-1))-.9)<1e-14,f'{label} selected mesh differs')
        need(abs(float(d.get('relative_error',-1))-1e-5)<1e-15 and abs(float(d.get('phase_d_requested_relative_accuracy',-1))-1e-3)<1e-15,f'{label} accuracy/input metadata differs')
        need(d.get('frequency_workers')==8 and d.get('solve_workers')==1 and d.get('phase_c_direct_adjoint') is True and d.get('phase_c_combined_reconstruction') is True,f'{label} route/worker metadata differs')
        need(d.get('phase_c_validate_reconstruction') is False and d.get('phase_d_common_validation') is False and d.get('phase_d_validate_batch') in (False,None),f'{label} timing path contains validation work')
        need(d.get('phase_c_zero_bin_skips') is True,f'{label} exact zero skip setting is not enabled')
        need(int(d.get('gravity_max_abs_nd_error',math.inf))>=0 and float(d['gravity_max_abs_nd_error'])<=1e-10,f'{label} gravity control failed')
        need(abs(float(d['time_shift_adjoint_identity'])-1)<=1e-10,f'{label} time-shift identity failed')
    need(scalar.get('phase_d_scalar_only') is True and scalar.get('phase_d_batched_angular') is False,'scalar native did not run scalar-only path')
    need(batch.get('phase_d_scalar_only') is False and batch.get('phase_d_batched_angular') is True,'batched native route flag is wrong')
    need(profile.get('phase_d_scalar_only') is False and profile.get('phase_d_batched_angular') is True,'batched profile route flag is wrong')
    need(accepted.get('status')=='pass' and accepted.get('requested_accuracy_met') is True,'accuracy-passing selected validation receipt absent')
    need(abs(float(accepted['candidate']['h'])-selected_h)<1e-14,'selected validation h differs')
    ymesh=yspec.get('backend')=='yspec' and yspec.get('degree')==200 and yspec.get('nq')==6 and abs(float(yspec.get('refinement',-1))-selected_h)<1e-14
    need(ymesh and yspec.get('phase_d_batched_angular') is False,'YSpec refresh did not use selected mesh/route')
    need(yspec.get('workers')==8 and yspec.get('frequency_workers')==1 and yspec.get('solve_workers')==8,'YSpec workers differ from the selected 8-worker comparator')
    need(yspec.get('frequency_count')==35 and yspec.get('first_bin')==7 and yspec.get('last_bin')==41,'YSpec frequency grid differs')
    scalar_batch=compare_science(a.scalar,a.batched);profile_batch=compare_science(a.profile,a.batched);selected_batch=compare_science(a.batched,a.accuracy)
    gate(scalar_batch,1e-11,'scalar-v-batched');gate(profile_batch,1e-11,'profile-v-native');gate(selected_batch,1e-11,'timed-v-validated')
    yspec_metrics=compare_science(a.yspec,a.batched)
    zero=rows(a.batched/'phase-c2-zero-bin-support.csv')
    need([int(x['bin']) for x in zero]==[7,8,41],'exact-zero support bins differ')
    need(all(int(x['taper_exact_zero'])==1 and int(x['measured_P_exact_zero'])==1 and int(x['unit_adjoint_skipped'])==1 for x in zero),'batched native exact-zero skip evidence differs')
    prof=obj(a.profile/'profile.json');rows_prof=prof.get('rows',[]);counters=prof.get('counters',[])
    angular=[x for x in rows_prof if x.get('scope_name','').endswith('phase_d_tiled_eigen_angular_reconstruction') and x.get('phase')=='initial']
    angular_summary={'calls':sum(int(x.get('calls',0)) for x in angular),'inclusive_worker_elapsed_s':sum(float(x.get('inclusive_s',0)) for x in angular),'exclusive_worker_elapsed_s':sum(float(x.get('exclusive_s',0)) for x in angular)}
    cvals={}
    for x in counters:cvals[x['name']]=cvals.get(x['name'],0)+int(x['value'])
    need(cvals.get('phase_c2_forward_spectrum_frequency_calls',0)==35,'profile must retain35 forward frequency calls')
    need(cvals.get('phase_c2_unit_adjoint_frequency_calls',0)==32,'profile must retain32 nonzero unit-adjoint frequency calls')
    need(cvals.get('phase_c2_exact_zero_unit_adjoint_skips',0)==3,'profile must skip exactly3 static-zero unit-adjoint bins')
    need(angular_summary['calls']==32,'batched angular profile should run32 nonzero bins')
    scalar_time=float(scalar['elapsed_s']);batch_time=float(batch['elapsed_s']);yspec_time=float(yspec['elapsed_s'])
    need(all(math.isfinite(x) and x>0 for x in (scalar_time,batch_time,yspec_time)),'invalid driver elapsed timing')
    report={'status':'pass','selected_mesh':{'degree':200,'nq':6,'h':selected_h,'mesh_scale':.9,'legacy_relative_error':1e-5,'requested_accuracy':1e-3,'frequency_bins':35,'native_sides':batch['native_side_count'],'along_angles':batch['along_angles'],'across_angles':batch['across_angles']},
      'native_timing_pair':{'same_binary_and_inputs':True,'scalar_driver_elapsed_s':scalar_time,'batched_driver_elapsed_s':batch_time,'scalar_over_batched_speedup':scalar_time/batch_time,'scalar_result_elapsed_s':scalar.get('elapsed_s'),'batched_result_elapsed_s':batch.get('elapsed_s')},
      'fresh_matched_yspec':{'backend':'yspec','driver_elapsed_s':yspec_time,'matched_case_mesh_bins_and_requested_output_grid':True,'yspec_over_batched_dspec_time_ratio':yspec_time/batch_time,'comparison_metrics_diagnostic_only':yspec_metrics,'worker_metadata':{'yspec_cli_workers':yspec['workers'],'yspec_frequency_workers':yspec['frequency_workers'],'yspec_solve_workers':yspec['solve_workers'],'dspec_frequency_workers':batch['frequency_workers'],'dspec_solve_workers':batch['solve_workers']}},
      'equivalence_gates':{'scalar_vs_batched_relative_l2':scalar_batch,'profile_vs_native_relative_l2':profile_batch,'timed_batch_vs_accepted_accuracy_run_relative_l2':selected_batch},
      'exact_zero_support':[7,8,41],'phase_d_batched_profile':angular_summary,'profile_counters':{k:cvals.get(k,0) for k in ('phase_c2_forward_spectrum_frequency_calls','phase_c2_unit_adjoint_frequency_calls','phase_c2_exact_zero_unit_adjoint_skips')},
      'derived_angular_temp_bound_bytes':batch['phase_d_actual_temporary_peak_bytes'],'actual_process_memory_and_wall':'recorded in terminal supervisor state.json; worker elapsed sums are not master wall time',
      'limits':['single matched native pair; no repeated timing distribution','YSpec output differences are diagnostics and are not a new accuracy/discrepancy gate','profile angular time is an overlapping worker-elapsed sum, not wall time']}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps(report,sort_keys=True))
if __name__=='__main__':main()
