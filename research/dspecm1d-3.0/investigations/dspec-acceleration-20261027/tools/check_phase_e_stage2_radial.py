#!/usr/bin/env python3
"""Audit the first Phase E Stage 2 radial substage receipts and parity."""
from __future__ import annotations
import argparse,csv,json,math,statistics,struct
from pathlib import Path


def need(ok,msg):
    if not ok: raise ValueError(msg)

def obj(path): return json.loads(Path(path).read_text())

def rows(path):
    with Path(path).open(newline='') as f: return list(csv.DictReader(f))

def rel(values_a,values_b):
    need(len(values_a)==len(values_b),'numeric vector lengths differ')
    err=norm=0.0
    for x,y in zip(values_a,values_b):
        x=float(x);y=float(y);need(math.isfinite(x) and math.isfinite(y),'nonfinite compared value')
        err+=(x-y)**2;norm+=y*y
    return math.sqrt(err/norm) if norm else (0.0 if err==0 else math.inf)

def compare_csv(a,b,grid,channels):
    x,y=rows(a),rows(b);need(len(x)==len(y),f'row count differs: {a} {b}')
    for i,(rx,ry) in enumerate(zip(x,y)):
        for key in grid: need(rx[key]==ry[key],f'grid/key differs at {a}:{i} {key}')
    return {key:rel([r[key] for r in x],[r[key] for r in y]) for key in channels}

def compare_bin(a,b):
    x=Path(a).read_bytes();y=Path(b).read_bytes();need(len(x)==len(y) and len(x)%8==0,f'binary shape differs: {a} {b}')
    return rel(struct.unpack('<'+str(len(x)//8)+'d',x),struct.unpack('<'+str(len(y)//8)+'d',y))

def compare_science(a,b):
    metrics={
      'trace':compare_csv(a/'trace.csv',b/'trace.csv',['time_s','window'],['velocity_m_s','rate_m_s2','adjoint_source']),
      'spectrum':compare_csv(a/'spectrum.csv',b/'spectrum.csv',['bin','frequency_hz'],['real_m','imag_m']),
      'components':compare_csv(a/'components.csv',b/'components.csv',['bin'],['up_re','up_im','north_re','north_im','east_re','east_im']),
    }
    for fn in ('nodes.csv','frequency-grid.csv','along-angles.csv','across-angles.csv','mesh-elements.csv'):
        pa,pb=a/fn,b/fn
        need(pa.exists()==pb.exists(),f'section/geometry file presence differs in {fn}')
        if pa.exists(): need(rows(pa)==rows(pb),f'exact geometry/frequency grid differs in {fn}')
    for fn in ('along.bin','across.bin'):
        pa,pb=a/fn,b/fn
        need(pa.exists()==pb.exists(),f'kernel section presence differs in {fn}')
        if pa.exists(): metrics['along_kernel' if fn=='along.bin' else 'across_kernel']=compare_bin(pa,pb)
    return metrics

def gate(metrics,limit,label):
    for group,values in metrics.items():
        channels=values.items() if isinstance(values,dict) else [(group,values)]
        for key,value in channels:
            need(math.isfinite(float(value)) and 0<=float(value)<=limit,
                 f'{label}.{key} relative L2 {value:.17g} exceeds {limit}')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--run-root',type=Path,required=True)
    ap.add_argument('--run-plan',type=Path,required=True)
    ap.add_argument('--supervisor-state',type=Path,required=True)
    ap.add_argument('--frozen-d',type=Path,required=True)
    ap.add_argument('--frozen-y',type=Path,required=True)
    ap.add_argument('--frozen-l3',type=Path,required=True)
    ap.add_argument('--l3-base',type=Path,required=True)
    ap.add_argument('--l3-plus',type=Path,required=True)
    ap.add_argument('--l3-minus',type=Path,required=True)
    ap.add_argument('--l3-fd-check',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args(); run=a.run_root
    def data(name): return run/name/'data'
    plan=obj(a.run_plan);need(plan.get('thread_cap')==20,'run supervisor thread cap differs')
    configured={j['id']:j for j in plan.get('jobs',[])}
    for jid in [f'dspec-generic-{i}' for i in range(1,4)]+[f'dspec-batched-{i}' for i in range(1,4)]+['dspec-generic-profile','dspec-batched-profile']:
        command=configured[jid]['command'];joined=' '.join(command)
        route='generic' if 'generic' in jid else 'batched'
        need(f'DSFAST_RADIAL_ROUTE={route}' in joined and 'DSFAST_RADIAL_ORACLE=1' not in joined,f'{jid} route/oracle controls differ')
        need('DSPEC_PHASE_D_VALIDATE_BATCHED=0' in joined and 'DSPEC_PHASE_C_VALIDATE_RECONSTRUCTION=0' in joined and 'DSPEC_PHASE_C_VALIDATE_DIRECT_ADJOINT=0' in joined,f'{jid} contains a timing oracle')
    audit_command=' '.join(configured['dspec-radial-oracle-audit']['command'])
    need('DSFAST_RADIAL_ROUTE=batched' in audit_command and 'DSFAST_RADIAL_ORACLE=1' in audit_command and 'DSPEC_PHASE_D_VALIDATE_BATCHED=1' in audit_command and 'DSPEC_PHASE_C_SKIP_ZERO_KERNEL_BINS=0' in audit_command,'all-bin audit command does not enable required validators')
    need('DSPEC_FREQUENCY_WORKERS=8' in audit_command and 'DSPEC_PHASE_C_DIRECT_ADJOINT=1' in audit_command and 'DSPEC_PHASE_D_VALIDATE_COMMON=0' in audit_command,'all-bin audit worker/scientific controls differ')
    l3_command=' '.join(configured['dspec-l3-fd-audit']['command'])
    need('DSFAST_RADIAL_ROUTE=batched' in l3_command and 'DSFAST_RADIAL_ORACLE=1' in l3_command and 'DSPEC_PHASE_D_L3_FD_CONTROL=1' in l3_command,'L3 command does not exercise the new route and frozen FD control')
    for i in range(1,4):
        ycmd=configured[f'yspec-native-{i}']['command']
        need(ycmd[:3]==['/usr/bin/env','-u','DSPEC_FREQUENCY_WORKERS'] and 'DSPEC_PHASE_E_YSPEC_PROFILE=0' in ycmd,'YSpec refresh worker/profile controls differ')
    report={"status":"pass","stage":"Phase E Stage 2 radial first substage"}
    baseline=[data(f'dspec-generic-{i}') for i in range(1,4)]
    optimized=[data(f'dspec-batched-{i}') for i in range(1,4)]
    base_profiles=data('dspec-generic-profile');fast_profiles=data('dspec-batched-profile')
    audit=data('dspec-radial-oracle-audit')
    yspec=[data(f'yspec-native-{i}') for i in range(1,4)]
    expected_h=0.02035382429229197
    for label,paths,route in [('generic',baseline+ [base_profiles],'generic-map'),('batched',optimized+[fast_profiles],'batched-direct-index')]:
        for p in paths:
            d=obj(p/'result.json')
            need(d.get('backend')=='dspec' and d.get('degree')==200 and d.get('nq')==6 and d.get('frequency_count')==35 and d.get('first_bin')==7 and d.get('last_bin')==41,f'{label} workload metadata differs')
            need(abs(float(d.get('refinement',-1))-expected_h)<1e-14 and d.get('native_radial_route')==route,f'{label} route or selected mesh differs')
            need(d.get('native_derivative_weights_mode')=='native-dref-over-jacobian',f'{label} derivative convention differs')
            need(d.get('native_radial_modal_scratch_bound_bytes',math.inf)+d.get('native_radial_derivative_weight_cache_bytes',math.inf)<=32*1024*1024,f'{label} radial cache/scratch bound differs')
            need((d.get('native_side_count'),d.get('along_angles'),d.get('across_angles'))==(208,135,19),f'{label} native side/angle inventory differs')
            need(d.get('frequency_workers')==8 and d.get('solve_workers')==1 and d.get('phase_c_direct_adjoint') is True,f'{label} worker/adjoint route differs')
            need(d.get('phase_c_zero_bin_skips') is True and d.get('phase_c_validate_reconstruction') is False and d.get('phase_c_validate_direct_adjoint') is False and d.get('phase_d_common_validation') is False,f'{label} timed run includes a validation oracle')
            need(abs(float(d.get('gravity_max_abs_nd_error',math.inf)))<=1e-10 and abs(float(d.get('time_shift_adjoint_identity',0))-1)<=1e-10,f'{label} gravity/time-shift control failed')
    audit_result=obj(audit/'result.json')
    need(audit_result.get('backend')=='dspec' and audit_result.get('degree')==200 and audit_result.get('nq')==6 and audit_result.get('frequency_count')==35 and audit_result.get('first_bin')==7 and audit_result.get('last_bin')==41,'radial audit workload/grid differs')
    need(abs(float(audit_result.get('refinement',-1))-expected_h)<1e-14 and (audit_result.get('native_side_count'),audit_result.get('along_angles'),audit_result.get('across_angles'))==(208,135,19),'radial audit selected mesh/inventory differs')
    need(audit_result.get('frequency_workers')==8 and audit_result.get('solve_workers')==1,'radial audit worker configuration differs')
    need(audit_result.get('phase_c_direct_adjoint') is True,'radial audit direct-adjoint route differs')
    need(abs(float(audit_result.get('gravity_max_abs_nd_error',math.inf)))<=1e-10 and abs(float(audit_result.get('time_shift_adjoint_identity',0))-1)<=1e-10,'radial audit gravity/time-shift control failed')
    need(math.isfinite(float(audit_result.get('elapsed_s',math.nan))) and float(audit_result['elapsed_s'])>0,'radial audit driver wall is invalid')
    need(audit_result.get('native_radial_route')=='batched-direct-index' and audit_result.get('native_derivative_weights_mode')=='native-dref-over-jacobian','radial audit did not use the reviewed fast path/default derivative')
    need(audit_result.get('native_radial_modal_scratch_bound_bytes',math.inf)+audit_result.get('native_radial_derivative_weight_cache_bytes',math.inf)<=32*1024*1024,'radial audit cache/scratch bound differs')
    need(audit_result.get('phase_c_zero_bin_skips') is False,'all-bin audit skipped an exact-zero bin')
    need((audit/'phase-d-batched-scalar-equivalence.csv').is_file(),'all-bin field/product oracle output is absent')
    need(obj(audit/'result.json').get('phase_d_batched_angular') is True,'selected NQ6 angular route changed')
    need(float(audit_result.get('phase_d_actual_temporary_peak_bytes',math.inf))<=32*1024*1024,'angular scratch bound changed')
    parity=compare_science(audit,a.frozen_d);gate(parity,1e-11,'NQ6-audit-v-frozen-PhaseD')
    for i,(g,f) in enumerate(zip(baseline,optimized),1):
        metrics=compare_science(g,f);gate(metrics,1e-11,f'generic-v-batched-repeat-{i}')
    for p in baseline+optimized+[base_profiles,fast_profiles]:
        metrics=compare_science(p,a.frozen_d);gate(metrics,1e-11,f'route-v-frozen-{p.name}')
    proof=[]
    proof_metrics=['all_radial_rel_l2','value_slots_rel_l2','derivative_slots_rel_l2','forward_rel_l2','adjoint_rel_l2','fluid_rel_l2','solid_rel_l2','strain_rel_l2']
    proof_rows=[]
    for bin_id in range(7,42):
        p=audit/f'radial-oracle.bin{bin_id}.csv'
        with p.open(newline='') as f: parsed=list(csv.reader(f))
        need(len(parsed)==1 and len(parsed[0])==16,f'radial oracle proof shape differs for bin {bin_id}')
        fields=parsed[0];r={'bin':fields[0],'compared_modes':fields[1],'conservative_capture_bytes':fields[2],
            'forward_modes':fields[3],'adjoint_modes':fields[4],'fluid_modes':fields[5],'solid_modes':fields[6],
            **dict(zip(proof_metrics,fields[7:15])),'route_marker':fields[15]}
        proof.append(r);need(int(r['bin'])==bin_id and r['route_marker']=='generic-and-batched-keyed-modes',f'radial proof bin/route key differs at {bin_id}')
    reference_inventory=None
    for r in proof:
        total=int(r['compared_modes']);need(total>0,f"empty radial proof bin {r['bin']}")
        need(int(r['forward_modes'])+int(r['adjoint_modes'])==total,f"forward/adjoint count mismatch bin {r['bin']}")
        need(int(r['fluid_modes'])+int(r['solid_modes'])==total,f"fluid/solid count mismatch bin {r['bin']}")
        need(all(int(r[k])>0 for k in ('forward_modes','adjoint_modes','fluid_modes','solid_modes')),f"radial category absent bin {r['bin']}")
        inventory=tuple(int(r[k]) for k in ('compared_modes','forward_modes','adjoint_modes','fluid_modes','solid_modes'))
        if reference_inventory is None: reference_inventory=inventory
        need(inventory==reference_inventory,f"all-bin radial inventory differs in bin {r['bin']}")
        need(int(r['conservative_capture_bytes'])<=256*1024*1024,f"oracle capture bound exceeded bin {r['bin']}")
        for key in proof_metrics:
            v=float(r[key]);need(math.isfinite(v) and v<=1e-11,f"{key} exceeds 1e-11 in bin {r['bin']}: {v}")
    scalar=rows(audit/'phase-d-batched-scalar-equivalence.csv')
    need(len(scalar)==70 and {(int(r['bin']),int(r['section'])) for r in scalar}=={(b,s) for b in range(7,42) for s in (0,1)},'same-workspace full-product proof does not cover all 35 bins and sections')
    scalar_max={k:0.0 for k in ('forward_field_rel_l2','adjoint_field_rel_l2','complex_product_rel_l2')}
    for r in scalar:
        for k in scalar_max:
            v=float(r[k]);need(math.isfinite(v) and 0<=v<=1e-11,f'same-workspace {k} gate failed at bin {r["bin"]}');scalar_max[k]=max(scalar_max[k],v)
    acc=plan.get('accuracy_anchor',{})
    need(acc.get('status')=='pass' and acc.get('requested_accuracy_met') is True and acc.get('accuracy_report_sha256'),'prior NQ10-bounded 0.1% accuracy evidence is not hash-bound')
    need(float(acc['max_solid_complex_field_relative_l2'])<1e-3 and float(acc['max_solid_full35_kernel_relative_l2'])<1e-3,'prior NQ10 support for 0.1% gate is absent')
    native_times={};supervisor=obj(a.supervisor_state);jobrows={x['id']:x for x in supervisor.get('jobs',[])}
    need(supervisor.get('status') in ('running','completed','failed'),'supervisor run state is invalid')
    expected_jobs={'dspec-radial-oracle-audit'}|{f'dspec-{route}-{i}' for route in ('generic','batched') for i in range(1,4)}|{'dspec-generic-profile','dspec-batched-profile'}|{f'yspec-native-{i}' for i in range(1,4)}|{'dspec-l3-fd-audit','phase-e-stage2-radial-fd-check','phase-e-stage2-radial-terminal-check'}
    need(expected_jobs.issubset(jobrows),f'supervisor state missing jobs: {sorted(expected_jobs-set(jobrows))}')
    prior_terminal=jobrows['phase-e-stage2-radial-terminal-check']
    prior_terminal_recovery=None
    if supervisor.get('status')=='failed':
        before=prior_terminal.get('inputs_sha256_before',{});after=prior_terminal.get('inputs_sha256_after',{})
        changed=[p for p,h in before.items() if after.get(p)!=h]
        expected_failure={'missing_outputs:data/phase-e-stage2-radial-report.json','input_hash_changed_during_job'}
        need(prior_terminal.get('status')=='failed' and prior_terminal.get('returncode')==1 and
             set(prior_terminal.get('failure_reasons',[]))==expected_failure,
             'prior terminal-check failure is not the preserved known output/state self-reference failure')
        need(len(changed)==1 and changed[0].endswith('/state.json') and
             before.get(changed[0])=='97d2c2420e8be075918f4e42cd8c497d20c2a233ba1d4108ff529a00acd5b59b' and
             after.get(changed[0])=='21ab218ea18fea149d0da940b130d45b166ffae0f16a6ca430051e325f648465',
             'prior terminal-check changed inputs differ from the recorded supervisor-state self-update')
        prior_terminal_recovery={'prior_task_id':supervisor.get('task_id'),
             'prior_terminal_status':'failed','prior_terminal_returncode':1,
             'prior_failure_reasons':prior_terminal['failure_reasons'],
             'only_changed_input':changed[0],
             'changed_input_sha256_before':before[changed[0]],
             'changed_input_sha256_after':after[changed[0]],
             'recovery_reason':'the terminal checker read the mutable state file of its own supervisor task; this retry reads the finalized immutable state snapshot, validates the 14 prior scientific/FD jobs, skips only the known failed terminal-check lifecycle, and relies on the retry supervisor receipt for the new terminal-check lifecycle'}
    else:
        need(prior_terminal.get('status') in ('running','completed'),
             'terminal-check lifecycle status in supervisor snapshot is invalid')
    for job_id,rec in jobrows.items():
        if job_id=='phase-e-stage2-radial-terminal-check': continue
        need(rec.get('returncode')==0 and rec.get('status')=='completed' and rec.get('inputs_unchanged') is True and not rec.get('failure_reasons') and rec.get('cleanup',{}).get('verified_gone') is True,f'{job_id} supervisor receipt is not clean')
        need(int(rec.get('peak_observed_threads',0))<=20 and float(rec.get('peak_observed_rss_gib_userspace',math.inf))<=16,f'{job_id} resource cap exceeded')
        need(math.isfinite(float(rec.get('elapsed_seconds',math.nan))) and float(rec['elapsed_seconds'])>0,f'{job_id} supervisor wall is invalid')
    def profile_check(path,route):
        profile=obj(path/'profile.json')
        phase='combined_forward_unit_adjoint_reconstruction_contraction'
        profile_rows=[r for r in profile.get('rows',[]) if r.get('phase')==phase]
        for row in profile_rows:
            inc=float(row.get('inclusive_s',math.nan));exc=float(row.get('exclusive_s',math.nan))
            need(math.isfinite(inc) and math.isfinite(exc) and 0<=exc<=inc+1e-9,
                 f'{route} raw profile scope invalid: {row.get("scope_name")}')
        by={r.get('scope_name','').removeprefix('frequency_worker.'):r for r in profile_rows}
        radial='phase_e_radial_mode_interpolation_and_strain_table_build'
        needed=['combined_frequency_workspace_forward_unit_adjoint_and_reconstruction',
                'frequency_workspace_prepare_and_modal_capture',
                'phase_d_stream_projected_radial_modes_to_strain_tables',
                'phase_e_angular_plan_and_weight_preparation',
                'phase_d_tiled_eigen_angular_reconstruction',radial,
                'phase_e_angular_tile_weight_packing',
                'phase_e_angular_coefficient_row_lookup_and_packing',
                'phase_e_angular_eigen_gemm',
                'phase_e_angular_tile_scatter',
                'phase_e_angular_material_evaluation_and_contraction']
        if route=='generic-map': needed.append('phase_e_angular_mode_lookup_materialization')
        for name in needed:
            row=by.get(name);need(row is not None and int(row.get('calls',0))>0 and float(row.get('inclusive_s',0))>0,f'{route} profile scope missing: {name}')
            inc=float(row['inclusive_s']);exc=float(row['exclusive_s']);need(math.isfinite(inc) and math.isfinite(exc) and 0<=exc<=inc+1e-9,f'{route} profile scope invalid: {name}')
        root=by['combined_frequency_workspace_forward_unit_adjoint_and_reconstruction'];need(int(root['calls'])==35,f'{route} profile forward-bin count differs')
        need(int(by[radial]['calls'])==32,f'{route} radial reconstruction call count differs')
        for name in needed[5:]:
            row=by[name]
            need(int(row['calls'])>0 and float(row['inclusive_s'])>0,
                 f'{route} fine angular profile scope is empty: {name}')
        counters={(r.get('phase'),r.get('name')):int(r.get('value',0)) for r in profile.get('counters',[])}
        need(counters.get((phase,'phase_c2_forward_spectrum_frequency_calls'))==35,f'{route} profile forward coverage differs')
        need(counters.get((phase,'phase_c2_unit_adjoint_frequency_calls'))==32,f'{route} profile adjoint coverage differs')
        need(counters.get((phase,'phase_c2_exact_zero_unit_adjoint_skips'))==3,f'{route} profile exact-zero skips differ')
        need(counters.get((phase,'phase_e_angular_eigen_gemm_calls'))==147840,
             f'{route} angular GEMM counter differs')
        need(int(by['phase_e_angular_eigen_gemm']['calls'])==147840,
             f'{route} angular GEMM scope call count differs')
        for counter in ('phase_e_angular_weight_tiles',
                        'phase_e_angular_coefficient_row_batches',
                        'phase_e_angular_material_contraction_batches'):
            need(counters.get((phase,counter),0)>0,f'{route} angular profile counter is empty: {counter}')
        wall=profile.get('wall_partition',{})
        for key in ('profile_pipeline_wall_s','master_exclusive_partition_s','residual_s','tolerance_s'):
            need(math.isfinite(float(wall.get(key,math.nan))),f'{route} profile wall partition missing {key}')
        need(abs(float(wall['residual_s']))<=float(wall['tolerance_s']),f'{route} profile master-wall partition fails')
        worker_exclusive=sum(float(r['exclusive_s']) for r in profile_rows if r.get('execution_domain')=='worker_elapsed_sum')
        root_inclusive=float(root['inclusive_s'])
        need(abs(worker_exclusive-root_inclusive)<=max(1e-6,1e-3*max(root_inclusive,1e-9)),f'{route} worker-sum closure fails')
        return {'radial_scope':by[radial],'profile_pipeline_wall_s':float(wall['profile_pipeline_wall_s']),
                'master_exclusive_partition_s':float(wall['master_exclusive_partition_s']),
                'worker_exclusive_sum_s':worker_exclusive,'combined_root_inclusive_worker_s':root_inclusive,
                'scope_names':sorted(by)}
    profile_summaries={'generic':profile_check(base_profiles,'generic-map'),'batched':profile_check(fast_profiles,'batched-direct-index')}
    for p,jid in ((base_profiles,'dspec-generic-profile'),(fast_profiles,'dspec-batched-profile')):
        value=float(obj(p/'result.json').get('elapsed_s',math.nan));need(math.isfinite(value) and value>0,f'{jid} driver elapsed is invalid')
        need(math.isfinite(float(jobrows[jid]['elapsed_seconds'])) and float(jobrows[jid]['elapsed_seconds'])>0,f'{jid} supervisor wall is invalid')
    def times(paths,prefix):
        vals=[float(obj(p/'result.json')['elapsed_s']) for p in paths]
        ids=[f'{prefix}-{i}' for i in range(1,4)]
        walls=[float(jobrows[i]['elapsed_seconds']) for i in ids]
        rss=[float(jobrows[i]['peak_observed_rss_gib_userspace']) for i in ids]
        threads=[int(jobrows[i]['peak_observed_threads']) for i in ids]
        for value in vals+walls: need(math.isfinite(value) and value>0,f'{prefix} driver/supervisor time is invalid')
        for job_id in ids:
            rec=jobrows[job_id]
            need(rec.get('returncode')==0 and rec.get('status')=='completed' and rec.get('inputs_unchanged') is True and not rec.get('failure_reasons') and rec.get('cleanup',{}).get('verified_gone') is True,f'{job_id} supervisor receipt is not clean')
        need(max(rss)<=16 and max(threads)<=20,f'{prefix} observed resources exceed configured numerical caps')
        return {'driver_elapsed_s':vals,'driver_median_s':statistics.median(vals),'driver_min_s':min(vals),'driver_max_s':max(vals),'supervisor_wall_s':walls,'supervisor_median_s':statistics.median(walls),'peak_observed_rss_gib_userspace':max(rss),'peak_observed_threads_including_supervisor':max(threads)}
    native_times['generic']=times(baseline,'dspec-generic');native_times['batched']=times(optimized,'dspec-batched')
    yinfo=[]
    for i,p in enumerate(yspec,1):
        d=obj(p/'result.json');need(d.get('backend')=='yspec' and d.get('degree')==200 and d.get('nq')==6 and abs(float(d.get('refinement',-1))-expected_h)<1e-14,'YSpec refresh workload differs')
        need(d.get('workers')==8 and d.get('frequency_workers')==1 and d.get('solve_workers')==8 and d.get('frequency_count')==35 and d.get('first_bin')==7 and d.get('last_bin')==41,'YSpec refresh worker/grid configuration differs')
        metrics=compare_science(p,a.frozen_y);gate(metrics,1e-11,'YSpec-v-Stage1-matched-reference')
        yinfo.append(float(d['elapsed_s']))
        j=jobrows[f'yspec-native-{i}'];need(math.isfinite(yinfo[-1]) and yinfo[-1]>0 and math.isfinite(float(j['elapsed_seconds'])) and float(j['elapsed_seconds'])>0,'YSpec timing is invalid')
    yjobs=[jobrows[f'yspec-native-{i}'] for i in range(1,4)]
    l3=data('dspec-l3-fd-audit');l3result=obj(l3/'result.json')
    need(l3result.get('degree')==3 and l3result.get('nq')==8 and l3result.get('native_radial_route')=='batched-direct-index','L3 audit did not execute the new generic-NQ native route')
    need(l3result.get('native_side_count')==157 and l3result.get('along_angles')==91 and l3result.get('frequency_count')==35,'L3 audit mesh/grid inventory differs')
    need(l3result.get('native_derivative_weights_mode')=='native-dref-over-jacobian' and l3result.get('phase_c_zero_bin_skips') is False,'L3 route/derivative/all-bin control differs')
    l3proof=[]
    for b in range(7,42):
        with (l3/f'radial-oracle.bin{b}.csv').open(newline='') as f: parsed=list(csv.reader(f))
        need(len(parsed)==1 and len(parsed[0])==16 and int(parsed[0][0])==b,f'L3 radial proof missing bin {b}')
        r=dict(zip(['bin','compared_modes','conservative_capture_bytes','forward_modes','adjoint_modes','fluid_modes','solid_modes']+proof_metrics+['route_marker'],parsed[0]))
        total=int(r['compared_modes'])
        need(total>0 and int(r['forward_modes'])+int(r['adjoint_modes'])==total and int(r['fluid_modes'])+int(r['solid_modes'])==total,'L3 radial proof inventory differs')
        need(r['route_marker']=='generic-and-batched-keyed-modes',f'L3 radial route marker differs at bin {b}')
        need(all(int(r[k])>0 for k in ('forward_modes','adjoint_modes','fluid_modes','solid_modes')),'L3 radial category is absent')
        need(int(r['conservative_capture_bytes'])<=256*1024*1024,'L3 radial proof bound differs')
        for k in proof_metrics:
            v=float(r[k]);need(math.isfinite(v) and 0<=v<=1e-11,f'L3 {k} fails at bin {r["bin"]}')
        l3proof.append(r)
    fd=obj(a.l3_fd_check);fd_error=float(fd.get('fd_relative_error',math.nan));need(fd.get('status')=='pass' and math.isfinite(fd_error) and 0<=fd_error<=2e-4,'frozen L3 central-difference control failed')
    l3parity=compare_science(l3,a.frozen_l3);gate(l3parity,1e-11,'L3-v-frozen-PhaseD')
    report.update({
      'prior_terminal_checker_recovery':prior_terminal_recovery,
      'selected_case':{'degree':200,'nq':6,'h':expected_h,'frequency_bins':35,'timed_zero_skips':[7,8,41],'audit_zero_skips':[],'native_sides':208,'angles':[135,19]},
      'radial_oracle':{'bins':len(proof),'compared_modes_total':sum(int(r['compared_modes']) for r in proof),'maximum_conservative_capture_bytes':max(int(r['conservative_capture_bytes']) for r in proof),'max_relative_l2_by_category':{k:max(float(r[k]) for r in proof) for k in proof_metrics}},
      'same_workspace_full_pipeline_scalar_parity_max_relative_l2':scalar_max,
      'timing':{'generic':native_times['generic'],'batched':native_times['batched'],'yspec_driver_elapsed_s':yinfo,'yspec_median_s':statistics.median(yinfo),'yspec_min_s':min(yinfo),'yspec_max_s':max(yinfo),'yspec_supervisor_wall_s':[float(j['elapsed_seconds']) for j in yjobs],'profile_driver_elapsed_s':{'generic':float(obj(base_profiles/'result.json')['elapsed_s']),'batched':float(obj(fast_profiles/'result.json')['elapsed_s'])},'profiles':profile_summaries},
      'l3_fd':{'status':'pass','fd_relative_error':fd['fd_relative_error'],'limit':2e-4,'route':l3result['native_radial_route'],'bins_with_radial_oracle':len(l3proof)},
      'parity_against_frozen_phase_d':parity,
      'l3_parity_against_frozen_phase_d':l3parity,
      'nq10_accuracy_evidence':{'report_path':acc['accuracy_report_path'],'report_sha256':acc['accuracy_report_sha256'],'max_solid_complex_field_relative_l2':acc['max_solid_complex_field_relative_l2'],'max_solid_full35_kernel_relative_l2':acc['max_solid_full35_kernel_relative_l2'],'both_below_1e-3':True},
      'resource_summary':{k:{'elapsed_s':v['elapsed_seconds'],'peak_rss_gib_userspace':v['peak_observed_rss_gib_userspace'],'peak_threads_including_supervisor':v['peak_observed_threads']} for k,v in jobrows.items() if k.startswith('dspec-') or k.startswith('yspec-')},
      'limits':['worker sums are not wall time','user-space RSS/thread watchdog is not kernel aggregate accounting','matched YSpec is a context refresh, not equivalent arithmetic'],
    })
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    Path(a.out).write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps(report,sort_keys=True))
if __name__=='__main__': main()
