#!/usr/bin/env python3
"""Check Stage 2 angular-route parity, controls, profiles, and matched timings."""
from __future__ import annotations
import argparse,csv,hashlib,json,math,statistics,struct
from pathlib import Path

def need(ok,msg):
    if not ok: raise ValueError(msg)
def obj(path): return json.loads(Path(path).read_text())
def rows(path):
    with Path(path).open(newline='') as f:return list(csv.DictReader(f))
def rel(a,b):
    need(len(a)==len(b),'numeric vectors differ in length');err=norm=0.0
    for x,y in zip(a,b):
        x=float(x);y=float(y);need(math.isfinite(x) and math.isfinite(y),'nonfinite value in comparison');err+=(x-y)**2;norm+=y*y
    return math.sqrt(err/norm) if norm else (0.0 if err==0 else math.inf)
def compare_csv(a,b,keys,channels):
    x,y=rows(a),rows(b);need(len(x)==len(y),f'row count differs: {a} {b}')
    for i,(u,v) in enumerate(zip(x,y)):
        for k in keys:need(u[k]==v[k],f'grid/key differs {a}:{i}:{k}')
    return {k:rel([r[k] for r in x],[r[k] for r in y]) for k in channels}
def compare_bin(a,b):
    x=Path(a).read_bytes();y=Path(b).read_bytes();need(len(x)==len(y) and len(x)%8==0,f'binary shape differs: {a} {b}')
    return rel(struct.unpack('<'+str(len(x)//8)+'d',x),struct.unpack('<'+str(len(y)//8)+'d',y))
def compare_science(a,b):
    m={'trace':compare_csv(a/'trace.csv',b/'trace.csv',['time_s','window'],['velocity_m_s','rate_m_s2','adjoint_source']),
       'spectrum':compare_csv(a/'spectrum.csv',b/'spectrum.csv',['bin','frequency_hz'],['real_m','imag_m']),
       'components':compare_csv(a/'components.csv',b/'components.csv',['bin'],['up_re','up_im','north_re','north_im','east_re','east_im'])}
    for fn in ('nodes.csv','frequency-grid.csv','along-angles.csv','across-angles.csv','mesh-elements.csv'):
        x,y=a/fn,b/fn;need(x.exists()==y.exists(),f'geometry inventory differs for {fn}')
        if x.exists():need(rows(x)==rows(y),f'geometry/frequency grid differs for {fn}')
    for fn in ('along.bin','across.bin'):
        x,y=a/fn,b/fn;need(x.exists()==y.exists(),f'section inventory differs for {fn}')
        if x.exists():m[fn]=compare_bin(x,y)
    return m
def gate(metrics,limit,label):
    for group,values in metrics.items():
        pairs=values.items() if isinstance(values,dict) else [(group,values)]
        for name,v in pairs:need(math.isfinite(float(v)) and 0<=float(v)<=limit,f'{label}.{name}={v} exceeds {limit}')
def route_check(plan,ids,radial,angular,oracle='0'):
    jobs={j['id']:j for j in plan['jobs']}
    for jid in ids:
        cmd=' '.join(jobs[jid]['command'])
        for setting in (f'DSFAST_RADIAL_ROUTE={radial}',f'DSFAST_ANGULAR_ROUTE={angular}',f'DSFAST_ANGULAR_ORACLE={oracle}'):
            need(setting in cmd,f'{jid} missing {setting}')
        need('DSPEC_PHASE_D_VALIDATE_BATCHED=0' in cmd and 'DSPEC_PHASE_C_VALIDATE_RECONSTRUCTION=0' in cmd and 'DSPEC_PHASE_C_VALIDATE_DIRECT_ADJOINT=0' in cmd,f'{jid} includes a timed oracle')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-root',type=Path,required=True);ap.add_argument('--run-plan',type=Path,required=True);ap.add_argument('--supervisor-state',type=Path,required=True)
    for n in ('frozen-d','frozen-y','frozen-l3','l3-base','l3-plus','l3-minus','l3-fd-check','out'):ap.add_argument('--'+n,type=Path,required=True)
    a=ap.parse_args();run=a.run_root;plan=obj(a.run_plan);need(plan.get('thread_cap')==20,'supervisor cap differs')
    jobs={j['id']:j for j in plan['jobs']};timed_legacy=[f'dspec-legacy-angle-{i}' for i in range(1,4)];timed_real=[f'dspec-real-angle-{i}' for i in range(1,4)]
    route_check(plan,timed_legacy+['dspec-legacy-angle-profile'],'batched','legacy-complex')
    route_check(plan,timed_real+['dspec-real-angle-profile'],'batched','real-legendre')
    auditcmd=' '.join(jobs['dspec-angular-all-bin-audit']['command'])
    need(all(x in auditcmd for x in ('DSFAST_RADIAL_ROUTE=batched','DSFAST_ANGULAR_ROUTE=real-legendre','DSFAST_ANGULAR_ORACLE=1','DSFAST_RADIAL_ORACLE=1','DSPEC_PHASE_D_VALIDATE_BATCHED=1','DSPEC_PHASE_C_SKIP_ZERO_KERNEL_BINS=0')),'all-bin angular/radial audit controls differ')
    l3cmd=' '.join(jobs['dspec-l3-fd-audit']['command'])
    need(all(x in l3cmd for x in ('DSFAST_RADIAL_ROUTE=batched','DSFAST_ANGULAR_ROUTE=real-legendre','DSFAST_ANGULAR_ORACLE=1','DSPEC_PHASE_D_L3_FD_CONTROL=1')),'L3 job does not exercise actual angular route and FD control')
    for i in range(1,4):
        ycmd=jobs[f'yspec-native-{i}']['command']
        need('DSPEC_PHASE_E_YSPEC_PROFILE=0' in ' '.join(ycmd),'YSpec refresh has wrong profile state')
        need(ycmd[:3]==['/usr/bin/env','-u','DSPEC_FREQUENCY_WORKERS'],'YSpec refresh must clear inherited frequency worker count')
    state=obj(a.supervisor_state);need(state.get('status')=='completed','producer supervisor is not completed')
    need(state.get('plan_sha256')==hashlib.sha256(a.run_plan.read_bytes()).hexdigest(),'producer state is not bound to the supplied run plan')
    jr={j['id']:j for j in state.get('jobs',[])};need(set(jobs)<=set(jr),'producer supervisor is missing declared jobs')
    for jid in jobs:
        r=jr[jid];need(r.get('status')=='completed' and r.get('returncode')==0 and r.get('inputs_unchanged') is True and not r.get('failure_reasons') and r.get('cleanup',{}).get('verified_gone') is True,f'{jid} producer receipt is not clean')
        need(int(r.get('peak_observed_threads',99))<=20 and float(r.get('peak_observed_rss_gib_userspace',99))<=16,f'{jid} observed resource cap exceeded')
        need(math.isfinite(float(r.get('elapsed_seconds',math.nan))) and float(r['elapsed_seconds'])>0,f'{jid} supervisor wall invalid')
    def data(jid):return run/jid/'data'
    audit=data('dspec-angular-all-bin-audit');ar=obj(audit/'result.json');h=0.02035382429229197
    def selected(d,label,angular_route,skip,oracle):
        need(d.get('backend')=='dspec' and d.get('degree')==200 and d.get('nq')==6 and d.get('frequency_count')==35 and d.get('first_bin')==7 and d.get('last_bin')==41,f'{label} workload differs')
        need(abs(float(d.get('refinement',-1))-h)<1e-14 and (d.get('native_side_count'),d.get('along_angles'),d.get('across_angles'))==(208,135,19),f'{label} mesh/angle inventory differs')
        need(d.get('native_radial_route')=='batched-direct-index' and d.get('native_angular_route')==angular_route and d.get('phase_e_angular_oracle') is oracle,f'{label} route/oracle metadata differs')
        need(d.get('native_derivative_weights_mode')=='native-dref-over-jacobian' and d.get('frequency_workers')==8 and d.get('solve_workers')==1 and d.get('phase_c_direct_adjoint') is True,f'{label} derivative/worker configuration differs')
        need(d.get('phase_c_validate_direct_adjoint') is False and d.get('phase_c_validate_reconstruction') is False and d.get('phase_d_common_validation') is False,f'{label} unexpected validation oracle')
        need(d.get('phase_c_zero_bin_skips') is skip and abs(float(d.get('gravity_max_abs_nd_error',math.inf)))<=1e-10 and abs(float(d.get('time_shift_adjoint_identity',0))-1)<=1e-10,f'{label} skip/gravity/time-shift control differs')
        need(math.isfinite(float(d.get('phase_d_actual_temporary_peak_bytes',math.inf))) and float(d['phase_d_actual_temporary_peak_bytes'])<=32*1024*1024,f'{label} angular scratch bound exceeded')
    selected(ar,'all-bin audit','real-legendre',False,True)
    for d,label in ((ar,'all-bin audit'),):
        need(d.get('backend')=='dspec' and d.get('degree')==200 and d.get('nq')==6 and d.get('frequency_count')==35 and d.get('first_bin')==7 and d.get('last_bin')==41,f'{label} workload differs')
        need(abs(float(d.get('refinement',-1))-h)<1e-14 and d.get('native_side_count')==208 and (d.get('along_angles'),d.get('across_angles'))==(135,19),f'{label} mesh/angle inventory differs')
        need(d.get('native_radial_route')=='batched-direct-index' and d.get('native_angular_route')=='real-legendre' and d.get('phase_e_angular_oracle') is True,f'{label} did not use selected routes/oracle')
        need(d.get('native_derivative_weights_mode')=='native-dref-over-jacobian' and d.get('frequency_workers')==8 and d.get('solve_workers')==1,f'{label} derivative/worker controls differ')
        need(d.get('phase_c_direct_adjoint') is True and d.get('phase_c_validate_direct_adjoint') is False and d.get('phase_c_validate_reconstruction') is False and d.get('phase_d_common_validation') is False,f'{label} direct-adjoint/validation controls differ')
        need(d.get('phase_c_zero_bin_skips') is False and abs(float(d.get('gravity_max_abs_nd_error',math.inf)))<=1e-10 and abs(float(d.get('time_shift_adjoint_identity',0))-1)<=1e-10,f'{label} zero-bin/gravity/time-shift control failed')
        need(float(d.get('phase_d_actual_temporary_peak_bytes',math.inf))<=32*1024*1024,f'{label} angular scratch bound exceeded')
    scalar=rows(audit/'phase-d-batched-scalar-equivalence.csv');need(len(scalar)==70 and {(int(r['bin']),int(r['section'])) for r in scalar}=={(b,s) for b in range(7,42) for s in (0,1)},'all-bin scalar field/product parity inventory differs')
    for r in scalar:
        for k in ('forward_field_rel_l2','adjoint_field_rel_l2','complex_product_rel_l2'):
            v=float(r[k]);need(math.isfinite(v) and 0<=v<=1e-11,f'scalar {k} fails at bin {r["bin"]}')
        for k in ('forward_scalar_norm2','adjoint_scalar_norm2','product_scalar_norm2'):
            v=float(r[k]);need(math.isfinite(v) and v>0,f'scalar reference norm {k} is not finite positive at bin {r["bin"]}')
    angular=rows(audit/'phase-e-angular-route-equivalence.csv');need(len(angular)==70 and {(int(r['bin']),int(r['section'])) for r in angular}=={(b,s) for b in range(7,42) for s in (0,1)},'same-workspace angular route proof inventory differs')
    angular_max={k:0.0 for k in ('forward_field_rel_l2','adjoint_field_rel_l2','complex_product_rel_l2')}
    angular_norms={'forward_field_rel_l2':'forward_legacy_norm2','adjoint_field_rel_l2':'adjoint_legacy_norm2','complex_product_rel_l2':'product_legacy_norm2'}
    for r in angular:
        for k in angular_max:
            v=float(r[k]);n=float(r[angular_norms[k]]);need(math.isfinite(v) and 0<=v<=1e-11 and math.isfinite(n) and n>0,f'angular proof {k}/norm fails at bin {r["bin"]}');angular_max[k]=max(angular_max[k],v)
    # Radial interpolation remains the independently checked Stage 2 first substage.
    radial=[];radial_inventory=None
    for b in range(7,42):
        with (audit/f'radial-oracle.bin{b}.csv').open(newline='') as f:q=list(csv.reader(f))
        need(len(q)==1 and len(q[0])==16 and int(q[0][0])==b,f'radial oracle missing bin {b}')
        r=q[0];need(r[-1]=='generic-and-batched-keyed-modes',f'radial oracle route marker differs at bin {b}')
        counts=(int(r[1]),int(r[3]),int(r[4]),int(r[5]),int(r[6]))
        need(counts[0]>0 and counts[1]+counts[2]==counts[0] and counts[3]+counts[4]==counts[0],f'radial oracle inventory differs at bin {b}')
        need(all(x>0 for x in counts[1:]),f'radial oracle category absent at bin {b}')
        need(radial_inventory is None or radial_inventory==counts, f'radial inventory changed at bin {b}');radial_inventory=counts
        need(int(r[2])<=256*1024*1024,f'radial oracle capture bound exceeded at bin {b}')
        radial.append({'bin':b,'row':r})
        for v in r[7:15]:need(math.isfinite(float(v)) and 0<=float(v)<=1e-11,f'radial oracle fails at bin {b}')
    need(compare_science(audit,a.frozen_d),'') if False else None
    audit_parity=compare_science(audit,a.frozen_d);gate(audit_parity,1e-11,'all-bin angular-v-radial-frozen')
    legacy=[data(x) for x in timed_legacy];real=[data(x) for x in timed_real]
    for i,(x,y) in enumerate(zip(legacy,real),1):gate(compare_science(x,y),1e-11,f'legacy-v-real-repeat-{i}')
    for p in legacy+real+[data('dspec-legacy-angle-profile'),data('dspec-real-angle-profile')]:gate(compare_science(p,a.frozen_d),1e-11,f'{p.name}-v-frozen-D')
    def profile(path,route):
        p=obj(path/'profile.json');phase='combined_forward_unit_adjoint_reconstruction_contraction';raw=[x for x in p.get('rows',[]) if x.get('phase')==phase]
        for r in raw:
            inc=float(r.get('inclusive_s',math.nan));exc=float(r.get('exclusive_s',math.nan));need(math.isfinite(inc) and math.isfinite(exc) and 0<=exc<=inc+1e-9,f'{route} raw profile row invalid')
        by={r.get('scope_name','').removeprefix('frequency_worker.'):r for r in raw}
        common=['combined_frequency_workspace_forward_unit_adjoint_and_reconstruction','frequency_workspace_prepare_and_modal_capture','phase_d_stream_projected_radial_modes_to_strain_tables','phase_e_angular_plan_and_weight_preparation','phase_d_tiled_eigen_angular_reconstruction','phase_e_radial_mode_interpolation_and_strain_table_build','phase_e_angular_coefficient_row_lookup_and_packing','phase_e_angular_tile_scatter','phase_e_angular_material_evaluation_and_contraction']
        gemm='phase_e_angular_eigen_gemm' if route=='legacy-complex' else 'phase_e_angular_real_legendre_gemm'
        for n in common+[gemm]:need(n in by and int(by[n].get('calls',0))>0 and float(by[n].get('inclusive_s',0))>0,f'{route} required profile scope missing: {n}')
        need(int(by[common[0]]['calls'])==35 and int(by[common[2]]['calls'])==32,f'{route} profile frequency/radial coverage differs')
        counters={(r.get('phase'),r.get('name')):int(r.get('value',0)) for r in p.get('counters',[])}
        need(counters.get((phase,'phase_c2_forward_spectrum_frequency_calls'))==35 and counters.get((phase,'phase_c2_unit_adjoint_frequency_calls'))==32 and counters.get((phase,'phase_c2_exact_zero_unit_adjoint_skips'))==3,f'{route} profile frequency counters differ')
        if route=='legacy-complex':
            need(int(by[gemm]['calls'])==147840 and counters.get((phase,'phase_e_angular_eigen_gemm_calls'))==147840,f'legacy complex GEMM count differs')
            need(counters.get((phase,'phase_e_angular_real_legendre_gemm_calls'),0)==0 and counters.get((phase,'phase_e_angular_route_real_legendre'),0)==0,'legacy route unexpectedly entered real-Legendre path')
        else:
            need(int(by[gemm]['calls'])==27552 and counters.get((phase,'phase_e_angular_real_legendre_gemm_calls'))==55104,f'real-Legendre GEMM call/product counts differ')
            need(counters.get((phase,'phase_e_angular_route_real_legendre'))==32,f'real-Legendre route marker count differs')
            need(int(by['phase_e_angular_coefficient_row_lookup_and_packing']['calls'])==1312 and int(by['phase_e_angular_material_evaluation_and_contraction']['calls'])==32,f'real-Legendre pack/material call counts differ')
            need(counters.get((phase,'phase_e_angular_eigen_gemm_calls'),0)==0 and counters.get((phase,'phase_e_angular_route_legacy_complex'),0)==0,'real-Legendre route unexpectedly entered legacy path')
        if route=='legacy-complex':need(counters.get((phase,'phase_e_angular_route_legacy_complex'))==32,f'legacy-complex route marker count differs')
        wall=p.get('wall_partition',{})
        for k in ('profile_pipeline_wall_s','master_exclusive_partition_s','residual_s','tolerance_s'):
            need(math.isfinite(float(wall.get(k,math.nan))),f'{route} wall partition lacks finite {k}')
        need(float(wall['profile_pipeline_wall_s'])>0 and float(wall['master_exclusive_partition_s'])>=0 and float(wall['tolerance_s'])>0 and abs(float(wall['residual_s']))<=float(wall['tolerance_s']),'profile wall partition fails')
        worker=sum(float(r['exclusive_s']) for r in raw if r.get('execution_domain')=='worker_elapsed_sum');root=by[common[0]]
        need(abs(worker-float(root['inclusive_s']))<=max(1e-6,1e-3*max(float(root['inclusive_s']),1e-9)),f'{route} worker profile closure fails')
        return {'gemm_scope':by[gemm],'worker_exclusive_sum_s':worker,'wall_partition':wall}
    profiles={'legacy_complex':profile(data('dspec-legacy-angle-profile'),'legacy-complex'),'real_legendre':profile(data('dspec-real-angle-profile'),'real-legendre')}
    selected(obj(data('dspec-legacy-angle-profile')/'result.json'),'legacy profile','legacy-complex',True,False)
    selected(obj(data('dspec-real-angle-profile')/'result.json'),'real profile','real-legendre',True,False)
    def timing(paths,prefix):
        ds=[obj(x/'result.json') for x in paths];vals=[float(x['elapsed_s']) for x in ds];ids=[f'{prefix}-{i}' for i in range(1,4)];walls=[float(jr[i]['elapsed_seconds']) for i in ids]
        route='legacy-complex' if 'legacy' in prefix else 'real-legendre'
        for d in ds:selected(d,prefix,route,True,False)
        for x in vals+walls:need(math.isfinite(x) and x>0,f'{prefix} driver/supervisor time invalid')
        return {'driver_elapsed_s':vals,'median_s':statistics.median(vals),'min_s':min(vals),'max_s':max(vals),'supervisor_wall_s':walls,'supervisor_median_s':statistics.median(walls),'peak_rss_gib_userspace':max(float(jr[i]['peak_observed_rss_gib_userspace']) for i in ids),'peak_threads':max(int(jr[i]['peak_observed_threads']) for i in ids)}
    lt=timing(legacy,'dspec-legacy-angle');rt=timing(real,'dspec-real-angle')
    ytimes=[]
    for i in range(1,4):
        p=data(f'yspec-native-{i}');d=obj(p/'result.json');need(d.get('backend')=='yspec' and d.get('degree')==200 and d.get('nq')==6 and d.get('frequency_count')==35,'YSpec workload differs')
        need(d.get('workers')==8 and d.get('frequency_workers')==1 and d.get('solve_workers')==8 and d.get('first_bin')==7 and d.get('last_bin')==41,'YSpec refresh worker/grid configuration differs')
        gate(compare_science(p,a.frozen_y),1e-11,'YSpec-v-frozen-Y');v=float(d['elapsed_s']);need(math.isfinite(v) and v>0,'YSpec time invalid');ytimes.append(v)
    l3=data('dspec-l3-fd-audit');ld=obj(l3/'result.json');need(ld.get('degree')==3 and ld.get('nq')==8 and ld.get('native_angular_route')=='real-legendre' and ld.get('native_radial_route')=='batched-direct-index' and ld.get('phase_e_angular_oracle') is True,'L3 did not exercise actual selected route/oracle')
    need(ld.get('native_side_count')==157 and ld.get('along_angles')==91 and ld.get('across_angles')==0 and ld.get('frequency_count')==35 and abs(float(ld.get('refinement',-1))-0.05)<1e-14 and ld.get('phase_c_zero_bin_skips') is False and ld.get('native_derivative_weights_mode')=='native-dref-over-jacobian','L3 grid/FD route differs')
    need(ld.get('frequency_workers')==1 and ld.get('solve_workers')==1 and ld.get('phase_c_direct_adjoint') is True and abs(float(ld.get('gravity_max_abs_nd_error',math.inf)))<=1e-10 and abs(float(ld.get('time_shift_adjoint_identity',0))-1)<=1e-10 and float(ld.get('phase_d_actual_temporary_peak_bytes',math.inf))<=32*1024*1024,'L3 angular control/resources differ')
    l3proof=rows(l3/'phase-e-angular-route-equivalence.csv');need(len(l3proof)==35 and {(int(r['bin']),int(r['section'])) for r in l3proof}=={(b,0) for b in range(7,42)},'L3 angular proof does not cover its 35-bin one-section inventory')
    for r in l3proof:
        for k in angular_max:
            v=float(r[k]);n=float(r[angular_norms[k]]);need(math.isfinite(v) and 0<=v<=1e-11 and math.isfinite(n) and n>0,f'L3 angular {k}/norm fails at bin {r["bin"]}')
    fd=obj(a.l3_fd_check);v=float(fd.get('fd_relative_error',math.nan));need(fd.get('status')=='pass' and math.isfinite(v) and 0<=v<=2e-4,'L3 finite-difference check failed')
    l3parity=compare_science(l3,a.frozen_l3);gate(l3parity,1e-11,'L3-v-frozen')
    acc=plan.get('accuracy_anchor',{});need(acc.get('status')=='pass' and acc.get('requested_accuracy_met') is True and acc.get('accuracy_report_sha256')=='013273941560685863015e48aeef42daeada129b3174ff47217637f21e605a95' and float(acc['max_solid_complex_field_relative_l2'])<1e-3 and float(acc['max_solid_full35_kernel_relative_l2'])<1e-3,'hash-bound NQ10 accuracy support is missing or changed')
    report={'status':'pass','stage':'Phase E Stage 2 angular reconstruction','timing':{'legacy_complex':lt,'real_legendre':rt,'yspec_elapsed_s':ytimes,'yspec_median_s':statistics.median(ytimes),'elapsed_reduction_fraction':(lt['median_s']-rt['median_s'])/lt['median_s'],'profiles':profiles},'same_workspace_angular_oracle':{'rows':len(angular),'max_relative_l2':angular_max},'radial_oracle_bin_count':len(radial),'all_bin_full_pipeline_parity_vs_final_radial':audit_parity,'scalar_pipeline_parity_rows':len(scalar),'l3':{'route':ld['native_angular_route'],'angular_oracle_rows':len(l3proof),'fd_relative_error':v,'parity_vs_frozen':l3parity},'nq10_accuracy_anchor':acc,'resource_caps':{'supervisor_threads':20,'job_rss_gib':16}}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps(report,sort_keys=True))
if __name__=='__main__':main()
