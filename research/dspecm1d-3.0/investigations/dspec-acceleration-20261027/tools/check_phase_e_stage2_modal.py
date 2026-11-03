#!/usr/bin/env python3
"""Check modal-operator parity, exact per-key solve evidence, controls and timings."""
from __future__ import annotations
import argparse,csv,hashlib,json,math,statistics,struct
from pathlib import Path

def need(ok,msg):
    if not ok: raise ValueError(msg)
def obj(p): return json.loads(Path(p).read_text())
def rows(p):
    with Path(p).open(newline='') as f:return list(csv.DictReader(f))
def rel(a,b):
    need(len(a)==len(b),'vectors differ in length');e=n=0.0
    for x,y in zip(a,b):
        x=float(x);y=float(y);need(math.isfinite(x) and math.isfinite(y),'nonfinite comparison value');e+=(x-y)**2;n+=y*y
    return math.sqrt(e/n) if n else (0.0 if e==0 else math.inf)
def compare_csv(a,b,keys,channels):
    x,y=rows(a),rows(b);need(len(x)==len(y),f'row count differs: {a} {b}')
    for i,(u,v) in enumerate(zip(x,y)):
        for k in keys:need(u[k]==v[k],f'key/grid differs: {a}:{i}:{k}')
    return {k:rel([r[k] for r in x],[r[k] for r in y]) for k in channels}
def compare_bin(a,b):
    x=Path(a).read_bytes();y=Path(b).read_bytes();need(len(x)==len(y) and len(x)%8==0,f'binary differs: {a} {b}')
    return rel(struct.unpack('<'+str(len(x)//8)+'d',x),struct.unpack('<'+str(len(y)//8)+'d',y))
def science(a,b):
    a=Path(a);b=Path(b)
    m={'trace':compare_csv(a/'trace.csv',b/'trace.csv',['time_s','window'],['velocity_m_s','rate_m_s2','adjoint_source']),
       'spectrum':compare_csv(a/'spectrum.csv',b/'spectrum.csv',['bin','frequency_hz'],['real_m','imag_m']),
       'components':compare_csv(a/'components.csv',b/'components.csv',['bin'],['up_re','up_im','north_re','north_im','east_re','east_im'])}
    for fn in ('nodes.csv','frequency-grid.csv','along-angles.csv','across-angles.csv','mesh-elements.csv'):
        x,y=a/fn,b/fn;need(x.exists()==y.exists(),f'geometry inventory differs: {fn}')
        if x.exists():need(rows(x)==rows(y),f'geometry differs: {fn}')
    for fn in ('along.bin','across.bin'):
        x,y=a/fn,b/fn;need(x.exists()==y.exists(),f'section inventory differs: {fn}')
        if x.exists():m[fn]=compare_bin(x,y)
    return m
def gate(m,limit,label):
    for group,vals in m.items():
        for k,v in (vals.items() if isinstance(vals,dict) else [(group,vals)]):
            need(math.isfinite(float(v)) and 0<=float(v)<=limit,f'{label}.{k}={v} exceeds {limit}')
def clean_state(state,plan_path,plan):
    need(state.get('status')=='completed','producer is not completed')
    need(state.get('plan_sha256')==hashlib.sha256(Path(plan_path).read_bytes()).hexdigest(),'state/run-plan hash mismatch')
    jobs={j['id']:j for j in plan['jobs']}; receipts={j['id']:j for j in state.get('jobs',[])}
    need(set(jobs)==set(receipts),'producer job inventory differs')
    for jid,r in receipts.items():
        need(r.get('status')=='completed' and r.get('returncode')==0 and r.get('inputs_unchanged') is True and not r.get('failure_reasons') and r.get('cleanup',{}).get('verified_gone') is True,f'{jid} receipt not clean')
        need(int(r.get('peak_observed_threads',99))<=20 and float(r.get('peak_observed_rss_gib_userspace',99))<=16,f'{jid} resource cap exceeded')
        need(math.isfinite(float(r.get('elapsed_seconds',math.nan))) and float(r['elapsed_seconds'])>0,f'{jid} supervisor wall invalid')
    return jobs,receipts
def modal_proof(path,expected_bins,degrees,basis_mode='absent'):
    with Path(path).open(newline='') as f: header=next(csv.reader(f))
    expected_header='bin,family,degree,rows,nnz,operator_ref_norm,operator_abs_l2,operator_rel_l2,operator_max_abs,forward_ref_norm,forward_abs_l2,forward_rel_l2,forward_zero,adjoint_ref_norm,adjoint_abs_l2,adjoint_rel_l2,adjoint_zero,basis_adjoint_ref_norm,basis_adjoint_abs_l2,basis_adjoint_rel_l2,basis_adjoint_zero,probe_forward_ref_norm,probe_forward_abs_l2,probe_forward_rel_l2,probe_adjoint_ref_norm,probe_adjoint_abs_l2,probe_adjoint_rel_l2,physical_new_forward_raw_residual,physical_new_forward_backward_error,physical_ref_forward_raw_residual,physical_ref_forward_backward_error,physical_new_adjoint_raw_residual,physical_new_adjoint_backward_error,physical_ref_adjoint_raw_residual,physical_ref_adjoint_backward_error,basis_new_adjoint_raw_residual,basis_new_adjoint_backward_error,basis_ref_adjoint_raw_residual,basis_ref_adjoint_backward_error,probe_new_forward_raw_residual,probe_new_forward_backward_error,probe_ref_forward_raw_residual,probe_ref_forward_backward_error,probe_new_adjoint_raw_residual,probe_new_adjoint_backward_error,probe_ref_adjoint_raw_residual,probe_ref_adjoint_backward_error,physical_new_forward_h_times_x_norm,physical_new_forward_rhs_norm,physical_new_forward_backward_denominator,physical_ref_forward_h_times_x_norm,physical_ref_forward_rhs_norm,physical_ref_forward_backward_denominator,physical_new_adjoint_h_times_x_norm,physical_new_adjoint_rhs_norm,physical_new_adjoint_backward_denominator,physical_ref_adjoint_h_times_x_norm,physical_ref_adjoint_rhs_norm,physical_ref_adjoint_backward_denominator,basis_new_adjoint_h_times_x_norm,basis_new_adjoint_rhs_norm,basis_new_adjoint_backward_denominator,basis_ref_adjoint_h_times_x_norm,basis_ref_adjoint_rhs_norm,basis_ref_adjoint_backward_denominator,probe_new_forward_h_times_x_norm,probe_new_forward_rhs_norm,probe_new_forward_backward_denominator,probe_ref_forward_h_times_x_norm,probe_ref_forward_rhs_norm,probe_ref_forward_backward_denominator,probe_new_adjoint_h_times_x_norm,probe_new_adjoint_rhs_norm,probe_new_adjoint_backward_denominator,probe_ref_adjoint_h_times_x_norm,probe_ref_adjoint_rhs_norm,probe_ref_adjoint_backward_denominator,new_dual_abs_error,new_dual_scalar_relative_diagnostic,ref_dual_abs_error,ref_dual_scalar_relative_diagnostic,new_dual_scalar_denominator,ref_dual_scalar_denominator,new_h_frobenius_norm,new_probe_forward_norm,new_probe_adjoint_norm,new_dual_normwise_denominator,ref_dual_normwise_denominator,new_dual_normwise_relative,ref_dual_normwise_relative'.split(',')
    need(len(header)==90 and len(set(header))==90 and header==expected_header,f'{path} header/schema differs from the reviewed 90-column modal proof')
    rr=rows(path);keys=set();family_norms={};zero_counts={'forward':0,'adjoint':0,'basis_adjoint':0}
    relcols=[x for x in rr[0] if x.endswith('_rel_l2') or x.endswith('_backward_error') or x.endswith('_relative_error') or x in ('new_dual_normwise_relative','ref_dual_normwise_relative')]
    rawcols=[x for x in rr[0] if x.endswith('_abs_l2') or x.endswith('_raw_residual') or x.endswith('_abs_error') or x.endswith('_max_abs')]
    need(relcols and rawcols,'modal proof lacks relative/raw diagnostics')
    for r in rr:
        b=int(r['bin']);fam=int(r['family']);d=int(r['degree']);rows_n=int(r['rows']);nnz=int(r['nnz'])
        need(b in expected_bins and rows_n>0 and nnz>0,'bad bin/operator shape in modal proof')
        need((fam==0 and d==0) or (fam in (1,2) and 1<=d<=degrees),f'invalid family/degree key {(fam,d)}')
        key=(b,fam,d);need(key not in keys,f'duplicate modal proof key {key}');keys.add(key)
        for c in header:
            if c in ('forward_zero','adjoint_zero','basis_adjoint_zero'):
                need(r[c] in ('0','1'),f'invalid zero classification {c} at {key}');continue
            v=float(r[c]);need(math.isfinite(v) and v>=0,f'nonfinite/negative {c} at {key}')
        for c in relcols:need(float(r[c])<=1e-11,f'{c} exceeds 1e-11 at {key}')
        for c in ('operator_ref_norm','probe_forward_ref_norm','probe_adjoint_ref_norm'):
            v=float(r[c]);need(math.isfinite(v) and v>0,f'{c} must be positive at {key}')
        for c in ('forward_zero','adjoint_zero','basis_adjoint_zero'):
            need(r[c] in ('0','1'),f'invalid zero classification {c} at {key}')
        for stem in ('forward','adjoint','basis_adjoint'):
            zero=r[stem+'_zero']=='1';ref=float(r[stem+'_ref_norm']);err=float(r[stem+'_abs_l2']);rrv=float(r[stem+'_rel_l2'])
            if zero:need(ref==0 and err==0 and rrv==0,f'{stem} exact-zero classification inconsistent at {key}');zero_counts[stem]+=1
            else:need(ref>0 and err/ref<=1e-11 and abs(rrv-err/ref)<=max(1e-15,1e-10*rrv),f'{stem} relative denominator/limit mismatch at {key}')
            family_norms.setdefault(fam,{}).setdefault(stem,0.0);family_norms[fam][stem]+=ref*ref
        for prefix,refcol,abscol,relcol in (('operator','operator_ref_norm','operator_abs_l2','operator_rel_l2'),('probe_forward','probe_forward_ref_norm','probe_forward_abs_l2','probe_forward_rel_l2'),('probe_adjoint','probe_adjoint_ref_norm','probe_adjoint_abs_l2','probe_adjoint_rel_l2')):
            ref=float(r[refcol]);err=float(r[abscol]);rv=float(r[relcol]);calc=err/ref
            need(calc<=1e-11 and abs(calc-rv)<=max(1e-15,1e-10*rv),f'{prefix} relative metric denominator/limit mismatch at {key}')
        # Independently recompute every backward error from its published raw terms.
        for prefix in ('physical_new_forward','physical_ref_forward','physical_new_adjoint','physical_ref_adjoint','basis_new_adjoint','basis_ref_adjoint','probe_new_forward','probe_ref_forward','probe_new_adjoint','probe_ref_adjoint'):
            raw=float(r[prefix+'_raw_residual']);be=float(r[prefix+'_backward_error'])
            base=prefix.removesuffix('_raw_residual')
            denom=float(r[base+'_backward_denominator'])
            recomputed=raw/denom if denom>0 else (0.0 if raw==0 else math.inf)
            need(math.isfinite(recomputed) and recomputed<=1e-11 and abs(recomputed-be)<=max(1e-15,1e-10*be),f'{prefix} backward-error denominator/limit mismatch at {key}')
            need(be<=1e-11,f'{prefix} backward error exceeds 1e-11 at {key}')
            hx=float(r[base+'_h_times_x_norm']);rhs=float(r[base+'_rhs_norm'])
            need(abs((hx+rhs)-denom)<=max(1e-14,1e-12*denom),f'{prefix} denominator terms do not sum at {key}')
        for stem in ('new_dual','ref_dual'):
            scalar_error=float(r[stem+'_abs_error']);scalar_den=float(r[stem+'_scalar_denominator']);scalar=float(r[stem+'_scalar_relative_diagnostic'])
            scalar_calc=scalar_error/scalar_den if scalar_den>0 else (0.0 if scalar_error==0 else math.inf)
            need(math.isfinite(scalar_calc) and abs(scalar_calc-scalar)<=max(1e-15,1e-10*scalar),f'{stem} scalar diagnostic is not reproducible at {key}')
        new_den=float(r['new_dual_normwise_denominator']);ref_den=float(r['ref_dual_normwise_denominator'])
        new_calc=float(r['new_dual_abs_error'])/new_den if new_den>0 else (0.0 if float(r['new_dual_abs_error'])==0 else math.inf)
        ref_calc=float(r['ref_dual_abs_error'])/ref_den if ref_den>0 else (0.0 if float(r['ref_dual_abs_error'])==0 else math.inf)
        need(abs(new_den-float(r['new_h_frobenius_norm'])*float(r['new_probe_forward_norm'])*float(r['new_probe_adjoint_norm']))<=max(1e-14,1e-12*new_den),f'new normwise dual denominator terms disagree at {key}')
        need(abs(ref_den-float(r['operator_ref_norm'])*float(r['probe_forward_ref_norm'])*float(r['probe_adjoint_ref_norm']))<=max(1e-14,1e-12*ref_den),f'reference normwise dual denominator terms disagree at {key}')
        for name,value in (('new',new_calc),('ref',ref_calc)):
            reported=float(r[name+'_dual_normwise_relative'])
            need(math.isfinite(value) and value<=1e-11 and abs(value-reported)<=max(1e-15,1e-10*reported),f'{name} normwise dual identity exceeds 1e-11 or is not reproducible at {key}')
    expected={(b,0,0) for b in expected_bins}|{(b,f,d) for b in expected_bins for f in (1,2) for d in range(1,degrees+1)}
    need(keys==expected,f'modal key inventory mismatch: got{len(keys)} expected{len(expected)}')
    need(all(fam in family_norms and family_norms[fam].get('forward',0)>0 for fam in (0,1,2)),f'physical forward coverage absent in a solve family')
    need(all(family_norms[fam].get('adjoint',0)>0 for fam in (1,2)),f'weighted physical adjoint coverage absent in active S/T families')
    if basis_mode=='absent':
        need(zero_counts['basis_adjoint']==len(rr),f'basis-adjoint zero classification differs from the direct-adjoint production route ({zero_counts["basis_adjoint"]}/{len(rr)})')
    elif basis_mode=='active-st':
        need(zero_counts['basis_adjoint']==len(expected_bins),f'L3 basis-adjoint exact-zero count differs from R0-only structural zeros ({zero_counts["basis_adjoint"]}/{len(expected_bins)})')
        need(all(family_norms[fam].get('basis_adjoint',0)>0 for fam in (1,2)),'L3 unit-basis adjoint coverage absent in active S/T families')
    else: raise ValueError(f'unknown modal basis proof mode {basis_mode}')
    return {'rows':len(rr),'unique_keys':len(keys),'physical_exact_zero_counts':zero_counts,'physical_reference_norm2_by_family':family_norms,'max_rel_or_backward_error':max(float(r[c]) for r in rr for c in relcols)}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run-root',type=Path,required=True);ap.add_argument('--run-plan',type=Path,required=True);ap.add_argument('--supervisor-state',type=Path,required=True)
    for n in ('frozen-d','frozen-y','frozen-l3','l3-fd-check','accuracy-report','out'):ap.add_argument('--'+n,type=Path,required=True)
    a=ap.parse_args();plan=obj(a.run_plan);need(plan.get('thread_cap')==20,'numerical thread cap differs');jobs,state=clean_state(obj(a.supervisor_state),a.run_plan,plan);run=a.run_root
    ids=[f'dspec-{r}-modal-{i}' for i in range(1,4) for r in ('legacy','cached')]
    for route,jids in (('legacy', [f'dspec-legacy-modal-{i}' for i in range(1,4)]+['dspec-legacy-modal-profile']),('cached-union',[f'dspec-cached-modal-{i}' for i in range(1,4)]+['dspec-cached-modal-profile'])):
        for jid in jids:
            cmd=' '.join(jobs[jid]['command'])
            for setting in ('DSFAST_RADIAL_ROUTE=batched','DSFAST_ANGULAR_ROUTE=real-legendre',f'DSPEC_MODAL_OPERATOR_ROUTE={route}','DSPEC_MODAL_OPERATOR_ORACLE=0'):
                need(setting in cmd,f'{jid} missing {setting}')
            need('DSFAST_ANGULAR_ORACLE=0' in cmd and 'DSPEC_PHASE_C_VALIDATE_DIRECT_ADJOINT=0' in cmd and 'DSPEC_PHASE_C_VALIDATE_RECONSTRUCTION=0' in cmd and 'DSPEC_PHASE_C_SKIP_ZERO_KERNEL_BINS=1' in cmd,f'{jid} timed route has oracle/skip mismatch')
    for jid in ('dspec-modal-all-bin-audit','dspec-l3-modal-fd-audit'):
        cmd=' '.join(jobs[jid]['command'])
        for setting in ('DSFAST_RADIAL_ROUTE=batched','DSFAST_ANGULAR_ROUTE=real-legendre','DSPEC_MODAL_OPERATOR_ROUTE=cached-union','DSPEC_MODAL_OPERATOR_ORACLE=1','DSFAST_ANGULAR_ORACLE=1','DSPEC_PHASE_C_SKIP_ZERO_KERNEL_BINS=0'):
            need(setting in cmd,f'{jid} missing proof control {setting}')
        need(any('modal-operator-oracle.csv' in x for x in jobs[jid]['outputs']),f'{jid} does not declare modal oracle output')
    def data(j):return run/j/'data'
    def selected(p,label,route,oracle,degree=200,nq=6,workers=8,skip=True,validators=False):
        d=obj(p/'result.json')
        need(d.get('backend')=='dspec' and d.get('degree')==degree and d.get('nq')==nq and d.get('frequency_count')==35 and d.get('first_bin')==7 and d.get('last_bin')==41,f'{label} workload/grid differs')
        h=.02035382429229197 if degree==200 else .05
        need(abs(float(d.get('refinement',-1))-h)<1e-14 and d.get('native_side_count')==(208 if degree==200 else 157),f'{label} mesh differs')
        need(d.get('native_radial_route')=='batched-direct-index' and d.get('native_angular_route')=='real-legendre' and d.get('modal_operator_route')==route,f'{label} route metadata differs')
        need(d.get('modal_operator_oracle') is oracle and d.get('native_derivative_weights_mode')=='native-dref-over-jacobian',f'{label} oracle/derivative metadata differs')
        need(d.get('phase_e_angular_oracle') is oracle,f'{label} angular oracle mode differs')
        need((d.get('along_angles'),d.get('across_angles'))==((135,19) if degree==200 else (91,0)),f'{label} angle inventory differs')
        scratch=float(d.get('phase_d_actual_temporary_peak_bytes',math.nan));need(math.isfinite(scratch) and 0<scratch<=32*1024*1024,f'{label} angular working-set bound differs')
        need(d.get('frequency_workers')==workers and d.get('solve_workers')==1 and d.get('phase_c_direct_adjoint') is True,f'{label} worker/direct-adjoint config differs')
        need(d.get('phase_c_zero_bin_skips') is skip and abs(float(d.get('gravity_max_abs_nd_error',math.inf)))<=1e-10 and abs(float(d.get('time_shift_adjoint_identity',0))-1)<=1e-10,f'{label} control check failed')
        need(d.get('phase_c_validate_direct_adjoint') is validators and d.get('phase_c_validate_reconstruction') is validators,f'{label} validator controls differ')
        if route=='cached-union':
            need(d.get('modal_operator_cache_entries')==(401 if degree==200 else 7),f'{label} cache entry count differs')
            need(d.get('modal_operator_shared_cap_bytes')==256*1024*1024 and d.get('modal_operator_worker_cap_bytes')==32*1024*1024 and 0<int(d.get('modal_operator_cache_bytes',0))<=256*1024*1024 and 0<int(d.get('modal_operator_worker_scratch_high_water_bytes',-1))<=32*1024*1024,f'{label} cache/scratch cap fails')
        else: need(d.get('modal_operator_cache_entries')==0 and d.get('modal_operator_cache_bytes')==0,f'{label} legacy baseline constructed modal cache')
        return d
    aud=data('dspec-modal-all-bin-audit');ad=selected(aud,'all-bin modal audit','cached-union',True,skip=False)
    need(ad.get('phase_c_validate_reconstruction') is False and ad.get('phase_c_validate_direct_adjoint') is False,'all-bin modal audit enables expensive generic validators')
    proof=modal_proof(aud/'modal-operator-oracle.csv',set(range(7,42)),200);need(proof['rows']==14035,'L200 modal proof row count differs')
    def table_proofs(directory,sections,label,expected_counts):
        scalar=rows(directory/'phase-d-batched-scalar-equivalence.csv')
        need(len(scalar)==35*sections and {(int(r['bin']),int(r['section'])) for r in scalar}=={(b,s) for b in range(7,42) for s in range(sections)},f'{label} scalar proof inventory differs')
        for r in scalar:
            for c in ('forward_field_rel_l2','adjoint_field_rel_l2','complex_product_rel_l2'):
                v=float(r[c]);need(math.isfinite(v) and 0<=v<=1e-11,f'{label} scalar {c} fails at bin {r["bin"]}')
            for c in ('forward_scalar_norm2','adjoint_scalar_norm2','product_scalar_norm2'):
                v=float(r[c]);need(math.isfinite(v) and v>0,f'{label} scalar reference norm {c} is not positive')
        angular=rows(directory/'phase-e-angular-route-equivalence.csv')
        need(len(angular)==35*sections and {(int(r['bin']),int(r['section'])) for r in angular}=={(b,s) for b in range(7,42) for s in range(sections)},f'{label} angular proof inventory differs')
        normmap={'forward_field_rel_l2':'forward_legacy_norm2','adjoint_field_rel_l2':'adjoint_legacy_norm2','complex_product_rel_l2':'product_legacy_norm2'}
        for r in angular:
            for c,nc in normmap.items():
                v=float(r[c]);n=float(r[nc]);need(math.isfinite(v) and 0<=v<=1e-11 and math.isfinite(n) and n>0,f'{label} angular proof {c}/norm invalid')
        radial=[]
        for b in range(7,42):
            with (directory/f'radial-oracle.bin{b}.csv').open(newline='') as f:q=list(csv.reader(f))
            need(len(q)==1 and len(q[0])==16 and int(q[0][0])==b and q[0][-1]=='generic-and-batched-keyed-modes',f'{label} radial proof missing or malformed at bin {b}')
            row=q[0];counts=(int(row[1]),int(row[3]),int(row[4]),int(row[5]),int(row[6]))
            need(counts[0]>0 and counts[1]+counts[2]==counts[0] and counts[3]+counts[4]==counts[0] and all(v>0 for v in counts[1:]),f'{label} radial key categories differ at bin {b}')
            need(int(row[2])<=256*1024*1024 and all(math.isfinite(float(v)) and 0<=float(v)<=1e-11 for v in row[7:15]),f'{label} radial proof bound/parity fails at bin {b}')
            radial.append(counts)
        need(all(x==radial[0] for x in radial) and radial[0]==expected_counts,f'{label} radial inventory changed or differs from frozen baseline')
        return {'scalar_rows':len(scalar),'angular_rows':len(angular),'radial_bins':len(radial)}
    l200proofs=table_proofs(aud,2,'L200 audit',(620440,387432,233008,44800,575640))
    # Native outputs must remain equal to the final validated angular candidate.
    audit_parity=science(aud,a.frozen_d);gate(audit_parity,1e-11,'all-bin modal-v-final-angular')
    legacy=[data(f'dspec-legacy-modal-{i}') for i in range(1,4)];cached=[data(f'dspec-cached-modal-{i}') for i in range(1,4)]
    for i,(x,y) in enumerate(zip(legacy,cached),1):
        selected(x,f'legacy timing {i}','legacy',False);selected(y,f'cached timing {i}','cached-union',False)
        gate(science(x,y),1e-11,f'legacy-v-cached-{i}')
        gate(science(x,a.frozen_d),1e-11,f'legacy-v-final-angular-{i}');gate(science(y,a.frozen_d),1e-11,f'cached-v-final-angular-{i}')
    def profile(jid,route):
        p=obj(data(jid)/'profile.json');phase='combined_forward_unit_adjoint_reconstruction_contraction';raw=[r for r in p.get('rows',[]) if r.get('phase')==phase]
        need(raw,'profile has no combined C2 phase rows')
        for r in raw:
            inc=float(r.get('inclusive_s',math.nan));exc=float(r.get('exclusive_s',math.nan));need(math.isfinite(inc) and math.isfinite(exc) and 0<=exc<=inc+1e-9,f'{jid} invalid profile row')
        by={r.get('scope_name','').removeprefix('frequency_worker.'):r for r in raw}
        for name in ('combined_frequency_workspace_forward_unit_adjoint_and_reconstruction','dynamic_H_assembly','solver_factorization_dispatch'):
            expected=35 if name=='combined_frequency_workspace_forward_unit_adjoint_and_reconstruction' else 14035
            need(name in by and int(by[name].get('calls',0))==expected,f'{jid} missing/incorrect {name} count')
        counters={(x.get('phase'),x.get('name')):int(x.get('value',0)) for x in p.get('counters',[])}
        key=(phase,'modal_operator_cached_union_calls' if route=='cached-union' else 'modal_operator_legacy_calls')
        need(counters.get(key)==14035,f'{jid} modal route count differs')
        cache_count=sum(int(x.get('value',0)) for x in p.get('counters',[]) if x.get('name')=='modal_component_cache_entries')
        need(cache_count==(401 if route=='cached-union' else 0),f'{jid} cache entry counter differs')
        all_rows=p.get('rows',[]);cold=[r for r in all_rows if r.get('scope_name','').removeprefix('frequency_worker.')=='modal_component_cache_cold_build']
        if route=='cached-union':need(len(cold)==1 and int(cold[0].get('calls',0))==1 and math.isfinite(float(cold[0].get('inclusive_s',math.nan))) and float(cold[0]['inclusive_s'])>0,'cold cache build is not separately timed')
        else:need(not cold,'legacy route paid cache construction')
        wall=p.get('wall_partition',{})
        for k in ('profile_pipeline_wall_s','master_exclusive_partition_s','residual_s','tolerance_s'):need(math.isfinite(float(wall.get(k,math.nan))),f'{jid} profile closure {k} invalid')
        need(wall.get('worker_elapsed_excluded') is True and float(wall['profile_pipeline_wall_s'])>0 and float(wall['master_exclusive_partition_s'])>=0 and float(wall['tolerance_s'])>0 and abs(float(wall['residual_s']))<=float(wall['tolerance_s']) and abs(float(wall['profile_pipeline_wall_s'])-float(wall['master_exclusive_partition_s']))<=float(wall['tolerance_s']) and abs((float(wall['master_exclusive_partition_s'])-float(wall['profile_pipeline_wall_s']))-float(wall['residual_s']))<=float(wall['tolerance_s']),f'{jid} master/profile wall closure fails')
        return {'modal_assembly':by['dynamic_H_assembly'],'factorization':by['solver_factorization_dispatch'],'wall_partition':wall,'cold_cache_build':cold[0] if cold else None}
    prof={'legacy':profile('dspec-legacy-modal-profile','legacy'),'cached_union':profile('dspec-cached-modal-profile','cached-union')}
    for r,jid in (('legacy','dspec-legacy-modal-profile'),('cached-union','dspec-cached-modal-profile')):
        selected(data(jid),jid,r,False)
        gate(science(data(jid),a.frozen_d),1e-11,f'{jid}-v-final-angular')
    def times(paths,prefix):
        ds=[obj(p/'result.json') for p in paths];v=[float(x['elapsed_s']) for x in ds];need(all(math.isfinite(x) and x>0 for x in v),f'{prefix} driver elapsed invalid')
        return {'driver_elapsed_s':v,'median_s':statistics.median(v),'min_s':min(v),'max_s':max(v)}
    lt=times(legacy,'legacy');ct=times(cached,'cached')
    ys=[]
    for i in range(1,4):
        d=obj(data(f'yspec-native-{i}')/'result.json');need(d.get('backend')=='yspec' and d.get('degree')==200 and d.get('nq')==6 and d.get('frequency_count')==35 and d.get('workers')==8 and d.get('frequency_workers')==1 and d.get('solve_workers')==8,f'YSpec workload {i} differs')
        gate(science(data(f'yspec-native-{i}'),a.frozen_y),1e-11,f'YSpec-{i}-v-frozen');ys.append(float(d['elapsed_s']))
    need(all(math.isfinite(v) and v>0 for v in ys),'YSpec time invalid')
    l3=data('dspec-l3-modal-fd-audit');l3d=selected(l3,'L3 modal audit','cached-union',True,degree=3,nq=8,workers=1,skip=False,validators=True)
    need(l3d.get('across_angles')==0 and l3d.get('along_angles')==91,'L3 angular inventory differs')
    l3proof=modal_proof(l3/'modal-operator-oracle.csv',set(range(7,42)),3,basis_mode='active-st');need(l3proof['rows']==245,'L3 modal proof row count differs')
    l3tables=table_proofs(l3,1,'L3 audit',(6760,3966,2794,504,6256))
    gate(science(l3,a.frozen_l3),1e-11,'L3-v-frozen-angular')
    fd=obj(a.l3_fd_check);fderr=float(fd.get('fd_relative_error',math.nan));need(fd.get('status')=='pass' and math.isfinite(fderr) and 0<=fderr<=2e-4,'L3 FD control failed')
    acc=plan.get('accuracy_anchor',{});anchor=obj(a.accuracy_report);anchor_sha=hashlib.sha256(a.accuracy_report.read_bytes()).hexdigest()
    need(acc.get('status')=='pass' and acc.get('requested_accuracy_met') is True and acc.get('accuracy_report_sha256')==anchor_sha=='013273941560685863015e48aeef42daeada129b3174ff47217637f21e605a95','hash-bound NQ10 accuracy evidence missing or changed')
    mesh=anchor.get('selected_mesh',{});accuracy=anchor.get('accuracy',{})
    need(anchor.get('status')=='pass' and mesh.get('nq')==6 and mesh.get('mesh_scale')==.9 and abs(float(mesh.get('h',-1))-.02035382429229197)<1e-14 and mesh.get('native_sides')==208 and mesh.get('frequency_bins')==35 and mesh.get('requested_accuracy_met') is True and float(mesh.get('requested_relative_accuracy',math.inf))==.001,'NQ10 accuracy anchor selected mesh/workload differs')
    need(accuracy.get('max_solid_complex_field_relative_l2',math.inf)<.001 and accuracy.get('max_solid_full35_kernel_relative_l2',math.inf)<.001,'existing NQ10-supported solid error budget exceeds 0.1 percent')
    report={'status':'pass','stage':'Phase E Stage 2 modal operator assembly','matched_timing':{'legacy':lt,'cached_union':ct,'yspec_elapsed_s':ys,'yspec_median_s':statistics.median(ys),'elapsed_reduction_fraction':(lt['median_s']-ct['median_s'])/lt['median_s'],'profiles':prof},'operator_oracle':proof,'table_and_radial_proofs':l200proofs,'L3_operator_oracle':l3proof,'L3_table_and_radial_proofs':l3tables,'L3_fd_relative_error':fderr,'all_bin_parity_vs_final_angular':audit_parity,'nq10_accuracy_anchor':{'sha256':anchor_sha,'selected_mesh':mesh,'solid_accuracy':{'max_solid_complex_field_relative_l2':accuracy['max_solid_complex_field_relative_l2'],'max_solid_full35_kernel_relative_l2':accuracy['max_solid_full35_kernel_relative_l2']}},'resource_caps':{'threads':20,'rss_gib':16}}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps(report,sort_keys=True))
if __name__=='__main__':main()
