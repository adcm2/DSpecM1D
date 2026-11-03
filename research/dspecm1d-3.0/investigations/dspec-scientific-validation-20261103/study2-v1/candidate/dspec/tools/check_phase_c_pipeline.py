#!/usr/bin/env python3
"""Compact Phase C equivalence, exact-skip, and native/profile parity gate."""
from __future__ import annotations
import argparse, csv, hashlib, json, math, struct
from pathlib import Path

EXPECTED_SKIP={7,8,41}

def require(ok,msg):
    if not ok: raise ValueError(msg)

def read_json(path):
    d=json.loads(Path(path).read_text())
    def walk(x):
        if isinstance(x,bool) or x is None or isinstance(x,str): return
        if isinstance(x,(int,float)): require(math.isfinite(float(x)),f"nonfinite JSON value at {path}")
        elif isinstance(x,list):
            for v in x:walk(v)
        elif isinstance(x,dict):
            for v in x.values():walk(v)
    walk(d);return d

def rows(path):
    with Path(path).open(newline="") as f:return list(csv.DictReader(f))

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def common_field_metrics(path):
    data=rows(path);groups={}
    require(data,"empty common field comparison")
    for r in data:
        key=(int(r['bin']),int(r['section']),r['location_class'],int(r['fluid']))
        g=groups.setdefault(key,[0.0]*4+[0])
        for role in ('forward','adjoint'):
            err=norm=0.0
            for c in range(6):
                x=complex(float(r[f'{role}_candidate{c}_re']),float(r[f'{role}_candidate{c}_im']))
                y=complex(float(r[f'{role}_ref{c}_re']),float(r[f'{role}_ref{c}_im']))
                require(math.isfinite(x.real) and math.isfinite(x.imag) and math.isfinite(y.real) and math.isfinite(y.imag),"nonfinite common field")
                err+=abs(x-y)**2;norm+=abs(y)**2
            i=0 if role=='forward' else 2;g[i]+=err;g[i+1]+=norm
        g[4]+=1
    return [{'bin':k[0],'section':k[1],'location_class':k[2],'fluid':k[3],
             'count':v[4],'forward_relative_l2':math.sqrt(v[0]/v[1]) if v[1] else (0.0 if v[0]==0 else math.inf),
             'adjoint_relative_l2':math.sqrt(v[2]/v[3]) if v[3] else (0.0 if v[2]==0 else math.inf)}
            for k,v in sorted(groups.items())]

def verify_case(d,label):
    for k,v in {'degree':200,'nq':10,'refinement':0.025,'frequency_count':35,'first_bin':7,'last_bin':41,
                'frequency_workers':8,'solve_workers':1}.items():
        require(d.get(k)==v,f"{label} {k} differs from reviewed workload")
    require(d.get('backend')=='dspec' and d.get('phase_c_direct_adjoint') is True and
            d.get('phase_c_combined_reconstruction') is True and
            d.get('kernel_adjoint_weight_source')=="this mesh's own waveform-derived adjoint source",
            f"{label} did not use the reviewed direct/common own-weight route")

def verify_skip(path,expected_skipped):
    data=rows(Path(path)/'phase-c-zero-bin-skip.csv')
    require({int(x['bin']) for x in data}==EXPECTED_SKIP,"zero-bin audit does not contain exactly bins 7,8,41")
    for x in data:
        value=complex(float(x['adjoint_weight_force_product_re']),float(x['adjoint_weight_force_product_im']))
        require(value==0j and int(x['exact_zero'])==1,"skip audit includes a nonzero measurement product")
        require(int(x['skipped'])==int(expected_skipped),"skip flag differs from requested on/off control")

def finite_binary(path):
    raw=Path(path).read_bytes();require(len(raw)%8==0,f"binary float file has invalid byte length: {path}")
    values=struct.unpack('<'+('d'*(len(raw)//8)),raw)
    require(all(math.isfinite(x) for x in values),f"nonfinite binary sample in {path}")

def equal_files(a,b,names):
    comparisons={}
    for name in names:
        pa=Path(a)/name;pb=Path(b)/name
        require(pa.is_file() and pb.is_file(),f"missing parity output {name}")
        if name.endswith('.bin'):finite_binary(pa);finite_binary(pb)
        ha,hb=digest(pa),digest(pb)
        require(ha==hb,f"native/profile or skip on/off bytes differ for {name}: {ha} != {hb}")
        comparisons[name]=ha
    return comparisons

def check_reconstruction(path):
    data=rows(Path(path)/'phase-c-reconstruction-equivalence.csv')
    require(len(data)==70,"reconstruction equivalence inventory must be 35 bins x 2 sections")
    keys=[(int(r['bin']),int(r['section'])) for r in data]
    require(len(set(keys))==70 and set(keys)=={(b,s) for b in range(7,42) for s in (0,1)},
            "reconstruction equivalence must contain each bin/section pair exactly once")
    for r in data:
        for k in ('forward_field_rel_l2','adjoint_field_rel_l2','complex_product_rel_l2'):
            x=float(r[k]);require(math.isfinite(x) and x<=1e-11,f"reconstruction {k} failed at bin {r['bin']} section {r['section']}: {x:.17g}")
        if int(r['bin']) in (9,25,40):
            for k in ('forward_oracle_norm','adjoint_oracle_norm','product_oracle_norm'):
                x=float(r[k]);require(math.isfinite(x) and x>0.0,
                    f"reconstruction {k} is not positive at bin {r['bin']} section {r['section']}")

def profile_calls(path):
    p=read_json(Path(path)/'profile.json')
    calls={};phase_calls={}
    for row in p.get('rows',[]):
        name=row.get('scope_name','').removeprefix('frequency_worker.')
        calls[name]=calls.get(name,0)+int(row.get('calls',0))
        key=(row.get('phase',''),name)
        phase_calls[key]=phase_calls.get(key,0)+int(row.get('calls',0))
    counters={}
    for row in p.get('counters',[]):
        name=row['name'].removeprefix('frequency_worker.')
        counters[name]=counters.get(name,0)+int(row['value'])
    return calls,counters,phase_calls

def check_skip_profile(path):
    d=read_json(Path(path)/'result.json')
    require(float(d.get('cache_bytes',-1))==0.0,"skip-on direct route retained response-cache bytes")
    calls,counters,phase_calls=profile_calls(path)
    for name,value in calls.items():
        forbidden=('minimal_unit_basis_HH_solve','receiver_basis_and_QR_reconstruction',
                   'receiver_basis_mix_and_unit_adjoint_projection',
                   'receiver_basis_assembly_and_qr','batch_cache','batch_response_cache')
        if any(token in name for token in forbidden):
            require(value==0,f"skip-on path performed forbidden unit-basis/QR/cache scope {name}: {value}")
    require(counters.get('solver_unit_basis_rhs_columns',0)==0,
            "skip-on path performed unit-basis adjoint work")
    require(calls.get('initial_frequency_workspace_forward_and_unit_adjoint_capture',0)==35 and
            sum(v for (phase,name),v in phase_calls.items()
                if phase=='initial_frequency_workspace_and_capture' and name=='frequency_workspace_prepare_and_modal_capture')==35,
            "skip-on path did not retain all 35 initial forward workspaces")
    require(sum(v for (phase,name),v in phase_calls.items()
                if phase=='kernel_frequency_parallel_pass_and_ordered_reduction_including_diagnostics' and
                name=='frequency_workspace_prepare_and_modal_capture')==32,
            "skip-on path did not execute exactly 32 nonzero kernel workspaces")
    require(counters.get('phase_c_exact_zero_kernel_bin_skips',0)==3,
            "skip-on path did not skip exactly three exact-zero bins")

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--validation',type=Path,required=True);ap.add_argument('--skip-off',type=Path,required=True);ap.add_argument('--skip-on',type=Path,required=True);ap.add_argument('--native',type=Path,required=True);ap.add_argument('--profile',type=Path,required=True);ap.add_argument('--a4-components',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    v=read_json(a.validation/'result.json');off=read_json(a.skip_off/'result.json');on=read_json(a.skip_on/'result.json');native=read_json(a.native/'result.json');profile=read_json(a.profile/'result.json')
    for d,n in ((v,'validation'),(off,'skip-off'),(on,'skip-on'),(native,'native'),(profile,'profile')):verify_case(d,n)
    require(v.get('phase_c_validate_direct_adjoint') is True and v.get('phase_c_validate_reconstruction') is True and v.get('phase_c_zero_bin_skips') is False,"full L200 equivalence control flags mismatch")
    check_reconstruction(a.validation)
    require(float(v.get('phase_c_max_reconstruction_validation_relative_l2',math.inf))<=1e-11,"L200 reconstruction scalar gate failed")
    validation_profile=read_json(a.validation/'profile.json')
    direct_counts={x['name'].removeprefix('frequency_worker.'):int(x['value'])
                   for x in validation_profile.get('counters',[])
                   if x['name'].removeprefix('frequency_worker.').startswith('phase_c_direct_')}
    direct_rhs_error=direct_counts.get('phase_c_direct_rhs_equivalence_scaled_1e12',0)
    require(direct_rhs_error==0 and 'phase_c_direct_solution_equivalence_scaled_1e12' in direct_counts and
            direct_counts.get('phase_c_direct_adjoint_rhs_columns',0)>0,
            "L200 direct-adjoint RHS/solution equivalence counters are absent or direct solves were not exercised")
    check_skip_profile(a.skip_on)
    require(off.get('phase_c_zero_bin_skips') is False and on.get('phase_c_zero_bin_skips') is True,"skip on/off metadata mismatch")
    verify_skip(a.skip_off,False);verify_skip(a.skip_on,True)
    parity_names=['trace.csv','spectrum.csv','components.csv','nodes.csv','frequency-grid.csv','mesh-elements.csv','along-angles.csv','across-angles.csv','along.bin','across.bin']
    skip_hashes=equal_files(a.skip_off,a.skip_on,parity_names)
    validation_hashes=equal_files(a.validation,a.skip_off,parity_names+['phase-c-zero-bin-skip.csv'])
    a4=rows(a.a4_components);candidate_components=rows(a.profile/'components.csv')
    require(len(a4)==len(candidate_components)==35,"Phase A4 forward spectrum row count changed")
    for x,y in zip(candidate_components,a4):require(int(x['bin'])==int(y['bin']),"Phase A4 forward bin grid changed")
    a4_component_rel={key:math.sqrt(sum((float(x[key])-float(y[key]))**2 for x,y in zip(candidate_components,a4))/
                                    max(sum(float(y[key])**2 for y in a4),1e-300))
                      for key in ('up_re','up_im','north_re','north_im','east_re','east_im')}
    require(all(math.isfinite(x) and x<=1e-11 for x in a4_component_rel.values()),
            f"Phase C versus accepted A4 forward component mismatch: {a4_component_rel}")
    common_names=[f'common-fields-bin{bin}.csv' for bin in (9,25,40)]+[f'kernel-workspace-receiver-bin{bin}.csv' for bin in (9,25,40)]
    native_profile_hashes=equal_files(a.native,a.profile,parity_names+['phase-c-zero-bin-skip.csv'])
    common=[]
    for bin in (9,25,40):
        file=f'common-fields-bin{bin}.csv'
        require((a.validation/file).is_file(),f"missing common reference packet comparison {file}")
        common.extend(common_field_metrics(a.validation/file))
    for label,d in (('native',native),('profile',profile)):
        require(abs(float(d.get('time_shift_adjoint_identity',0))-1)<=1e-10,f"{label} time-shift identity failed")
        require(float(d.get('gravity_max_abs_nd_error',math.inf))<=1e-10,f"{label} gravity check failed")
    report={'status':'pass','case':{'degree':200,'nq':10,'refinement':0.025,'bins':35,'frequency_workers':8,'inner_workers':1},
      'l200_direct_adjoint_and_common_reconstruction_gate':True,'reconstruction_rows':len(rows(a.validation/'phase-c-reconstruction-equivalence.csv')),
      'direct_adjoint_quantized_equivalence_counters':{
          'rhs_scaled_1e12':direct_rhs_error,
          'solution_scaled_1e12':direct_counts['phase_c_direct_solution_equivalence_scaled_1e12'],
          'rhs_columns':direct_counts['phase_c_direct_adjoint_rhs_columns']},
      'validation_vs_timed_skip_off_exact_output_hashes':validation_hashes,
      'skip_off_on_exact_output_hashes':skip_hashes,'native_profile_exact_output_hashes':native_profile_hashes,
      'common_reference_unmasked_group_metrics':common,'zero_skip_bins':sorted(EXPECTED_SKIP),
      'a4_forward_component_relative_l2':a4_component_rel,
      'time_shift_identity':{'native':native['time_shift_adjoint_identity'],'profile':profile['time_shift_adjoint_identity']},
      'gravity_max_abs_nd_error':{'native':native['gravity_max_abs_nd_error'],'profile':profile['gravity_max_abs_nd_error']}}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps(report,sort_keys=True))
if __name__=='__main__':main()
