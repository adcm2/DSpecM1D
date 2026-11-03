#!/usr/bin/env python3
"""Bounded C2 one-pass, exact-skip, and complete-pipeline parity checker."""
from __future__ import annotations
import argparse,csv,hashlib,json,math,struct
from pathlib import Path

BINS={7,8,41}
SCIENCE_FILES=('trace.csv','spectrum.csv','components.csv','nodes.csv','frequency-grid.csv','mesh-elements.csv','along-angles.csv','across-angles.csv','along.bin','across.bin')

def need(ok,msg):
    if not ok: raise ValueError(msg)
def read_json(path):
    d=json.loads(Path(path).read_text())
    def walk(x):
        if isinstance(x,bool) or x is None or isinstance(x,str): return
        if isinstance(x,(int,float)): need(math.isfinite(float(x)),f'nonfinite JSON number in {path}')
        elif isinstance(x,list):
            for v in x: walk(v)
        elif isinstance(x,dict):
            for v in x.values(): walk(v)
    walk(d);return d
def rows(path):
    with Path(path).open(newline='') as f: return list(csv.DictReader(f))
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def rel_l2(a,b):
    need(len(a)==len(b),'array lengths differ')
    err=norm=0.0
    for x,y in zip(a,b):
        need(math.isfinite(x) and math.isfinite(y),'nonfinite array sample')
        err+=(x-y)*(x-y);norm+=y*y
    return math.sqrt(err/norm) if norm else (0.0 if err==0 else math.inf)
def doubles(path):
    raw=Path(path).read_bytes();need(len(raw)%8==0,f'invalid binary byte length: {path}')
    vals=struct.unpack('<'+'d'*(len(raw)//8),raw)
    need(all(math.isfinite(x) for x in vals),f'nonfinite binary value: {path}')
    return vals
def verify_case(d,label):
    for k,v in {'degree':200,'nq':10,'refinement':0.025,'frequency_count':35,'first_bin':7,'last_bin':41,'frequency_workers':8,'solve_workers':1}.items():
        need(d.get(k)==v,f'{label} {k} mismatch')
    need(d.get('backend')=='dspec' and d.get('phase_c2_preweighted_complex_contraction') is True and d.get('phase_c_direct_adjoint') is True and d.get('phase_c_combined_reconstruction') is True,f'{label} did not use C2 combined direct/common route')
    need(d.get('phase_c2_product_storage_bytes')==25268320,f'{label} C2 product payload is not 25,268,320 bytes')
    need(d.get('kernel_adjoint_weight_source')=="this mesh's own waveform-derived adjoint source",f'{label} does not use its own measured P')
    need(abs(float(d.get('time_shift_adjoint_identity',0))-1)<=1e-10,f'{label} time-shift identity failed')
    need(float(d.get('gravity_max_abs_nd_error',math.inf))<=1e-10,f'{label} gravity check failed')
def profile(path):
    p=read_json(Path(path)/'profile.json');calls={};counters={};phases={}
    for r in p.get('rows',[]):
        n=r.get('scope_name','').removeprefix('frequency_worker.')
        calls[n]=calls.get(n,0)+int(r.get('calls',0));key=(r.get('phase',''),n);phases[key]=phases.get(key,0)+int(r.get('calls',0))
    for r in p.get('counters',[]):
        n=r['name'].removeprefix('frequency_worker.');counters[n]=counters.get(n,0)+int(r['value'])
    return calls,counters,phases
def check_zero(path,skip):
    data=rows(Path(path)/'phase-c2-zero-bin-support.csv')
    need({int(x['bin']) for x in data}==BINS,'C2 exact-zero support is not bins7/8/41')
    for x in data:
        need(x['taper_exact_zero']=='1' and x['measured_P_exact_zero']=='1' and x['forward_solve_retained']=='1','C2 static zero proof fields failed')
        need(int(x['unit_adjoint_skipped'])==int(skip),'C2 skip on/off marker mismatch')
def check_profile(path,skip):
    calls,counters,phases=profile(path)
    need(counters.get('phase_c2_forward_spectrum_frequency_calls',0)==35,'C2 did not retain 35 forward frequency calls')
    expected=32 if skip else 35
    need(counters.get('phase_c2_unit_adjoint_frequency_calls',0)==expected,'C2 active unit-adjoint frequency calls mismatch')
    need(counters.get('phase_c2_exact_zero_unit_adjoint_skips',0)==(3 if skip else 0),'C2 exact zero skip counter mismatch')
    if skip:
        need(float(read_json(Path(path)/'result.json').get('cache_bytes',-1))==0,'C2 retained response-cache bytes')
        for n,v in calls.items():
            if any(x in n for x in ('minimal_unit_basis_HH_solve','receiver_basis_and_QR_reconstruction','batch_cache','batch_response_cache')):
                need(v==0,f'C2 timing route performed forbidden basis/cache work: {n}={v}')
def check_common(path):
    data=rows(Path(path)/'phase-c2-common-reconstruction-equivalence.csv')
    need(len(data)==70,'C2 reconstruction proof inventory must be 35 bins x 2 sections')
    keys={(int(x['bin']),int(x['section'])) for x in data}
    need(keys=={(b,s) for b in range(7,42) for s in (0,1)},'C2 reconstruction proof key inventory differs')
    for x in data:
        b=int(x['bin'])
        for k in ('forward_field_rel_l2','adjoint_field_rel_l2','complex_product_rel_l2'):
            v=float(x[k]);need(math.isfinite(v) and v<=1e-11,f'C2 common reconstruction {k} failed at bin {b}, section {x["section"]}: {v:.17g}')
        if b in (9,25,40):
            for k in ('forward_oracle_norm2','adjoint_oracle_norm2','product_oracle_norm2'):
                v=float(x[k]);need(math.isfinite(v) and v>0,f'C2 common oracle {k} missing at bin {b}')
    return {k:max(float(x[k]) for x in data) for k in ('forward_field_rel_l2','adjoint_field_rel_l2','complex_product_rel_l2')}
def check_linearity(path):
    data=rows(Path(path)/'phase-c2-weight-linearity.csv')
    need(len(data)==6,'C2 unit/weighted linearity audit must contain 3 bins x 2 sections')
    need({(int(x['bin']),int(x['section'])) for x in data}=={(b,s) for b in (9,25,40) for s in (0,1)},'C2 linearity inventory differs')
    for r in data:
        for k in ('adjoint_input_vector_scale_rel','assembled_direct_rhs_scale_rel','forward_radial_mode_scale_rel','adjoint_radial_mode_scale_rel','forward_field_scale_rel','adjoint_field_scale_rel','stored_unit_product_rel','complex_product_scale_rel'):
            v=float(r[k]);need(math.isfinite(v) and v<=1e-11,f'C2 linearity {k} failed at bin {r["bin"]}, section {r["section"]}: {v:.17g}')
def check_direct_profile(path):
    _,counters,_=profile(path)
    rhs=counters.get('phase_c_direct_rhs_equivalence_scaled_1e12',0)
    need(rhs==0,'C2 direct RHS equivalence counter is nonzero')
    need('phase_c_direct_solution_equivalence_scaled_1e12' in counters,'C2 direct solution equivalence counter is absent')
    need(counters.get('phase_c_direct_adjoint_rhs_columns',0)>0,'C2 direct-adjoint validation did not exercise RHS columns')
def check_common_packets(path):
    summaries=[]
    for b in (9,25,40):
        data=rows(Path(path)/f'common-fields-bin{b}.csv');need(len(data)==293*(135+19),f'common packet inventory differs at bin {b}')
        groups={}
        for r in data:
            key=(int(r['section']),r['location_class'],int(r['fluid']))
            g=groups.setdefault(key,[0.0]*4+[0])
            for role in ('forward','adjoint'):
                err=norm=0.0
                for c in range(6):
                    x=complex(float(r[f'{role}_candidate{c}_re']),float(r[f'{role}_candidate{c}_im']))
                    y=complex(float(r[f'{role}_ref{c}_re']),float(r[f'{role}_ref{c}_im']))
                    need(math.isfinite(x.real) and math.isfinite(x.imag) and math.isfinite(y.real) and math.isfinite(y.imag),'nonfinite C2 packet field')
                    err+=abs(x-y)**2;norm+=abs(y)**2
                off=0 if role=='forward' else 2;g[off]+=err;g[off+1]+=norm
            g[4]+=1
        for key,g in sorted(groups.items()):
            fr=math.sqrt(g[0]/g[1]) if g[1] else (0.0 if g[0]==0 else math.inf)
            ar=math.sqrt(g[2]/g[3]) if g[3] else (0.0 if g[2]==0 else math.inf)
            summaries.append({'bin':b,'section':key[0],'location_class':key[1],'fluid':key[2],'count':g[4],'forward_relative_l2':fr,'adjoint_relative_l2':ar})
    return summaries

def compare_csv_channels(candidate, reference, label):
    a,b=rows(candidate),rows(reference)
    need(len(a)==len(b),f'{label} row count differs')
    import csv as _csv
    with Path(candidate).open(newline='') as f: ah=_csv.DictReader(f).fieldnames
    with Path(reference).open(newline='') as f: bh=_csv.DictReader(f).fieldnames
    need(ah==bh,f'{label} column inventory differs')
    exact={'bin','time_s','frequency_hz','frequency_bin','sample','window'}
    metrics={}
    for key in ah:
        av=[r[key] for r in a];bv=[r[key] for r in b]
        if key in exact:
            need(av==bv,f'{label} exact grid column {key} differs')
            continue
        try:
            x=[float(v) for v in av];y=[float(v) for v in bv]
        except (TypeError,ValueError):
            need(av==bv,f'{label} nonnumeric column {key} differs')
            continue
        e=rel_l2(x,y);need(e<=1e-11,f'{label} channel {key} differs: {e:.17g}')
        metrics[key]=e
    return metrics

def scalar_rel(a,b):
    a=float(a);b=float(b)
    need(math.isfinite(a) and math.isfinite(b),'nonfinite measurement-rate norm')
    return abs(a-b)/max(abs(b),1e-300)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--validation',type=Path,required=True);ap.add_argument('--skip-off',type=Path,required=True);ap.add_argument('--skip-on',type=Path,required=True);ap.add_argument('--native',type=Path,required=True);ap.add_argument('--profile',type=Path,required=True);ap.add_argument('--phase-c-oracle',type=Path,required=True);ap.add_argument('--f8a-reference',type=Path,required=True);ap.add_argument('--a4-components',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    dirs=(('validation',a.validation),('skip-off',a.skip_off),('skip-on',a.skip_on),('native',a.native),('profile',a.profile))
    results={label:read_json(path/'result.json') for label,path in dirs}
    for label,d in results.items():verify_case(d,label)
    need(results['validation'].get('phase_c_validate_direct_adjoint') is True and results['validation'].get('phase_c_validate_reconstruction') is True,'C2 validation flags are not enabled')
    need(results['skip-off'].get('phase_c_zero_bin_skips') is False and results['skip-on'].get('phase_c_zero_bin_skips') is True,'C2 skip control metadata mismatch')
    check_zero(a.skip_off,False);check_zero(a.skip_on,True);check_profile(a.skip_off,False);check_profile(a.skip_on,True);common_reconstruction_max=check_common(a.validation);check_linearity(a.validation);check_direct_profile(a.validation)
    packet_metrics=check_common_packets(a.validation)
    pc=read_json(a.phase_c_oracle/'result.json');f8=read_json(a.f8a_reference/'result.json')
    for label,d in (('two-pass Phase C',pc),('accepted F8a',f8)):
        for k in ('degree','nq','refinement','frequency_count','first_bin','last_bin','native_side_count','along_angles','across_angles'):
            need(d.get(k)==results['native'].get(k),f'{label} {k} grid mismatch')
    comparisons={}
    # All C2 jobs must emit byte-identical science products, including their own
    # retained-grid and mesh inventories. Historical F8a predates those two
    # inventory files, so compare its saved radial/angle geometry instead.
    for name in SCIENCE_FILES:
        for _,path in dirs:
            need((path/name).is_file(),f'missing C2 scientific output {path/name}')
        hashes={label:sha(path/name) for label,path in dirs}
        need(len(set(hashes.values()))==1,f'C2 skip/native/profile science outputs differ for {name}')
        comparisons[name]={'candidate_sha256':hashes['native'],'profile_sha256':hashes['profile']}
    historical_geometry={}
    for name in ('nodes.csv','along-angles.csv','across-angles.csv'):
        for path in (a.phase_c_oracle,a.f8a_reference):
            need((path/name).is_file(),f'missing historical geometry output {path/name}')
        hashes={label:sha(path/name) for label,path in (('C2',a.native),('two_pass_C',a.phase_c_oracle),('accepted_F8a',a.f8a_reference))}
        need(len(set(hashes.values()))==1,f'historical geometry {name} differs: {hashes}')
        historical_geometry[name]=hashes
    historical_channels={}
    for name in ('trace.csv','spectrum.csv','components.csv'):
        need((a.phase_c_oracle/name).is_file(),f'missing two-pass waveform output {a.phase_c_oracle/name}')
        historical_channels[name]={'two_pass_C':compare_csv_channels(a.native/name,a.phase_c_oracle/name,'C2 vs two-pass '+name)}
    need((a.f8a_reference/'components.csv').is_file(),'missing accepted F8a forward components')
    historical_channels['components.csv']['accepted_F8a_forward']=compare_csv_channels(a.native/'components.csv',a.f8a_reference/'components.csv','C2 vs accepted F8a components')
    rate_norm_errors={label:scalar_rel(results['native']['measurement_rate_norm'],ref['measurement_rate_norm']) for label,ref in (('two_pass_C',pc),('accepted_F8a',f8))}
    need(all(v<=1e-11 for v in rate_norm_errors.values()),f'C2 measurement_rate_norm parity failed: {rate_norm_errors}')
    # C2 weighted contraction may change summation roundoff; preserve the established common-route 1e-11 check.
    kernel_errors={}
    for name in ('along.bin','across.bin'):
        values=doubles(a.native/name);old=doubles(a.phase_c_oracle/name);accepted=doubles(a.f8a_reference/name)
        e_old=rel_l2(values,old);e_accepted=rel_l2(values,accepted)
        need(e_old<=1e-11 and e_accepted<=1e-11,f'C2 {name} differs from preserved complete-pipeline reference: C={e_old:.17g}, F8a={e_accepted:.17g}')
        kernel_errors[name]={'two_pass_phase_c_relative_l2':e_old,'accepted_f8a_relative_l2':e_accepted}
    a4=rows(a.a4_components);comp=rows(a.native/'components.csv');need(len(a4)==len(comp)==35,'A4/C2 component inventory differs')
    component_errors={}
    for key in ('up_re','up_im','north_re','north_im','east_re','east_im'):
        err=norm=0.0
        for x,y in zip(comp,a4):
            need(int(x['bin'])==int(y['bin']),'A4/C2 component bin mismatch')
            dx=float(x[key]);dy=float(y[key]);need(math.isfinite(dx) and math.isfinite(dy),'nonfinite waveform component')
            err+=(dx-dy)**2;norm+=dy*dy
        v=math.sqrt(err/norm) if norm else (0.0 if err==0 else math.inf);need(v<=1e-11,f'A4 waveform component {key} differs: {v:.17g}');component_errors[key]=v
    timing={k:results['native'].get(k) for k in ('combined_forward_unit_adjoint_reconstruction_contraction_s','measurement_and_actual_adjoint_spectrum_s','weighted_complex_product_reduction_s','elapsed_s','cache_bytes','phase_c2_product_storage_bytes')}
    fluid_diagnostics=[{k:float(r[k]) for k in ('bin','section','fluid_forward_rel','fluid_adjoint_rel','fluid_complex_product_rel')} for r in rows(a.validation/'phase-c2-weight-linearity.csv')]
    report={'status':'pass','case':{'degree':200,'nq':10,'refinement':0.025,'bins':35,'frequency_workers':8,'inner_workers':1},'exact_zero_bins':sorted(BINS),'skip_on_off_work_proof':{'forward_calls':35,'unit_adjoint_calls':{'skip_off':35,'skip_on':32},'unit_adjoint_skips':{'skip_off':0,'skip_on':3}},'c2_common_tensor_reconstruction_rows':70,'common_reconstruction_relative_l2_maxima_from_csv':common_reconstruction_max,'legacy_result_scalar_used':False,'legacy_result_scalar_note':'result.json phase_c_max_reconstruction_validation_relative_l2 is a stale placeholder; actual maxima are derived from phase-c2-common-reconstruction-equivalence.csv','weighted_unweighted_linearity_rows':6,'corrected_phase_a_packet_side_metrics':packet_metrics,'fluid_linearity_diagnostics':fluid_diagnostics,'a4_forward_component_relative_l2':component_errors,'full35_kernel_relative_l2':kernel_errors,'scientific_file_hashes':comparisons,'historical_geometry_hashes':historical_geometry,'historical_trace_spectrum_component_relative_l2':historical_channels,'measurement_rate_norm_relative_error':rate_norm_errors,'native_timing_fields':timing}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps(report,sort_keys=True))
if __name__=='__main__':main()
