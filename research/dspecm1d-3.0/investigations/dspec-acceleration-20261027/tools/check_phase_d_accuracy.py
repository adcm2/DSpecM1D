#!/usr/bin/env python3
"""Bounded Phase D fixed-point/full35 accuracy and batching proof checker."""
from __future__ import annotations
import argparse, csv, json, math, pathlib, struct

def rows(path):
    with pathlib.Path(path).open(newline="") as f:
        return list(csv.DictReader(f))

def read_json(path):
    return json.loads(pathlib.Path(path).read_text())

def rel_l2(err, norm):
    return math.sqrt(err / norm) if norm else (0.0 if err == 0 else math.inf)

def finite(x, label):
    if not math.isfinite(x):
        raise ValueError(f"non-finite {label}: {x}")

def scalar_channels(candidate, reference, cols, exact=()):
    a, b = rows(pathlib.Path(candidate)), rows(pathlib.Path(reference))
    if len(a) != len(b): raise ValueError("row count mismatch")
    metrics = {}
    for col in exact:
        if [r[col] for r in a] != [r[col] for r in b]: raise ValueError(f"grid mismatch: {col}")
    for col in cols:
        err = norm = 0.0
        for x, y in zip(a, b):
            xv, yv = float(x[col]), float(y[col]); finite(xv, col); finite(yv, col)
            err += (xv-yv)**2; norm += yv*yv
        metrics[col] = rel_l2(err, norm)
    return metrics

def inventory(path):
    result=[]
    for r in rows(path):
        result.append((int(r['section']),int(r['reference_element']),int(r['reference_node']),
          int(r['section_point']),int(r['layer']),float(r['radius_m']),int(r['fluid']),r['location_class']))
    return result

def load_p(path):
    out={}
    for r in rows(path):
        b=int(r['bin']); p=complex(float(r['P_re']),float(r['P_im'])); finite(p.real,f"P[{b}].re"); finite(p.imag,f"P[{b}].im")
        if b in out: raise ValueError(f"duplicate P bin {b}")
        if int(r['exact_zero']) != (1 if b in (7,8,41) else 0): raise ValueError(f"P exact-zero support differs at bin {b}")
        out[b]=p
    return out

def compare_common_fields(candidate, reference, inv):
    cp, rp = load_p(candidate/'phase-d-own-P.csv'), load_p(reference/'phase-d-own-P.csv')
    if set(cp)!=set(range(7,42)) or set(rp)!=set(range(7,42)): raise ValueError("own-P inventory must be all35 bins")
    groups={}
    for b in range(7,42):
        cf=(candidate/f"phase-d-common-unit-fields-bin{b}.bin").read_bytes()
        rf=(reference/f"phase-d-common-unit-fields-bin{b}.bin").read_bytes()
        need=len(inv)*12*16
        if len(cf)!=need or len(rf)!=need: raise ValueError(f"common unit-field byte size differs at bin {b}")
        pc,pr=cp[b],rp[b]
        for i,key in enumerate(inv):
            section,_,_,_,_,_,fluid,loc=key
            g=groups.setdefault((b,section,loc,fluid),[0.0]*8)
            off=i*192
            cv=struct.unpack_from('<24d',cf,off);rv=struct.unpack_from('<24d',rf,off)
            # Six forward components followed by six unit-adjoint components.
            for c in range(6):
                cfr=complex(cv[2*c],cv[2*c+1]);rfr=complex(rv[2*c],rv[2*c+1])
                car=complex(cv[12+2*c],cv[13+2*c])*pc.conjugate()
                rar=complex(rv[12+2*c],rv[13+2*c])*pr.conjugate()
                for z,label in ((cfr,'candidate forward'),(rfr,'reference forward'),(car,'candidate weighted adjoint'),(rar,'reference weighted adjoint')):
                    finite(z.real,label);finite(z.imag,label)
                g[0]+=abs(cfr-rfr)**2;g[1]+=abs(rfr)**2
                g[2]+=abs(car-rar)**2;g[3]+=abs(rar)**2
                if fluid:
                    g[4]+=abs(cfr-rfr)**2;g[5]+=abs(rfr)**2
                    g[6]+=abs(car-rar)**2;g[7]+=abs(rar)**2
    summary=[]
    for (b,section,loc,fluid),g in sorted(groups.items()):
        f=rel_l2(g[0],g[1]);a=rel_l2(g[2],g[3])
        item={'bin':b,'section':section,'location_class':loc,'fluid':fluid,
          'forward_relative_l2':f,'forward_reference_norm2':g[1],
          'forward_absolute_l2':math.sqrt(g[0]),
          'adjoint_relative_l2':a,'adjoint_reference_norm2':g[3],
          'adjoint_absolute_l2':math.sqrt(g[2])}
        if fluid:
            item['fluid_forward_diagnostic_l2']=rel_l2(g[4],g[5]);item['fluid_adjoint_diagnostic_l2']=rel_l2(g[6],g[7])
        summary.append(item)
    return summary

def compare_common_kernel(candidate, reference, inv):
    a,b=rows(candidate/'phase-d-common-full35-kernel.csv'),rows(reference/'phase-d-common-full35-kernel.csv')
    if len(a)!=len(inv) or len(b)!=len(inv): raise ValueError("common full35 kernel row count mismatch")
    groups={}
    for i,(x,y,key) in enumerate(zip(a,b,inv)):
        xk=(int(x['section']),int(x['reference_element']),int(x['reference_node']),int(x['section_point']),int(x['layer']),float(x['radius_m']),int(x['fluid']),x['location_class'])
        yk=(int(y['section']),int(y['reference_element']),int(y['reference_node']),int(y['section_point']),int(y['layer']),float(y['radius_m']),int(y['fluid']),y['location_class'])
        if xk!=key or yk!=key: raise ValueError(f"common kernel keyed sample mismatch at row {i}")
        fluid=key[6];g=groups.setdefault((key[0],key[7],fluid),[0.0]*2)
        xv,yv=float(x['kernel_vs']),float(y['kernel_vs']);finite(xv,'candidate common kernel');finite(yv,'reference common kernel')
        g[0]+=(xv-yv)**2;g[1]+=yv*yv
    return [{'section':s,'location_class':c,'fluid':f,'relative_l2':rel_l2(v[0],v[1]),
      'absolute_l2':math.sqrt(v[0]),'reference_norm2':v[1]} for (s,c,f),v in sorted(groups.items())]

def check_zero(path):
    data=rows(path/'phase-c2-zero-bin-support.csv')
    if [int(r['bin']) for r in data] != [7,8,41]: raise ValueError("exact zero support must be bins7,8,41")
    for r in data:
        if any(int(r[k])!=1 for k in ('taper_exact_zero','measured_P_exact_zero','forward_solve_retained')) or int(r['unit_adjoint_skipped'])!=0:
            raise ValueError("validation exact-zero support provenance or skip-off setting mismatch")

def proof(path, filename, label):
    data=rows(path/filename)
    expected={(b,s) for b in range(7,42) for s in (0,1)}
    if len(data)!=70 or {(int(x['bin']),int(x['section'])) for x in data}!=expected:
        raise ValueError(f"{label} proof inventory must be exact35x2 bin/section rows")
    maxima={k:0.0 for k in ('forward_field_rel_l2','adjoint_field_rel_l2','complex_product_rel_l2')}
    for x in data:
        for k in maxima:
            v=float(x[k]);finite(v,f"{label}.{k}");maxima[k]=max(maxima[k],v)
            if v>1e-11: raise ValueError(f"{label} {k} exceeds1e-11 at bin{x['bin']} section{x['section']}: {v:.17g}")
        if int(x['bin']) in (9,25,40):
            normkeys=('forward_oracle_norm2','adjoint_oracle_norm2','product_oracle_norm2') if 'reconstruction-equivalence' in filename else ('forward_scalar_norm2','adjoint_scalar_norm2','product_scalar_norm2')
            if any(not math.isfinite(float(x[k])) or float(x[k])<=0 for k in normkeys):
                raise ValueError(f"{label} positive active oracle norms missing at bin{x['bin']}")
    return maxima

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--reference',type=pathlib.Path,required=True);ap.add_argument('--candidate',type=pathlib.Path,required=True)
    ap.add_argument('--accepted',type=pathlib.Path,required=True);ap.add_argument('--accuracy',type=float,default=1e-3)
    ap.add_argument('--out',type=pathlib.Path,required=True);a=ap.parse_args()
    r,c,f=(read_json(p/'result.json') for p in (a.reference,a.candidate,a.accepted))
    grids={}
    for label,d,path in (('reference',r,a.reference),('candidate',c,a.candidate)):
        if d.get('backend')!='dspec' or d.get('phase_d_batched_angular') is not True or d.get('phase_d_common_validation') is not True:
            raise ValueError(f"{label} must be batched DSpec common-field validation")
        if d.get('phase_c_validate_reconstruction') is not True or d.get('phase_c_zero_bin_skips') is not False:
            raise ValueError(f"{label} common/tensor oracle or zero-skip settings invalid")
        if d.get('frequency_count')!=35 or d.get('frequency_workers')!=8 or d.get('solve_workers')!=1:
            raise ValueError(f"{label} frequency/workers inventory mismatch")
        if d.get('first_bin')!=7 or d.get('last_bin')!=41: raise ValueError(f"{label} grid endpoints differ")
        if not math.isfinite(float(d['gravity_max_abs_nd_error'])) or float(d['gravity_max_abs_nd_error'])>1e-10: raise ValueError(f"{label} gravity check failed")
        if not math.isfinite(float(d['time_shift_adjoint_identity'])) or abs(float(d['time_shift_adjoint_identity'])-1)>1e-10: raise ValueError(f"{label} time-shift identity failed")
        if int(d['phase_d_actual_temporary_peak_bytes'])>32*1024*1024: raise ValueError(f"{label} angular bound exceeded")
        grid=rows(path/'frequency-grid.csv');grids[label]=grid
        if [int(x['bin']) for x in grid]!=list(range(7,42)): raise ValueError(f"{label} frequency-grid inventory differs")
        for x in grid:
            finite(float(x['z_real']),f"{label} z.real");finite(float(x['z_imag']),f"{label} z.imag")
    if [(x['bin'],x['z_real'],x['z_imag']) for x in grids['reference']] != [(x['bin'],x['z_real'],x['z_imag']) for x in grids['candidate']]:
        raise ValueError("NQ10/NQ6 frequency grid differs")
    if r['nq']!=10 or abs(float(r['refinement'])-.025)>1e-15 or r['relative_error']!=1e-5: raise ValueError("Phase D NQ10 reference contract differs")
    if f.get('degree')!=200 or f.get('nq')!=10 or f.get('frequency_count')!=35 or abs(float(f.get('refinement',-1))-.025)>1e-15:
        raise ValueError("accepted C2/F8a reference contract differs")
    if c['nq']!=6 or c['degree']!=200 or abs(float(c['relative_error'])-1e-5)>1e-15: raise ValueError("Phase D NQ6 candidate contract differs")
    if abs(float(c['phase_d_requested_relative_accuracy'])-a.accuracy)>1e-15: raise ValueError("requested accuracy metadata differs")
    check_zero(a.reference);check_zero(a.candidate)
    invr=inventory(a.reference/'phase-d-common-point-inventory.csv');invc=inventory(a.candidate/'phase-d-common-point-inventory.csv')
    if invr!=invc: raise ValueError("reference/candidate common physical-side inventories differ")
    identity=[(x[0],x[1],x[2],x[3]) for x in invr]
    if len(invr)!=45122 or len(set(identity))!=45122: raise ValueError("common physical-side inventory must contain 45,122 unique keyed locations")
    for name in ('nodes.csv','along-angles.csv','across-angles.csv'):
        geometry=rows(a.reference/name)
        scalar_channels(a.reference/name,a.accepted/name,(),exact=tuple(geometry[0].keys()))
    f8_components=scalar_channels(a.reference/'components.csv',a.accepted/'components.csv',
        ('up_re','up_im','north_re','north_im','east_re','east_im'),exact=('bin',))
    f8_trace=scalar_channels(a.reference/'trace.csv',a.accepted/'trace.csv',('velocity_m_s','rate_m_s2','adjoint_source'),exact=('time_s','window'))
    f8_spectrum=scalar_channels(a.reference/'spectrum.csv',a.accepted/'spectrum.csv',('real_m','imag_m'),exact=('bin','frequency_hz'))
    f8_measurement_norm=abs(float(r['measurement_rate_norm'])-float(f['measurement_rate_norm']))/max(abs(float(f['measurement_rate_norm'])),1e-300)
    f8_kernels={}
    for name in ('along.bin','across.bin'):
        x=(a.reference/name).read_bytes();y=(a.accepted/name).read_bytes()
        if len(x)!=len(y) or len(x)%8: raise ValueError(f"accepted kernel binary shape mismatch for {name}")
        xv=struct.unpack('<'+str(len(x)//8)+'d',x);yv=struct.unpack('<'+str(len(y)//8)+'d',y)
        for u,v in zip(xv,yv): finite(u,f"NQ10 {name}");finite(v,f"F8a {name}")
        e=sum((u-v)**2 for u,v in zip(xv,yv));n=sum(v*v for v in yv);f8_kernels[name]=rel_l2(e,n)
    for v in [*f8_components.values(),*f8_trace.values(),*f8_spectrum.values(),*f8_kernels.values(),f8_measurement_norm]: finite(v,'NQ10/F8a parity')
    if max([*f8_components.values(),*f8_trace.values(),*f8_spectrum.values(),*f8_kernels.values(),f8_measurement_norm])>1e-11:
        raise ValueError("NQ10 DSpec reference does not match accepted C2/F8a within1e-11")
    reference_errors={}
    for name,cols in (('components.csv',('up_re','up_im','north_re','north_im','east_re','east_im')),
                      ('spectrum.csv',('real_m','imag_m')),
                      ('trace.csv',('velocity_m_s','rate_m_s2'))):
        exact=('bin','frequency_hz') if name=='spectrum.csv' else ('bin',) if name=='components.csv' else ('time_s','window')
        reference_errors[name]=scalar_channels(a.candidate/name,a.reference/name,cols,exact=exact)
    refnorm=abs(float(c['measurement_rate_norm'])-float(r['measurement_rate_norm']))/max(abs(float(r['measurement_rate_norm'])),1e-300)
    inv=invc
    field_errors=compare_common_fields(a.candidate,a.reference,inv)
    kernel_errors=compare_common_kernel(a.candidate,a.reference,inv)
    batch_max=proof(a.candidate,'phase-d-batched-scalar-equivalence.csv','same-data scalar/batched')
    ref_batch_max=proof(a.reference,'phase-d-batched-scalar-equivalence.csv','NQ10 same-data scalar/batched')
    tensor_max=proof(a.candidate,'phase-c2-common-reconstruction-equivalence.csv','common/tensor')
    ref_tensor_max=proof(a.reference,'phase-c2-common-reconstruction-equivalence.csv','NQ10 common/tensor')
    accuracy_components={k:max(v.values()) for k,v in reference_errors.items()}
    solids_f=[x for x in field_errors if not x['fluid']];solids_k=[x for x in kernel_errors if not x['fluid']]
    max_solid_field=max(max(x['forward_relative_l2'],x['adjoint_relative_l2']) for x in solids_f)
    max_solid_kernel=max(x['relative_l2'] for x in solids_k)
    max_wave=max(accuracy_components.values())
    requested_met=max(max_solid_field,max_solid_kernel,max_wave,refnorm)<=a.accuracy
    report={'status':'pass' if requested_met else 'accuracy_not_met_try_next_mesh',
      'requested_relative_accuracy':a.accuracy,'relative_error_input':float(c['relative_error']),
      'reference':{'nq':r['nq'],'h':r['refinement'],'accepted_F8a_components_relative_l2':f8_components,
        'accepted_F8a_spectrum_relative_l2':f8_spectrum,'accepted_F8a_trace_relative_l2':f8_trace,
        'accepted_F8a_measurement_norm_relative_error':f8_measurement_norm,'accepted_F8a_kernel_relative_l2':f8_kernels,
        'scalar_vs_batched_max_relative_l2':ref_batch_max,'common_vs_tensor_max_relative_l2':ref_tensor_max},
      'candidate':{'nq':c['nq'],'h':c['refinement'],'mesh_scale':c['phase_d_mesh_scale'],
        'measurement_rate_norm_relative_error':refnorm,'waveform_relative_l2':reference_errors,
        'common_solid_field_by_bin_section_location':field_errors,'common_full35_kernel_by_section_location':kernel_errors,
        'maximum_solid_field_relative_l2':max_solid_field,'maximum_solid_full35_kernel_relative_l2':max_solid_kernel,
        'maximum_waveform_relative_l2':max_wave,'common_points':len(inv),
        'scalar_vs_batched_max_relative_l2':batch_max,'common_vs_tensor_max_relative_l2':tensor_max,
        'angular_workspace_bound_bytes':c['phase_d_batched_cache_max_bytes'],
        'derived_angular_temporary_peak_bytes':c['phase_d_actual_temporary_peak_bytes']},
      'requested_accuracy_met':requested_met,
      'next_mesh_scale':None if requested_met else (0.9 if float(c['phase_d_mesh_scale'])==1.0 else 0.8 if float(c['phase_d_mesh_scale'])==0.9 else None),
      'limits':['fluid field/kernel groups are diagnostics separate from solid acceptance',
        'angular bytes are derived per-worker workspace accounting; supervisor RSS is the measured process resource',
        'only requested 1e-3 tolerance is validated in this phase']}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n');print(json.dumps(report,sort_keys=True))
if __name__=='__main__': main()
