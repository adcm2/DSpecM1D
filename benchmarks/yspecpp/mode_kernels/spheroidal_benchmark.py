"""Optional full-gravity 0S2 comparison using production DSpecM1D and YSpecpp kernels."""
from __future__ import annotations
import csv, json, math, shutil, sys, time
from pathlib import Path
import numpy as np
import run_benchmark as common
from spheroidal_reference import (A0, C0, F0, L0, N0, R, RHO, VS, VP,
                                  compare_orders, solve_reference)

BAND=(0.30,0.36)
DT=60.0
NT=8192
DURATION=72*3600.0
EPS=5.0/DURATION
DF=1.0/(NT*DT)
FIRST,END=146,178
BINS=np.arange(FIRST,END)
H=1e-4
WINDOW_SAMPLES=4321
MAXSTEP=.02
NQ=8
TOL=1e-7
MESH=1000.0
PARAMETERS=('A','C','F','L','N','rho')
GATE={'metadata_rel':5e-11,'baseline_frequency_rel':1e-5,'raw_spectrum_rel_l2':1e-5,
      'processed_spectrum_rel_l2':1e-5,'frequency_derivative_rel':1e-3,
      'cross_code_derivative_rel':1e-3,'profile_weighted_l2':1e-3,
      'centered_spectrum_rel_l2':1e-3,'fit_q_identity':1e-8,
      'reference_order_frequency_rel':1e-5,'reference_displacement_overlap_min':1-1e-6,
      'reference_boundary_residual':1e-7,'toroidal_leakage':1e-12}


def _read(path):
    with open(path,newline='') as f: return list(csv.DictReader(f))


def _curve_error(solver,path,reference):
    rows=_read(path)
    radius=np.array([float(r['radius_m']) for r in rows])
    x=radius/R
    if solver=='dspec': weight=np.array([float(r['weight_m']) for r in rows])/R
    else: weight=np.array([float(r['weight_m3'])/max(float(r['radius_m'])**2,1e-300) for r in rows])/R
    result={}
    for name in PARAMETERS:
        observed=np.array([float(r['radial_'+name]) for r in rows])*R
        ref=reference['curves'][name]
        result[name]=float(np.sqrt(np.sum((observed-ref)**2*weight)/np.sum(ref**2*weight)))
    return result


def _write_reference_profiles(out):
    x=np.linspace(0.0,1.0,501)
    curves=solve_reference(10,curve_x=x)['curves']
    with (out/'reference_radial_profiles.csv').open('w',newline='') as f:
        writer=csv.writer(f)
        writer.writerow(['x','radius_m',*(f'domega_dh_{name}_per_x_rad_s' for name in PARAMETERS)])
        for i,xi in enumerate(x):
            writer.writerow([xi,xi*R,*(curves[name][i] for name in PARAMETERS)])


def run_spheroidal(args):
    out=args.output.resolve()
    build=args.build_root.resolve()/'spheroidal'
    common.output_guard(out)
    deps_record={name:common.pinned_dependency((common.DEPS/name if name.endswith('-src') else
        common.DEPS/'planetmodel-src' if name=='planetmodel-dspec' else common.PLANET),commit)
        for name,commit in common.DEP_PINS.items()}
    env=common.env_bounded()
    model=common.HERE/'homogeneous_sphere.deck'
    ref8=solve_reference(8)
    ref10=solve_reference(10)
    refcheck=compare_orders(ref8,ref10)
    protocol={'reference':'independent DT full-gravity weak-energy Rayleigh-Ritz; order 10, order 8 cross-check',
      'fixture':{'radius_m':R,'density_kg_m3':RHO,'Vp_m_s':VP,'Vs_m_s':VS,'G_m3_kg_s2':6.6723e-11},
      'mode':{'name':'0S2','degree':2,'reference_frequency_mHz':ref10['frequency_mHz'],
              'reference_omega_rad_s':ref10['omega'],'regular_displacement_basis':'U,V = ell,k times r^(ell-1) coupled leading term plus independent r^(ell+1) Legendre(r^2) terms',
              'potential_basis':'P = r^ell Legendre(r^2)','basis_order_primary':10,'basis_order_check':8},
      'geometry':{'source_lat_lon_deg':[0,0],'source_depth_m':50000,'moment_tensor_Nm':[1e19,0,0,0,0,0],
                  'receiver_lat_lon_deg':[0,90],'receiver_depth_m':3000,'component':'Up'},
      'spectrum':{'band_mHz':list(BAND),'duration_s':DURATION,'dt_s':DT,'nt':NT,'df_Hz':DF,
                  'damping_s-1':EPS,'common_bins_inclusive':[FIRST,END-1],'window':'Hann, 4321 samples'},
      'reference_energy':'DT 8.208 elastic energy plus full-gravity mixed functional; eliminate P by stationarity with regular-center and exterior-decay boundary term',
      'density_kernel':'[8πGρr²U²−2grUt+2r²UP\u2032+2krVP−ω²r²(U²+V²)−8πGr²∫r^Rρ(s)U(s)t(s)/s ds]/(2ωI), t=2U−kV',
      'modulus_kernels':{'A':'t²/(2ωI)','C':'(rU\u2032)²/(2ωI)','F':'2rU\u2032t/(2ωI)','L':'s²/(2ωI)','N':'[(ell(ell+1)−2)V²−t²]/(2ωI)','s':'rV\u2032−V+kU'},
      'perturbations':{'h':H,'directions':list(PARAMETERS),'type':'separate centered uniform fractional changes; density recomputes g(r), fixed moduli'},
      'gates':GATE,'status':'predeclared'}
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    for solver,source,planet in [('dspec',args.dspec,common.DEPS/'planetmodel-src'),('yspecpp',args.yspec,common.PLANET)]:
        bdir=build/solver
        cm=['cmake','-S',common.HERE,'-B',bdir,'-DSOLVER='+('dspec' if solver=='dspec' else 'yspecpp'),
            '-DDSPEC_SOURCE_DIR='+str(args.dspec),'-DYSPEC_SOURCE_DIR='+str(args.yspec),
            '-DDEPENDENCY_ROOT='+str(common.DEPS),'-DPLANETMODEL_SOURCE_DIR='+str(planet),'-DCMAKE_BUILD_TYPE=Release']
        common.run(cm,env)
        common.run(['cmake','--build',bdir,'--target','kernel_adapter','-j1'],env)
    exeD=build/'dspec'/'kernel_adapter'
    exeY=build/'yspecpp'/'kernel_adapter'
    fixture=out/'homogeneous_sphere.deck'
    shutil.copy2(model,fixture)
    empty=out/'empty_q.csv'
    empty.write_text('basis,bin,real,imag\n')
    dirs={}
    start=time.monotonic()
    def dspec(label,typ,delta=0.0,prop='none',qfile=empty):
        dest=out/label
        common.run([exeD,'spectrum',fixture,dest,MAXSTEP,NQ,*BAND,delta,typ,prop,qfile,'up','spheroidal'],env)
        return dest
    dirs['d_s']=dspec('dspec_spheroidal_baseline',3)
    dirs['d_t']=out/'dspec_toroidal_preflight'
    common.run([exeD,'spectrum',fixture,dirs['d_t'],MAXSTEP,NQ,*BAND,0,2,'none',empty,'up','spheroidal'],env)
    dirs['y']=out/'yspecpp_spheroidal_baseline'
    common.run([exeY,'spectrum',fixture,dirs['y'],*BAND,0,TOL,empty,MESH,'spheroidal','none'],env)
    dsraw,dsproc,_=common.read_spectrum(dirs['d_s']/'spectrum.csv')
    traw,tproc,_=common.read_spectrum(dirs['d_t']/'spectrum.csv')
    yraw,yproc,_=common.read_spectrum(dirs['y']/'spectrum.csv')
    leakage=float(np.linalg.norm(tproc[BINS])/np.linalg.norm(dsproc[BINS]))
    raw_l2=float(np.linalg.norm(dsraw[BINS]-yraw[BINS])/np.linalg.norm(yraw[BINS]))
    proc_l2=float(np.linalg.norm(dsproc[BINS]-yproc[BINS])/np.linalg.norm(yproc[BINS]))
    fits={}
    cfg={'NT':NT,'DT':DT,'WINDOW_SAMPLES':WINDOW_SAMPLES,'EPS':EPS,'BINS':BINS,'DF':DF}
    for solver,proc in [('dspec',dsproc),('yspecpp',yproc)]:
        fits[solver]=common.fit_model(proc,cfg,omega_initial=ref10['omega'])
        common.write_q(out/f'q_{solver}.csv',fits[solver],BINS)
    kernel_runtime={}
    kernel_rows={}
    finite={}
    spectra={}
    for solver,exe,base_dir in [('dspec',exeD,dirs['d_s']),('yspecpp',exeY,dirs['y'])]:
        qfile=out/f'q_{solver}.csv'
        kd=out/f'{solver}_kernels'
        if solver=='dspec': cmd=[exe,'kernel',fixture,kd,MAXSTEP,NQ,*BAND,0,3,'none',qfile,'up','spheroidal']
        else: cmd=[exe,'kernel',fixture,kd,*BAND,0,TOL,qfile,MESH,'spheroidal','none']
        t0=time.monotonic()
        common.run(cmd,env)
        kernel_runtime[solver]=time.monotonic()-t0
        kernel_rows[solver]={r['basis']:r for r in _read(kd/'kernel_summary.csv')}
        for param in PARAMETERS:
            pair=[]
            for sign,label in ((1,'plus'),(-1,'minus')):
                dest=out/f'{solver}_{param}_{label}'
                if solver=='dspec': cmd=[exe,'spectrum',fixture,dest,MAXSTEP,NQ,*BAND,sign*H,3,param,empty,'up','spheroidal']
                else: cmd=[exe,'spectrum',fixture,dest,*BAND,sign*H,TOL,empty,MESH,'spheroidal',param]
                common.run(cmd,env)
                pair.append(common.read_spectrum(dest/'spectrum.csv')[1])
            coords=[float(kernel_rows[solver][f'q{i}']['distributed_'+param])+
                    float(kernel_rows[solver][f'q{i}'].get('source_'+param,0.0)) for i in range(7)]
            finite[solver,param]=common.fit_change(fits[solver],pair[0][BINS],pair[1][BINS],coords,H)
            spectra[solver,param]=(pair[0],pair[1])
    # Independently compare each q6 production radial curve with DT fractional-parameter profiles.
    reference_curves={}
    profile_errors={}
    for solver in ('dspec','yspecpp'):
        rows=_read(out/f'{solver}_kernels'/'kernel_q6.csv')
        radii=np.array([float(r['radius_m']) for r in rows])
        reference_curves[solver]=solve_reference(10,curve_x=radii/R)
        profile_errors[solver]=_curve_error(solver,out/f'{solver}_kernels'/'kernel_q6.csv',reference_curves[solver])
    parameter_results={}
    for p in PARAMETERS:
        ref=ref10['derivatives'][p]
        parameter_results[p]={'reference_domega_dh_rad_s':ref,'production':{}}
        for solver in ('dspec','yspecpp'):
            row=kernel_rows[solver]['q6']
            distributed=float(row['distributed_'+p])
            source=float(row.get('source_'+p,0.0))
            total=distributed+source
            parameter_results[p]['production'][solver]={'distributed':distributed,'source_point':source,'total':total,
                'relative_reference_error':abs(total-ref)/abs(ref),'profile_weighted_l2':profile_errors[solver][p],
                'centered_spectrum_residual':finite[solver,p]['residual_l2'],
                'finite_fit_omega_shift_rad_s':(common.fit_model(spectra[solver,p][0][BINS],cfg,ref10['omega'])['omega']-
                    common.fit_model(spectra[solver,p][1][BINS],cfg,ref10['omega'])['omega'])/2,
                'kernel_runtime_s':kernel_runtime[solver]}
    metaD=json.loads((dirs['d_s']/'meta.json').read_text())
    metaY=json.loads((dirs['y']/'meta.json').read_text())
    meta_errors={}
    for key in ('nt','dt_s','df_hz','damping_s','source_radius_m'):
        meta_errors[key]=abs(float(metaD[key])-float(metaY[key]))/max(abs(float(metaD[key])),abs(float(metaY[key])),1e-300)
    meta_errors['first_bin']=0.0 if int(metaD['first'])==int(metaY['first']) else math.inf
    meta_errors['end_bin_exclusive']=0.0 if int(metaD['end'])==int(metaY['last'])+1 else math.inf
    gates={
      'metadata_grid':{'value':max(meta_errors.values()),'limit':GATE['metadata_rel']},
      'baseline_frequency_dspec':{'value':abs(fits['dspec']['fit_frequency_rel_bias']),'limit':GATE['baseline_frequency_rel']},
      'baseline_frequency_yspecpp':{'value':abs(fits['yspecpp']['fit_frequency_rel_bias']),'limit':GATE['baseline_frequency_rel']},
      'baseline_raw_spectrum':{'value':raw_l2,'limit':GATE['raw_spectrum_rel_l2']},
      'baseline_processed_spectrum':{'value':proc_l2,'limit':GATE['processed_spectrum_rel_l2']},
      'toroidal_leakage':{'value':leakage,'limit':GATE['toroidal_leakage']},
      'reference_order_frequency':{'value':refcheck['frequency_relative_change'],'limit':GATE['reference_order_frequency_rel']},
      'reference_displacement_overlap':{'value':1-refcheck['displacement_overlap'],'limit':1-GATE['reference_displacement_overlap_min']},
      'reference_boundary':{'value':ref10['boundary_residual'],'limit':GATE['reference_boundary_residual']},
      'fit_q_identity_dspec':{'value':fits['dspec']['q_identity'],'limit':GATE['fit_q_identity']},
      'fit_q_identity_yspecpp':{'value':fits['yspecpp']['q_identity'],'limit':GATE['fit_q_identity']}}
    for p in PARAMETERS:
        vals=[parameter_results[p]['production'][s]['relative_reference_error'] for s in ('dspec','yspecpp')]
        cross=abs(parameter_results[p]['production']['dspec']['total']-parameter_results[p]['production']['yspecpp']['total'])/abs(ref10['derivatives'][p])
        gates[f'{p}_reference_derivative']={'value':max(vals),'limit':GATE['frequency_derivative_rel']}
        gates[f'{p}_cross_code_derivative']={'value':cross,'limit':GATE['cross_code_derivative_rel']}
        for solver in ('dspec','yspecpp'):
            gates[f'{solver}_{p}_profile']={'value':profile_errors[solver][p],'limit':GATE['profile_weighted_l2']}
            gates[f'{solver}_{p}_centered_spectrum']={'value':finite[solver,p]['residual_l2'],'limit':GATE['centered_spectrum_rel_l2']}
    for record in gates.values(): record['pass']=bool(record['value']<=record['limit'])
    result={'protocol':protocol,'reference_checks':{'order8_frequency_mHz':ref8['frequency_mHz'],'order10_frequency_mHz':ref10['frequency_mHz'],**refcheck,
       'order10_boundary_residual':ref10['boundary_residual']},'metadata':{'dspec':metaD,'yspecpp':metaY,'relative_errors':meta_errors},
       'preflight':{'toroidal_over_spheroidal_processed_l2':leakage,'raw_spectrum_relative_l2':raw_l2,'processed_spectrum_relative_l2':proc_l2},
       'fits':{s:{k:v for k,v in f.items() if k not in ('matrix','qrows','fit_target')} for s,f in fits.items()},
       'analytic_reference':{'frequency_mHz':ref10['frequency_mHz'],'omega_rad_s':ref10['omega'],'derivatives_rad_s':ref10['derivatives']},
       'parameter_results':parameter_results,'finite_changes':{f'{s}_{p}':v for (s,p),v in finite.items()},'gates':gates,
       'elapsed_numeric_s':time.monotonic()-start,
       'provenance':{'DSpecM1D':common.identity(args.dspec),'YSpecpp':common.identity(args.yspec),
          'pinned_dependencies':deps_record,'source_sha256':{p.name:common.sha(p) for p in (common.HERE/'dspec_adapter.cpp',common.HERE/'yspecpp_adapter.cpp',common.HERE/'run_benchmark.py',common.HERE/'spheroidal_reference.py',common.HERE/'spheroidal_benchmark.py')}}}
    _write_reference_profiles(out)
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    _write_report(out,result)
    print(json.dumps({'status':'PASS' if all(g['pass'] for g in gates.values()) else 'FAIL','gates':gates},indent=2))
    if not all(g['pass'] for g in gates.values()): raise SystemExit(1)


def _write_report(out,result):
    import html
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    profile_rows=_read(out/'reference_radial_profiles.csv')
    ref_x=np.array([float(row['x']) for row in profile_rows])
    fig,axes=plt.subplots(3,2,figsize=(11,11),sharex=True)
    colors={'dspec':'tab:blue','yspecpp':'tab:orange'}
    for ax,name in zip(axes.flat,PARAMETERS):
        ref=np.array([float(row[f'domega_dh_{name}_per_x_rad_s']) for row in profile_rows])
        ax.plot(ref_x,ref,color='black',lw=2,label='DT Ritz reference')
        for solver in ('dspec','yspecpp'):
            rows=_read(out/f'{solver}_kernels'/'kernel_q6.csv')
            x=np.array([float(row['radius_m'])/R for row in rows])
            values=np.array([float(row[f'radial_{name}'])*R for row in rows])
            ax.plot(x,values,color=colors[solver],lw=1,label=solver)
        ax.set_title(name)
        ax.set_ylabel(r'$d\omega/dh$ per $d(r/R)$ (rad s$^{-1}$)')
        ax.grid(alpha=.25)
    for ax in axes[-1,:]: ax.set_xlabel('normalized radius r/R')
    axes.flat[0].legend(fontsize=8)
    fig.suptitle('0S2 fractional-parameter derivative profiles',y=1.01)
    fig.tight_layout()
    fig.savefig(out/'summary.png',dpi=160,bbox_inches='tight')
    plt.close(fig)

    frequencies= BINS * DF * 1e3
    spectral_data={}
    for name in PARAMETERS:
        rows=[]
        for solver in ('dspec','yspecpp'):
            record=result['finite_changes'][f'{solver}_{name}']
            actual=np.asarray(record['finite_delta'])
            predicted=np.asarray(record['predicted_delta'])
            n=len(BINS)
            spectral_data[solver,name]=(actual,predicted)
        for i,frequency in enumerate(frequencies):
            row=[frequency]
            for solver in ('dspec','yspecpp'):
                actual,predicted=spectral_data[solver,name]
                row.extend((actual[i],actual[n+i],predicted[i],predicted[n+i]))
            rows.append(row)
        with (out/f'centered_spectrum_{name}.csv').open('w',newline='') as f:
            writer=csv.writer(f)
            writer.writerow(['frequency_mHz','dspec_actual_real','dspec_actual_imag',
                'dspec_kernel_real','dspec_kernel_imag','yspecpp_actual_real',
                'yspecpp_actual_imag','yspecpp_kernel_real','yspecpp_kernel_imag'])
            writer.writerows(rows)

    for solver in ('dspec','yspecpp'):
        fig,axes=plt.subplots(3,2,figsize=(11,10),sharex=True)
        for ax,name in zip(axes.flat,PARAMETERS):
            actual,predicted=spectral_data[solver,name]
            n=len(BINS)
            ax.plot(frequencies,actual[:n],color='tab:blue',label='centered change, real')
            ax.plot(frequencies,actual[n:],color='tab:blue',ls='--',label='centered change, imaginary')
            ax.plot(frequencies,predicted[:n],color='tab:orange',label='kernel prediction, real')
            ax.plot(frequencies,predicted[n:],color='tab:orange',ls='--',label='kernel prediction, imaginary')
            ax.axhline(0,color='black',lw=.5)
            ax.set_title(name)
            ax.grid(alpha=.25)
        for ax in axes[-1,:]: ax.set_xlabel('frequency (mHz)')
        axes.flat[0].set_ylabel('centered complex ΔS (m s)')
        axes.flat[2].set_ylabel('centered complex ΔS (m s)')
        axes.flat[4].set_ylabel('centered complex ΔS (m s)')
        axes.flat[0].legend(fontsize=7)
        fig.suptitle(f'{solver}: centered spectrum and q-kernel prediction',y=1.01)
        fig.tight_layout()
        fig.savefig(out/f'{solver}_centered_spectrum.png',dpi=160,bbox_inches='tight')
        plt.close(fig)
    params=''.join(f"<tr><td>{p}</td><td>{result['analytic_reference']['derivatives_rad_s'][p]:.9g}</td>"+
      ''.join(f"<td>{result['parameter_results'][p]['production'][s]['total']:.9g}</td><td>{result['parameter_results'][p]['production'][s]['relative_reference_error']:.3g}</td>" for s in ('dspec','yspecpp'))+
      f"<td>{result['parameter_results'][p]['production']['dspec']['profile_weighted_l2']:.3g}</td><td>{result['parameter_results'][p]['production']['yspecpp']['profile_weighted_l2']:.3g}</td></tr>" for p in PARAMETERS)
    gates=''.join(f"<tr><td>{html.escape(k)}</td><td>{v['value']:.4g}</td><td>{v['limit']:.4g}</td><td>{'PASS' if v['pass'] else 'FAIL'}</td></tr>" for k,v in result['gates'].items())
    status='PASS' if all(g['pass'] for g in result['gates'].values()) else 'FAIL'
    body=f'''<!doctype html><html><head><meta charset="utf-8"><title>Full-gravity spheroidal modal kernel benchmark</title>
<style>body{{font:15px system-ui,sans-serif;max-width:1200px;margin:2rem auto;padding:0 1rem;color:#222}}table{{border-collapse:collapse;width:100%;margin:1rem 0}}th,td{{border:1px solid #bbb;padding:.4rem;text-align:left}}code{{background:#eee;padding:.1rem .25rem}}</style></head><body>
<h1>Full-gravity spheroidal modal kernel benchmark: {status}</h1>
<p>The independent homogeneous-sphere DT weak-energy Rayleigh–Ritz reference isolates 0S2 (ℓ=2) with regular displacement bases U,V = coupled leading terms proportional to r^(ℓ−1) plus independent r^(ℓ+1) Legendre terms, and potential basis P = r^ℓ Legendre(r²). The order-10 Ritz frequency is {result['analytic_reference']['frequency_mHz']:.10g} mHz; order 8 is the independent basis-size check. DT analytic kernels are evaluated on this numerical Ritz eigenmode, so the reported frequency is not a closed-form result.</p>
<p>Fixture: R=6371 km, ρ=5515 kg m⁻³, Vp=8000 m s⁻¹, Vs=4500 m s⁻¹, G=6.6723×10⁻¹¹ m³ kg⁻¹ s⁻². A pure Mrr source at (0°,0°), depth 50 km, and vertical receiver at (0°,90°), depth 3 km, select the spheroidal contribution. DSpec toroidal-only leakage is {result['preflight']['toroidal_over_spheroidal_processed_l2']:.3g}.</p>
<p>The reference uses the DT-normalized V convention with k=√(ℓ(ℓ+1)), t=2U−kV, s=rV′−V+kU, and I=∫ρ(U²+V²)r²dr. Its elastic weak-energy terms are C(rU′)² + 2F rU′t + (A−N)t² + Ls² + (ℓ(ℓ+1)−2)NV². The mixed gravitational functional includes 4πGρ²r²U², −2ρgrUt, 2ρ(r²UP′+krVP), the potential-gradient terms, and the exterior-decay boundary term (ℓ+1)RP(R)²/(4πG); P is eliminated by stationarity.</p>
<p>These weak-energy expressions follow Dahlen and Tromp (1998), <i>Theoretical Global Seismology</i>, Eqs. (8.201) and (8.208).</p>
<p>For A,C,F,L,N, the radial kernels are respectively t², (rU′)², 2rU′t, s², and [(ℓ(ℓ+1)−2)V²−t²], each divided by 2ωI. The fixed-modulus density kernel is [8πGρr²U²−2grUt+2r²UP′+2krVP−ω²r²(U²+V²)−8πGr²∫<sub>r</sub><sup>R</sup>ρ(s)U(s)t(s)/s ds]/(2ωI). It includes inertia and explicit gravity, while the integral is the nonlocal δg contribution. Density changes g(r); moduli remain fixed. The q6 source-point terms for moduli are scaled by the source-point modulus in pascals, and density has no point term.</p>
<p>For each of A,C,F,L,N,ρ, separate centered uniform fractional perturbations h=10⁻⁴ are propagated through both production solvers. The spectrum uses a 72-hour record sampled every 60 seconds, 0.30–0.36 mHz (bins 146–177), a 4321-sample Hann window, and damping 5/259200 s⁻¹. Six nuisance coordinates describe complex amplitude and smooth background changes; the seventh is the fitted angular-frequency tangent. Each parameter’s complete complex centered change is compared with its six production kernel q-coordinate derivatives and modal-frequency derivative. Radial profiles are plotted as fractional-frequency derivative per unit r/R.</p>
<table><thead><tr><th>Parameter</th><th>DT δω/dh (rad/s)</th><th>DSpec derivative</th><th>Relative error</th><th>YSpec derivative</th><th>Relative error</th><th>DSpec profile L2</th><th>YSpec profile L2</th></tr></thead><tbody>{params}</tbody></table>
<img src="summary.png" alt="DT, DSpecM1D and YSpecpp radial modal-kernel profiles for all six parameters" style="width:100%;height:auto">
<h2>Predeclared gates</h2><table><thead><tr><th>Gate</th><th>Value</th><th>Limit</th><th>Result</th></tr></thead><tbody>{gates}</tbody></table>
<p>Centered spectral changes are compared as complete complex vectors; the panels below show real and imaginary parts separately.</p>
<img src="dspec_centered_spectrum.png" alt="DSpecM1D centered spectra and kernel predictions for all six parameters" style="width:100%;height:auto">
<img src="yspecpp_centered_spectrum.png" alt="YSpecpp centered spectra and kernel predictions for all six parameters" style="width:100%;height:auto">
<p>Saved evidence includes <a href="results.json">full-precision results and provenance</a>, <a href="protocol.json">protocol and frozen gates</a>, <a href="summary.png">all six radial profiles</a>, centered-spectrum plots and per-direction CSVs, <a href="reference_radial_profiles.csv">sampled DT reference profiles</a>, solver spectra, q weights, and radial-kernel CSVs. The run is bounded to one homogeneous sphere, degree 2, basis orders 8 and 10, and a 1 km YSpecpp radial mesh. It does not claim general model or mesh convergence. Reproduce in a new empty output directory with <code>python3 benchmarks/yspecpp/mode_kernels/run_benchmark.py --mode spheroidal --output build/mode-frequency-kernels-20261024/spheroidal-reproduction-1 --build-root /tmp/modal-band-build --dspec "$PWD" --yspec "$(dirname "$PWD")/YSpecpp"</code>. Exact source and pinned-dependency identities are in results.json.</p>
</body></html>'''
    (out/'report.html').write_text(body)
