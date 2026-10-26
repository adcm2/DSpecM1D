#!/usr/bin/env python3
"""Bounded homogeneous-sphere modal-frequency kernel benchmark."""
from __future__ import annotations
import argparse, csv, hashlib, html, json, math, os, shutil, subprocess, time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import spherical_jn

HERE = Path(__file__).resolve().parent
DSPEC = HERE.parents[2] if len(HERE.parents)>2 else Path('/home/adcm2/Documents/Research/Projects/DSpecM1D')
YSPEC = DSPEC.parent / "YSpecpp"
DEPS = DSPEC / "build/phase4-20261021/_deps"
PLANET = YSPEC / "build/library/_deps/planetmodel-src"
DEP_PINS = {"numericconcepts-src":"e03a2889d43a03660d0879f46dda1406d90c9517",
            "interpolation-src":"d5a1711ca6a89fc8ad8148fe13212398fe569bc0",
            "fftwpp-src":"cb152c0e85f3822493c2d9e85a56307ac7b4f1c0",
            "gaussquad-src":"7d4786d192f68077cd99a34dcdf1379ed1eeac9a",
            "gshtrans-src":"8dd116df8df4222fe9b10a40c7213d34f2a274fb",
            "planetmodel-dspec":"35be4ee54b3d84df1483e647878c5446a886c90b",
            "planetmodel-yspec":"fa6e230fde25152d40f27d9839de44c4f3182a24"}
BAND = (0.24, 0.32)                 # mHz
DT, DURATION, NT = 60.0, 72*3600.0, 8192
EPS = 5.0 / DURATION
DF = 1.0 / (NT*DT)
FIRST, END, BINS = 116, 159, np.arange(116,159)
WINDOW_SAMPLES, H = 4321, 1e-4
R, RHO, VS, VP = 6_371_000.0, 5515.0, 4500.0, 8000.0
GATE = {"metadata_rel":5e-11,"baseline_frequency_rel":1e-5,
        "raw_spectrum_rel_l2":1e-5,"processed_spectrum_rel_l2":1e-5,
        "centered_linear_prediction_rel_l2":1e-3,"analytic_omega_rel":1e-3,
        "cross_code_omega_rel":1e-3,"radial_LN_rel":1e-3,
        "sensitivity_identity":1e-8,"spheroidal_leakage":1e-12}


def run(cmd, env):
    print("+", " ".join(map(str,cmd)), flush=True)
    subprocess.run(list(map(str,cmd)), check=True, env=env)

def readcsv(path):
    with open(path,newline="") as f: return list(csv.DictReader(f))

def complex_rows(path, field="processed"):
    rows=readcsv(path)
    return np.array([complex(float(r[field+"_real"]),float(r[field+"_imag"])) for r in rows])

def output_guard(out):
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise SystemExit(f"output directory must be empty: {out}")
    out.mkdir(parents=True,exist_ok=True)

def env_bounded():
    e=os.environ.copy(); e.update(OMP_NUM_THREADS="2",OPENBLAS_NUM_THREADS="2",MKL_NUM_THREADS="2")
    return e

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def git(root,*args): return subprocess.check_output(["git","-C",str(root),*args],text=True).strip()
def identity(root):
    status=subprocess.check_output(["git","-C",str(root),"status","--porcelain","--untracked-files=all"],text=True).rstrip()
    dirty={}
    for line in status.splitlines():
        p=root/line[3:]
        if p.is_file(): dirty[line[3:]]=sha(p)
    return {"path":str(root),"commit":git(root,"rev-parse","HEAD"),"status_porcelain":status,"dirty_file_sha256":dirty}
def dependency_identity(path): return {"path":str(path),"commit":git(path,"rev-parse","HEAD"),"status":git(path,"status","--porcelain","--untracked-files=all")}
def pinned_dependency(path, expected):
    record=dependency_identity(path)
    if record['commit']!=expected or record['status']:
        raise SystemExit(f"dependency pin or cleanliness mismatch: {record}")
    return record

def write_report(out,result):
    os.environ.setdefault('MPLCONFIGDIR',str(out/'.mplconfig'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4.6),constrained_layout=True)
    colors={'dspec':'#1769aa','yspecpp':'#c45500'}
    for solver in ('dspec','yspecpp'):
        rows=readcsv(out/f'{solver}_spectrum_change.csv')
        f=np.array([float(r['frequency_hz'])*1e3 for r in rows])
        actual=np.array([complex(float(r['centered_delta_real']),float(r['centered_delta_imag'])) for r in rows])
        predicted=np.array([complex(float(r['kernel_prediction_real']),float(r['kernel_prediction_imag'])) for r in rows])
        axes[0].plot(f,np.abs(actual),color=colors[solver],lw=1.5,label=f'{solver}: centered forward')
        axes[0].plot(f,np.abs(predicted),color=colors[solver],lw=1,ls='--',label=f'{solver}: SK prediction')
    axes[0].set(xlabel='Frequency (mHz)',ylabel='|processed ΔS| (m s)',title='Centered h=10⁻⁴ spectrum change')
    axes[0].legend(fontsize=7); axes[0].grid(alpha=.25)
    analytic=readcsv(out/'dt_analytic_radial.csv')
    ra=np.array([float(r['radius_m'])/1e6 for r in analytic])
    for param,style in [('L','-'),('N','--')]:
        key='domega_d'+param+'_per_m'
        ya=np.array([float(r[key])*1e6 for r in analytic])
        axes[1].plot(ra,ya,color='black',ls=style,lw=1.5,label=f'DT {param}')
        for solver in ('dspec','yspecpp'):
            rows=readcsv(out/f'{solver}_kernels/kernel_q6.csv')
            stride=max(1,len(rows)//350)
            col='radial_'+param
            radii=np.array([float(r['radius_m'])/1e6 for r in rows[::stride]])
            values=np.array([float(r[col])*1e6 for r in rows[::stride]])
            axes[1].plot(radii,values,color=colors[solver],ls=style,lw=.8,alpha=.8,
                         label=f'{solver} {param}')
    axes[1].set(xlabel='Radius (Mm)',ylabel='dω/dh radial integrand (rad s⁻¹ Mm⁻¹)',title='Uniform L and N frequency-kernel integrands')
    axes[1].legend(fontsize=7,ncol=2); axes[1].grid(alpha=.25)
    fig.savefig(out/'summary.png',dpi=160); plt.close(fig)
    gates=''.join(f"<tr><td>{html.escape(k)}</td><td>{v['value']:.4g}</td><td>{v.get('limit',v.get('minimum_detectable',''))}</td><td>{'PASS' if v['pass'] else 'FAIL'}</td></tr>" for k,v in result['gates'].items())
    fits=''.join(f"<tr><td>{s}</td><td>{result['fits'][s]['omega']:.15g}</td><td>{result['fits'][s]['fit_frequency_rel_bias']:.4g}</td><td>{result['fits'][s]['residual_l2']:.4g}</td><td>{result['kernel_results'][s]['domega_dh']:.15g}</td><td>{result['kernel_results'][s]['kernel_runtime_s']:.3f}</td></tr>" for s in ('dspec','yspecpp'))
    shifts=''.join(f"<tr><td>{s}</td><td>{result['centered_finite_difference'][s]['kernel_predicted_omega_shift']:.12g}</td><td>{result['centered_finite_difference'][s]['finite_fit_omega_shift']:.12g}</td><td>{result['centered_finite_difference'][s]['finite_exact_omega_shift_expected']:.12g}</td><td>{result['centered_finite_difference'][s]['residual_l2']:.4g}</td></tr>" for s in ('dspec','yspecpp'))
    curves=''.join(f"<tr><td>{s}</td><td>{p}</td><td>{result['radial_curve_comparison'][s][p]['weighted_L2']:.4g}</td><td>{result['radial_curve_comparison'][s][p]['peak_scaled_max']:.4g}</td></tr>" for s in ('dspec','yspecpp') for p in ('L','N'))
    status='PASS' if all(v['pass'] for v in result['gates'].values()) else 'FAIL'
    body=f'''<!doctype html><html><head><meta charset="utf-8"><title>Modal-frequency kernel benchmark</title>
<style>body{{font:15px system-ui,sans-serif;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#222}}table{{border-collapse:collapse;width:100%;margin:1rem 0}}th,td{{border:1px solid #bbb;padding:.4rem;text-align:left}}img{{max-width:100%}}code{{background:#eee;padding:.1rem .25rem}}</style></head><body>
<h1>Narrow-band modal-frequency kernel benchmark: {status}</h1>
<p>Homogeneous isotropic sphere; 0T2, ℓ=2, North displacement at an equatorial receiver 90° from the equatorial source. DSpecM1D uses its toroidal family. YSpecpp computes both families; the spheroidal North angular factor is zero at this geometry, verified by DSpec S-only leakage of {result['preflight']['spheroidal_processed_l2_over_toroidal']:.3g}.</p>
<p>R=6371 km; ρ=5515 kg m⁻³; Vp=8000 m s⁻¹; Vs=4500 m s⁻¹. Source: depth 50 km, Mrt=10¹⁹ N m. Receiver: depth 3 km. Band 0.24–0.32 mHz, 72 h, dt=60 s, nt=8192, common bins 116–158, damping 5/259200 s⁻¹. A 4321-sample Hann window is applied after inverse FFT and damping removal. All other positive-frequency bins are zero in the processed template.</p>
<p>The processed complex spectrum is fit with the isolated pole <code>a/[i z (ω²−z²)]</code> plus a linear complex background in <code>z</code>. The real-linear processing operator is applied separately to each real and imaginary nuisance column. Seven tangent coordinates comprise six nuisance coordinates and fitted ω. Production SKs are evaluated for each coordinate using q weights from the fitted real-linear inverse. The centered material perturbation is h=10⁻⁴ with L=N=μ scaled uniformly and A,C,F,ρ fixed.</p>
<img src="summary.png" alt="Processed spectral change and radial frequency kernels">
<h2>Frequency result</h2><p>Independent Bessel root/quadrature: x={result['analytic_reference']['x']:.15g}, ω={result['analytic_reference']['omega']:.15g} rad s⁻¹ ({result['analytic_reference']['frequency_mHz']:.12g} mHz), dω/dh=ω/2={result['analytic_reference']['omega_over_2']:.15g} rad s⁻¹. The exact centered shift is reported separately in <a href="results.json">results.json</a>.</p>
<p>The reference uses W=j₂(xr/R), I=∫ρW²r²dr, κL=r²(W′−W/r)²/(2ωI), and κN=4W²/(2ωI), with δω=∫[κLδL+κNδN]dr. Plotted profiles are LκL and NκN (here L=N=μ), so integrating them against fractional changes gives δω; the saved totals also include source-point terms.</p>
<table><thead><tr><th>Solver</th><th>Fitted baseline ω (rad s⁻¹)</th><th>Relative baseline bias</th><th>Complex fit residual</th><th>Production dω/dh (rad s⁻¹)</th><th>Kernel runtime (s)</th></tr></thead><tbody>{fits}</tbody></table>
<h2>Predeclared gates</h2><table><thead><tr><th>Gate</th><th>Value</th><th>Limit</th><th>Result</th></tr></thead><tbody>{gates}</tbody></table>
<h2>Centered finite change versus linear prediction</h2><p>At h=10⁻⁴ the tangent prediction uses h(dω/dh) as its linear shift. The finite-model column uses the fitted half difference between +h and −h spectra; its exact uniform-scaling reference is ω[√(1+h)−√(1−h)]/2. The complex-spectrum residual includes every measured bin, including the portion outside the seven-coordinate tangent subspace.</p>
<table><thead><tr><th>Solver</th><th>SK linear h·dω/dh (rad s⁻¹)</th><th>Fitted centered finite shift (rad s⁻¹)</th><th>Exact centered reference (rad s⁻¹)</th><th>Full complex ΔS residual L2</th></tr></thead><tbody>{shifts}</tbody></table>
<h2>Radial curve checks</h2><table><thead><tr><th>Solver</th><th>Modulus</th><th>Weighted profile L2</th><th>Peak scaled error</th></tr></thead><tbody>{curves}</tbody></table>
<p>Separate distributed and source-point L/N integrals are in <a href="results.json">results.json</a>; per-bin spectral changes are in <code>dspec_spectrum_change.csv</code> and <code>yspecpp_spectrum_change.csv</code>. Seven q weights per solver are saved in <code>q_dspec.csv</code> and <code>q_yspecpp.csv</code>. YSpecpp's public calculation API has no family selector, so the geometric S cancellation is specific to this component/geometry.</p>
<p>Bound: one homogeneous elastic toroidal mode only; this does not test spheroidal/density/attenuation kernels or model/mesh convergence. Reproduce into a new empty directory from the repository root with <code>python3 benchmarks/yspecpp/mode_kernels/run_benchmark.py --output build/mode-frequency-kernels-20261024/reproduction-1 --build-root /tmp/mode-frequency-kernels-build --dspec "$PWD" --yspec "$(dirname "$PWD")/YSpecpp"</code>.</p>
<p>Saved evidence: <a href="results.json">results.json</a>, <a href="protocol.json">protocol.json</a>, and source/dependency identity records inside the results file.</p>
</body></html>'''
    (out/'report.html').write_text(body)

def read_spectrum(path):
    rows=readcsv(path)
    raw=np.zeros(NT//2+1,complex); proc=np.zeros_like(raw)
    for r in rows:
        k=int(r['bin']); raw[k]=complex(float(r['raw_real']),float(r['raw_imag']))
        proc[k]=complex(float(r['processed_real']),float(r['processed_imag']))
    return raw,proc,rows

def P(raw, cfg=None):
    cfg=cfg or {}; nt=cfg.get('NT',NT); dt=cfg.get('DT',DT); nw=cfg.get('WINDOW_SAMPLES',WINDOW_SAMPLES); eps=cfg.get('EPS',EPS)
    n=np.arange(nt)*dt
    win=np.zeros(nt); win[:nw]=.5-.5*np.cos(2*np.pi*np.arange(nw)/(nw-1))
    return np.fft.rfft(np.fft.irfft(raw)*np.exp(eps*n)*win)

def fit_model(proc, cfg=None, omega_initial=None):
    cfg=cfg or {}; nt=cfg.get('NT',NT); dt=cfg.get('DT',DT); nw=cfg.get('WINDOW_SAMPLES',WINDOW_SAMPLES)
    eps=cfg.get('EPS',EPS); bins=np.asarray(cfg.get('BINS',BINS)); df=cfg.get('DF',DF)
    # Pole and two smooth complex background columns all share the physical 1/(i z) factor.
    freq=bins*df; z=2*np.pi*freq-1j*eps
    observed=proc[bins] if len(proc)>bins[-1] else proc
    target=np.r_[observed.real,observed.imag]
    def model(omega, with_columns=False):
        bases=[1/(1j*z*(omega**2-z**2)),1/(1j*z),(z-omega)/(1j*z)]
        columns=[]
        for basis in bases:
            full=np.zeros(nt//2+1,complex); full[bins]=basis
            p=P(full,cfg)[bins]
            full_i=np.zeros_like(full); full_i[bins]=1j*basis
            pi=P(full_i,cfg)[bins]
            columns.extend((np.r_[p.real,p.imag],np.r_[pi.real,pi.imag]))
        A=np.column_stack(columns)
        scale=np.linalg.norm(A,axis=0); As=A/scale
        xscaled=np.linalg.lstsq(As,target,rcond=None)[0]
        x=xscaled/scale
        residual=target-A@x
        if with_columns: return float(np.vdot(residual,residual).real),x,A
        return float(np.vdot(residual,residual).real)
    # analytic root is only the initial bracket; the fit remains independent of that value.
    if omega_initial is None:
        x=brentq(lambda v:v*spherical_jn(2,v,derivative=True)-spherical_jn(2,v),2.4,2.6)
        omega0=x*VS/R
    else: omega0=float(omega_initial)
    fit=minimize_scalar(model,bounds=(omega0*.97,omega0*1.03),method="bounded",
                        options={"xatol":1e-15,"maxiter":250})
    residual,coef,A=model(fit.x,True)
    # Real-linear coordinate tangent: six nuisance columns plus fitted-frequency tangent.
    a=complex(coef[0],coef[1])
    deriv=-2*fit.x/(1j*z*(fit.x**2-z**2)**2)
    # Because P is real-linear, evaluate complex coefficient parts separately.
    # a*g' = a.real*g' + a.imag*(i*g')
    gfull=np.zeros(nt//2+1,complex); gfull[bins]=deriv
    igfull=np.zeros_like(gfull); igfull[bins]=1j*deriv
    pg=P(gfull,cfg)[bins]; pig=P(igfull,cfg)[bins]
    t=a.real*np.r_[pg.real,pg.imag]+a.imag*np.r_[pig.real,pig.imag]
    basis_norms=np.linalg.norm(A,axis=0); t_norm=np.linalg.norm(t)
    nuisance=A/basis_norms[None,:]
    D=np.column_stack((nuisance,t))
    Dnorm=np.column_stack((nuisance,t/t_norm))
    qscaled=np.linalg.pinv(Dnorm,rcond=1e-13)
    qrows=qscaled.copy(); qrows[6,:]/=t_norm
    result={"omega":float(fit.x),"omega_analytic":float(omega0),"residual_l2":math.sqrt(residual/np.vdot(target,target).real),
            "coefficients":coef.tolist(),"basis_condition":float(np.linalg.cond(Dnorm)),
            "q_identity":float(np.max(np.abs(qrows@D-np.eye(7)))),
            "matrix":D.tolist(),"qrows":qrows.tolist(),"fit_target":target.tolist()}
    result["amplitude"]=[float(coef[0]),float(coef[1])]
    result["fit_frequency_rel_bias"]=(float(fit.x)-omega0)/omega0
    return result

def write_q(path, fit, bins=None):
    bins=np.asarray(BINS if bins is None else bins)
    q=np.asarray(fit['qrows'])
    with open(path,'w',newline='') as f:
        w=csv.writer(f); w.writerow(['basis','bin','real','imag'])
        for j in range(7):
            for ix,k in enumerate(bins):
                w.writerow([f'q{j}',int(k),float(q[j,ix]),float(q[j,len(bins)+ix])])

def fit_change(fit, plus, minus, kernel_coordinates, h=H):
    delta=np.r_[((plus-minus)/2).real,((plus-minus)/2).imag]
    # qrows project the finite signal change; production SKs independently predict the coordinates.
    observed=np.asarray(fit['qrows'])@delta
    coords=h*np.asarray(kernel_coordinates)
    D=np.asarray(fit['matrix']); predicted=D@coords
    subspace=D@np.linalg.lstsq(D,delta,rcond=None)[0]
    return {"projected_finite_coordinates":observed.tolist(),"kernel_predicted_coordinates":coords.tolist(),
            "finite_delta":delta.tolist(),"predicted_delta":predicted.tolist(),
            "residual_l2":float(np.linalg.norm(delta-predicted)/np.linalg.norm(delta)),
            "outside_subspace_l2":float(np.linalg.norm(delta-subspace)/np.linalg.norm(delta)),
            "projected_coordinate_rel_l2":float(np.linalg.norm(observed-coords)/np.linalg.norm(observed))}

def analytic_reference():
    x=brentq(lambda v:v*spherical_jn(2,v,derivative=True)-spherical_jn(2,v),2.4,2.6)
    omega=x*VS/R
    # Gaussian quadrature for the independent DT radial integrals.
    nodes,weights=np.polynomial.legendre.leggauss(600)
    r=(nodes+1)*R/2; wr=weights*R/2
    W=spherical_jn(2,x*r/R); dW=(x/R)*spherical_jn(2,x*r/R,derivative=True)
    I=np.sum(wr*RHO*W**2*r**2)
    kL=r**2*(dW-W/r)**2/(2*omega*I)
    kN=(2*(2+1)-2)*W**2/(2*omega*I)
    dL=np.sum(wr*kL*RHO*VS**2)
    dN=np.sum(wr*kN*RHO*VS**2)
    radial={"radius_m":r.tolist(),"domega_dL_per_m":(kL*RHO*VS**2).tolist(),
            "domega_dN_per_m":(kN*RHO*VS**2).tolist()}
    return {"x":float(x),"omega":float(omega),"frequency_mHz":float(omega/(2*np.pi)*1e3),
            "domega_duniform_L_fraction":float(dL),"domega_duniform_N_fraction":float(dN),
            "domega_duniform_LN_fraction":float(dL+dN),"omega_over_2":float(omega/2),
            "quadrature_identity_rel":float(abs(dL+dN-omega/2)/(omega/2)),
            "radius_m":R,"density_kg_m3":RHO,"Vs_m_s":VS,"Vp_m_s":VP,"radial":radial,"I_kg_m2":float(I)}

def analytic_profile(radius, parameter, reference):
    r=np.asarray(radius,dtype=float); x=reference['x']; omega=reference['omega']; I=reference['I_kg_m2']
    W=spherical_jn(2,x*r/R); dW=(x/R)*spherical_jn(2,x*r/R,derivative=True)
    if parameter=='L': value=r**2*(dW-W/np.where(r==0,1,r))**2/(2*omega*I)
    else: value=4*W**2/(2*omega*I)
    value=np.where(r==0,0,value)
    return value*RHO*VS**2

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',type=Path,required=True); ap.add_argument('--build-root',type=Path,default=Path('/tmp/modal-band-build'))
    ap.add_argument('--dspec',type=Path,default=DSPEC); ap.add_argument('--yspec',type=Path,default=YSPEC)
    ap.add_argument('--mode',choices=('toroidal','spheroidal'),default='toroidal'); args=ap.parse_args()
    if args.mode=='spheroidal':
        from spheroidal_benchmark import run_spheroidal
        return run_spheroidal(args)
    out=args.output.resolve(); build=args.build_root.resolve(); output_guard(out)
    deps_record={name:pinned_dependency((DEPS/name if name.endswith('-src') else
        DEPS/'planetmodel-src' if name=='planetmodel-dspec' else PLANET),commit)
        for name,commit in DEP_PINS.items()}
    env=env_bounded(); model=HERE/'homogeneous_sphere.deck'; maxstep=.02; nq=8; tol=1e-7; mesh=1000
    dsrc=build/'dspec'; ysrc=build/'yspec'
    for solver,bdir,source,planet in [('dspec',dsrc,args.dspec,DEPS/'planetmodel-src'),('yspecpp',ysrc,args.yspec,PLANET)]:
        cm=['cmake','-S',HERE,'-B',bdir,'-DSOLVER='+('dspec' if solver=='dspec' else 'yspecpp'),
            '-DDSPEC_SOURCE_DIR='+str(args.dspec),'-DYSPEC_SOURCE_DIR='+str(args.yspec),'-DDEPENDENCY_ROOT='+str(DEPS),
            '-DPLANETMODEL_SOURCE_DIR='+str(planet),'-DCMAKE_BUILD_TYPE=Release']
        run(cm,env); run(['cmake','--build',bdir,'--target','kernel_adapter','-j1'],env)
    exeD=dsrc/'kernel_adapter'; exeY=ysrc/'kernel_adapter'; fixture=out/'homogeneous_sphere.deck'; shutil.copy2(model,fixture)
    # Record the predeclared protocol before producing any solver comparison output.
    protocol={"fixture":"homogeneous isotropic elastic sphere","geometry":{"source_lat_lon_deg":[0,0],"source_depth_m":50000,"moment_tensor_Nm":[0,1e19,0,0,0,0],"receiver_lat_lon_deg":[0,90],"receiver_depth_m":3000,"component":"North","degree":2},
      "spectrum":{"band_mHz":list(BAND),"duration_s":DURATION,"dt_s":DT,"nt":NT,"df_Hz":DF,"damping_s-1":EPS,"common_bins_inclusive":[FIRST,END-1],"window":"Hann, 4321 samples"},
      "kernel":{"h":H,"perturbation":"uniform fractional L=N=mu; A,C,F,rho fixed","linear_frequency_shift_rad_s":H*0.5*brentq(lambda v:v*spherical_jn(2,v,derivative=True)-spherical_jn(2,v),2.4,2.6)*VS/R,
        "phase_shift_over_record_rad":H*0.5*brentq(lambda v:v*spherical_jn(2,v,derivative=True)-spherical_jn(2,v),2.4,2.6)*VS/R*DURATION,
        "centered_shape_curvature_estimate":(H*0.5*brentq(lambda v:v*spherical_jn(2,v,derivative=True)-spherical_jn(2,v),2.4,2.6)*VS/R*DURATION)**2/6,
        "fit":"real-linear processed complex pole plus linear complex background after common 1/(i z); six nuisance coordinates plus fitted omega"},
      "gates":GATE,"status":"predeclared"}
    (out/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    # DSpec family preflight: toroidal signal and spheroidal cancellation at this geometry.
    dirs={}
    def dspec_run(label,typ,delta=0):
        dest=out/label; run([exeD,'spectrum',fixture,dest,maxstep,nq,*BAND,delta,typ,out/'empty_q.csv'],env); return dest
    (out/'empty_q.csv').write_text('basis,bin,real,imag\n')
    dirs['dspec_t']=dspec_run('dspec_baseline',2)
    dirs['dspec_s']=dspec_run('dspec_spheroidal_preflight',3)
    dirs['yspec']=out/'yspec_baseline'; run([exeY,'spectrum',fixture,dirs['yspec'],*BAND,0,tol,out/'empty_q.csv',mesh],env)
    ds_raw,ds_proc,ds_rows=read_spectrum(dirs['dspec_t']/'spectrum.csv')
    s_raw,s_proc,_=read_spectrum(dirs['dspec_s']/'spectrum.csv')
    y_raw,y_proc,y_rows=read_spectrum(dirs['yspec']/'spectrum.csv')
    idx=BINS
    s_leak=float(np.linalg.norm(s_proc[idx])/np.linalg.norm(ds_proc[idx]))
    raw_l2=float(np.linalg.norm(ds_raw[idx]-y_raw[idx])/np.linalg.norm(y_raw[idx]))
    proc_l2=float(np.linalg.norm(ds_proc[idx]-y_proc[idx])/np.linalg.norm(y_proc[idx]))
    metaD=json.loads((dirs['dspec_t']/'meta.json').read_text()); metaY=json.loads((dirs['yspec']/'meta.json').read_text())
    if s_leak>GATE['spheroidal_leakage']: raise SystemExit(f'spheroidal leakage gate failed: {s_leak}')
    if raw_l2>GATE['raw_spectrum_rel_l2'] or proc_l2>GATE['processed_spectrum_rel_l2']: raise SystemExit(f'baseline spectrum gate failed: raw={raw_l2}, processed={proc_l2}')
    fits={}
    for solver,proc in [('dspec',ds_proc),('yspecpp',y_proc)]:
        fits[solver]=fit_model(proc)
        write_q(out/f'q_{solver}.csv',fits[solver])
    # Centered actual forward perturbations. They are diagnostic, never substitutes for production kernels.
    kernel_runtime={}
    for solver,exe in [('dspec',exeD),('yspecpp',exeY)]:
        pair=[]
        for sign,label in [(1,'plus'),(-1,'minus')]:
            dest=out/f'{solver}_{label}';
            if solver=='dspec': cmd=[exe,'spectrum',fixture,dest,maxstep,nq,*BAND,sign*H,2,out/'empty_q.csv']
            else: cmd=[exe,'spectrum',fixture,dest,*BAND,sign*H,tol,out/'empty_q.csv',mesh]
            run(cmd,env); pair.append(read_spectrum(dest/'spectrum.csv')[1])
        # Production kernels are evaluated at baseline model; seven q rows share one Y solve.
        dest=out/f'{solver}_kernels'; qfile=out/f'q_{solver}.csv'
        if solver=='dspec': cmd=[exe,'kernel',fixture,dest,maxstep,nq,*BAND,0,2,qfile]
        else: cmd=[exe,'kernel',fixture,dest,*BAND,0,tol,qfile,mesh]
        start=time.monotonic(); run(cmd,env); kernel_runtime[solver]=time.monotonic()-start
    # Compare each kernel's seven-coordinate omega derivative with analytic and centered projections.
    ref=analytic_reference()
    kernel_results={}
    for solver in ('dspec','yspecpp'):
        rows=readcsv(out/f'{solver}_kernels/kernel_summary.csv')
        by={r['basis']:r for r in rows}
        dw=float(by['q6']['total'])
        dL=float(by['q6']['distributed_L'])+float(by['q6']['source_L'])
        dN=float(by['q6']['distributed_N'])+float(by['q6']['source_N'])
        kernel_results[solver]={"domega_dh":dw,"omega_over_2":ref['omega_over_2'],
           "relative_analytic_error":abs(dw-ref['omega_over_2'])/ref['omega_over_2'],
           "domega_duniform_L_fraction":dL,"domega_duniform_N_fraction":dN,
           "L_relative_analytic_error":abs(dL-ref['domega_duniform_L_fraction'])/abs(ref['domega_duniform_L_fraction']),
           "N_relative_analytic_error":abs(dN-ref['domega_duniform_N_fraction'])/abs(ref['domega_duniform_N_fraction']),
           "kernel_runtime_s":kernel_runtime[solver],
           "q_coordinate_derivatives":{name:float(row['total']) for name,row in by.items()},
           "q_omega_nuisance_orthogonality_max":fits[solver]['q_identity']}
    radial_curve_comparison={}
    for solver in ('dspec','yspecpp'):
        rows=readcsv(out/f'{solver}_kernels/kernel_q6.csv')
        radius=np.array([float(r['radius_m']) for r in rows])
        if solver=='dspec': weights=np.array([float(r['weight_m']) for r in rows])
        else:
            weights=np.array([float(r['weight_m3'])/max(float(r['radius_m'])**2,1e-300) for r in rows])
        radial_curve_comparison[solver]={}
        for param in ('L','N'):
            observed=np.array([float(r['radial_'+param]) for r in rows])
            analytic=analytic_profile(radius,param,ref)
            radial_curve_comparison[solver][param]={
                "weighted_L2":float(np.sqrt(np.sum((observed-analytic)**2*weights)/np.sum(analytic**2*weights))),
                "peak_scaled_max":float(np.max(np.abs(observed-analytic))/np.max(np.abs(analytic)))}
    # Full centered complex processed spectrum change vs its seven-dimensional kernel projection.
    finite={}
    for solver in ('dspec','yspecpp'):
        plus=read_spectrum(out/f'{solver}_plus/spectrum.csv')[1][idx]
        minus=read_spectrum(out/f'{solver}_minus/spectrum.csv')[1][idx]
        coords=[kernel_results[solver]['q_coordinate_derivatives'][f'q{i}'] for i in range(7)]
        finite[solver]=fit_change(fits[solver],plus,minus,coords)
        finite[solver]['finite_fit_omega_shift']=(fit_model(plus)['omega']-fit_model(minus)['omega'])/2
        finite[solver]['kernel_predicted_omega_shift']=H*kernel_results[solver]['domega_dh']
        finite[solver]['linear_omega_shift_expected']=H*ref['omega_over_2']
        finite[solver]['finite_exact_omega_shift_expected']=ref['omega']*(math.sqrt(1+H)-math.sqrt(1-H))/2
        d=np.asarray(fits[solver]['matrix'])@(H*np.asarray(coords)); n=BINS.size
        with (out/f'{solver}_spectrum_change.csv').open('w',newline='') as f:
            w=csv.writer(f); w.writerow(['bin','frequency_hz','centered_delta_real','centered_delta_imag','kernel_prediction_real','kernel_prediction_imag','residual_real','residual_imag'])
            delta=(plus-minus)/2; predicted=d[:n]+1j*d[n:]
            for j,k in enumerate(BINS): w.writerow([int(k),k*DF,delta[j].real,delta[j].imag,predicted[j].real,predicted[j].imag,(delta-predicted)[j].real,(delta-predicted)[j].imag])
    with (out/'dt_analytic_radial.csv').open('w',newline='') as f:
        w=csv.writer(f); w.writerow(['radius_m','domega_dL_per_m','domega_dN_per_m'])
        for values in zip(ref['radial']['radius_m'],ref['radial']['domega_dL_per_m'],ref['radial']['domega_dN_per_m']): w.writerow(values)
    metadata_errors={}
    for key in ('nt','dt_s','df_hz','damping_s','source_radius_m'):
        metadata_errors[key]=abs(float(metaD[key])-float(metaY[key]))/max(abs(float(metaD[key])),abs(float(metaY[key])),1e-300)
    metadata_errors['first_bin']=0.0 if int(metaD['first'])==int(metaY['first']) else math.inf
    metadata_errors['end_bin_exclusive']=0.0 if int(metaD['end'])==int(metaY['last'])+1 else math.inf
    metadata_max=max(metadata_errors.values())
    negative_control=float(np.linalg.norm(1.01*ds_proc[idx]-y_proc[idx])/np.linalg.norm(y_proc[idx]))
    # Verification gates and minimal self-contained JSON output.
    gate_results={
      "metadata_grid":{"value":metadata_max,"limit":GATE['metadata_rel'],"pass":metadata_max<=GATE['metadata_rel']},
      "spheroidal_leakage":{"value":s_leak,"limit":GATE['spheroidal_leakage'],"pass":s_leak<=GATE['spheroidal_leakage']},
      "baseline_raw_spectrum":{"value":raw_l2,"limit":GATE['raw_spectrum_rel_l2'],"pass":raw_l2<=GATE['raw_spectrum_rel_l2']},
      "baseline_processed_spectrum":{"value":proc_l2,"limit":GATE['processed_spectrum_rel_l2'],"pass":proc_l2<=GATE['processed_spectrum_rel_l2']},
      "corrupted_spectrum_negative_control":{"value":negative_control,"minimum_detectable":GATE['processed_spectrum_rel_l2'],"pass":negative_control>GATE['processed_spectrum_rel_l2']},
      "analytic_integral_identity":{"value":ref['quadrature_identity_rel'],"limit":GATE['analytic_omega_rel'],"pass":ref['quadrature_identity_rel']<=GATE['analytic_omega_rel']}}
    gate_results['cross_code_omega']={"value":abs(kernel_results['dspec']['domega_dh']-kernel_results['yspecpp']['domega_dh'])/ref['omega_over_2'],"limit":GATE['cross_code_omega_rel'],"pass":abs(kernel_results['dspec']['domega_dh']-kernel_results['yspecpp']['domega_dh'])/ref['omega_over_2']<=GATE['cross_code_omega_rel']}
    for param in ('L','N'):
        key=param+'_relative_analytic_error'
        value=max(kernel_results[s][key] for s in ('dspec','yspecpp'))
        gate_results['radial_'+param]={"value":value,"limit":GATE['radial_LN_rel'],"pass":value<=GATE['radial_LN_rel']}
        cross=(abs(kernel_results['dspec']['domega_duniform_'+param+'_fraction']-
                    kernel_results['yspecpp']['domega_duniform_'+param+'_fraction'])/
               abs(ref['domega_duniform_'+param+'_fraction']))
        gate_results['cross_code_radial_'+param]={"value":cross,"limit":GATE['radial_LN_rel'],"pass":cross<=GATE['radial_LN_rel']}
    for s in ('dspec','yspecpp'):
        gate_results[s+'_baseline_frequency']={"value":abs(fits[s]['fit_frequency_rel_bias']),"limit":GATE['baseline_frequency_rel'],"pass":abs(fits[s]['fit_frequency_rel_bias'])<=GATE['baseline_frequency_rel']}
        gate_results[s+'_omega_kernel']={"value":kernel_results[s]['relative_analytic_error'],"limit":GATE['analytic_omega_rel'],"pass":kernel_results[s]['relative_analytic_error']<=GATE['analytic_omega_rel']}
        gate_results[s+'_linear_prediction']={"value":finite[s]['residual_l2'],"limit":GATE['centered_linear_prediction_rel_l2'],"pass":finite[s]['residual_l2']<=GATE['centered_linear_prediction_rel_l2']}
        gate_results[s+'_q_orthogonality']={"value":fits[s]['q_identity'],"limit":GATE['sensitivity_identity'],"pass":fits[s]['q_identity']<=GATE['sensitivity_identity']}
    result={"protocol":protocol,"metadata":{"dspec":metaD,"yspecpp":metaY},"preflight":{"spheroidal_processed_l2_over_toroidal":s_leak,"raw_spectrum_relative_l2":raw_l2,"processed_spectrum_relative_l2":proc_l2},
      "analytic_reference":{k:v for k,v in ref.items() if k not in ('radial','I_kg_m2')},"metadata_relative_errors":metadata_errors,"radial_curve_comparison":radial_curve_comparison,"fits":{s:{k:v for k,v in f.items() if k not in ('matrix','qrows','fit_target')} for s,f in fits.items()},
      "kernel_results":kernel_results,"centered_finite_difference":finite,"gates":gate_results,
      "provenance":{"DSpecM1D":identity(args.dspec),"YSpecpp":identity(args.yspec),"pinned_dependencies":deps_record,"adapter_sha256":{"dspec_adapter.cpp":sha(HERE/'dspec_adapter.cpp'),"yspecpp_adapter.cpp":sha(HERE/'yspecpp_adapter.cpp'),"run_benchmark.py":sha(HERE/'run_benchmark.py'),"fixture":sha(model)}}}
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    write_report(out,result)
    passed=all(v['pass'] for v in gate_results.values())
    print(json.dumps({"status":"PASS" if passed else "FAIL","gates":gate_results},indent=2))
    if not passed: raise SystemExit(1)
if __name__=='__main__': main()
