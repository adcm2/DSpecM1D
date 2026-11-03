from pathlib import Path
import json,statistics
R=Path(__file__).resolve().parent;rows=json.loads((R/'final-timings-v6/summary.json').read_text());spectral=json.loads((R/'spectral/isolated-repetitions-v2/summary.json').read_text());ab=json.loads((R/'final-ablation-v5/summary.json').read_text());gates=json.loads((R/'final-timings-v6/frozen_comparisons.json').read_text());assert all(x['passed'] for x in gates)
def get(l,nf,m,k):return next(x for x in rows if (x['degree'],x['bins'],x['objectives'],x['kind'])==(l,nf,m,k))
macros={}
tr=[]
for l,nf,m in [(3,3,1),(6,3,1),(12,3,1),(12,3,2),(12,12,1),(12,48,1)]:
 o,f,y=[get(l,nf,m,k)['median_s'] for k in ['old','fast','yspec']];tr.append(f'{l} & {nf} & {m} & {o:.3f} & {f:.3f} & {y:.3f} & {o/f:.2f}$\\times$'+r'\\')
macros['HeadlineRows']='\n'.join(tr)
o,f,y=[get(12,3,1,k)['median_s'] for k in ['old','fast','yspec']];o2,f2,y2=[get(12,3,2,k)['median_s'] for k in ['old','fast','yspec']]
macros['HeadlineConclusion']=f'The degree-12 three-bin prototype improves the accepted DSpec calculation by {o/f:.2f} times for one objective and {o2/f2:.2f} times for two. Its measured times are {f:.3f}/{f2:.3f} s versus YSpec\'s {y:.3f}/{y2:.3f} s. The smaller-degree case remains slower than YSpec. These are measured single-thread prototypes, not changes to the production default.'
tr=[]
for case in dict.fromkeys(x['case'] for x in spectral):
 aa=[x['aggregate'] for x in spectral if x['case']==case];l=aa[0]['degree'];nf=aa[0]['frequency_count'];med=lambda fn:statistics.median(fn(x) for x in aa)
 vals=[med(lambda x:x['direct_synthesis_seconds_sum']),med(lambda x:x['covariance_formation_seconds_sum']+x['covariance_evaluation_seconds_sum']),med(lambda x:x['recoupling_seconds_sum']+x['kernel_assembly_seconds_sum']+x['fft_synthesis_seconds_sum']),med(lambda x:x['global_plan_prepare_seconds']),med(lambda x:x['svd_seconds_sum'])]
 tr.append(f'{l} & {nf} & '+' & '.join(f'{1000*v:.2f}' for v in vals)+r'\\')
macros['SpectralRows']='\n'.join(tr)
macros['SpectralPlanTwelve']=f"{statistics.median(x['aggregate']['global_plan_prepare_seconds'] for x in spectral if x['case']=='l12-f48'):.2f}"
tr=[]
for x in ab[:9]:tr.append(x['variant'].replace('-',' ')+' & '+f"{x['median_api_s']:.3f}"+' & '+('pass' if x['passed_volume_frozen'] else 'fail')+r'\\')
macros['AblationText']=r'''\begin{center}
\begin{tabular}{lrl}
\toprule
Actual L12/F3/full-grid route & Median (s) & Frozen volume gates\\
\midrule
'''+ '\n'.join(tr)+r'''
\bottomrule
\end{tabular}
\end{center}
The complex tensor layout is faster here than the real-degree split, despite its higher nominal arithmetic count. Small packed Eigen/OpenBLAS calls lose to the packed scalar route; matrix-matrix work over many radii pays. The nodal r-strain bridge loses on this grid. Strict native-architecture compilation (without fast-math, FMA contraction disabled) fails the processed-spectrum ($1.42\times10^{-12}$) and objective ($6.46\times10^{-12}$) gates and is rejected. The final headline table therefore uses portable O3 and the complex tensor layout.

The v6 forward angular cache adds actual reuse across objectives: it retains the volume forward strains and receiver displacement by exact grid/frequency identity. Hits skip source tensor synthesis and receiver suffix sampling; workspace preparation and representative coefficient gathering still run. The separate radial response cache saves repeated solves. At L12/F3 the forward field payload is 21.8 MB. Zero-budget fallback and contact-bearing reuse with radial response caching disabled both pass frozen gates. This scheduling/cache prototype is distinct from a general multi-objective API.
'''
checks=[v for g in gates for ob in g['objectives'] for k,v in ob['checks'].items() if k.startswith('kernel_')];maxl=max(x['relative_L2'] for x in checks);maxp=max(x['max_scaled_pointwise'] for x in checks)
macros['ValidationText']=r'''The final v6 headline outputs pass all 18 independent process/configuration comparisons against the unchanged production evaluator: returned raw/processed spectra, objectives and all six stored modulus channels, at exact matched coordinates. Repeated outputs are deterministic. Maximum kernel relative L2 is '''+('$'+f'{maxl:.2e}'.split('e')[0]+'\\times10^{'+str(int(f'{maxl:.2e}'.split('e')[1]))+'}$')+r''' and scaled pointwise '''+('$'+f'{maxp:.2e}'.split('e')[0]+'\\times10^{'+str(int(f'{maxp:.2e}'.split('e')[1]))+'}$')+r'''. Expanded tests separately establish 336 public field and 240 source-contact comparisons per bridge setting against the saved production arrays. An additional complete permuted PREM tensor product agrees with the standard field reconstruction in solid strains/fluid divergence below $1.13\times10^{-15}$; it is not the source of those historical contact arrays.

Fresh v5 optimized tensor Born checks pass the independent Cartesian regular and source-local volume/contact/total gates, with the existing factorization and radial-response cache disabled. The v6 forward angular cache is checked separately against frozen full outputs; it was not used in that Born adapter run. The final forward-cache tests pass all eight frozen comparisons at default and zero budgets, covering two/three objectives and complex/real layouts. Specifically, contact-bearing L3/three-objective and L12/two-objective cases use the complex layout and radial-response cache off; L12/three-objective complex/real cases use radial-response cache on and contacts off. The response-cache contacts-on and native L12 failures remain explicitly rejected; all gates were retained.
'''
mem=[]
for k,label in [('old','Phase 1'),('fast','prototype'),('yspec','YSpec')]:
 a=get(12,3,1,k);b=get(12,3,2,k);mem.append(f"{label}: {a['peak_process_hwm_bytes']/2**20:.1f}/{b['peak_process_hwm_bytes']/2**20:.1f} MiB")
macros['MemoryText']='For degree 12, one/two objectives, maximum sampled process high-water memory is '+ '; '.join(mem)+r'''. These values are polled after executable entry every 10 ms and include setup, outputs and, in the prototype driver, its subsequent same-engine reference calculation. They are process-level diagnostics, not isolated API allocator peaks or aggregate limits. The runtime advantage has a clear memory cost. The tensor coefficient multiplication is tiled, but full volume strain/output arrays and plan maps are still materialized; production integration should bound those intermediates more tightly.
'''
(R/'report/numbers.tex').write_text('\n'.join('\\newcommand{\\'+k+'}{'+v+'}' for k,v in macros.items())+'\n')
print('Report numbers generated',flush=True)
