from pathlib import Path
import json,shutil
R=Path(__file__).resolve().parent; M=R.parent
S=json.loads((M/'summary.json').read_text())
def get(solver,l=12,f=3,m=1,c=0):return next(s for s in S if (s['solver'],s['degree'],s['bins'],s['objectives'],s['contacts'])==(solver,l,f,m,c))
def ms(x):return f'{1000*x:,.2f}'
def esc(x):return x.replace('_',r'\_')
Dlabels={
'measurement':'Measurement adjoint $\\mathcal P^*$','operator':'Dynamic $H_\\ell(z)$ assembly','equilibration':'Matrix equilibration and scaling','symbolic':'Sparse symbolic analysis','factorization':'Sparse numerical factorization','factor_overhead':'Factor dispatch and wrapper overhead','forward_substitution':'Forward substitutions, including scaling','adjoint_substitution':'Hermitian-adjoint substitutions, including scaling','source_receiver_recovery':'Source lifts, receiver QR, modal recovery','response_cache':'Radial-response cache lookup and miss overhead','radial_strain':'Radial interpolation and strain numerators','representative_receiver_sampling':'Representative/receiver sampling and record capture','grid_key_lookup':'Full-volume semantic key construction/lookup','grid_cold_layout':'Cold tensor grid layout and maps','angular_plan_basis':'Angular geometry, Wigner tables and basis plans','result_allocation':'Volume result allocation and material tagging','angular_gemm':'Actual complex BLAS matrix products','gemm_buffers_copies':'BLAS layout conversion, buffers and copies','coefficient_packing':'Radial coefficient keys, packing and tile scatter','rotation_scatter':'Volume spin rotations and scatter','l0_scatter':'Degree-zero scatter','forward_cache':'Forward-field cache lookup, store and copy','receiver_scatter':'Receiver spectra scatter','kernel_products':'Pointwise modulus products and scatter','source_contacts':'Contact branch (disabled in main comparison)','Pforward':'Forward processing $\\mathcal P$, FFT/window','objective':'Final objective inner product','tensor_misc_cleanup':'Tensor bookkeeping, allocation and cleanup remainder','api_misc_cleanup':'API bookkeeping, copies and cleanup remainder'}
Ylabels={
'setup_api':'API geometry, sampling and measurement setup','measurement':'Measurement/basis frequency weights','ode_up_integration':'Initial conditions and upward ODE integration','ode_down_integration':'Initial conditions and downward ODE integration','b_vector_transport':'Spheroidal source-jump $b$-vector transport','matching_transport':'Determinants, jumps, matching and transport','radial_recovery':'Physical radial variable recovery','radial_strain':'Combine source/receiver coefficients and strains','radial_driver_overhead':'Radial context and family driver overhead','source_receiver_weights':'Forward source and receiver unit-force weights','angular_tables':'Wigner and generalized-Legendre tables','angular_source':'Source angular harmonic sums','angular_adjoint':'Adjoint angular harmonic sums','rotation':'Relative-frame spin rotations','energy_pairs':'Five-modulus strain-energy pair products','parameter_products':'Objective mixing and eight parameter products','modulus_weights':'Frequency/dispersion parameter weights','frequency_accumulation':'Block-to-output frequency accumulation','output_scaling':'Output transforms, scaling and construction','angular_block_misc_cleanup':'Angular-block allocation, loop and cleanup remainder','api_misc_cleanup':'API allocation, bookkeeping and cleanup remainder'}
def table_groups(s,labels):
 g=s['representative_profile']['groups'];total=sum(g.values())
 rows='\n'.join(f'{labels[k]} & {ms(g[k])} & {100*g[k]/total:.2f}\\\\' for k in labels if k in g)
 assert set(g)<=set(labels),set(g)-set(labels)
 return r'\begin{tabular}{p{.72\textwidth}rr}\toprule Component & ms & \%\\\midrule'+'\n'+rows+'\n'+r'\midrule Total instrumented API & '+ms(total)+r' & 100.00\\\bottomrule\end{tabular}'
def native_rows():
 rows=[]
 for l,f,m in [(3,3,1),(6,3,1),(12,3,1),(12,3,2),(12,12,1),(12,48,1)]:
  d=get('dspec',l,f,m);y=get('yspec',l,f,m)
  rows.append(f'{l} & {f} & {m} & {ms(d["native_api_median_s"])} & {ms(y["native_api_median_s"])} & {d["native_api_median_s"]/y["native_api_median_s"]:.3f}\\\\')
 return '\n'.join(rows)
d=get('dspec');y=get('yspec');d48=get('dspec',f=48);y48=get('yspec',f=48);dm=get('dspec',m=2);ym=get('yspec',m=2);dc=get('dspec',c=1)
pages=[]
pages.append(r'''\textbf{Finding.} The final modulus products are a small part of the cost. Both methods spend substantial time constructing radial responses and synthesizing angular fields. DSpec's tested tensor prototype moves much of its angular arithmetic into BLAS, but still pays for coefficient packing, copied BLAS buffers, map-based grid addressing, rotations and output movement. YSpec uses direct degree-contiguous angular sums and an ODE boundary-response method. Its fine phase timers are more intrusive, so native totals and instrumented components are reported separately.
\section{Native comparison and experiment}
\begin{center}\begin{tabular}{rrrrrr}\toprule
$L$ & bins & objectives & DSpec (ms) & YSpec (ms) & DSpec/YSpec\\\midrule
'''+native_rows()+r'''
\bottomrule\end{tabular}\end{center}
Each native entry is the median of three fresh processes, one numerical thread on one fixed allowed CPU, alternating solver/build order. These are complete numerical API times, including first-use plans and return cleanup, excluding frontend input/setup and export. DSpec sums its objective calls; YSpec handles objectives in one batch. Volume-only comparisons disable source contacts. A separate contact-enabled DSpec case follows.

The homogeneous shell fixture uses full gravity, no attenuation, SEM order parameter NN=8, the original 50 km-deep source and 3 km-deep receiver, and 56 one-sided element/radius traces. $L=3$ has 192 angles (10,752 points); $L=6,12$ have 1,352 angles (75,712 points). The selected direct tensor route retains the original arbitrary-radius compatibility grid; this is not a test of making native SEM output the default. Maximum frequencies for 3, 12 and 48 bins are respectively 1.171875, 4.6875 and 18.75 mHz.\footnote{Bins are $k/(128\times20\,\mathrm{s})$, $k=1,\ldots,F$, with damping $1/360\,\mathrm{s}^{-1}$. A bounded degree/frequency fixture is not a convergence result for a full Earth calculation.}

DSpec additionally returns raw/processed receiver spectra and the objective. YSpec's sampled API returns volume kernels and computes eight elastic/Q channels; only the five common elastic moduli are exported. These differing products are costed explicitly, rather than hidden behind an equal-work claim.
''')
pages.append(r'''\clearpage\section{What the mathematics asks each solver to do}
For each active frequency and spherical-harmonic family, DSpec forms
\begin{equation}
H_\ell(z)=K_\ell-z^2M_\ell,\qquad
H_\ell u_\ell=b_\ell,\qquad
H_\ell^{H}\lambda_\ell=R_\ell^{H}\widehat q,
\quad \widehat q=\mathcal P^*q.
\end{equation}
The prototype shares one equilibrated sparse LU factorization between forward and Hermitian-adjoint substitutions. ``One forward and one adjoint field'' does not mean two scalar solver calls: at $L=12$ each bin has 12 spheroidal, 12 toroidal and one radial operators, hence 75 distinct dynamic matrices over three bins. The measured one-objective route factors each once.

YSpec instead integrates radial ODE response systems,
\begin{equation}
\frac{dy}{dr}=A_\ell(r,z)y,
\end{equation}
constructing upward/downward boundary solutions, determinant and source-jump matching, spheroidal $b$-vector transport, and physical radial responses at requested nodes. Forward source and receiver unit-force responses share these fundamentals. It has no SEM matrix assembly or LU stage to time.

Both methods reconstruct six spin-strain components using radial coefficients:
\begin{equation}
\varepsilon_s(r,\Omega;z)=
\sum_{\ell,m} c_{s\ell m}(r,z)\,
 d^\ell_{s m}(\theta')e^{im\phi'},\qquad s=0,0,1,-1,2,-2.
\end{equation}
The tensor prototype evaluates these sums as $E=BC$ across 16-radius tiles, then rotates and scatters to output points. Angular geometry and basis values are frequency-independent and reused. YSpec builds tables in 64-angle blocks at each frequency, reusing each table over every radial node; its degree dots use complex coefficients and real Legendre values.

Finally, both form local bilinear elastic energy products and sum their correctly weighted real parts over frequency:
\begin{equation}
K_p(r,\Omega)=\sum_k w_k\,\Re\,B_p
\bigl(\varepsilon^F_k,\varepsilon^A_k\bigr).
\end{equation}
The dual/conjugation convention is carried by each implementation's coefficients. Summing each field over frequency before multiplication introduces unwanted cross-frequency products. Exact covariance contractions can reorder the frequency sum of coefficient \emph{products}; that alternative was explored in the previous audit and is not the direct tensor route profiled here.

The one-objective DSpec path still solves the full weighted global adjoint RHS plus a supplemental minimal receiver basis, even with contacts disabled. The two-objective cache path uses only the minimal unit-force basis and mixes weighted adjoints algebraically. Their different first-call costs must not be attributed solely to cache hits.
''')
pages.append(r'\clearpage\section{DSpec: complete three-bin component accounting}'+'\n'+f'At $L=12$, one objective, native median is {ms(d["native_api_median_s"])} ms. The representative instrumented run totals {ms(d["representative_profile"]["api_s"])} ms. Entries below are disjoint exclusive times; their unrounded sum closes the API root.\n'+r'\begin{center}\small'+'\n'+table_groups(d,Dlabels)+r'\end{center}'+'\n'+r'''Nested solve scopes provide a useful further split without being added again: full weighted adjoint solves and minimal supplemental basis solves are shown in the machine-readable stage table. Radial strain time includes interpolation-plan use and strain-numerator construction; it is not a pure interpolation-only timer. No interpolated raw-strain approximation is introduced.''')
pages.append(r'\clearpage\section{YSpec: complete three-bin component accounting}'+'\n'+f'At $L=12$, one objective, native median is {ms(y["native_api_median_s"])} ms. The representative instrumented run totals {ms(y["representative_profile"]["api_s"])} ms. These are instrumented durations, not a synthetic rescaling of native time.\n'+r'\begin{center}\small'+'\n'+table_groups(y,Ylabels)+r'\end{center}'+'\n'+r'''The energy-pair row computes the five bilinear modulus terms. The parameter-products row includes measurement/basis mixing and eight output channels, including three Q channels even for this non-attenuating calculation. It therefore covers more than five final multiplications. The angular-block remainder includes allocation, loop control, temporary energy-vector handling, cleanup and timer boundary overhead; it is retained explicitly.''')
# Rest of report generated after all results are ready, with exact values below.
(R/'findings.tex').write_text('\n'.join(pages)+'\n')
# Vector plot uses disjoint API-parent partitions, not sums of nested stages.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
fig,axes=plt.subplots(1,2,figsize=(10,3.2))
for ax,f in zip(axes,[3,48]):
 dd=get('dspec',f=f);yy=get('yspec',f=f);dp=dd['representative_profile'];yp=yy['coarse_representative_profile']
 dr=dp['stages']['frequency_workspace_and_modal_solve']['inclusive_s']
 df=dp['stages']['tensor_grid_semantic_plan_and_sampling']['inclusive_s']+dp['groups']['kernel_products']
 yr=sum(v for k,v in yp['groups'].items() if k.startswith(('ode_','radial_')) or k in ['b_vector_transport','matching_transport','source_receiver_weights'])
 yf=yp['stages']['kernel.angular_block_batch']['inclusive_s']
 vals=np.array([[dr,df,dp['api_s']-dr-df],[yr,yf,yp['api_s']-yr-yf]])
 left=np.zeros(2)
 for i,(label,color) in enumerate(zip(['Radial responses','Field reconstruction and products','Other API work'],['#344c64','#4e8c83','#c6aa75'])):
  ax.barh([0,1],vals[:,i],left=left,label=label,color=color);left+=vals[:,i]
 ax.set_yticks([0,1],['DSpec profile','YSpec coarse']);ax.invert_yaxis();ax.set_xlabel('Instrumented seconds');ax.set_title(f'{f} bins, degree 12');ax.grid(axis='x',alpha=.2);ax.set_axisbelow(True)
 ax.spines[['top','right']].set_visible(False)
handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3,frameon=False,fontsize=9);fig.tight_layout(rect=[0,.12,1,1]);fig.savefig(R/'stage-partitions.pdf');plt.close(fig)
# Scaling table exposes actual operations as well as time.
def keygroups(s,names):return sum(s['representative_profile']['groups'].get(k,0) for k in names)
def coarseval(s,key):return s['coarse_representative_profile']['groups'].get(key,0)
workrows=[]
for name,a,b in [
 ('DSpec numerical LU',d['representative_profile']['groups']['factorization'],d48['representative_profile']['groups']['factorization']),
 ('DSpec forward + adjoint substitutions',keygroups(d,['forward_substitution','adjoint_substitution']),keygroups(d48,['forward_substitution','adjoint_substitution'])),
 ('DSpec actual complex GEMM',d['representative_profile']['groups']['angular_gemm'],d48['representative_profile']['groups']['angular_gemm']),
 ('DSpec coefficient packing + BLAS buffers',keygroups(d,['coefficient_packing','gemm_buffers_copies']),keygroups(d48,['coefficient_packing','gemm_buffers_copies'])),
 ('DSpec volume rotations/scatter',d['representative_profile']['groups']['rotation_scatter'],d48['representative_profile']['groups']['rotation_scatter']),
 ('DSpec final modulus products',d['representative_profile']['groups']['kernel_products'],d48['representative_profile']['groups']['kernel_products']),
 ('YSpec up/down ODE integration (coarse)',coarseval(y,'ode_up_integration')+coarseval(y,'ode_down_integration'),coarseval(y48,'ode_up_integration')+coarseval(y48,'ode_down_integration')),
 ('YSpec source $b$ transport (coarse)',coarseval(y,'b_vector_transport'),coarseval(y48,'b_vector_transport')),
 ('YSpec whole angular/product block (coarse)',y['coarse_representative_profile']['stages']['kernel.angular_block_batch']['inclusive_s'],y48['coarse_representative_profile']['stages']['kernel.angular_block_batch']['inclusive_s'])]:workrows.append(f'{name} & {ms(a)} & {ms(b)}\\\\')
pages.append(r'''\clearpage\section{Why the wider frequency band separates the methods}
\begin{center}\includegraphics[width=\textwidth]{stage-partitions.pdf}\end{center}
\noindent\small The bars partition instrumented API roots exactly. DSpec's field block includes radial strain evaluation and the final point-kernel pass; YSpec forms radial strain coefficients in the radial block and fuses angular sums with products. These are comparable pipeline regions, not identical kernels. YSpec's coarse build has no per-point phase clocks or fine-child commits.\normalsize
\begin{center}\small\begin{tabular}{p{.68\textwidth}rr}\toprule
Measured work & 3 bins (ms) & 48 bins (ms)\\\midrule
'''+ '\n'.join(workrows)+r'''
\bottomrule\end{tabular}\end{center}
YSpec attempts 78,751 Runge--Kutta steps and evaluates 481,137 RHS callbacks over three bins. At 48 bins these rise to 3,998,468 and 24,129,174: about 51 and 50 times the work for 16 times as many bins. Its higher-frequency ODE integration/transport cost explains much of the growing gap. DSpec factors 75 versus 1,200 operators of the same bounded SEM sizes.

Actual DSpec complex-GEMM terms are 125,473,464 versus 2,007,575,424; YSpec performs 141,732,864 versus 2,267,725,824 complex-coefficient/real-Legendre degree terms. These are different operation types, not equal FLOP counts. Counting only $75{,}712\times6$ output strain slots misses the degree/order sums, two fields and frequency loop.
''')
# Reuse and non-API work, with measurement medians kept independent.
objp=dm['representative_per_objective'];contactstage=dc['representative_profile']['stages']['source_contact_separate_case']['inclusive_s']
frontendrows=[]
for label,sv1,sv2 in [
 ('Frontend setup, native median',d['setup_native_median_s'],y['setup_native_median_s']),
 ('Model/input to returned arrays, native median',d['model_to_arrays_native_median_s'],y['model_to_arrays_native_median_s']),
 ('Export, native median (binary versus CSV)',d['binary_export_native_median_s'],y['CSV_export_native_median_s'])]:frontendrows.append(f'{label} & {ms(sv1)} & {ms(sv2)}\\\\')
pages.append(r'''\clearpage\section{Repeated objectives, source contacts and frontend costs}
At $L=12$, three bins, DSpec's two-objective native total is '''+ms(dm['native_api_median_s'])+r''' ms, against YSpec's '''+ms(ym['native_api_median_s'])+r''' ms. DSpec's native median call durations are '''+ms(dm['native_per_objective_median_s'][0])+r''' and '''+ms(dm['native_per_objective_median_s'][1])+r''' ms. Their individual medians need not sum to the median of process totals.

The representative instrumented DSpec calls take '''+ms(objp[0]['api_s'])+r''' and '''+ms(objp[1]['api_s'])+r''' ms. The first factors all 75 operators; the second factors none, hits 75 radial-response entries and all three forward angular-field entries. It still performs weighted-adjoint reconstruction, source/modal bookkeeping, representative radial coefficient sampling, adjoint angular synthesis, point products and result allocation. Full-source angular synthesis and receiver-suffix sampling are skipped on the cache hit; workspace preparation and coefficient gathering are not.

The one-objective nonbatch path solves 219 forward RHS columns, 1,011 full weighted global-adjoint columns and 108 supplemental unit-basis columns. The two-objective batch path actually solves 219 forward plus 111 unit-basis columns once, then reconstructs the 2,022 requested weighted columns across both calls. Submitted cache RHS widths are not counts of substitutions on cache hits.

YSpec's two-objective batch retains the same 75 ODE response-family solves and the same 481,137 RHS evaluations. A second independent receiver direction adds coefficient combinations and angular/bilinear work rather than another fundamental radial integration. These are two different ways of sharing response information.

\textbf{Source contacts.} In a separate, nonbatch, contact-enabled DSpec run, the contact branch costs '''+ms(contactstage)+r''' ms across three bins. Native total is '''+ms(dc['native_api_median_s'])+r''' ms; the small difference from the contacts-off total is within the observed run spread, so subtracting the totals is not a valid contact timer. All six contact-enabled clone/output comparisons pass with zero contact error. YSpec's sampled API omits this source-shell term; no contact speed ratio is claimed.

\begin{center}\small\begin{tabular}{p{.68\textwidth}rr}\toprule
Three-bin boundary & DSpec (ms) & YSpec (ms)\\\midrule
'''+ '\n'.join(frontendrows)+r'''
\bottomrule\end{tabular}\end{center}
DSpec's representative setup includes '''+ms(d['setup_profile'][0]['stages']['sample_CSV_read']['exclusive_s'])+r''' ms sample-CSV parsing, '''+ms(d['setup_profile'][0]['stages']['SEM_mesh_and_matrix_setup']['exclusive_s'])+r''' ms SEM construction and '''+ms(d['setup_profile'][0]['stages']['result_stream_creation_and_coordinate_export']['exclusive_s'])+r''' ms stream/coordinate export. YSpec frontend setup combines model/configuration/sample preparation. Output formats differ; their export times are not solver times.

The DSpec diagnostic executable also runs a separate reference calculation after every candidate call: native medians are '''+ms(d['reference_benchmark_only_native_median_s'])+r''' ms (three bins) and '''+ms(d48['reference_benchmark_only_native_median_s'])+r''' ms (48 bins), plus roughly 3 ms validation and 11 ms binary export. This harness verification is excluded from solution-only API and model-to-array times. The benchmark wrapper's wall time additionally includes parent-side hashing/compression/comparison and must not be presented as child solver time.
''')
# Calibration, trustworthy accounting, and next work bounded by the measured scopes.
overrows=[]
for l,f,m in [(3,3,1),(6,3,1),(12,3,1),(12,3,2),(12,12,1),(12,48,1)]:
 dd=get('dspec',l,f,m);yy=get('yspec',l,f,m)
 overrows.append(f'{l}/{f}/{m} & {dd["profile_over_native_percent"]:+.1f} & {yy["coarse_over_native_percent"]:+.1f} & {yy["profile_over_native_percent"]:+.1f}\\\\')
pages.append(r'''\clearpage\section{Calibration, completeness and validation}
\begin{center}\small\begin{tabular}{lrrr}\toprule
$L$/bins/objectives & DSpec profile (\%) & YSpec coarse (\%) & YSpec fine (\%)\\\midrule
'''+ '\n'.join(overrows)+r'''
\bottomrule\end{tabular}\end{center}
Percentages compare instrumented median API time with the separate native median. They include hook perturbation and run-to-run variability; negative values do not mean instrumentation accelerates the solver. At the main $L=12$/three-bin point, native ranges are '''+ms(d['native_api_range_s'][0])+r'''--'''+ms(d['native_api_range_s'][1])+r''' ms for DSpec and '''+ms(y['native_api_range_s'][0])+r'''--'''+ms(y['native_api_range_s'][1])+r''' ms for YSpec.

YSpec fine timing accumulates raw phase intervals inside its original fused node/place/basis loop, with no new angular work buffers or changed operation order. It makes 2,725,632 fine timer reads at three bins, so small phase durations include appreciable timer boundary cost. Coarse timing compiles those phase hooks and child commits out while retaining outer scopes and counters. The difference is combined fine-hook perturbation, not a pure clock latency measurement. Its legacy \texttt{angular\_phase\_timer\_reads} counter counts disabled call sites, not actual clock reads. Neither detailed table should be relabeled as exact native per-component time, and no global rescaling is used.

All API components, including retained bookkeeping/cleanup remainders, are assigned once by exclusive time. For every root,
\begin{equation}
t_{\mathrm{API}}=\sum_s t_s^{\mathrm{exclusive}}.
\end{equation}
Closure errors are zero to floating-point summation roundoff (maximum below $10^{-8}$ s). Representative tables use the process with median instrumented API duration, so their unrounded rows close their own root; they are not sums of independently chosen component medians. Rounded table rows may differ in the last printed digit.

The final matrix comprises 84 processes, plus 18 coarse-control and seven contact processes: 109 final processes. All 96 clone-versus-original output gates pass; arrays are byte-identical in all these comparisons. DSpec also retains its original per-process reference comparison. Setup and export have independent accounting roots. Diagnostic packet export is disabled. No timings overlap compilation or another campaign numerical run; each process is pinned to one fixed allowed CPU with one numerical thread, a 180 s CPU limit and 3 GiB address-space limit. Hardware performance counters were unavailable, so these are wall-time scopes and exact software work counters, not PMU measurements.\footnote{Portable GCC Release/O3 without fast-math; original equilibrated DSpec solver and cached YSpec dependency library retained. Process address-space limits are not an aggregate kernel-enforced campaign memory budget.}

Frozen source review restored YSpec's original early-exit scan and DSpec's direct initialization forms before final runs. All 366 production DSpec tracked-file hashes and seven Phase0 artifact hashes are unchanged; its existing dirty state is preserved. YSpec remains clean at \texttt{280eab0}; DSpec remains at \texttt{a778b13}. All 14 protected original prototype files and 66 profiling source/build/runtime files verify against pre-run hashes. Instrumentation equivalence inherits the prototype's previous scientific checks; it adds no convergence claim or new whole-planet receiver accuracy claim.
''')
pages.append(r'''\clearpage\section{What the profile says to improve next}
\begin{enumerate}
\item \textbf{Remove avoidable data preparation in the direct tensor path.} Canonical coefficient indices and retained BLAS-ready bases can eliminate repeated map construction, tile reordering and basis copies. Current packing plus BLAS buffers cost '''+ms(keygroups(d,['coefficient_packing','gemm_buffers_copies']))+r''' ms at three bins and '''+ms(keygroups(d48,['coefficient_packing','gemm_buffers_copies']))+r''' ms at 48 bins, versus '''+ms(d['representative_profile']['groups']['angular_gemm'])+r''' and '''+ms(d48['representative_profile']['groups']['angular_gemm'])+r''' ms in actual matrix-product scopes. These are measured costs to target, not guaranteed removable savings.
\item \textbf{Make the volume-only minimal-adjoint path explicit.} The current one-objective branch computes both the full weighted global adjoint and a supplemental receiver basis. The nested full-weighted solve scope alone costs '''+ms(d['representative_profile']['stages']['full_weighted_HH_solve']['inclusive_s'])+r''' ms at three bins. The batch route already demonstrates basis reuse, but source contacts require their separate validated path and cache failure/admission rules still need integration work. Reusing symbolic sparsity analysis is another smaller target; numerical factors still change with frequency.
\item \textbf{Fuse field movement, relative rotation and energy products.} Volume rotations/scatter cost '''+ms(d['representative_profile']['groups']['rotation_scatter'])+r''' ms, result allocation/tagging '''+ms(d['representative_profile']['groups']['result_allocation'])+r''' ms, and separate degree-zero scatter '''+ms(d['representative_profile']['groups']['l0_scatter'])+r''' ms at three bins. YSpec shows the alternative common-frame fused product organization. Validate spin/dual signs, TI invariance, fluid behavior and source-contact separation before integrating this change.
\item \textbf{Fix reuse admission and grid addressing before chasing tiny products.} At 48 bins DSpec builds two tensor frame plans but records 104 basis misses, compared with 14 at three bins. Retained plan payload reaches 33,528,992 bytes of a 33,554,432-byte admission budget. The counters are consistent with two bases being rebuilt at each additional frequency because the remaining payload cannot admit them. Compact grid addressing and avoiding duplicated basis/storage payload should be evaluated before increasing the budget. The cold grid map is '''+ms(d['representative_profile']['groups']['grid_cold_layout'])+r''' ms; repeated full-volume semantic keys cost '''+ms(d48['representative_profile']['groups']['grid_key_lookup'])+r''' ms across 48 bins.
\item \textbf{Keep covariance and native SEM output as separate decisions.} Exact frequency covariance could remove repeated field synthesis and some buffers, but the current BLAS arithmetic is already fast and radial work remains substantial. A useful integrated covariance implementation must eliminate packing/movement as well as contractions. Native SEM output is a useful representation feature; the previous audit's strict high-degree native-node discrepancy must be resolved before adopting it as the default. Radial interpolation/strain construction here is only '''+ms(d['representative_profile']['groups']['radial_strain'])+r''' ms, so it does not explain the present runtime.
\end{enumerate}
These are proposed follow-on changes. This profiling follow-up changes only isolated instrumentation; it does not promote a new numerical route to production.

\section{Evidence and reproducibility}
The follow-up evidence root is \texttt{math-profile/}, beside the previous \texttt{algorithm-audit/} under this chat's \texttt{investigations/dspec-sk-profile-20261026/}. The report tables and vector plot are generated directly from \texttt{summary.json}. Exact native ranges, medians, selected run indices, every stage/counter and every gate remain in the saved JSON.
\begin{center}\small\begin{tabular}{p{.45\textwidth}p{.49\textwidth}}\toprule
Evidence & Location relative to the follow-up root\\\midrule
All 2,781 component timings & \texttt{all-stage-times.csv}\\
Exhaustive 105-stage mathematical map & \texttt{complete-component-map.json}\\
Native/fine full matrix and arrays & \texttt{runs-01/}\\
No-fine-hook YSpec control & \texttt{coarse-runs-01/}\\
Separate source-contact case & \texttt{contacts-01/}\\
Exact configuration/resource protocol & \texttt{protocol.json}\\
Instrumented source/build freeze & \texttt{final-before-runs.sha256.json}\\
Production/prototype preservation & \texttt{preservation-verification.json}\\
Fresh independent reviews & \texttt{review/}\\
\bottomrule\end{tabular}\end{center}
''')
(R/'findings.tex').write_text('\n'.join(pages)+'\n')
print('Generated',len(pages),'content pages')
