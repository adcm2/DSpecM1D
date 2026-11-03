#!/usr/bin/env python3
"""Build the Phase A comparison from saved, verified receipts; no solver execution."""
from pathlib import Path
import csv, hashlib, json, shutil
A=Path(__file__).resolve().parents[1]; N=A.parent/'tunya-profile-20261027'
P=A/'receipts/post-analysis-completed-20261027T1834Z/aggregate-and-resolution/analysis'
R=A/'reports/acceleration-comparison-20261027'
paths={'layout':A/'analyses/layout-all-components.json','native_v4':A/'analyses/final-direct-all-components.json','native_v5':A/'analyses/native-final-v5-all-components.json','global':P/'global_v8.json','serial_global':P/'serial_global_preserved.json','paired':P/'paired_repeats.json','resolution':P/'resolution.json','post_summary':P/'summary.json','prior_yspec':N/'analyses/profile06-complete-material-aware-20261027T123112Z.json','independent_review':A/'review/phase-a-final-evidence-data.json'}
d={k:json.loads(p.read_text()) for k,p in paths.items()}
report={'scope':'Phase A acceleration, fixed physical case; no figures or requested high-degree full179sideglobalvolume','inputs':{k:{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for k,p in paths.items()},'data':d}
(A/'analyses/phase-a-final-comparison.json').write_text(json.dumps(report,indent=2)+'\n')
records=[]
for route in ['layout','native_v4','native_v5','global','serial_global']:
 for name,j in d[route]['jobs'].items():
  for stage,v in j.get('profile',{}).get('stages',{}).items():records.append([route,name,j.get('result',{}).get('degree'),stage,v['inclusive_s'],v['exclusive_s'],v['calls']])
for name,j in d['paired'].items():
 for stage,v in j['profile']['stages'].items():records.append(['paired',name,j['result']['degree'],stage,v['inclusive_s'],v['exclusive_s'],v['calls']])
for name,j in d['prior_yspec']['profiled_stage_records'].items():
 for stage,v in j['stages'].items():records.append(['prior_profile06',name,d['prior_yspec']['telemetry'].get(name,{}).get('actual_degree_from_result',''),stage,v['inclusive_s'],v['exclusive_s'],v['calls']])
with (A/'analyses/phase-a-all-profile-components.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['route','job','degree','stage','inclusive_s','exclusive_s','calls']);w.writerows(records)

def esc(s):
 return str(s).replace('\\',r'\textbackslash{}').replace('_',r'\_').replace('&',r'\&').replace('%',r'\%').replace('#',r'\#')
def fmt(x):return f'{x:.3f}' if abs(x)>=.001 else (f'{x:.2e}' if x else '0')
def table(headers,rows,layout=None):
 out=[r'\begin{center}',r'\small',r'\begin{tabular}{'+(layout or ('l'+'r'*(len(headers)-1)))+'}',r'\toprule',' & '.join(headers)+r' \\',r'\midrule']
 out += [' & '.join(map(str,row))+r' \\' for row in rows]
 return '\n'.join(out+[r'\bottomrule',r'\end{tabular}',r'\end{center}'])
def prof(route,job):return d[route]['jobs'][job]['profile']
def ex(route,job,key):return prof(route,job)['stages'].get(key,{}).get('exclusive_s',0)
tex=[r'''\documentclass[11pt,a4paper]{article}
\usepackage{alex-report}
\usepackage{longtable,array}
\ReportSetup{DSpecM1D / YSpecpp}{Acceleration phase A}{27 October 2026}
\hypersetup{pdftitle={DSpecM1D acceleration and YSpecpp comparison}}
\begin{document}
\ReportTitle{Sensitivity-kernel acceleration}{Timings, validation and remaining limits on Tunya}
The approved acceleration campaign is complete. The latest native candidate reduces the degree-1536 unprofiled runtime from 536.442 to 313.361 seconds (41.6\%), with byte-identical scientific outputs. The previous accelerated candidate remains faster at degrees 128 and 384. Global synthesis/storage/threading checks pass internally, but its frozen comparison against the original synthesis route still fails. No sensitivity-kernel figures or requested high-degree complete 179-side global volume were generated.
\section{Fixed experiment and comparison domains}
The observable is windowed vertical-displacement energy,
\[
E=\frac{\Delta t}{2}\sum_j H_{1000\text{--}1200}(t_j)\,u_z(t_j)^2.
\]
The PREM500 model uses a source at $10^\circ$N, $20^\circ$E, depth 50 km, and receiver at $45^\circ$N, $75^\circ$E, depth 3 km. The nominal band is 0--10 mHz, flat to 8 mHz and tapered to 10 mHz. Sampling is 5 s, record duration 3600 s, FFT length 1024 and period 5120 s; damping is $5/3600\ \mathrm{s}^{-1}$. There are 52 physical bins including DC, plus one zero-weight guard bin in response capture. The SEM uses six GLL points per element and radial refinement 0.05, with 179 positive native element/node sides. Both sides of interfaces remain distinct. Contacts are off; source-shell contact contributions are omitted.

At high degree the native DSpec/YSpec comparison uses the same \textbf{64-angle, 179-side diagnostic}, not a complete global volume. Global diagnostics use the complete product-bandwidth grid, with $n_\theta=2L+1$ and $n_\phi=\operatorname{FastFFTSize}(4L+1)$, on eight sides at low degree, two at $L=128,384$, and one fluid side at $L=1536$. Global timings must not be compared to YSpec's 64-angle timing as if their output domains matched.\footnote{Timings are individual observed runs on the same Tunya host, not statistical estimates or comparisons at matched local-kernel accuracy. Inclusive profiler scopes overlap; only exclusive entries may be summed. All stages at every measured degree are retained in the accompanying JSON and CSV.}
\section{Native timings and individually measured changes}
''']
rows=[]
for L in [128,384,1536]:
 old=d['layout']['jobs'][f'reference-l{L}']['process_wall_s']
 v4=d['paired'][f'paired-v4-l{L}']['process_wall_s'] if L<1536 else d['native_v4']['jobs'][f'final-symbolic-l{L}']['process_wall_s']
 v5=d['paired'][f'paired-v5-l{L}']['process_wall_s'] if L<1536 else d['native_v5']['jobs'][f'final-v5-l{L}']['process_wall_s']
 rows.append([L,fmt(old),fmt(v4),fmt(v5),f'{100*(1-v5/old):.1f}\\%'])
tex += [table(['$L$','Original (s)','Candidate v4 (s)','Candidate v5 (s)','v5 reduction'],rows),r'''All entries above are instrumented process wall times. Degrees 128 and 384 use the isolated paired reruns for v4/v5; degree 1536 uses the accepted full component runs. At 1536 the unprofiled v4/v5 runs take 369.677/313.361 s and agree byte-for-byte with their instrumented builds. Native outputs also match the original production route through a verified hash chain, despite the matched-response control having a different baseline capture setting.

At 128 and 384, v5 is respectively 11.0\% and 6.5\% slower than v4 in the isolated repeats. The regression is concentrated in dynamic assembly, scaling and numerical LU, with identical work counters. Its underlying compiler/allocation/cache cause has not been established. The low-degree result cannot be attributed solely to earlier overlapping work.

The retained changes are shared coefficient-label layouts/permutations, nonradial response reuse, exact-pattern symbolic reuse with fresh numerical LU, and removal of two unused full-order zero adjoint matrices. The radial $L=0$ weighted solve is retained to preserve the frozen numerical result. The separate native-node recovery prototype was rejected: at 1536 it took 414.657 s versus 400.374 s for packing alone, and increased the recovery component from 64.631 to 78.965 s.

\newpage
The last two-line allocation removal avoids about 2.15 TB of cumulative zero-fill traffic in the pole-only cached route at degree 1536. The source/solver exclusive remainders fall by 64.529 s; the observed net wall saving is 56.021 s. It is retained for its measured high-degree benefit, with its low-degree regression disclosed.
''']
keys=[('Coefficient packing/reorder','tensor_coefficient_packing_and_key_reorder'),('Layout deduplication','tensor_coefficient_layout_dedup_and_permutation_build'),('Dynamic matrix assembly','dynamic_H_assembly'),('Numerical LU','sparse_numerical_factorization'),('Symbolic analysis','sparse_symbolic_analysis'),('Equilibration/scaling','sparse_equilibration_and_scaling'),('Forward back-solves','sparse_forward_solve_including_rhs_scaling'),('Adjoint back-solves','sparse_adjoint_solve_including_rhs_scaling'),('Native strain/derivative recovery','radial_interpolation_and_strain_numerator_build'),('Angular Eigen multiplication','tensor_angular_gemm_complex_and_packing'),('Receiver basis/QR','receiver_basis_and_QR_reconstruction'),('Spheroidal source/solver remainder','spheroidal_source_lift_projection_and_solver'),('Toroidal source/solver remainder','toroidal_source_lift_projection_and_solver'),('Final pointwise products','pointwise_strain_rotation_kernel_contraction')]
rows=[[label,fmt(ex('layout','reference-l1536',k)),fmt(ex('native_v4','final-symbolic-l1536',k)),fmt(ex('native_v5','final-v5-l1536',k))] for label,k in keys]
tex += [table(['Exclusive component, $L=1536$','Original (s)','v4 (s)','v5 (s)'],rows,'p{8.2cm}rrr'),r'''The historical recovery timer still contains ``interpolation'' in its name. At native SEM nodes, the sampling plan selects the node and applies cached $D_\mathrm{ref}/J$ derivative weights to the current modal displacement fields. It does not require arbitrary off-grid radial output. Gathering modal values, applying these weights, constructing strains and allocating small coefficient records remain real work. The rejected stack-allocation prototype did not implement a successful fully batched derivative/recovery algorithm.
''']
r=d['native_v5']['jobs']['final-v5-l1536']['result'];c=prof('native_v5','final-v5-l1536')['counters']
tex += [f"At degree 1536, factorisations fall from 322665 to {c['sparse_factorization_calls']}; only {r['sparse_symbolic_analysis_calls']} symbolic analyses are performed, with {r['sparse_symbolic_analysis_reuses']} reuses. Nonradial response capture has {r['response_cache_misses']} misses and {r['response_cache_hits']} later hits; 105 radial bypass calls retain the original solve. Logical retained response storage is {r['response_cache_bytes']/2**30:.3f} GiB within the 32 GiB budget, and observed profiled-native peak RSS is 8.736 GiB. There is no cache fallback.\n"]
tex += [r'\section{YSpec comparison at the same native diagnostic points}']
rows=[]
for L in [128,384,1536]:
 t=d['prior_yspec']['telemetry'];rows.append([L,fmt(t[f'dsv6-native-l{L}']['wall_elapsed_s']),fmt(t[f'ysv6-native-wave-l{L}']['wall_elapsed_s']),fmt(t[f'ysv6-native-l{L}']['wall_elapsed_s'])])
tex += [table(['$L$','Original DSpec (s)','YSpec waveform (s)','YSpec kernel (s)'],rows),r'''These are the saved unprofiled profile-06 results. YSpec's independent waveform run is separate; its kernel diagnostic uses the frozen shared DSpec measurement baseline. Their times are not added and presented as a measured end-to-end YSpec run. DSpec's kernel timing includes its own baseline. The latest unprofiled DSpec degree-1536 result is 313.361 s, compared with the saved YSpec kernel diagnostic's 3338.354 s, subject to these methodology and scientific limitations.
''']
y=d['prior_yspec']['profiled_stage_records']['ysv6-profiled-l1536'];rows=[]
for label,key in [('Upward integration','radial.initial_conditions_and_upward_integration'),('Downward integration','radial.initial_conditions_and_downward_integration'),('Spheroidal b-vector transport','radial.spheroidal_b_vector_transport'),('Matching/response transport','radial.determinant_jumps_matching_and_response_transport'),('Coefficient/strain recovery','radial.coefficient_combine_and_strain_recovery'),('Source angular synthesis','kernel.angular_harmonic_synthesis.source'),('Adjoint angular synthesis','kernel.angular_harmonic_synthesis.adjoint'),('Elastic contraction','kernel.elastic_energy_pair_contraction'),('Parameter/objective products','kernel.objective_weight_and_parameter_products')]:
 if key in y['stages']:rows.append([label,fmt(y['stages'][key]['exclusive_s'])])
tex += [table(['YSpec exclusive component, $L=1536$','Seconds'],rows,'p{11cm}r'),r'''YSpec's cost is dominated by radial integration and b-vector transport, while the accelerated DSpec cost is shared between sparse solves, native strain recovery and angular synthesis. Final contractions are small in both native diagnostics. This is a comparison of the configured algorithms and tolerances, not proof of equal local-kernel accuracy.
\section{Global synthesis: complete angular grids on selected sides}
''']
rows=[]
for L in [128,384,1536]:
 for mode in ['generated','stored']:
  j=d['global']['jobs'][f'global-{mode}-t32-l{L}'];rows.append([L,mode,2 if L<1536 else 1,fmt(j['process_wall_s']),f"{j['peak_rss_gib_userspace']:.3f}"])
tex += [table(['$L$','Wigner mode','Sides','Wall (s)','Peak RSS (GiB)'],rows,'r l r r r'),r'''All runs request one outer worker and 32 inner GSH threads, with other numerical-library thread limits controlled. The maximum observed resident count is 34 including the supervisor. The source releases the main-thread OpenMP pool after grid construction and before the worker team is formed; runtime cleanup succeeds. The global cap is 40 and the waiting follow-up cap is 4, keeping the configured aggregate below 50. Supervision is userspace observation/enforcement, not kernel aggregate accounting.

Threading is harmful on the tiny cases: Generated degree-3/12 runs rise from 3.502/5.242 s serial to 21.939/27.510 s with 32 threads. Keep those cases serial. At degree 384, the serial Generated/Stored runs took 259.723/173.746 s; the 32-thread versions take 104.868/100.433 s. There is no completed serial degree-1536 reference: its earlier unfinished run was deliberately stopped and preserved at 1952.070 s. Do not infer a 32-fold speedup or turn that partial duration into a completed benchmark.

The degree-1536 Stored table is about 270.44 GiB; observed total peak RSS is 276.787 GiB, versus 6.860 GiB Generated. Stored saves only 25.565 s (5.0\%) in this one-fluid-side end-to-end case. Generated is the sensible default for the initial bounded outputs; many-side cache amortisation remains a separate measured decision. Both modes spool 11.293 GiB of modal records for all native sides, while the one-side kernel output is 0.846 GiB.
''']
rows=[]
for mode in ['generated','stored']:
 j=d['global']['jobs'][f'global-{mode}-t32-l1536'];t=j['global_manifest']['timings'];s=t['worker_stage_elapsed_sums'];worker=t['worker_elapsed_sum_s'];setup=t['serial_grid_and_context_setup_wall_s'];total=t['total_wall_s'];ww=t['worker_execution_wall_s']
 vals=[setup,s['spool_reads'],s['coefficient_rotation_layout'],s['gsh_inverse_fft'],s['point_kernel_contractions'],s['frequency_accumulation'],s['output_writes'],worker-sum(s.values()),ww-worker,total-setup-ww,total,j['process_wall_s']-total,j['process_wall_s']]
 if mode=='generated':gvals=vals
 else:svals=vals
labels=['Grid/context setup','Spool reads','Coefficient rotation/layout','Inverse GSH/FFT','Pointwise contractions','Frequency accumulation','Output writes','Worker residual','Worker scheduling/wall difference','Consumer preflight/cleanup/other residual','Total consumer','Process outside consumer','Process wall']
tex += [table(['Degree-1536 wall decomposition','Generated (s)','Stored (s)'],[[k,fmt(g),fmt(s)] for k,g,s in zip(labels,gvals,svals)],'p{10cm}rr'),r'''The manual consumer timings are a separate timing layer from MATHPROFILE and must not be added to its already-inclusive consumer coverage. Residuals are explicit: roughly 4.5--4.7 s outside setup/worker intervals and 1.3 s inside each worker remain unclassified. MATHPROFILE root exclusive sums close, but that does not imply every manual sub-operation is individually instrumented. The GSH backend here is the library's Loop plus FFT implementation, not a new external-BLAS GEMM; the native tensor route uses Eigen internal multiplication.

The consumer accumulates matching-frequency products side-by-side and writes each completed side once. It skips the exact-zero guard mixing and ten of twelve structurally zero fluid transforms. It does not sum the two wavefields independently over frequency before multiplying them: that would change the kernel. Dense cross-degree precontraction instead introduces $O(L^2)$ coupling and loses the low-rank frequency advantage when $L$ greatly exceeds the 52 active bins. Approximate compression would require a new error budget.
\section{Validation, unresolved discrepancies and resolution}
All 17 final-native jobs and 11 native comparison reports pass; all 19 global jobs and nine internal exact comparisons pass; the post-analysis queue's four additional paired exact checks pass. Native and global scientific-array comparisons are byte-exact where reported, work counters are unchanged by the final allocation fix, recorded inputs remain unchanged and process cleanup is verified. At degree 1536 the energy is $8.5734206967\times10^{-10}\ \mathrm{m}^2\mathrm{s}$ and the $q\cdot p$ versus $2E$ relative residual is $6.03\times10^{-16}$.

\textbf{Global scientific promotion remains held.} The original-direct versus global comparison at degrees 3 and 12 exceeds the frozen $2\times10^{-13}$ threshold in the L and fluid-bulk channels (roughly $3.1$--$4.3\times10^{-13}$ scaled differences). Internal threading/storage parity does not clear this separate gate, and no tolerance was relaxed. These differences are distinct from the inherited DSpec/YSpec discrepancies: at 1536 the saved solid-channel scaled differences include about 42.8\% in C, 34.1\% in F and 15.6\% in L on the diagnostic points. The fixed-degree waveform agreement is much better (relative $L^2\approx4.64\times10^{-5}$, energy difference $\approx2.68\times10^{-5}$). Neither proves local-kernel agreement.\footnote{The high-degree native samples occupy only the first 64 angular points of the relevant rule; these discrepancies are descriptive, not a representative global error norm. Native/interface identities are matched without radial interpolation.}

The bounded resolution comparison uses a 720-point source-to-receiver great-circle arc, a 5-degree angular exclusion at each endpoint, and approximate bilinear angular interpolation on two native sides. There is no radial interpolation. Between degrees 128 and 384, max-scaled differences are 8.22\% for the selected fluid bulk trace and 8.29--17.65\% across active solid channels; norm-scaled differences are 8.66\% and 8.93--16.63\%. This does not establish convergence. The earlier degree-384 to 1536 energy difference is about 2.93\% relative to the degree-1536 value (2.84\% relative to degree 384), so those degrees must remain labelled exploratory.
\section{Decision and next approved boundary}
Use the validated native route for initial exploratory sections and depth slices, retaining raw values, one-sided interfaces, contact omissions and explicit scales. At low degree v4 is currently faster; v5 has the demonstrated large-degree gain. The global route remains experimental until its original-route numerical gate is resolved. Its complete 179-side, degree-1536 output would occupy approximately 151.47 GiB; that requested high-degree full volume has not been materialised in this phase.

The proposed figure configuration keeps the same source, receiver, window and band, with kernel-bandwidth-derived section sampling ($4L+1$ points around the full great circle) and selected native depths near 96.446, 600 and 1019.368 km. Degree 384/1536 panels must be exploratory. Receiver lifts, Sobolev/Riesz gradients and receiver-point peak convergence remain deferred, as requested; they are not prerequisites for the initial figures.

After figures, the largest remaining native targets are genuine batched SEM recovery, angular synthesis and numerical LU. A vertical-receiver adjoint specialisation is a plausible further reduction, but has not been implemented or validated here. The low-degree assembly/scaling/LU regression needs a bounded compiler/memory investigation if further acceleration is requested. No additional optimisation campaign or figure generation starts automatically.
\clearpage
\appendix
\section{Complete degree-1536 native exclusive component table}
All entries below are exclusive seconds; zero denotes no entry under that scope name, not proof that corresponding work costs zero. Historical scope names are retained. Calls are for v5. Inclusive timings and all other measured degrees remain in the accompanying machine-readable files.
\begingroup\scriptsize
\begin{longtable}{>{\raggedright\arraybackslash}p{10cm}rrrr}
\toprule Scope & Original & v4 & v5 & Calls\\\midrule\endhead
''']
pa=prof('layout','reference-l1536')['stages'];pb=prof('native_v4','final-symbolic-l1536')['stages'];pc=prof('native_v5','final-v5-l1536')['stages']
for key in sorted(set(pa)|set(pb)|set(pc)):
 tex.append(esc(key)+' & '+' & '.join(fmt(s.get(key,{}).get('exclusive_s',0)) for s in [pa,pb,pc])+' & '+str(pc.get(key,{}).get('calls','--'))+r' \\')
tex += [r'\bottomrule\end{longtable}\endgroup',r'\clearpage\section{Complete degree-1536 YSpec kernel components}',r'\begingroup\scriptsize\begin{longtable}{>{\raggedright\arraybackslash}p{12cm}rr}\toprule Scope & Exclusive (s) & Calls\\\midrule\endhead']
for key,v in sorted(y['stages'].items()):tex.append(esc(key)+' & '+fmt(v['exclusive_s'])+' & '+str(v['calls'])+r' \\')
tex += [r'''\bottomrule\end{longtable}\endgroup
\section{Saved evidence and provenance}
The report data are derived from saved receipts, not rerun solvers. The independent final review qualifies execution/internal parity as passing and retains the separate global scientific failure. Its local rehash covers the collected small files; large kernel/spool payloads remain on Tunya, bound through verified remote artifact digests and comparison receipts. No production research checkout, existing evidence, commit or remote branch was changed by the acceleration phase.

Companion files: \texttt{analyses/phase-a-final-comparison.json} contains all source datasets and input SHA-256 values; \texttt{analyses/phase-a-all-profile-components.csv} retains every inclusive/exclusive profiler row at all measured degrees. \texttt{review/phase-a-final-evidence-review.md} states the independent audit limits. Proposed figure settings are in \texttt{configs/proposed-figure-case.json}. All paths are relative to the acceleration investigation folder.
\end{document}
''']
(R/'main.tex').write_text('\n'.join(tex))
shutil.copy2('/home/adcm2/.codex/skills/standard-latex-report/assets/alex-report.sty',R/'alex-report.sty')
print(json.dumps({'profile_rows':len(records),'report_source':str(R/'main.tex'),'comparison_json':str(A/'analyses/phase-a-final-comparison.json')}))
