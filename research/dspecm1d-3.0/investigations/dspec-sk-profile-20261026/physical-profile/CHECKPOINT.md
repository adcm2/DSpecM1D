# DSpecM1D 2.0: shutdown / resume checkpoint

Saved 2026-10-26 22:11 UTC. This is a resumable checkpoint, not a completed physical full-kernel result.

## Authorized scope

Resolve explicit DC, waveform/energy convergence, then matched component profiling for PREM vertical-displacement windowed energy (Hann 1000–1200 s; raw-contour 8–10 mHz cosine taper). User explicitly requested **materialize the complete global kernel**, all positive native SEM radial sides and complete kernel-bandwidth angles. Partial profiling blocks do not fulfill materialization. Production stays read-only; no commit/push. Follow Luna High implementation → fresh Luna High review → short Sol High review and log every code change. No high-degree angular-library repair is authorized without renewed approval. Sources and earlier saved evidence are preserved.

## Confirmed results

- Explicit DC is at damped contour z=-i epsilon, epsilon=5/3600 per second, not undamped static zero. dt=5 s, N=1024, period5120 s, 52 bins including DC below10mHz. Qex2 doubles period/N at fixed epsilon.
- Dropping DC changes weighted waveform by1.4334% and energy0.8702% atL1536. Exclusion is a different observable.
- L1536 SEM n6→n8/n10: waveform~4.3e-5; halve radial mesh maxstep:1.03e-6. Time quadrature5→0.625s changes energy9.16e-7. Qex1→2 changes waveform2.647e-4, energy3.896e-4. These pass fixed1e-3 checks.
- Degree1536→1792 waveform1.23955e-3 fails; energy7.949e-4 passes. 1792→1827 waveform5.43067e-5 and energy3.47982e-5 pass, but only35degrees and1827 is current field-library safety ceiling. No broad higher-degree tail check exists; do not overclaim convergence.
- DSpec L1827 n8 vsYSpec L1827: waveform1.04547e-5, energy6.98193e-7. Energy~8.56667e-10 m²s. Independent SI spectrum reconstruction matches DSpec traces to~2e-16.
- Small complete-global L3 volume kernels: DSpec reference/pole/profile/contact-enabled byte-identical. YSpec v5 streamed output exactly equals v3 monolithic after layout transpose; native/profile output exact. YSpec v6 count64 probe equals first64angles of v5 full result.
- Global GSH exact rotation prototype v2: all12 complexchannels match tensor atL3/6/12 solid and L12 fluid/interface-side controls near roundoff. V1 fluid representation failure retained; v2 projects nonradial fluid rr+H into channel0 and preserves l0 separately. Luna+Sol reviews saved.
- Full global atL1827, n6:3655theta×7350phi×179radial sides=4,808,700,750points. Five-double result192.35GB; DSpec six-double230.82GB; YSpec scratch+result500.10GB. No such high-degree full file exists yet.
- Preservation check:438 protected files rehashed, zero changed/missing (preservation-check.json).

## Critical new algorithmic finding — validate next

Frozen DSpec v4 `Windowed3DFields.h::receiverAtPole` repeatedly calls `harmonic(l,m,n,0,0)`. Each call constructs full Wigner(degree,degree,nMax2) tables. Repeated over every degree, family, receiver and frequency this is O(F L³) setup even though pole values are exactly sqrt((2l+1)/(4pi))*delta_nm.

L384/64angle/all179radii profile: physical_case_total460.444s; frequency_workspace_prepare_and_modal_capture exclusive293.548s, numericalfactorization30.576s, radial interpolation/strain26.393s, tensor coefficientpacking41.259s, Eigencomplexsynthesis19.398s. Concurrent accuracy jobs make these component diagnostics, not controlled speed-ratio headlines.

Luna is compiling isolated **DSpec v5** from frozenv4 with ONLY closed-form receiverAtPole, corrected actual-small-m angular-cache byte estimate, and narrow setup/mixing timers. source-v5-frozen/, source manifest, freeze-v5-build.json exist; verify binary build completion/manifests. No v5 numerical gate yet. Preserve v4. Then fresh Luna review, exact/small-tolerance fullL3 parity, L128/384/1536 physical probe and serial sameCPU native/profile controls. Do not state speedup until measured.

Optional `underflow-probe/libdspec_profile_ftz-v1.so` setsMXCSR FTZ+DAZ onlywithDSPEC_PROFILE_FTZ=1. L128 samefrozenbinary waveform/spectrum/kernel bit-identical; timing under concurrency is inconclusive. HighL test pending. It is a broad arithmetic perturbation, not accepted production behavior. Luna staticreview noteslegacyDAZsupportcheckcaveat. No Sol finalreviewyet.

## Exact full-grid strategy still needed

Read physical-algorithm-options.md and global-gsh-feasibility.md. Promising exact options: baseline forward+unitreceiveradjoint solve once, spool compactpolecoefficients and apply objectiveq later; bounded paired-theta GSH matrix synthesis using Eigeninternal products (stock GSHMatrix uses external FortranBLAS, not authorized preference); exploit vertical receiver m=0 spheroidal symmetry; batchfrequencies; accumulateeachoutputtile overfreqthenwriteonce. Frequency-first cross-degree coefficient products are generally more expensive; static-tail subtraction/radialcutoffs need separate proof/error budget. Fullglobalmaterialization is still outstanding and must not be replaced by forecast or selectedblocks.

## Build / review locations

- dspec/source-v{3,4,5}-frozen and matching physical_driver binaries/manifests.
- global-gsh/source-v2-frozen and native/profiledv2 binaries. v1 highLsolidtiming stillrunning or near timeout.
- yspec/v5: bounded positionalI/O full kernel API, reviewedLuna+Sol.
- yspec/v6: kernel-probe count64, allradialrows, reviewedLuna conditionalPASS. CLI stoull acceptsnegative counts: use explicitpositive64 only; do not advertise generalboundedCLI untilfixed. Fullmode unchanged.
- review/*.md records independent verdicts and limits.
- results.json from collect_results.py holds completed cases only and excludes partial outputs as results.
- run_case.py now streams hashes and records child CPU/RSS and preload environment; old active runners retain older code. Run unique names, oneCPU/thread,8GiB AS limits, timeouts. Compiler gate algorithm-audit/run_compile.py at most2jobs.

## Active numerical runs at checkpoint

- `dspec-v4-probe-l1536-nq6-a64`: elapsed 1790 s, wall guard 1830 s. Outputs are partial until run.json records success.
- `dspec-v4-probe-l1536-nq6-a64-ftz`: elapsed 643 s, wall guard 1830 s. Outputs are partial until run.json records success.
- `yspec-l1792-t9-step0125-q1-v5`: elapsed 694 s, wall guard 1830 s. Outputs are partial until run.json records success.
- `yspec-v6-probe-l1536-nq6-a64`: elapsed 308 s, wall guard 1830 s. Outputs are partial until run.json records success.

Known unified-exec sessions (poll, do not duplicate):
-93937 DSpecv4L1536normal64angleprobe (1800s).
-6243 sameFTZprobe (1800s).
-2619 GlobalGSHv1L1536onefreq/oneRfullangles (1200s).
-98364 YSpecv5L1792tol1e-9step.125wave (1800s).
-77059 YSpecv6L1536first64angles/allradii probe (1800s; preceded two completedL3gates).
Other numerical sweeps completed. DSpecv5 compiler sessions owned by physical_dspec_luna_impl; ask status before recompile. Full PDF report not generated yet; PDF artifact-start marker ran once, so do not repeat it for this same planned PDF. Standard LaTex report skill and PDF skill read.

## Resume order

1. Inspect completed run.json/logs and v5 compiler manifests; preserve any timeout/failure/partialbin. Never treat partial .bin as complete.
2. Complete v5 review and numerical equivalence gates; profile the removed pole setup separately; assess optional FTZ only from controlled paired evidence.
3. Finish YSpec ODE tolerance check, all-component matched probes and serial sameCPU native/profile overhead checks, with exact exclusive-timer closure. Review v6 after Luna (shortSol pending).
4. Decide broad angular tail resolution honestly under1827ceiling and retain known near-receiver/source localstrain limits. Waveformconvergence does not prove fullkernelconvergence.
5. Present concrete full-materialization cost/algorithm result; if new high-degree library repair/major algorithm phase required, obtain approval under userphase rule. Fullglobaloutput remains authorized target but unfinished.
6. Generate scientific report from results.json with standardLaTexstyle, plots and PDFvisualinspection; rehashproduction/evidence preservation; complete fresh independentreview. No commits/push.
