# Physical PREM profiling phase

User authorized DC resolution, waveform/energy convergence and full matched component profiling. Initial protocol saved before implementation. Production and earlier evidence remain read-only.

Added parent-owned run_case.py: one CPU/thread, per-process CPU/address-space limits, unique output directories, command/binary/output hashes and raw logs preserved on failure. This is per-process accounting, not an aggregate kernel memory limit.

Added analyze_waveforms.py to independently recompute .5*dt*sum(Hann*u_z^2) and true-relative weighted/unweighted window errors from SI trace CSV. Fixed1e-3 refinement gates are recorded, with unweighted/peak diagnostics.

Mathematical policy clarified after read-only Sol analysis: explicit DC at finite damping; firstoutofband guard haszero taper; fixed-epsilon padding separate from contour-filter damping sensitivity; continuous-timeYSpec energyadjoint samples excludeanextra dt. Saved in protocol.json and dc-math-review.md.

Runner creates the empty per-run data directory explicitly so both adapters share the same output guard.

Initial YSpec run stopped at its output guard because the runner pre-created data/. Preserved failed run yspec-initial-l32-t7-q1. Removed data/ precreation: adapters now own creating the absent data directory; parent still guards the unique run directory.

Parent adapter audit found YSpec v1 trace reconstruction used an SI taper on a nondimensional grid. Preserved its successful-but-rejected waveform run and requested an immutable v2 repair. Added independent NumPy SI-spectrum reconstruction check to detect filter/normalization mismatches; no tolerance relaxation.

Parent v3 full-native radial/complete L3 grid controls executed: pole/reference/profile/contact-enabled volume outputs compared in kernel-gates-v3.json; finite/informative kernels and profile closure checked. These are controlled low-degree algorithm gates, not a resolved physical Earth kernel. User explicitly chose full global materialization; large-degree disk and runtime assessment continues after waveform convergence.

Added explicit --cpu-index to run_case.py for independent untimed accuracy/probe scheduling; all final paired headline timings will use one identical CPU and run serially, with no compiler overlap. Each numerical process still has one thread.

2026-10-26: Derived validation_inputs/dspec-l3-baseline-721.csv from the frozen v3 L3 complete-native-grid trace, retaining header and first 721 rows through 3600 s for matched YSpec streaming-vs-monolithic checks. Original trace preserved. Exact global GSH v1 low-L rotation/field gates executed on complete grids (L3/6/12 at native side118; L12 side175); numerical results retained under runs/gsh-v1-gate-*.

2026-10-26: Added analyze_time_quadrature.py to reconstruct the identical saved contour spectrum at 5/2.5/1.25/0.625 s spacing, with fixed epsilon, FFT period and frequency bins. This isolates Hann-window time-quadrature convergence from padding/frequency convergence and requires no duplicated Earth solves. YSpec v5 streamed-vs-v3 monolithic low-L validation passes exact values after explicit layout transpose; native/profile output values also exact (yspec-stream-gates-v5.json). Full high-degree materialization remains pending degree convergence and measured runtime assessment.

2026-10-26: Changed runner hashes to bounded 4 MiB streaming reads so validating a full global output cannot allocate the entire file. Runner now records optional DSPEC_PROFILE_FTZ/LD_PRELOAD settings and preload hashes for the isolated gradual-underflow diagnostic; numerical executable arguments and one-thread/per-process guards are unchanged. Older running processes retain their original runner code and hashes.

2026-10-26: Saved runner v2, then added its own source hash plus child user/system CPU time and exact Linux child peak-RSS accounting via getrusage after completion. Sampled VmHWM alone underestimated sub-0.1-s diagnostic peaks. CPU/wall evidence helps distinguish concurrent accuracy-job contention from the final serial timing controls. Existing numerical runs and their metadata are preserved.

2026-10-26: Added collect_results.py to assemble completed saved cases and independent trace-based energy/refinement/cross-code checks into a single results JSON. Pending/failed runs are explicit, concurrent accuracy timings are not labelled controlled comparisons, and high-degree full-kernel materialization remains explicit as pending. Added full-grid-sizes.json using the cached GSHTrans FastFFTSize factor rule, with output/scratch size estimates distinguished from existing files.

2026-10-26: User requested ETA before shutdown. Saved CHECKPOINT.md with current verified evidence, outstanding full-global materialization, exact hotspot/v5-build state, active bounded run sessions, reviews and resume order. Rehashed438 protected files: zero changes. No production code changed.
