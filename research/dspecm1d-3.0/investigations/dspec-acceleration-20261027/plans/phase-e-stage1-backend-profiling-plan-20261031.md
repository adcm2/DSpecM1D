# Phase E Stage 1: matched backend profiling plan

Status: source and configuration are implemented in the isolated Phase E candidate and are awaiting bounded independent review. No build or numerical run has started.

## Scope and frozen workload

This stage profiles the validated Phase D selected case only, then stops. Preserve the selected NQ6 scale-0.9 mesh (`h=0.02035382429229197`), degree 200, 35 frequency bins, 208 physical one-sided sides, 135 along-track and 19 across-track angles, and the existing full-gravity/measurement configuration. Use the exact Phase D case and sampling inputs recorded in `configs/phase-d-selected-grid-timing-plan.json` and `receipts/phase-d-selected-grid-timing-run-v1/state.json`.

The requested common observable is relative `lnVs`: DSpec contracts this parameter directly, whereas the pinned YSpec case requests `ParameterSet::Velocities` and returns a full velocity-parameter vector before the driver selects `vs`. Record the actual returned parameter names/channels and counts. This is the same requested case/output, not identical backend arithmetic or computational work. Also report the unchanged YSpec radial-solver defaults and available source/receiver-depth and integration/response dimensions; the Phase D SEM mesh is not assumed to describe YSpec's internal radial discretization.

The DSpec route uses eight frequency workers and one inner numerical-library thread. The YSpec route uses CLI workers 8, explicitly unsets `DSPEC_FREQUENCY_WORKERS`, and pins numerical-library inner threads to 1. Run all timing jobs sequentially on the reserved CPU set, with no competing numerical or build job.

## Source and binary provenance

- Preserve the successful D v2 source/build and all Phase D evidence. Seed a new isolated candidate from the exact D source files and manifest in `phase-d-nq6-batched-v2-inputs.tar.gz` / `phase-d-nq6-batched-v2-inputs.json`, verified against `receipts/phase-d-nq6-batched-build-v2/state.json`. Do not edit the Phase D candidate or binaries.
- The native control executable is the frozen Phase D v2 native driver at `/space/adcm2/research/builds/phase-d-nq6-batched-v2/phase-d-native/phase-d-native-driver` (SHA-256 `c780867c4fdd84639e75cb29516b70868d02c7ad54cb4437ab87f5c95566b223`). Its YSpec dependency is `/space/adcm2/research/builds/documentation-wedges-20261027-01/yspec-library/libyspec-doc.so` (SHA-256 `f891e44ed9d9ca8aaa871fcbb3d9287526aee90dc33bf4e58867abc14e9a6d8f`).
- Pin YSpec to commit `280eab0200495391539c0f97e62f978083c307bc`. The recorded Stage 2 preflight has a clean checkout, `src` tree SHA-256 `26f29edb53e94fda4859f3bb962b49b33cd3c8da8e1d3433a3037a6c7608e692`, `YSpec` tree SHA-256 `eb023ebfd0953a9b60d5be3b4d7c652ff320f573755b4cbd4913ca3118615853`, and the library SHA above. Verify current tree and library hashes again before staging. Build the instrumented YSpec library from a fresh source archive at that exact commit into a new build root. Do not edit the live checkout or describe an earlier unrelated profiler as this profile.
- Build fresh variants from the same isolated candidate/source and compiler options: a D driver with timers compiled out, a D driver with D scopes enabled, a YSpec DSO built from the exact pinned source with its recorder compiled out, and the same YSpec source with scopes enabled. Native controls use the newly built timer-disabled artifacts so instrumentation overhead is not confounded by comparing old and new builds. Use the frozen D v2 executable and original YSpec library only as immutable output/provenance anchors; require parity before interpreting timings. DSpec profile runs pair the D instrumented driver with timer-disabled YSpec; YSpec profile runs pair the timer-disabled D driver with instrumented YSpec, isolating each backend's instrumentation. Bind every executable, DSO, case, manifest, plan, and checker by SHA-256 in the build/run receipts.

## Instrumentation boundaries

Keep instrumentation observational: no numerical expression, traversal order, reduction order, scheduling policy, or worker count changes. Use steady-clock scopes around a frequency, degree, radial solve, sampling batch, or contraction phase; do not time scalar inner-loop iterations.

For DSpec, retain existing modal and pipeline scopes and add narrowly split scopes in `phase_d_batched_reconstruction.hpp` for angular-plan/map lookup, phase/rotation weight-tile packing, Eigen GEMM, tile scatter/accumulation, and projected radial interpolation/strain table preparation. Keep material evaluation/weak-form contraction and ordered reduction as separate scopes. Add call/sample counters at batch or degree boundaries. The existing D profile uses thread-local maps and merges worker data after joins; extend that mechanism without adding a global lock per event.

For YSpec, instrument the pinned source at these boundaries:

1. `DepthSolves` setup and `Record` degree work, split into toroidal/radial/spheroidal radial solve/integration scopes.
2. `StrainOf`/`StrainTable` coefficient-to-strain interpolation and construction, counted per degree/node batch.
3. `AngularFunctions` preparation/evaluation and the sampled-kernel angular synthesis blocks.
4. Material-parameter contraction/`Transform` and accumulation into sampled outputs.
5. Per-frequency ordered/local reductions, plus top-level setup, output I/O, and measurement/transform work.

Use a YSpec-owned recorder in the instrumented shared library. Reset it immediately before the outer YSpec API call. Each thread accumulates into its own fixed set of per-scope inclusive/exclusive nanoseconds and call counts. At the end of each existing frequency work item, after the scopes have closed, merge that thread's record once into the run aggregate; merge remaining main-thread scopes once after the worker join. The OpenMP master thread participates in frequency work, so each WorkItem preserves its preexisting caller entries and reports a separate master-thread work-item time/count bridge; caller API wall, worker-summed time, and that bridge are overlapping domains, not additive wall partitions. Finalize/write only after the outer API returns. The D driver calls explicit reset/finalize entry points exported by the DSO. Timer-disabled `Scope`/`WorkItem` compile away inline in both YSpec and D drivers; profile-enabled YSpec has a separate D driver linked to the enabled DSO. This avoids assuming inline TLS symbols are shared across the executable/DSO boundary and avoids a mutex on each scope event. The profile output reports bin coverage, calls, inclusive/exclusive time, caller totals, worker totals, and the master-thread bridge; supervisor/driver wall remains separate.

Before any profile run is trusted, require finite inclusive/exclusive totals, positive expected-scope calls, and explicit work coverage: 35 forward bins, plus kernel coverage of 35 work bins or 32 work bins and the three exact-zero skips, with the resulting count exactly 35. Keep D unit-adjoint/skipped-bin counts separate. Require positive D angular packing/GEMM/scatter counts. Record actual observed threads from the supervisor alongside configured DSpec 8/1 or YSpec CLI-8/inner-1 settings; metadata alone is not thread-scheduling proof.

## Run matrix and order

Use fresh output directories for every invocation and execute sequentially:

| Order | Backend | Build | Repeats | Purpose |
|---|---|---|---:|---|
| 1–3 | DSpec | fresh timer-disabled D driver + timer-disabled YSpec DSO | 3 | Same-workload native timing distribution and parity controls |
| 4 | DSpec | fresh instrumented D driver + timer-disabled YSpec DSO | 1 | DSpec phase profile and instrumentation-overhead context |
| 5–7 | YSpec | same fresh timer-disabled D driver + timer-disabled YSpec DSO | 3 | Same-workload native timing distribution and parity controls |
| 8 | YSpec | same timer-disabled D driver + instrumented YSpec DSO | 1 | YSpec-only phase profile and instrumentation-overhead context |
| 9 | checker | immutable staged checker | 1 | Read-only receipt, output, timing, profile, and parity audit |

Keep the CLI case, mesh selector inputs, sample grids, model, measurement, bin set, and output observables identical to the Phase D selected-grid run. Each timing invocation writes a new result tree; no output directory is reused. Timed runs keep validation/oracle CSVs disabled, as in the frozen Phase D timing run.

## Acceptance and interpretation

- Require each backend's three native outputs to be finite and grid-identical, and compare their physical outputs to the corresponding existing Phase D selected-grid native output with per-channel relative L2 at `1e-11` (exact grids/window/bin ordering). Compare the instrumented run to the same backend's native reference at the same gate. Do not gate YSpec against DSpec: their raw sampled arrays already differ measurably; report that comparison only as context.
- Report each invocation's supervisor wall time, driver wall time, peak RSS, and peak observed threads. For profiles, report master wall and summed worker durations separately. For each nested scope give calls, inclusive and exclusive worker seconds; inclusive rows overlap children and worker sums overlap across concurrently scheduled threads. Do not subtract summed worker time from wall or treat it as wall-time share.
- Summarize native repeat median/min/max/range (and relative range), then describe the instrumented run separately as one observation versus the native median; do not claim a speed change from a single profile run.
- Report configured caps (build: 8 threads/32 GiB; numerical: 20 threads/16 GiB; aggregate ceiling: 50) and the user-space watchdog limitation. Run no more than 8 solve workers with one inner numerical-library thread and recheck resources before each detached launch.

Stage 1 ends with the reviewed side-by-side profile, provenance, repeat spread, instrumentation overhead context, and any parity limitation. Do not begin optimization, production integration, discrepancy work, or Stage 2 in this task.

## Review gate

Sol must review the actual candidate scopes, YSpec recorder/export boundary, build configuration, run plan, and checker before any build is staged. If the YSpec profile export cannot prove 35-bin coverage without changing scheduling or adding per-event synchronization, stop for a focused design correction before building.
