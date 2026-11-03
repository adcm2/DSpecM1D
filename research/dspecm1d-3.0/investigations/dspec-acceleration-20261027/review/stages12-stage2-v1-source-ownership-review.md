# Stage 2 frequency scheduler: implementation and ownership review

Date: 2026-10-31. Candidate: `candidates/stages12-stage2-frequency-v1`.

## Frozen provenance and allowed delta

The candidate was copied byte-for-byte from the accepted archive extraction at `receipts/stages12-stage1-v1-frozen-inputs/candidate`, archive SHA256 `52fccc457361e1fc8b1e865f05738f93a04e933f2725ac24169f59770de3ad63`. Its 20-file seed manifest is `analyses/stages12-stage2-v1-baseline-candidate-sha256.json`, SHA256 `f9e29601f520c253f618f2d7ca386adcf3263e785401ccc4e1ec2bc7267fe7c4`.

Candidate source changes are confined to `doc_driver.cpp` and `profile.hpp`; `frequency_scheduler.hpp` and its small test are additive. `stage-to-math-map.json` and `implementation_status.md` are reporting/provenance changes. The Stage 2 checker requires that exact file inventory and rejects any other baseline-file changes. It pins the frozen `measure()` implementation SHA256 `7b5aa2a57fb90d4ea858159bdcda5312524bb9d687127683e673a185183cce34` and unchanged `kernelVs()` formula SHA256 `99687637e94756b9342cb758570a5e165d3fea087bd674dd53532e753412e295`.

Fresh Stage 2 preflight is recorded in `receipts/stages12-stage2-v1-preflight.json` (SHA256 `f12feb49fff452f1f560524f0f0a9b2dfd6067dd74630e2584c98cefc136951a`). It matched all 12 pinned source/dependency tree hashes from the successful historical build. DSpec HEAD was `a778b13b22dca9709e1f96604b62925707654e5d`, preserving five pre-existing tracked modifications and untracked `validation/phase0/`; YSpec HEAD was `280eab0200495391539c0f97e62f978083c307bc` and clean. The linked YSpec library hash is `f891e44ed9d9ca8aaa871fcbb3d9287526aee90dc33bf4e58867abc14e9a6d8f`. No production checkout or dependencies are edited.

## Scheduling and ownership

- `DSPEC_FREQUENCY_WORKERS` opts into `dspec` or `dspec-forward`, accepts only a fully parsed integer 1–8, and requires CLI solve workers 1. With the variable unset, the legacy degree route is unchanged. Eigen/BLAS/OpenMP inner pools remain one thread.
- Both passes construct a fresh `WorkerCapture` and fresh `WorkerData` slots for each outer OpenMP worker. The context factory is `noexcept` and compile-time checked; context teardown is `noexcept`, restores the thread-local diagnostic backend, and the capture destructor checks local scope depth.
- Each frequency task writes only its indexed raw/component or contribution/error/validation slots. `SEM`, source, parameters, model, grid and the single response cache are shared at stable addresses; the same response cache is used in both passes. Both workspace calls use `solveWorkers==1`.
- Kernel tasks use independent, empty `SamplingPlanCache` values for every `(worker, section)`, constructed once per pass and kept across that worker's statically assigned frequencies. No cache copy or shared plan pointer is used.
- The contribution allocation uses checked size arithmetic and is capped at 12,634,160 bytes. Workers assign the unchanged `kernelVs()` value into one frequency slot. After every task joins successfully, the parent applies the existing `+=` in increasing `k`, then section, angle, and native-node order. No OpenMP reduction or shared `+=` exists on the frequency route.
- Task failures are captured by bin, the remaining scheduled work joins, and the parent rethrows the lowest failed bin with pass context. Result serialization occurs only after successful completion.
- The scheduler test exercises indexed out-of-order completion, deterministic lowest-bin selection across multiple errors after all work finishes, worker-context lifetime, and an injected callback capacity exception. The capacity test is an injected callback failure; it does not replace or modify the actual `BatchResponseCache`. The actual shared cache implementation remains unchanged and is covered by scientific cache/factor/RHS invariants in controls and full measurements.

## Profiling and output semantics

Worker TLS rows are merged only after a pass with `frequency_worker.`. The emitter strips exactly one worker prefix for exact scope-name lookup; unknown names remain unclassified and raw names remain intact. The map is seeded from the 101-scope accepted Stage 1 audited map plus six explicit Stage 2 scope names.

Per-bin diagnostic snapshots, durations, and overhead stay private until the join. The writer removes diagnostic deltas only from frequency-worker production rows and counters; master-wall rows are never adjusted by summed worker elapsed. The raw inclusive parallel-pass umbrella retains diagnostics and is labeled as such. Diagnostic rows/counters appear once in a distinct overlay. The profile checker enforces `diagnostic=true` only in that overlay and rejects worker rows mislabeled as master wall.

The frequency route reports `kernel_s` as kernel setup + parallel pass including diagnostics + parent ordered reduction. It also writes `kernel_parallel_pass_wall_s`, `reduction_wall_s`, `diagnostic_worker_elapsed_s`, and timing semantics. No sum of worker diagnostic time is subtracted from master wall. The legacy degree route's prior timing semantics remain labelled separately.

## Tests and limits

Locally, all 13 Stage 2 checker synthetic tests pass, covering byte parity, receipt integrity, conditional across outputs, worker metadata, wall closure, diagnostic overlay duplication, frozen candidate scope, and build-receipt binding. The scheduler helper executable compiles and passes. Candidate source-contract/hash checks pass. These are implementation gates only; no C++ driver build, race detector, numerical control, or L200 run has yet been made from this Stage 2 candidate.

The pinned read-only call-tree audit in `plans/stages12-stage2-thread-safety-audit-draft.txt` found no lazy shared mutation on the selected SEM, YSpec model, interpolation, or double-precision Wigner paths. It explicitly does not establish dynamic race freedom or cover alternate routes. Its source hashes and limits remain part of the release record. Runtime build and numerical checks remain mandatory before acceptance.
