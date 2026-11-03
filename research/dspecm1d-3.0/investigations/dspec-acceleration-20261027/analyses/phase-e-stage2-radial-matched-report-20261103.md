# Phase E Stage 2: radial reconstruction result

## Scope

This stage changes only native radial-mode gathering, derivative construction, and strain-table assembly in the isolated DSpec candidate. The generic visitor remains selectable as the same-source baseline and numerical oracle. Angular synthesis formulas, GEMM, tile sizes, and the modal solver are unchanged; the fast route also uses direct-index access to the existing angular coefficient tables.

The matched case uses L=200, NQ=6, h=0.02035382429229197, 35 bins (7–41), 208 native sides, and 135/19 along/across angles. Timed jobs use eight frequency workers, one inner solver thread, and the existing exact-zero skips at bins 7, 8, and 41. The all-bin oracle audit disables those skips.

## Correctness

The 35-bin keyed radial oracle compared 21,715,400 keyed radial records covering node/order/family/field-role entries; this count is records, not independent model solves. The six-component radial packets (three values and three derivatives), per-mode strain contributions, forward/adjoint inventories, and fluid/solid subsets each had maximum relative L2 error 0. This is packet and contribution parity, not a separate comparison of accumulated strain tables. The largest conservative diagnostic capture estimate was 210,949,600 bytes per worker, below its 256 MiB audit bound.

The same-workspace end-to-end field and complex-product proof covered all 35 bins and both sections. Maximum relative L2 errors were 8.631317934446982e-16 for forward fields, 7.389022962721731e-16 for adjoint fields, and 1.0882120385497865e-15 for complex products. All selected-grid output channels and grids matched the frozen Phase D result exactly.

The prior hash-bound Phase D final report links the NQ10 accuracy evidence used for the 0.1% support bound: 2.1952474713321954e-4 for solid complex fields and 1.450918716573929e-4 for the solid full-35-bin kernel. That report SHA-256 is `013273941560685863015e48aeef42daeada129b3174ff47217637f21e605a95`; it is the report hash, not a raw NQ10 evidence-file hash.

The actual L3/NQ8 batched path also ran the frozen finite-difference control. Its relative error was 6.440261900046464e-9 against a 2e-4 limit; all 35 radial proof bins passed, and the L3 outputs matched the frozen Phase D outputs exactly.

## Matched timings

Three sequential same-case repeats gave these driver elapsed times:

| Route | Repeats (s) | Median (s) | Supervisor wall median (s) |
|---|---:|---:|---:|
| Generic baseline | 7.800024017, 7.781501814, 7.788696483 | 7.788696483 | 8.057737327 |
| Batched radial | 7.576651356, 7.589408306, 7.585382455 | 7.585382455 | 7.878134965 |

The batched route reduced median driver elapsed by 0.203314028 s, or 2.61037297%, in this run. This is a modest full-pipeline reduction; the batched DSpec driver median remained 1.7615 times the YSpec context median. One profile run per route measured 4.965903364 versus 4.017802243 worker-seconds in the radial reconstruction scope (19.09% lower); those worker sums are not wall time and should not be subtracted from the pipeline wall. The route also changes angular coefficient table access: coefficient row lookup/packing was 1.972437253 versus 1.728386235 worker-seconds, and generic lookup materialization was 0.004497936 worker-seconds while the batched route removes that scope. Angular GEMM was essentially unchanged at 19.213565570 versus 19.224154517 worker-seconds. Thus the full-pipeline difference cannot be attributed to the radial scope alone. The radial scope still includes per-order callbacks and the existing strain formula; these profiles do not split out finer radial substeps. Profile-only driver elapsed was 7.814170005 s generic and 7.597191790 s batched.

The three YSpec context-refresh driver times were 4.298078351, 4.306145970, and 4.306802959 s (median 4.306145970 s). This is context only; YSpec and DSpec use different arithmetic paths.

## Build, resources, and receipt history

Both timer-disabled and timer-enabled builds succeeded. Each took about 75 s; the sequential build plan took 151.07 s overall, peaking at 6 observed threads and 2.2914 GiB RSS. Tracked source and dependency hashes were unchanged before and after each variant. The independent comparison found all 12 shared external dependency paths identical to the reviewed Stage 1 v2 build; only the candidate source tree differed. The pinned Stage 1 YSpec DSO hash remained `9224a655ba4fcbe9dc153a5f0c233f46c741741f0fa1d08629117e477df5e738`.

All 14 scientific and finite-difference jobs completed with return code 0 and verified cleanup. The first terminal-check attempt failed after the scientific checks because the supervisor updated its own state file during the check and the checker’s output directory did not exist. The preserved record shows exactly one changed input (`state.json`); the other 275 input hashes were unchanged. A reviewed checker-only retry read the finalized immutable state, validated all 14 prior jobs, and completed with return code 0. No scientific jobs were rerun.

Numerical jobs stayed within the 20-thread and 16 GiB per-job caps; the highest reported job peak was 9 threads and 2.1153 GiB RSS. The supervisor uses a user-space watchdog; this is not kernel aggregate accounting.

## Evidence

- Machine-readable terminal report: `analyses/phase-e-stage2-radial-report-20261103.json` (SHA-256 `6cadfe132104162ad831d11bffd80caf1ba172f9ed2bba6732ed035452b5bcc8`).
- Compact raw evidence archive and manifest: `receipts/phase-e-stage2-radial-evidence-canonical-20261103.tar.gz` and `.json` (archive SHA-256 `c8ec48d8a6a64ce1fb4ae526a47f304957c46a34248eed871365472556cac9d9`). They include both build command records, all profile JSON, all 35 L200 and 35 L3 proof rows, the 70-row scalar parity CSV, all scientific result JSON, the finite-difference result, supervisor states, and failed/retry checker logs.
- Build dependency detail: `receipts/phase-e-stage2-radial-v1-build-evidence-20261103.json`.
- Independent source/build review: `reviews/phase-e-stage2-radial-source-build-independent-review-20261103.txt`.
- Independent run/checker review: `reviews/phase-e-stage2-radial-run-checker-independent-review-20261103.txt`.

This completes only the radial reconstruction substage. No subsequent angular or modal-solver optimization substage, production integration, or discrepancy investigation was performed.
