# Phase E Stage 2: angular reconstruction result

## Scope

This stage adds a real-Legendre angular reconstruction route to the isolated DSpec candidate. The same-source baseline retains the legacy complex angular route, and **both timing routes use the already validated batched radial route**. The new route sums ordinary real and imaginary degree products against real Legendre weights, then applies the degree-independent complex phase/rotation. It batches coefficient packing and angular tiles, shares spin-0 weight work, and specializes the proven degree-zero weighted-adjoint radial frame with its existing Legendre normalization. The modal solver and production sources were not changed.

The matched L=200, NQ=6 workload used mesh scale 0.9, h=0.02035382429229197, 35 frequency bins (7–41), 208 native sides, and 135/19 along/across angles, with eight frequency workers, one inner solver thread, and exact-zero skips at bins 7, 8, and 41 for timed jobs. The 35-bin numerical audit turned zero-bin skips off. A conservative derived real-angular scratch bound was 23,663,080 bytes per worker (<32 MiB); the legacy route bound was 17,940,200 bytes and the L3 real-route bound was 5,177,576 bytes. These are derived working-set bounds, not measured allocations or RSS. The actual L3/NQ8/h=0.05 finite-difference control ran the new real-Legendre route.

## Correctness

The same-workspace angular oracle compared 70 section/bin rows across all 35 bins. Maximum relative L2 differences were 8.378e-16 for forward fields, 7.168e-16 for adjoint fields, and 1.086e-15 for complex products. The independent scalar pipeline proof also covered 70 rows; maximum relative L2 error was below 4.15e-16. All 35 radial-oracle bins passed.

Compared with the final validated radial-stage result on the same physical mesh and grids, the all-bin output differences were 2.232e-15 for `along.bin`, 1.468e-15 for `across.bin`, and zero for components, spectrum, and trace channels. The unchanged mesh and full-pipeline parity retain the existing 0.1% NQ10 accuracy support: maximum solid complex-field error 2.1952474713321954e-4 and solid full-35-bin kernel error 1.450918716573929e-4, linked by the hash-bound Phase D final report (`013273941560685863015e48aeef42daeada129b3174ff47217637f21e605a95`).

The L3 real-Legendre finite-difference audit covered 35 one-section angular-oracle rows. Its relative FD error was 6.440261900046e-09; the greatest difference from the frozen L3 result was 3.186e-16 in `along.bin`, with zero component, spectrum, and trace differences. This is below the frozen 2e-4 FD limit.

## Matched timings

Driver elapsed times come from completed `data/result.json` files; supervisor task walls are reported separately. Each timing pair uses the same workload and environment, with legacy complex followed by real-Legendre in the reviewed interleaved sequence.

| Route | Driver elapsed repeats (s) | Median (s) | Supervisor wall repeats (s) | Median (s) |
|---|---|---:|---|---:|
| Legacy complex angular, batched radial | 7.608922249, 7.601138393, 7.581567730 | 7.601138393 | 7.892535, 7.992918, 7.894200 | 7.894200 |
| Real-Legendre angular, batched radial | 4.899015515, 4.927772044, 4.967048291 | 4.927772044 | 5.304443, 5.222636, 5.363008 | 5.304443 |

The real-Legendre route reduced median DSpec driver elapsed by 35.17% (2.673366 s). The three fresh YSpec driver runs were 4.358693065, 4.312707284, and 4.308166653 s (median 4.312707284 s). DSpec driver elapsed remains 1.1426 times the YSpec context median; these use different arithmetic paths. For a separate supervisor-wall view, the legacy DSpec, real-Legendre DSpec, and YSpec medians were 7.894200, 5.304443, and 4.555655 s, respectively; that comparison includes task orchestration overhead and is distinct from driver elapsed.

One profile per route observed 21.746719850 to 4.662851660 worker-seconds in the tiled angular reconstruction scope. Its worker sum overlaps the frequency-pass wall and must not be subtracted from driver elapsed. The GEMM profile counted 147,840 legacy complex GEMM invocations versus 55,104 real GEMM invocations in 27,552 real/imaginary pair scopes on the new route. Coefficient row lookup/packing was 1.699598032 to 0.219113381 worker-seconds. Tile scatter rose from 0.094675903 to 0.357326292 worker-seconds with the post-sum phase work. The same profile showed 24.3717 of 33.5702 total worker-exclusive seconds in modal scopes (72.6%). These are single-profile worker-time observations, not repeated wall-time measurements or values to subtract from driver elapsed.

## Build, resources, and evidence

Both timer-disabled and timer-enabled isolated builds passed. Their supervisor build task took 156.71 s, peaking at six observed threads and 2.369 GiB RSS; shared external dependency hashes matched the reviewed finalized radial build. The pinned YSpec context used the unchanged Stage 1 native DSO.

All 14 numerical and finite-difference producer jobs completed with return code 0, unchanged inputs, and verified cleanup. A separate terminal checker completed with return code 0 against an immutable producer-state snapshot. Maximum numerical-task peak was nine threads and 2.1825 GiB RSS under the 20-thread and 16-GiB per-job caps. The supervisor uses a user-space watchdog; it is not kernel aggregate accounting.

- Machine-readable terminal result: [`phase-e-stage2-angular-report-20261103.json`](phase-e-stage2-angular-report-20261103.json), SHA-256 `ab867a3caf8da0a4b63eb65c9d000f3981d3818a7ac2513881e608d71c938468`.
- Compact raw result and proof bundle: [`phase-e-stage2-angular-compact-evidence-20261103.tar.gz`](../receipts/phase-e-stage2-angular-compact-evidence-20261103.tar.gz), SHA-256 `7b4e6faeb4e7cf3d8eed40366c994b3ea304607b653c84cecc1132ae8e93f674`. It contains the 70-row angular and scalar proofs, 35+35 radial proof rows, profiles, result JSONs, FD result, immutable producer state, and the manifest binds all 242 declared outputs.
- Compact-manifest format-correction sidecar: [`phase-e-stage2-angular-compact-evidence-20261103-manifest-canonical-bound.json`](../receipts/phase-e-stage2-angular-compact-evidence-20261103-manifest-canonical-bound.json); it preserves the original malformed bytes and binds their hash to the parseable equivalent content.
- Build evidence: [`phase-e-stage2-angular-v1-build-evidence-20261103.json`](../receipts/phase-e-stage2-angular-v1-build-evidence-20261103.json).
- Independent source/build review: [`phase-e-stage2-angular-source-build-independent-review-20261103.txt`](../reviews/phase-e-stage2-angular-source-build-independent-review-20261103.txt).
- Independent run/checker review: [`phase-e-stage2-angular-run-checker-independent-review-20261103.txt`](../reviews/phase-e-stage2-angular-run-checker-independent-review-20261103.txt).

This completes the angular reconstruction substage. No subsequent modal-solver or other optimization substage, production integration, or discrepancy investigation was performed.
