# Phase A fixed-mesh reconstruction report

**Status:** Phase A accepted; stop here pending a separate decision on Phase B. This is an isolated validation candidate, not a production integration or a complete-pipeline speed claim.

## Scope and result

The frozen case is the attenuated PREM, full-gravity transverse-velocity traveltime calculation at 1000–1200 s: L200, NQ10, refinement 0.025, 35 bins 7–41, dt=2 s, 2048 FFT samples over a 4096 s period, damping 5 e-folds over 2400 s, and filter corners 2/4/8/10 mHz. The captured packet uses active bins 9, 25, and 40. This is distinct from the earlier vertical-displacement physical-energy campaign. The source/receiver packet uses the accepted Stage 2 waveform/adjoint weights and the same model, mesh, angular frequencies, frames, and receiver. No NQ5/6 change, truncation adoption, cache-policy change, or Phase B work was performed. A2 consumed unchanged DS modal radial bytes captured in A1 and reconstructed the six-spin degree-inner fields; it did not rerun the YSpec ODE. The existing Hermitian-dual `H^H`, kernel `Vs`, and one-sided lift conventions were retained. No new-backend equivalence claim is made.

All five planned controls passed: (A1) accepted radial packet capture; (A2) same-packet reconstruction and field/kernel comparison; (A3) L3 finite-difference control; (A4) true-forward-only route control; and (A5) standard forward route in its untruncated and existing truncated submodes. A1–A4 ran in the preserved v3 campaign; A5 and a receipt-backed A4 recheck ran in v4. The v3 campaign is marked failed because of a checker defect after A4 had completed successfully; its original state was preserved. V4 requalified that A4 receipt, ran A5, and completed successfully. No new A1–A4 numerical processes were run during recovery.

| Check | Result | Timing boundary |
|---|---:|---|
| A1 packet capture | 13.238 s | Profiled, three active bins |
| A2 post-radial replay | old DS/Eigen 3.871 s; common degree-inner 3.351 s | Native harness replay; includes map lookup/dispatch/collection, excludes solve, capture and full pipeline |
| A3 L3 finite difference | relative difference 6.4403e-9 | Native, existing 2e-4 gate |
| A4 true-forward-only | driver 6.738 s; supervisor 7.032 s | Profiled; 35 bins, 8 frequency workers |
| A5 standard forward | untruncated 5.811 s; truncated 2.534 s | Profiled; 8 degree workers; two distinct submodes |

A2 passed its approved aggregate field gates across solid and fluid rows for each `(bin, section, role)`: maximum relative L2 was 8.5059e-16 for packed-reference/common and 9.1189e-16 for same-packet DS/common, against 1e-11. Kernel-product comparisons were also at roundoff scale: maximum relative L2 2.0545e-15 (packed/common) and 1.7309e-15 (DS/common). The saved per-fluid diagnostics remain visible. In particular, along-section fluid/adjoint rows show relative L2 values 0.0310849, 0.00266990, and 0.00605437 for bins 9, 25, and 40, with maximum absolute differences 0.0588585, 3.76996, and 0.150429 in the stored tagged Hermitian-dual strain convention. These per-fluid subset ratios are diagnostics, not separate thresholds: fluid rows remain included in the approved pooled solid-plus-fluid acceptance groups. They are not physical SI-strain errors; the fluid kernel `Vs` contribution remains zero. No inference about the cause of these subset residuals is made here.

A3 reproduced the frozen L3 result: common reconstruction `-106.02359024964330`, frozen DS `-106.02359024964328`, and finite-difference estimate `-106.02358956682346`. The finite-difference relative difference was 6.4403e-9, below the pre-existing 2e-4 gate.

A4 recorded 35 true-forward route calls, 42,035 forward RHS columns, 14,035 factorizations, zero unit/weighted adjoint RHS work, and zero receiver QR scopes. A5 matched the A4 receiver to 3.8523e-14 relative L2 in its untruncated submode. The existing truncated policy differed by 6.5837e-11; it remains a control-only result and was not adopted. Both A5 submodes used the prepared single-SEM policy and `nskip=10`; their symbolic-analysis counts can differ from the multi-SEM Stage 2 history.

## Timing interpretation

A2's 3.871 s versus 3.351 s is a post-radial harness comparison, not a production tensor-kernel benchmark or end-to-end speedup. It includes symmetric map lookup, dispatch, and output collection.

A4 and A5 timings are instrumented. A4 uses frequency-parallel scheduling; A5 retains the approved standard degree-parallel 8-worker decomposition, so their elapsed times are not a like-for-like speed comparison. A5 submode elapsed includes cut-index capture/allocation and excludes later CSV serialization. Its master region-wall times were 0.05845/0.36680/5.38377 s (radial/toroidal/spheroidal) for untruncated and 0.02738/0.19831/2.30777 s for truncated. The corresponding worker-region elapsed sums across workers were 0.37916/2.91175/43.06968 s and 0.21865/1.58615/18.46170 s. These are summed elapsed durations across workers, not CPU seconds or master wall time. The largest worker-exclusive scopes were numerical factorization (23.464 s untruncated; 9.864 s truncated), solve (6.484; 2.526 s), spheroidal dynamic-operator assembly (6.159; 2.605 s), and equilibration/scaling (5.989; 2.601 s). Inclusive child rows overlap; their inclusive values must not be summed into a total.

Peak recorded numerical-job resource use was 0.5659 GiB RSS and 9 threads for A4; A1 0.4541 GiB, A2 0.3808 GiB, A3 0.0419 GiB, and A5 0.1419 GiB. The three builds completed with return code 0 in 151.1765 s; peak build use was 2.2575 GiB and 6 threads. The candidate implementation-status file received a post-build completion-log append; the archived build receipt is authoritative for the compiled input tree. All Tunya jobs were detached with `nohup`; campaign caps were 20 threads per numerical job, 8 per build job, and at most 50 aggregate threads. Parallel forward controls used 8 workers; A1–A3 used serial inner work, with Eigen/library teams set to 1 throughout. Resource enforcement used `USER_SPACE_SESSION_WATCHDOG`; this is a user-space watchdog, not kernel aggregate accounting.

## Provenance and evidence

- Frozen build-input candidate tree SHA256 (pre-status-log append): `4c94847f0a46ed96fe166165cc993ee349e9ec0ad0531fd564f89bd58bd49b4e`.
- Key runtime source SHA256s: `phase_a_driver.cpp` `301b9bbc9825a9f6b21a30fcfea48b24c7b5c857fc5d98af44d1e9bbb6626ce9`; `phase_a_reconstruction.hpp` `3644be9b75ce439e7f2e94c06850cdf81e96f0d832bfe77aba376c88c7911fbb`; `phase_a_standard_forward.hpp` `2c33797e91eac442feedbc54ae37adb82ebeb5ebb816b262cd5d3094c1b3296c`; contract test `32d02e53dd7d9e034ee774d8fe469f981957b2c9c3cde2faa09b402122405beb`; copied diagnostic hook `Windowed3DFields.h` `5ad3ea5197d0b2880b56c9a4292bcdec67983d3e412638e18004bc2c2c7745e7`.
- Native/profile/tiny binary SHA256: `514bef2e13bff9d84a08e540b821babe261fb3eaa3d24a7e0421d0755f77db85`, `f6d52ca5439a560cedf0083a2241993784b477c37eeb7e191791593a6c3a7020`, `2b7953f08c70ff419622688d29799485b6a6dc0ee146d9972bae2424f58ca0db`.
- V4 terminal state SHA256: `15d723f4e25ff1ede30b32b185abd6ec353b061ce8d361d3040236ee370a9d72`; completed remote task: `/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya/phase-a-v4-controls/continuation/phase-a-v4-controls`.
- V4 full remote verification covered 10 output artifacts and 3 logs, all three jobs return code 0, unchanged inputs, and verified cleanup. The local V4 snapshot intentionally omits `A5-profiled-standard-forward-two-submodes/data/parameters.txt`; its remote hash is present in `receipts/phase-a-v4-controls-completed/snapshot-verification.json`. V3's local snapshot has 13 output and 11 log hashes; 10 larger artifacts remain remote-only and were independently verified.
- Build receipt: `receipts/phase-a-v1-build-completed` (state SHA256 `03f9e98c0977ca38ddf675097f11239d8f4c0a1757252caa50b9dd679ab5806b`). Preserved failed setup/checker attempts and their exact states remain under `receipts/phase-a-v2-controls-failed` and `receipts/phase-a-v3-controls-failed`.
- Detailed JSON, raw A1–A5 component/profile/cut-index artifacts, and their hashes are linked from the companion [validation JSON](phase-a-reconstruction-20261031-report.json). A2's full field-comparison CSV was not included in the local snapshot; its accepted group diagnostics are preserved in the hashed summary, while raw remote artifacts remain under the Tunya task receipt.

Phase A is complete. No Phase B implementation or execution is included.
