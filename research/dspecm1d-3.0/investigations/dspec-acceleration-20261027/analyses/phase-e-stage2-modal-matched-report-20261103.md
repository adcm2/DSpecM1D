# Phase E Stage 2: modal operator cache result

## Scope

This substage caches exact SEM stiffness, mass, and attenuation components and their sparse assembly structure by solve family and degree. The original modal assembly remains selectable as the same-source baseline. Both routes in both timing blocks use the already validated batched radial and real-Legendre angular paths. The accepted COLAMD solver is unchanged, and each family/degree/frequency still receives a fresh numerical factorization. The cache cold build is inside the measured full-driver path.

The L=200, NQ=6 case used mesh scale 0.9, h=0.02035382429229197, 35 frequency bins (7–41), 208 native sides, and 135/19 along/across angles. Timed jobs used eight frequency workers, one inner solver thread, and the existing exact-zero skips at bins 7, 8, and 41; the modal audit disabled those skips. The L3 finite-difference control used the actual cached route at L=3, NQ=8, h=0.05, one frequency worker, and the frozen 157-side grid.

## Numerical result

The all-bin modal oracle covered 14,035 unique family/degree/frequency keys; L3 covered 245 keys. Cached and legacy operator, full physical solution, and independent probe solution equivalence errors were exactly zero across all keys. The maximum backward error was 8.931e-17 in the complete L200 proof and 6.982e-17 for L3. Full-pipeline output differences against the final angular result were exactly zero in trace, spectrum, components, and both section arrays.

**The dual-identity audit definition changed after a preserved failure.** At bin 33, R0, degree 0, the original scalar relative quotient was 1.014258591984e-11 for both cached and legacy routes, with identical raw error 7.640312286e-06 and scalar denominator 753290.3686. The full final CSV shows a maximum scalar diagnostic quotient of 1.111113330628e-11 at bin 40, R0, degree 0 for both routes; the L3 maximum was 3.085887646994e-11 at bin 24, R0, degree 0. These quotients are cancellation-sensitive relative to the matrix and solution norm scale, not absolutely near-zero scalar products. The original scalar criterion therefore did not pass for either route, and this report does not claim that it did.

The accepted same-workspace dual check is separately recomputable for each route: `|yᴴHx − (Hᴴy)ᴴx| / (||H||F ||x||₂ ||y||₂) ≤ 1e-11`, without a denominator floor; a zero denominator requires an exactly zero numerator. The maximum cached/reference normwise ratio was 1.717e-17 (bin 8, family 2, degree 168). The separate operator, solution, and backward-residual limits remain 1e-11.

The actual L3 finite-difference relative error was 6.440262e-09, below its 2e-4 limit. The unchanged mesh and full-pipeline parity retain the existing hash-bound NQ10 solid-accuracy support: maximum solid field error 2.195247e-04 and full-35-bin kernel error 1.450919e-04. The linked Phase D report hash is `013273941560685863015e48aeef42daeada129b3174ff47217637f21e605a95`.

## Matched timings

The two three-pair blocks are kept separate. Driver elapsed values come from each completed `data/result.json`; supervisor task walls are shown separately. Pairing used the same DSpec candidate workload and binaries. The second block reversed route order to check timing drift. The ranges overlap substantially; no pooled median or selected best run is used.

| Block | Route | Driver repeats (s) | Driver median (s) | Driver range (s) | Supervisor-wall median (s) |
|---|---|---|---:|---:|---:|
| Original | Legacy modal assembly | 6.000284, 5.554909, 5.961739 | 5.961739 | 5.554909–6.000284 | 6.337082 |
| Original | Cached components | 6.140654, 5.880600, 5.034904 | 5.880600 | 5.034904–6.140654 | 6.321484 |
| Original | Fresh YSpec context | 6.284758, 5.321028, 5.362524 | 5.362524 | 5.321028–6.284758 | 5.629956 |
| Confirmation | Cached components | 5.033527, 5.982223, 5.371863 | 5.371863 | 5.033527–5.982223 | 5.830511 |
| Confirmation | Legacy modal assembly | 6.514198, 5.228265, 5.689460 | 5.689460 | 5.228265–6.514198 | 6.058552 |
| Confirmation | Fresh YSpec context | 4.930904, 5.264622, 5.879470 | 5.264622 | 4.930904–5.879470 | 5.510886 |

In the first block the cached median elapsed was 1.36% lower (0.0811 s); in the reverse-order confirmation it was 5.58% lower (0.3176 s). Both are modest, and the repeat ranges overlap. Both block medians were lower, but pair results and spread do not establish a reliable speedup. Confirmation pair reductions were +22.73%, −14.42%, and +5.58%, so the paired signs were inconsistent; repeat ranges also overlap. DSpec cached median elapsed was 1.0966× the YSpec context median in the first block and 1.0204× in the confirmation; those use different arithmetic paths. A read-only comparison of the six fresh YSpec `along.bin`/`across.bin` arrays against the first-block YSpec run found they were not byte-identical; the largest relative L2 difference was 3.08e-16 and the largest absolute difference was 4.73e-30. These YSpec arrays are timing context, not a DSpec acceptance gate.

One profile per route observed these worker-sum scopes:

| Profile scope | Legacy (s) | Cached (s) |
|---|---:|---:|
| Dynamic H assembly | 7.1483 | 2.4827 |
| Receiver-map assembly | 1.7283 | 0.03736 |
| Radial strain construction | 4.4412 | 3.9345 |
| Sparse numerical factorization | 10.5683 | 8.7804 |
| Sparse equilibration and scaling | 3.1310 | 2.3559 |
| Sparse forward solve | 3.1171 | 2.0708 |
| Sparse adjoint solve | 2.4341 | 1.6804 |
| Real-Legendre angular GEMM | 4.9802 | 3.5319 |

The cached route separately spent 0.1473 master seconds building 401 cache entries. Its component-cache payload was 136,612,316 bytes under a 268,435,456-byte (256 MiB) shared-cache cap; worker scratch high-water was 305,052 bytes (32 MiB cap). These are single-profile worker/master observations; worker sums overlap and must not be subtracted from driver elapsed. Other measured scopes and route-level variation also shifted in the same profiles. The full-pipeline timing result remains the end-to-end comparison.

## Build, resources, and evidence

The final timer-disabled and profile-enabled build variants succeeded. An earlier compile attempt failed and was preserved. The build supervisor completed in 161.61 s with a peak of six observed threads and 2.5072 GiB RSS; the independent build review confirmed identical external dependencies to the preceding validated candidate. The numerical producer and the separate checker-only retry completed cleanly. The first checker attempt failed only because its L3 basis-adjoint zero classification assumed the L200 direct-adjoint mode; the corrected checker distinguished the already-selected L200 and L3 validation modes and used the original immutable producer state. No numerical job was rerun for that checker correction.

Maximum numerical producer peaks were nine threads and 2.308 GiB RSS in the initial block, and nine threads and 1.2839 GiB in the timing confirmation. The supervisor uses a user-space watchdog, not kernel aggregate accounting.

- Machine-readable report: [`phase-e-stage2-modal-report-20261103.json`](phase-e-stage2-modal-report-20261103.json), SHA-256 `f8cfd6f589b1637036810b4ee5ec656ca63d3d2fc4e8bcc097162a92bd6fca7d`.
- Final modal validation report: [modal-stage-report.json](../receipts/phase-e-stage2-modal-v1-rebuild3-terminal-v2-snapshot-20261103/phase-e-stage2-modal-terminal-check/data/modal-stage-report.json), SHA-256 `324365882e2c8f86e4a9aab4c82b84c084e14eb6d979ab31a8996a7542c5b28f`.
- Initial producer snapshot: [`../receipts/phase-e-stage2-modal-v1-rebuild3-run-snapshot-20261103/`](../receipts/phase-e-stage2-modal-v1-rebuild3-run-snapshot-20261103/), archive SHA-256 `5c5c9d45a32fc38d54003386d0bed2a3cfb8530131c976f2439a813d7a654108`.
- Timing confirmation snapshot: [`../receipts/phase-e-stage2-modal-v1-rebuild3-timing-confirmation-snapshot-20261103/`](../receipts/phase-e-stage2-modal-v1-rebuild3-timing-confirmation-snapshot-20261103/), archive SHA-256 `24b635848dd4808a4ae78d8bc8ea55d193e49a9df54ce1fc00545b9750883fd1`.
- Build evidence: [`../receipts/phase-e-stage2-modal-v1-rebuild3-build-evidence-20261103.json`](../receipts/phase-e-stage2-modal-v1-rebuild3-build-evidence-20261103.json), SHA-256 `93441e93b1dfd0bccb684bbc24d63365907e123da86fb54022c266bbf83442dc`.
- YSpec repeated section-array comparison (read-only): [`../receipts/phase-e-stage2-modal-v1-rebuild3-yspec-repeat-array-comparison-20261103.json`](../receipts/phase-e-stage2-modal-v1-rebuild3-yspec-repeat-array-comparison-20261103.json), SHA-256 `0b541cb58b42ddbd89e8ceb22a89c9e404f9f38bb66736874393077b7766ef28`.
- Independent source/build review: [`../reviews/phase-e-stage2-modal-rebuild3-normwise-independent-review-20261103.txt`](../reviews/phase-e-stage2-modal-rebuild3-normwise-independent-review-20261103.txt), SHA-256 `7f80a600b27d55e91006adffe5dfdec33a7c268c60f5bffcc9f2d3b7508cd779`.
- Independent build-evidence review: [`../reviews/phase-e-stage2-modal-rebuild3-build-evidence-independent-review-20261103.txt`](../reviews/phase-e-stage2-modal-rebuild3-build-evidence-independent-review-20261103.txt), SHA-256 `40cabce86a555b91db5c014b2b77c8ce1a2fdd47e73548206cb1dc2c79f0759c`.
- Independent producer/checker review: [`../reviews/phase-e-stage2-modal-rebuild3-run-checker-independent-review-20261103.txt`](../reviews/phase-e-stage2-modal-rebuild3-run-checker-independent-review-20261103.txt), SHA-256 `cb00ef34c7ec339344595efd47efa7c8324c549fa7095bc449c299897000dd96`.
- Independent checker-only retry review: [`../reviews/phase-e-stage2-modal-checker-v2-independent-review-20261103.txt`](../reviews/phase-e-stage2-modal-checker-v2-independent-review-20261103.txt), SHA-256 `588938adf24202b70fea502b2cc2b27ff531894efa09c935e0825fc3f234dd59`.
- Independent timing-plan review: [`../reviews/phase-e-stage2-modal-timing-confirmation-independent-review-20261103.txt`](../reviews/phase-e-stage2-modal-timing-confirmation-independent-review-20261103.txt), SHA-256 `8872142cebdc00159ed62e231ea09689ab22e9aec061eaeb2bfe0363ff496aaa`.
- Independent numerical review: [`../reviews/phase-e-stage2-modal-numerical-independent-review-20261103.txt`](../reviews/phase-e-stage2-modal-numerical-independent-review-20261103.txt), SHA-256 `eae3a62224aadfa905ffd67ce27379105087c46fd3b33118ade0e8ef98e0706a`.

This completes the authorized modal-solver substage. No further optimization, production integration, or discrepancy investigation was performed.
