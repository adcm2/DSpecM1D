# Angular-v1 handover — 2026-11-03

## State

- No angular-v1 numerical jobs have been staged or launched. Preserve this state and all existing artifacts.
- The isolated build completed and independently passed. Build evidence review: `reviews/angular-build-v1-evidence-independent-review-20261103.txt`, SHA-256 `e0842d70542f6e3240105a8e6f357a2abe86735ff057b70879496c55fc3d2290`.
- The full compiled-binary SHA-256 is recorded in the build evidence and root restart receipt. Build source/config gate: `reviews/angular-source-build-v1-independent-review-20261103.txt`, SHA-256 `2d3f1cef310ffe4239838e6cc47046e0a7945e3f9084ad73c00e44d5110259d9`.
- Current DSpec sources remain frozen at the reviewed hashes: `phase_d_driver.cpp` `6f30aa9818d82dfaa50d4566a6f1e88444532ec7ccff4d552c1d73a2fd687988`; `phase_d_batched_reconstruction.hpp` `2b0f36e65fb7d81e41ee686013f91b10b6ff11ebd08cda322002d6179fb5fdbb`.
- Build root and run resource admission are root-owned. The shared 40 GB resource reservation was admitted; do not create a new reservation or run work outside the reviewed angular service.

## Current local artifacts

- Producer plan: `configs/angular-run-v1-plan.json`, SHA-256 `33ebc8c0e9012a3239d78c7d9271bf985488b7d23ac75de03f4d8baf4e1d0207`; it declares the 20-thread supervisor cap and eight jobs: `angular-l400-dspec-nq6`, `angular-l400-yspec-nq6`, `angular-l600-dspec-nq6`, `angular-l600-yspec-nq6`, `angular-l800-dspec-nq6`, `angular-l800-yspec-nq6`, `angular-l800-dspec-nq8`, `angular-l800-dspec-nq10`.
- Checker: `tools/check_angular_v1.py`, current SHA-256 `9fff18710bc7b221f389d422f045a5d887f28ce92b34e515f7c2d5b13995ffe7`.
- Standalone checker plan v2: `configs/angular-terminal-check-v2-plan.json`, SHA-256 `584f9c947da3aa9a26dd1285c9664e44725b2a0bdeca41884565fad322535f62`, 241 canonical inputs; differs only in the two corrected build-input paths.
- Reference manifest: `reference/reference-manifest.json`, SHA-256 `e615d81d1d40db1ca6db23f36f56018759b2599739e3e8b4661765fb2f201428`.
- Fixed probes: `data/study1-fixed-probes-v2.csv` and its JSON receipt; exact 168-point physical probe set.
- Reviewed design: `reviews/angular-plan-feasibility-independent-review-20261103.txt`, SHA-256 `1bc147d24b519ba20db8e129375ee3ec6de2e9fa165dd571051e0e29c688c618`.
- Every checker edit through change 16 is recorded immediately in candidate, angular-v1, and campaign `implementation_status.md` files.

## Current gates and payloads

- Run/checker gate: `reviews/angular-run-checker-v1-independent-review-20261103.txt`, SHA-256 `774cceeeb4d7ce541a0c5a9b229b1a22570ae3976bf4afa0ae6b541580f9ba9c`; it binds the exact producer plan, terminal-v2 plan, checker, helper, and reference-manifest hashes above.
- Producer payload-v2: `receipts/angular-producer-v2.tar.gz` (90 members, SHA-256 `ee40cc3e235ddbe1b7bf977c4cb05eb3537df10878b748f75925ca1269e63966`); manifest `receipts/angular-producer-v2-manifest.json` SHA-256 `0dac99e59b78237bc9400d65853ef2ad53f5836857649a54202f5969bc994a35`.
- Validation payload-v4: `receipts/angular-validation-v4.tar.gz` (3 members, SHA-256 `96e1e652ec16421b88a920404cdb19003c3a75cf16f31d1db486a5268846d97a`); manifest `receipts/angular-validation-v4-manifest.json` SHA-256 `13a699eb846f3ae6b2e01cf8818baeb221a06fa6a0d125b72f1cbb071a80735b`.
- Sol verified producer-v2 and the terminal-v2 path delta. A root-owned stage preflight detected two incorrect build inputs; it stopped before numerical launch, both campaign control roots remain absent, and v1/v2/v3 payloads are preserved. The corrected terminal-v2 plan is sealed; root is reviewing the recovery-stage action before any launch.

## Checker work completed locally

- Comparisons use the higher L cutoff as reference, with complete case objects passed to waveform/P metrics.
- Full frequency/time/taper and section-grid checks include NQ8/10. Taper arrays must match exactly.
- The unweighted 99-sample window stays diagnostic; only the tapered window metric enters the predeclared gate.
- All three L800 NQ6/8/10 radial comparisons are retained in the report without early exit. DSpec certification requires those controls; YSpec angular metrics remain reported separately.
- DSpec and YSpec matched-degree tables/window/P are diagnostics using each backend's own P.
- Node metadata must match `nodes.csv` for every backend; DSpec geometry bytes are independently derived from the reviewed formula. JSON output disallows NaN/Infinity.
- Unit-adjoint bins 7, 8, and 41 are explicitly skipped placeholders; only the 32 active bins enter unit-product tail metrics.
- `python3 -m py_compile` passes. A local function exercise using saved L300 data passed complete-case wiring, higher-reference denominator, and exact taper mismatch rejection. These are not angular convergence results.
- Read-only result collector is prepared at `operations/retrieve_angular_producer.py`, SHA-256 `00082ab11f4bc9fa737ba2fd76474a9109325a0d1fe3d7d65252a9a914f57f80`; it requires root's completed finalization receipt and has not been run.

## Pending review and next steps

1. Await Sol's validation-v4 byte check and root's recovery-stage operator review. Root owns the staging/launch operators and kernel service hook; do not create overlapping operators.
2. Reuse the verified producer-v2 extraction and stage only validation-v4 into the already-correct checker stage, using root's reviewed exclusive recovery action. Launch only the sealed eight-job plan after fresh quiet/resource guards.
3. Proposed roots: producer candidate stage `/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya/scientific-inputs/dspec-science-angular-v1`; validation stage `.../scientific-inputs/dspec-science-angular-v1-validation`; campaign control `.../dspec-science-angular-v1-control`; task `.../dspec-science-angular-v1-control/continuation/dspec-science-angular-v1-run-v1`; service `dspec-science-angular-v1-run.service`. Checker stage/control/task: `.../scientific-inputs/dspec-science-angular-v1-check`, `.../dspec-science-angular-v1-check-control`, `.../dspec-science-angular-v1-check-control/continuation/dspec-science-angular-v1-check-v1`.
4. Run only the sealed eight-job plan within the admitted angular service and resource cap. Finalize an immutable producer snapshot, then run the standalone checker as a separate supervised task. Preserve any failed attempt and do not broaden beyond L800 or the frozen case.
5. Sol's exact checker/run gate and root's build/service evidence are required before any numeric launch. No cutoff or 0.1% claim is established yet.

## FINAL — prior pending steps superseded

Eight physical jobs plus standalone checker completed and retrieved; AFTER-results282595ff and report-claims15cfbfed independently sealed. Final report e38a87cc delivered with manifest; no0.1% cutoff certified, conditional extension false. All services inactive, logical40GB released474455aa. Do not restage/relaunch/check or extend. No further source implementation authorized in this stage.
