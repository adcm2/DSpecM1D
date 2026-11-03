# Study execution configuration status

2026-11-03: Created `configs/study1-build-v1.json`. Native/profile solver build configuration cloned from final modal build; only fresh candidate/stage/build roots change. All twelve shared external dependency paths and compiler flags retained. SHA-256 `b5372290f075a451dfecf97f9edeef25441b82870f8c3354dc9eec001417fda8`. Pending independent source/build review.

2026-11-03: Created `configs/study1-build-v1-plan.json`. Detached sequential build supervisor plan cloned to a fresh study task/control root, 8 threads and 32 GiB. No remote action has occurred. SHA-256 `31eb960a56cb3cadb1f11e90e864848b8d1f99c2ed90e85c4ac3c6c2a9fc9ef2`. Pending independent source/build review.

2026-11-03: Added build_stage.py for exact include-closure packaging, pinned dependency preflight, exclusive staging and detached build launch. Every network action uses the established SSH helper; scripts compile before transfer, failed operations preserve receipts, no remote action yet. SHA-256 `168fcbb72b96944f295d8e204d3bfd4d2a853c8197a84ceed1ccd08460fa7365`.

2026-11-03: Tightened staging order: the complete generated extraction/verification script now compiles before creating the remote stage or transferring payloads. No remote action has occurred. Operator SHA-256 `6d86f30da2649753dfcc345199f73e57e0e15843cecd7e5d925aaf61421c6446`.

2026-11-03: Build operation `preflight` completed; receipt saved. No solver job launched by this operation.

2026-11-03: Exact source-hash package guard rejected a metadata-only driver change made after the initial source seal. It stopped before archive/manifest creation and before remote actions. Guard receipt preserved; awaiting updated focused source seal.

2026-11-03: Prepared exact reviewed build allowlist, 23 files; archive SHA-256 9377828fa005d45e9aa9a9d73deba04fe674c13349742bf979a27f9d1353b302. No remote action.

2026-11-03: Build operation `stage` completed; receipt saved. No solver job launched by this operation.

2026-11-03: Build operation `launch` completed; receipt saved. Detached build launched.

2026-11-03: Retrieved successful terminal build evidence, compiler command/input maps, observed binary hashes and supervisor cleanup/resource state. Evidence SHA-256 `976951d9afbd9afc7aaeca73c7247d9bd51967b7e15b19e5e0a92ad648f408bf`. Independent build-evidence review pending.

2026-11-03 Angular-v1 source change 1: fresh candidate `angular-v1/candidate/dspec/phase_d_batched_reconstruction.hpp` now has explicit memory-policy plumbing with unchanged 16/32 MiB defaults. No Tunya action or numerical run; full-angle 32/64 MiB policy remains to be wired only to the reviewed angular-v1 workload. Candidate status and angular-v1 status were updated immediately.

2026-11-03 Angular-v1 source change 2: completed explicit policy propagation through plan preparation and visitor checks; added the known Wigner recursion-table transient to the peak formula while leaving default limits/math unchanged. Candidate and angular-v1 logs updated immediately; no build/run or remote action.

2026-11-03 Angular-v1 source change 3: fresh `phase_d_driver.cpp` now recognizes only the explicit angular-v1 mode and exact D/Y L400/L600/L800 and DSpec L800 NQ8/NQ10 workloads, while leaving prior degree guardrails unchanged by default. Candidate and angular-v1 logs updated immediately; no build/run or remote action.

2026-11-03 Angular-v1 source change 4: propagated the explicit policy through the batched sample wrapper with legacy defaults. Candidate and angular-v1 logs updated immediately; no build/run or remote action.

2026-11-03 Angular-v1 source change 5: new angular-v1 mode now selects 32/64 MiB caps only at its explicit visitor call sites; existing paths retain 16/32 MiB. Candidate and angular-v1 logs updated immediately; metadata/oracle guard work remains, no build/run or remote action.

2026-11-03 Angular-v1 source change 6: result metadata exposes angular-v1 mode, individual per-worker caps, Wigner recursion transient, and calculated cache/temp peaks. Candidate and angular-v1 logs updated immediately; no build/run or remote action.

2026-11-03 Angular-v1 source change 7: corrected the newly added JSON output stream quote; no scientific behavior changed. Candidate and angular-v1 logs updated immediately; no build/run or remote action.

2026-11-03 Angular-v1 source change 8: expanded geometry/temp policy is now DSpec-only; new DSpec jobs require explicit real-Legendre, batched direct-index radial, and Dref routes, while YSpec remains unchanged. Candidate and angular-v1 logs updated immediately; no build/run.

2026-11-03 Angular-v1 source change 9: result metadata now reports expanded-cap use explicitly and emits zero geometry/temp caps for YSpec. Candidate and angular-v1 logs updated immediately; no build/run.

2026-11-03 Angular-v1 source change 10: expanded-cap DSpec runs fail closed if extra reconstruction/common-field oracle buffers are requested. Candidate and angular-v1 logs updated immediately; no build/run.

2026-11-03 Angular-v1 source change 11: result JSON now names the exact filter/bins/window/taper/component/weight source. Candidate and angular-v1 logs updated immediately; no build/run.

2026-11-03 Angular-v1 source change 12: metadata now identifies the 64 MiB cap and observed peak as cache-plus-temporary totals, with legacy fields retained. Candidate and angular-v1 logs updated immediately; no build/run.

2026-11-03: Source change 13 corrects the angular-v1 route guard to accept the actual `DSFAST_RADIAL_ROUTE=batched` control and adds pre-allocation dimension checks (L<=800, exactly 68+10 angles, <=350 radial sides) so the expanded-mode memory arithmetic stays inside reviewed bounds. No build or numerical run.

2026-11-03: Source change 14 makes the retained `phase_d_batched_cache_max_bytes` output alias report the policy combined cache-plus-temporary cap (32 MiB default, 64 MiB angular-v1) instead of a hardcoded value; the separate geometry-cap field remains explicit. No build or numerical run.

2026-11-03: Source change 15 fail-closes the expanded angular-v1 DSpec route if radial/modal oracles are enabled or the accepted legacy modal route is overridden. This keeps unbudgeted oracle buffers and alternate modal allocation out of the declared memory envelope. No build or numerical run.

2026-11-03: Prepared hash-bound angular-v1 inputs by copying the approved 168-point probe set, generator, tight YSpec case and the sealed L200/L300 DSpec/YSpec reference outputs into the isolated angular-v1 tree. Reference manifest `reference/reference-manifest.json` SHA-256 `360dffb1c6abfa756cf82a816fadba0fbb8e2ebf9401ee045bdda53e57b19c36` records 68 copied files and binds the Study2 after-results review. No source or numerical changes.

2026-11-03: Prepared the isolated YSpec tight case input with its model path bound to the angular-v1 staged PREM deck; all physical fields and 1e-7 solver tolerance are copied unchanged. Case SHA-256 `f7bdca24fc9e657b7b4b4b9dfea8441e28ad8e6fad2097cb4444177fda2030ad`. No source or numerical changes.

2026-11-03: Drafted the eight-job angular producer plan with D/Y L400, L600, L800 and DSpec L800 NQ8/NQ10 controls. It pins the fixed absolute h, 68+10 output angles, real-Legendre/batched routes, own-waveform weights, per-job resource caps and exact service campaign root; stale inherited inputs were removed before review. Current plan SHA-256 `7ecc63b1807b3ef5644122b1544c458590f31558909b6559ceabdb60e07cf28d`. This is a local draft pending source/build and independent run/checker review; no launch.

2026-11-03: Run-plan input closure was tightened after local path inspection: all eight jobs now reference only their exact angular-v1 case/model/probe, angular build binary/command, runtime utilities, and the pinned YSpec DSO/command. Removed stale Study1/Study2 candidate paths and recorded the queue as an unsealed draft. Current plan SHA-256 `09ddf6ccd34431eb6126c582d0ef9bc94809c2f9a53e2e2e5a6e3515d88a7f41`. No remote action.

2026-11-03: Added `tools/check_angular_v1.py` for the exact eight-job queue: state/output integrity, L400/L600/L800 waveform/normalization/P/kernel comparisons against saved L200/L300 and L800, far-solid depth strata with exact-zero handling, full section-grid and high-L radial controls, expanded memory guards, and DSpec P800 counterfactuals. It explicitly excludes source/contact pointwise gates and the unsupported full 3D high-L action. Checker is draft pending static review and source/build/run gates; no remote execution.

2026-11-03: Checker change 2 adds an explicit 168×35 per-probe/bin unit-product reader with exact zero-bin validation and own-P binding; it no longer calls a helper that does not expose that API. Python syntax compilation passes; no scientific runs.

2026-11-03: Added the previously reviewed `check_study1_helpers.py` beside the angular checker for the minimal validation package; SHA-256 `610271198745973b654f2154395bde7aac7b034ca39cd47b0238837d3a982c17`, byte-identical to the sealed Study1/Study2 helper. No source or run changes.

2026-11-03: Checker/reference change 3 binds the L200 DSpec comparison to Study2’s hash-qualified sparse-grid bridge (68+10), not the older full-grid L200 section capture. Both retain the same fixed probes and own waveform/kernel outputs; L300 is the 68+10 section-grid baseline. Added the bridge outputs to the reference manifest without replacing prior copied references. No solver run.

2026-11-03: Updated the reference manifest after the YSpec case file was rebound to the angular-v1 staged PREM path, refreshing exact size/hash entries for the copied study inputs. Manifest SHA-256 is now `8aa9ce749223ae28549b019f13ff93ee13134c47a56df310da2c05c8b317050f`. No numerical files were modified.

2026-11-03: Checker change 4 binds angular-v1 mode, DSpec-vs-YSpec absolute-h/batched-angular metadata and the 64/32 MiB cap alias; cutoff certification now requires two successive increments for L400, while L600 is explicitly uncertifiable with only one higher increment. Required kernel cohorts gate both L2 and max-absolute ratios. Python syntax compilation passes; no scientific run.

2026-11-03: Checker change 5 adds per-bin receiver-spectrum and fixed-probe unit-product tail records plus raw along/across section-vector difference diagnostics. The report labels spectrum increments additive by cutoff and unit-product differences as complete products with cross terms, not isolated degree blocks. No science run.

2026-11-03: Checker change 6 requires NQ6 L400/L600/L800 outputs to preserve the exact L300 fixed-h radial nodes and element-boundary files before section-vector comparisons. The checker retains output data paths internally for tail/section diagnostics; no new solver data is introduced.

2026-11-03: Checker change 7 emits frequency-resolved receiver-spectrum and unit-product increments and along/across section diagnostics for adjacent cutoffs and each level versus L800. It now fails the checker if L800 NQ8 or NQ10 exceeds the predeclared 0.1% radial control on required waveform, normalization, P, far-solid and depth-stratum kernel measures. Python syntax compilation passes; no numerical execution.

2026-11-03: Checker change 8 independently recomputes exact generic/real visitor workspace and full-angle product-byte formulas from result degree, side and angle counts; the reported combined peak and product buffer must match the formula before cap acceptance. This verifies bounded arithmetic rather than trusting cap metadata alone. Python syntax compilation passes; no C++ or numerical build.

2026-11-03: Added the hash-bound baseline angle and mesh tables required to independently recheck fixed-probe geometry provenance. The reference manifest now includes both CSVs and has SHA-256 `e615d81d1d40db1ca6db23f36f56018759b2599739e3e8b4661765fb2f201428`. These are copied from the sealed Study2 inputs; no numerical content changed.

2026-11-03: Checker change 9 adds explicit YSpec case input and validates the 1e-7-only solver-tolerance difference, fixed PREM hash, and probe generator/mesh/path-reference hashes against the frozen receipt. This binds physical geometry provenance independently of producer output hashes. Python syntax compilation passes; no solver run.

2026-11-03: Producer plan change 2 binds the exact reference-manifest SHA-256 `e615d81d1d40db1ca6db23f36f56018759b2599739e3e8b4661765fb2f201428` into the run plan. Plan SHA-256 is now `52cde2cab449ae7b9396fefd9b78d261f2244f2e346f02b6359f587295326084`. No remote stage or launch.

2026-11-03: Checker change 10 validates exact per-job argv/environment/output/input closure and binds the reference manifest plus reviewed angular build plan into the producer plan. The checker and plan pass Python syntax/JSON parsing; hashes are checker `2c9f7ac1a049aa6e563c3f826e7a0a0eb2b1170f433a21c035dad478e3b5bde8`, run plan `52cde2cab449ae7b9396fefd9b78d261f2244f2e346f02b6359f587295326084`, and reference manifest `e615d81d1d40db1ca6db23f36f56018759b2599739e3e8b4661765fb2f201428`. All remain drafts pending independent run/checker review.

2026-11-03: Local pre-run closure review: all eight serialized job commands pass the checker’s exact ID/backend/L/NQ/case/h/workers/routes/argv/output/input invariants; the reference manifest verifies all 88 files and is bound to the sealed Study2 after-results review. This is local plan/evidence integrity only, not a solver run.

2026-11-03: Local checker-data closure check on the hash-bound L200 bridge/L300 DSpec and YSpec references passed: 168 physical IDs, 35 own-P bins, DSpec 5,880 keyed probe/bin products with field/product consistency, and 28 explicitly unavailable YSpec source-trace rows. This inspects saved Study2 artifacts only, not new solves.

2026-11-03: Checker change 11 adds a whole-trace unweighted velocity comparison as a separate diagnostic alongside the two frozen measurement-window metrics; it does not affect 0.1% acceptance. Python syntax compilation passes.

2026-11-03: Checker review packet updated after whole-trace diagnostic addition. Current hashes: checker `5987c0f9aafde1bc8518ad0d9d38dd10f44d2624b399563c40c68202abdd0d10`; run plan `52cde2cab449ae7b9396fefd9b78d261f2244f2e346f02b6359f587295326084`; reference manifest `e615d81d1d40db1ca6db23f36f56018759b2599739e3e8b4661765fb2f201428`. Build has launched under the reviewed admitted service; numerical jobs remain held pending build and run/checker gates.

2026-11-03: Local fixed-probe receipt check confirms the DSpec case/probe CSV hashes and the YSpec case differs only by the staged PREM path plus ODE tolerance1e-7. No candidate or numerical input was changed.
2026-11-03: Angular-v1 checker change 5 independently validates node counts and geometry-byte arithmetic for all DSpec quadrature cases and marks the unweighted window measure diagnostic-only. No source, plan, build, or numerical changes.
2026-11-03: Angular-v1 checker change 6 extends frozen frequency/time/taper invariants to the NQ8/NQ10 controls and validates their section-angle grid. No source, plan, build, or numerical changes.
2026-11-03: Angular-v1 checker change 7 fixes measurement-call wiring and uses the higher cutoff as the reference for all convergence ratios. No source, plan, build, or numerical changes.
2026-11-03: Angular-v1 checker change 8 retains all three higher-reference radial-order comparisons even when a gate fails, removes the unweighted-window diagnostic from acceptance, and adds clearly labeled D/Y and skipped-bin diagnostics. No source, plan, build, or numerical changes.
2026-11-03: Angular-v1 checker change 9 requires every L800 NQ6/8/10 radial comparison to pass before any angular cutoff can be certified; the report still retains all diagnostics on failure. No source, plan, build, or numerical changes.
2026-11-03: Angular-v1 checker change 10 validates node inventory metadata for DSpec and YSpec and corrects diagnostic tail norm denominators to the higher cutoff. No source, plan, build, or numerical changes.
2026-11-03: Angular-v1 checker change 11 verifies the exact 32 active/three skipped unit-product bin inventory for each DSpec result and hash-bound reference. No source, plan, build, or numerical changes.
2026-11-03: Angular-v1 checker change 12 prevents skipped zero-adjoint bins from entering unit-product tail metrics and labels those entries as uncomputed placeholders. No source, plan, build, or numerical changes.
2026-11-03: Angular-v1 checker change 13 fixes reference initialization order in the NQ8/NQ10 invariant checks. No source, plan, build, or numerical changes.
2026-11-03: Angular-v1 checker change 14 removes a redundant reference assignment left by the preceding ordering fix. No source, plan, build, or numerical changes.
2026-11-03: Angular-v1 producer plan now declares thread_cap=20 to match its reviewed supervisor limit; the resource service separately enforces the admitted task/CPU/memory slice. No numerical run.
2026-11-03: Angular-v1 checker change 15 rejects producer plans that omit the exact 20-thread supervisor cap. No source, build, or numerical changes.
2026-11-03: Angular-v1 checker change 16 separates angular cutoff results by backend from the DSpec-specific radial gate. No source, plan, build, or numerical changes.
2026-11-03: Drafted the angular-v1 standalone terminal checker plan with an immutable finalized producer snapshot, explicit completed producer/build dependencies, and the full declared output/reference/input closure. It remains unsealed pending independent review; no remote staging or numerical run.
2026-11-03: Prepared minimal angular producer and validation allowlists using the root-agreed archive/manifest names. Local archive/member/current-byte verification passes; remote staging remains gated on Sol byte review and root operator review.
2026-11-03: Replaced the pending v1 staging layout with preserved v2 allowlists matching the operator roots: 90 producer-stage files and 3 checker-stage files. Both archive manifests and member hashes pass local verification. Root operator constants must use the v2 archive/manifest names. No remote staging or numerical launch.
2026-11-03: Angular-v1 terminal-plan change 1 canonicalizes implicit probe-generator/geometry inputs from `reference/../data/` to the same producer-stage `data/` files, avoiding duplicate noncanonical paths in the input manifest. No remote staging or numerical run.
2026-11-03: Validation-v3 package now binds the canonicalized terminal plan and updated Sol review; v1/v2 payloads remain preserved. Archive SHA-256 `e24830a81345f6ec3e9c232caeb556e6aff28c1a972fcbfbdae679fadd43bed0`; manifest SHA-256 `2e8b449eedf6ab6f629982859a118a33c5c9a606fcbef379340f295c4128af1f`. Await root operator review; no remote staging or numerical run.
2026-11-03: Created versioned terminal-plan v2 after the staging preflight identified two wrong declared build-input paths. Only those two path strings change; the stopped preflight occurred before numerical execution and both control roots remain absent.
2026-11-03: Validation-v4 packages terminal-plan v2 and the unchanged checker/helper for the root-owned recovery staging path. Exact 3-member archive and manifest hashes are recorded in receipts; no remote numerical work.
2026-11-03: The reviewed read-only collector retrieved and verified all 132 angular-v1 output files, 8 logs, and finalized state; archive hash is recorded in the angular-v1 retrieval receipt. No remote write or rerun.
2026-11-03: Saved a compact preliminary angular comparison table from the finalized, hash-verified raw producer bundle. It uses the higher cutoff as reference and is explicitly non-final pending Sol's independent scientific review; no certification claim is recorded.
