# DSpecM1D 3.0 — isolated backup checkpoint

This branch backs up the current development state before a separate, deliberate PR review. It is not a release or production integration. The user's 3 November 2026 instruction authorizes this commit and push only; historical no-commit/no-push statements in the preserved campaign records describe the earlier execution boundaries.

## What is preserved

The branch starts from `a778b13b22dca9709e1f96604b62925707654e5d` on `feature/radial-sensitivity-kernels-k1`. The five previously modified tracked files and seven existing `validation/phase0` files were copied verbatim from the original DSpecM1D working tree. They were already present before this backup; no new optimization was integrated into the library during backup preparation.

The `investigations/` folder here preserves the recent isolated acceleration and scientific-validation candidates, their plans/checkers, independent reviews, execution receipts, saved numerical outputs, reports and figures. It also contains selected earlier profiling source, scripts and metadata. `backup-manifest.json` records every copied file's original location, byte count and SHA-256, plus the explicitly excluded files. Large earlier profiling arrays and transfer bundles are still at their original local/Tunya locations: this Git branch does not back those up. Existing provenance records retain the available hashes for those external artifacts. Compiler outputs, caches and temporary PDF renders were excluded.

## Where to begin the later review

- Acceleration history and status: `investigations/dspec-acceleration-20261027/implementation_status.md`.
- Final modal-optimization candidate: `investigations/dspec-acceleration-20261027/candidates/phase-e-stage2-modal-v1/dspec/`.
- Scientific Studies 1 and 2: `investigations/dspec-scientific-validation-20261103/plan.md`, its `implementation_status.md`, reviews and reports.
- Latest high-degree candidate: `investigations/dspec-scientific-validation-20261103/angular-v1/candidate/dspec/`; the compiled driver entry was `phase_d_driver.cpp`.
- Angular results: `investigations/dspec-scientific-validation-20261103/angular-v1/output/pdf/dspec-angular-convergence-20261103.pdf` and `reports/delivery-manifest.json` under that same `angular-v1` directory.
- Independent numerical findings: `angular-v1/reviews/angular-after-results-independent-review-20261103.txt`; final report claim seal: `angular-v1/reviews/angular-report-claims-independent-review-20261103.txt`.

The experiments retain the batched radial/strain and angular reconstruction and the reviewed modal-solver work. Candidate implementations are still under research directories. Their historical build configurations point to pinned Tunya source/dependency/build locations; this is not a standalone binary or dependency distribution. Read the saved build receipts before attempting a new build, rather than treating the historical orchestration state as a live queue.

## Scientific status to retain

For the 2/4/8/10 mHz four-corner filter, 1000–1200 s measurement window and 40 s taper, the L400/L600/L800 ladder did not establish the 0.1% angular convergence target. Both codes show substantial shallow-kernel cutoff dependence. Fixed-L800 DSpec–YSpec agreement improves from about 0.788% at NQ6 to 0.00104% at NQ10 on the selected 89-probe kernel vector; that is not angular convergence. NQ8→NQ10 still fails the required 50 km stratum at about 0.187%, so the conditional extra ladder was not run.

The user accepts the current limitation for now. A limitation of the truncated method is a working interpretation, not a proof that every implementation issue has been ruled out. No new H1 formulation or singularity regularization was implemented. Distinguish the actual transverse-velocity-traveltime fixture from the older displacement-energy fixture.

All authorized numerical jobs finished. The logical 40 GB Tunya reservation was released through CCDSM, with DSpecM3D informed. Preserved empty slice registration is not admission for a future campaign. No new numerical work, PR, merge, deployment or production integration is part of this backup.

## Backup verification

`backup-manifest.json` binds the copied source/evidence; `backup-verification.json` records the final local integrity checks. The original working copy and investigations were not edited. This operation verifies snapshot and Git integrity; it does not rerun or extend the scientific acceptance tests.
