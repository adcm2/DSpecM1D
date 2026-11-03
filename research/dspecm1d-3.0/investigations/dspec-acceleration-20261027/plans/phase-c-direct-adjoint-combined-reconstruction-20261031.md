# Phase C plan: direct adjoint and combined reconstruction

## Gate and scope

Phase B numerical and result review passed. Use the accepted NQ10 mesh with its original `h=0.025` and full radial domain. Phase B's principal restriction remains rejected for use: its paired controls are descriptive and no accuracy budget was specified. All edits stay in `candidates/phase-c-direct-adjoint-v1`; all builds and numerical work run on Tunya through the existing detached supervisor. Preserve Phase B source, outputs, and failed checker attempt.

The approved sequence is C source implementation and independent review, fresh profile build, focused equivalence proof, then native/profile complete-pipeline matched validation. Do not integrate into production or investigate field discrepancies.

## C implementation

1. Add a direct measurement-direction Hermitian-adjoint path. At each degree/frequency, form the same receiver frame and `localWeights = conj(receiverAtPole) * conj(fz) * frozenWeights` used by the accepted route. Solve the receiver-adjoint RHS directly, and compare its solutions and complex synthesized modes against the existing QR/unit-basis projection at the same operator. Preserve existing source/receiver coefficient conventions and isolated radial `l=0` adjoint treatment.
2. Integrate the accepted Phase A degree-inner reconstruction into the isolated combined candidate driver. Stream `visitProjectedRadialModes` into a worker-local contiguous `phasea::StrainTable`, use unchanged `phasea::strainOf` / slot assignment, `LegendreTable` / `evaluate`, frame rotation, and inherited fluid projection on the declared native one-sided radii. Keep the current tensor reconstruction available as an equivalence oracle; compare complex fields and products before enabling the combined path in the pipeline.
   The angular plans are immutable and worker-local, keyed by role, receiver, and angle. Reuse is allowed only while degree, input angles, and exact source/receiver frame axes match; radial side data is not retained in the angular cache.
3. Use worker-bounded two-pass reconstruction/recomputation. Keep all 35 forward bins. Skip kernel work only for bins 7, 8, and 41 when the measurement-adjoint weight/forcing product is exactly zero before kernel work. Assert that the exact-zero bins are precisely 7, 8, and 41; use no threshold or a posteriori kernel cancellation test. Include an on/off equivalence proof for those skips; do not skip any other bins.
4. Keep the accepted solver/backend and existing 8 total frequency workers, one inner thread, current aggregate thread and memory caps. Profile and native variants are isolated under fresh Tunya build roots.

## Validation sequence

- Source review by Sol before build.
- Run the existing frozen L3 audit as `L=3`, `NQ=8`, `h=0.05`, audit quadrature, with one outer frequency worker, one inner solve thread, direct measurement adjoint enabled, and common reconstruction enabled. The candidate checker compares candidate baseline trace/spectrum/components with the exact frozen Stage 2 base and checks the candidate integrated kernel against the existing frozen plus/minus traces under the existing `2e-4` central-FD limit. No new perturbed solves are performed.
- Profile build and focused equivalence run on Tunya. Require direct RHS/solution/mode equivalence, combined reconstruction fields/products against the retained tensor route, and zero-bin skip on/off equality for every saved nonzero output.
- After that review, run native and profile complete-pipeline controls on the frozen Stage 2 workload, with 8 total workers and inner thread count 1. Compare against the accepted Stage 2 native median of 16.8176 seconds and the accepted isolated Phase A references. Report all component times, total wall, output hashes, raw unmasked sides and source/receiver/bulk limits, time-shift and gravity controls, finite-difference control where supported, RSS, threads, and exact skips. No YSpec rerun.
- If any equivalence or resource gate fails, preserve the artifact and stop the affected stage. Use full-domain NQ10 throughout.

## Evidence and interpretation

Retain the compact Phase B metrics and remote hash receipts. The 49 MB checker JSON and raw CSV remain remote. The B common-side relative differences are descriptive because no accuracy tolerance was supplied; no reduced mesh or radial truncation is accepted. Report implementation, saved-artifact validation, and scientific conclusions separately.
