# Phase B fixed-mesh v3 run-plan review

Independent review: Sol (`phase_b_plan_review`), 2026-10-31.

Disposition: PASS for the isolated B run queue and checker after the final path, profile-schema, per-solve cut-radius, fluid-family classification, and reference-file corrections.

Scope reviewed: four NQ5/NQ6 nominal and h/2 meshes; sequential own-weight waveform controls, full-domain frozen-weight controls, and paired principal-restriction controls; exact A1 common-side packets; frozen NQ10 A4 components and A2 packed-reference kernel products; direct polynomial radial derivative path (`DSFAST_RSTRAIN_BRIDGE=0`); eight frequency workers, inner thread limits one, 16 GiB solve cap, 20 aggregate thread cap; checker requires actual restricted H^H work and omitted DOFs and classifies common samples affected by relevant family/degree cut radii.

Interpretation limits: no continuous-accuracy threshold is imposed. Full controls use the inherited response cache while restricted controls bypass it, so no truncation speedup claim is supported. The `dspec-forward` own-weight control also retains legacy unit-adjoint capture; it is a waveform/measurement control, not a pure-forward timing. The restriction boundary is a principal matrix with eliminated cross-cut entries and zero-extension output, not a physical boundary condition.
