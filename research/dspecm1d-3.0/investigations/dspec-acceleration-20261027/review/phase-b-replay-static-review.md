# Phase B parent static review

Scope: source–receiver great-circle sections L384/L1536 at4L+1angles and179positive native SEM sides; depth slices on sides147,130,114. No YSpec, full-order Stored Wigner table, radial solves or full global volume.

The isolated adapter includes the unchanged global-v8 physical driver for SEM setup. It verifies native element/node/radius/material identities against each saved spool. Section synthesis uses the same local-frame coordinates, real Wigner degree bases, exp(+im phi), Eigen internal multiplication and spin rotations as native tensor synthesis. It adds separate L0 terms and applies the saved complex adjoint mixing once, then calls the unchanged pointKernel and sums matching-frequency products. The exact-zero guard bin is skipped; DC retained. Narrow-order section bases depend only on geometry, never frequency/radius. Fluid pointKernel uses divergence; it is equivalent to the native collapsed spin-zero form.

Numerical replay is gated by saved native comparisons atL3/L1536 and byte-exact Generated slice replay atL384. The new roundoff sanity ceiling is2e-12; original2e-13 flags remain visible and historical gates unchanged. No degree-convergence assertion.

Resources: one sequential supervisor, cap40 including parent, actual12section workers with Eigen/OpenMP1; slice1outerworker/GSH32 and existing checked hard OpenMP pool release. Long jobs nohup with unchanged identity-aware userspace watchdog. No overlapping campaign queues. Section output<=128MiB, slice output<=4GiB, working budget64GiB. Rendering uses24GiB RSS watchdog, native-sample overview stride, saved unmasked raw data and explicit normalization/masks. Supervisor manifests hash all declared past/future inputs at runtime, before/after.

Build01 compile overload failure preserved; v2 build completed cleanly75s. No production checkout edits. Full arithmetic controls required before production launch.
