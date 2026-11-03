# Short final Sol review: layout build gate

2026-10-27. **PASS for build only.** Static review of `packing-shared`, `nodal-shared`, `configs/layout-build.json`, and `configs/layout-build-plan.json` against `investigations/tunya-profile-20261027/dspec/source-v6-tunya` and the prior `configs-v3/build.json` profiled DSpec command. No SSH, compilation, tests, or scientific execution performed; no implementation edited.

- Both candidate trees preserve every baseline file byte-for-byte except `engine_headers/Windowed3DFields.h`; their added status notes and `.orig` evidence remain intact. Packing hashes select buckets only: exact label-vector equality prevents collision aliasing. The internal producer emits unique labels in ordered-map `m` order and ascending `l`, supporting the monotone permutation walk. Sorted/unique union keys reproduce the old set order, and absolute-spin thresholds reproduce the old masks. Complex/real coefficient placement and fluid projection agree with the previous mappings; GEMM, phase multiplication, and accumulation loops retain their order.
- Nodal plans cache the same family-specific LTG indices at the same explicit element. Source lift selection and addition, interpolation weights, ascending node summation, and field mapping remain unchanged. Toroidal index construction and use are restricted to the solid elastic range. Element/radius cache keys retain one-sided interface identity; optional r-strain node plans retain that explicit element. Added index payload is included in cache accounting, with existing local/stack fallback preserved. No source-level blocker found for compilation. Exact output equality still requires execution evidence.
- Both manifest compiler argument lists match prior `dsv6-profiled` exactly after substituting only the candidate `physical_driver.cpp` path. Includes, library linkage, flags, dependency paths, and environment are unchanged. Each variant invokes g++ directly; there are no nested build/parallel/LTO commands. The existing reviewed builder serializes variants, requires a fresh empty build root, and fingerprints inputs before/after. The new build root and `acceleration-inputs/layout-v1` paths avoid the saved old input/build locations. The plan contains one bounded build job only, using the existing receipt wrapper and builder; their remote identity/immutability is accepted from the parent handoff, not reverified here.

This PASS authorizes no science run and makes no numerical, parity, timing, or completion claim. The fresh Luna checker HOLD remains outside this build-only decision and must be resolved separately before comparison acceptance; preserve the fixed `2e-13` gate and required byte equality for unchanged-arithmetic outputs.

Reviewed SHA-256:

| File | SHA-256 |
|---|---|
| packing-shared `Windowed3DFields.h` | `e720a9ee3f926e7e5e86c8985ded28692622cf97b46d4752173bee33579550ac` |
| nodal-shared `Windowed3DFields.h` | `87ac1416c66adfdef0fcf641343410847932ccf5a3eca7891342bdf973ad6b45` |
| `layout-build.json` | `fb382c8e1e44221c22956a3f958e09ff930be771201e70b68147395100afcf49` |
| `layout-build-plan.json` | `0420c008ee4bc3e531d62442d41a43b36f14f763900be85ccbb654fafa7f9d11` |
