# Independent v1→v2 review: frozen global GSH adapter

Scope: reviewed only the frozen v1-to-v2 source delta, v2 build/freeze records, and the parent-saved L=12 gate artifacts. I did not compile or run anything. The v2 source and both binary manifests validate; all profiled gates use the profiled v2 binary hash, and the native side-175 gate uses the native v2 hash.

## Delta review

The v2 change is narrow and addresses both v1 fluid issues:

- It records whether each captured source/adjoint pole record belongs to a fluid element, then projects non-radial fluid coefficients before rotation: output component 0 receives `rr + horizontalTrace`; components 1–5 are zero. The reserved radial l=0 contribution remains separate and unchanged. This matches the accepted tensor route's fluid coefficient packing (`Windowed3DFields.h:1761-1763`) and radial l=0 scatter (`:1824-1844`).
- Relative-L2 output now writes JSON `null` for non-finite values and adds `channel_relative_l2_defined`, avoiding v1's invalid `inf` tokens. The fluid projection is described in result metadata.

I found no code-level blocker in this delta. The side-element fluid classification is taken from the captured record's element, which is the selected native one-sided SEM side.

## Saved L=12 gates

Every listed run returned 0, its `result.json` parses as strict JSON, and the grid node/weight parity flag is true. The solid-side all-channel relative-L2 errors are near machine precision. Fluid-side outputs now match the accepted representation: projected trace in component 0 and zero in components 1–5.

| Side(s) | Material / interface context | Channel relative-L2 result | Pointwise-kernel result |
|---|---|---|---|
| 118 | solid interior | all 12 channels 2.17e-15–4.15e-15 | max abs error 4.84e-56 |
| 175 | fluid interior | forward trace 3.23e-15; adjoint trace 1.89e-15; other channels 0 | all elastic kernel channels zero; kappa error 1.53e-58 |
| 172 / 173 | solid/fluid sides at shared 6,368,000 m interface | solid 2.46e-15–4.28e-15; fluid trace 3.74e-15 forward and 2.50e-15 adjoint, others 0 | max abs error 1.39e-54 solid; kappa error 3.68e-58 fluid |
| 154 / 155 | adjacent solid sides at shared 6,321,000 m interface | all 12 channels 2.24e-15–3.97e-15 | max abs error 6.93e-55 |
| 22 / 23 | solid/fluid sides at shared 1,221,500 m interface | solid 2.10e-16–1.74e-15; fluid trace 6.93e-16 forward and 4.80e-16 adjoint, others 0 | max abs error 2.07e-59 solid; kappa error 2.64e-59 fluid |

The native and profiled v2 side-175 runs have identical GSH-field and reference-file hashes, supporting profile/native parity for that fluid gate. These are pointwise one-frequency/one-radial-side comparisons on the tested complete angular grid. They do not establish full physical convergence or performance benefit.

## Metric interpretation caveat

For fluid channels whose reference norm is exactly zero and whose candidate is also exactly zero, v2 emits relative error 0 and marks it defined. Numerically the field arrays agree, and the recorded absolute errors are exactly zero, so this does not undermine the parity result. Mathematically, however, a relative-L2 ratio with zero reference norm is undefined even when the numerator is also zero. `channel_relative_l2_defined` currently means the computed quotient is finite, rather than that its denominator is nonzero; either rename/document that field or set it false for zero-reference channels. The strict-JSON correction itself works for all inspected artifacts.
