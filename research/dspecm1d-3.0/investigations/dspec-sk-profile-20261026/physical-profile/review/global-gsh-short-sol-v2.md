# Short independent Sol check: frozen global GSH v2

26 October 2026. Read-only mathematical/source delta and saved-binary check. I compared frozen v1 and v2 `physical_global_gsh.cpp`, the accepted frozen v4 tensor projection, and all nine parent-saved L=12 native-grid gates. I did not compile or run a solver or change the prototype.

**Decision:** the narrow v2 change fixes v1's fluid **field-representation** mismatch on the tested one-sided native radial sides. The Wigner coefficient rotation itself is unchanged, and the saved solid result is byte-for-byte unchanged. The complete-grid L=12 complex-field and pointwise-kernel gates pass at solid interiors, fluid interiors, and three duplicated-radius interface pairs. This is an angular-equivalence result at one frequency and low degree, not a physical convergence or high-degree performance result.

The v1→v2 numerical delta is exactly where needed: each captured nonradial source/adjoint pole record is tagged by `sem->elementIsFluid(record.element)`; in fluid, coefficient channel 0 receives `raw rr + raw horizontalTrace` and channels 1–5 receive zero **before** the existing unnormalized-small-d rotation. Solid coefficients keep their original six values. The separate radial `l=0` rr/trace scatter is unchanged. This matches frozen v4 `sampleForKernelsTensor`'s pressure/divergence packing before synthesis and its separate l=0 addition. Thus the change preserves the correct `d^l_{m_local,m_global}` orientation and global spin transform established in the v1 solid gates. The v2 solid side-118 candidate and reference binary files are each byte-identical to their v1 counterparts.

I independently read each strict-JSON result and both channel-major complex128 files, recomputed all 12 channel errors at every one of the 1,225 native angles per side, and checked the saved pointwise kernel maxima. The nine runs cover solid side 118, fluid side 175 in both native and profiled builds, and three exact duplicated-radius pairs: 22/23 at 1,221,500 m (solid/fluid), 154/155 at 6,321,000 m (solid/solid), and 172/173 at 6,368,000 m (solid/fluid). The pairs retain distinct element/node identities at the same radius; matching each side separately is material to the one-sided interface contract.

| Side(s) | Material | Largest nonzero-reference channel relative L2 | Largest complex field absolute error | Largest pointwise kernel absolute error |
|---:|---|---:|---:|---:|
| 118 | solid | 4.15e-15 | 1.42e-26 | 4.84e-56 |
| 175, native and profiled | fluid | 3.23e-15 | 8.09e-28 | 1.53e-58 |
| 22 / 23 | solid / fluid | 1.74e-15 / 6.93e-16 | 3.18e-29 / 2.92e-29 | 2.07e-59 / 2.64e-59 |
| 154 / 155 | solid / solid | 3.92e-15 / 3.97e-15 | 1.23e-25 / 1.18e-25 | 6.93e-55 / 6.53e-55 |
| 172 / 173 | solid / fluid | 4.28e-15 / 3.74e-15 | 1.15e-25 / 1.25e-27 | 1.39e-54 / 3.68e-58 |

For every fluid gate, the eight forward/adjoint spin ±1/±2 channels have **exactly zero** candidate and reference arrays and exactly zero absolute error. The two spin-zero channels retain the separate radial-l0 contribution; the nonradial pressure projection gives finite near-roundoff relative errors. At side 175, native and profiled `gsh_fields.bin` and `reference.bin` are byte-identical. The v2 `result.json` files all parse strictly; the v1 bare `inf` failure is removed.

One metric name remains imprecise: `channel_relative_l2_defined` is true for the zero/zero fluid shear channels because the implementation floors the reference denominator before division, producing finite numeric zero. The arrays agree exactly, but the mathematical relative norm is undefined when the reference norm is zero. Read those entries as **absolute-zero parity**, not defined relative errors; a future metadata-only clarification could mark them false or describe the denominator-floor convention. This does not change the sampled fields or kernels.

No inference about high-L runtime follows from these low-L gates. In particular, a suspected denormal/underflow contribution to high-degree cost is a separate diagnostic hypothesis; v2 neither measures it nor changes floating-point behavior. Full 179-side/multi-frequency output, RSS/I/O feasibility, and waveform/kernel physical convergence remain separate parent-owned gates.

Source/evidence anchors: `global-gsh/source-v1-frozen/physical_global_gsh.cpp`, `global-gsh/source-v2-frozen/physical_global_gsh.cpp` (fluid tag and coefficient projection near lines 135–180, metadata/JSON near 401–415); frozen DSpec v4 `Windowed3DFields.h` fluid tensor packing near 1757–1794 and radial l0 scatter near 1821–1844; `runs/gsh-v2-gate-l12-r*-*/data/{result.json,gsh_fields.bin,reference.bin}`.
