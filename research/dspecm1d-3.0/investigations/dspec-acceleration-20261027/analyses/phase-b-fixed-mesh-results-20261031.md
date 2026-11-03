# Phase B fixed-mesh results

**Decision:** Sol's independent result review passed. NQ5/NQ6 fixed-mesh comparisons are descriptive because no numerical accuracy budget was specified. The paired principal radial restriction was executed and measured, then rejected for use. Phase C uses the accepted NQ10 mesh (`h=0.025`) and full radial domain.

## Mesh and waveform controls

| Mesh | Elements / native sides | 35-bin own waveform vs paired refinement | 35-bin waveform vs frozen NQ10 |
|---|---:|---:|---:|
| NQ5 nominal | 70 / 199 | 9.62e-6 | 9.42e-6 |
| NQ5 finer | 132 / 363 | 9.62e-6 | 3.55e-7 |
| NQ6 nominal | 53 / 189 | 3.16e-6 | 3.20e-6 |
| NQ6 finer | 96 / 322 | 3.16e-6 | 1.80e-7 |

All rows contain 35 bins and 1,201 waveform samples. The time-shift adjoint identity ranged from `0.9999999999999997` to `1.0000000000000004`. Frozen NQ10 weights and each row's own waveform weights were kept distinct.

Across common one-sided complex fields/products at bins 9, 25, and 40, worst nominal-to-finer group relative L2 values were: NQ5 forward `6.48e-4`, adjoint `1.11e-3`, product `3.12e-3`; NQ6 forward `5.20e-4`, adjoint `4.37e-4`, product `5.12e-4`. These are descriptive errors, with no pass threshold. Pooled field-relative metrics include small-denominator groups and should not be read as absolute errors.

## Restricted H/Hᴴ control

Each restricted row executed 14,000 kernel-pass Hᴴ solve scopes, including solves that fell back to the full domain. Actual principal restrictions and omitted DOFs were nonzero. The checker classified all 45,122 requested one-sided locations per active bin against exact family/degree cut radii. The largest pooled restricted-minus-full relative differences were forward `8.10e-5`, product `9.87e-5`, and adjoint `0.7354`. The latter was the NQ6-finer bin 9 group `material_interface|fluid=1|outside_all_relevant_cut_radii` (154 samples); its pooled absolute complex error norm was `1.12159` against a full-field norm of `1.52508`. This measured approximation difference has no accepted tolerance, so truncation is not adopted.

Restricted active-bin receiver spectra agreed with their full controls to at most `3.48e-10` relative L2. Residual/duality counters are integer values scaled by `1e-15` and floor/cap their readings; the largest observed adjoint residual counter was `1,884,781` and H/Hᴴ identity counter `2,167,754`. They are quantized audit bounds, not raw full-precision residuals.

The full frozen-weight jobs took `11.35–17.91 s`; restricted jobs took `13.80–22.01 s`. These timings do not measure a truncation speedup: full controls used cached responses, while restricted controls bypassed the cache to guarantee actual restricted solves.

## Provenance and limits

- Numerical queue: `phase-b-fixed-mesh-v1-fixed-mesh-03`, 12 jobs, all successful; 8 frequency workers, one inner thread, 16 GiB/job and 20-thread caps. Maximum observed RSS was about 5.95 GiB and 9 threads.
- Profile build: Tunya binary SHA256 `3de02d1c95545a52223f6f664ba4f9ee1f89b4614e56192f1080c0b2b7ca9b67`.
- Checker-only retry: `phase-b-fixed-mesh-check-01`, completed successfully; 133.15 s, 1.291 GiB peak RSS, 2 threads. The full remote summary is 49,300,939 bytes, SHA256 `a603ca8d79d9e34668aa56cc805b5db764fb0d885f4c08d58d613b6cbe3dd13e`; raw CSVs remain remote.
- The original checker `KeyError` and all attempts remain preserved. The review record is [phase-b-final-independent-review-20261031.txt](../reviews/phase-b-final-independent-review-20261031.txt).
- This report does not accept a lower mesh or truncation, make a speedup claim, or infer field-error causes.
