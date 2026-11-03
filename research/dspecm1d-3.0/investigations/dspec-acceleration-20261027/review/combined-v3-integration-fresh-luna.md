# Combined-v3 integration review: fresh Luna

**Result: HOLD pending one spool-metadata correction.** Static integration review only; no SSH, build, or science execution. Scope is `combined-v3` against the previously reviewed response/symbolic, packing-shared, nodal-shared, and global-v3 candidates. The global module review used as an input is hash-bound below.

## Findings

- **Pass — three-way engine header.** `combined-v3/engine_headers/Windowed3DFields.h` contains the response cache and basis-residual hooks alongside the nodal `PlannedValues` fixed-size recovery and per-element cached `ltg*` DOF layout. The packing merge retains per-trace label-layout deduplication, spin-specific key unions, and precomputed permutations. Dof caches remain keyed by one-sided `(element,radius)` sampling identity; fluid/inactive toroidal recovery is skipped before reading its unused slots. The response solver, factorization wrapper, and batch-cache headers are byte-identical to `responses-symbolic`.
- **Pass — global module integration.** The combined consumer and spool headers are byte-identical to `global-v3`. The combined `global_gsh_tile.hpp` differs from reviewed global-v3 only in its relative include path (`engine_headers/...` to `../engine_headers/...`). Its contents otherwise match the reviewed copy.
- **Pass — single baseline/capture control flow.** With `responseSpoolUsed`, the driver creates the spool and calls `baseline()` once with the record sink. The later fallback branch is guarded by `else if(!responseSpoolUsed)`, so it does not run another baseline in spool mode. The spool branch also bypasses `runKernel()`. With spool off, direct-cache enabled takes the cached baseline once; cache-off or preflight fallback takes the ordinary baseline once. The explicit runtime-capacity fallback intentionally discards the cache and reruns the baseline, with fallback metadata set.
- **Pass — DC, frequency, weights, and basis audits.** Kernel modes reject DC exclusion. `endBin` includes bin 0 and the first out-of-band guard, so the fixed case stores 53 frequency records when `frequency_bins` is 52. The capture loop visits `k=0..endBin-1`; unit receiver weights select vertical component 0, and post-energy spool audits store the physical adjoint mix `energyAudit->weights(0,k)` plus each bin's measured receiver-basis residual. Those audits are passed to the spool and written to its frequency-record manifest.
- **Pass — CLI controls and reported actuals.** Worker count is capped at 48; angle tile, output cap, working-memory cap, transform scratch reserve, spool/cache caps, and at-most-eight-side selection are validated and passed into the consumer. `result.json` reports configured caps and actual consumer frequency/side/output counts. The consumer side manifest records the actual native-side indices and identities; the driver marks whether selection was explicit. This covers the actual identities across the consumer manifest and summary result.
- **Hold — spool basis-residual summary is reported as zero.** In spool mode `responseCache` is null, so the `result.json` field `response_basis_reconstruction_max_relative` is hardwired to `0.0`, although the driver collected nonzero-capable per-bin residuals in `responseBasisResiduals` and writes them to the spool manifest. The per-frequency evidence is present, but the top-level result field can read as a measured zero. Set this summary to the maximum captured residual when spooling (or rename/scope the current field explicitly as cache-only and add a spool-residual summary) before treating result metadata as an honest summary.

## Hash bindings

SHA-256 values observed during this review:

- `combined-v3/engine_headers/Windowed3DFields.h`: `74d7c1bc07e27a247b59f875b3e4ce8ccb5630c0533b800e8c16460cd7cb48e0`
- `combined-v3/physical_driver.cpp`: `afaab77fda831ef5140e8a1fc929e71487042f63170ca1b66a81f183bb02ccf4`
- `combined-v3/engine_headers/SensitivityKernels.h`: `3c4dc037841e2923ad2c3064ddda4c818266379e47bf74ca79efacca5b1307a8` (same as `responses-symbolic`)
- `combined-v3/profile_equilibrated_sparse_lu.hpp`: `8dd2b4c7be5eb42d35d49d462eb0d0b04a36e5264660df4126abcdac7ee03047` (same as `responses-symbolic`)
- `combined-v3/solvers/batch_response_cache.hpp`: `e4ae716d3794b2ab32d5455802178bc1c837a0a53e46e2d47f65e64b5a76ffc0` (same as `responses-symbolic`)
- `combined-v3/global/global_kernel_consumer.hpp`: `ab741d8fda380b31c473a75034d100fd314c422f66031703df4b29c1301ec994` (same as `global-v3`)
- `combined-v3/global/pole_response_spool.hpp`: `39a8e1bd6de5d52257fdfcb365f68fb09c9609439696dbba3e9c81b7d3605ea4` (same as `global-v3`)
- `combined-v3/global/global_gsh_tile.hpp`: `513500d769c2b70da05a0d7a82911c5a63690e0643b4b24b26a68d58a0308445`; reviewed `global-v3` source: `8cd1503bc22181d4e711f1b1de0841d3a471ac9cfcf77570284d3ff63b7f8fdf` (only the relative include differs)
- `packing-shared/engine_headers/Windowed3DFields.h`: `e720a9ee3f926e7e5e86c8985ded28692622cf97b46d4752173bee33579550ac`
- `nodal-shared/engine_headers/Windowed3DFields.h`: `87ac1416c66adfdef0fcf641343410847932ccf5a3eca7891342bdf973ad6b45`
- `review/global-spool-v3-fresh-luna.md`: `939793120f1b310a61d142181c3e332abdeaff850cf4e46fa14ee8694a03b92f`

## Gate boundary

This review does not establish compilation, runtime numerical correctness, resource peak, or speedup. The reviewed global-v3 module's prior PASS remains separate; the new HOLD concerns the combined driver's top-level residual summary.
