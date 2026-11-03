# Combined-v3 integration review: fresh Luna v2

**Result: PASS for the reviewed source-static integration gate.** No SSH, build, or science run. This v2 review supersedes the residual-metadata HOLD in `combined-v3-integration-fresh-luna.md`; that original review is preserved. Other integration findings are carried forward from that hash-bound review.

## Residual aggregation correction

The current driver computes `response_basis_reconstruction_max_relative` from the maximum captured per-frequency residual when `responseSpoolUsed` is true, and from the response cache's measured weighted-basis reconstruction maximum otherwise. In spool mode, kernel/DC validation guarantees bin 0 is included, `endBin` is positive, and `baseline()` fills each entry in the `responseBasisResiduals` vector before the maximum is taken. The same per-bin residuals continue to be written into the spool audit records. This resolves the prior top-level zero-summary mismatch without changing the cache-mode measurement.

The previous integration checks remain valid: the three-way field header retains response hooks, fixed-size nodal recovery, cached one-sided LTG layout, and packed coefficient layouts; global-v3 consumer/spool modules remain byte-identical, with only the relative include path changed in the GSH header copy; spool mode uses one baseline/capture pass and skips the direct kernel path; DC-inclusive frequency, vertical unit receiver weights, q-mix, residual audit, resource/side controls, and consumer metadata flow are unchanged.

## Hash bindings

- `combined-v3/physical_driver.cpp`: `357e0038c5e56188d24b47de99697c79fe4a1b3bd7e3b2b340cfb274b928d5cf`
- `combined-v3/engine_headers/Windowed3DFields.h`: `74d7c1bc07e27a247b59f875b3e4ce8ccb5630c0533b800e8c16460cd7cb48e0`
- `combined-v3/engine_headers/SensitivityKernels.h`: `3c4dc037841e2923ad2c3064ddda4c818266379e47bf74ca79efacca5b1307a8`
- `combined-v3/profile_equilibrated_sparse_lu.hpp`: `8dd2b4c7be5eb42d35d49d462eb0d0b04a36e5264660df4126abcdac7ee03047`
- `combined-v3/solvers/batch_response_cache.hpp`: `e4ae716d3794b2ab32d5455802178bc1c837a0a53e46e2d47f65e64b5a76ffc0`
- `combined-v3/global/global_gsh_tile.hpp`: `513500d769c2b70da05a0d7a82911c5a63690e0643b4b24b26a68d58a0308445` (reviewed global-v3 source differs only by relative include path)
- `combined-v3/global/global_kernel_consumer.hpp`: `ab741d8fda380b31c473a75034d100fd314c422f66031703df4b29c1301ec994`
- `combined-v3/global/pole_response_spool.hpp`: `39a8e1bd6de5d52257fdfcb365f68fb09c9609439696dbba3e9c81b7d3605ea4`
- Prior global-v3 Luna review: `939793120f1b310a61d142181c3e332abdeaff850cf4e46fa14ee8694a03b92f`
- Preserved original combined-v3 integration HOLD review: `fce4dfb9ce7721b73d7e8a5532eb52ce003ec933c17fae0c09440bd2e30828e7`

## Gate boundary

This PASS supports the next independent build-only review. It does not establish compilation, runtime numerical correctness, resource peak, or scientific completion, and it does not authorize a science run.
