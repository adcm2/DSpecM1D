# Combined-final-v1 integration review: fresh Luna

**Result: PASS for the bounded pre-build source-integration gate.** The final candidate retains the reviewed cache-v2 integration and replaces only `engine_headers/Windowed3DFields.h` with the exact combined-v5 header. The only additional difference from cache-v2 is the implementation-status note recording that assembly. No build, SSH, or science run was performed.

## Review findings

- Recursive source comparison confirms `combined-final-v1` matches `combined-cache-v2` byte-for-byte except `engine_headers/Windowed3DFields.h`; the status file adds the final-candidate provenance note. The physical driver, batch cache, global consumer/spool modules, copied GSH header, and other engine headers are unchanged from reviewed cache-v2.
- The final field header is byte-identical to combined-v5. Comparing v4→v5 shows the reverse nodal-shared delta: it removes the cached SEM-DOF arrays and `PlannedValues`, restores the prior `Eigen::MatrixXcd` recovery path and per-node `sem.ltg*` lookups, removes precomputed nodal bridge plans, and restores on-demand plan construction in the optional bridge. These are precisely the changed recovery/packing sections in the nodal-shared versus packing-shared comparison; surrounding response/cache/global integration is retained. This is a source comparison, not a claim of runtime bitwise output.
- The response integration remains present in the transplanted header: `BatchResponseFootprint` and `estimateBatchResponseFootprint()` remain available, `solveFrequencyModes()` resumes private scope, and the pole block/source-frame/cache hooks remain in the solve and reconstruction paths.
- The previously reviewed cache-v2 source PASS applies to the byte-identical remaining source set. Its gate explicitly excluded compilation, runtime numerical agreement, and scientific claims; those limits carry forward.

## Hash bindings

- `combined-final-v1/engine_headers/Windowed3DFields.h` (identical to v5): `fcf8c2ec7d428bd3593b0030e75cfcf524382e2d4d5f6496f9aa01824149cb35`
- `combined-final-v1/physical_driver.cpp`: `8b763d85bf1552809bdf44b6ed1a5e16af6c218c1f52ad07f30e6fe5c237ae74`
- `combined-final-v1/solvers/batch_response_cache.hpp`: `e4ae716d3794b2ab32d5455802178bc1c837a0a53e46e2d47f65e64b5a76ffc0`
- `combined-final-v1/global/global_kernel_consumer.hpp`: `3a84e1cfdf02e48ff48f8ffd08146e122894335117126ff3081781e8a1845b90`
- `combined-final-v1/global/pole_response_spool.hpp`: `5e8868439771cb502cb4e3005039f6986e9d0d6ad23a58b168a7a61c2c9152dc`
- `combined-final-v1/global/global_gsh_tile.hpp`: `c19f438af8ba56dadb84617a4f03302e52075d36fa6939b88146f7c392fe9eac`
- `combined-final-v1/engine_headers/Windowed3DKernels.h`: `a852dead24ba03c1c4139ef926b807a88551ce2c2d77f443ddf149cf150d5506`
- `combined-final-v1/implementation_status.md`: `f89a782f11744043a198c3a0a63ef2c525f1e75d30c5d3c36a76eb85d06a5480`
- Prior cache-v2 fresh Luna review: `global-cache-v2-fresh-luna.md`, SHA-256 `96487bb573d0a465caf0f9e2e3a7be9a22e12ea9bec05874a5ddb934c2de116c`

## Gate boundary

This PASS is source-static and permits the next authorized build-only gate. It does not establish compilation, runtime behavior, bitwise output, scientific correctness, or speedup, and does not authorize a science run.
