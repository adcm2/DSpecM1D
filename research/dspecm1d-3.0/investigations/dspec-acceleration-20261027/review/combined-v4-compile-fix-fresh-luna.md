# Combined-v4 compile-fix review: fresh Luna

**Result: PASS for the bounded source-static compile-fix gate.** The three targeted header edits preserve the prior combined-v3 source-integration PASS. No build, SSH, or science run was performed.

The spool writer accessor is now non-const, matching `std::fstream::tellp()` on the open-stream path; the closed-stream file-size path and caller behavior are unchanged. `BatchResponseFootprint` and `estimateBatchResponseFootprint()` are now public for the intended preflight call, with `private:` restored immediately before `solveFrequencyModes()`. The consumer replaces calls to engine-private angular-memory helpers with local checked helpers.

The new consumer helpers retain the engine formulas: for each `n=-nMax..nMax`, `nMax=min(2,maxDegree)`, they sum `GSHIndices<All>(maxDegree,min(mMax,maxDegree),n).Size()` and charge each value as `sizeof(double)`, plus `(2*nMax+1)*sizeof(size_t)` index bytes. Coefficient storage remains `(2*nMax+1)*5*(maxDegree+1)*sizeof(double)`. The consumer adds checked size arithmetic; the existing `lmax<=1536` guard also keeps the remaining small integer expressions in range. Include resolution is unchanged through `global_kernel_consumer.hpp` → `pole_response_spool.hpp` → `global_gsh_tile.hpp`, which includes the engine GSH declarations. The copied module/header surface was not otherwise changed in this diff.

The direct-cache footprint estimate remains internally consistent: with `values=n*(forwardCols+basisCols)`, the two complex RHS matrices occupy `16*values` bytes, the two complex response matrices occupy another `16*values`, and the two bit signatures occupy `(2*values+4)*8 = 16*values+32` bytes. Total: `48*values+32`, exactly the estimator. The reviewed example `n=10`, forward width `5`, basis width `3` gives `3,872` bytes.

## Hash bindings

- `combined-v4/engine_headers/Windowed3DFields.h`: `a66f74690d5655335598948d2b804c83d19eca6606481a26214389e0dff9feb0`
- `combined-v4/global/global_kernel_consumer.hpp`: `7bdb242b22635e53684f0c2e5e8fc4319e0c235b39023cc655409449373d24dd`
- `combined-v4/global/pole_response_spool.hpp`: `5e8868439771cb502cb4e3005039f6986e9d0d6ad23a58b168a7a61c2c9152dc`
- `combined-v4/physical_driver.cpp` (unchanged from v3): `357e0038c5e56188d24b47de99697c79fe4a1b3bd7e3b2b340cfb274b928d5cf`
- Prior combined-v3 integration-v2 PASS: `f316d2d458c916823edf71fa5a694d0cc2083fd5968d2be6297e891622e780f5`

## Gate boundary

This PASS covers only the targeted source edits and their access, formula, and include relationships. It does not establish compilation or runtime correctness; it does not authorize a science run. The prior source-integration PASS remains in force for the reviewed unchanged integration paths.
