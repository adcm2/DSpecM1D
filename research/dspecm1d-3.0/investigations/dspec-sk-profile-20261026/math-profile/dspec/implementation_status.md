# DSpec mathematical-stage profiling clone

This directory is an isolated instrumentation copy of frozen v6. Its mathematical algorithm is unchanged. The original source and binary remain under `algorithm-audit/fast/v6/`.

## Change log

- 2026-10-26: Copied frozen v6, including solver helpers, and preserved its origin freeze manifest and build command.
- 2026-10-26: Added candidate-only nested profiling, count families and tensor dimensions, API/setup/reference/export timing metadata, a copied EquilibratedSparseLU timing shadow, and a candidate-stage map. See `../implementation_status.md` and `stage-to-math-map.json` for scope and boundaries.
- 2026-10-26: Compiled `fast_driver-profiled` and `fast_driver-profile-disabled` successfully through the shared two-slot compiler gate. Both use the same instrumented source; the second defines `MATHPROFILE_DISABLED`. No numerical run or output validation was performed here.
- 2026-10-26: Added a separate setup root that closes before candidate roots, with named stages for model/frequency input, SEM mesh/matrix setup, CSV and native-grid work, output streams/coordinates, and cache bookkeeping. Built `fast_driver-profiled-v2` and `fast_driver-profile-disabled-v2` through the shared two-slot compiler gate. Both use the same setup-instrumented source; the disabled build defines `MATHPROFILE_DISABLED`. No numerical run or validation was performed.

## Parent-owned next step

Parent owns fresh-process timed runs, native v6 comparisons, numerical output comparisons, overhead/accounting review, the independent Luna and Sol reviews, and final reporting. Use one OpenMP and one OpenBLAS thread. Keep `DSFAST_DISABLE_PROFILE` enabled so the older hot-path audit timer stays off.
- 2026-10-26: Added v3 tensor-stage scopes for semantic grid key build/lookup versus cold layout construction, representative/record capture, result allocation/tagging, angular-plan key lookup versus cold frame/Wigner setup, complex/real basis key lookup versus cold fill, radial l=0 scatter, forward-cache copying, and the receiver suffix route. Added the sum of actual BLAS GEMM m*n*k products as a counter. These remain coarse scopes outside per-point loops; v2 artifacts remain unchanged. Built `fast_driver-profiled-v3` and `fast_driver-profile-disabled-v3` through the shared two-slot compiler gate from identical source, with `MATHPROFILE_DISABLED` on the latter. No numerical run.
- 2026-10-26: Restored direct initialization for dynamic H assembly, candidate evaluation results, and tensor result vectors while retaining their profile scopes through return-value lambdas. Replaced v3 default-constructed/assigned basis keys with direct vector initialization; split key-construction and cache-lookup timing to preserve object-construction semantics. No numerical run.
- 2026-10-26: Built `fast_driver-profiled-v4` and `fast_driver-profile-disabled-v4` through the shared two-slot compiler gate from the same direct-initialization-preserving source; the disabled build defines `MATHPROFILE_DISABLED`. No numerical run. Source and all v1-v4 build hash entries verify.
