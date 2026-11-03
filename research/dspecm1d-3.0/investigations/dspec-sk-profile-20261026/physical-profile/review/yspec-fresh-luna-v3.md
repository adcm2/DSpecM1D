# Fresh independent review: YSpec physical adapter v3

Review date: 2026-10-26  
Scope: review of the v3 corrections against the v2 findings and the requested v3 checkpoint/snapshot provenance. Read-only; no compile or numerical execution. The mutable v4 work was not inspected.

## Decision

**Conditional code-review PASS for the parent’s v3 waveform sweep only.** The three v2 code/build blockers are corrected in the frozen v3 source/build: the YSpec library translation units and driver use the same isolated `TimeSeries.hpp`, the nested sampled-kernel profiler root is removed, and the sampled-kernel preflight now budgets substantially more than the output arrays. The CLI degree ceiling is 1536, and the waveform grid override still fixes epsilon while applying explicit N padding.

This is not a pass on physical convergence or on full-global sampled-kernel materialization. V3 is explicitly compile-only. The parent’s 512–1536 waveform results, their fixed-gate comparisons, and native/instrumented timing and output checks remain the evidence needed after execution. The full-global L>384 materialized-kernel case remains outside this v3 approval; review the separate v4 streaming implementation when it is ready.

## Re-review of v2 findings

- **FrequencyGrid ODR/header mismatch: corrected in the v3 build.** Both native and instrumented command manifests compile 14 library translation units from the isolated `physical-profile/yspec/source/src` tree plus the driver. Each compile command includes that directory, the local `TimeSeries.cpp` is compiled, and the link contains the 14 local objects with no `libYSpec.a`. Thus the physical driver, `InverseTransform`, `MakeGrid`, and all YSpec library TUs use the same copied `TimeSeries.hpp` definition. The v3 manifest reports `archive_linked: false` and 14 compiled library TUs.
- **Nested profiler Root: corrected.** The v3 `physical_yspec.cpp` owns the sole `mathprofile::Root`; the isolated `SampledKernels.cpp` has no `Root` construction. Its API scope can therefore nest under the driver's active root without invoking the root constructor’s terminate-on-nesting path.
- **Sampled-kernel memory preflight: materially corrected for materialized kernel runs.** Before allocating the full `places` vector, the code estimates 24 doubles per radial-angular point, 2048 bytes per radial-degree pair, 1 MiB per degree for angular workspace, and a 64 MiB reserve. It defaults to a 6144 MiB budget, reports the estimate, and keeps the cap kernel-only so waveform mode is unaffected. This covers the simultaneous accumulation and returned arrays with headroom and explicitly accounts for angular workspace. It remains a conservative model, not measured peak-RSS validation; the parent should use recorded RSS for any kernel run that fits the cap.

## Provenance and waveform readiness

- The v3 checkpoint manifest’s 141 listed files all match their recorded hashes, including both v3 executables, command manifests, copied source tree, and frozen v2 snapshot paths. The snapshot remains separated under `v2-checkpoint/`; the v2 snapshot’s original manifest points at live absolute paths, so the v3 manifest’s snapshot-file hashes are the useful freeze evidence here.
- Native and instrumented builds each contain 15 compile commands: the same 14 isolated library translation units and the physical driver. Their link object lists select the corresponding mode’s objects and do not reference the YSpec archive.
- The v3 grid constructor still receives `PhysicalIncludeDC` and `1024*padding`; the driver checks `N=1024*padding`, the corresponding `df`, and SI epsilon `5/3600 s^-1`. The new degree ceiling changes the allowed range without changing the waveform calculation path. These source checks support the requested fixed-epsilon 512–1536 sweep.
- The review found no v3 change to the previously reviewed v2 spectrum processing, vertical receiver selection, source/coordinate conversion, or adjoint weighting. This is a source review; those properties still need to be associated with the actual output provenance for each saved run.

## Evidence boundary

No v3 numerical case or matched native/instrumented pair is archived in the checkpoint. Therefore this conditional pass does not assert that degree 1536 converges, that any degree pair meets the protocol’s waveform and energy gates, that component timers close within the allowed overhead, or that a v3 profile is equivalent to its native run. The v2 corrected L=32 trace/spectrum and cross-code checks remain one-point consistency evidence only. Full-global L>384 kernels require a separate review of v4’s streamed implementation and its radial-domain/interface coverage.
