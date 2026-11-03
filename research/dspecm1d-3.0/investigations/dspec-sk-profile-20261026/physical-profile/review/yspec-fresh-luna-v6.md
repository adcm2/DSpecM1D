# Fresh independent review: YSpec kernel probe v6

Review date: 2026-10-26  
Scope: v6 checkpoint/source only, focused on probe-vs-full sampling and metadata, unchanged full-kernel behavior relative to v5, allocation bounds, and profiler scope names. No compilation or high-cost execution was performed by this reviewer.

## Decision

**Conditional code-review PASS for the explicit positive 64-angle probe command; one probe-boundary input validation issue remains.** With `--probe-angle-count 64`, the code constructs only the first 64 place coordinates, retains every supplied radial node, and reports the full grid separately. The v6 checkpoint contains 91 hash entries and all 91 recorded hashes and sizes verified. The parent-owned L=3 native/profile probe outputs are byte-identical and each matches the first 64-angle prefix of the full v5 native and instrumented output files.

The new CLI parser uses `std::stoull` and checks only complete token consumption. That conversion accepts a leading minus: `--probe-angle-count -1` becomes a very large unsigned value, which `min(requested, targetAngleCount)` turns into the complete target angle count. It can therefore make `kernel-probe` allocate and process the full geometry. The supplied explicit `64` command is unaffected, but reject a leading minus (or parse with an unsigned `from_chars` path that rejects signs) before describing probe mode as bounded for arbitrary CLI input.

## Probe subset and metadata

- The full target grid dimensions are still computed as `(2L+1) x FastFFTSize(4L+1)`, but `angleCount` is the requested count clamped to the target. Preflight, `places.reserve`, and the theta-major/longitude-major generation loop all use this actual count (`v6/physical_yspec.cpp:404-459`). In the valid 64-angle path there is no full `KernelPlace` or `PointGeometry` allocation; only the `O(L)` colatitude vector and requested places are held. All radial rows remain in the sample input.
- Metadata distinguishes the probe subset from full target: for L=3 and 179 radial rows, the saved result reports 64 actual versus 91 target angles, 11,456 actual versus 16,289 full-target points, 458,240 actual output bytes versus 651,560 full-target output bytes, and `kernel_probe_subset=true`. The actual file is exactly 458,240 bytes. The estimated scratch-plus-output is 1,191,424 bytes, matching 13 doubles per actual radial-angle point. The separate name `kernel_probe.bin` and sampling-scope text explicitly identify a diagnostic subset (`:491-559`, `:618-647`).
- I compared both parent-saved L=3 probe files with the corresponding full v5 output. Each probe file equals the first 458,240 bytes of the 651,560-byte full file; native and instrumented probe files are also byte-identical. This is direct layout-prefix evidence for the tested L=3 fixture, not proof of high-L convergence or full-kernel correctness.

## Full-kernel path relative to v5

For `kernel` mode, `maxAngles` remains zero, so `angleCount` equals the full target count. The full place order, full preflight sizing, output destination `kernels.bin`, frequency accumulation, transform, and binary storage path are unchanged from v5. The v6 driver adds probe metadata fields and labels the unchanged binary representation as v6; those are metadata/version-label changes, not a data-layout change. The only sampled-kernel implementation change from v5 is the profiler scope string (`api.ComputeSampledKernelsImpl`); no numerical expressions or accumulation order changed.

## Profiling and evidence boundary

- The two scopes are now distinct: the outer driver timer is `api.ComputeSampledKernels`, and the implementation timer is `api.ComputeSampledKernelsImpl`. The saved L=3 instrumented profile records one call to each, removing v5's same-name aggregation. They remain nested inclusive timers, so their inclusive durations must not be summed; use exclusive values or the single driver wall timer for closure.
- The v6 checkpoint is compile-only by declaration; the parent-owned L=3 outputs provide the limited probe checks above. No conclusion about L=1536 cost, physical convergence, full-grid memory feasibility, or high-degree profile closure follows from this source review. Those depend on the separately scheduled parent probe and its matched run evidence.
