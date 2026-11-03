# Fresh independent review: YSpec bounded streaming adapter v5

Review date: 2026-10-26  
Scope: frozen v5 source, build provenance, disk/RAM safeguards, v3 output-layout parity, and profile interpretation. Read-only review; no compilation or numerical execution by this reviewer.

## Decision

**Conditional code-review PASS for proceeding with parent-owned bounded high-degree profiling.** I found no source-level blocker in the v5 streaming accumulation or output conversion. The parent-provided low-degree gate records exact v3/v5 values after layout transpose, exact native/instrumented v5 bytes, and finite output. This establishes the tested L=3 adapter-path equivalence, not high-degree physical convergence or successful full-domain resource behavior.

The v5 checkpoint manifest has 90 hash entries; all 90 files verified against the recorded hashes. Both build command records compile 14 isolated local YSpec library translation units plus the driver and link the corresponding local objects, with no YSpec archive. The frozen v4 snapshot remains separately identified by the v5 checkpoint.

## Correctness and output layout

- The disk branch is limited to exactly one measurement, one sampling, absolute Love parameters. It creates a zero-filled scratch file ordered `[angle][node][8 elastic channels]`, then performs each block read/modify/write in the same frequency and angular-block order as the v3 in-memory sum. The local contraction is transposed into the scratch order without changing its per-angle arithmetic (`v5/source/src/SampledKernels.cpp:975-1018`).
- Final conversion reads at most 64 angles, fills only the five Love moduli, applies the same `modulus * volumeUnit` normalization and `Transform(..., ParameterSet::Love)` as v3, checks every result for finiteness, and writes native-endian float64 `[angle][node][A,C,F,L,N]` (`:1033-1094`). This intentionally differs from v3's in-memory parameter/node/place storage; consumers must apply the declared transpose when comparing or reading old-format files.
- The direct angular coordinate constructor uses `2L+1` Gauss-Legendre colatitudes and `FastFFTSize(4L+1)` longitudes in theta-major, ascending-longitude order, and avoids a global Wigner table. The supplied exact L=3 output comparison provides a low-degree end-to-end parity point. High-degree coordinate and full-domain behavior remain unverified by this review.
- The new disk path is separate from the retained v3 in-memory path. The v5 diff adds the accumulator/output helpers and streamed branch; the monolithic contraction and transform path remains intact for comparison.

## RAM, disk, and failure behavior

- Before reserving the global place vector, the driver estimates 80 bytes per angle for the 32-byte `KernelPlace` plus 40-byte `PointGeometry`, strain and bounded angular workspaces, both the contraction and scratch-update tiles, the simultaneous 8-channel/5-channel transform buffers, node materials, and a 512 MiB fixed reserve. It rejects cases above the configured memory budget (`v5/physical_yspec.cpp:388-430`). This is a conservative source estimate, not a measured high-degree RSS guarantee; retain peak-RSS evidence for the parent runs.
- The disk estimate is 13 doubles per radial-angle point (8 scratch plus 5 output), checked against the configured cap and currently available space plus a 1 GiB reserve before angular allocation. Scratch and output sizes have checked integer bounds; positioned I/O handles EINTR/short transfers and checks every requested range (`v5/source/src/SampledKernels.cpp:390-535`).
- Scratch/output use exclusive creation, fsync, and preserve scratch and any partial `kernels.bin` after an API/conversion failure. Scratch is unlinked only after the output fsync succeeds (`:429-445`, `:475-529`, `:1050-1094`). The preflight cannot guarantee later writes against filesystem quotas, concurrent disk use, or allocation failure; a failed run must be retained as failed evidence, and downstream readers should require successful `result.json` plus the expected byte count rather than treating a present `kernels.bin` as complete. Retries should use a fresh output directory or explicitly preserve and move prior failure artifacts.
- No induced full-disk/quota failure or full-domain high-L RSS run is in the v5 checkpoint. Those remain parent validation items; this is not a disk-capacity pass.

## Build, dependencies, and profile interpretation

- Native and instrumented v5 link commands have the same external dependency list and order as their corresponding v2 adapter link commands. Neither list contains BLAS, LAPACK, OpenBLAS, or MKL. Both binaries are directly linked from their mode-matched v5 objects; no `libYSpec.a` is linked. The command manifests and checkpoint hashes provide the source/object/binary provenance.
- The parent-provided `yspec-stream-gates-v5.json` records `v3_vs_v5_exact_values: true`, zero relative L2 error, finite output, and byte-identical native/profile output. The instrumented L=3 run records 16.9 MB peak RSS, which is a small-case measurement only.
- The profile contains two nested scopes both named `api.ComputeSampledKernels`: the driver wrapper (`v5/physical_yspec.cpp:597-604`) and the implementation (`v5/source/src/SampledKernels.cpp:539-545`). The saved profile therefore sums both inclusive intervals: 5.6346 s across two calls, exceeding the 2.8191 s enclosing `physical_case_total` and the 2.9105 s process wall time. This is a duplicate-name/inclusive-accounting caveat, not evidence of incorrect kernel execution. For timing closure use exclusive entries and the driver's single wall timer; do not report the summed API inclusive value as elapsed API time. The new scratch/output component scopes are separately present.

## Evidence boundary

The reviewed sources and L=3 comparison support a bounded-streaming code-path pass. They do not establish angular-degree convergence, high-L resource feasibility, disk-exhaustion behavior, full-domain numerical equivalence at production degree, or physical validation. Those claims remain tied to the parent's saved runs and explicit convergence/resource gates.
