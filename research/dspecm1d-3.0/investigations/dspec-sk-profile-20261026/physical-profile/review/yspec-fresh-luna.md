# Fresh independent review: YSpec physical adapter

Review date: 2026-10-26  
Scope: read-only source, protocol, build-manifest, and saved-run review. No compilation, numerical execution, or implementation edits were made.

## Decision

**Hold before matched profiling and sampled-kernel profiling.** The waveform adapter's corrected v2 trace and spectrum are internally consistent with one another and its cross-code L=32 comparison passes the stated 1e-3 gates. The build still has a `FrequencyGrid` definition mismatch across translation units, and the instrumented kernel path has nested profiler roots that terminate. The sampled-kernel memory preflight also understates peak allocations. These are concrete blockers to the requested native/instrumented and kernel profiling evidence.

The saved single-degree comparison is not a physical convergence result. No YSpec degree, FFT-padding, radial-tolerance, or full-grid kernel convergence sequence is present here, so no physical convergence or global-kernel conclusion is supported.

## Confirmed correct in the inspected adapter

- The v2 waveform uses raw `SeismogramSpectra` SI m·s values, applies the frequency taper through `MakeTaper` in `InverseTransform`, then the library's inverse transform undamps by `exp(epsilon*t)` and applies the frequency spacing. The time-domain energy then uses the 1000–1200 s Hann factor and `0.5*dt`. The receiver is `Vertical` and the kernel measurement direction is `.up=1`.
- The YSpec source and receiver latitude inputs are converted from the protocol's geocentric latitudes to geographic latitudes with `atan(tan(phi)/0.993277)`. The tensor entries, 50 km source depth, and 3 km receiver depth are explicit and match the protocol's intended SI setup. Gravity is full, attenuation is off, and the requested quantity is displacement.
- The isolated grid override actually supplies `includeDC` and the explicit `1024*padding` sample count to `FrequencyGrid`. Its decay and duration both scale with padding, leaving SI epsilon fixed at `5/3600 s^-1`; the adapter checks N, df, and epsilon at runtime. The corrected archive used `--include-dc --fft-padding 1`, reporting N=1024, df=1/5120 Hz, and epsilon=5/3600 s^-1.
- The copied `Kernels.cpp` uses `df` at DC and Nyquist and `2df` for interior positive frequencies. The frozen adjoint input is `Hann*u_z` without an extra dt, and the returned linear measurement is recorded as `2E`.
- Sampled kernels use the complete `ForBand(lmax,2,2.0)` angular grid. The radial input retains row order and repeated radii; at model interfaces, the implemented endpoint mapping is q=0 to `Side::Above` and q=NN-1 to `Side::Below`. The YSpec sampled-kernel source performs radial solves outside the per-64-angle geometry blocks, so the ODE solve is not duplicated per angular batch.
- Both `checkpoint-manifest-v2.json` (45 entries) and `math-profile/yspec/source-build-hashes.json` (35 entries) match the current files. The recorded corrected waveform run executable hash matches the v2 native binary. Native and instrumented command manifests link the corresponding native/instrumented math objects; the v2 source/object/binary hashes are present.

## Blockers

1. **`FrequencyGrid` has differing class definitions across the executable (formal ODR violation).** The physical adapter's copied `source/src/TimeSeries.hpp` adds `includeDC` and `sampleOverride` constructor arguments and changes constructor logic. The linked YSpec library and the linked math-profile hook objects use the unmodified production/math-profile `TimeSeries.hpp`. The data members currently remain the same, and the saved v2 waveform passes its independent spectrum-to-trace check, but the class definition is not token-identical across translation units. This fails the requested no-ODR/mismatched-header criterion. Compile every translation unit that sees `FrequencyGrid` against one identical isolated header/source set, including the transform implementation and hook objects, or redesign the override so the class definition remains identical.

2. **Instrumented sampled-kernel execution will terminate on nested profiler roots.** The physical driver opens `mathprofile::Root` around both modes. The instrumented `ComputeSampledKernels` implementation also constructs a `mathprofile::Root`; the shared root constructor calls `std::terminate()` if a root is already active. Therefore the instrumented kernel mode cannot produce the requested component profile in its current form. Keep one root owner for that call path, then compare native/profile outputs and timer closure/overhead before accepting profile results.

3. **Sampled-kernel memory preflight understates peak memory.** `readSamples` gates on eight doubles per radial-angular point. `ComputeSampledKernels` keeps a five-double-per-point `total` array while constructing a separate five-double-per-point returned `SampledKernels` array, so the principal arrays alone peak at about ten doubles per point. The preflight also omits the full angular `places`/geometry storage and other solver data. The generalized-Legendre tables themselves are built in 64-angle blocks, so they are not a full-grid table, but the current estimate still cannot serve as the advertised pre-allocation memory cap. Raise and document a conservative peak estimate or stream/chunk the output before scheduling large complete grids, especially L>=192 under the 8 GiB limit.

## Saved evidence and limits

- `runs/yspec-corrected-l32-t7-q1-v2/spectrum-integrity.json` reports a 2.23e-16 trace/spectrum relative L2 discrepancy and matching independently recomputed energy. Its `cross-code.json` compares one DSpec L=32 trace and energy and passes both 1e-3 thresholds. This is a corrected single-point comparison, not a convergence sequence.
- The earlier `runs/yspec-initial-l32-t7-q1-v2/spectrum-integrity.json` is rejected evidence: it fails with relative L2 0.954 and an energy mismatch. Preserve it as rejected; do not pool it with v2 results. I did not find a directory literally named `v1rejectedresults` under `physical-profile`.
- The corrected run is native-only (`run.json` has `profiles: []`; `result.json` has empty `profile_stages`). There is no archived instrumented physical-profile run here, so native/instrumented timing equivalence, component-timer closure, and profiling overhead bounds remain unverified. The current native waveform run is not evidence that the instrumented kernel mode is executable.
- No YSpec waveform refinement series or sampled-kernel result is archived under this physical-profile run set. In particular, code support for N overrides and a complete angular grid does not demonstrate convergence or full-global radial coverage.
