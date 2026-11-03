# Fresh Luna High review: frozen DSpec physical adapter v2

**Scope:** Independent source/freeze review only. No compilation or numerical run was performed. The review is limited to the frozen v2 source, build metadata, adapter protocol, and DC/measurement math note. It does not establish PREM waveform or kernel convergence.

## Decision

The code implements most requested adapter mechanics, but the kernel campaign is **not ready for a full physical pass** because the required baseline identity `q·P(B raw)=2E` is neither checked nor recorded, and the kernel loop does not use the accepted complex tensor baseline path. The profiling map also needs to be reconciled with the physical adapter's stages before timing claims are made. Physical correctness and profiled/native parity remain runtime questions for the parent-owned runs.

## Findings

1. **Blocking audit gap — no `q·P(B raw)=2E` check.** `physicalEnergyAdjoint()` constructs `B raw`, computes `P(B raw)`, forms `q=df_SI eta P(B raw)`, applies `P*`, and applies `B` on the adjoint path. This is the required algebraic ordering. However, it returns only `qhat`: it does not compute or emit the real inner product, compare it to `2E`, or fail outside a tolerance. The driver separately computes `E`, so the protocol's mandatory identity cannot be audited from run artifacts. See `physical_driver.cpp:197-210` and energy evaluation at `:329-333`.

2. **Profiling provenance gap — map and roots do not match the physical adapter.** The source creates one `mathprofile::Root("dspec", -1)` for the entire executable and a nested `dspec_setup` scope; it does not create the separate setup and per-candidate roots described in `stage-to-math-map.json`. The map contains no physical-driver stage entries such as `baseline_forward_receiver_only`, `measurement_filter_and_inverse_fft`, `physical_energy_adjoint_all_output_bins`, `frequency_workspace_modal_solve_and_adjoint_construction`, or the streamed physical kernel I/O/contraction stages. Exclusive timing is emitted by the generic profiler, but these physical costs are not fully mapped to the report taxonomy. See `physical_driver.cpp:316-325` and `source-v2/stage-to-math-map.json:1-67`.

3. **DC/endpoint audit is partial.** The baseline explicitly solves bin 0 at `z=w[0]-i epsilon`, multiplies by `-i/z`, and records the absolute imaginary part of component 0. The spectrum and inverse-FFT paths project DC and Nyquist to real values. There is no explicit finite-value check on the raw endpoint or normalized DC imaginary-residual gate, and the raw Nyquist endpoint/residual is not saved. `dc_imaginary_residual_raw` alone is an absolute number whose significance depends on signal scale. See `physical_driver.cpp:135-159`, `:162-169`, and `:340-353`.

4. **External baseline match is not established by source inspection.** The adapter fixes the same stated source/receiver coordinates, step moment tensor entries, finite-damping convention, and vertical component in its parameter text, and exports the raw SI spectrum and filtered SI trace. The frozen DSpec code does not itself compare against the unchanged YSpec baseline, so matching must be established from the parent-owned paired evidence; source shape alone is insufficient.

5. **Blocking performance/methodology mismatch — kernel loop uses the generic packed sampler.** The angle-tile loop calls `workspace.sampleForKernelsPacked(...)`, which routes through `samplePoleFrames`; it does not call the accepted v6 complex tensor sampler. Thus a profiled physical kernel run would measure a different contraction route from the optimized prototype baseline, and its performance cannot be reported as the accepted tensor-kernel method. Preserve this distinction in the v2 evidence; a future implementation phase should route the campaign through the agreed tensor path and verify numerical parity before using its timings. See `physical_driver.cpp:272-273` and `Windowed3DFields.h:1410-1419` versus `:1430-1444`.

6. **Full-grid storage has overflow protection but no disk-capacity preflight.** Kernel mode computes the full `pointCount` and `resultBytes` and checks `size_t` overflow, then creates/zeros a file of that size and performs per-frequency read/modify/write accumulation. It does not check available filesystem space or estimate the combined kernel and optional contact output before starting. A high-L all-radial-sides run can therefore fail late after expensive solves. See `physical_driver.cpp:230-244` and `:306-309`.

7. **Pole-only still pays avoidable setup work.** In spheroidal and toroidal loops the code still allocates an `rhs` with `2l+1` columns, constructs source coefficients and receiver maps, then leaves the weighted RHS columns zero when `poleOnly` is true. It avoids populated nonradial weighted RHS solves and filters their blocks from capture, but the remaining preparation costs should be included in performance interpretation. See `Windowed3DFields.h:723-755` and `:767-799`.

## Criteria supported by source/build inspection

- The vertical observable uses component 0 and reports it as outward/up radial displacement. Raw spectrum export applies `lengthScaleMetres * timeNorm`; the inverse-FFT trace applies the length scale after `freq2time`. Energy is `0.5 * dt * sum(Hann * u_z^2)`, with square-root Hann supplied to the linear map.
- The band is applied to the damped-contour spectrum before inverse FFT. The energy adjoint forms all output bins before `P*`, then applies the band again. Its endpoint projection uses the real FFT subspace and the one-sided factor is `df` at DC/Nyquist and `2df` at interior bins.
- `qex` is passed into the frequency-grid constructor while `dt` and requested duration remain explicit. The resulting `epsilon` and FFT size are written to run metadata, allowing the fixed-epsilon padding schedule to be audited.
- Native radial element/node sides are enumerated, retaining one-sided interface values and excluding the zero-radius center. `--radial-traces 0` selects all available sides; the default selects one side, so a full radial-side kernel run must set the flag explicitly. The angular grid has `2L+1` Gauss-Legendre colatitudes and `FastFFTSize(4L+1)` longitudes, and points are streamed in tiles of 64 angles. Output metadata explicitly says the selected-side product is not a full Earth volume integral.
- Each frequency workspace is prepared once outside its angle-tile loop; kernels and optional source-shell contacts accumulate over frequency. The contact file is separate. `--contacts on` with pole-only enabled is rejected, forcing the global modal route for contacts.
- Pole-only mode skips the per-order weighted nonradial adjoints while retaining the unit receiver-basis adjoints and radial l=0 branch used by the pole sampler. This route has not been numerically compared to the full route in this review.
- Kernel sampling passes backend value `2`, which is `DotBackend::Eigen`; the frozen build command has neither `DSFAST_HAVE_CBLAS` nor `EIGEN_USE_BLAS`. LAPACK/BLAS linkage is recorded for existing solver symbols only.
- The profiled and profile-disabled frozen executables are built from the same frozen source tree with profiling macros as the recorded difference. Both executable hashes and every listed source-v2 file hash passed `sha256sum -c`.
- Cache budgets are explicit (32 MiB) and sampling uses bounded radial/angle tiles. This supports bounded cache/tile design by inspection; actual peak memory and runtime are not established here.

## Required parent-owned evidence

Record the `q·P(B raw)/(2E)` identity and a scale-aware DC/endpoint residual check; compare the external YSpec baseline without renormalization; verify full-route versus pole-only equivalence at a controlled case; run profiled/native parity on the accepted tensor path and capture complete solve-cost timing with the physical stage names mapped; preflight the combined output size against available storage; and keep the configured radial subset, degree, padding, convergence failures, and resource ceiling explicit. Do not infer full physical convergence from this code review.

## Freeze verification

`freeze-v2.sha256` verified both binaries (`physical_driver-profiled-v2-frozen`, `physical_driver-profile-disabled-v2-frozen`). `source-v2.sha256` verified the frozen driver, headers, solver helpers, and stage map. The review follows the provenance note that v1's complete source tree was not retained; no v1 reproducibility or scientific result is claimed.
