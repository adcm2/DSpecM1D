# Fresh independent review: native grid and physical measurement setup

**Status: PASS for the bounded grid, Eigen-dispatch/equivalence, and synthetic measurement checks only. No physical PREM kernel or convergence PASS.**

This is a read-only source and evidence review of the grid generator, isolated Eigen candidate, and synthetic measurement probe. I did not compile or run numerical code. Runtime artifacts were produced by the parent. The PREM/windowed-energy fixture is a proposed next benchmark setup, not a solved Earth model or physical kernel result.

## Findings

1. **The native grid call matches the requested product-resolving grid.** The generator instantiates `GaussLegendreGrid<double, All, All>` and calls `ForBand(L, 2, 2.0, Estimate)`. In the inspected GSHTrans header, `ForBand` resolves `ceil(2L)`; the grid then uses `lGrid + 1` Gauss–Legendre colatitudes and a fast longitude count at least `2*lGrid + 1`. It samples the complete degree-`2L` product band and writes the Cartesian product with GSHTrans weights. `resolved_band` comes from `grid.MaxDegree()`.

2. **Trace metadata and one-sided radial identities are preserved in source and generated fixtures.** Target `J_Z,A` rows are selected before remaining fields are parsed; trailing CRs are removed from header and rows. Trace identity is `(element, exact radius bits)`, preserving same-radius interface samples on opposite elements. Split trace blocks and duplicate angles are rejected. Metadata, `dr_m`, ordered angular blocks, and input weights are checked across traces. Output keeps the radial metadata and replaces only angles and `dOmega`; the kernel value stays blank and `deltaX_Pa` stays zero.

3. **The dispatch correction and bounded Eigen output gates pass.** The original candidate's implicit enum values (`Scalar=0`, `Eigen=1`, `OpenBlas=2`) disagreed with the driver's IDs; its failed no-CBLAS run is preserved as rejected-dispatch evidence. The isolated candidate now declares `Scalar=1`, `Eigen=2`, and `OpenBlas=3`. With ID 2, `synthesizeComplex` reaches the explicit Eigen `basis * coeff` branch before any CBLAS branch. The corrected binary hash is `ab60c650e58824b94ea0d3cb269c5cb56e1c714ec45052dddd43e6c6c6006320`. Its command omits `DDSFAST_HAVE_CBLAS`, `-lopenblas`, and `EIGEN_USE_BLAS`; `-llapack` remains for frozen solver routines and may transitively load OpenBLAS. This LAPACK dependency is separate from the selected Eigen dense multiplication.

   Parent's strict comparisons pass at L3 and L12 on the newly generated native grids, both against the candidate's internal reference and against frozen explicit CBLAS. Across the Eigen-versus-CBLAS comparisons, kernel L2 errors stay below `7.4e-16`, scaled pointwise errors below `3.3e-14`, and raw/processed receiver spectra and objectives match exactly in the saved comparisons. Contacts were disabled in these bounded runs. No backend latency advantage is claimed.

4. **The generated L3/L6/L12 grids pass their bounded structural checks.** Resolved bands are 6, 12, and 24; dimensions are 7×13, 13×25, and 25×49; each grid has 56 one-sided traces. Output row counts are 5,096, 18,200, and 68,600. Weight sums agree with `4π` to roundoff, radial trace metadata is exact, and zonal/longitude Gram errors are below `1.6e-14`. These are native-grid fixture checks, not a high-degree PREM run.

5. **The requested measurement math passes on the synthetic probe.** The code applies the 0–10 mHz band filter to raw bins, then applies the production window map using square-root Hann over 1000–1200 s. It forms baseline `p=P(BU)`, `q=df*eta*p`, and raw gradient `B P* q`, with one-sided DC and Nyquist weights. The v2 probe reports Parseval relative error `2.97e-15`, adjoint relative error `5.55e-15`, and central-FD relative errors `4.58e-15`, `3.43e-14`, and `2.04e-13` for `h=1e-2`, `1e-3`, and `1e-4`. Its energy is `0.05374319150371884`; `q·p=0.10748638300743736`, twice the energy to rounding. This validates the identities on the synthetic signal.

6. **The selected PREM geometry places likely S phases inside the window; the actual waveform remains unverified.** The updated receiver is 45°N, 75°E at 3 km depth. TauP PREM ray arrivals put S at 1070.84 s, sS at 1093.50 s, and ScS at 1175.59 s. The earlier 45°N, 90°E receiver put S at 1198.80 s, where the Hann window is zero. TauP arrival times do not establish a DSpec/YSpec receiver waveform or its energy derivative.

7. **DC, the upper bracket, and physical convergence remain open for any physical-kernel run.** The desired band includes DC while `FreqFull::i1()` starts at 1; its upper index uses an exclusive bracket and the first bin above 10 mHz is 0.01015625 Hz. The diagnostic candidate remains a 128-sample / 20 s fixed-weight control, not the proposed 1024-sample, 5 s-step physical energy case. Windowing mixes frequencies, so the exact `P*` may require bins beyond the band-limited signal support. The estimated `L≈114` is dimensional scale only. Degree/radial/time-padding convergence, finite-band baseline waveform and energy-derivative checks, and DSpec–YSpec agreement remain unperformed.

## Evidence reviewed

- `evidence/grid-checks.json` and generated L3/L6/L12 CSVs.
- `evidence/eigen-grid-gates.json`, corrected candidate source and binary, and `dispatch_erratum.json` (including rejected first dispatch attempt).
- `evidence/measurement-v2-stdout.log`, the v2 probe source, and `physical-check/fixture.json`.
- `evidence/prem-arrival-check.json` and the updated `checked-decisions.json`.

This PASS applies only to those bounded checks and artifacts. It does not authorize or imply a high-degree PREM solve, convergence campaign, or physical kernel performance claim.
