# Read-only analysis of the saved physical degree sequence

This examines the saved matched DSpecM1D and YSpecpp PREM vertical waveforms at `Lmax=32,64,96,128,192,256,384` (5 s sampling, `N=1024`, `epsilon=5/3600 s^-1`, 0–10 mHz contour-bin taper, 1000–1200 s Hann energy). It is diagnostic analysis, not a new Earth solve or a convergence pass. The calculated values and exact run names are in `angular-convergence-analysis.json` in this directory. The parent is extending YSpec to higher degrees; no result above 384 was used here.

## What the saved spectra show

I reconstructed each DSpec window trace directly from its saved *raw* spectrum after applying the specified 8–10 mHz cosine taper:

`u_n = df exp(epsilon t_n) [Re X_0 + 2 sum_{k=1}^{51} Re{X_k exp(i 2 pi k n/1024)}]`.

The maximum absolute difference from any saved DSpec window trace is below `1.6e-20 m`. Thus the frequency decomposition below describes the actual saved waveform, including the same damping removal and filter order. DSpec and YSpec saved window traces agree at each degree to weighted relative L2 between `5.8e-6` and `1.6e-5`, with energy differences below `8e-6`; that confirms their close *same-truncation* agreement, not convergence in degree.

| Lmax | Full E (10^-9 m² s) | Full trace rel. L384 | No-DC trace rel. L384 | k>=4 trace rel. L384 | k>=16 trace rel. L384 |
|---:|---:|---:|---:|---:|---:|
| 32 | 0.284 | 0.456 | 0.430 | 0.412 | 0.389 |
| 64 | 1.083 | 0.132 | 0.130 | 0.130 | 0.121 |
| 96 | 1.239 | 0.371 | 0.134 | 0.0255 | 0.0169 |
| 128 | 1.115 | 0.268 | 0.0963 | 0.0183 | 0.0120 |
| 192 | 0.835 | 0.631 | 0.223 | 0.0425 | 0.0288 |
| 256 | 1.501 | 0.552 | 0.196 | 0.0389 | 0.0251 |
| 384 | 0.882 | 0 | 0 | 0 | 0 |

Every relative trace entry compares the **same modified observable** at that degree with L384. “No DC” removes only `k=0`; `k>=4` removes 0–3 (below 0.78125 mHz); `k>=16` removes 0–15 (below 3.125 mHz). These removals are diagnostic decompositions and do **not** redefine the authorized physical measurement. The full energies are the saved protocol energies, independently reproduced from the spectra. In particular, removing only DC leaves L192 and L256 waveform differences of 22.3% and 19.6% from L384. Removing the first four bins makes those 4.25% and 3.89%, and even keeping only `k>=16` leaves 2.88% and 2.51%. DC plus the adjacent low frequencies dominate the large oscillations, while higher bins also remain well beyond the `1e-3` waveform gate. A temporarily close energy alone is misleading: L192 full energy differs by only 5.43% from L384 while its windowed trace differs by 63.1%.

The saved raw DC values (m s) alternate with degree: `+5.2074e-4` at L32, `+9.0378e-5` at L64, `-1.1710e-3` at L96, `-8.6469e-4` at L128, `+1.8528e-3` at L192, `-1.7323e-3` at L256, and `-6.1606e-5` at L384. The `k=1` frequency is 0.1953125 mHz, close to the contour damping frequency `epsilon/(2 pi)=0.22105 mHz`; treating DC alone as an anomalous static term misses the neighbouring quasi-static response. Both solvers' close agreement at fixed L and the directly reconstructed traces support a common angular-truncation problem. They do not exclude a shared convention error or establish the final infinite-degree answer.

## Why the 10 mHz traveling-wave scale is insufficient

The protocol's `ell≈2 pi fR/c≈114` is a dimensional scale for a 10 mHz **propagating** phase. The point moment tensor has a step time dependence: DSpec displacement multiplies its impulsive source by `-i/z=1/(i z)`, while YSpec `SourceWeights::Factor` contains `1/(i omega)`. At `z=-i epsilon`, this is finite but enhances slow, quasi-static response. Neither source is spatially smoothed. In YSpec `SourceWeights.hpp`, the source coefficients contain `zeta=sqrt(ell(ell+1))` and `sqrt((ell-1)(ell+2))`; DSpec `SEMForceSpheroidal.h` has the corresponding `sqrt(ell(ell+1))` and `sqrt((ell+2)(ell-1))` factors. High degrees can therefore enter the near-field/evanescent Green function even when no 10 mHz wave propagates at those degrees.

The source and receiver radii are 6321 and 6368 km, just 47 km apart radially. A local static spherical-harmonic Green-function term has an illustrative radial envelope `(r_source/r_receiver)^ell≈exp(-ell/135)`: it is about 0.150 at L256 and 0.058 at L384, so radius separation alone offers no small L384 bound. A moment derivative can add powers of ell, but the complete modal response also has inverse-degree factors and cancellations. Consequently **the code and this envelope do not prove a peak at L200–400**; they make it plausible that important cumulative contributions extend through that range. At 58.52° angular separation, a Legendre asymptotic phase oscillates with a period of about 6.15 in degree index, so sparse cumulative L checkpoints can alternate and falsely appear stable. Dense degree-block contributions and a higher-L reference are needed for a quantitative tail bound.

## Consequence for the complete-2L kernel grid

The saved current DSpec radial-node file has **179 one-sided rows** from 35.9 km to the surface, including 28 repeated-radius rows. This is a materially larger radial domain than the earlier 56-row, 4200–6200 km fixture. `ForBand(L,2,2.0)` needs at least `(2L+1)(4L+1)` angular samples per row before longitude rounding. Lower bounds for the 179-row grid are:

| Lmax | Minimum points | 8-double result alone | 5-double binary alone |
|---:|---:|---:|---:|
| 256 | 94.1 million | 5.61 GiB | 3.51 GiB |
| 384 | 211.6 million | 12.61 GiB | 7.88 GiB |
| 512 | 375.9 million | 22.41 GiB | 14.00 GiB |
| 768 | 845.5 million | 50.39 GiB | 31.50 GiB |
| 1024 | 1.503 billion | 89.57 GiB | 55.98 GiB |
| 1536 | 3.380 billion | 201.47 GiB | 125.92 GiB |

These exclude points, intermediate fields, basis/cache data, the input grid, and additional output. The YSpec adapter estimates eight doubles per point before calling an API that materializes the sampled result; its output cap correctly prevents an unbounded allocation. Under the current 8 GiB run cap, a complete 179-side L384 grid is already infeasible as one in-memory result, and any higher converged L makes this more severe. A partial radial or angular sample may still benchmark a precisely labelled *sampled-domain* kernel, but cannot be reported as a complete Earth-kernel measurement.

The immediate sequence is to let the parent finish the bounded higher-L **waveform** checks, compare full window traces and energies against an actually higher reference, and retain degree-block/bin-group diagnostics if the oscillation persists. Only after that should a separate approved kernel design choose bounded streaming/chunking, a different output contract, or a mathematically derived acceleration. A pole-only receiver-adjoint route can reduce solve work but does not reduce the cardinality of a complete grid; analytical static-tail subtraction or radius-dependent angular truncation would need a derivation and validation against the full observable. This note authorizes none of those changes.
