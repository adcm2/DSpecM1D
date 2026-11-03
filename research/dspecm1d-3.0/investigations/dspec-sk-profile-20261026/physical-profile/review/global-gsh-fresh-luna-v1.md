# Independent review: frozen global GSH v1

Reviewed the exact `source-v1-frozen/physical_global_gsh.cpp` against frozen DSpec v4 headers, GSHTrans Wigner/grid APIs, and the parent-saved low-degree gates. I did not compile or run binaries. Both frozen source hashes and both frozen binary hashes validate against their manifests. The saved runs report the same profiled binary hash as the frozen manifest.

## Rotation concern resolved by the saved gates

My first static derivation incorrectly treated the coefficient rotation as transposed. The saved non-polar complete-grid gates compare the frozen adapter against the accepted tensor sampler for source (10°N, 20°E) and receiver (45°N, 75°E) frames. At solid side 118, all 12 complex channels agree to near floating-point precision:

| Degree | Channel relative-L2 error range | Grid points | Kernel max absolute error range |
|---:|---:|---:|---:|
| 3 | 7.62e-16–1.27e-15 | 91 | 0–4.46e-57 |
| 6 | 9.29e-16–1.56e-15 | 325 | 0–1.66e-56 |
| 12 | 2.17e-15–4.15e-15 | 1,225 | 0–4.84e-56 |

The source and receiver frames are both asymmetric and away from a pole. These gates therefore resolve the convention concern for the sampled model, frequency, source/receiver setup, solid side, and tested degrees. The Wigner access `rotations[ri][mprime][l][m]` intentionally stores the local source order in the upper-index slot (bounded by 2 for this source) and the output global order in the lower-index slot; the tested results agree with the accepted route. The table's upper-index cap of 2 is intentional. The source result metadata still says `polar_conventions_runtime_validated:false` (`physical_global_gsh.cpp:399`); that field is now stale relative to these saved external gates and should not be read as an up-to-date run verdict.

## Fluid-side comparison caveat

At L=12, side 175 is fluid. Its accepted tensor route projects fluid coefficient data to `rr+horizontalTrace` and zeros the other channels before synthesis (`Windowed3DFields.h:1761-1763`), while this adapter synthesizes its captured channel coefficients independently. Consequently, channel comparison is not like-for-like: the run reports forward `rr` and trace relative errors 73.7 and 112.2, and infinite relative errors for zero-reference shear channels. The pointwise fluid kernel uses only the compressional combination and its absolute errors are near zero, but this does not establish six-channel fluid parity. Also, the fluid run's `result.json` contains the token `inf` in the relative-error array, which is not valid JSON and fails strict parsing. Fluid channel metrics need projection-aware comparison and non-finite JSON serialization before they can be consumed as machine-readable parity evidence.

## Other static checks

- **Six strain channels and spins:** channel ordering `rr`, horizontal trace, radial ±, horizontal ± and spins `{0,0,+1,-1,+2,-2}` match `SpinStrainField` and the accepted strain rotation (`Windowed3DFields.h:34-36, 508-511`). Kernel contraction consumes the same channels and pairings (`Windowed3DKernels.h:534-548`). Solid-side numerical parity covers all twelve forward/adjoint channels.
- **Radial degree zero:** the visitor creates explicit `(m,l)=(0,0)` radial records carrying only `rr` and horizontal trace (`Windowed3DFields.h:1929-1949`). The adapter captures both records and adds them only to spin-zero channel fields after inverse transforms (`physical_global_gsh.cpp:283-287, 346-348`). This matches the native degree-zero radial accumulation (`Windowed3DFields.h:880-908`); the solid-side gates include these channels.
- **Angular grid and comparison shape:** `ForBand(L,2,2.0)` resolves through degree `2L`, giving `2L+1` Gauss-Legendre colatitudes and `FastFFTSize(4L+1)` longitudes (`GaussLegendreGrid.hpp:109-123`; `SphericalGrid.hpp:104-115`). The adapter constructs the same GL nodes via `acos(-x)`, supplies that longitude count, and uses the same `i*2π/nPhi` longitude formula (`physical_global_gsh.cpp:112-126, 208-209, 246-261`; `SphericalGrid.hpp:236-240`). For L≤64 it checks theta nodes, weights, and dimensions at runtime; longitude equality follows from the shared formula rather than a separate runtime comparison. The saved solid gates compare every channel at every angular point for one selected one-sided radial trace and compare pointwise kernels on the same points (`physical_global_gsh.cpp:301-313, 349-362, 370-385`). These are same-grid parity checks, not full physical-convergence evidence.
- **BLAS boundary:** the frozen build metadata leaves `GSHTRANS_HAVE_BLAS`, `DSFAST_HAVE_CBLAS`, and `EIGEN_USE_BLAS` undefined and records Eigen internal math. Both binaries still link `-lblas -llapack`, attributed in the manifest to inherited DSpec Fortran solver symbols. The GSH synthesis/contraction path is not configured to call external BLAS, but the executable is not literally BLAS-free.
- **Resources:** complete-grid field output scales as O(L²), and the known-array estimate includes two Wigner rotation tables, coefficient storage, transform buffers, and optional point/reference arrays (`physical_global_gsh.cpp:208-229`). Metadata explicitly excludes solver and transform scratch and disclaims a peak-RSS guarantee (`:398`). I found no O(L² Ngrid) angular basis allocation in this adapter. This estimate is not a hard memory bound.

## Evidence boundary

The saved solid-side gates establish channel and kernel parity for this model, source/receiver setup, selected native radial trace, and L=3, 6, 12. The fluid-side run establishes pointwise kernel parity for that sample but not six-channel parity; its relative-error artifact is also invalid JSON. Neither set establishes physical convergence or a performance benefit.
