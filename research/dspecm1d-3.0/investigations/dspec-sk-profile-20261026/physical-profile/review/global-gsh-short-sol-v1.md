# Short mathematical and saved-array review: frozen global GSH v1

26 October 2026. I read the frozen `global-gsh/source-v1-frozen/physical_global_gsh.cpp`, the frozen v4 DSpec tensor route, GSHTrans's Wigner convention/test, and the parent-saved complex binary arrays. I did not compile, run an Earth solve, or modify the prototype.

## Coefficient orientation

The v1 rotation `A^N_{lm}=exp(-im lambda) sum_{m'} d^l_{m',m}(theta_0) a^N_{lm'}` is the correct orientation for the implemented `Y^N_{lm}=sqrt((2l+1)/(4pi)) d^l_{N,m}(theta) exp(+im phi)` convention and the DSpec source/receiver frame `R_z(lambda)R_y(theta_0)`. The `rotations[ri][mprime][l][m]` accessor places **local order** in Wigner's first/upper slot and **global order** in its second/order slot intentionally. GSHTrans's accessor is normalized by `sqrt((2l+1)/(4pi))`; v1 divides by that factor, so the rotation uses unnormalized small-d. Only `|m'|<=2` is needed by the local source (and `<=1` for the receiver adjoint), hence the rotation table's upper-index cap 2 is sufficient even though output `|m|` reaches L.

A direct l=1 scalar sign check rules out the proposed transpose. Tilt a local zonal `Y^0_{10}=C cos(theta')` by positive y angle beta. Since `cos(theta')=cos(beta)cos(theta)+sin(beta)sin(theta)cos(phi)`, and GSHTrans has `Y^0_{1,+1}=-C sin(theta)e^{+i phi}/sqrt(2)` and `Y^0_{1,-1}=+C sin(theta)e^{-i phi}/sqrt(2)`, the global coefficients are `a_{+1}=-sin(beta)/sqrt(2)` and `a_{-1}=+sin(beta)/sqrt(2)`. These are `d^1_{0,+1}` and `d^1_{0,-1}`; the transposed `d^1_{m,0}` would reverse both signs. A subsequent positive z rotation gives `exp(-im lambda)` because the local azimuth is `phi-lambda`. For spin ±1/±2, DSpec's tangent factor `(a-ib)^N` supplies the spin-frame phase in the same rotation representation; the saved full-channel tests below exercise that convention rather than only the scalar case.

## Independent check of saved complex arrays

I read both channel-major `complex128` files independently of their JSON summaries, reshaped them as 12 channels by the reported native angular point count, and computed per-channel relative L2 norms and maximum complex absolute differences against the accepted frozen tensor sampler. Source/receiver locations are asymmetric and nonpolar: (10°N,20°E) and (45°N,75°E). Frequency and model are common to candidate and reference, so this is an **angular synthesis equivalence** gate, not an independent modal/physical prediction.

| Native radial side | L | Angular points | Max channel relative L2 | Max complex absolute field difference |
|---:|---:|---:|---:|---:|
| 118, solid, r=5,600,000 m | 3 | 91 | 1.27e-15 | 8.10e-28 |
| 118, solid, r=5,600,000 m | 6 | 325 | 1.56e-15 | 1.67e-27 |
| 118, solid, r=5,600,000 m | 12 | 1,225 | 4.15e-15 | 1.42e-26 |

All 12 complex forward/adjoint channels were compared at **every** angular point; the L=12 maximum point error divided by that channel's maximum reference magnitude is 6.75e-15. The saved pointwise solid-kernel maximum absolute differences are `4.46e-57`, `1.66e-56`, and `4.84e-56` respectively. These nonpolar results resolve the frozen v1 Wigner-index concern for the tested solid cases. They do not establish high-degree numerical stability, timing improvement, all radial sides, or physical degree/FFT convergence.

## Fluid representation is the remaining parity failure

The completed L=12 run on native side 175 (`r=6,369,072.1527 m`, element 29, fluid) fails 12-channel parity: forward rr/trace relative L2 errors are 73.7/112.2; the accepted reference's spin ±1/±2 channels are exactly zero while v1's are finite. This has a specific source-level cause. In frozen v4 `sampleForKernelsTensor`, fluid nonradial pole coefficients are converted **before angular synthesis** to component 0 = `rawStrain[0]+rawStrain[1]` (pressure/divergence representation), and components 1–5 = zero; the radial l=0 rr/trace record is subsequently added separately. V1 captures the same pole records but rotates/synthesizes all six unprojected raw components. Comparing those six outputs to the projected tensor reference is therefore not like-for-like. To match the accepted field representation, apply that same fluid coefficient projection before the GSH rotation/synthesis and retain the existing separate l=0 addition. This is an implementation recommendation, not a change made in this review.

The **physical fluid contraction already agrees** because `Windowed3DKernelEngine::pointKernel` uses only `rr+horizontalTrace` for fluid `kappa`. Independent binary calculation gives candidate-versus-reference divergence maximum absolute differences `4.83e-27` forward and `8.56e-31` adjoint (relative L2 `3.23e-14` and `2.61e-14`). Saved fluid `kappa` pointwise maximum absolute error is `5.57e-57` against maximum reference `2.92e-44`. Thus fluid channel representation is a parity/contract problem, while this saved one-frequency fluid kernel happens to pass. The fluid `result.json` also writes unquoted `inf` for zero-reference relative errors, so it is not strict JSON despite exit code zero; finite/null serialization is needed for a machine-readable gate.

After the fluid projection and JSON fix, rerun an all-12-channel fluid gate and a duplicated interface-side pair. Continue to treat the `polar_conventions_runtime_validated:false` metadata in v1 as stale relative to these **nonpolar** saved gates; it is not evidence about exact coordinate-pole behavior. A complete 179-side, multi-frequency physical kernel still requires separate resource, output, and convergence evidence.

Key source anchors: frozen global GSH v1 `physical_global_gsh.cpp` rotation and inverse sections (~147–178, 292–365); frozen v4 `Windowed3DFields.h` fluid projection (~1757–1763 and 1788–1794), radial l=0 scatter (~1821–1840), and `Windowed3DKernels.h` fluid `pointKernel` (~253–267); cached GSHTrans `Wigner.hpp` and `tests/CheckWignerConvention.hpp`.
