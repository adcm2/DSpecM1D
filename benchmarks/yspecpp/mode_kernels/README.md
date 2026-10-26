# Narrow-band modal-frequency kernel benchmarks

These are optional, reproducible benchmarks and are not registered with CI. Run `run_benchmark.py --mode toroidal` for the historical 0T2 case (the default), or add `--mode spheroidal` for the full-gravity 0S2 case. Each run requires a new empty output directory; the runner rejects nonempty destinations before writing. It caps OpenMP and BLAS at two threads and builds one adapter at a time.

Both cases use DSpecM1D and YSpecpp production sensitivity kernels, a real-linear fit to a damped complex pole with smooth complex background, and seven q measurements: six nuisance coordinates and the fitted angular frequency. Each production kernel is applied to every q weight. The frozen centered spectral tangent is compared against separate centered forward runs with fractional perturbation `h=1e-4`.

## Cases and reference

The shared homogeneous sphere has radius 6371 km, density 5515 kg m⁻³, Vp=8000 m s⁻¹, Vs=4500 m s⁻¹, and G=6.6723×10⁻¹¹ m³ kg⁻¹ s⁻². The 0T2 case uses degree 2, a pure Mrt source at (0°,0°) at 50 km depth, and a North receiver at (0°,90°) at 3 km depth. Its analytic spherical-Bessel reference satisfies `x j₂′(x)−j₂(x)=0`; the receiver geometry cancels the spheroidal component. The 0S2 case uses a pure Mrr source and vertical receiver at the same locations, isolating the spheroidal component. The DSpec family preflights use the same source and receiver geometry as the measured mode; YSpecpp computes both mode families and the vertical toroidal angular factor vanishes for 0S2.

The spheroidal reference is an independent global Rayleigh–Ritz solve of the DT weak-energy formulation, not either production solver’s radial eigensolver or matrix. The formulation follows Dahlen and Tromp (1998), *Theoretical Global Seismology*, Eqs. (8.201) and (8.208). Its regular displacement basis has the coupled leading terms U=ℓr^(ℓ−1), V=√(ℓ(ℓ+1))r^(ℓ−1), followed by independent r^(ℓ+1) Legendre terms. The potential basis is P=r^ℓ times Legendre polynomials in r². Order 10 is the primary numerical eigenmode and order 8 is the bounded basis-size check. DT analytic kernels evaluated on this Ritz mode provide the reference derivatives; its frequency is numerical, not closed form.

Using the DT-normalized V convention, let `k=√(ℓ(ℓ+1))`, `t=2U−kV`, `s=rV′−V+kU`, and `I=∫ρ(U²+V²)r²dr`. The elastic weak-energy integrand is

`C(rU′)² + 2F rU′t + (A−N)t² + Ls² + (ℓ(ℓ+1)−2)NV²`.

The mixed gravity integrand contains `4πGρ²r²U² − 2ρgrUt + 2ρ(r²UP′+krVP) + (r²P′²+ℓ(ℓ+1)P²)/(4πG)`. Its integral includes the exterior-decay boundary term `(ℓ+1)RP(R)²/(4πG)`; the potential is eliminated by stationarity. With each numerator below divided by `2ωI`, the five modulus kernels are `A:t²`, `C:(rU′)²`, `F:2rU′t`, `L:s²`, and `N:[(ℓ(ℓ+1)−2)V²−t²]`. The fixed-modulus density kernel is

`[8πGρr²U²−2grUt+2r²UP′+2krVP−ω²r²(U²+V²)−8πGr²∫ᵣᴿρ(s)U(s)t(s)/s ds]/(2ωI)`.

The density expression includes inertia, explicit gravity, and the nonlocal `δg(r)` term. A density change recomputes `g(r)` while holding elastic moduli fixed. Modulus source-point terms are reported as fractional derivatives by multiplying the objective-per-pascal source kernel by the local source modulus; density has no point term.

## Frozen processing and gates

The toroidal band is 0.24–0.32 mHz; the spheroidal band is 0.30–0.36 mHz. Both use a 72-hour record, 60-second sampling, `nt=8192`, damping `5/259200 s⁻¹`, and a 4321-sample Hann window. The positive-frequency spectrum is zero outside the computed band before inverse transform, damping removal, windowing, and forward transform. The fit applies the processing map separately to each complex pole/background basis and its imaginary multiple, preserving the map’s real-linearity. YSpecpp uses a 1 km radial mesh to represent the receiver depth break.

The spheroidal run applies separate uniform fractional changes to A,C,F,L,N and density. Each direction must pass the predeclared frequency-derivative, cross-code, radial-profile, and complete complex-spectrum tangent limits (each 1e-3), as well as fit/q identity (1e-8), baseline spectrum and fitted-frequency limits (1e-5), metadata-grid limit (5e-11), order-8/10 frequency and displacement checks, and normalized boundary-residual limit (1e-7). The toroidal run retains its original L=N reference and frozen gates. The cases cover one homogeneous sphere, one degree, one band each, and fixed solver meshes; they do not establish general model or mesh convergence.

## Dated results

All results below are from 2026-10-24. Relative derivative errors compare production `dω/dh` with the stated independent reference. The spheroidal column reports the largest error across all six directions; its separate results and every gate are in the saved JSON.

| Mode and evidence | Reference frequency | Reference derivative(s) | Max DSpec relative derivative error | Max YSpec relative derivative error | Max centered-spectrum residual | Gates |
|---|---:|---:|---:|---:|---:|---|
| 0T2, analytic spherical-Bessel; fresh regression after adapter extension | 0.281165325600 mHz | `dω/dh=ω/2=8.83306921350e-4 rad s⁻¹` for uniform L=N | 1.57e-10 | 2.96e-8 | 3.30e-5 | PASS |
| 0S2, DT full-gravity Ritz order 10; order 8 frequency change 2.96e-13, displacement overlap 1.000000000000 | 0.328735516651 mHz | A 3.15909105433e-4; C 1.56391382569e-4; F −6.99966312630e-5; L 2.61739478274e-4; N 2.76427341966e-4; ρ −8.48188269874e-4 rad s⁻¹ | 1.36e-7 | 3.77e-7 | 3.05e-5 | PASS |

For the toroidal case, the saved historical output is `build/mode-frequency-kernels-20261024/reviewed/`; the fresh post-extension default-path regression is `build/mode-frequency-kernels-20261024/toroidal-regression/`. The final spheroidal raw spectra, centered changes, q weights, reference curves, six-profile plot, actual-versus-predicted spectral plots, gates, and exact provenance are saved in `build/mode-frequency-kernels-20261024/spheroidal-accepted/`.

## Run

From the DSpecM1D repository root, pass a new empty directory for each run:

```sh
python3 benchmarks/yspecpp/mode_kernels/run_benchmark.py \
  --mode toroidal \
  --output build/mode-frequency-kernels-20261024/toroidal-reproduction-1 \
  --build-root /tmp/modal-frequency-kernels-build \
  --dspec "$PWD" --yspec "$(dirname "$PWD")/YSpecpp"
```

```sh
python3 benchmarks/yspecpp/mode_kernels/run_benchmark.py \
  --mode spheroidal \
  --output build/mode-frequency-kernels-20261024/spheroidal-reproduction-1 \
  --build-root /tmp/modal-frequency-kernels-build \
  --dspec "$PWD" --yspec "$(dirname "$PWD")/YSpecpp"
```

Each output contains `results.json`, `protocol.json`, solver spectra, q weights, radial kernels, and provenance; the spheroidal output also contains `report.html`, `summary.png`, two centered-spectrum comparison figures, six per-direction centered-spectrum CSVs, and sampled reference profiles. A required gate failure is recorded and exits nonzero. Dependency and source identities are recorded in the results JSON.
