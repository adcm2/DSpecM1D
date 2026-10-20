# Optional YSpecpp benchmark

This standalone benchmark builds DSpecM1D and YSpecpp in separate CMake
processes and build trees. It does not add YSpecpp to DSpecM1D's public build.
Use the harness runner with explicit paths to the DSpecM1D checkout and
immutable YSpecpp source snapshot. Generated run products belong under the
ignored `build/phase3/benchmark-results/` directory.

Run the matrix with:

```sh
python3 benchmarks/yspecpp/run_benchmark.py \
  --dspec-source "$PWD" \
  --yspec-source /path/to/yspecpp-snapshot \
  --dependency-root build/phase12-fix/fetched-pinned/_deps --jobs 1
```

It requires CMake, a C++23
compiler, Python 3, and the pinned dependency checkouts; matplotlib is optional.
The runner verifies each dependency and the YSpecpp source file-tree digest
before configuring the isolated adapters.

The harness was developed against DSpecM1D base
`65b47ccba6eae5f42776974de73974b6e6fa60cb`; the measured pre-3C baseline is
checkpoint `00e08d8f0c2ac47271d9435cb00e9fc1242705d9`. Runtime provenance is
authoritative for each run. YSpecpp is pinned to
`84be5a6cbe442324c164e5193a66b0911fca5860`. Its dependency
revisions are NumericConcepts `e03a2889d43a03660d0879f46dda1406d90c9517`,
Interpolation `d5a1711ca6a89fc8ad8148fe13212398fe569bc0`, FFTWpp
`cb152c0e85f3822493c2d9e85a56307ac7b4f1c0`, GaussQuad
`7d4786d192f68077cd99a34dcdf1379ed1eeac9a`, GSHTrans
`8dd116df8df4222fe9b10a40c7213d34f2a274fb`, and planetmodel
`fa6e230fde25152d40f27d9839de44c4f3182a24`. The three version tags declared
by YSpecpp were resolved to these exact SHAs on 2026-10-20. Keep each solver's
FetchContent closure in its own build tree; verify cached repository `HEAD`s
before reuse. The supplied YSpecpp archive SHA-256 is
`6d79288770cd6b73f8679e532d8810ddc6bbc79e8fe2a06cba04e8d0ae4183c3`; its
76-file content tree matches SHA-256
`5e77d93a8b2dbe469813ff3d190a664752158d7c66008a5b4c61a20ced38979e`.

The shared case JSON defines geocentric source and receiver latitudes, source
moment components, model, frequency samples, output quantity, gravity,
attenuation, and solver controls. The runner converts these latitudes to the
geographic convention expected by YSpecpp while preserving the local
spherical moment-tensor components. Start with homogeneous or
piecewise-constant fixtures to avoid the known within-layer interpolation
difference: YSpecpp splines elastic moduli, while DSpecM1D splines velocities
before constructing moduli. Report physical frequency, complex angular
frequency, radii, normalizations, source, and observable conventions with each
spectrum. Raw complex spectra are converted to m/s per Hz; time traces are in
m/s. The runner independently checks each native time trace against a direct
inverse Fourier sum from the saved physical complex spectrum. Raw spectra are
compared before filtering using matching integer bins after checking FFT
length, frequency spacing, complex angular frequency, and damping. The real
angular-frequency check allows `1e-10` relative to account for DSpecM1D's
existing truncated pi constant; other grid checks use `1e-12`. DSpecM1D's
gravity setting must be explicitly routed through
its `forceCowling` argument or cutoff; a parsed selector alone does not prove
which gravity operator ran. Both solvers use damped complex frequency and the
same linear-Q attenuation form, but their reference-period and finite-damping
implementation details must be recorded and checked in the measured run.

The reviewed pre-3C fixture bounds are 5e-6 for solid full-gravity cases and
1e-5 for ocean full-gravity and both Cowling fixtures. For every component,
degree 0-4 and all-degree sum, both raw complex spectra and native SI time
traces must satisfy relative L2 and peak-scaled maximum error at the same
limit. Bounds cover attenuation-off cases and their paired full-gravity
attenuation-on regression cases after Phase 3C; the bounds do not widen. The
attenuation-on pre-3C runs are diagnostic and are expected to exceed them.
These fixture-specific bounds do not change the existing matrix 1e-12, thread
1e-10, off-version 1e-9, or helper 1e-12 gates. The separate direct-inverse
transform check remains 1e-12 for both relative L2 and peak-scaled maximum.
Exact zero channels are checked for zero absolute error. This compact scope does
not establish fluid-source support or paper-scale validation. Historical
Fortran YSpec data are not YSpecpp results.

For a small convergence sweep, run `run_benchmark.py` for the same case with
`--refinement coarse`, `medium`, and `fine`, each to its own `--output`
directory, then pass those three directories to `analyze_convergence.py`.
The resolved coarse/medium/fine DSpecM1D `nq`/`maxstep` profiles are
3/0.2, 5/0.0125, and 7/0.00625; YSpecpp tolerances are 1e-3, 1e-7, and
1e-9. Medium and fine meshes separate the 100 km source from the solid
10 km-deep receiver. Settings and source hashes are saved with each result;
the analyzer compares matching physical cases, provenance, FFT bins, and time
samples.
