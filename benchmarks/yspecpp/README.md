# Optional YSpecpp benchmark

## Published-style Phase 3 suite

`run_paper_suite.py --case fig1|fig3|fig4` runs the new matched-input
comparisons; `run_fig5_record.py` runs the high-frequency, 91-receiver record
section. Each comparison records its own settings, source/model hashes,
adapter logs, binaries, and harness source snapshots under the ignored
`build/phase3/paper-suite-20261020/` directory. `--pilot` selects a bounded
top-frequency diagnostic and is never labeled as a paper result. Figure 1's
raw Z spectrum is primary; its ex1-style pre-window/filter/final-window
spectrum is a separately reported secondary product. Figures 3 and 4 report
acceleration in SI units, with full-window and 400–1300 s metrics. Figure 5
adopts acceleration from ex4 because the paper does not state the plotted
quantity; it reports record-section products without a cross-solver score.
Figure 2 is the separately completed `run_paper2.py` comparison. The outputs
are new published-style comparisons, not a pixel-for-pixel reproduction of
the paper's historical solver figures; current DSpecM1D and YSpecpp use
different continuous PREM interpolation representations.

The reviewed case matrix is tracked at
`cases/paper-suite-protocol.json`; this is the default for the Figure 1, 3, 4,
and 5 runners, so a clean checkout does not depend on ignored build data.
Figure 5 remains blocked: DSpecM1D refuses the required degree 7500 because
its current Wigner recursion is only safe through degree 1827. A separate,
explicitly non-paper three-receiver pilot capped at degree 1827 exited with
signal 11 before producing result files, so no supported-degree runtime
estimate is claimed. Both failure records are preserved under the ignored
suite output; no production code or angular guard was changed.

Figure 6's internal convergence sweep uses `run_convergence_benchmark.py`
(`--pilot` for the bounded cost check, `--full` for the sequential order/mesh
sweep). Figure 7's isolated degree-2 tidal-resolution study uses
`run_tidal_benchmark.py`; it preserves the old ex6 example and reports the
normalized response without claiming an absolute forcing amplitude. Full
Figure 1–4 comparisons are complete, Figure 5 is blocked as described above,
Figure 6 has a bounded high-frequency pilot while its multi-hour convergence
sweep awaits an execution preference, and Figure 7 has completed isolated
degree-2 tidal runs. The suite and Figure 5 runners refuse to overwrite
nonempty run directories; choose a new output root for a repeat measurement.

## Paper Figure 2 comparison

`run_paper2.py` writes an isolated comparison under
`build/phase3/paper2-20261020/`. Its default is a narrow-band diagnostic pilot
with degrees 0–300 and the paper's five-hour grid, nq=5, and 1e-5 DSpecM1D
accuracy. Run with `--full` only to calculate the published 1.5–21 mHz band.
The outputs include raw physical spectra, shared-filter SI traces, per-channel
misfit metrics, provenance, and an overlay/residual plot. The DSpecM1D and
YSpecpp PREM interpolants differ because they interpolate different material
quantities; this comparison records that model difference.

The runner reads the tracked `cases/paper2-protocol.json`, records the current
DSpecM1D revision in provenance, and does not require the ignored historical
protocol file under `build/`. It still needs an external immutable YSpecpp
checkout and a verified dependency cache; pass `--yspec-source` and
`--dependency-root` when their default locations are unavailable. Use
`--protocol` to select an alternate protocol; if that file pins `source_commit`,
the runner also enforces that exact DSpecM1D revision.

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
