# Optional YSpecpp benchmark

This harness compares DSpecM1D and YSpecpp velocity spectra and native time
traces for the six compact solid and ocean fixtures in `cases/cases.json`.
Each fixture runs degrees 0–4 separately and as a sum. It checks raw complex
spectra before filtering and independently verifies each native time trace.

```sh
python3 benchmarks/yspecpp/run_benchmark.py \
  --dspec-source "$PWD" \
  --yspec-source /path/to/yspecpp-snapshot \
  --dependency-root /path/to/verified-dependencies \
  --jobs 1
```

The runner requires CMake, a C++23 compiler, Python 3, and a clean immutable
YSpecpp checkout. It verifies the 76-file YSpecpp tree and dependency revisions
before configuring separate solver build trees. Results go to
`build/phase3/benchmark-results/`; pass `--output` for another empty directory.
Use `--case ID` to select a fixture. Use `--refinement fine` for the measured
acceptance gates; coarse and medium runs can serve as convergence diagnostics.
The runner records the resolved refinement profile.

YSpecpp is pinned to `84be5a6cbe442324c164e5193a66b0911fca5860`, with tree
SHA-256 `5e77d93a8b2dbe469813ff3d190a664752158d7c66008a5b4c61a20ced38979e`.
The dependency pins are recorded in `run_benchmark.py`; the YSpecpp archive
SHA-256 is `6d79288770cd6b73f8679e532d8810ddc6bbc79e8fe2a06cba04e8d0ae4183c3`.

The reviewed limits for each component and degree are relative L2 and
peak-scaled maximum error of `5e-6` for solid full-gravity fixtures and
`1e-5` for ocean full-gravity and Cowling fixtures. They cover attenuation-off
cases and their paired attenuation-on full-gravity regressions. Report which
refinement was measured; coarse and medium runs can serve as convergence
diagnostics. Direct inverse-transform checks use `1e-12` for both metrics;
exact-zero channels require zero absolute error. Frequency bins, complex
angular frequencies, damping, and time grids are checked separately.

The homogeneous fixtures avoid the known difference between YSpecpp's modulus
spline and DSpecM1D's velocity spline. This compact comparison does not claim
paper-scale validation or establish fluid-source support.

Each benchmark run also writes `report.html` beside the JSON, CSV, and plots.
It summarizes the resolved physical and solver setup, measured degree and
summed comparisons, inverse-transform checks, measured limits, and source
identities. Use the command above as the single benchmark entry point; the
report records the resolved command and selected refinement.

### Reproducing paper comparisons

Replay one saved case without rerunning either solver:

```sh
python3 benchmarks/yspecpp/run_paper.py --saved-run /path/to/old-run --case fig1 --output /path/to/replay
```

Omit `--case` to replay Figures 1–4 into `OUTPUT/fig1` through `OUTPUT/fig4`.
For a new full run, give one explicit case and all source/dependency paths:

```sh
python3 benchmarks/yspecpp/run_paper.py \
  --case fig1 \
  --dspec-source /path/to/DSpecM1D \
  --yspec-source /path/to/YSpecpp \
  --dependency-root /path/to/dependencies \
  --output /path/to/new-run
```

`--pilot` is available only for new runs. Each case writes `report.html`,
`summary.json`, `provenance.json`, raw and processed spectrum plots, and
comparison plots; a display-window plot is included when the case defines one.
The report labels saved-data replay and links parent provenance and raw-spectrum
hashes. Raw spectra and
native traces are saved under `results/`. New summaries use the consolidated
runner's format. Existing web artifacts retain their original run derivation;
the corrected Figure 1 secondary spectrum is a separately derived product from
saved raw spectra, with its own input hashes and derivation record.

For the focused matched full-gravity elastic modulus-kernel comparison, see
[`kernels/README.md`](kernels/README.md). It documents the bounded homogeneous
fixture, scalar/window mapping, units, gates, and exact reproduction command.

For optional toroidal and spheroidal modal-frequency kernel benchmarks against
the independent DT references (not CI), see
[`mode_kernels/README.md`](mode_kernels/README.md).
