# DSpecM1D

A C++23 spectral element library for computing synthetic seismograms in a
spherically symmetric (1D) Earth model. The solver works in the
frequency domain, assembling sparse stiffness and inertia matrices for each
angular degree _ℓ_ and solving the resulting linear system at each frequency, before performing an inverse Fourier-Laplace transform.

---

## Features

- Spheroidal, toroidal, and radial mode computation
- Sparse frequency-domain SEM solver with OpenMP parallelisation
- Centroid-moment-tensor (CMT) source representation
- Arbitrary receiver geometry (Z / T / R components)
- Hann-window bandpass filtering pipeline
- Non-dimensionalisation against a user-supplied reference model (PREM by default)
- Clean dependency management via CMake `FetchContent`
- Modular self-contained unit/component tests for regression catching
- Protected paper-reproduction examples with external-code comparisons intact

---

## Dependencies

| Library | Fetched automatically |
|---|---|
| [Eigen 3.4](https://gitlab.com/libeigen/eigen) | ✓ |
| [FFTWpp](https://github.com/da380/FFTWpp) | ✓ |
| [GaussQuad](https://github.com/da380/GaussQuad) | ✓ |
| [GSHTrans](https://github.com/da380/GSHTrans) | ✓ |
| [Interpolation](https://github.com/da380/Interpolation) | ✓ |
| [FFTW3](https://www.fftw.org/) | **System** |
| BLAS / LAPACK | **System** |
| OpenMP | **System** |

FFTW3 plus a BLAS/LAPACK implementation must be installed on your system before building. On Ubuntu/Debian:

```bash
sudo apt install build-essential cmake git libfftw3-dev libopenblas-dev liblapack-dev
```

On macOS with Homebrew:

```bash
brew install fftw openblas
```

NetCDF is not required by DSpecM1D. OpenMP support is provided by the compiler
toolchain; on Ubuntu/Debian, GCC normally provides it through `libgomp`.

`PlanetaryModel` and `EarthMesh` are now maintained in-house within the
DSpecM1D source tree under `DSpecM1D/src/model_info/`.
The frequency-axis and FFT/post-processing utilities formerly provided via
`SpectraSolver` are now maintained in-house under
`DSpecM1D/src/frequency_info/`.

---

## Building

Quick start:

```bash
git clone https://github.com/adcm2/DSpecM1D.git
cd DSpecM1D
cmake -S . -B build
cmake --build build --parallel
```

CMake will automatically download all header-only dependencies during the
configure step. An internet connection is required on first build.

If you switch CMake generators for an existing build tree, use a new build
directory or remove the old one first. CMake build trees are generator-specific.

### Presets

The repository includes three configure/build/test presets:

- `dev`: `RelWithDebInfo` with tests, tutorials, benchmarks, and smoke checks
- `release`: optimized release build with tests, tutorials, benchmarks, and smoke checks
- `debug`: debug build with tests, tutorials, benchmarks, and smoke checks

Example preset usage:

```bash
cmake --preset dev
cmake --build --preset dev
ctest --preset dev -L "unit|component"
ctest --preset dev -L dspecm1d_smoke
```

For a release or debug tree, replace `dev` with `release` or `debug`.

### Build options

| Option | Default | Description |
|---|---|---|
| `DSPECM1D_BUILD_BENCHMARKS` | `ON` | Build the benchmark executables |
| `DSPECM1D_BUILD_TUTORIALS` | `ON` | Build the tutorial executables |
| `DSPECM1D_BUILD_TESTS` | `OFF` | Build the self-contained unit/component test suite |
| `DSPECM1D_ENABLE_SMOKE_TESTS` | `OFF` | Register runtime smoke checks such as `t1` |
| `DSPECM1D_ENABLE_MIGRATION_TESTS` | `OFF` | Register temporary migration-only checks such as the `ex2` tolerance comparison |
| `DSPECM1D_ENABLE_PAPER_VALIDATION` | `OFF` | Register optional paper-reproduction validation checks |
| `DSPECM1D_BUILD_DOCS` | `OFF` | Enable Doxygen and static website build targets |

For example, to build only the tutorial executables:

```bash
cmake -S . -B build/tutorials-only \
  -DDSPECM1D_BUILD_BENCHMARKS=OFF \
  -DDSPECM1D_BUILD_TUTORIALS=ON
cmake --build build/tutorials-only --parallel
```

---

## Quick Start

After building, the tutorial binaries are placed in `build/bin/`. Tutorial 1
(`t1`) is the self-contained quickstart and runtime smoke example:

```bash
cd build
./bin/t1
```

The parameter file is expected at `build/data/params/t1.txt`. See
[`tutorials/t1.cpp`](tutorials/t1.cpp) for a fully documented walkthrough of
the pipeline.

---

## Testing

DSpecM1D separates small regression-catching tests from heavier
paper-reproduction workflows.

- `tests/`: modular, self-contained unit/component checks for parsers,
  helpers, readers, filters, writers, and small SEM invariants
- `t1`: self-contained runtime smoke example
- `ex1` through `ex7`: protected paper-reproduction examples that intentionally
  keep their YSpec / MinEOS / SpecNM comparisons

Manual development build:

```bash
cmake -S . -B build/dev \
  -DDSPECM1D_BUILD_TESTS=ON \
  -DDSPECM1D_BUILD_TUTORIALS=ON \
  -DDSPECM1D_BUILD_BENCHMARKS=ON \
  -DDSPECM1D_ENABLE_SMOKE_TESTS=ON
cmake --build build/dev --parallel
ctest --test-dir build/dev --output-on-failure -L "unit|component"
ctest --test-dir build/dev --output-on-failure -L dspecm1d_smoke
```

Equivalent preset workflow:

```bash
cmake --preset dev
cmake --build --preset dev
ctest --preset dev -L "unit|component"
ctest --preset dev -L dspecm1d_smoke
```

Temporary migration-only validation, if enabled:

```bash
ctest --test-dir build/dev --output-on-failure -L migration
```

For manual tutorial or benchmark runs, cap OpenMP explicitly if you do not want
the solver to use all available cores:

```bash
export OMP_NUM_THREADS=2
```

---

## Parameter File Format

The input file controls all aspects of the simulation. The current parser reads
an ordered sequence of value lines, not a keyed `name value` format. Comments
and blank lines are ignored, but the remaining values must appear in the
expected order. Full-line comments beginning with `#` are supported; trailing
inline comments on value lines are not.

```text
# prefix for output files
"./output/t1.out"

# Earth model
"models/prem.200.noatten.txt"

# which type of modes to include
4

# attenuation switch: 1 = on, 0 = off
0

# gravitation: 0 = none, 1 = cowling, 2 = self
2

# output: 0 = displacement, 1 = velocity, 2 = acceleration
0

# potential and tilt corrections: 0 = no, 1 = yes
0

# relative error
1e-3

# lmin
0

# lmax
150

# fmin (mHz)
1.5

# fmax (mHz)
10.0

# length of time series (min)
240

# time step (sec)
1.0

# f11 filter (mHz)
1.9

# f12 filter (mHz)
2.0

# f21 filter (mHz)
9.0

# f22 filter (mHz)
10.0

# source depth (km)
647.1

# source latitude (deg)
-13.82

# source longitude (deg)
-67.25

...
```

See [data/params/t1.txt](/home/adcm2/Documents/c++/DSpecM1D_Draft/data/params/t1.txt) for a fully commented example in the exact parser order.

---

## Examples

- `t1` is the recommended first executable and the default runtime smoke test.
- `ex1` through `ex7` are paper-reproduction examples. Their external
  comparison inputs are intentional and are preserved so the workflows continue
  to reproduce the results discussed in the accompanying paper.
- The canonical YSpec, MinEOS, and SpecNM comparison files used by the
  protected examples are stored in `data/reference/` and copied into the build
  tree with the rest of the project data.

## Repository Layout

```
DSpecM1D_Draft/
├── DSpecM1D/
│   ├── ModelInput             # Public Earth model entry point
│   ├── EarthMesh              # Public radial-mesh entry point
│   └── src/
│       ├── SpectraMaster.h    # Core SEM orchestration
│       ├── InputParser.h      # Parameter file reader
│       ├── InputParametersNew.h
│       ├── SignalFiltering.h
│       ├── frequency_info/
│       │   ├── FreqFull.h
│       │   └── PostprocessFunctions.h
│       ├── model_info/
│       │   ├── ModelInput.h
│       │   ├── ModelConcepts.h
│       │   ├── MeshModel.h
│       │   └── earthmesh/
│       └── ...
├── benchmarks/                # Benchmark executables
│   └── ex1.cpp ... ex7.cpp    # Paper-reproduction examples
├── tutorials/
│   └── t1.cpp                 # Tutorial 1: basic seismogram synthesis
├── docs/                      # Static website and Doxygen configuration
├── data/
│   ├── params/                # Input parameter files
│   └── models/                # Reference Earth model files (PREM variants)
├── cmake/                     # CMake helper modules
└── CMakeLists.txt
```

---

## Normalisation

All quantities are non-dimensionalised against three reference scales:

| Scale | Value |
|---|---|
| Length | 1 km = 1000 m |
| Density | ρ_ref = 5515 kg m⁻³ |
| Time | `1 / √(π G ρ_ref)` |

These are encapsulated by the normalization helper in `DSpecM1D/src/NormClass.h`.

---

## Documentation

Enable documentation targets with:

```bash
cmake -S . -B build/docs -DDSPECM1D_BUILD_DOCS=ON
cmake --build build/docs --target website
```

This builds a static website into `build/docs/site/` and, when Doxygen is
available, generates API reference pages under `build/docs/site/api/`.

To preview the generated site locally:

```bash
cmake --build build/docs --target preview-website
```

The preview server runs at `http://127.0.0.1:8000`.

The repository also includes a dedicated GitHub Pages workflow in
`.github/workflows/pages.yml` that builds the same `website` target and deploys
the generated `build/pages/site/` artifact. To activate it, set the repository
Pages source to `GitHub Actions` in `Settings > Pages`.

---

## Repository

- Repository: `https://github.com/adcm2/DSpecM1D`
- Issues: `https://github.com/adcm2/DSpecM1D/issues`
- Contact: use the issue tracker for release feedback and usage questions

## License

DSpecM1D is released under the GNU General Public License v3.0.
See [`LICENSE`](LICENSE).
