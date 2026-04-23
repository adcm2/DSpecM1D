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
| [PlanetaryModel](https://github.com/da380/PlanetaryModel) | ✓ |
| [Interpolation](https://github.com/da380/Interpolation) | ✓ |
| [EarthMesh](https://github.com/adcm2/EarthMesh) | ✓ |
| [SpectraSolver](https://github.com/adcm2/SpectraSolver) | ✓ |
| [FFTW3](https://www.fftw.org/) | **System** |
| BLAS / LAPACK | **System** |
| OpenMP | **System** |

FFTW3 plus a BLAS/LAPACK implementation must be installed on your system before building. On Ubuntu/Debian:

```bash
sudo apt install libfftw3-dev libopenblas-dev liblapack-dev
```

On macOS with Homebrew:

```bash
brew install fftw openblas
```

---

## Building

```bash
git clone https://github.com/adcm2/DSpecM1D.git
cd DSpecM1D
mkdir build && cd build
cmake ..
cmake --build . --parallel
```

CMake will automatically download all header-only dependencies during the
configure step. An internet connection is required on first build.

### Build options

| Option | Default | Description |
|---|---|---|
| `DSPECM1D_BUILD_BENCHMARKS` | `ON` | Build the benchmark executables |
| `DSPECM1D_BUILD_TUTORIALS` | `ON` | Build the tutorial executables |
| `DSPECM1D_BUILD_TESTS` | `OFF` | Build the self-contained unit/component test suite |
| `DSPECM1D_ENABLE_SMOKE_TESTS` | `OFF` | Register runtime smoke checks such as `t1` |
| `DSPECM1D_ENABLE_PAPER_VALIDATION` | `OFF` | Register optional paper-reproduction validation checks |
| `DSPECM1D_BUILD_DOCS` | `OFF` | Enable Doxygen and static website build targets |

```bash
cmake .. -DDSPECM1D_BUILD_BENCHMARKS=OFF -DDSPECM1D_BUILD_TUTORIALS=ON
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

Example development build:

```bash
cmake -S . -B build/dev -G Ninja \
  -DDSPECM1D_BUILD_TESTS=ON \
  -DDSPECM1D_BUILD_TUTORIALS=ON \
  -DDSPECM1D_BUILD_BENCHMARKS=ON \
  -DDSPECM1D_ENABLE_SMOKE_TESTS=ON
cmake --build build/dev --parallel
ctest --test-dir build/dev --output-on-failure -L "unit|component"
ctest --test-dir build/dev --output-on-failure -L smoke
```

---

## Parameter File Format

The input file controls all aspects of the simulation. The current parser reads
an ordered sequence of value lines, not a keyed `name value` format. Comments
and blank lines are ignored, but the remaining values must appear in the
expected order:

```
"./output/yspec.lf.out"
"models/prem.200.noatten.txt"
4
0
2
0
0
1e-5
0
100
0.1
7.0
6000
1.0
0.1
0.2
4.9
5.0
647.1
-13.82
-67.25
...
```

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
│   └── src/
│       ├── SpectraMaster.h    # Core SEM orchestration
│       ├── MeshModel.h        # Material parameter interpolation
│       ├── InputParser.h      # Parameter file reader
│       ├── InputParametersNew.h
│       ├── SignalFiltering.h
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
cmake -S . -B build/docs -G Ninja -DDSPECM1D_BUILD_DOCS=ON
cmake --build build/docs --target website
```

This builds a static website into `build/docs/site/` and, when Doxygen is
available, generates API reference pages under `build/docs/site/api/`.

To preview the generated site locally:

```bash
cmake --build build/docs --target preview-website
```

The preview server runs at `http://127.0.0.1:8000`.

---

## Repository

- Repository: `https://github.com/adcm2/DSpecM1D`
- Issues: `https://github.com/adcm2/DSpecM1D/issues`
- Contact: use the issue tracker for release feedback and usage questions

## License

DSpecM1D is released under the GNU General Public License v3.0.
See [`LICENSE`](LICENSE).
