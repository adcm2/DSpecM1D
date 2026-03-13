# DSpecM1D

NOTE: The current release should be viewed as a pre-release version of the software. API changes and a general cleaning up of the library will be complete soon.

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

---

## Dependencies

| Library | Fetched automatically |
|---|---|
| [Eigen 3.4](https://gitlab.com/libeigen/eigen) | ✓ |
| [FFTWpp](https://github.com/da380/FFTWpp) | ✓ |
| [GaussQuad](https://github.com/da380/GaussQuad) | ✓ |
| [GSHTrans](https://github.com/da380/GSHTrans) | ✓ |
| [PlanetaryModel](https://github.com/da380/PlanetaryModel) | ✓ |
| [TomographyModels](https://github.com/adcm2/TomographyModels) | ✓ |
| [Interpolation](https://github.com/da380/Interpolation) | ✓ |
| [EarthMesh](https://github.com/adcm2/EarthMesh) | ✓ |
| [SpectraSolver](https://github.com/adcm2/SpectraSolver) | ✓ |
| [FFTW3](https://www.fftw.org/) | **System** |
| OpenMP | **System** |

FFTW3 must be installed on your system before building. On Ubuntu/Debian:

```bash
sudo apt install libfftw3-dev
```

On macOS with Homebrew:

```bash
brew install fftw
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
| `BUILD_BENCHMARKS` | `ON` | Build the benchmark executables |
| `BUILD_TUTORIALS` | `ON` | Build the tutorial executables |

```bash
cmake .. -DBUILD_BENCHMARKS=OFF -DBUILD_TUTORIALS=ON
```

---

## Quick Start

After building, the tutorial binaries are placed in `build/bin/`. Tutorial 1
runs a complete seismogram synthesis for a single earthquake and receiver
geometry defined in the parameter file:

```bash
cd build
./bin/t1
```

The parameter file is expected at `build/data/params/ex1.txt`. See
[`tutorials/t1.cpp`](tutorials/t1.cpp) for a fully documented walkthrough of
the pipeline.

---

## Parameter File Format

The input file controls all aspects of the simulation. Key fields are:

```
earth_model    prem.txt        # path relative to data/
lmax           200             # maximum angular degree
NQ             5               # GLL quadrature points per element
f1  0.3  f2  10.0              # passband corners (mHz)
f11 0.2  f12 0.5               # lower taper corners (mHz)
f21 9.0  f22 11.0              # upper taper corners (mHz)
dt  1.0                        # time step (s)
t_out  60.0                    # output duration (minutes)
output_type  2                 # 0=displacement, 1=velocity, 2=acceleration
```

---

## Repository Layout

```
DSpecM1D_Draft/
├── DSpecM1D/
│   └── src/
│       ├── spectra_master.h   # Core SEM class (specsem)
│       ├── mesh_model.h       # Material parameter interpolation
│       ├── input_parser.h     # Parameter file reader
│       └── ...
├── benchmarks/                # Benchmark executables
├── tutorials/
│   └── t1.cpp                 # Tutorial 1: basic seismogram synthesis
├── data/
│   ├── params/                # Input parameter files
│   └── prem.txt               # Reference Earth model
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

These are encapsulated in the `prem_norm<FLOAT>` helper class defined in each
executable.

---

## License

To be added.

## Contact

To be added.
