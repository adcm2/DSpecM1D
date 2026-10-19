# DSpecM1D

[![DOI](https://zenodo.org/badge/1170173312.svg)](https://doi.org/10.5281/zenodo.19793472)

DSpecM1D is a C++23 spectral-element library for synthetic seismograms in a spherically symmetric Earth model. The quickest route is to build `generate_seismogram`, run the bundled one-receiver input, and plot its output.

- First calculation: [docs/site/first-seismogram.html](https://adcm2.github.io/DSpecM1D/first-seismogram.html)
- Install and build: [docs/site/installation.html](https://adcm2.github.io/DSpecM1D/installation.html)
- C++ library example: [docs/site/library-example.html](https://adcm2.github.io/DSpecM1D/library-example.html)
- Parameter reference: [docs/site/parameter-files.html](https://adcm2.github.io/DSpecM1D/parameter-files.html)
- Generated docs site: [adcm2.github.io/DSpecM1D](https://adcm2.github.io/DSpecM1D/)

## Linux quick start

The primary documented build is Linux. Install a C++23/OpenMP-capable compiler, CMake 3.24 or newer, Git, and FFTW3 development files. The first CMake configure fetches the remaining dependencies from pinned revisions and needs network access.

```bash
sudo apt update
sudo apt install build-essential cmake git libfftw3-dev

git clone https://github.com/adcm2/DSpecM1D.git
cd DSpecM1D
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release \
  -DDSPECM1D_BUILD_TOOLS=ON -DDSPECM1D_BUILD_TUTORIALS=ON
cmake --build build/release --target generate_seismogram --parallel 2
OMP_NUM_THREADS=2 build/release/bin/generate_seismogram \
  build/release/data/params/first_seismogram.txt
```

The CLI prints the generated file path: `build/release/output/first_seismogram.out_t.out`. The example uses attenuation off, one receiver, and a modest degree/frequency range. It demonstrates the workflow; it does not establish numerical convergence. Use `OMP_NUM_THREADS` to choose a modest runtime thread count independently of `--parallel` build jobs.

## Build choices

Use `build/release` for the documented user build and `build/dev` for development with tests and examples enabled. Keep separate CMake binary directories for distinct settings. Manual configuration requires CMake 3.24 or newer; the supplied preset file uses schema version 6 and requires CMake 3.25 or newer.

```bash
# Development configuration
cmake -S . -B build/dev \
  -DDSPECM1D_BUILD_TESTS=ON -DDSPECM1D_BUILD_TUTORIALS=ON \
  -DDSPECM1D_BUILD_BENCHMARKS=ON -DDSPECM1D_BUILD_TOOLS=ON \
  -DDSPECM1D_ENABLE_SMOKE_TESTS=ON
cmake --build build/dev --parallel 2
```

For test options, configured targets, and the current CTest inventory, see [the testing guide](https://adcm2.github.io/DSpecM1D/library-testing.html) and run `ctest --test-dir build/dev -N`.

## Platform and dependencies

Linux is the primary documented platform. macOS is best-effort with a C++23 compiler, CMake, FFTW, and a matching OpenMP runtime such as `libomp`; check that CMake finds OpenMP. Windows has not been verified. NetCDF and LAPACK are not required. BLAS is not a required direct project dependency with the pinned GSHTrans configuration.

## License

GNU General Public License v3.0; see [LICENSE](LICENSE).
