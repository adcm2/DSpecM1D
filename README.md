# DSpecM1D

DSpecM1D is a C++23 spectral element library for computing synthetic
seismograms in a spherically symmetric (1D) Earth model.

The README is a quick-start guide for the first public release. Full
documentation, testing notes, debugging guidance, and release-facing API
reference live on the website:

- Live docs: `https://adcm2.github.io/DSpecM1D/`
- Repository: `https://github.com/adcm2/DSpecM1D`
- Issues: `https://github.com/adcm2/DSpecM1D/issues`

## Features

- Sparse frequency-domain SEM solver for radial, spheroidal, and toroidal modes
- CMT earthquake source representation with arbitrary receiver geometry
- Filtered time-domain seismogram generation from standard parameter files
- Self-contained unit/component tests plus protected paper-reproduction examples
- Static website and Doxygen API reference

## Dependencies

Required system packages for the first release:

- FFTW3
- OpenMP-capable compiler toolchain
- CMake 3.24 or newer

Remaining third-party C++ dependencies are fetched automatically from pinned
revisions during the first configure.

On Ubuntu or Debian:

```bash
sudo apt install build-essential cmake git libfftw3-dev
```

On macOS with Homebrew:

```bash
brew install cmake fftw
```

This first release does not require NetCDF, BLAS, or LAPACK.

## Building

The default single-config build now configures in `Release` mode when
`CMAKE_BUILD_TYPE` is unset, so a fresh clone builds in release mode by
default:

```bash
git clone https://github.com/adcm2/DSpecM1D.git
cd DSpecM1D
cmake -S . -B build
cmake --build build --parallel
```

If you prefer presets:

```bash
cmake --preset release
cmake --build --preset release
```

## Quick Start

The main user-facing executable for the first release is
`generate_seismogram`. It reads a standard DSpecM1D parameter file and writes a
filtered time-domain seismogram to the output prefix specified inside that
parameter file.

Using the bundled tutorial parameter file:

```bash
cmake -S . -B build
cmake --build build --target generate_seismogram
cd build
export OMP_NUM_THREADS=2
./bin/generate_seismogram data/params/t1.txt
```

That writes:

```text
build/output/t1.out_t.out
```

The output location comes from the parameter-file `output_prefix`. The CLI does
not take a separate output-path override for this first release.

The output file is a semicolon-delimited table with:

- column 1: time in seconds
- remaining columns: seismogram rows in solver order

For one receiver, those rows correspond to the three output components. For
multiple receivers, the component triplets are written in parameter-file order.

`t1` remains the self-contained tutorial and runtime smoke example, while
`ex1` through `ex7` remain protected paper-reproduction workflows.

## Parameter File Format

The parser reads an ordered sequence of value lines, not a keyed `name value`
format. Blank lines and full-line comments beginning with `#` are ignored, but
the remaining values must appear in the expected order.

See:

- [data/params/t1.txt](data/params/t1.txt)
- website parameter-file guide: `https://adcm2.github.io/DSpecM1D/parameter-files.html`

## Normalisation

DSpecM1D uses three reference scales:

- length: `1 km = 1000 m`
- density: `5515 kg m^-3`
- time: `1 / sqrt(pi G rho_ref)`

The corresponding helper is documented in `DSpecM1D/src/NormClass.h` and in the
API reference on the docs site.

## First Release Notes

This is the first public release. The README intentionally stays short and
focuses on getting the library compiled and running. For testing, debugging,
API details, and the paper-reproduction workflows, use the website.

Planned next steps after the first release:

- highest priority: Love numbers
- next: optional LAPACK banded solver support
- next: further cleanup and small API improvements

## Documentation

To build the static site locally:

```bash
cmake -S . -B build/docs -DDSPECM1D_BUILD_DOCS=ON
cmake --build build/docs --target website
cmake --build build/docs --target preview-website
```

The local preview server runs at `http://127.0.0.1:8000`.

The repository also includes `.github/workflows/pages.yml`, which publishes the
same site through GitHub Pages from `main`.

## Repository

- Repository: `https://github.com/adcm2/DSpecM1D`
- Issues: `https://github.com/adcm2/DSpecM1D/issues`
- Contact: use the issue tracker for questions, bug reports, and release feedback

## License

DSpecM1D is released under the GNU General Public License v3.0.
See [LICENSE](LICENSE).
