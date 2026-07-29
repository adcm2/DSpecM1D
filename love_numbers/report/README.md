# Reproducing the Love-number technical report

The editable report is `love_numbers_report.tex`; the compiled deliverable is
`love_numbers_report.pdf`. The workflow requires a C++23 compiler, CMake,
Python 3 with Matplotlib, and a standard LaTeX installation. Paper validation
also requires Fortran, BLAS/LAPACK, the pinned gia3D sources, and the official
ELLN supporting-information archive.

## 1. Configure and test

From the repository root:

```bash
ROOT=$(pwd)
NORMAL_BUILD=/tmp/dspecm1d-stage20-normal
PAPER_BUILD=/tmp/dspecm1d-stage20-paper
ELLN_ARCHIVE=/tmp/ggx444_supp.zip

cmake -S "$ROOT" -B "$NORMAL_BUILD" \
  -DDSPECM1D_BUILD_LOVE_NUMBERS=ON \
  -DDSPECM1D_BUILD_TESTS=ON \
  -DDSPECM1D_ENABLE_PAPER_VALIDATION=OFF
cmake --build "$NORMAL_BUILD" --parallel 4
ctest --test-dir "$NORMAL_BUILD" --output-on-failure

cmake -S "$ROOT" -B "$PAPER_BUILD" \
  -DDSPECM1D_BUILD_LOVE_NUMBERS=ON \
  -DDSPECM1D_BUILD_TESTS=ON \
  -DDSPECM1D_ENABLE_PAPER_VALIDATION=ON \
  -DDSPECM1D_ELLN_ARCHIVE="$ELLN_ARCHIVE"
cmake --build "$PAPER_BUILD" --parallel 4
ctest --test-dir "$PAPER_BUILD" --output-on-failure
ctest --test-dir "$PAPER_BUILD" --show-only=json-v1
```

`DSPECM1D_GIA3D_SOURCE_DIR` and `DSPECM1D_GIA3D_CORE_SOURCE_DIR` may name
local checkouts at the pinned revisions recorded in the report. Existing
FetchContent source overrides may likewise be supplied for an offline build.

## 2. Benchmark datasets

The following single command reproduces every CSV in `data/` by running the
existing validation programs and parsing their labelled output:

```bash
python3 "$ROOT/love_numbers/report/collect_benchmarks.py" \
  --source-root "$ROOT" \
  --build-dir "$PAPER_BUILD" \
  --output-dir "$ROOT/love_numbers/report/data" \
  --elln-archive "$ELLN_ARCHIVE"
```

An extracted official package can be used with `--elln-source-dir` instead.
The collector performs no Love-number or gravity calculation itself. It runs
these validation paths:

| Dataset | Reproducing validation command |
| --- | --- |
| `background_gravity_sfs.csv` | `compare_controlled_sphere.py ... solid_fluid_solid.dspec solid_fluid_solid.gia3d solid-fluid-solid` |
| `controlled_sfs_convergence.csv` | The same controlled solid--fluid--solid command, including `N=8,16,32,64,128,256` |
| `controlled_finest_summary.csv` | `compare_controlled_sphere.py` separately for the homogeneous, central-fluid, and solid--fluid--solid paired decks |
| `isotropic_prem_knot_convergence.csv` | `compare_isotropic_prem.py $PAPER_BUILD/bin/dspecm1d_controlled_love_numbers $PAPER_BUILD/bin/dspecm1d_gia3d_controlled_reference $PAPER_BUILD/bin/dspecm1d_gia3d_export_isotropic_prem $PAPER_BUILD/love_numbers/validation/isotropic-prem-models` |
| `isotropic_prem_finest.csv` | The same isotropized-PREM command, filtered to 12.5 km knots and radial step 0.0025 |
| `isotropic_prem_diagnostics.csv` | The solve diagnostics from the same finest isotropized-PREM runs |
| `elln_gravitational_constant.csv` | `compare_gravitational_constants.py` using the production and validation-only ELLN-`G` executables, 6.25 km TI deck, manifest, Table S7, and official archive |
| `elln_sensitivity_summary.csv` | The sensitivity rows emitted by that gravitational-constant comparison |

Each CSV begins with comment lines recording its model, quantity definition,
resolution, and generating validation script. The labelled column row follows
those comments.

## 3. Figures

Generate the four vector PDF figures from the CSV files:

```bash
MPLCONFIGDIR=/tmp/dspecm1d-stage20-matplotlib \
  python3 "$ROOT/love_numbers/report/generate_figures.py"
```

The script reads the stored datasets only. It does not call production code or
reimplement any scientific calculation.

## 4. Report

Build the appendix from `../TESTS.md` and compile the report:

```bash
"$ROOT/love_numbers/report/build_report.sh"
```

The script prefers:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error \
  love_numbers_report.tex
```

If `latexmk` is unavailable it uses repeated `pdflatex`, with `bibtex` when
available. No shell escape or unusual LaTeX package is required.

For final presentation checks:

```bash
pdfinfo "$ROOT/love_numbers/report/love_numbers_report.pdf"
pdftotext "$ROOT/love_numbers/report/love_numbers_report.pdf" -
pdftoppm -png -r 120 \
  "$ROOT/love_numbers/report/love_numbers_report.pdf" \
  /tmp/dspecm1d-love-report-page
```
