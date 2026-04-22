# DSpecM1D Release Readiness Implementation Status

Last updated: 2026-04-22 (Europe/London)

This file is a durable handoff log for the release-readiness work currently in
progress. It exists so work can continue safely even if the session stops
unexpectedly.

## Goal

Implement the release-readiness plan in safe phases while preserving the
current behavior of:

- `benchmarks/ex1.cpp` through `benchmarks/ex7.cpp` as paper-reproduction
  examples
- `tutorials/t1.cpp` as the self-contained tutorial/smoke example

The main regression safety net should be small, modular, self-contained tests,
not full seismogram comparisons.

## Phase Plan

### Phase 1: Lock The Baseline

Intent:

- Preserve `ex1`-`ex7` as paper-reproduction workflows
- Preserve `t1` as the self-contained tutorial/smoke example
- Document external dependencies of the paper examples

Status: partially implemented

Completed:

- README now distinguishes `t1` from `ex1`-`ex7`
- Static docs site includes a dedicated paper-reproduction page
- `benchmarks/CMakeLists.txt` now defines a `paper_examples` target
- Optional paper-validation CTest registration was added behind
  `DSPECM1D_ENABLE_PAPER_VALIDATION`

Still to verify:

- Rebuild `ex1`-`ex7` and `t1` after all current edits
- Confirm no workflow drift in the paper examples

### Phase 2: Build A Real Small-Test Suite

Intent:

- Replace placeholder tests with modular, self-contained unit/component tests
- Keep tests synthetic and local

Status: verified

Completed:

- Deleted placeholder `tests/test_math.cpp`
- Added test configuration helper:
  - `tests/config.h.in`
  - `tests/test_utils.h`
- Added modular test files:
  - `tests/test_input_parser.cpp`
  - `tests/test_input_parameters_new.cpp`
  - `tests/test_output_writers.cpp`
  - `tests/test_reference_series_io.cpp`
  - `tests/test_signal_filtering.cpp`
  - `tests/test_spec_helpers.cpp`
  - `tests/test_sem_component.cpp`
- Updated `tests/CMakeLists.txt` to build the new suite and register tests
  with `unit;component` labels

Verified on 2026-04-22:

- the modular test suite compiles cleanly
- all 20 unit/component tests pass with:

```bash
ctest --test-dir build/dev --output-on-failure -L "unit|component"
```

Implementation note:

- the original single `dspecm1d_tests` executable caused ODR/linker clashes
  because several dependency headers pull in non-`inline` definitions
- the suite was restructured into multiple narrow test executables:
  - `dspecm1d_parser_tests`
  - `dspecm1d_input_context_tests`
  - `dspecm1d_output_writer_tests`
  - `dspecm1d_signal_filter_tests`
  - `dspecm1d_sem_component_tests`

### Phase 3: Strengthen Smoke And Validation Workflows

Intent:

- Keep `t1` as the mandatory runtime smoke test
- Distinguish unit/component tests, smoke, and paper validation

Status: partially verified

Completed:

- `tutorials/CMakeLists.txt` labels `smoke_tutorial_t1` with `smoke;tutorial`
- `benchmarks/CMakeLists.txt` labels optional paper validation with
  `validation;paper`
- `.github/workflows/ci.yml` now runs:
  1. configure
  2. build
  3. unit/component tests
  4. smoke tests
  5. docs site build

Still to verify:

- CI logic remains valid after local build/test verification
- optional paper-validation registration still needs an explicit run in the
  supported reproduction environment

Verified on 2026-04-22:

- smoke tests pass with:

```bash
ctest --test-dir build/dev --output-on-failure -L smoke
```

- `smoke_build_benchmark_ex1` passed
- `smoke_tutorial_t1` passed

### Phase 4: Public API Documentation And Safe Cleanup

Intent:

- Add Doxygen comments to release-facing headers
- Clarify preferred vs legacy APIs

Status: verified

Completed:

- Added/expanded Doxygen-style comments in:
  - `DSpecM1D/src/NormClass.h`
  - `DSpecM1D/src/InputParser.h`
  - `DSpecM1D/src/InputParametersNew.h`
  - `DSpecM1D/src/SignalFiltering.h`
  - `DSpecM1D/src/OutputWriters.h`
  - `DSpecM1D/src/FullSpec.h`
  - `DSpecM1D/src/SEM/SEM.h`
- Documented the preferred release-facing API in README and site content:
  - `InputParametersNew`
  - `SPARSESPEC::SparseFSpec`
  - `DSpecM::FilterOptions`
  - `DSpecM::applyFilter`
  - output helpers
- Added documentation or compatibility notes for:
  - `DSpecM1D/src/ReferenceSeriesIO.h`
  - `DSpecM1D/src/SpecHelpers.h`
  - `DSpecM1D/src/SpectraRunContext.h`
  - `DSpecM1D/src/SpectraMaster.h`
- Cleaned one harmless implementation signature mismatch in
  `DSpecM1D/src/SEM/SEMConstructor.h` to remove a Doxygen warning
- Adjusted Doxygen configuration to degrade cleanly when Graphviz `dot` is not
  available and to exclude noisy external-style docs from `BiCGSTABT.h`

Verified on 2026-04-22:

- fresh configure/build succeeded
- all 20 modular unit/component tests still passed
- `cmake --build build/dev --target website` completed successfully
- Doxygen now builds cleanly enough for the current release-facing scope

### Phase 5: Packaging And Reproducible Release Setup

Intent:

- Add install/export rules
- Generate a usable package config
- Pin dependency revisions

Status: verified

Completed:

- Top-level `CMakeLists.txt` now:
  - includes `GNUInstallDirs` and `CMakePackageConfigHelpers`
  - installs `DSpecM1D`
  - configures `DSpecM1DConfig.cmake`
  - installs `DSpecM1DConfig.cmake` and version file
  - installs `FindFFTW.cmake`
- Added `cmake/DSpecM1DConfig.cmake.in`
- Pinned FetchContent dependencies in `CMakeLists.txt` to fixed revisions
- Replaced the old global `macro(install)` workaround with a safer
  `CMAKE_SKIP_INSTALL_RULES`-based helper while fetching dependencies
- Cleaned test dependency packaging so project install no longer pulls
  `googletest` artifacts into the install prefix
- Added a minimal downstream package-consumer project in:
  - `work/phase5_consumer/CMakeLists.txt`
  - `work/phase5_consumer/main.cpp`
- Hardened `tests/test_utils.h` temporary-directory creation so parallel test
  runs do not collide on `std::rand()`-derived temp paths

Verified on 2026-04-22:

- fresh configure/build succeeded in `build/phase5_verify`
- modular test executables still passed after packaging changes
- clean install succeeded to `build/install_phase5_verify`
- installed package contents are clean:
  - headers under `include/DSpecM1D/`
  - package files under `lib/cmake/DSpecM1D/`
  - no stray `googletest` install artifacts
- fresh downstream `find_package(DSpecM1D)` configure succeeded against the
  installed package
- downstream consumer compiled and linked successfully, and the resulting
  executable ran with exit code `0`
- local repo hygiene was tightened after verification:
  - `.gitignore` now ignores `work/`
  - `.gitignore` now ignores `compile_commands.json`
  - `.gitignore` now ignores `CMakeUserPresets.json`

Implementation note:

- The installed package config recreates pinned dependency targets with
  `FetchContent` for downstream consumers. This keeps the header-only package
  usable without exporting third-party targets directly, but means a first
  downstream configure may require network access unless those dependencies are
  already cached locally.

### Phase 6: Website And Release-Facing Docs

Intent:

- Build a full static documentation website
- Fix README and release-facing docs

Status: implemented in code, not fully verified yet

Completed:

- Added docs build support:
  - `docs/CMakeLists.txt`
  - `docs/Doxyfile.in`
- Added static website files:
  - `docs/site/index.html`
  - `docs/site/installation.html`
  - `docs/site/parameter-files.html`
  - `docs/site/paper-reproduction.html`
  - `docs/site/library-testing.html`
  - `docs/site/release-notes.html`
  - `docs/site/assets/style.css`
- README now:
  - documents the real ordered parameter-file format
  - distinguishes tests vs paper examples
  - documents docs build target

Still to do or verify:

- Build `website` target
- If Doxygen is available, confirm API docs copy into `site/api/`
- License/contact are still placeholders and need real project decisions

### Phase 7: Final Release Validation

Intent:

- Run the full release checklist on a clean build

Status: not done yet

Needs:

- clean configure/build
- unit/component tests
- smoke test
- paper-reproduction validation in supported environment
- install/package verification
- docs/site verification

## Files Changed So Far

### Build and packaging

- `CMakeLists.txt`
- `cmake/DSpecM1DConfig.cmake.in`

### Docs and website

- `docs/CMakeLists.txt`
- `docs/Doxyfile.in`
- `docs/site/assets/style.css`
- `docs/site/index.html`
- `docs/site/installation.html`
- `docs/site/parameter-files.html`
- `docs/site/paper-reproduction.html`
- `docs/site/library-testing.html`
- `docs/site/release-notes.html`
- `README.md`

### CI and workflow classification

- `.github/workflows/ci.yml`
- `benchmarks/CMakeLists.txt`
- `tutorials/CMakeLists.txt`

### Public header docs / release-facing cleanup

- `DSpecM1D/src/NormClass.h`
- `DSpecM1D/src/InputParser.h`
- `DSpecM1D/src/InputParametersNew.h`
- `DSpecM1D/src/SignalFiltering.h`
- `DSpecM1D/src/OutputWriters.h`
- `DSpecM1D/src/FullSpec.h`
- `DSpecM1D/src/SEM/SEM.h`
- `DSpecM1D/src/ReferenceSeriesIO.h`
- `DSpecM1D/src/SpecHelpers.h`
- `DSpecM1D/src/SpectraRunContext.h`
- `DSpecM1D/src/SpectraMaster.h`
- `DSpecM1D/src/SEM/SEMConstructor.h`

### Tests

- `tests/CMakeLists.txt`
- `tests/config.h.in`
- `tests/test_utils.h`
- `tests/test_input_parser.cpp`
- `tests/test_input_parameters_new.cpp`
- `tests/test_output_writers.cpp`
- `tests/test_reference_series_io.cpp`
- `tests/test_signal_filtering.cpp`
- `tests/test_spec_helpers.cpp`
- `tests/test_sem_component.cpp`
- removed `tests/test_math.cpp`

## Last Known Verification State

The following was verified before the session interruption:

- CMake reconfigure succeeded for:

```bash
cmake -S . -B build/dev -G Ninja \
  -DDSPECM1D_BUILD_TESTS=ON \
  -DDSPECM1D_BUILD_BENCHMARKS=ON \
  -DDSPECM1D_BUILD_TUTORIALS=ON \
  -DDSPECM1D_ENABLE_SMOKE_TESTS=ON \
  -DDSPECM1D_BUILD_DOCS=ON
```

- An earlier build attempt exposed a real header dependency bug:
  `DSpecM1D/src/InputParametersNew.h` was missing a direct include of
  `NormClass.h`
- That issue was fixed by adding:

```cpp
#include "NormClass.h"
```

- The test helper and test files were also adjusted after the first compile
  pass:
  - `tests/test_utils.h` now includes `<cstdlib>`
  - `tests/test_output_writers.cpp` now includes `NormClass.h`
  - `tests/test_signal_filtering.cpp` now includes `NormClass.h`

After these fixes, the following has now been verified:

- `cmake --build build/dev --parallel` succeeds
- all Phase 2 unit/component tests pass
- smoke tests pass, including `t1`

## Next Commands To Run

Run these from the repository root:

```bash
cmake --install build/dev --prefix build/install
```

If packaging verification is reached, then also run:

```bash
cmake --install build/dev --prefix build/install
```

Then create a tiny downstream consumer project and test:

```cmake
find_package(DSpecM1D REQUIRED)
add_executable(consumer main.cpp)
target_link_libraries(consumer PRIVATE DSpecM1D::DSpecM1D)
```

## Known Open Risks

- `cmake/DSpecM1DConfig.cmake.in` is implemented but not yet validated in a
  downstream consumer
- docs/site build is implemented but not yet verified locally
- downstream package consumption is not yet verified
- License and Contact sections still need real project decisions and should not
  be considered complete
- The optional paper-validation path still needs an explicit supported-environment
  verification run

## Resume Rule

If work resumes later, continue in this order:

1. verify install/package flow
2. run optional paper-validation checks in the supported environment

Do not refactor `ex1`-`ex7` away from their comparison workflows. Those are
protected paper-reproduction examples.
