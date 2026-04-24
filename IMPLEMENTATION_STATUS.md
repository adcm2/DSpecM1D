# DSpecM1D Release Readiness Implementation Status

Last updated: 2026-04-24 (Europe/London)

This file is a durable handoff log for the release-readiness work currently in
progress. It exists so work can continue safely even if the session stops
unexpectedly.

Follow-on hardening work is tracked separately in `API_HARDENING_PLAN.md`.

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

Status: verified

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
- Added a local documentation preview target:
  - `preview-website`
- Added repository/support metadata to the README and website:
  - repository URL
  - issue tracker URL
  - contact via repository issues
- Added the project `LICENSE` file using the GNU General Public License v3.0

Verified on 2026-04-22:

- fresh docs configure succeeded in `build/phase6_verify`
- `cmake --build build/phase6_verify --target website` completed successfully
- the generated site includes the expected release-facing pages and API-doc
  integration path
- the `preview-website` target is available for local inspection

### Phase 7: Final Release Validation

Intent:

- Run the full release checklist on a clean build

Status: verified

Completed so far on 2026-04-22:

- Copied the canonical paper-comparison files into the repository under:
  - `data/reference/yspec/`
  - `data/reference/specnm/`
  - `data/reference/mineos/bolivia/`
  - `data/reference/mineos/noheader/`
- Added `data/reference/README.md` to document provenance and example usage
- Updated `benchmarks/ex1.cpp`, `benchmarks/ex2.cpp`, `benchmarks/ex3.cpp`,
  and `benchmarks/ex7.cpp` to read the in-repo comparison files via
  build-directory-relative paths
- Verified a fresh benchmark-only configure/build in `build/phase7_verify`
- Verified the copied reference files appear in
  `build/phase7_verify/data/reference/`
- Verified a clean release-validation configure/build in `build/phase7_release`
- Verified all 20 modular unit/component tests pass in `build/phase7_release`
- Verified the docs site builds successfully in `build/phase7_release`
- Tightened the `t1` smoke workflow so it is safe on a developer machine:
  - smoke test remains serial
  - OpenMP thread count is capped
  - timeout remains enforced
- Reworked `t1` to use a dedicated lighter input file:
  - added `data/params/t1.txt`
  - updated `tutorials/t1.cpp` to read `t1.txt` instead of `ex1.txt`
- Set single-config builds to default to `Release` when
  `CMAKE_BUILD_TYPE` is unset
- Verified the new default with a fresh tree:
  - `build/phase7_smoke_release` configured with
    `CMAKE_BUILD_TYPE:STRING=Release`
- Verified guarded smoke now passes in the fresh release-default tree:

```bash
ctest --test-dir build/phase7_smoke_release --output-on-failure -L smoke
```

  Results:
  - `smoke_build_benchmark_ex1` passed
  - `smoke_tutorial_t1` passed in about 15.5 seconds

Important current status:

- The main software-engineering release gates now pass locally:
  - clean configure/build
  - modular tests
  - guarded smoke
  - docs build
  - protected example rebuilds with in-repo reference data
- The final install/downstream package recheck from the latest tree now also
  passes:
  - clean install succeeded to `build/install_phase7_final`
  - fresh downstream `find_package(DSpecM1D)` configure succeeded against that
    install
  - downstream consumer built successfully in
    `work/phase5_consumer/build_phase7_final`
  - the downstream consumer executable ran with exit code `0`
- Paper-reproduction validation in the supported environment has now been
  confirmed by the user.

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
- The optional paper-validation path still needs an explicit supported-environment
  verification run

## Resume Rule

If work resumes later, continue in this order:

1. verify install/package flow
2. run optional paper-validation checks in the supported environment

Do not refactor `ex1`-`ex7` away from their comparison workflows. Those are
protected paper-reproduction examples.

## 2026-04-24 Dependency And In-Housing Update

This section is the durable handoff point for the follow-on dependency work.
It records both what was already changed and the agreed direction for bringing
 selected external libraries in house.

### Current dependency status

Verified on 2026-04-24:

- build/doc/test hygiene was tightened:
  - tracked generated root `CMakeFiles/` artifacts were removed
  - tracked Python `__pycache__` bytecode files were removed
  - `full_codebase.txt` was removed
  - `.gitignore` was expanded to ignore common generated CMake and Python files
- README and static docs were updated so:
  - Ninja is optional rather than required
  - NetCDF is explicitly documented as not required by DSpecM1D
  - manual runs document `OMP_NUM_THREADS`
  - smoke-test examples now use the project-specific `dspecm1d_smoke` label
- FetchContent build/subbuild state was moved from source-level `.deps/` into
  per-build-tree `_deps/` directories in:
  - `CMakeLists.txt`
  - `cmake/DSpecM1DConfig.cmake.in`
- This fixed the generator-mismatch failure caused by sharing one `.deps`
  cache between Ninja and Unix Makefiles build trees

Verified build/test results on 2026-04-24:

- a clean default-generator tree in `build/clean` configured and built
- `ctest --test-dir build/clean -N` showed exactly 46 project tests:
  - 44 unit/component
  - 2 smoke
  - no leaked dependency tests such as `fastmath`
- in the main `build/` tree:
  - `cmake --build build/` succeeded
  - `ctest --test-dir build --output-on-failure -L "unit|component"` passed
  - `ctest --test-dir build --output-on-failure -L dspecm1d_smoke` passed

### Public dependency interface correction

The public interface target was audited against the actual public headers.

Important finding:

- `GSHTrans` is a real public dependency because DSpecM1D public headers include
  `<GSHTrans/Core>`
- `GaussQuad` and `FFTWpp` are not directly included by DSpecM1D public headers
  and should not be named as public interface dependencies unless needed for a
  real downstream contract

Changes made and verified on 2026-04-24:

- `CMakeLists.txt` now links `DSpecM1D` publicly against:
  - `Eigen3::Eigen`
  - `Interpolation`
  - `GSHTrans`
  - `PlanetaryModel`
  - `EarthMesh`
  - `SpectraSolver`
  - `OpenMP::OpenMP_CXX`
- `GaussQuad` and `FFTWpp` were removed from the DSpecM1D public interface list
- `cmake/DSpecM1DConfig.cmake.in` was updated to mirror the same installed
  interface

Verification outcome:

- the trimmed public interface still builds and passes:
  - all 44 unit/component tests
  - both smoke tests

Implementation note:

- the previous build happened to work without explicitly linking `GSHTrans`
  because the fetched `GSHTrans` project currently leaks its include directory
  via a typo in its own CMake setup
- DSpecM1D should not rely on that accidental transitive include path

### Agreed in-housing order

The agreed order for bringing selected libraries in house is:

1. `PlanetaryModel`
2. `EarthMesh`
3. `SpectraSolver`

The following libraries should remain external:

- `GSHTrans`
- `Interpolation`
- `GaussQuad`
- `FFTWpp`

### Agreed migration approach

The agreed migration rule is:

- avoid breaking working builds
- first get imported code into a usable, stable in-repo state
- only then do cleanup, API tidying, or structural simplification

That means each library migration should be split into two phases:

1. Stabilization phase:
   - copy or vendor the code into the repository
   - preserve current behaviour and include paths as much as possible
   - keep tests and smoke workflows passing
   - prefer compatibility shims over eager refactoring
2. Cleanup phase:
   - reduce legacy wrappers
   - rationalize include paths and naming
   - remove redundant compatibility layers
   - tighten DSpecM1D-facing boundaries once behaviour is proven stable

### Next planned work for in-housing

Start with `PlanetaryModel`.

Initial tasks for the `PlanetaryModel` in-housing pass:

- audit exactly which headers/types DSpecM1D currently consumes from
  `PlanetaryModel`
- copy the library into an in-repo home without changing behaviour
- switch DSpecM1D includes/build logic to the in-repo copy
- keep benchmarks, tutorial, unit/component tests, and smoke tests passing
- defer cleanup/refactoring until the in-repo copy is verified stable

Do not start `EarthMesh` or `SpectraSolver` until `PlanetaryModel` is stable
and the resulting build/test state has been rechecked.

### `PlanetaryModel` consumption map

Audit completed on 2026-04-24.

Current include shape:

- DSpecM1D code only includes `PlanetaryModel` through `<PlanetaryModel/All>`
- Current include sites:
  - `DSpecM1D/src/InputParametersNew.h`
  - `DSpecM1D/src/FullSpec.h`
  - `tutorials/t1.cpp`
  - `benchmarks/ex1.cpp` through `benchmarks/ex7.cpp`
  - `tests/test_sem_component.cpp`

PlanetaryModel umbrella contents actually used by DSpecM1D:

- `PlanetaryModel/src/EARTHMODELS.h`
- `PlanetaryModel/src/Concepts.h`
- `PlanetaryModel/src/Concepts2.h`

Concrete `PlanetaryModel` class definitions currently used directly:

- `EarthModels::ModelInput<FLOAT, INTEGRAL>`

Important correction:

- `prem_norm<double>` is not provided by `PlanetaryModel`
- DSpecM1D defines it locally in `DSpecM1D/src/NormClass.h`
- the current `PlanetaryModel` loading path passes that DSpecM1D-owned
  normalization class into `EarthModels::ModelInput(...)`

`InputParametersNew` currently depends on:

- `decltype(EarthModels::ModelInput(path, prem_norm<double>{}, "true"))`
  as its concrete stored `ModelType`
- `EarthModels::ModelInput(...)` for earth-model loading
- DSpecM1D-owned `prem_norm<double>` for normalization during model loading
- returned model methods:
  - `TimeNorm()`
  - `TREF()`
  - `LengthNorm()`

Template/generic model interface DSpecM1D relies on downstream:

- normalization:
  - `LengthNorm()`
  - `MassNorm()`
  - `TimeNorm()`
  - `DensityNorm()`
- geometry:
  - `NumberOfLayers()`
  - `LowerRadius(i)`
  - `UpperRadius(i)`
  - `OuterRadius()`
- material / state:
  - `Density(i)(r)`
  - `IsSolid(i)`
  - `IsFluid(i)`
- wave speeds / elastic quantities:
  - `VP(i)(r)`
  - `VPV(i)(r)`
  - `VPH(i)(r)`
  - `VS(i)(r)`
  - `VSV(i)(r)`
  - `VSH(i)(r)`
  - `A(i)(r)`
  - `C(i)(r)`
  - `F(i)(r)`
  - `L(i)(r)`
  - `N(i)(r)`
  - `Kappa(i)(r)`
  - `Mu(i)(r)`
- metadata:
  - `TREF()`

Important scoping note:

- DSpecM1D does not currently appear to use the wider example/demo model zoo in
  `PlanetaryModel`
- the first in-house pass should focus on preserving:
  - `EarthModels::ModelInput`
  - only the supporting code inside `EARTHMODELS.h` that `ModelInput` itself
    needs
- classes such as `PREM`, `PERTPREM`, `HOMOSPHERE`, `HOMOBOUND*`, and
  `TomographyZeroModel` are not directly referenced by DSpecM1D code today
- `PlanetaryModel/src/Concepts.h` and `PlanetaryModel/src/Concepts2.h` are not
  currently required for the immediate migration, but should be retained in the
  repository for now and reconsidered during a later API-hardening pass
- do not try to tidy or reduce the imported model code during the first copy
  step; preserve behavior first

### `PlanetaryModel` in-housing progress

Progress completed on 2026-04-24.

What was moved in-repo:

- a local `ModelInput` implementation now lives at:
  - `DSpecM1D/src/model_info/ModelInput.h`
- a DSpecM1D-facing umbrella header was added at:
  - `DSpecM1D/ModelInput`
- the DSpecM1D public umbrella now exposes that header through:
  - `DSpecM1D/All`

Compatibility layer added to avoid breaking `EarthMesh`:

- a local compatibility tree now exists at:
  - `PlanetaryModel/All`
  - `PlanetaryModel/src/EARTHMODELS.h`
  - `PlanetaryModel/src/Concepts.h`
  - `PlanetaryModel/src/Concepts2.h`
  - `PlanetaryModel/CMakeLists.txt`
- this keeps the old `<PlanetaryModel/...>` include shape valid while routing
  `EARTHMODELS.h` to the vendored `DSpecM1D/src/model_info/ModelInput.h`
- the compatibility `CMakeLists.txt` provides an interface target named
  `PlanetaryModel` so fetched dependencies that still expect that target
  continue to configure

Code switched to local headers:

- direct DSpecM1D/example/test includes of `<PlanetaryModel/All>` were changed
  to `<DSpecM1D/ModelInput>` in:
  - `DSpecM1D/src/InputParametersNew.h`
  - `DSpecM1D/src/FullSpec.h`
  - `tutorials/t1.cpp`
  - `benchmarks/ex1.cpp` through `benchmarks/ex7.cpp`
  - `tests/test_sem_component.cpp`

Build-system changes:

- top-level `CMakeLists.txt` no longer declares `PlanetaryModel` as a direct
  `FetchContent` dependency of DSpecM1D
- `DSpecM1D` no longer links `PlanetaryModel` on its public interface
- `cmake/DSpecM1DConfig.cmake.in` was updated to mirror that installed
  interface change
- the install now includes the local `PlanetaryModel/` compatibility tree
- `FETCHCONTENT_SOURCE_DIR_PLANETARYMODEL` is set to the vendored local
  compatibility directory in both the top-level build and installed package
  config so `EarthMesh` resolves `PlanetaryModel` locally instead of trying to
  fetch the external repository

Important compatibility finding:

- the first build after moving `ModelInput` failed because `EarthMesh` still
  included and fetched external `PlanetaryModel`, causing duplicate definitions
  of `EarthModels::ModelInput`
- the local `PlanetaryModel/` compatibility shim resolved that conflict
- `EarthMesh` also tries to fetch `PlanetaryModel` from an SSH Git URL; the
  `FETCHCONTENT_SOURCE_DIR_PLANETARYMODEL` override is what prevents smoke-test
  rebuilds from attempting that network access

Verification outcome after the vendoring pass:

- `cmake --build build/` succeeded
- `ctest --test-dir build --output-on-failure -L "unit|component"` passed
  `44/44`
- `ctest --test-dir build --output-on-failure -L dspecm1d_smoke` passed `2/2`

Regression found and corrected on 2026-04-24:

- the first vendored `ModelInput` implementation changed the behavior of the
  three-argument constructor
- in the original external `PlanetaryModel`, the constructor
  `ModelInput(path, norm, bool)` does not use that boolean to force isotropic
  column parsing; it still follows the file header `ifanis`
- the vendored copy incorrectly used that argument to switch anisotropic files
  onto the six-column isotropic read path
- because DSpecM1D call sites commonly pass `"true"` as the third argument,
  that bug flattened `VPH`, `VSH`, and `eta` for anisotropic PREM inputs and
  changed seismogram output even though builds and existing tests still passed
- the constructor semantics were restored to match the original implementation:
  - one- and two-argument constructors keep the original parsing behavior
  - the three-argument constructor now again follows `ifanis` from the file and
    only preserves the original informational print behavior

Regression coverage added:

- `tests/test_input_parameters_new.cpp` now includes a check that the current
  DSpecM1D loading path preserves anisotropic columns from
  `data/models/prem.200.noatten.txt`
- this test specifically guards against silently collapsing an anisotropic model
  into isotropic `VPH/VSH/eta` values during future in-housing edits

Current state:

- `PlanetaryModel` is now effectively in-house for DSpecM1D through the local
  vendored `ModelInput` and compatibility layer
- the migration is in the intended stabilization phase, not cleanup phase
- some docs may still mention `PlanetaryModel` as an external dependency and
  should be cleaned up later

Next recommended step:

- begin the same stabilization-first audit for `EarthMesh`

### `EarthMesh` in-housing progress

Progress completed on 2026-04-24.

Scope and approach:

- the full `EarthMesh` header implementation was copied into the repository
  under:
  - `DSpecM1D/src/model_info/earthmesh/RadialMesh.h`
  - `DSpecM1D/src/model_info/earthmesh/RadialMeshDeclaration.h`
  - `DSpecM1D/src/model_info/earthmesh/RadialMeshDefinition.h`
- compatibility wrappers were added at:
  - `EarthMesh/All`
  - `EarthMesh/Declaration`
  - `EarthMesh/src/RadialMesh.h`
  - `EarthMesh/src/RadialMeshDeclaration.h`
  - `EarthMesh/src/RadialMeshDefinition.h`
  - `EarthMesh/CMakeLists.txt`
- as with `PlanetaryModel`, the compatibility tree preserves the existing
  `<EarthMesh/...>` include shape while routing the owned implementation to
  `DSpecM1D/src/model_info/earthmesh`

Build-system changes:

- top-level CMake now sets `FETCHCONTENT_SOURCE_DIR_EARTHMESH` to the vendored
  local `EarthMesh/` directory
- the installed package config mirrors that setting so downstream builds also
  resolve `EarthMesh` locally
- the project install now includes the vendored `EarthMesh/` compatibility tree
- the top-level fetch still creates the `EarthMesh` target, but it now resolves
  to the local vendored source rather than the external repository

Compatibility findings during the move:

- the first wrapper version reused the same include guards as the copied
  `EarthMesh` headers, which suppressed the wrapped definitions and caused
  `EarthMesh::RadialMesh` to disappear from compile units
- the wrapper guards were changed to unique compatibility-only names
- the first EarthMesh-local CMake pass also exposed that `EarthMesh` expects a
  real `PlanetaryModel` target name during CMake configure/link propagation
- restoring the top-level `PlanetaryModel` fetch declaration with the vendored
  local source override fixed that target expectation without reintroducing the
  external dependency

Verification outcome after the EarthMesh vendoring pass:

- `cmake --build build/` succeeded
- `ctest --test-dir build --output-on-failure -R "InputParametersNewTests|SEMComponentTests|PreferredSolverApiTests"`
  passed `13/13`
- `ctest --test-dir build --output-on-failure -L dspecm1d_smoke` passed `2/2`

Current state after EarthMesh move:

- `EarthMesh` is now effectively in-house via the vendored implementation in
  `DSpecM1D/src/model_info/earthmesh` plus the local compatibility tree
- both `PlanetaryModel` and `EarthMesh` are now in the stabilization phase of
  migration, with build and smoke verification passing

Immediate next steps agreed in discussion:

- add a temporary migration-phase seismogram reference check around `ex2`
  (useful during in-housing work, not intended as a final release test)
- then do a documentation cleanup pass for `PlanetaryModel` and `EarthMesh`
- after that, review whether DSpecM1D-owned earth/model-on-mesh code should be
  reorganized into the same `model_info` area

### Post-migration cleanup completed

Cleanup completed on 2026-04-24 after the initial vendoring passes.

Final architecture after cleanup:

- `ModelInput` now lives directly in:
  - `DSpecM1D/src/model_info/ModelInput.h`
- `EarthMesh` now lives directly in:
  - `DSpecM1D/src/model_info/earthmesh/RadialMesh.h`
  - `DSpecM1D/src/model_info/earthmesh/RadialMeshDeclaration.h`
  - `DSpecM1D/src/model_info/earthmesh/RadialMeshDefinition.h`
- local model concepts now live in:
  - `DSpecM1D/src/model_info/ModelConcepts.h`
- public DSpecM1D entry points were added at:
  - `DSpecM1D/ModelInput`
  - `DSpecM1D/EarthMesh`
  - `DSpecM1D/ModelConcepts`

Compatibility-layer cleanup performed:

- the temporary `PlanetaryModel/` and `EarthMesh/` compatibility trees were
  removed after internal code was switched over to DSpecM1D-owned headers
- internal includes now use:
  - `<DSpecM1D/ModelInput>`
  - `<DSpecM1D/EarthMesh>`
  - `<DSpecM1D/ModelConcepts>`
- there are no remaining active library/test/example includes of
  `<PlanetaryModel/...>` or `<EarthMesh/...>`

`ModelInput` constructor cleanup:

- the temporary three-argument constructor path was removed
- the model loader now exposes only the direct constructor forms without the
  legacy boolean override
- model loading now always follows the file metadata to determine isotropic vs
  anisotropic parsing
- DSpecM1D, tutorials, benchmarks, and tests were updated to stop passing the
  old `"true"` argument

Build-system cleanup:

- `PlanetaryModel` and `EarthMesh` are no longer fetched as external projects
  for DSpecM1D
- `GaussQuad` is now linked explicitly on the DSpecM1D public interface because
  the vendored EarthMesh headers depend on it directly
- the installed package config was updated to mirror the cleaned dependency
  surface

Temporary migration validation added:

- a migration-only `ex2` reference check is now available behind
  `DSPECM1D_ENABLE_MIGRATION_TESTS`
- this check runs `ex2`, writes outputs into the build tree, and compares the
  DSpecM1D `Z/N/E` traces against the checked-in `plotting/outputs/ex2_t.out`
  using numeric tolerances rather than exact file identity
- this is intentionally a migration-phase check, not a final release test

Verification after cleanup:

- `cmake --build build/` succeeded
- `ctest --test-dir build --output-on-failure -L "unit|component"` passed
  `45/45`
- `ctest --test-dir build --output-on-failure -L dspecm1d_smoke` passed `2/2`
- `ctest --test-dir build --output-on-failure -L migration` passed `1/1`

Current state after cleanup:

- `PlanetaryModel` and `EarthMesh` are now fully in-house in the active library
  layout, not just behind compatibility wrappers
- the library is back in a clean, verified state before documentation updates
- the next step should now be the documentation pass for the updated dependency
  story and model-info layout

Additional model-info cleanup:

- `MeshModel.h` was moved into `DSpecM1D/src/model_info/` so the mesh-derived
  model helpers now live alongside `ModelInput`, `ModelConcepts`, and the
  vendored `earthmesh` headers
- the temporary forwarding header at `DSpecM1D/src/MeshModel.h` has now been
  removed, so the tree uses only the `model_info/MeshModel.h` path

### SpectraSolver audit

Current DSpecM1D usage of `SpectraSolver` is concentrated in the frequency-axis
helper and FFT/post-processing utilities rather than its internal matrix solver
machinery.

Direct in-tree include pattern:

- all active code included `SpectraSolver` through `<SpectraSolver/FF>` at the
  time of the audit
- no active DSpecM1D file includes `SpectraSolver/ODES`, `SpectraSolver/All`,
  or the old-code headers

Confirmed `FreqFull` usage surface in DSpecM1D:

- construction via the full fourteen-argument constructor in the main library
  and examples
- construction via the short nine-argument constructor in unit tests
- methods used in active code:
  - `w()`
  - `w(int)`
  - `dt()`
  - `df()`
  - `ep()`
  - `f11()`
  - `f22()`
  - `i1()`
  - `i2()`
  - `nt()`
  - `nt0()`
  - `t2()`

Confirmed post-processing usage surface in DSpecM1D:

- `processfunctions::freq2time`
- `processfunctions::filtfreq2time`
- `processfunctions::fulltime2freq`

Not currently used by active DSpecM1D code:

- `FreqFull::Spectra_Raw*` solver member functions
- `filters::` or `filterclass::` APIs directly
- `MatrixReplaceFT` directly
- `BlockMatFreePreconditioner` directly

Implication for an in-house move:

- the minimal practical subset is much smaller than the full upstream
  repository
- the active dependency appears to be:
  - `SpectraSolver/FF`
  - `FrequencyFull.h`
  - `FrequencyFullConstructor.h`
  - `FrequencyFullDefinitions.h`
  - `postprocessfunctions.h`
  - `filter_header.h`
- however, if we copy those files verbatim, `FrequencyFullDefinitions.h` still
  pulls in `BlockPreconditioner.h` and `matrix_rcopy.h` because upstream keeps
  the unused `Spectra_Raw*` methods inline in the same header
- that means the cleanest in-house migration is likely not a blind copy of the
  upstream umbrella header, but a DSpecM1D-owned trimmed frequency/FFT layer
  that preserves `FreqFull` behaviour and the three `processfunctions::*`
  helpers we actually use

### `SpectraSolver` in-housing progress

Completed on 2026-04-24.

What was moved in-repo:

- a local public frequency/FFT umbrella now exists at:
  - `DSpecM1D/FrequencyTools`
- the owned implementation now lives in:
  - `DSpecM1D/src/frequency_info/FreqFull.h`
  - `DSpecM1D/src/frequency_info/PostprocessFunctions.h`

What was preserved intentionally:

- the namespace `SpectraSolver`
- the class name `SpectraSolver::FreqFull`
- the namespace `processfunctions`
- the active DSpecM1D-facing behavior of:
  - `FreqFull` construction
  - `FreqFull` accessors used in the library/examples/tests
  - `processfunctions::freq2time`
  - `processfunctions::filtfreq2time`
  - `processfunctions::fulltime2freq`

What was intentionally not imported:

- the inline `FreqFull::Spectra_Raw*` solver methods
- `MatrixReplaceFT`
- `BlockMatFreePreconditioner`
- the upstream `OLD_CODE` and `ODES` surfaces

Codebase changes:

- active includes of `<SpectraSolver/FF>` were replaced with
  `<DSpecM1D/FrequencyTools>` across the library, benchmarks, tutorial, and
  tests
- top-level `CMakeLists.txt` no longer fetches `SpectraSolver`
- `DSpecM1D` no longer links `SpectraSolver` on its public interface
- `cmake/DSpecM1DConfig.cmake.in` was updated so installed consumers no longer
  fetch or link the external `SpectraSolver` target

Current state:

- `SpectraSolver` is no longer an external build dependency of DSpecM1D
- the active library now owns the frequency-axis and FFT/post-processing layer
- the retained external dependencies remain:
  - `Eigen3`
  - `FFTWpp`
  - `GaussQuad`
  - `GSHTrans`
  - `Interpolation`

Verification after in-housing:

- `cmake --build build/` succeeded
- focused tests covering the moved functionality passed:
  - `SignalFilteringTests`
  - `OutputWriterTests`
  - `InputParametersNewTests`
  - `SEMComponentTests`
  - `PreferredSolverApiTests`
- `ctest --test-dir build --output-on-failure -L dspecm1d_smoke` passed `2/2`
- `ctest --test-dir build --output-on-failure -L migration` passed `1/1`
