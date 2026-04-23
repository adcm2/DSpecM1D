# DSpecM1D API Hardening Status

Last updated: 2026-04-23 (Europe/London)

This file records progress against `API_HARDENING_PLAN.md` so the hardening
work can resume safely after interruptions.

## Goal

Strengthen modular API coverage and simplify the paper-comparison examples
without changing their intended scientific behaviour.

## Phase Status

### Phase A: Example Refactor Audit

Status: completed

Intent:

- Identify repeated logic in `benchmarks/ex1.cpp` through `benchmarks/ex7.cpp`
- Identify patterns that already exist as library functionality so we avoid
  duplicate definitions and repeated code
- Decide which logic should be owned by the library, which should be shared in
  benchmark-local helpers, and which should remain example-specific

#### High-level findings

- `ex1` is already the clearest example of the preferred public path:
  `InputParametersNew` + `SparseFSpec::spectra(InputParametersNew&)` +
  `FilterOptions`/`applyFilter` + reference loaders + output writers.
- `ex2`, `ex3`, and `ex7` still duplicate a large amount of logic that already
  exists in the library:
  - manual YSpec loading
  - manual Mineos loading
  - manual normalization selection
  - manual frequency/time comparison writing
- `ex4` shares some of the same setup duplication, but its record-section output
  is specific enough that it likely needs either a benchmark-local helper or a
  new dedicated writer helper.
- `ex5` and `ex6` are more specialised numerical examples. They should be
  cleaned up where obvious duplication exists, but they should remain more
  explicit than the other paper examples so their purpose stays readable.

#### Pattern mapping

1. Full workflow context setup

- Repeated in:
  - `ex2`, `ex3`, `ex4`, `ex5`, `ex7`
- Existing library owner:
  - `InputParametersNew`
- Notes:
  - `InputParametersNew` already bundles:
    - `InputParameters`
    - Earth model loading
    - `SourceInfo::EarthquakeCMT`
    - `SRInfo`
    - `FreqFull`
    - norm-factor calculation metadata
  - This is a strong reuse target for `ex2`, `ex3`, `ex5`, and `ex7`
  - `ex4` is only a partial fit because it currently uses `t2 = tout + 1.0`,
    while `InputParametersNew` currently fixes the end of the time window to
    `t_out()/60.0`
  - `ex6` is intentionally lower-level and should stay local

Recommendation:

- Prefer `InputParametersNew` wherever the example does not need unusual manual
  `FreqFull` configuration
- If `ex4` benefits materially from the preferred path, consider a small,
  non-breaking extension to `InputParametersNew` to allow a custom end-window;
  otherwise keep `ex4` on the legacy setup

2. Output-type normalization selection

- Repeated in:
  - `ex2`, `ex3`, `ex5`, `ex7`
  - related acceleration-only scaling in `ex4`
- Existing library owner:
  - `InputParametersNew::normFactor()`
- Missing reusable path:
  - a small standalone helper for legacy `InputParameters` workflows

Recommendation:

- When an example moves to `InputParametersNew`, use `paramsNew.normFactor()`
- For legacy examples that should remain legacy, add a small library helper for
  normalization-factor selection rather than repeating the same `output_type`
  switch

3. Standard filtering pipeline

- Repeated in:
  - `ex1`, `ex3`, `ex7`
  - partially in `ex2`, `ex4`, and `ex5`
- Existing library owner:
  - `DSpecM::FilterOptions`
  - `DSpecM::applyFilter(...)`
- Notes:
  - `ex3` and `ex7` can largely switch to `applyFilter`
  - `ex4` can likely use `applyFilter` and then keep only the time-domain
    output it needs
  - `ex5` can likely use `applyFilter` inside the convergence loop
  - `ex2` is only a partial fit because it inserts a source-time-function
    convolution between filter stages and applies a special Mineos correction

Recommendation:

- Reuse `applyFilter` wherever the standard pipeline is intended
- Keep the source-time-function-specific part of `ex2` local unless another
  example needs the same workflow

4. Filter option construction

- Repeated in:
  - `ex1`, `ex2`, `ex3`, `ex7`
- Existing library owner:
  - `DSpecM::FilterOptions`
- Missing helper:
  - none required in the public API

Recommendation:

- Avoid adding a new public helper just for default construction
- If repetition still feels noisy after refactoring, use a benchmark-local
  helper or a short local factory in the examples/support code

5. YSpec loading

- Repeated manually in:
  - `ex2`, `ex3`, `ex7`
- Existing library owner:
  - `DSpecM::loadYSpecTimeSeries(...)`
- Notes:
  - `ex1` already uses the helper
  - `ex3` and `ex7` should clearly migrate to the helper
  - `ex2` can also use the helper before its STF-specific processing

Recommendation:

- Replace the manual `YSPECREADER::DataColumns` blocks with
  `loadYSpecTimeSeries(...)`

6. Mineos loading

- Repeated manually in:
  - `ex2`, `ex3`, `ex7`
- Existing library owner:
  - `DSpecM::loadMineosTimeSeries(...)`
  - `DSpecM::loadReferenceTimeSeries(...)` for bundled YSpec + Mineos loading
- Notes:
  - `ex2`, `ex3`, and `ex7` can all use `loadMineosTimeSeries(...)`
  - `loadReferenceTimeSeries(...)` is a candidate only where a paired YSpec +
    Mineos bundle genuinely improves readability
  - `ex2` still needs local STF-specific postprocessing after loading

Recommendation:

- Replace manual Mineos loading with `loadMineosTimeSeries(...)`
- Only use `loadReferenceTimeSeries(...)` where it truly simplifies the code

7. SpecNM loading and filtering

- Repeated in:
  - `ex1`, `ex2`, `ex3`, `ex7`
- Existing library owner:
  - `DSpecM::loadSpecnmTimeSeries(...)`
  - `DSpecM::applyFilter(...)`

Recommendation:

- Keep using the existing helpers; this path already matches the intended
  library surface

8. Frequency-domain comparison writing

- Existing library owner:
  - `DSpecM::writeFrequencyComparison(...)`
- Current coverage:
  - `ex1` uses the helper already
- Missing reusable path:
  - a helper for the four-way comparison layout used by `ex2`, `ex3`, and
    `ex7` (primary + YSpec + Mineos + SpecNM)
- Keep local:
  - `ex5` convergence-matrix output
  - `ex6` tidal/kinetic-energy output

Recommendation:

- Reuse `writeFrequencyComparison(...)` where two-way or three-way output is
  sufficient
- Add a new helper for the repeated four-way paper-comparison layout used by
  `ex2`, `ex3`, and `ex7`
- Do not try to force `ex5` or `ex6` into the same writer abstraction

9. Time-domain comparison writing

- Existing library owner:
  - `DSpecM::writeTimeComparison(...)`
- Current coverage:
  - `ex1` uses the helper already
- Missing reusable path:
  - a helper for the repeated four-way comparison layout used by `ex2`, `ex3`,
    and `ex7`
  - possibly a record-section writer for `ex4`
- Keep local:
  - `ex5` convergence time-matrix output
  - `ex6` radial/tidal outputs

Recommendation:

- Add a reusable helper for the repeated multi-reference time output shared by
  `ex2`, `ex3`, and `ex7`
- Consider a dedicated record-section helper for `ex4` only if it improves
  readability without obscuring the example

10. YSpec path construction from `output_prefix`

- Repeated in:
  - `ex3`, `ex7`
- Existing library owner:
  - none
- Recommended owner:
  - benchmark-local support helper, not public API

Recommendation:

- Add a small benchmark-local helper that converts `output_prefix` into the
  in-repo YSpec reference path, rather than duplicating the
  `std::filesystem::path(...).filename().string()` logic

11. Source-time-function-specific postprocessing

- Repeated in:
  - `ex2` only
- Existing library owner:
  - none
- Recommended owner:
  - remain example-local for now

Recommendation:

- Keep this local unless another example or public workflow needs the same STF
  machinery

12. Convergence-study / low-level SEM outputs

- Repeated in:
  - `ex5`, `ex6`
- Existing library owner:
  - none
- Recommended owner:
  - remain example-local for now

Recommendation:

- Keep these outputs explicit in the examples
- Avoid over-abstracting them into the public library unless a second concrete
  consumer appears

### Phase B: Simplify The Paper Examples

Status: completed

Summary:

- Added [benchmarks/PaperExampleSupport.h](/home/adcm2/Documents/c++/DSpecM1D_Draft/benchmarks/PaperExampleSupport.h)
  to hold repeated paper-example glue that did not already belong in the
  public library.
- Simplified `ex2`, `ex3`, `ex4`, `ex5`, and `ex7` to reuse existing library
  functionality where that made the paper workflows shorter and clearer.
- Kept `ex6` intentionally low-level because its purpose is to expose the
  explicit SEM/matrix workflow.
- Cleaned up and then reintroduced structured comments across `ex1` through
  `ex7` so the flow is easier to follow without bringing back noisy commented
  code.

Verification:

- `cmake --build build/phase7_final --target paper_examples`
- protected example validation completed successfully on the user side

### Phase C: Output Writer Hardening

Status: completed

Summary:

- Expanded [tests/test_output_writers.cpp](/home/adcm2/Documents/c++/DSpecM1D_Draft/tests/test_output_writers.cpp)
  to directly cover:
  - low-level two-way and three-way frequency writers
  - `InputParametersNew` two-way and three-way convenience frequency writers
  - low-level and `InputParametersNew` time writers
  - invalid output-path failures
  - truncation to the shortest available input
  - explicit delimiter/layout expectations

Verification:

- `cmake --build build/phase7_final --target dspecm1d_output_writer_tests`
- `build/phase7_final/bin/dspecm1d_output_writer_tests`

### Phase D: Reference Loader And Filter Hardening

Status: completed

Summary:

- Extended [tests/test_reference_series_io.cpp](/home/adcm2/Documents/c++/DSpecM1D_Draft/tests/test_reference_series_io.cpp)
  with direct coverage for:
  - `loadReferenceTimeSeries(...)`
  - SpecNM padding on short inputs
  - Mineos truncation to the shortest component
  - zero-column requests across all three loader families
- Expanded [tests/test_signal_filtering.cpp](/home/adcm2/Documents/c++/DSpecM1D_Draft/tests/test_signal_filtering.cpp)
  to cover both `applyFilter(...)` overloads with correctly sized synthetic
  inputs, multiple passes, `enforceRealSignal`, and finite-output checks.
- Hardened [DSpecM1D/src/SignalFiltering.h](/home/adcm2/Documents/c++/DSpecM1D_Draft/DSpecM1D/src/SignalFiltering.h)
  so the wrapper now rejects invalid input matrix widths explicitly instead of
  letting mismatched dimensions fall through into the FFT/filter routines.

Notes:

- During implementation, the first draft of the new filter tests exposed that
  the underlying `SpectraSolver` routines require exact `FreqFull`-compatible
  input widths:
  - frequency input must have `freq.nt() / 2 + 1` columns
  - time input must have `freq.nt()` columns
- Phase D now codifies that contract directly at the `DSpecM::applyFilter(...)`
  boundary.

Verification:

- `cmake --build build/phase7_final --target dspecm1d_parser_tests dspecm1d_signal_filter_tests`
- `build/phase7_final/bin/dspecm1d_parser_tests`
- `build/phase7_final/bin/dspecm1d_signal_filter_tests`
- `cmake --build build/phase7_final --target t1 paper_examples`
- `ctest --output-on-failure -R '^smoke_tutorial_t1$'` in
  `build/phase7_final`

### Phase E: Preferred Solver API Hardening

Status: completed

Summary:

- Added [tests/test_solver_api.cpp](/home/adcm2/Documents/c++/DSpecM1D_Draft/tests/test_solver_api.cpp)
  to cover the preferred release-facing solver path directly.
- Added a small configurable parameter fixture in
  [tests/test_utils.h](/home/adcm2/Documents/c++/DSpecM1D_Draft/tests/test_utils.h)
  so the preferred-path tests can run with a deliberately tiny radial-only
  setup rather than reusing a heavier full-mode configuration.
- Added direct modular coverage for:
  - `SPARSESPEC::SpectraRunContext`
  - `Full1D::SEM(const InputParametersNew &)`
  - `SPARSESPEC::SparseFSpec::spectra(InputParametersNew &)`
  - `SPARSESPEC::SparseFSpec::spectra(InputParametersNew &, Full1D::SEM &)`
  - `SPARSESPEC::SparseFSpec::spectra(Full1D::SEM &, InputParametersNew &)`
- Kept the solver assertions structural:
  - the request object exposes the expected workflow objects
  - the preferred `SEM` constructor honors the `InputParametersNew` settings
  - all three preferred `spectra(...)` overloads produce finite outputs with
    the expected shape
  - the two SEM-reuse convenience overloads agree exactly for the same reused
    `SEM`

Notes:

- The first attempt put the new tests into the existing
  `dspecm1d_input_context_tests` binary. That failed at link time because some
  upstream headers used by this surface expose non-`inline` definitions.
- The stable fix was to keep the new preferred-solver coverage in its own test
  executable: `dspecm1d_solver_api_tests`.

Verification:

- `cmake --build build/phase7_final --target dspecm1d_input_context_tests dspecm1d_solver_api_tests`
- `build/phase7_final/bin/dspecm1d_input_context_tests`
- `build/phase7_final/bin/dspecm1d_solver_api_tests`
- `cmake --build build/phase7_final --target paper_examples t1`
- `ctest --output-on-failure -R '^smoke_tutorial_t1$'` in
  `build/phase7_final`

### Phase F: SEM And Helper Contract Hardening

Status: completed

Summary:

- Expanded [tests/test_sem_component.cpp](/home/adcm2/Documents/c++/DSpecM1D_Draft/tests/test_sem_component.cpp)
  so the SEM component suite now covers more than the original single
  monotonic local-to-global check.
- Added direct low-level contract coverage for:
  - receiver/source element lookup staying within mesh bounds
  - additional local-to-global ordering invariants for spheroidal and radial
    unknowns
  - consistent square matrix shapes for `hS`, `pS`, `hTk`, `pTk`, `hR`, and
    `pR`
  - expected dimensions and non-zero support for `rvZR`, `rvRedZR`,
    `rvBaseFull`, and `rvBaseFullT`
- Kept the tests local and structural rather than trying to reproduce the
  numerics of `ex5` or `ex6`.

Notes:

- The new fixtures reuse the same tiny synthetic parameter-file approach as
  Phase E so the low-level SEM tests remain fast and diagnostic.
- This phase is intentionally about contracts that would fail close to the bug
  source, not about asserting full waveform equivalence.

Verification:

- `cmake --build build/phase7_final --target dspecm1d_sem_component_tests`
- `build/phase7_final/bin/dspecm1d_sem_component_tests`
- `cmake --build build/phase7_final --target paper_examples t1`
- `ctest --output-on-failure -R '^smoke_tutorial_t1$'` in
  `build/phase7_final`

### Phase G: Coverage Closeout

Status: completed

Summary:

- Added [API_COVERAGE_MATRIX.md](/home/adcm2/Documents/c++/DSpecM1D_Draft/API_COVERAGE_MATRIX.md)
  as the explicit hardening closeout record for the release-facing surface.
- Classified the main public surface by coverage source:
  - direct modular tests
  - indirect modular tests
  - paper examples
  - smoke
  - legacy compatibility surface
- Recorded the intentional residual risks for:
  - remaining low-level SEM helper methods
  - legacy `SparseFSpec` overloads
  - specialist low-level access paths such as `meshModel()`
- Recorded the docs follow-on decision that the website should present testing
  material as multiple linked pages rather than one long page.
- Implemented that docs split in the static site with:
  - `library-testing.html` as the testing hub
  - `testing-smoke.html`
  - `testing-inputs.html`
  - `testing-io.html`
  - `testing-solver.html`
  - `testing-coverage.html`
- Added one more direct pass on the legacy/public API surface by covering:
  - the legacy single-SEM `SparseFSpec` overload
  - the legacy multi-SEM `SparseFSpec` overload
  - a small direct `meshModel()` contract check in the SEM suite

Outcome:

- The hardening workstream now has a stable stop condition:
  - the preferred release-facing API is directly covered
  - lower-level solver and SEM contracts have targeted modular coverage
  - the paper examples remain the protection for scientific workflow fidelity
  - the remaining uncovered items are intentional and documented

#### Example-by-example refactor direction

- `ex1`
  - Keep as the canonical concise example of the preferred public API
  - Use it as the model for how readable a paper example should be

- `ex2`
  - Strong candidate for helper reuse in:
    - context setup
    - YSpec loading
    - Mineos loading
  - Keep STF-specific logic local
  - Likely needs the new four-way writer helper before it becomes much shorter

- `ex3`
  - Strongest candidate for simplification
  - Should move toward:
    - `InputParametersNew`
    - helper-based YSpec/Mineos loading
    - `applyFilter`
    - shared four-way comparison writers

- `ex4`
  - Moderate simplification target
  - Record-section writing likely needs either a dedicated helper or explicit
    local code
  - `InputParametersNew` is only a partial fit because of the custom time
    window

- `ex5`
  - Moderate simplification target
  - Setup and normalization can likely reuse existing helpers
  - Convergence and error-output logic should remain explicit

- `ex6`
  - Keep intentionally low-level
  - No major abstraction push recommended
  - Only obvious duplication/cleanup should be touched

- `ex7`
  - Same direction as `ex3`
  - Strong candidate for helper-based cleanup using the existing library
    surface plus the new multi-reference writers

#### Phase A implementation checklist for Phase B

- Introduce or reuse `InputParametersNew` where it clearly reduces duplicate
  setup without changing behaviour
- Replace manual YSpec loading with `loadYSpecTimeSeries(...)`
- Replace manual Mineos loading with `loadMineosTimeSeries(...)`
- Replace standard filter pipelines with `applyFilter(...)` where the workflow
  is actually the standard one
- Add benchmark-local helper(s) for:
  - deriving in-repo YSpec paths from `output_prefix`
  - constructing repeated comparison-data paths if that makes the examples
    easier to read
- Add library or benchmark helper(s) for:
  - repeated four-way frequency comparison writing
  - repeated four-way time comparison writing
- Do not abstract away:
  - `ex2` STF-specific logic
  - `ex5` convergence/error logic
  - `ex6` low-level SEM/tidal logic

Checkpoint:

- Phase A is documentation/audit only
- No behavioural changes were made to the examples in this phase
- `cmake --build build/phase7_final --target paper_examples` completed
  successfully on 2026-04-23
- The next phase can proceed with a concrete “reuse vs add helper vs keep local”
  checklist

### Phase B: Simplify The Paper Examples

Status: implemented, pending protected-example validation

Intent:

- Make `ex1` through `ex7` shorter and easier to read
- Reuse existing library functionality where it already exists
- Keep specialised paper-specific logic explicit where abstraction would make
  the examples harder to understand

Completed on 2026-04-23:

- Added benchmark-local shared helper header:
  - `benchmarks/PaperExampleSupport.h`
- Centralised repeated benchmark glue there:
  - legacy normalization-factor selection for legacy `InputParameters` flows
  - standard paper-example `FilterOptions` construction
  - YSpec path derivation from `output_prefix`
  - bundled YSpec + Mineos loading from a Mineos base path
  - repeated four-way frequency comparison writing
  - repeated four-way time comparison writing
- Simplified `ex1` to reuse the shared filter-option helper
- Simplified `ex2` by reusing:
  - shared normalization helper
  - shared filter-option helper
  - bundled YSpec + Mineos loading helper
- Simplified `ex3` by reusing:
  - shared normalization helper
  - `DSpecM::applyFilter(...)` for the standard filtering pipeline
  - shared YSpec path helper
  - shared bundled YSpec + Mineos loading helper
  - shared four-way comparison writers
- Simplified `ex4` by reusing `DSpecM::applyFilter(...)` for the standard
  record-section filtering pipeline while keeping the record-section output
  explicit
- Simplified `ex5` by reusing:
  - shared normalization helper
  - `DSpecM::applyFilter(...)` inside the convergence loop
- Simplified `ex7` in the same style as `ex3`
- Left `ex6` intentionally low-level and unchanged in structure

Deliberately kept local:

- `ex2` source-time-function-specific frequency-domain processing
- `ex5` convergence/error-study output structure
- `ex6` tidal/low-level `SEM` workflow

Local verification completed:

- `cmake --build build/phase7_final --target paper_examples` completed
  successfully
- targeted modular test binaries still pass:
  - `build/phase7_final/bin/dspecm1d_output_writer_tests`
  - `build/phase7_final/bin/dspecm1d_signal_filter_tests`
  - `build/phase7_final/bin/dspecm1d_input_context_tests`

Still to verify before Phase B is treated as closed:

- protected paper-example validation in the supported environment to confirm
  there is no workflow drift in `ex1` through `ex7`

### Phase C: Output Writer Hardening

Status: completed

Intent:

- Add direct modular coverage for the release-facing writing helpers

Completed on 2026-04-23:

- Expanded `tests/test_output_writers.cpp` from the original two basic tests to
  cover:
  - low-level two-way frequency writing
  - low-level three-way frequency writing
  - `InputParametersNew` two-way frequency convenience writing
  - `InputParametersNew` three-way frequency convenience writing
  - low-level time writing
  - `InputParametersNew` time convenience writing
  - invalid output-path failures for both frequency and time writers
  - explicit delimiter/layout checks
  - truncation to the shortest available input

Local verification completed:

- `cmake --build build/phase7_final --target dspecm1d_output_writer_tests`
  completed successfully
- `build/phase7_final/bin/dspecm1d_output_writer_tests` now passes all 9 tests

Outcome:

- The public `OutputWriters.h` overload surface now has direct modular
  protection rather than being exercised only indirectly by the paper examples
