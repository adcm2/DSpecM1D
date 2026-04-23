# DSpecM1D API Hardening Plan

Last updated: 2026-04-23 (Europe/London)

This file is the follow-on plan after the release-readiness work. Its goal is
to strengthen API coverage and simplify the paper-comparison examples without
changing their intended scientific behaviour.

## Goal

- Add modular, self-contained tests for release-facing APIs that are only
  lightly exercised today
- Make `benchmarks/ex1.cpp` through `benchmarks/ex7.cpp` as concise and clear
  as possible while preserving their paper-reproduction role
- Reuse existing library functionality wherever it already exists, so we do not
  keep duplicate logic in the examples or introduce double definitions

## Constraints

- `ex1` through `ex7` remain protected paper-comparison examples
- Scientific behaviour of the paper examples must not drift
- New tests should be small, local, and synthetic unless a paper example is the
  only realistic place a workflow can be validated
- After each phase there must be a stable checkpoint so the examples and tests
  can be rechecked before moving on

## Current Audit Summary

Covered reasonably well already:

- `InputParser`
- `InputParametersNew` basics and norm-factor selection
- selected `SpecHelpers`
- selected `SEM` local-to-global invariants
- individual reference-data loaders
- basic `FilterOptions` validation
- basic low-level output-writer truncation/time-window behaviour

Covered mainly by examples rather than modular tests:

- `InputParametersNew` convenience writing path in `ex1`
- preferred `SparseFSpec::spectra(InputParametersNew &)` workflow
- legacy lower-level `SparseFSpec` workflows used by `ex2`, `ex3`, `ex5`,
  and `ex7`
- specialised low-level `SEM`/matrix workflow used by `ex6`

Main gaps to harden:

- output-writer overload coverage
- reference-series bundle coverage
- deeper filter-behaviour coverage
- preferred solver overload coverage
- `SpectraRunContext` coverage
- `SEM(const InputParametersNew &)` coverage
- receiver/source helper and matrix-contract tests

## Phase A: Example Refactor Audit

Intent:

- Identify duplicated logic in `ex1` through `ex7`
- Identify patterns that are already implemented as library functions so we do
  not keep duplicate definitions or repeated code in the examples
- Decide which repeated logic should move to shared helpers and which should
  remain example-local because it is genuinely paper-specific

Work:

- Review `ex1` through `ex7` for repeated blocks:
  - normalization selection
  - filter option construction
  - reference-data loading
  - frequency/time output writing
  - repeated `FreqFull` setup patterns
- Cross-check those repeated blocks against the existing library surface:
  - `InputParametersNew`
  - `DSpecM::FilterOptions`
  - `DSpecM::applyFilter`
  - `DSpecM::ReferenceSeriesIO` loaders
  - `DSpecM::OutputWriters`
  - `SPARSESPEC::SparseFSpec`
  - `SEM` and related helper paths
- Produce a refactor checklist that labels each repeated pattern as:
  - already covered by an existing helper and should be reused
  - missing from the library and worth introducing as a helper
  - intentionally example-specific and should remain local

Checkpoint:

- A written mapping exists from repeated example logic to the helper/API that
  should own it
- No code behaviour changes yet
- `ex1` through `ex7` still build unchanged

## Phase B: Simplify The Paper Examples

Intent:

- Make `ex1` through `ex7` shorter, more readable, and more obviously aligned
  with the library API

Work:

- Refactor examples to use existing helper functions wherever that reduces
  duplicate code without obscuring the workflow
- Prefer code that reads as:
  1. load parameters/context
  2. run solver
  3. normalize/filter
  4. load comparison data
  5. write outputs
- Replace manual writer code with `OutputWriters` helpers where appropriate
- Replace manual reference-data loading with `ReferenceSeriesIO` helpers where
  appropriate
- Introduce small shared helpers only where the library does not already have a
  good home for the repeated logic
- Keep `ex5` and `ex6` more explicit if abstraction would make their numerical
  purpose harder to follow

Checkpoint:

- `ex1` through `ex7` still build
- Protected example outputs/workflows are rechecked
- Examples are materially shorter or clearer without hiding the paper logic

## Phase C: Output Writer Hardening

Intent:

- Add direct modular coverage for the release-facing writing helpers

Work:

- Add tests for all `writeFrequencyComparison` overloads:
  - low-level two-way writer
  - low-level three-way writer
  - `InputParametersNew` two-way convenience writer
  - `InputParametersNew` three-way convenience writer
- Add tests for both `writeTimeComparison` overloads:
  - low-level writer
  - `InputParametersNew` convenience writer
- Add failure-path tests for invalid output paths
- Add shape/contract tests for mismatched column counts
- Add checks for expected delimiter/column layout on tiny synthetic inputs

Checkpoint:

- Writer tests pass
- `t1` and `ex1` through `ex7` still build
- No paper-example output drift from refactoring

## Phase D: Reference Loader And Filter Hardening

Intent:

- Strengthen the small utility APIs that prepare comparison data and filtered
  outputs

Work:

- Add tests for `loadReferenceTimeSeries(...)`
- Extend edge-case coverage for:
  - `loadYSpecTimeSeries(...)`
  - `loadSpecnmTimeSeries(...)`
  - `loadMineosTimeSeries(...)`
- Add focused tests for both `applyFilter(...)` overloads
- Verify:
  - dimension preservation
  - `passes > 1`
  - `enforceRealSignal`
  - stable behaviour on tiny synthetic inputs

Checkpoint:

- Loader/filter tests pass
- Smoke tests still pass
- Example builds remain unchanged

## Phase E: Preferred Solver API Hardening

Intent:

- Add direct modular coverage for the preferred release-facing solver path

Work:

- Add minimal tests for:
  - `SparseFSpec::spectra(InputParametersNew &)`
  - `SparseFSpec::spectra(InputParametersNew &, Full1D::SEM &)`
  - `SparseFSpec::spectra(Full1D::SEM &, InputParametersNew &)`
- Add a focused `SpectraRunContext` construction/accessor test
- Add a focused `SEM(const InputParametersNew &)` test
- Keep assertions structural and local:
  - successful execution
  - expected dimensions
  - consistent reuse behaviour

Checkpoint:

- Preferred API tests pass
- `t1` smoke still passes
- Paper examples still build and behave as expected

## Phase F: SEM And Helper Contract Hardening

Intent:

- Add narrow tests for the lower-level contracts that underlie the more
  specialised examples

Work:

- Extend `SEM` component tests for:
  - receiver/source element lookup behaviour
  - additional local-to-global mapping invariants
  - selected matrix dimension/shape checks
- Extend `SpecHelpers` tests for any still-uncovered helper contracts used by
  the solver
- Keep these tests local and diagnostic rather than trying to reproduce
  full-example numerics

Checkpoint:

- SEM/helper tests pass
- `ex5` and `ex6` still build and remain readable

## Phase G: Coverage Closeout

Intent:

- Finish with an explicit record of what is and is not covered

Work:

- Build an API coverage matrix that marks each release-facing function/class as:
  - covered by modular tests
  - covered mainly by examples
  - legacy/documented but not a hardening target for now
- Record any deliberate residual risks
- Update docs/status so the hardened surface is easy to maintain
- For the website/docs follow-on, present the test inventory as a small set of
  linked pages rather than one long page:
  - overview/testing strategy hub
  - smoke and integration checks
  - parser and workflow-context tests
  - reference I/O and output-writer tests
  - filtering, solver API, and SEM/component tests

Checkpoint:

- Coverage matrix is complete
- Remaining uncovered items are intentional and documented
- The hardening workstream can be declared complete

## Recommended Order

1. Phase A: Example refactor audit
2. Phase B: Simplify the paper examples
3. Phase C: Output writer hardening
4. Phase D: Reference loader and filter hardening
5. Phase E: Preferred solver API hardening
6. Phase F: SEM and helper contract hardening
7. Phase G: Coverage closeout

## Immediate Next Step

Start with Phase A and produce the mapping from repeated example logic to the
existing library helper or API that should own it. That prevents duplicated
definitions and keeps the later test work aligned with the real public surface.
