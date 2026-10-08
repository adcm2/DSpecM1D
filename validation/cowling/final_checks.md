# Parent verification before the final review

This records the original review-time checks; no validation was rerun for the
checkpoint. Named logs and binary outputs remain local and are excluded from
version control. Reproduction instructions are in `README.md`.

Branch: `dev/cowling`; feature HEAD: `8f7e7720bb2201bba035ae4abc8e919857a43e32`.
Immutable baseline: `889452b6790d2c60afacc156b4c7512990758aee`.
Test-source SHA-256: `ff21acaa2921f5c219ad507a31b520fd875bd0fe4a7124597da389e95303ab9b`.

Fresh builds use this repository as their CMake source directory, GNU C++ 13.3,
and cached dependency sources under `build/release/_deps`. Builds use `-j2`;
all checks below use `OMP_NUM_THREADS=1`.

- Full Release CTest: 63/63 active tests passed. The two opt-in diagnostics
  appear as disabled among 65 discovered tests; each was run explicitly in its
  numbered stage. Log: `final_release_ctest.log`.
- Debug solver API: 9/9 passed. Log: `final_debug_solver_tests.log`.
- Debug Cowling SEM component tests: 8/8 passed. Log: `final_debug_sem_tests.log`.
- `git diff --check`: passed.
- Original-baseline full-gravity regression: byte-identical serialized full
  `hS`, `hSa`, `pS`, forces and receivers at l=1,2,7, and single-/multi-SEM
  spectra at l=1–4 over requested 2–8 mHz, two receivers, attenuation off/on.
  Both 9,055,120-byte files have SHA-256
  `ec838d8173387a019879565053f52fa7d7dcbe9cc449c8bae81e8370317e9c2b`.
  The original baseline-header archive and binaries remain in the ignored
  directory `build/cowling_final_audit`. The probe source is retained verbatim
  as `full_regression.cpp`; it can reproduce the comparison using the recipe
  in `README.md`. Both compilations used the same dependency headers,
  input/model files and compiler flags. Production source has not changed
  during Stages 8–9.

Commands for the saved test logs:

```sh
OMP_NUM_THREADS=1 ctest --test-dir build/cowling_final_release --output-on-failure
OMP_NUM_THREADS=1 build/cowling_final_debug/bin/dspecm1d_solver_api_tests --gtest_filter='PreferredSolverApiTests.*'
OMP_NUM_THREADS=1 build/cowling_final_debug/bin/dspecm1d_sem_component_tests --gtest_filter='SEMComponentTests.Cowling*'
```

These checks establish the stated tested behavior. The limited Stage 8
scientific campaign and Stage 9 timing scope remain as described in `README.md`.
