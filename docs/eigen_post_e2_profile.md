# Post-E2 Eigen SparseLU profiling reassessment

This reassessment measures the committed adaptive multi-SEM Eigen path after
E1b NaturalOrdering and E2 fixed-pattern/value-only frequency updates.  It is
profiling-only: ordinary builds retain the existing `compute()`/
`factorize()` cadence and no production optimization is included here.

## Current-path audit

`FullSpecMultiSem.h` creates one SEM per adaptive frequency chunk.  Radial
systems have a fixed full active range and rebuild one worker-private
`FixedPatternFrequencyMatrix` lazily, then call `compute()` at every frequency.
Toroidal and spheroidal paths compute `ridx` with their existing
`allIndices*()` helpers.  A worker rebuilds its fixed compressed union pattern
when the `ridx` changes, and otherwise updates only `valuePtr()` values.

The existing multi-SEM cadence remains `nskip = max(1, (i2-i1)/20)`.  At each
chunk position where the cadence requires recomputation, the Eigen solver does
`compute()` (symbolic analysis plus numerical factorization); intermediate
positions call `factorize()` only.  Thus E2 pattern rebuilding and Eigen
symbolic analysis are separate events: a pattern rebuild is needed only when a
worker sees a new `ridx`, whereas `compute()`/symbolic analysis follows the
existing chunk cadence and can repeat for an unchanged structure.  Each
OpenMP worker still owns independent `SparseLU` and fixed-pattern state; no
symbolic state is shared between workers, degrees, or chunks.

## Instrumentation semantics

With `DSPECM1D_ENABLE_PROFILING`, the profile executable expands only the
profiling build's Eigen `compute()` equivalent into:

```text
analyzePattern(matrix)
factorize(matrix)
```

at the same cadence.  Non-profiling builds continue to call `compute()`.
`analyze_pattern` and `factorization` scopes therefore measure separate worker
times without changing production cadence.  `frequency_matrix_construction`
contains E2 active-pattern creation and the per-frequency A(w) value update;
`source_receiver_projection`, `solve`, and the existing categories retain
their prior meanings.  Counts include symbolic analyses, all numerical
factorizations, solves, fixed-pattern rebuilds, and the union of observed
`ridx` values.  Category times inside OpenMP regions are summed worker-thread
seconds and are not wall-time-additive.

## Measured point

Only the required primary point was run: one untimed warm-up followed by one
measured complete `spectra()` call at `fmax=30 mHz`, `OMP_NUM_THREADS=32`.
The physics is the ex3/PREM-style type-4 all-mode attenuating workload.  The
fixed environment was `OMP_DYNAMIC=FALSE`, `OPENBLAS_NUM_THREADS=1`,
`MKL_NUM_THREADS=1`, and `BLIS_NUM_THREADS=1`, using GCC 15.2 and the Release
profiling executable recorded in the metadata.

| quantity | measured value |
|---|---:|
| complete wall time | 7.045931729 s |
| worker time | 224.966921316 s |
| frequency systems / solves | 1,475,483 / 1,475,483 |
| analyzePattern calls | 32,483 |
| numerical factorizations | 1,475,483 |
| intermediate factorize calls | 1,443,000 |
| fixed-pattern rebuilds | 9,878 |
| distinct ridx values | 138 |
| systems / degrees / SEMs | 1,475,483 / 1,500 / 3 |

| worker category | seconds | share of worker time |
|---|---:|---:|
| numerical factorization | 150.951972259 | 67.0996% |
| solve | 38.611247045 | 17.1631% |
| active-pattern/A(w) preparation | 20.739476044 | 9.2189% |
| start/truncation extraction | 5.799151879 | 2.5777% |
| analyzePattern | 0.358393603 | 0.1593% |
| source/receiver projection | 0.997481 | 0.4434% |
| other categorized work | 1.033157876 | 0.4592% |
| unclassified remainder | 6.497387563 | 2.8881% |

The principal mode split was:

| mode | A(w) preparation | analyzePattern | factorization | solve | worker total |
|---|---:|---:|---:|---:|---:|
| radial | 0.240638 | 0.055015 | 0.450948 | 0.049761 | 0.917398 |
| toroidal | 1.803023 | 0.046262 | 12.909937 | 3.128733 | 20.897561 |
| spheroidal | 18.695815 | 0.257117 | 137.591087 | 35.432753 | 203.151962 |

All measured outputs were finite with shape `6 x 16385`; the norm was
`0.0060603090464708366`.

## Decision

Stage E3 analyzePattern caching is not justified by this measurement.  The
measured symbolic-analysis work is 0.358394 s.  If it were removed entirely,
the strict wall-time ceiling (even if every analysis lay on the critical path)
is 0.358394 / 7.045932 = **5.09%**.  Under the observed 32-worker work-sum
accounting, the realistic balanced ceiling is approximately
0.358394 / 32 / 7.045932 = **0.1590%**, consistent with the 0.1593% worker
share.  The repeated analyses account for 32,483 calls, of which 22,605 were
on unchanged fixed patterns and 9,878 coincided with a rebuild. Numerical
factorization is the dominant remaining kernel at 67.0996% of worker time,
with solve at 17.1631% and the E2 active-pattern / A(w) preparation at 9.2189%.
E2 has therefore made further matrix-formation
work secondary to factorization.

The factorization/solve worker-time ratio is **3.9095**. Compared with the
retained pre-E2 f35 profiling artifact, where frequency-matrix construction
was about 19.2% of worker time, the current 30 mHz E2 profile is 9.22%.
This is contextual only: the grids, matrix sizes, build, and instrumentation
differ, so it is not an exact speedup claim. The Eigen solve is qualitatively
a different kernel from the custom-band backsolve; retained LAPACK/custom
results are not directly comparable to this worker profile.

The reproduction parameter is generated from the retained repository
`data/params/ex3.txt` template into the build's `data/params` directory. This
generated path is intentional because `InputParametersNew` resolves the model
relative to that build data tree. Template and generated-parameter hashes are
recorded in metadata.

The recommendation is **B: skip E3 and investigate numerical factorization**
in a separately scoped future stage.  This profile does not implement that
stage, alter ordering/cadence, or claim direct wall-time speedups against the
historical E1/E2 or LAPACK runs; those builds and instrumentation differ.

## Reproduction artifacts

- Runner: `scripts/eigen_post_e2_profile.py`
- Raw profile: `results/eigen_post_e2_profile/eigen_post_e2_profile_raw.tsv`
- Summary: `results/eigen_post_e2_profile/eigen_post_e2_profile_summary.json`
- Metadata: `results/eigen_post_e2_profile/eigen_post_e2_profile_metadata.json`
- Build: `/tmp/dspecm1d-backend-on-gcc15-final`
- Raw SHA-256: `b2a2214eacab64236ba4d3a5b8410e5532278e8f27c12e695feec460ad8089d2`
- Summary SHA-256: `6435688e75d20e4574d19671e45b0d7c0beb28f18e28084e81394b7d01615711`
- Metadata SHA-256: `2f8edf0f368d12baeb1cd85907855ac55c82606dede87290f3324d0b7d472a5a`
