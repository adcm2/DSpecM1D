# DSpecM1D Phase 0 baseline

The baseline for the dependency and documentation phases is remote `master` at
`fa265f584496f204589839fa6c9f10168d290b57`, observed on 19 October 2026.
PR [17](https://github.com/adcm2/DSpecM1D/pull/17) was squash merged at
09:57:47 Europe/London that day. The tracked source tree is
`8b57f310e8ce2ea4fb4f7a10d3a687be48bb99a5`.

The local checkout remains on `dev/cowling` at `6627fd2`. Its tracked tree is
identical to the merged baseline. Phase 0 used an isolated archive of that
tree and clean dependency archives; existing build folders were preserved.

## Branch situation

| Reference | Local or cached state at entry | Live remote state | Interpretation |
| --- | --- | --- | --- |
| Default branch | `origin/HEAD -> origin/master` | `master` | There is no remote `main` branch. |
| `master` | Local and cached remote `889452b` | `fa265f5` | Local refs predate the Cowling squash merge. Local master can later fast-forward by one commit. |
| `dev/cowling` | Local and cached remote `6627fd2` | `6627fd2` | Same tracked tree as current master. Git counts one master-only and five feature-only commits because the PR was squashed. |
| `dev/lapack-banded-pr1` | No local branch or cached remote ref | `adc2744` | 31 unmerged commits from `889452b`; current master contributes one separate Cowling squash commit. No open PR was returned by the live API. |
| `dev/love-numbers` | Cached remote `45e8b2b` | Absent | PR 15 was squash merged. Its tree is identical to local master at `889452b`; the cached remote ref is stale. |
| `archive/toroidal-modes-verbose` | Local and cached remote `4a0d46b` | `4a0d46b` | This commit is an ancestor of current master. Preserve unless cleanup is explicitly approved. |
| `v1.0.0` | `9b13dc5` | `9b13dc5` | Historical release tag, rather than the current baseline. |

Only one registered worktree exists in the inspected local repository. Its
working tree was clean at entry. Remote refs were read through an isolated
bare clone and GitHub; local refs were not fetched, pruned, switched, reset,
merged, renamed or deleted. No commits, pushes or GitHub edits were made.

The LAPACK branch modifies solver paths also changed by Cowling, including
`FullSpecSingleSem.h`, `FullSpecMultiSem.h`, `InputParametersNew.h`,
`StartElement.h` and their tests. Keep its integration separate from Phase 1.
The branch audit is not an assessment of that branch's numerical correctness.

After Phase 0 approval, a normal fetch and fast-forward of local master can
refresh the checkout. Retain the feature branches until their disposition is
explicitly decided. An ancestry-only deletion check does not recognize the
Cowling and Love-number squash merges as equivalent feature histories.

## Dependency baseline and cache mismatch

| Dependency | Revision used for this baseline |
| --- | --- |
| FFTWpp | `06cc1fb04c4398407839e637cb5414b179b0e9c3` |
| GaussQuad | `fbe37c7eef93695317dbf9c48715110bf64c800f` |
| GSHTrans | `b880ac9aeacdb54a6480fca3aee70a20a727c678` |
| Interpolation | `4508a7b3d2205f5647e87cdbaa4496cc0283c389` |
| NumericConcepts | `888126b44a979fda5dcc9d0b11d80cb8b3f61b52` |
| GoogleTest | `f8d7d77c06936315286eb55f8de22cd23c188571` (`v1.14.0`) |

Eigen is the installed 3.4.0 package; FFTW is the installed 3.3.10 package.
Compiler: GNU C++ 13.3.0. The fresh Release build used Ninja, two build workers,
and fully disconnected dependency resolution. NumericConcepts was fixed to
the recorded cached revision for this baseline; the production dependency
declarations still allow a moving transitive revision.

The existing `build/release/_deps/interpolation-src` is at `1557aad`, although
DSpecM1D declares `4508a7b`. The historical Cowling final builds point to that
Release cache. The separate `build/_deps` contains the declared revision.
All inspected dependency checkouts were clean.

Changing only the Interpolation include source back to `1557aad`, while
compiling the unchanged current regression driver and DSpecM1D source,
reproduced the original historical output byte-for-byte. Its SHA-256 remains
`ec838d8173387a019879565053f52fa7d7dcbe9cc449c8bae81e8370317e9c2b`.
The historical evidence retains that dependency provenance.

## Fresh validation

| Check | Result |
| --- | --- |
| Complete Release build | Passed; main library tools, examples, tests and optional Love-number module built |
| Normal CTest | 125 active tests passed: 64 main, 58 Love-number and 3 smoke checks |
| Optional diagnostics | Two expensive Cowling diagnostics remain disabled; they were not counted as passes |
| Four-thread solver API | 10 tests passed |
| Installed consumer without language selection | Failed: the imported target does not propagate the C++23 requirement |
| Same installed consumer with explicit C++23 | Configured, built and ran; model reported 13 layers and normalized radius 6371 |
| Source and dependency snapshot hashes | Unchanged after the checks |
| CI for merged master | CI and Pages workflows both completed successfully |

The consumer language requirement is a bounded packaging issue for Phase 1
and a setup-guide requirement for Phase 2. Phase 0 did not correct it.
Love-number installation/export integration remains outside this consumer
check, which exercises the main public umbrella header.

The unchanged `validation/cowling/full_regression.cpp` saved 34 matrices and
spectra: stiffness, attenuation stiffness, mass, forcing and receiver arrays
at degrees 1, 2 and 7, plus single- and multi-SEM spectra, with attenuation
off and on. Serial repetition was byte-identical. The fresh declared-pin
reference SHA-256 is
`a82d5e78edc22829164f178385e0d4aca524ad614854e38985b0bfc79b6db778`.
Serial versus four-thread differences were at most `1.61e-16` in relative
Frobenius norm.

The unchanged `multi_sem_cleanup.cpp` saved six spectra covering full gravity,
Cowling and a crossing cutoff, with attenuation off/on. It used degrees
1–24, three mesh chunks, 157 active bins at 4.6875–35.15625 mHz, and six output
rows. Serial versus four-thread differences were at most `1.97e-16`.
All saved numerical arrays were finite. Incidental driver timings establish
no performance conclusion.

The fresh declared-pin reference differs from the historical reference by up
to `1.44e-6` in relative Frobenius norm across the matrix comparisons and
`2.31e-8` across the selected spectra. The controlled legacy-Interpolation reproduction explains that
difference; it is not attributed to a Phase 0 source change.

These are regression references for the stated cases. They do not establish
universal Cowling accuracy, a preferred cutoff, broad convergence, or new
agreement with independent seismic solvers. Existing scientific evidence and
its limitations remain in `validation/cowling/` and `love_numbers/`.

## Evidence and reproduction

`baseline.json` records revisions, key input hashes, packages, checks and the
documentation requirements. `source_manifest.json` hashes all 255 baseline
tracked files. `numerical_summary.json` identifies each saved matrix and its
comparison. `branch_inventory.json` retains the branch commands and output.
`artifact_manifest.json` hashes the retained local files and numerical arrays.

Raw evidence is local and ignored under `build/phase0-2026-10-19/`: configure,
build and test logs; exact commands; test discovery; the consumer source;
the numerical drivers' compressed binary outputs; and the CLI output.
Binary arrays use the existing driver's native-long headers and column-major
real/imaginary doubles. Decompress `.bin.gz` before reading that format.
Hashes establish provenance; comparisons across toolchains should use
numerical tolerances rather than assume byte identity.

To reproduce, archive commit `6627fd2` or the identical tree at `fa265f5`, then
archive each dependency at the revision above into separate directories.
Use the retained `configure_command.json`, substituting the new source, build
and dependency directory paths. Do not reuse the old Release dependency cache.
Build with two workers and run:

```sh
OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE ctest --test-dir BUILD --output-on-failure --parallel 1
OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE BUILD/bin/dspecm1d_solver_api_tests --gtest_filter='PreferredSolverApiTests.*'
```

CTest smoke checks set their own two-thread environment. The retained
`run_checks.py` records the exact driver compilation and execution procedure;
`check_install.py` and `install_cpp23_commands.json` preserve both consumer
language-selection cases. The scripts use the directory containing them as
their run root, with `source`, `dependencies` and `release` subdirectories.
Restore those directories before running them.

## Documentation phase requirements

Prioritize a short path through downloading, installing and using the library.
Develop small examples that explain each step and show an expected result;
separate the introductory path from advanced reference material. Review the
current quick start for both missing guidance and excessive detail.

Visually assess the website for readability, navigation, usability and
professional presentation. Compare suitable professional scientific-software
websites, then make a focused improvement with simple styling. Avoid a large
redesign campaign. This work remains in Phase 2.

Phase 0 is complete as an observed baseline with the packaging finding above.
Phase 1, Phase 2 and branch cleanup await approval.
