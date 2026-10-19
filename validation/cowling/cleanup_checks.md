# Final Cowling PR cleanup checks — 2026-10-08

This checks only the approved chunk-local accumulation restoration, comments,
and removal of unrelated historical reports. No weak form, routing condition,
truncation/cadence, source/receiver, attenuation or output scaling changed.
Pre-cleanup headers: `cba9c34a607290043af41e1f1291816808bbb441`.
Immutable full-gravity baseline: `889452b6790d2c60afacc156b4c7512990758aee`.

## Numerical checks

- Release rebuild and CTest: 64/64 active tests passed, two diagnostics disabled.
- Focused solver tests: 10/10 passed separately with one/four threads, in both
  Release and Debug. The new multi-degree test compares all three routing cases
  against a matched-mesh single-SEM reference, attenuation off/on, and checks
  zero output outside the frequency band. Existing tests retain exact boundary
  routing and the full-only/crossing/Cowling-only three-chunk fixture.
- `multi_sem_cleanup.cpp` was compiled unchanged against pre-cleanup/current
  headers using identical Release flags/dependencies. Six outputs (full-only,
  Cowling-only, mixed; attenuation off/on) were compared as complete complex
  matrices, including columns outside the band. Serial outputs were byte-equal;
  four-thread before/after maximum relative Frobenius difference was
  `1.695e-16`. Current one/four-thread maximum difference was `1.258e-16`.
- The unchanged `full_regression.cpp` was freshly compiled against current
  headers and compared to the exact original-baseline executable/output from
  the independent PR review. Serial output remained byte-identical: 9,055,120
  bytes, SHA-256 `ec838d8173387a019879565053f52fa7d7dcbe9cc449c8bae81e8370317e9c2b`.
  Fresh four-thread runs of both executables differed by at most `3.009e-16`
  relative Frobenius norm across 34 matrices. Parallel degree addition order
  remains nondeterministic; bitwise parallel equality is not promised.

## One fixed multi-SEM benchmark

Whole `spectra(params)` calls include mesh/both-layout assembly, solves,
accumulation, rotation and scaling. PREM without attenuation, degrees 1–24,
two receivers, `nq=5`, relative error `0.001`, requested 5–35 mHz, 60-minute
window, 5-second sample spacing. Actual band: 157 bins at 4.6875–35.15625 mHz;
cutoff 19.921875 mHz, three chunks, derived `nskip=7`. Adaptive maximum steps
were 0.0341603, 0.0202827 and 0.0144232. Both versions used GCC 13.3, `-O3
-DNDEBUG -std=gnu++23 -fopenmp`, the same cached dependencies and Intel Core
Ultra 5 235U machine. Environment: `OMP_NUM_THREADS=4`, `OMP_DYNAMIC=FALSE`,
`OMP_PROC_BIND=close`, `OMP_PLACES=cores`.

After one warm-up, five calls per version were measured, restored version
first, pre-cleanup version second. Builds/tests were finished before these
measurements. An earlier overlapping run is excluded from timing conclusions.

| Seconds per complete call | Pre-cleanup | Restored |
|---|---:|---:|
| Five observations | 0.462197, 0.469856, 0.465767, 0.465389, 0.467114 | 0.468903, 0.472857, 0.477978, 0.470966, 0.470802 |
| Median | 0.465767 | 0.470966 |

The median was 1.1% slower in this fixed case; no runtime improvement is
demonstrated. Spheroidal locking nevertheless changes from 157 acquisitions
per degree to three (one per nonempty chunk). The earlier Stage 9 measured
Cowling-versus-full speed-up is a different comparison. No broad performance
claim or large campaign follows from this cleanup check.

## Reproduction and provenance

Configure/build the Release test target as in [README.md](README.md). Run
`PreferredSolverApiTests.*` with `OMP_NUM_THREADS=1` and `4`, dynamic teams off.
For the fixed driver, compile with the same generated test flags and put the
old headers first on the include path. From the repository root:

```sh
mkdir -p build/cowling_cleanup/before_headers
git archive cba9c34a607290043af41e1f1291816808bbb441 DSpecM1D | \
  tar -x -C build/cowling_cleanup/before_headers
python3 - <<'PY'
import pathlib, shlex, subprocess
root = pathlib.Path.cwd()
flags = {}
for line in (root / 'build/cowling_validation/tests/CMakeFiles/dspecm1d_solver_api_tests.dir/flags.make').read_text().splitlines():
    if line.startswith('CXX_'):
        key, value = line.split('=', 1)
        flags[key.strip()] = shlex.split(value)
out = root / 'build/cowling_cleanup'
for name, headers in [('before', out / 'before_headers'), ('after', root)]:
    subprocess.run(['/usr/bin/c++', *flags['CXX_FLAGS'], '-I' + str(headers),
                    '-I' + str(root / 'tests'), *flags['CXX_INCLUDES'],
                    str(root / 'validation/cowling/multi_sem_cleanup.cpp'),
                    '-o', str(out / name), '-lfftw3', '-lfftw3f', '-lfftw3l'], check=True)
PY
OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores \
  build/cowling_cleanup/after build/cowling_cleanup/after.bin
OMP_NUM_THREADS=4 OMP_DYNAMIC=FALSE OMP_PROC_BIND=close OMP_PLACES=cores \
  build/cowling_cleanup/before build/cowling_cleanup/before.bin
```

Repeat with one thread for deterministic byte comparison using `cmp`. Binary
format is six matrices in attenuation-off/on order, full/Cowling/mixed within
each: native `long` rows/columns followed by column-major real/imaginary doubles.
For parallel comparisons use relative Frobenius norm, tolerance `1e-10`;
platform-independent binary hashes are not promised. The baseline probe recipe
is unchanged in README.md. Actual local runs used `build/cowling_final_release`
and `build/cowling_final_debug`; raw logs/binaries are local only under
`/tmp/dspecm1d-cowling-cleanup/`.

- Driver SHA-256: `95d4f7137bf40c1dac8b383d3f255374f9171c40089bc9441511546684b88ff5`.
- Current solver test SHA-256: `96a39997144743c2eb58753fdeca95e37fa770a34c6fb05968d9f84d831da044`.
- Serial before/after driver output SHA-256:
  `1580a17088ad07cb909bb4effe914b314f635ed38db9b79fa3fe92fb86e44a62`.
- Both PREM model files and all four preserved historical benchmark files
  retain their pre-cleanup bytes. The historical paths are untracked with exact
  local `.git/info/exclude` entries; they are absent from the baseline PR diff.

Fresh Luna-high cleanup review: `PASS`, with no required corrections. It
independently inspected the complete cleanup and reran all 10 focused tests
with one/four threads in Release/Debug. Fresh Sol-6.1-high review returned
`PASS WITH NON-BLOCKING NOTES`, with no required corrections. It reviewed the
complete baseline-to-current feature, independently passed the 10 solver and
eight Cowling component tests in Release/Debug with one/four threads, and
recomputed the numerical comparisons and benchmark medians.
The earlier scientific coverage and low-level member-function-pointer
compatibility limitations remain unchanged. Historical original-source paths
in the preserved benchmark provenance were already unavailable at cleanup
entry; preservation checks cover the four repository-local copies.
