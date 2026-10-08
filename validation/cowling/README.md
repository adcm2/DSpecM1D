# Cowling Stage 8 validation

This is a focused diagnostic comparison of the full and Cowling spheroidal
solver on the repository PREM model. It reuses the existing solver API test
fixture as an explicitly disabled GTest case; ordinary CTest does not run the
campaign. No Cowling interpolation/blending is applied.

## Checkpoint packaging

The checkpoint retains the diagnostic C++ source in `tests/test_solver_api.cpp`,
its parameter generator in `tests/test_utils.h`, the existing PREM model,
this methodology/results/provenance report, the three independent review
reports, `final_checks.md`, and the unchanged baseline probe source
`full_regression.cpp`. No diagnostic was rerun during packaging.

Logs, the four derived result CSVs, `campaign.time`, and archived
`input_parameters_l*.txt` files remain local and ignored. The input archives
contain machine-specific absolute model paths; they are generated snapshots,
not fixtures read by a test. The harness recreates its inputs using the current
repository path. All output names and hashes below describe the original local
runs, not files required to be present in a checkout. Review reports describe
the pre-packaging snapshot and retain their original statistics.

## Reproduce

From the repository root, build the solver API test target and run only the
opt-in campaign:

```sh
cmake -S . -B build/cowling_validation \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DDSPECM1D_BUILD_TESTS=ON \
  -DDSPECM1D_BUILD_BENCHMARKS=OFF \
  -DDSPECM1D_BUILD_TUTORIALS=OFF \
  -DDSPECM1D_BUILD_TOOLS=OFF \
  -DDSPECM1D_BUILD_DOCS=OFF
cmake --build build/cowling_validation --target dspecm1d_solver_api_tests -j2
OMP_NUM_THREADS=1 build/cowling_validation/bin/dspecm1d_solver_api_tests \
  --gtest_also_run_disabled_tests \
  --gtest_filter=PreferredSolverApiTests.DISABLED_CowlingValidationCampaignOnPrem
OMP_NUM_THREADS=1 build/cowling_validation/bin/dspecm1d_solver_api_tests \
  --gtest_filter=PreferredSolverApiTests.*
```

The build uses the dependencies and pinned versions declared in the repository
CMake files. An offline build may supply existing dependency source directories
through CMake's `FETCHCONTENT_SOURCE_DIR_*` overrides; cached build trees are
not checkpoint inputs. The commands print all measured rows to stdout. The
`COWLING_*_CSV_BEGIN/END` markers identify comparison, resolution and cutoff
sections; `COWLING_RESOLUTION_META` rows contain mesh/DOF information. The
comparison section ends when the resolution section begins. CSV files are
optional derived views, not inputs to the diagnostics.

The original local run's captured stdout is `campaign.log`; tables were split into
`full_vs_cowling.csv`, `resolution.csv`, `resolution_meshes.csv`, and
`cutoff_continuity.csv`. `campaign.time` records elapsed time and maximum RSS.
The ordinary focused solver API run passed 9/9 tests; its output is
`preferred_solver_tests.log`.

The fixture uses `data/models/prem.200.noatten.txt`, one receiver at 45 degrees
latitude / 90 degrees longitude, output type 0, attenuation off, `nskip=1`,
`dt=5 s`, a 5 minute output window, and frequency limits 5–80 mHz. It evaluates
l=1, 20, and 100. The campaign mesh uses `nq=5`, `maxstep=0.005`; the selected
resolution check compares that mesh with `nq=5`, `maxstep=0.0025` at l=1 and
100. The input parameter templates are generated for each degree by
`makeCowlingSolveParams`; the source and methodology retain their settings.
Mesh settings and `nskip` are supplied through `InputParametersNew` in the
harness.

The frequency-bin spacing is 1.5625 mHz. The complex frequency also includes
the existing imaginary cyclic damping of about 2.6526 mHz (`5/(300*2*pi)`),
with attenuation disabled. These settings limit interpretation near narrow
resonances and do not validate time-domain ringing. The five requested target
frequencies map to these nearest computed bins:
4.6875, 9.375, 20.3125, 40.625, and 79.6875 mHz. The finite window therefore
does not evaluate exactly 5, 10, 20, 40, or 80 mHz. Relative complex and
amplitude differences use a per-output-component floor equal to `1e-12` times
the largest magnitude in that component's full and Cowling spectra. Phase is
reported only when both values exceed this floor; `NA` means the component is
below it.

## Findings

Across these matched PREM solves, Cowling's relative complex response error
usually falls as frequency rises, while degree and component matter. The
largest component-wise relative complex errors at l=1 are 2.59e-4, 8.88e-5,
1.49e-5, 4.73e-6, and 1.13e-6 at the five computed bins. At l=20 they are
2.41e-2, 1.12e-2, 4.44e-3, 1.53e-3, and 1.67e-4. At l=100 they are
2.47e-2, 7.85e-3, 1.13e-2, 1.01e-2, and 3.59e-4. The l=100 sequence is not
monotonic between 10 and 40 mHz, but is much smaller at the high-frequency end.
Component amplitude and phase differences are retained row-by-row in the CSV.

The selected mesh refinement changed l=100 response components at 79.6875 mHz
by at most 0.785% between 207 and 406 elements; l=1 changed by at most 0.190%.
At the same bins, the full-versus-Cowling relative errors were nearly unchanged
across those two meshes. The comparison supports the reported formulation
difference at these sampled bins, but it does not establish general mesh
convergence or validate all frequencies/degrees.

For the hard cutoff at 20.3125 mHz, the adjacent bins are 18.75 mHz (full) and
21.875 mHz (Cowling). For each degree, the mixed result matches the standalone
full result below the cutoff and standalone Cowling at and above it exactly.
The maximum switch excess over the three output components, measured against
the adjacent full-response scale, is 1.19e-5 at l=1, 0.353% at l=20, and
0.748% at l=100. Corresponding maximum amplitude excesses are 1.18e-5, 0.0413%,
and 0.722%; maximum phase excesses are 0.000122, 0.253, and 0.418 degrees.
The natural full/full complex changes across those bins reach 19.9%, 20.8%,
and 33.8%, respectively. On these sampled cases, the cutoff adds a small
change compared with the natural adjacent-bin variation, which supports a hard
switch for this check. It does not establish that every source, receiver,
degree, or cutoff is similarly smooth; no waveform ringing check was run.

No matched MINEOS spectrum using both full gravity and an equivalent Cowling
cutoff was available through a straightforward repository workflow. The
repository contains MINEOS model readers and saved seismogram traces, but not a
matched full/Cowling reference for this campaign; no new MINEOS infrastructure
was introduced.

## Provenance

- Feature source baseline: `8f7e7720bb2201bba035ae4abc8e919857a43e32`.
- Immutable original implementation baseline: `889452b6790d2c60afacc156b4c7512990758aee`.
- PREM model SHA-256: `35071296c7fccf52730c3cdc476eb0ac5b2d42550c77dec484bc33dcbbac32f8`.
- `tests/test_solver_api.cpp` SHA-256 at final campaign run:
  `afcf98eaefc81bb6f7d5fb80a3c0c01106d322c65cee3d95696dc7003cabd919`.
- Input template SHA-256 values:
  - l=1: `b139c588af9972a7d64ec3a4ad1fd33b317361b8b1cfe44d78a0f2dc90d7f1f8`
  - l=20: `ee6383dbabe75d8e1ea359f8f6a692519fec84f7533700b210be801cb11d9f34`
  - l=100: `2474a6a4c9e509f32f864f662f64ee2a6b7d6d36d665e4a9e867befc36a3cfef`
- Captured campaign stdout SHA-256:
  `94c6b9db2f8944d1e5175d61ee303b288c7b9abb932ae113d8e6e4a8c49046d0`.
- Focused solver API output SHA-256:
  `c4eec3b2a51f17a3eecf01dfeca3a3c0c17209aa19b1f39eafc96cf1185fcfa7`.

The independent parent full-gravity regression also passed: full matrices,
forces, receivers, and single/multi-SEM spectra serialized byte-identically
against the immutable baseline. Its shared output hash was
`ec838d8173387a019879565053f52fa7d7dcbe9cc449c8bae81e8370317e9c2b`.

## Stage 9 performance diagnostic

Run the disabled diagnostic from the repository root after building the Release
solver API test target:

```sh
cmake --build build/cowling_validation \
  --target dspecm1d_solver_api_tests -j2
OMP_NUM_THREADS=1 \
  build/cowling_validation/bin/dspecm1d_solver_api_tests \
  --gtest_also_run_disabled_tests \
  --gtest_filter=PreferredSolverApiTests.DISABLED_CowlingPerformanceOnPrem
```

The source uses PREM without attenuation, `nq=5`, `maxstep=0.005`, degree 20,
one receiver, and a requested 5–80 mHz spectrum. The solver evaluates 13 bins
from 6.25 to 81.25 mHz for this one-minute window. The representative matrix
frequency is the 18.75 mHz bin nearest the 20 mHz target. The mesh has 2490
full-gravity and 1661 Cowling spheroidal DOFs (ratio 1.4991). Production
frequency truncation starts at global rows 84 and 56, leaving actual reduced
operators of 2406 and 1605 rows. Their stored nonzero counts are 35224 and
19213; lower/upper half-bandwidths are 15/15 and 10/10.

The `stage9_performance.log` timings are medians of five warm single-thread
measurements after one untimed warm-up. Fresh `SparseLU::compute` timings are
2.138 ms full and 0.609 ms Cowling; these include COLAMD symbolic analysis and
numeric factorization. Four-RHS solves with those systems take 0.341 ms and
0.143 ms. The prepared-SEM spectrum calls over the 13-bin band take 41.645 ms
and 16.112 ms. These are local ratios for this build and case, not general
speedup claims. `spectra(request, sem, ...)` includes per-call degree/frequency
matrix, source, and receiver setup but excludes model and SEM construction.
The shared SEM construction measurement is 8.777 ms and assembles both matrix
layouts, so it is not a full-gravity-only setup comparison. All four-RHS
solutions passed finite-value checks and relative residuals below `1e-8`.

The accumulated production diff against the immutable baseline was inspected
for duplicated full/Cowling loops, one-use helpers, wrappers, unnecessary API,
temporary production diagnostics, narrating comments, and over-generalized
abstractions. Existing solver loops share the same frequency traversal and
select the matching matrix/map/source/receiver layout at the point of use;
the Cowling matrix accessors mirror the existing full-gravity accessors. No
clearly beneficial production simplification was justified. The Stage 8
review's optional nearest-bin helper cleanup was not needed for correctness,
and the disabled diagnostics remain useful validation evidence. No production
source was changed in Stage 9.

### Stage 9 provenance and checks

- Feature/base commit: `8f7e7720bb2201bba035ae4abc8e919857a43e32`;
  immutable full-gravity baseline: `889452b6790d2c60afacc156b4c7512990758aee`.
- Performance test source SHA-256:
  `ff21acaa2921f5c219ad507a31b520fd875bd0fe4a7124597da389e95303ab9b`.
- PREM model SHA-256: `35071296c7fccf52730c3cdc476eb0ac5b2d42550c77dec484bc33dcbbac32f8`.
- Release executable SHA-256:
  `f4c3e0654ef0290c00dffecd74c3e77dc912556da25d0777b2a1107a4a8828be`.
- Build: GCC 13.3.0, CMake Release flags `-O3 -DNDEBUG`, `-j2`; both
  measurements used `OMP_NUM_THREADS=1`.
- Performance log SHA-256:
  `b9281ae0afca98f2f3b78e7ff7367ffa992172bfc58a9c5892dba16671af2606`.
- After adding Stage 9, the complete Stage 8 CSV section was rerun and matched
  `campaign.log` byte-for-byte (89 lines). Its log SHA-256 is
  `6783e47b999ef0ea7efa4496c0b498bf1d3bef9621346e46929f231c2e7afe00`.
- The ordinary `PreferredSolverApiTests.*` filter passed 9/9 tests; its output
  is `stage9_focused_tests.log`. `git diff --check` passed.

## Reproduce the recorded full-gravity baseline comparison

`full_regression.cpp` is copied byte-for-byte from the probe used for the final
review (SHA-256 `0363215a8bd92bd7e8868754d62facbd2474c7ac77b36267cfb289fcb92480e6`).
It serializes the same full matrices, forces, receivers and spectra described
in `final_checks.md`. Use identical dependencies/compiler flags for both
header sets. After configuring and building the Release test target above,
the following commands compile the existing probe; they do not add a CMake
production target. All archives, executables and binary results stay in `build/`.

```sh
mkdir -p build/cowling_baseline/headers
git archive 889452b6790d2c60afacc156b4c7512990758aee DSpecM1D | \
  tar -x -C build/cowling_baseline/headers
python3 - <<'PYTHON'
import json, pathlib, shlex, subprocess
root = pathlib.Path.cwd()
entries = json.loads((root / "build/cowling_validation/compile_commands.json").read_text())
entry = next(e for e in entries if e["file"].endswith("/tests/test_solver_api.cpp"))
args = entry.get("arguments") or shlex.split(entry["command"])
args = list(args)
i = args.index("-o")
del args[i:i + 2]
args.remove("-c")
args.remove(entry["file"])
for name, extra in [("original", ["-I" + str(root / "build/cowling_baseline/headers")]),
                    ("feature", [])]:
    subprocess.run([args[0], *extra, *args[1:], "-I" + str(root / "tests"),
                    str(root / "validation/cowling/full_regression.cpp"),
                    "-o", str(root / "build/cowling_baseline" / name),
                    "-lfftw3", "-lfftw3f", "-lfftw3l"],
                   cwd=entry["directory"], check=True)
PYTHON
OMP_NUM_THREADS=1 build/cowling_baseline/original build/cowling_baseline/original.bin
OMP_NUM_THREADS=1 build/cowling_baseline/feature build/cowling_baseline/feature.bin
cmp build/cowling_baseline/original.bin build/cowling_baseline/feature.bin
```

The archived binary hash records the original machine/build. Byte equality
between the two runs is the regression comparison; cross-platform binary
serialization and absolute hash equality are not promised.
