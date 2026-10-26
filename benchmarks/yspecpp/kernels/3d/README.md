# Phase D: bounded 3D elastic-kernel comparison

The operational fixture is `fixture.json`; it preserves the accepted physical setup, objective map, perturbation profiles, finite-difference rule, gates, and dependency pins. The original reviewed planning record remains in `protocol.json`, and its SHA-256 is carried into the fixture. DSpecM1D now provides opt-in windowed 3D field and material-kernel APIs; this benchmark contains the independent validation and bounded evidence harness. The default 1D workflow remains unchanged.

From a DSpecM1D checkout with YSpecpp alongside it, run the matched baseline stage into an absent or empty output directory:

```sh
python3 benchmarks/yspecpp/kernels/3d/run_phase_d.py \
  --output build/3d-kernels/phase-d/baseline-01 --baseline-only
```

The output guard runs before any write or build and rejects a nonempty destination. Baseline-only produces matched spectra, regular-volume YSpec comparisons, an independent Cartesian modal baseline, and regular-volume linear Born predictions. The oracle checks raw and P-processed spectra against DSpec component by component over the full grid at the frozen `1e-9` bridge limit. The runner compares regular-volume Born derivatives with DSpec kernel integrals at the `1e-6` relative gate and `1e-6 S` informative floor. It also compares the independent source-local volume, contact, and total Born terms against DSpec production integrals for J_Z/J_N and C/F/L, including the source-owning upper-element lift. Baseline-only runs no ±h finite differences and give no Richardson result, so they do not complete Phase D. Failed gates and their outputs are retained for diagnosis.

Run the complete bounded campaign with one command and a fresh output directory:

```sh
python3 benchmarks/yspecpp/kernels/3d/run_phase_d.py \
  --output build/3d-kernels/phase-d/full-01 --full
```

The full runner builds and saves the matched baseline first, freezes its scales and source/executable identities, then invokes the independent Cartesian solver for all six prescribed ±h pairs. It saves centered derivatives, Richardson estimates and the selected plateau, source-volume/contact splits, coupling evidence, true solve residuals, and resource estimates beside the HTML report. `--baseline-reference PATH` is optional and adds a same-physics S/Kscale replay check against a previously accepted baseline packet. The report uses readable three-significant-figure tables and plots; JSON packets retain full precision.

By default, the runner finds `../YSpecpp` relative to the DSpecM1D checkout and uses `build/phase-d-build-cache` as a mutable build cache outside run outputs. Set `PHASE_D_YSPEC_SOURCE_DIR`, `PHASE_D_DEPENDENCY_ROOT`, `PHASE_D_DSPEC_PLANETMODEL_SOURCE_DIR`, `PHASE_D_YSPEC_PLANETMODEL_SOURCE_DIR`, or `PHASE_D_BUILD_CACHE_ROOT` to use other locations. Dependency sources are checked against the pins in `fixture.json`; absent caches may be populated by the pinned CMake FetchContent revisions. The DSpec and YSpec planetmodel pins differ and can be supplied separately.

Each run records source HEADs, dirty-file hashes, all Phase D source/config hashes, resolved dependency revisions and dirty hashes, exact commands, working directories, the reproducibility-relevant environment plus a complete environment fingerprint, exit statuses, CMake cache fingerprints, and copied executable hashes. It archives the dirty DSpecM1D and YSpecpp source files under `source_snapshot/`, with a hash-checked manifest tied to the run provenance. `report.html` is a responsive summary with scientific PNG plots; `baseline_metrics.json` and both oracle packets retain full precision. `fixture.json`, logs, solver CSVs, and `run.json` preserve the inputs and provenance.

The fixture retains l≤3, the one source/receiver pair, active raw bins 1–3, six fixed h steps, one job, OMP=1, BLAS=2, and no degree or mesh sweep. The frozen campaign resource cap is 1 GiB for sparse matrices and factorization estimates. YSpec sampled-volume kernels omit source-shell contact terms. Compare YSpec with DSpec only on the regular [4200,6200] km volume profile, whose perturbation vanishes at both endpoints; do not treat that comparison as source-contact evidence. Source-local C/F/L volume, contact, and total terms are a separate Phase D comparison against the independent Cartesian SEM oracle, using the finite-h source projection and source-owning upper-element lift. The standalone timing mode below measures the volume API only and excludes source contacts.

## Matched volume-kernel timing

`benchmark_kernel_runtime.py` times the same saved tutorial sample grid, model, full-gravity setup, raw bins 1–3, and one or two real windowed SI objectives in both solvers. It requires Release builds with explicit `-O3`, the same C++ compiler path and version, and BLAS disabled. Exact common dependency revisions and both planetmodel pins are in [`fixture.json`](fixture.json). The existing default Phase D cache layout is `build/phase-d-build-cache/build-yspecpp/_deps` for common dependencies and YSpec's planetmodel, plus `build/phase-d-build-cache/build-dspec/_deps/planetmodel-src` for DSpec's planetmodel. These are the concrete defaults below. If that cache is elsewhere, set `PINNED_DEPS`, `DSPEC_PLANETMODEL`, and `YSPEC_PLANETMODEL` to verified checkouts at the fixture revisions. The YSpecpp source itself is the live checkout used for this run, not a dependency pin; its HEAD, dirty state, and source hashes are recorded. DSpecM1D and YSpecpp deliberately use different planetmodel revisions, also recorded in `build_info.json`.

Set paths and the one compiler used by all builds:

```sh
export DSPEC_SRC=/path/to/DSpecM1D
export YSPEC_SRC=/path/to/YSpecpp
export PINNED_DEPS="$DSPEC_SRC/build/phase-d-build-cache/build-yspecpp/_deps"
export DSPEC_PLANETMODEL="$DSPEC_SRC/build/phase-d-build-cache/build-dspec/_deps/planetmodel-src"
export YSPEC_PLANETMODEL="$DSPEC_SRC/build/phase-d-build-cache/build-yspecpp/_deps/planetmodel-src"
export CXX=/usr/bin/g++
export PREFIX="$DSPEC_SRC/build/kernel-timing-install"
export D_BUILD="$DSPEC_SRC/build/kernel-timing-dspec"
export Y_BUILD="$DSPEC_SRC/build/kernel-timing-yspec"
```

Install DSpecM1D with its pinned dependency sources, then build the installed consumer and timing adapter. `FETCHCONTENT_SOURCE_DIR_*` values avoid an unplanned dependency fetch; they must point to the fixture-pinned checkouts. The commands use `FETCHCONTENT_FULLY_DISCONNECTED=ON` and therefore expect the shown cache to exist. If it does not, first obtain the exact revisions listed in `fixture.json`, or remove the disconnected flag and allow the CMake wrappers to fetch their hard-coded fixture pins. The standalone wrapper accepts `-DDEPENDENCY_ROOT=...` only as a directory containing the named `*-src` dependency directories; it does not infer the DSpecM1D checkout from that path. `DSPEC_SOURCE_DIR` is an empty cache setting and must be supplied for `-DSOLVER=dspec`; `YSPEC_SOURCE_DIR` must be supplied for `-DSOLVER=yspecpp`. The commands below pass the latter and use the installed DSpecM1D package for its timing target.

```sh
cmake -S "$DSPEC_SRC" -B "$DSPEC_SRC/build/kernel-timing-core" \
  -DCMAKE_CXX_COMPILER="$CXX" -DCMAKE_BUILD_TYPE=Release \
  '-DCMAKE_CXX_FLAGS_RELEASE=-O3 -DNDEBUG' -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DGSHTRANS_WITH_BLAS=OFF \
  -DDSPECM1D_BUILD_BENCHMARKS=OFF -DDSPECM1D_BUILD_TUTORIALS=OFF \
  -DDSPECM1D_BUILD_TOOLS=OFF \
  -DFETCHCONTENT_FULLY_DISCONNECTED=ON \
  -DFETCHCONTENT_SOURCE_DIR_NUMERICCONCEPTS="$PINNED_DEPS/numericconcepts-src" \
  -DFETCHCONTENT_SOURCE_DIR_INTERPOLATION="$PINNED_DEPS/interpolation-src" \
  -DFETCHCONTENT_SOURCE_DIR_FFTWPP="$PINNED_DEPS/fftwpp-src" \
  -DFETCHCONTENT_SOURCE_DIR_GAUSSQUAD="$PINNED_DEPS/gaussquad-src" \
  -DFETCHCONTENT_SOURCE_DIR_GSHTRANS="$PINNED_DEPS/gshtrans-src" \
  -DFETCHCONTENT_SOURCE_DIR_PLANETMODEL="$DSPEC_PLANETMODEL"
cmake --build "$DSPEC_SRC/build/kernel-timing-core" --parallel 1
cmake --install "$DSPEC_SRC/build/kernel-timing-core"

cmake -S "$DSPEC_SRC/tutorials/sensitivity_kernels_3d" -B "$D_BUILD" \
  -DCMAKE_CXX_COMPILER="$CXX" -DCMAKE_BUILD_TYPE=Release \
  '-DCMAKE_CXX_FLAGS_RELEASE=-O3 -DNDEBUG' -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DGSHTRANS_WITH_BLAS=OFF -DCMAKE_PREFIX_PATH="$PREFIX" \
  -DFETCHCONTENT_FULLY_DISCONNECTED=ON \
  -DFETCHCONTENT_SOURCE_DIR_NUMERICCONCEPTS="$PINNED_DEPS/numericconcepts-src" \
  -DFETCHCONTENT_SOURCE_DIR_INTERPOLATION="$PINNED_DEPS/interpolation-src" \
  -DFETCHCONTENT_SOURCE_DIR_FFTWPP="$PINNED_DEPS/fftwpp-src" \
  -DFETCHCONTENT_SOURCE_DIR_GAUSSQUAD="$PINNED_DEPS/gaussquad-src" \
  -DFETCHCONTENT_SOURCE_DIR_GSHTRANS="$PINNED_DEPS/gshtrans-src" \
  -DFETCHCONTENT_SOURCE_DIR_PLANETMODEL="$DSPEC_PLANETMODEL"
cmake --build "$D_BUILD" --target dspec_3d_kernel_consumer dspec_3d_kernel_timing --parallel 1

cmake -S "$DSPEC_SRC/benchmarks/yspecpp/kernels/3d" -B "$Y_BUILD" \
  -DCMAKE_CXX_COMPILER="$CXX" -DCMAKE_BUILD_TYPE=Release \
  '-DCMAKE_CXX_FLAGS_RELEASE=-O3 -DNDEBUG' -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DSOLVER=yspecpp -DYSPEC_SOURCE_DIR="$YSPEC_SRC" \
  -DDEPENDENCY_ROOT="$PINNED_DEPS" -DYSPEC_PLANETMODEL_SOURCE_DIR="$YSPEC_PLANETMODEL" \
  -DGSHTRANS_WITH_BLAS=OFF -DFETCHCONTENT_FULLY_DISCONNECTED=ON
cmake --build "$Y_BUILD" --target YSpec yspec_kernel_timing --parallel 1
```

Generate the lmax-3 tutorial values once; retain this output for both the timed comparison and later replay:

```sh
export SMALL="$DSPEC_SRC/build/kernel-examples-phase12/small-timing-01"
python3 "$DSPEC_SRC/tutorials/sensitivity_kernels_3d/run.py" \
  --consumer "$D_BUILD/dspec_3d_kernel_consumer" \
  --model "$DSPEC_SRC/tutorials/sensitivity_kernels_3d/homogeneous_shell_3d.deck" \
  --lmax 3 --output "$SMALL"
```

Run the primary scalar comparison with five alternating process-order repeats. Use a separate empty output directory for each run. The optional two-measurement run evaluates `J_Z` and `J_N`: YSpecpp batches them in one API call while DSpecM1D makes two scalar API calls, and the runner reports their summed API time as the matched workload.

```sh
export TIMED="$DSPEC_SRC/build/kernel-runtime/one-objective-01"
python3 "$DSPEC_SRC/benchmarks/yspecpp/kernels/3d/benchmark_kernel_runtime.py" \
  --dspec "$D_BUILD/dspec_3d_kernel_timing" --yspec "$Y_BUILD/yspec_kernel_timing" \
  --dspec-source "$DSPEC_SRC" --yspec-source "$YSPEC_SRC" \
  --dspec-build "$D_BUILD" --yspec-build "$Y_BUILD" \
  --model "$DSPEC_SRC/tutorials/sensitivity_kernels_3d/homogeneous_shell_3d.deck" \
  --samples "$SMALL/samples.csv" --lmax 3 --measurements 1 --repeats 5 \
  --accuracy-tolerance 1e-3 \
  --output "$TIMED"

export BATCHED="$DSPEC_SRC/build/kernel-runtime/two-objectives-01"
python3 "$DSPEC_SRC/benchmarks/yspecpp/kernels/3d/benchmark_kernel_runtime.py" \
  --dspec "$D_BUILD/dspec_3d_kernel_timing" --yspec "$Y_BUILD/yspec_kernel_timing" \
  --dspec-source "$DSPEC_SRC" --yspec-source "$YSPEC_SRC" \
  --dspec-build "$D_BUILD" --yspec-build "$Y_BUILD" \
  --model "$DSPEC_SRC/tutorials/sensitivity_kernels_3d/homogeneous_shell_3d.deck" \
  --samples "$SMALL/samples.csv" --lmax 3 --measurements 2 --repeats 3 \
  --accuracy-tolerance 1e-3 \
  --output "$BATCHED"
```

Replay saved timings without solver calls into a new directory:

```sh
python3 "$DSPEC_SRC/benchmarks/yspecpp/kernels/3d/benchmark_kernel_runtime.py" \
  --replay "$TIMED" --output "$TIMED/replay-01"
```

The default environment sets OMP and common BLAS thread counts to one. Each process performs a warmup before the measured calls; for two objectives DSpec warms and times two scalar calls, while YSpecpp warms and times one batched call. Executable process order alternates by replicate. API time excludes model loading, setup, warmup, and CSV export. The report separately retains first model-to-arrays time (setup plus first cold/warmup call), setup, setup-plus-API calculation time, and end-to-end process time. The raw `getrusage` high-water RSS field is retained as a diagnostic; it may inherit launcher pre-exec high-water and is not a solver-memory comparison. Each child call has a 300-second timeout, with a 900-second limit on the timed child-call loop; final report generation and replay are outside that limit. Actual timed arrays must pass the accuracy gate against one another and against the saved DSpec tutorial values before timing results are reported. These timings cover the recorded builds and grid only; they are not a universal speed ranking or a convergence study. The lmax-12 mode is a bounded pilot, not a convergence claim. See [`tutorials/sensitivity_kernels_3d/README.md`](../../../../tutorials/sensitivity_kernels_3d/README.md) for the installed consumer and replayable scientific plots.

## Saved timing result — 25 October 2026

On an Intel Core Ultra 5 235U, GCC 13.3, Release `-O3 -DNDEBUG`, one thread and BLAS off, the following whole-workload API times were measured. The five moduli are returned together; points count spatial locations per modulus and objective.

| Degrees | Points | Objectives | Repeats | DSpec median [min–max], s | YSpec median [min–max], s |
| --- | ---: | --- | ---: | --- | --- |
| 0–3 | 10,752 | J_Z | 5 | 0.716 [0.682–0.826] | 0.0248 [0.0234–0.0274] |
| 0–3 | 10,752 | J_Z + J_N | 5 | 1.54 [1.35–1.59] | 0.0336 [0.0289–0.0520] |
| 0–12 | 75,712 | J_Z | 3 | 60.4 [57.4–71.3] | 0.246 [0.221–0.416] |

All actual timed arrays pass the unchanged `1e-3` informative profile/integral gates. Maximum weighted profile/integral relative differences are `1.24e-9`/`4.68e-9` at lmax=3 and `2.32e-9`/`1.66e-8` at lmax=12. DSpec's timed J_Z arrays also match the installed-example values. Numerical inputs, raw timing stdout, hashes, plots and readable reports are saved locally under `build/kernel-examples-phase12/`; generated artifacts are not committed. These measurements expose a substantial performance difference for the stated native APIs and workloads, without establishing a general speed ratio or angular/mesh convergence.
