# Aggregate adaptive multi-SEM profiling

`stage_profile_solver` enables `DSPECM1D_ENABLE_PROFILING` only for this
dedicated executable. It calls the preferred adaptive `spectra(paramsNew)` API,
including one warm-up and the requested measured calls. The production path
still forms Eigen sparse frequency matrices; the LAPACK build then packs each
active sparse matrix into general-band storage before `zgbtrf`/`zgbtrs`.

The profile records use summed thread-seconds for work performed inside OpenMP
regions and a wall-clock total for the complete `spectra()` call. Categories
are kept separate for SEM construction, base operators, start selection,
frequency matrices, sparse compression, LAPACK packing, factorization, solve,
and source/receiver projection. The existing Eigen statement fuses trailing
sparse-block extraction with frequency-dependent arithmetic; that combined
cost is reported as frequency-matrix construction rather than changing the
production expression to manufacture separate measurements. The active RHS
`fRed` is an Eigen block view in the toroidal loop; the spheroidal loop's
existing dense copy is timed with start/truncation/extraction work. Counts
include SEMs, degrees, systems, compute/factorize/solve calls, RHS, dimensions,
nonzeros, and bandwidth ranges.

Each mode also records summed OpenMP worker-region seconds. Categorized worker
seconds plus a derived nonnegative unclassified remainder reconcile to that
worker total; the remainder includes loop/counter scans, waits, and other work
outside nested scopes. Initial serial setup is also reported as unclassified,
while SEM construction and final output projection have explicit serial
categories. Thread-work seconds are workload measures and do not add to the
separately reported complete-call wall time.

Build both Release configurations with GCC 15.2 and run the intentionally
small pilot first:

```bash
cmake --build /tmp/dspecm1d-stage-g-off-verify --target stage_profile_solver -j2
cmake --build /tmp/dspecm1d-stage-e-verify --target stage_profile_solver -j2
python3 scripts/stage_profile.py \
  --eigen-build /tmp/dspecm1d-stage-g-off-verify \
  --lapack-build /tmp/dspecm1d-stage-e-verify \
  --output-dir results/profiling --fmax 1 --lmax 20 \
  --threads 2 --repetitions 1 --tag stage_profile_pilot
```

The full run is deliberately prepared for a normal SSH shell outside Codex;
it is not launched by the implementation stage:

```bash
cd /home/adcm2/space/CleanBuilds/DSpecM1D && {
  nohup env \
    OMP_NUM_THREADS=2 \
    OMP_DYNAMIC=FALSE \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    BLIS_NUM_THREADS=1 \
    python3 scripts/stage_profile.py \
      --eigen-build /tmp/dspecm1d-stage-g-off-verify \
      --lapack-build /tmp/dspecm1d-stage-e-verify \
      --output-dir results/profiling \
      --fmax 20 35 55 --lmax 750 --threads 2 --repetitions 3 \
      --tag stage_profile_full \
    > results/profiling/stage_profile_full.log 2>&1 < /dev/null &
  echo $! > results/profiling/stage_profile_full.pid
}
```

Monitor it from another SSH session:

```bash
cd /home/adcm2/space/CleanBuilds/DSpecM1D
cat results/profiling/stage_profile_full.pid
pid=$(cat results/profiling/stage_profile_full.pid)
kill -0 "$pid" 2>/dev/null && ps -p "$pid" -o pid,etime,stat,cmd || echo DONE
pgrep -af 'stage_profile.py|stage_profile_solver'
tail -n 50 -f results/profiling/stage_profile_full.log
ls -lh results/profiling
python3 - <<'PY'
import json
from pathlib import Path
p = Path('results/profiling/stage_profile_full_metadata.json')
if not p.exists(): print('metadata not written')
else:
    d = json.loads(p.read_text())
    print(d['status'], d['completed_points'])
PY
```

The raw JSONL file is streamed and flushed after every measured call. Metadata
is atomically updated after every completed backend and frequency. It records
the base git revision and short worktree status, CPU model/topology, GCC
version, effective Release/C++ flags and LAPACK option from each build cache,
and SHA-256 hashes of the uncommitted profiling sources and harness. Thus a
detached result identifies the exact instrumentation even though this stage is
not yet committed. Expected full-run files:
`stage_profile_full_profiles.jsonl`, `stage_profile_full_metadata.json`, the
log, and PID file under `results/profiling/`.

The PID file is transient runtime state and is not a retained campaign
artifact. `stage_profile.py` appends to its raw JSONL output, so a completed
campaign must not be rerun with the same existing tag/output file. Use a new
tag, or explicitly move/remove the previous raw and metadata files first.

The retained tiny pilot was generated before this expanded metadata block was
added; its timing records remain valid and were deliberately not rerun or
rewritten. The full detached command above captures the expanded metadata.

## Completed detached campaign

The full campaign completed on `earth-tunya` (AMD EPYC 9334) at git revision
`04d5708c684101cce4efb8d3a99ff6a3e2ad597b`. Both builds used GCC 15.2,
Release `-O3 -DNDEBUG`, C++23, two OpenMP threads, dynamic OpenMP disabled, and
single-threaded OpenBLAS/MKL/BLIS. Each backend/frequency block performed one
untimed warm-up followed by three complete preferred adaptive
`spectra(paramsNew)` measurements. Backend launch order alternated by
frequency. Physics was ex3/PREM, all modes, attenuation enabled, `lmax=750`.

All 18 measured records are present: three repetitions for each backend at
20, 35, and 55 mHz. Every result is finite and has shape `6 x 16385`.

### Complete-call wall time

| fmax (mHz) | Eigen min/median/max (s) | LAPACK min/median/max (s) | Eigen/LAPACK |
|---:|---:|---:|---:|
| 20 | 65.537 / 65.614 / 65.662 | 64.253 / 64.309 / 64.436 | 1.0203x |
| 35 | 181.654 / 181.680 / 181.773 | 175.085 / 175.171 / 175.553 | 1.0372x |
| 55 | 492.166 / 492.464 / 493.635 | 471.234 / 471.712 / 471.915 | 1.0440x |

The spread is small relative to each median. The profiled result reproduces
the Stage I conclusion: LAPACK band LU gives a modest advantage that grows
with frequency.

### Detailed worker-time breakdown

The following values are median aggregate worker time in thread-seconds, with
share of the corresponding worker total in parentheses. They are not
wall-additive. `A(omega)` includes the current fused Eigen trailing-block
extraction and frequency-dependent sparse arithmetic.

| fmax/backend | Worker total | Base ops | Start/trunc. | A(omega) | Compress | Band pack | Factor | Solve | Projection | Unclassified |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 35 Eigen | 363.33 | 1.44 (0.40%) | 6.90 (1.90%) | 69.90 (19.24%) | 0.08 (0.02%) | 0 | 216.60 (59.62%) | 61.41 (16.90%) | 1.21 (0.33%) | 5.79 (1.59%) |
| 35 LAPACK | 350.31 | 1.41 (0.40%) | 6.55 (1.87%) | 68.00 (19.41%) | 0.08 (0.02%) | 58.14 (16.60%) | 125.45 (35.81%) | 77.72 (22.19%) | 1.20 (0.34%) | 11.76 (3.36%) |
| 55 Eigen | 984.87 | 3.05 (0.31%) | 10.63 (1.08%) | 188.43 (19.13%) | 0.17 (0.02%) | 0 | 597.28 (60.65%) | 169.10 (17.17%) | 2.01 (0.20%) | 14.19 (1.44%) |
| 55 LAPACK | 943.36 | 3.00 (0.32%) | 10.09 (1.07%) | 183.69 (19.47%) | 0.14 (0.01%) | 157.10 (16.65%) | 345.82 (36.66%) | 211.25 (22.39%) | 1.98 (0.21%) | 30.29 (3.21%) |

Per-call raw worker accounting reconciles exactly as `total = categorized +
unclassified`. For the stacked median table/figure, independently selected
category medians include tiny serial radial preparation/output scopes; the
displayed nonnegative unclassified residual closes the median stack to the
median worker total. Initial serial setup, SEM construction, and final output
projection remain separately available in the raw records.

Spheroidal work dominates: 91.0--92.1% of worker time at 35 mHz and
91.8--92.8% at 55 mHz, depending on backend. Toroidal work contributes about
7.0--8.8%; radial work is below 0.3%. A future optimization should therefore
be validated primarily on the spheroidal path while retaining all-mode
numerical checks.

### Workload and numerical validation

At 20, 35, and 55 mHz the solves contain respectively 983,155; 1,720,146; and
2,704,802 frequency systems. Eigen compute plus reused-factorize counts equal
the system count at every endpoint. On LAPACK, factorization count, band-pack
count, solve count, and system count are identical at every endpoint. Matrix
dimensions range from 17 to 930, 1554, and 2418 respectively; bandwidth stays
within `kl,ku = 4..15`.

Median Eigen/LAPACK output-norm relative differences are:

| fmax (mHz) | Relative norm difference |
|---:|---:|
| 20 | 2.20e-12 |
| 35 | 4.51e-14 |
| 55 | 7.17e-14 |

Compared with the same-machine Stage I medians, profiling increased complete
wall time by approximately 1.35%/2.15% (Eigen/LAPACK) at 20 mHz,
1.60%/1.78% at 35 mHz, and 0.88%/1.34% at 55 mHz. This is a cross-run overhead
estimate, not a paired microbenchmark, but it shows no large perturbation.

### Optimization conclusion

The evidence favors direct construction of each dynamic matrix into
preallocated LAPACK band storage as the next optimization experiment. On the
55 mHz LAPACK path, current sparse `A(omega)` formation consumes 183.69
thread-s (19.47%) and sparse-to-band preparation consumes another 157.10
thread-s (16.65%): 340.79 thread-s, or 36.1% combined. The same combined share
is about 36.0% at 35 mHz. Direct band construction can target the entire
form-sparse-then-pack boundary and can naturally represent trailing principal
systems by offsets/views.

Band-native storage of the base SEM operators is not the first priority:
measured base preparation is only about 0.3--0.4%, and sparse compression is
about 0.02%. A standalone change to start-element logic is also unsupported;
explicit start/truncation work is only about 1.1--1.9%. Trailing sparse block
extraction remains fused into the 19% dynamic-matrix category, so its benefit
should be evaluated as part of direct band construction rather than inferred
from the small explicit start-selection number.

Factorization and solve still form the largest LAPACK block (59.1% at 55 mHz,
58.0% at 35 mHz), so direct assembly will not remove the main intrinsic linear
algebra cost. After eliminating packing/temporary sparse construction, solver
layout, RHS conversion, and solve cost are the next measured areas to examine.

### Retained analysis artifacts

- `stage_profile_full_profiles.jsonl`: immutable raw profile records.
- `stage_profile_full_metadata.json`: immutable campaign metadata.
- `stage_profile_summary.json`: validated derived metrics and shares.
- `stage_profile_summary.csv`: compact machine-readable table.
- `stage_profile_wall_time.png`: complete wall-time medians and observed range.
- `stage_profile_worker_breakdown.png`: detailed worker-time categories.

Regenerate only the derived artifacts with:

```bash
MPLCONFIGDIR=/tmp/dspecm1d-matplotlib \
python3 scripts/stage_profile_analyze.py
```
