# Tunya profiling campaign tools

These portable Python 3 utilities prepare fresh native/profiled DSpecM1D v4/v5
and YSpec++ v6 builds and supervise a separately reviewed sequential job plan.
They do not contain project source or infer compiler/dependency flags. The
parent campaign owner must fill in and review absolute Tunya source and pinned
dependency paths, adapter file lists, exact compiler arguments and run commands.
The frozen transferred evidence stays read-only.

## Builder manifest

`build_campaign.py CONFIG.json` accepts this structure:

```json
{
  "build_root": "/space/adcm2/research/builds/kernel-profile-20261027",
  "environment": {"CXX": "/usr/bin/g++"},
  "variants": [
    {
      "id": "dsv4-native",
      "cwd": "/space/adcm2/research/projects/DSpecM1D",
      "output_name": "physical-driver",
      "command": ["/usr/bin/g++", "-O3", "-std=gnu++23", "-I{BUILD_DIR}",
                  "/space/.../dspec/source-v4-frozen/physical_driver.cpp",
                  "-o", "{OUTPUT}"],
      "sources": ["/space/.../dspec/source-v4-frozen/physical_driver.cpp"],
      "dependencies": ["/space/.../pinned-dependency"]
    }
  ]
}
```

The example uses abbreviated paths and is illustrative only. Add all compiler
flags, macros, source translation units, include/library paths and link flags
to the `command` arrays exactly as reviewed. `{BUILD_DIR}` and `{OUTPUT}` are
the only substitutions. Each command must invoke a supported compiler directly
(`g++`, `gcc`, `c++`, `clang++` or `clang`); shell/build wrappers, parallel job
flags and LTO flags are rejected. This keeps each variant to one compiler
driver invocation, with variants built sequentially. Use separate native and
profiled variant commands, explicit outputs, and include every compiled source
plus every pinned dependency in the hash lists. `dependencies` may name files
or directories; directory fingerprints include relative names, file contents
and symlink targets.
The builder sets OpenMP and common BLAS thread variables to one, writes each
exact argument array and input hash before compilation, streams stdout/stderr
to per-variant files, re-hashes inputs afterward, and stops at the first failure.
It refuses a nonempty build root and existing variant output directories.

## Run plan and supervisor

`supervise_campaign.py PLAN.json [--task-id ID]` runs jobs sequentially. Use a
new unique ID for every attempt. A run creates
`CAMPAIGN_ROOT/continuation/ID/JOB_ID/`, records state atomically, and refuses
existing task/job folders. Each job has this shape:

```json
{
  "campaign_root": "/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya",
  "thread_cap": 50,
  "jobs": [
    {
      "id": "profile-case-a",
      "cwd": "/space/.../working-directory",
      "command": ["/space/.../program", "--case", "a", "{JOB_DIR}/result.json"],
      "inputs": ["/space/.../program", "/space/.../model.dat", "/space/.../baseline.bin"],
      "outputs": ["result.json"],
      "timeout_seconds": 7200,
      "memory_gib": 24
    }
  ]
}
```

Each job must declare input files that define its result, including the
executable, model, scripts and baseline data as applicable. These files and
file-valued command arguments are hashed before and after execution; resume
requires current input hashes to match the saved hashes. Declared output paths
are excluded from the input set. Final records include peak observed threads
and userspace RSS across polls, elapsed wall time, and child user/system CPU
time deltas. `memory_gib` and `timeout_seconds` are optional. Output paths are
relative to that job's new directory and must be files by successful completion. Output
hashes are streamed after completion, and combined stdout/stderr is streamed
to `job.log`. The supervisor forces OpenMP/BLAS thread variables to one and
rejects caps above 50. It holds one nonblocking campaign-wide `flock`, polls
Linux `/proc` for the supervisor plus active job descendants, and on timeout,
thread-budget breach or memory-budget breach signals only recorded `(PID,
start-time)` identities, escalating from TERM to KILL. No process-group-only
cleanup is used. Failed output remains in place. The single-instance lock
file and all state/log artifacts remain in the campaign tree.

For detached execution, choose a fresh task ID and matching log path. Use the
shell's no-clobber mode to exclusively create an empty log, then append to that
file:

```sh
task_id=reviewed-unique-id
supervisor_log="/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya/supervisor-${task_id}.log"
( set -o noclobber; : > "$supervisor_log" ) || exit 1
nohup python3 /path/to/supervise_campaign.py /path/to/plan.json --task-id "$task_id" >> "$supervisor_log" 2>&1 < /dev/null &
```

Record the returned supervisor PID in the operator log. The Python supervisor
itself detaches each child into a new session. Do not start the builder or this
runner until the campaign plan and exact source/dependency map have passed the
independent review sequence.

The campaign owner reports that detached `nohup` execution over SSH succeeded
in the earlier CCDSM campaign. Tunya's user service manager has linger disabled,
so do not use a systemd user service for this campaign. This prior evidence does
not replace checking the actual campaign log and process state after launch.

`--resume --task-id ID` requires an existing, readable schema-1 state file
whose task ID and plan hash match before it creates directories or opens the
campaign lock. It only skips jobs whose saved command hash, log hash and
declared output hashes still verify. If a job was started but did not reach
verified completion, the supervisor refuses that task ID; keep its evidence
and make any retry under a new task ID. A changed plan hash also blocks resume.

## Accounting boundary

The supervisor enables and verifies Linux `PR_SET_CHILD_SUBREAPER` before
creating campaign work. It retains every observed `(PID,start-time)` identity,
including descendants adopted after reparenting, and cleanup polls, escalates,
reaps adopted children and verifies that identities are gone. Cleanup failures
are recorded with status `unclean`; they do not count as successful cleanup.
This remains a userspace watchdog based on periodic `/proc` snapshots, not
kernel aggregate accounting, and a very short-lived process/thread may start
and exit between polls. RSS is the sum of currently observed descendant
`VmRSS` values; shared pages may be counted more than once. Use site-provided
kernel/cgroup limits separately if guaranteed aggregate enforcement is needed.

## Completed-job qualification checker

`qualification_checks.py TASK_DIR OUTPUT_JSON` reads six completed L=3 job
directories: `dsv4-native-l3`, `dsv4-profiled-l3`, `dsv6-native-l3`,
`dsv6-profiled-l3`, `ysv6-native-l3`, and `ysv6-profiled-l3`. Each has
`data/result.json` and `job.log`. DSpec jobs also have `data/kernel.bin`,
`kernel_domain.json`, `kernel_radial_nodes.csv`, `trace.csv`, and
`spectrum.csv`; YSpec kernel jobs have `data/kernels.bin` and a readable
shared baseline at `kernel_baseline_trace` in `result.json`. YSpec radial
nodes are read from `data/radial_nodes.csv` if present, otherwise from the
recorded `kernel_sampling_input` path.

It enforces complete L=3 angular grids, all 179 radial nodes, finite numeric
outputs, DSpec q dot P(B raw)/(2E) error <= 1e-10, byte-exact DSpec
native/profile parity and stable result/domain metadata, DSpec v4/v6
physical-option/layout consistency plus per-kernel-channel and per-CSV-column
scaled maximum error <= 2e-13, and exact YSpec kernel and stable metadata
parity. DSpec trace/spectrum tables are compared per numeric column; YSpec's
shared baseline is checked for finiteness only, not treated as a computed
YSpec trace. The JSON report gives descriptive cross-code A/C/F/L/N
differences and fluid-shell DSpec kappa vs YSpec A+C+F over 4200-6200 km.
These comparisons are not global accuracy gates, and the checker makes no
near-source accuracy claim. Existing report paths are refused. The checker
reads completed data; it does not run a build or solver.

## Detached prerequisite gate

`wait_and_launch.py ROOT INPUT_ROOT EXPECTED_CHECKER_SHA256 EXPECTED_PLAN_SHA256 [TASKID]`
waits up to four hours for `ROOT/continuation/build-03/state.json` to report
`completed` and `ROOT/transfer/transfer-verification.json` to report
`verified`. It logs only phase changes, fails on terminal build errors, verifies
the exact plan/checker hashes, checks the build receipt and all seven current
binary hashes under the build root bound to the reviewed build config, confirms
the receipt hash in completed supervisor state, and requires those exact paths
as taskset-wrapped executables in the reviewed plan. It refuses an existing
`ROOT/continuation/TASKID`. Partial transfer-marker JSON remains pending during
the bounded wait. Once all checks pass it replaces itself with the
reviewed supervisor and plan. The default task ID is `profile-01`. Launch this
waiter under `nohup` with an exclusive fresh log; this helper creates no task
folders or logs before the final supervisor exec.
