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
  "compiler_jobs": 2,
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
to the `command` arrays exactly as reviewed. `{BUILD_DIR}`, `{OUTPUT}` and
`{JOBS}` are the only substitutions; `{JOBS}` is the validated setting 1 or 2.
Use separate native and profiled variant commands, explicit outputs, and
include every compiled source plus every pinned dependency in the hash lists.
`dependencies` may name files or directories; directory fingerprints include
relative names, file contents and symlink targets. Variants build sequentially.
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
      "outputs": ["result.json"],
      "timeout_seconds": 7200,
      "memory_gib": 24
    }
  ]
}
```

`memory_gib` and `timeout_seconds` are optional. Output paths are relative to
that job's new directory and must be files by successful completion. Output
hashes are streamed after completion, and combined stdout/stderr is streamed
to `job.log`. The supervisor forces OpenMP/BLAS thread variables to one and
rejects caps above 50. It holds one nonblocking campaign-wide `flock`, polls
Linux `/proc` for the supervisor plus active job descendants, and on timeout,
thread-budget breach or memory-budget breach signals only recorded `(PID,
start-time)` identities, escalating from TERM to KILL. No process-group-only
cleanup is used. Failed output remains in place. The single-instance lock
file and all state/log artifacts remain in the campaign tree.

For detached execution, the campaign owner starts the reviewed plan with
`nohup python3 /path/to/supervise_campaign.py /path/to/plan.json > supervisor.log 2>&1 < /dev/null &`
and records the returned supervisor PID in the operator log. The Python
supervisor itself detaches each child into a new session. Do not start the
builder or this runner until the campaign plan and exact source/dependency map
have passed the independent review sequence.

`--resume --task-id ID` only skips jobs whose saved command hash, log hash and
declared output hashes still verify. If a job was started but did not reach
verified completion, the supervisor refuses that task ID; keep its evidence
and make any retry under a new task ID. A changed plan hash also blocks resume.

## Accounting boundary

This is a userspace watchdog based on periodic `/proc` snapshots. It is not
kernel aggregate accounting, and a very short-lived process/thread may start
and exit between polls. RSS is the sum of currently observed descendant
`VmRSS` values; shared pages may be counted more than once. Use site-provided
kernel/cgroup limits separately if guaranteed aggregate enforcement is needed.
