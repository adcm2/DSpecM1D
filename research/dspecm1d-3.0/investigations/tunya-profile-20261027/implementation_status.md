# Tunya profiling campaign preparation

## Scope and current status — 27 October 2026

### Review iteration v2 — Luna HOLD findings

Before editing, preserved the reviewed v1 utilities byte-for-byte in
`tools-v1/`. The active `tools/` copy is the v2 candidate; the five reported
blockers are being addressed there. This snapshot operation did not touch the
transferred campaign or source material.

1. **Builder parallelism claim constrained.** Removed the generic
   `compiler_jobs` setting and `{JOBS}` substitution. Each variant must now be
   one direct supported C/C++ compiler invocation; shell/build wrappers,
   parallel-job flags and LTO flags are rejected. Variant commands still run
   sequentially. This is the promised first v2 code change.

2. **Resume state preflight added.** `--resume` now requires an existing,
   readable schema-1 state file whose task ID and plan hash match before it
   creates campaign directories or opens the campaign lock. Missing, malformed,
   mismatched or unsupported state fails closed. This is the second v2 code
   change.

3. **Subreaper ownership added.** The supervisor now enables and verifies
   Linux `PR_SET_CHILD_SUBREAPER` before creating campaign work. Each polling
   pass discovers descendants from both recorded live identities and the
   supervisor identity, retaining start-time identities after reparenting.
   `/proc` inspection errors are surfaced and trigger cleanup.

4. **Cleanup verification and reaping added.** Limit breaches and handled
   SIGTERM/SIGINT now send identity-checked TERM, escalate to KILL, reap
   adopted descendants without stealing the leader's Popen wait, and poll
   until all recorded identities are verified gone or a bounded ten-second
   deadline expires. Signal failures and unresolved identities are recorded;
   unresolved cleanup yields terminal `unclean` status.

5. **Detached log creation made exclusive.** The README now uses a fresh
   task-specific supervisor log and shell noclobber for exclusive creation,
   then appends under `nohup` so a repeated launch cannot truncate prior logs.

6. **Builder parallelism statement corrected.** The README now documents one
   direct compiler invocation per variant and no parallel-job or LTO flags;
   variants run sequentially. This matches the schema guard and avoids a
   `compiler_jobs` setting that did not enforce the command.

7. **Synthetic cleanup probe exposed and fixed a supervisor identity bug.**
   The first setsid-child timeout probe showed that the supervisor identity was
   included as a cleanup target, which could signal the supervisor itself.
   Descendant discovery now uses that identity only as a root and removes it
   from the set of job identities eligible for signaling. The failed probe
   stopped before terminal state was written; its synthetic child was already
   gone when checked. The negative-state preflight probe passed.

8. **Second orphan probe improved leader reaping and signal evidence.** The
   setsid child was discovered and reaped, but the first bounded cleanup loop
   still saw the direct leader as a zombie because `Popen.wait()` ran only
   afterward. Cleanup now polls `Popen` during verification and reaps adopted
   children separately. TERM/KILL outcomes are recorded once per identity and
   signal stage instead of being duplicated every polling interval.

9. **PID reuse race closed for signaling.** Cleanup now opens a Linux pidfd,
   rechecks the recorded start time, and sends TERM/KILL through that stable
   process handle. If pidfd signaling is unavailable, the signal operation is
   recorded as an error and final identity verification controls the `unclean`
   result.

10. **Process-table enumeration failures made explicit.** `/proc` directory
    enumeration failure is now returned as a watchdog inspection error so the
    active job enters cleanup instead of the monitor silently losing process
    visibility.

11. **Normal-exit poll adapted to the v2 process-table return type.** A fresh
    normal-completion check found one call site still treating `proc_table()`
    as the old single-value API. It now unpacks the table and inspection
    issues, allowing completed leaders to finish their worker check correctly.

12. **Input fingerprint support added.** Each planned job must now supply an
    explicit nonempty `inputs` file list. The supervisor hashes these files
    plus command-line executable/file arguments before launch, enabling resume
    to bind old results to the actual binary and other declared inputs.

13. **Job metrics added.** The runner now retains peak observed thread count
    and peak userspace RSS across polling samples, wall elapsed time, and
    `RUSAGE_CHILDREN` user/system CPU deltas for each sequential job.

14. **Declared outputs excluded from input fingerprints.** A successful
    result path may appear as a command argument and exist during resume; the
    input hasher now excludes paths explicitly listed as job outputs while
    retaining executable and other file-valued command arguments.

14a. **RSS accounting aligned with thread accounting.** Peak and threshold RSS
     samples now include the supervisor process along with observed job
     descendants, using the same process set as the thread total.

15. **Run-plan documentation updated.** The example now requires `inputs` and
    documents hash-bound resume plus peak thread/RSS, elapsed, and child CPU
    fields so the operator plan exposes the measurement contract.

15a. **Detached-session note added.** The README records the parent-provided
     prior CCDSM evidence for `nohup` over SSH and Tunya's disabled user-service
     linger; it directs the operator away from a systemd user service and
     keeps post-launch observation as a separate deployment check.

16. **Resume verification moved ahead of state writes.** Completed records,
    command/input hashes, logs and outputs are now checked before the runner
    rewrites task state. Incomplete records are rejected without erasing their
    saved evidence; verified completed records remain intact when later jobs
    continue.

17. **Synthetic v2 verification passed.** A successful synthetic task recorded
    peak threads/RSS, elapsed time and child user/system CPU, then resumed only
    after input/command/log/output hash checks. Negative cases confirmed that a
    missing state file causes no lock/state/job side effects, a changed input
    binary/script blocks resume without rewriting state, and wrapper, parallel
    and LTO compiler commands are rejected. Timeout and thread-limit probes
    both stopped the job; the timeout probe used a `setsid()` child that ignored
    TERM and verified its reaping and exit. All checks ran in temporary
    directories under `/tmp`, with no scientific executables or SSH.

Approved scope is limited to portable campaign build, run-supervision and
documentation tools in `tools/`. The transferred original campaign and its
source/evidence remain untouched. No production repository, original source,
frozen metadata or remote system has been edited or accessed. No campaign
build or scientific job has been launched.

### Tooling implementation

Added `tools/build_campaign.py`, which consumes reviewed explicit build
commands, refuses nonempty output roots, accepts only one direct compiler
invocation per variant (rejecting wrappers, parallel-job flags and LTO), forces
common numerical thread controls to one, hashes declared source and dependency
inputs before and after each sequential build, saves exact commands and logs,
hashes binaries, and stops at the first error. The reviewed manifest must
supply the exact Tunya source and dependency mapping; paths and flags are not
guessed by the tool.

Added `tools/supervise_campaign.py`, which creates unique task/job directories,
uses a campaign-wide nonblocking `flock`, records atomic state checkpoints,
runs jobs sequentially, streams combined logs, records per-job peak observed
threads/RSS, elapsed wall time and child user/system CPU, and enforces a
maximum 50-thread userspace descendant budget plus optional timeout/RSS
thresholds. It enables Linux child-subreaper behavior, retains observed
`(PID, /proc start-time)` identities across reparenting, signals via pidfd,
reaps adopted children, and verifies cleanup or records `unclean`. Explicit
inputs and file-valued command arguments are hashed before and after; resume
validates state, command, inputs, log and declared output hashes before writing
state. Failed task folders remain intact.

Added `tools/README.md` with the JSON interfaces, `nohup` invocation, review
boundary and accounting limitations. The watchdog is a periodic userspace
`/proc` monitor, not kernel aggregate accounting; short-lived tasks may escape a
poll, and summed RSS can double-count shared pages.

### Validation

Both scripts passed Python 3 syntax compilation. Temporary synthetic checks
covered successful completion and hash-verified resume, missing-state
no-side-effect failure, changed-input rejection without state rewrite, per-job
resource/CPU metrics, compiler wrapper/parallel/LTO rejection, thread-budget
cleanup, and timeout cleanup/reaping for a reparented `setsid()` child that
ignored TERM. The build tool was not run against project sources, and no
scientific executable or SSH was used. The Tunya source/dependency map remains
a parent-owned input.

The v2 candidate is ready for a fresh independent tools review. It is not
campaign-ready until that review and the parent-owned actual Tunya build/run
manifest review are complete. Tunya's detached-session behavior and any
kernel/cgroup task limit remain deployment gates; this implementation does not
claim kernel aggregate thread accounting.

### Remaining campaign inputs and gates

The parent campaign owner must supply and independently review exact pinned
DSpecM1D/YSpec++ source and dependency paths, per-variant compiler commands,
run jobs, output declarations and numerical/resource acceptance gates. Tools
are preparation only; no profiling result or Wigner optimization validation is
claimed here.
# Transfer preparation

2026-10-27: Added transfer/stage.py to archive every original campaign file, including failed/partial evidence and frozen binaries, plus exact clean dependency source snapshots. It hashes file payloads while archiving, detects concurrent source changes, and retains the archive plus bounded transfer parts. Original files are read only. Long transfer preparation will run under nohup.

2026-10-27: Staging completed: 5,418 campaign/dependency files, 6,480,238,934 compressed bytes in nine parts. Added directory metadata including empty original evidence directories. Added a four-connection nohup upload helper with unique per-artifact logs/receipts and no deletions; failed transfers remain available. Numerical execution has not begun.

2026-10-27: Added transfer/extract_on_tunya.py: independently hashes every uploaded part, archive and extracted payload; refuses existing roots and unexpected paths/types; preserves modes, timestamps and empty evidence directories. Synthetic fresh extraction and overwrite/path rejection checks passed. Intended to run detached under nohup on Tunya.

2026-10-27: Added a 779-file, hash-verified bootstrap source/dependency copy for fresh builds while the complete evidence transfer continues. It uses separate bootstrap paths and the same extractor; original evidence paths remain reserved for the complete archive.

2026-10-27: Added explicit Tunya build and sequential numerical plans plus a build receipt wrapper. Seven fresh compiler variants use pinned dependency snapshots and Eigen internal matrix multiplication. Numerical jobs pin CPU0, one library thread, watchdog maximum50; full L3 parity gates precede 64-angle diagnostic component profiles at128/384/1536. YSpec forward and kernel runs are separate. Original first generated configs preserved as configs-v1-pre-review; explicit input paths corrected before review. These probes do not fulfill high-degree global materialization.

2026-10-27: Added `tools/qualification_checks.py` for data-only checks of the six planned DSpec v4/v6 and YSpec v6 native/profile L=3 jobs. It validates complete angular grids, 179 radial nodes, binary shape and byte counts, finite numeric outputs, DSpec q dot P(B raw)/(2E), DSpec trace/spectrum parity, DSpec native/profile and v4/v6 parity, exact YSpec native/profile kernel and stable metadata parity, and MATHPROFILE exclusive/inclusive accounting. It reports cross-code parameter errors and fluid-shell DSpec kappa versus YSpec A+C+F descriptively without a global accuracy gate. YSpec's shared baseline input is checked for finite values only, not treated as a computed YSpec trace. Existing output JSON is never overwritten.

Validation used a temporary synthetic six-job fixture under `/tmp`: the complete fixture passed, and changing DSpec sampled-point metadata to an incomplete grid produced a failing JSON report. Python syntax compilation passed. No campaign data was inspected; no build, scientific job, or SSH was started. The checker is ready for independent review and makes no global scientific accuracy claim.

2026-10-27: Preserved the pre-channel-fix qualification checker as `tools-v1/qualification_checks.py` before tightening parity semantics, so the reviewed previous implementation remains recoverable.

2026-10-27: Static fresh Luna and short Sol reviews passed build-only launch. Frozen build bundle contains exact v6 source, reviewed receiver gate, explicit builder/supervisor and build config. Numerical plan remains held for stricter qualification-checker review. Added SSH transport helper using existing read-only SSH config; failed remote diagnostics now retain stdout/stderr.

2026-10-27: Tightened `tools/qualification_checks.py` parity semantics after independent review: DSpec native/profile kernel, trace, and spectrum parity is now byte/value exact; DSpec v4/v6 checks use separate relative scales for each kernel parameter channel and each numeric trace/spectrum column. YSpec native/profile remains exact. Kernel comparisons still stream bounded chunks and now retain per-channel error metrics. The original checker was preserved under `tools-v1/qualification_checks.py` before this edit.

2026-10-27: Added CSV header names to the per-column parity metrics so any failing trace/spectrum column is reported by its actual field name.

2026-10-27: Added frozen-driver binary layout gates (`kernel_domain.json.binary_layout` for DSpec and `result.json.kernel_binary_layout` for YSpec), exact DSpec native/profile stable result/domain metadata comparison (excluding elapsed time, paths, and disk-capacity settings), and v4/v6 consistency for stable physical options/domain layout plus per-field energy checks. DSpec trace/spectrum headers must match the frozen column contract.

2026-10-27: Tightening validation passed Python syntax compilation and the synthetic six-job qualification fixture. A synthetic mixed-scale kernel with a large first channel and a 0.5% relative error in a tiny second channel fails per-channel gating although the former whole-array scale would have hidden it; the analogous tiny waveform-column error next to a 3600-second time column fails by `vertical_m`. A non-byte-exact kernel fails the exact-parity check. These checks used temporary files under `/tmp` only.

2026-10-27: Build-01 failed before compilation because another migration moved CCDSM. Preserved its source bundle, task state, logs and empty build attempt. Live verification found 530 Eigen headers byte-identical at /space/adcm2/research/projects/CCDSM/deps/toolchain/eigen/include/eigen3 and confirmed system Python3 NumPy1.26.4. Configs-v2 changes only relocated inputs root, fresh retry build root, actual Eigen location and analysis Python; code/physics/compiler flags unchanged.

2026-10-27: Added bounded detached transfer waiter and launched it under nohup on Tunya. It waits for complete exact-size parts, then invokes the hash-pinned safe extractor. All existing files remain. Build-02 reached compiler and exposed five absolute laptop include directives in the frozen DSpec header; failed build/logs preserved. Preparing separate portable source copies changing only those include directives before fresh retry.

2026-10-27: Prepared build-only v3 bundle with separately reviewed portable v4/v6 copies; each differs only in five include directives. Configs-v3 points to those copies and a new build-root03/input-rootv3. Original frozen code, earlier deployed inputs, failed builds and logs remain. Fresh Luna portable/config delta PASS; short Sol delta review pending. No numerical source equations changed.

2026-10-27: Preserved the per-channel/per-column checker as `tools-v1/qualification_checks-pre-metadata.py` before adding driver-layout and DSpec metadata gates.

2026-10-27: Metadata-gate tests passed on a temporary synthetic six-job fixture. Deliberately changing the DSpec v6 profiled `qex` setting failed stable native/profile metadata parity; changing YSpec's declared binary layout failed the driver-layout gate. No source or campaign data was changed and no campaign run was launched.

2026-10-27: Final qualification checker passed fresh independent Luna and short Sol reviews, closing exact instrumentation parity, per-channel scaling and metadata/layout gates. Numerical plan remains conditional on completed/hashverified build03 and full evidence verification. Preparing detached prerequisite-wait/launch helper so the campaign can continue independently of SSH and the laptop after upload.

2026-10-27: Added `tools/wait_and_launch.py`, a bounded prerequisite gate for detached continuation. It waits up to four hours for build-03 state `completed` and immutable transfer status `verified`, logs phase changes, verifies the reviewed checker/plan hashes, validates all seven current binary hashes against `build-03/build/build-result.json` using the receipt's actual build root, refuses an existing task directory, then `execv`s the supervisor with the reviewed plan and task ID. It creates no campaign outputs or logs and launches no build or solver before the gates.

2026-10-27: Waiter syntax compilation and a temporary mock prerequisite test passed. The mock good path verified seven report-bound binaries and the reviewed plan/checker hashes without invoking `execv`; wrong checker hash and pre-existing task directory were rejected. The waiter was not run against campaign paths and no remote process was started.

2026-10-27: Tightened the binary-receipt check to require each output to reside directly in its actual report-bound variant directory, preventing a receipt from redirecting a variant hash check to an unrelated path.

2026-10-27: Preserved the pre-provenance prerequisite waiter as `tools-v1/wait_and_launch-pre-provenance.py` before binding its build receipt to the reviewed config, task plan and build state.

2026-10-27: Bound the waiter to the reviewed build config root/output names, the completed build job's receipt hash in supervisor state, and successful/unchanged per-variant command receipts. It now requires all seven binaries as absolute taskset CPU-0 executable arguments in the reviewed plan and rejects any other plan argument/input under the build root. A partially written transfer marker remains pending until complete verification or timeout. README documents the strengthened chain.

2026-10-27: Provenance-waiter mock checks passed: valid plan/config/state/receipt chain verified all seven binaries without `execv`; a partial transfer marker remained pending; a changed plan executable and a receipt build root differing from reviewed config were rejected; malformed build state failed closed. Python syntax compilation passed. Tests used temporary files under `/tmp`; no campaign path or remote process was touched.

2026-10-27: Restricted the waiter and its validation API to the reviewed task ID `profile-01`, matching the fixed paths embedded in the reviewed plan.

2026-10-27: The `profile-01` restriction passed syntax and synthetic checks: only the reviewed task ID succeeds; `profile-02` is rejected; the valid chain still verifies all seven binaries; the partial marker remains pending; and a mismatched plan executable is rejected without invoking `execv`.

2026-10-27: Build03 completed all seven variants in432.5401s with peak6observedthreads,2.3957GiB userspace RSS and unchanged input hashes. All binaries and53linked libraries recorded in remote build-runtime-provenance.json; no directCBLAS symbols in any variant. Deployed final checker and CPU0 sequential plan after fresh Luna+shortSol PASS. Final hash/receipt-bound waiter also passed reviews and launched undernohup asPID1985242; it waits for full evidence transfer/extraction before runningprofile-01. Remote post-SSH process verification is being recorded.

2026-10-27: Full archive transfer verified all5,418payloads and1,327originaldirectories. First detached numerical attemptprofile-01 stopped atreceiver audit, preserving its log/state. Allnew closedform outputs matchedtheanalyticoracle exactly; oldGSH recurrence endpointnormalization showed relative error7.7113e-12 atL1536 (zero/supportchecks passed). IndependentSolmathreview classifies this asoldreference recurrence drift, notnewclosedform error. Preparing validation-v2 only with independentlongdouble oracle, unchangednew8ULP/exactzero requirements and explicit retainedoldamplitude diagnostics. Configs-v4 rebuild onlythatgate and retryprofile-02 withall sixexisting solverbinaries andL3paritytolerances unchanged. No optimization/physics source edits.

2026-10-27: Added `tools-v3/prefix_trace.py`, a data-only adapter for the
validated 1024-row/5 s DSpec trace format. It checks the exact CSV header,
finite values, the anchored grid, Hann weights over 1000–1200 s, and exact zero
weights for removed rows 721–1023; it then saves the original header and first
721 rows byte-for-byte to an exclusive fresh path. It independently verifies
positive equal full/prefix Hann-weighted energies with `math.fsum`, rechecks
the input hash after writing, and records paths, hashes, row/time ranges and
energy in an exclusive receipt. Added six synthetic-only tests and documented
the CLI in `tools-v3/README.md`; all passed locally. No solver/build/science,
remote activity, or campaign data modification occurred. The adapter notes
are isolated in `tools-v3/prefix_trace-README.md`; the supervisor snapshot's
`tools-v3/README.md` is unchanged.

2026-10-27: Corrected the adapter's Hann-weighted energy metric to
`0.5*dt*sum(weight*vertical_m^2)` and labeled receipt/console values with the
`m2_s` unit suffix. Matching tests and adapter documentation were updated.
This is data-only receipt accounting; the trace rows and physical calculations
remain unchanged.

2026-10-27: Preserved the pre-evidence-retention adapter as
`tools-v3/prefix_trace-pre-preservation.py` (SHA-256
`932b6856936441a20499e9c0828ea1befc15b68374ddeda748390ee4567fb224`). The
active adapter now leaves any output/receipt created before a later failure in
place rather than deleting it. Added a synthetic forced input-change test that
verifies the byte-exact partial prefix remains and no receipt is fabricated.
The shared supervisor README remains untouched.

2026-10-27: Preserved the currently deployed supervisor unchanged and added `tools-v3/supervise_campaign.py` as a separate snapshot. Fixed its nonblocking pipe reader to distinguish `EAGAIN` from a zero-length EOF read; the supervisor now waits for true pipe EOF, a finished leader, and no remaining owned threads before closing and hashing `job.log`. A 100-case synthetic test exercised fast failures writing short and 131,071-byte final tails to alternating stdout/stderr; all exact log contents, SHA-256 values, failed statuses, and exit code 23 were verified. Python syntax compilation passed. Tests used temporary files under `/tmp`; no build, scientific job, SSH action, or deployed supervisor change occurred.

2026-10-27: Preserved the first v3 candidate as `tools-v3-pre-cleanup` and extended only `tools-v3/supervise_campaign.py` to retain output on watchdog cleanup exits. After cleanup, it drains through true EOF; if EOF is not observed within a bounded interval, it records `incomplete_pipe_drain` with byte count and cleanup verification state before hashing the available log. Python syntax compilation and a synthetic timeout test passed: a SIGTERM handler's exact final stderr tail and SHA-256 were retained, the job remained failed for timeout, and cleanup verified all owned identities gone. The earlier 100-case normal-exit regression test remains applicable and passed before this cleanup-only extension. No SSH, build, science job, or source change occurred.

2026-10-27: Preserved the reviewed v3 cleanup candidate as `tools-v3-pre-deadline` before a narrow deadline correction. The post-cleanup drain now checks its deadline on every iteration, including continuous-readable-pipe iterations. A continuous synthetic writer test with a 0.15-second bound returned at the bound, drained bytes, and recorded the `incomplete_pipe_drain` condition; only that test's writer PID was terminated afterward. Re-ran the 100 normal-exit tail cases and timeout/SIGTERM-tail test against the final v3; all passed. Python syntax compilation passed. No remote, build, scientific, or source-tree action occurred.

2026-10-27: Versioned pole gatev2 passedTunya; fourDSfullL3runs completed withidenticalenergy1.2917153256316188e-11m2s andqdotrelativeerror6.25365e-16. Profile02 then stoppedatYSbaselineformat check. Directnohupdiagnostic in newfolder preservedactualstderr: baseline must721rows5sthrough3600s; DSfullFFTtrace has1024rows. Preparingdata-onlybyteprefixadapter withzero-tailHann andexactenergy guards, andseparatenewsupervisor snapshot fixinglostpipe-taildiagnostics. Configsv5 retryprofile03 addsfourtimedbaselineformatjobs andusesprefixonlyforYSkernelinput; allseven scientificbinaries, completeL3checks andphysicalobservable remainunchanged. Priorattempts/logs/data preserved.

2026-10-27: Added an isolated `dspec/source-v4-analytic-pole-control/` source
copy and `dspec/validation-v3-control/` standalone receiver gate to test whether
old v4 scalar Wigner pole rounding explains the L=3 v4/v6 cross-channel
residual. Exactly one helper receives a degree>=1 exact-pole closed-form early
branch; the original Wigner fallback, v4 receiver aggregation, degree-zero
radial callsite, v6 implementation, and fixed numerical comparison threshold
remain unchanged. Twelve-file hashes, exact diff, and provenance are saved.
Generator static checks passed for exact body extraction, call graph, radial
L=0 preservation, complete includes, and no L=0 receiver call. No build,
science, or remote activity occurred; independent review is required before
the parent-owned Tunya control run.

2026-10-27: Added `tools-v4/supervise_campaign.py` as a new snapshot copied from the final v3 supervisor; prior snapshots remain unchanged. Its sole code change removes precreation of parent directories for declared outputs while retaining path validation and rejection of already-existing output files. Python syntax compilation and a synthetic producer test passed: the declared `data/` directory was absent at process start, the producer created it and wrote `result.json`, and the supervisor recorded a successful output hash. No SSH, build, scientific run, or source change occurred.

2026-10-27: Added configs-v6 for the isolated analytic-pole control build and ten low-L3 validation jobs. The plan uses the reviewed output-directory lifecycle fix, retains the six unchanged build-03 scientific binaries, and adds a separately named native v4 analytic control plus independent receiver gate. All scientific jobs are serial and pinned to CPU 0; every plan has thread_cap=50 and numerical library thread counts remain one. This pre-assessment plan deliberately contains no high-degree dispatch. Staged 22 hash-verified files in the fresh Tunya continuation-inputs/v6 directory; all older inputs, failed attempts and scientific evidence remain intact. Fresh Luna static review passed; final Sol review precedes execution. The original strict v4/v6 C/F comparison remains a recorded failure, not reclassified as a pass.

2026-10-27: Added `tools-v4/qualification_checks_control.py` without editing the frozen six-job checker. The adapter imports and reuses that checker for all six original jobs, including exact native/profile parity, full L=3 grids, all 179 radial nodes, finite values, q·P/(2E), metadata, and descriptive cross-code comparisons. It requires a completed independent receiver-gate predecessor, reads the seventh `dsv4-analytic-control-l3` output, and gates control-to-v6 native/profile kernel channels plus trace/spectrum columns at the unchanged 2e-13 tolerance. Original v4-v6 error values and per-channel/column maxima, scales, and locations remain in the report; any original mismatch is marked FAIL, and a passing control yields only `PASS_WITH_ANALYTIC_REFERENCE_CONTROL` with `original_parity_passed: false`. Control mismatch reports FAIL. Python syntax and synthetic seven-job tests passed for both a retained original mismatch with passing control and a control mismatch. All test data stayed under `/tmp`; no SSH, solver run, production source change, or original evidence change occurred.

2026-10-27: Preserved the first control adapter as `tools-v4/qualification_checks_control-pre-predecessor.py` before tightening evidence requirements. The current adapter requires both `receiver-gate` and `analytic-control-receiver-gate` completed before the analytic-control job and verifies each saved log hash against supervisor state. It now checks control-v6 result/domain energy with the frozen checker's unchanged energy comparator, rejects missing physical-option/domain fields, and reports an old-v4/control diagnostic alongside old-v4/v6 and control/v6. `original_parity_passed` is set from the historical comparison regardless of control outcome. Synthetic tests passed for original mismatch plus passing control, control mismatch with the original-failure flag retained, and rejection when the analytic-control receiver gate is absent. Python syntax passed; frozen checker remains unchanged; tests used `/tmp` only.

2026-10-27: Preserved the adapter as `tools-v4/qualification_checks_control-pre-state-verification.py` before closing the provenance gap. The adapter now requires the exact 10 reviewed profile04 job IDs, unique state records, terminal completed status/return code zero/unchanged inputs/verified cleanup for every job, exact required output lists, and verified log/output SHA-256 values from the original state. It hashes the full original state, all outputs/logs, the imported checker and adapter, YSpec baseline and sampling references, and optional radial-node files both before and after assessment; an absent optional YSpec radial file is tracked as absent and may not appear mid-assessment. A present optional file must be byte-identical to the state-bound shared sampling input. Tests passed for clean synthetic data, missing predecessor, failed job state, a pre-assessment output tamper, state mutation during assessment, and optional-file appearance during assessment. No original checker/source, remote path, or scientific result was changed.

2026-10-27: After fresh Luna and final Sol approval, launched the isolated two-variant analytic-control build via nohup in a new session. It completed with unchanged inputs and verified cleanup; observed peak was five threads. New binary hashes: control physical driver 0d05d8a5f9625f41db9373412855fe42fb9bf437c4a599bacec0fc7c6696ca1c, independent control receiver gate 9c0a15cea2a2d9869cb73700cc7e5404cbb9a8194bdcc03f1903bb74efa8f43e. Verified all existing build-03 binary hashes before use. The ten reviewed profile-04 jobs then completed under nohup: both receiver gates, five DSpec full-L3 cases, the energy-preserving baseline prefix, and two YSpec full-L3 cases. YSpec directory creation now succeeds. Successful execution is not yet a parity acceptance claim; the read-only assessment and its provenance checks remain in review. No high-degree dispatch yet, and no full high-degree global kernel has been materialised.

2026-10-27: Read-only profile-04 SHA checks show the analytic-control v4 kernel.bin, trace.csv and spectrum.csv are byte-for-byte identical to v6 native; both DSpec v4 and v6 native/profiled kernel/trace/spectrum pairs are exact, and YSpec native/profiled kernels.bin are exact. The original recurrence reference residual remains preserved. Added configs-v7/qualified-component-profile-plan.json for a separately reviewed read-only profile-04 assessment followed by the existing 25 longer component jobs. These jobs are 64-angle diagnostics at L=128/384/1536, with separate YSpec waveform jobs; they do not materialise a full high-degree global kernel. Fresh review put dispatch on hold until the assessment verifies the completed jobs' entire output/log manifests before and after reading. This provenance-only correction is being implemented in a new adapter snapshot; all numerical thresholds remain unchanged.

2026-10-27: Final provenance adapter checks now also bind YSpec baseline and shared radial-sampling paths to each YSpec job's recorded input hashes, reject any path differing from the reviewed baseline-l3/dsv6-native-l3 inputs, and perform post-assessment provenance verification even when an assessment raises. The absent optional YSpec radial-file case is now explicit in the manifest and any appearance during assessment fails. Repeated synthetic tests passed for all 10 job receipts, clean assessment, historical/control outcomes, failed state, tampered output, mid-assessment state mutation, and optional-file appearance; both adapters compile and the frozen checker remains unchanged. No remote or scientific action occurred.

2026-10-27: Froze the provenance-closed read-only adapter (SHA-256 b835467dd487ae13c650ea69dd88154677a823c8dc7cdd6944e4c561d67c986b), unchanged original checker (3af8a351c66aadd1bcc86a2cc04e303f0cdde27a4aba4b634cede871139901a7), and component plan (76610964e9ff92c6026e449a2977cadf36d8cbe93a1d6416fe984121e2f90aa6) in a three-file hash-verified Tunya continuation-inputs/v7 package. Fresh Luna provenance addendum passed, preserving its earlier HOLD report. Final Sol static review also passed. The original six scientific binaries remain unchanged. The first profile-05 job must qualify completed profile-04 before any high-degree job is released; failure stops the serial queue. Every subsequent high job is exactly the previous reviewed v5 high job with profile-05 task paths. All run plans keep max 50 threads, one numerical thread, and CPU 0 for comparable measurements. Full global high-degree output remains outside these diagnostic measurements.

2026-10-27: Added an isolated `tools-v5` qualification snapshot for the reporting-only fluid-shell correction. `tools-v5/qualification_checks.py` is copied from the unchanged frozen checker and changes only the descriptive cross-code fluid-kappa shell from the mantle range 4200–6200 km to regular outer-core nodes at 2200–3200 km, with matching metric/error labels. `tools-v5/qualification_checks_control.py` copies the final v4 adapter and imports its sibling v5 checker; no parity, energy, finite, full-grid, or provenance gate changed. Static diff confirms only the intended shell/comment/error/label edits in the copied checker and the import path change in the adapter. Synthetic outer-core fluid-node selection and the existing synthetic provenance suite passed; no SSH or scientific runs were made.

2026-10-27 09:11 UTC: Corrected only the versioned descriptive fluid comparison to regular outer-core radii 2200-3200 km (22 actual saved native fluid radial sides), with fresh Luna and final Sol approval. Original checker and failed profile-05 remain unchanged. The new tools-v5 adapter differs from its reviewed predecessor only in loading its sibling checker; all numerical acceptance thresholds remain fixed. Staged hash-verified continuation-inputs/v8 and launched profile-06 under nohup, PID/session 2045526. After a separate SSH connection, verified PPID=1, session=2045526, status running. The complete read-only qualification passed as PASS_WITH_ANALYTIC_REFERENCE_CONTROL: original recurrence comparison still FAIL, analytic-control/v6 parity PASS, exact native/profile pairs, full L3/179 radial sides, qdot, finiteness and complete before/after evidence manifests all passed. First L128 numerical job is actively running; observed campaign threads=2 and thread_cap=50. All 25 component jobs are now locally controlled by the remote supervisor and need no laptop process or connection. Native/profile high-degree kernels are 64-angle diagnostics; full high-degree global materialisation remains outstanding. Copying the full small continuation audit (sources, tools, configs, all reviews and receipts) to a fresh remote notes directory for sleep readiness.

2026-10-27 09:13 UTC: SLEEP READY. Verified 238 continuation audit files on Tunya; all 5418 original campaign/dependency files were already verified. Remote profile-06 remains running under nohup, PID/session 2045526, PPID 1, max 50 threads, one numerical thread. Qualification completed, native L128 completed and instrumented L128 is progressing. No laptop process or SSH connection is required. Final sleep receipt saved both locally and in remote continuation/sleep-ready-20261027T0913.json. Full high-degree global materialisation remains outstanding.

2026-10-27 09:47 UTC continuation: Tunya SSH reconnected; profile-06 survived laptop disconnection and remains within cap 50 (two observed campaign threads). All L128/L384 native and profiled cases completed; L1536 started. Copied 100 files from completed jobs into receipts/profile06-snapshot-20261027T094727Z, verifying every copied log/output against the original supervisor SHA before transfer and again locally. No running output was copied or changed. Native L384 wall times: v4 295.4885 s, v6 122.8253 s, YSpec kernel diagnostic 321.0469 s plus its separately measured waveform 85.6767 s. These are 64-angle/179-side diagnostic runs, not full global grids. Raw v6 L384 profile: baseline forward 31.9520 s inclusive; kernel workspace 34.4556 s inclusive; field sampling 55.8778 s inclusive; final pointwise contraction 0.0094013 s. Exclusive packing is 21.3024 s; radial coefficient recovery 15.4865 s; angular GEMM 16.0909 s. No overlapping inclusive scopes are added to an exclusive budget. Factorization count remains 80,745 = 105*(2*384+1), reflecting 53 baseline plus 52 kernel passes. A data-only collection/finalisation tool is being prepared for independent review; production numerical sources remain unchanged and no new algorithmic phase has begun.

2026-10-27: Added `tools-v6/collect_component_profile.py` and `tools-v6/finalize_profile.py` as isolated data-only helpers; prior tool snapshots are unchanged. The collector validates every copied state/log/output digest against the snapshot receipt and completed-job state, scans numeric JSON/CSV/binary outputs for non-finite values, checks byte-exact native/profiled kernel/trace/spectrum pairs, retains full raw profiler stage/counter records, and reports exclusive-stage root closure and userspace telemetry. It reports first-64-angle cross-code channel differences and shared-prefix waveform L2/energy comparisons descriptively, preserving explicit limits on angle coverage, degree convergence, and single-pass timing. The finalizer waits for the existing `profile-06` queue with a 12-hour default bound and at most 30-second polling, snapshots only completed-job evidence with verified hashes, and invokes the collector with numerical-library thread variables set to one; it never launches or changes campaign jobs. Synthetic provenance/tamper, profile-closure, and finalizer CLI checks passed. The actual verified L128/L384 snapshot produced 12 exact parity passes, eight profiler records with exclusive/root closure within 1e-6 s, and nine explicit pending L1536 jobs. Cross-code waveform relative energy differences were 2.62e-6 (L128) and 7.45e-6 (L384), recorded descriptively against the existing 1e-3 physical scale with no new threshold. No SSH, scientific run, build, production source edit, or remote finalization occurred; fresh Luna and Sol reviews remain required before remote use.

2026-10-27: Tightened the v6 collector/finalizer before review: profile identity and objective are checked against the DSpec, YSpec kernel-probe, and YSpec waveform scopes; the waveform relative energy denominator is the saved DSpec reference energy; the snapshot includes the state hash; and the finalizer rejects path traversal, duplicate destination writes, incomplete output hash maps, or changed state during copying. It logs the terminal snapshot/collection phase and preserves verified completed-job evidence on terminal queue failure. Re-ran Python syntax compilation, synthetic provenance/tamper/profile/finalizer CLI checks, and collection of the saved profile06 snapshot; all passed with the same 12 completed-pair exact parity passes and nine L1536 jobs pending. No remote action occurred.

2026-10-27: Preserved the pre-review tools as `tools-v6-pre-parent-review/` before applying the parent review corrections in `tools-v6/`. The expected plan is now the exact 26 IDs: DSpec v4 only at L128/L384, DSpec v6 and YSpec native/profiled jobs at all three L values, with baseline jobs at each L. Profile accounting requires the unique `physical_case_total` root, nonnegative stage times/calls, and the existing `max(1e-8 s, 2e-12*root)` closure tolerance. Waveform checks verify 5-second grids, 721/1024 sample extents, independent DSpec/YSpec Hann weights within 1e-12, zero DSpec weights beyond the YSpec prefix, Hann formula, reference-normalized weighted L2 and YSpec-referenced relative energy at the existing 1e-3 gates, and consistency with both saved energy fields. Finalization requires all exact 26 jobs, verified cleanup/input invariance, observed thread counts, and a passing complete summary; failed queues still write a failure manifest and preserve completed-job snapshots. Summary explicitly scopes contacts-off volume-only kernels, omitted separate source-shell contacts, checks DSpec qdot/(2E) <=1e-10, reports actual degree from result metadata, and binds YSpec radial sampling to the exact state-hashed DSpec file without interpolation. Synthetic tests and Python syntax pass. The newer state-bound snapshot `profile06-snapshot-20261027T095752Z` yields status INCOMPLETE with six L1536 jobs pending, eight completed exact parity passes, eight profiler root closures within tolerance, nine completed DSpec qdot gates passing, exact L128/L384 radial input identities, and both shared-prefix waveform gates passing. No remote execution or scientific run occurred; fresh Luna and Sol review remains required.

2026-10-27: Copied all three v6 tools into a new `tools-v7/` snapshot, leaving v6 and its pre-review copy unchanged. Changed only the descriptive cross-code kernel comparison to respect material channels: compare A/C/F/L/N on solid radial sides and DSpec kappa against YSpec A+C+F on fluid sides; separately report regular solid 4200–6200 km and outer-core 2200–3200 km. Each maximum includes absolute/scaled difference and its angle, radial side, element, node, radius, fluid flag, and values. The exact 179-side radial table is hash-bound through the YSpec sampling input and shapes each angle-major kernel; reports make no kernel tolerance or pass claim. A synthetic material-mapping test confirms fluid kappa matches the YSpec sum despite DSpec fluid A–N zeros, while a deliberately wrong solid A pairing remains detected. Python syntax, all synthetic checks, and collection from the state-bound 100543 snapshot passed. That snapshot reports status INCOMPLETE with four YSpec L1536 jobs pending; L128/L384 material-aware metrics were produced from 64 angles x 179 radial sides. No SSH or scientific run occurred; v7 needs its own fresh review before use.

2026-10-27 10:22 UTC: Completed immutable snapshots at 09:57:52, 10:05:43 and 10:22:20, each binding the captured state plus all copied completed-job artifacts to SHA-256; earlier snapshots remain untouched. DSpec v6 L1536 native/profiled kernel, waveform and spectrum are byte-exact. Native wall 536.4424 s, instrumented wall 533.0479 s, native peak 1.0914 GiB. Instrumented nonoverlapping top-level budgets: baseline 115.5141 s, kernel frequency workspace 135.0675 s, tensor sampling 281.4451 s, pointwise contraction 0.009594 s. Exclusive coefficient packing is 134.8233 s, angular GEMM 64.6127 s, radial coefficient recovery 64.5504 s. Sparse factorization count is 322,665 = 105*(2*1536+1); current angular/radial cache ceilings cause misses at this degree. The independent YSpec L1536 waveform completed in 1061.7258 s. Relative windowed waveform L2 is 4.6411e-5 and energy difference 2.6797e-5, both below the unchanged 1e-3 cross-code gates. This is agreement at fixed L, not proof of degree-tail or pointwise kernel convergence.

2026-10-27: Parent inspection caught a material-convention problem in the first descriptive report: fluid DSpec A/C/F/L/N are zero with the fluid response stored in kappa, whereas YSpec uses parameter components whose fluid bulk response is A+C+F. Preserved tools-v6, its first partial report, and the Sol HOLD; staged input v9 was never executed and is annotated DO_NOT_RUN on Tunya. New tools-v7 changes only descriptive comparison/matching metadata and tests: solid A/C/F/L/N, fluid kappa versus YSpec A+C+F, plus explicit regular solid/outer-core subregions, scales, absolute errors and one-sided node locations. No physical threshold was relaxed. Fresh Luna and final Sol reviews and data-only synthetic tests passed. Source/runtime hashes are recorded in component-analysis-inputs-v10.json. The updated partial material-aware report is analyses/profile06-partial-material-aware-20261027T102220Z.json, status INCOMPLETE with three remaining YSpec jobs. Solid kernel discrepancies remain unresolved, notably C/F at the receiver interface; waveform agreement is not global kernel validation.

2026-10-27 10:20-10:22 UTC: After both reviews, launched the one-off analysis finalizer under nohup on Tunya, PID/session 2120868, CPU 1, one numerical thread, max wait 43,200 s. Separate SSH verified PPID 1 and active waiting. It never changes or launches physics jobs: it waits for profile-06, then copies a fresh verified terminal snapshot and writes analysis/profile06-final-v1/component-profile-summary.json plus a receipt. Failure/incomplete outcomes remain explicit. Before launch, verified live I8 plan SHA against the approved 1ecfc12a... digest, the same digest in running state, all 23 saved command/cwd records against that plan, and all deployed analysis source hashes. Main queue continues unchanged with cap 50 and two observed campaign threads; the waiting analysis controller adds one thread. No production numerical code, branches, commits or original evidence changed.

Next numerical phase remains unstarted. Current evidence prioritises reuse of the duplicate forward pass, fixed coefficient packing/index layouts, batched native-node recovery, and a frequency-accumulating global GSH path that writes each output tile once. At L1536 the six-channel full volume is 151.474 GiB; the current 52 file read/write passes imply 15.384 TiB logical accumulation traffic. All-theta five-spin real Wigner payload is about 270.44 GiB before other arrays; Tunya currently has about 1468 GiB available, so large-cache and bounded-tile alternatives should be measured rather than inheriting laptop cache limits. These are execution-plan estimates, not implemented speedups or authorisation for another phase. Full high-degree global materialisation and its physical kernel-convergence proof remain outstanding.

2026-10-27 12:31 UTC: The profile-06 queue completed all 26 jobs, with no failed jobs, unchanged recorded inputs, verified cleanup, and at most two observed campaign threads (configured cap 50). The detached finalizer completed at 12:31:23 UTC; remote summary SHA-256 is ecb50b6892fea55d723c76c86976ba4f956b56a595d2d02d89059194634c046f. Retrieved a fresh immutable snapshot, receipts/profile06-snapshot-20261027T123112Z, with 134 state/log/output digests verified remotely and locally. The unchanged reviewed tools-v7 collector produced analyses/profile06-complete-material-aware-20261027T123112Z.json, status PASS for the scoped campaign gates: all expected jobs, byte-exact native/profile output pairs, exclusive profile closure, qdot/(2E), radial sampling identity, and fixed-L independent waveform agreement. This PASS does not certify pointwise kernel agreement or global convergence. No scientific source change, build, commit, push, or further numerical campaign occurred.

Completed native timings at L128/L384/L1536: DSpec v6 (including its baseline) 41.8656/122.8253/536.4424 s; YSpec independent waveform 16.8427/85.6767/1061.7258 s; YSpec kernel diagnostic using the shared DSpec baseline 77.1774/321.0469/3338.3542 s. At L1536 the instrumented YSpec kernel root is 3356.5473 s, including exclusive upward integration 1005.2696 s, downward integration 57.9293 s, spheroidal b-vector transport 2195.8723 s, strain recovery 11.7059 s, source/adjoint angular synthesis 37.4087/22.4493 s, elastic/objective contractions 0.08135/0.09870 s. All raw exclusive and inclusive stages are retained in the complete JSON. These are single numerical thread comparisons at configured SEM discretisation and YSpec integration tolerance, not equal local-kernel-accuracy comparisons or full global throughput measurements.

The final L1536 material-aware kernel comparison remains descriptive and unresolved: max absolute difference divided by the larger channel maximum over all solid sample sides is A 2.493%, C 42.781%, F 34.069%, L 15.572%, N 2.700%; fluid kappa versus YSpec A+C+F is 8.893% over all fluid sides, with its largest absolute error in the shallow ocean. The regular outer-core region is 0.02165%; regular-solid region errors range from 4.79% to 14.68%. These are first-64-angle diagnostics, with small absolute amplitudes recorded separately; the sampled angular positions differ with L. DSpec energy changes by 2.927% between L384 and L1536 relative to L1536, so fixed-L cross-code agreement is not a degree-convergence certificate. Proposed next phase remains unstarted: establish representative local field/kernel correctness, replace repeated map/set packing with a fixed contiguous layout, reuse frequency responses/factors across baseline and kernel passes, batch native-nodal recovery, then validate the frequency-accumulating global GSH execution path and materialise the full requested global volume under the existing evidence and resource constraints.

2026-10-27 plan update requested by Alex: Added acceleration-and-figures-plan.md. The proposed next implementation phase now includes fixed-layout packing, shared baseline/kernel frequency responses, batched native-node recovery, reusable angular/global-GSH plans, tile-local frequency accumulation, compatible symbolic reuse and controlled parallelism. Receiver lifts/enrichment, receiver-point convergence and Sobolev-gradient work are explicitly deferred; existing source/interface semantics remain preservation gates. A bounded away-from-endpoint validation accompanies the accelerations, followed by a separate approved figure/global-volume phase with early sections, transparent display masks and complete-bandwidth sampling. Full global materialisation remains required. This is a documentation-only update: no numerical source edit, new run, build, commit or remote change occurred; implementation awaits approval under the staged-change instruction.
