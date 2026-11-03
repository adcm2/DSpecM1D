# Sol delta review: Tunya build path retry

27 October 2026. Static review only; no SSH, build, scientific execution or source edits. This supplements `tunya-launch-sol.md` and leaves that review's numerical-plan HOLD in place.

**Build-only retry: PASS for the reviewed path delta, conditional on parent-owned transfer/host verification and completion of the moved Eigen header comparison.** The unfinished numerical checker does not block this separate build-only task. No numerical launch or numerical result is accepted here.

Compared preserved `prepare_configs-v1.py` with active `prepare_configs.py`, and parsed all three old/new config JSON files. Each new JSON is exactly the old structure after only these four substitutions:

- continuation input root `continuation-inputs/v1` to `continuation-inputs/v2`;
- fresh build root `kernel-profile-20261027` to `kernel-profile-20261027-02`;
- Eigen include root `/space/adcm2/CCDSM/deps/toolchain/eigen/include/eigen3` to `/space/adcm2/research/projects/CCDSM/deps/toolchain/eigen/include/eigen3`;
- analysis Python from the former CCDSM virtual environment to `/usr/bin/python3`.

The generator additionally writes local output into `configs-v2` instead of `configs`; no compiler flag, source selection, variant count, thread/resource cap, job order or numerical argument changed. Seven direct sequential compiler variants and one build-plan job with `thread_cap=50` remain. The build-plan interpreter was already `/usr/bin/python3`; the analysis-Python substitution affects only the future checker job.

The remote package must install the **contents** of local `configs-v2` as `continuation-inputs/v2/configs/`: the reviewed build plan intentionally invokes `continuation-inputs/v2/configs/build.json`, not a remote `configs-v2/build.json`. Confirm that mapping and its hashes before detached launch. Use a new supervisor task ID and preserve the failed build-01 task and original partially populated build root; do not resume the failed task or repurpose its output directory.

Reviewed local config digests:

```text
06828df1795216e524c8d97efd8b432d1227aebabd86211f6034778167890720  configs-v2/build.json
9bc3ad0c3538a99b5c78e7789bc556ae8660f7fab9852ddb90641063337da52c  configs-v2/build-plan.json
7e49a8f47e291f248668cccc9417236100813532932ed573a299a1b13a2c5ddc  configs-v2/qualification-and-profile-plan.json
```

The parent reports build-01 stopped during source fingerprinting before compiler execution, failed evidence is preserved, and system Python has NumPy 1.26.4. Those are parent-provided remote facts, not independently reverified by this static review. Moved Eigen content identity remains a parent-owned prerequisite; a path's existence alone does not establish byte identity. The builder will still fingerprint source/dependency trees before and after compilation.

Launch only the build plan under the already reviewed detached supervisor after those checks. Retain userspace-watchdog accounting limitations and the 50-thread maximum. The numerical plan remains on HOLD until the final checker edits pass fresh Luna and Sol review; the path-only PASS does not override that gate.
