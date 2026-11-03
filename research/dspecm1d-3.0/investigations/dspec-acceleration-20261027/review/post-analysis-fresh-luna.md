# Fresh static review: post-analysis queue and summary helper

**Result: PASS for the bounded post-analysis plan and summary helper.** No queue was launched and no science/build/SSH was run.

The plan contains 10 jobs under a 4-thread cap. Its first job waits up to eight hours for the global-thread queue to reach a terminal state, and requires every recorded job cleanup to be verified before writing a dependency receipt. The state argument is now the dependency task directory, with `state.json` constructed inside the wait script; this avoids treating the mutable state file as an argv-file input. The wait job also declares the immutable staged plan itself as a nonempty input, satisfying the supervisor’s input requirement. The supervisor processes jobs in plan order, so this wait is a gate ahead of the paired reruns.

The four paired reruns are v4/v5 L128/L384 64-angle kernel probes on one pinned CPU with matched physical settings. Each has a corresponding `check_pair.py --exact` comparison against the existing symbolic direct reference. The aggregate job reads the intended native-v5, global-v8, and preserved serial-global snapshots; its resolution arguments use the `accel-layout-run-01` layout reference and the preserved serial-global L128/L384 jobs. Its output root is the fresh `accel-post-analysis-01` job directory. It records the global and resolution statuses but allows a descriptive helper failure to be summarized; the summary states that bookkeeping is not scientific acceptance and that the original-global gate remains failed.

`check_phase_a_summary.py` is byte-identical to `summarize_snapshots.py`. It reports saved state/job metadata and checks hashes for existing reported artifacts without rerunning the solvers. Missing artifacts are not upgraded into a pass or acceptance claim by the aggregator.

## SHA-256

- `configs/post-analysis-plan.json`: `d54a68ed32389fa6430f9b43b23ce892df7793e9b2b4557e35601ed4e3ed53b2`
- `tools/check_phase_a_summary.py`: `b042a7e8ca4a3fd95a02f25676723b28f6d8cfdd8bd639e4162e3c0b1d88eefa`
- `tools/summarize_snapshots.py`: `b042a7e8ca4a3fd95a02f25676723b28f6d8cfdd8bd639e4162e3c0b1d88eefa`
- `tools/check_global_resolution.py`: `025e52a38b66bf4aec365d6e30ac4bc35260545f322997c5f964c3c8ce684020`
