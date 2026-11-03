# Final short Sol review: prefix adapter, logger v3 and profile-03

27 October 2026. Static read-only review after `pipeline-v3-fresh-luna.md` PASS, including its final preservation-only follow-up. Inspected current adapter, supervisor diff, actual frozen YSpec baseline reader and configs-v5 wiring. No source edits, builds, solver/scientific runs or SSH performed. Reported synthetic test results were not independently repeated.

**PASS for the reviewed pipeline correction and fresh profile-03 launch contract.** This is source/config acceptance, not a successful numerical-gate result. Preserve profile-01/profile-02, their failed evidence and all prior scripts.

Reviewed digests:

```text
ebf8a91282e3653c258415369e8a1241b721d900b548721387c2a5399837baf5  tools-v3/prefix_trace.py
3961cb872581217537a3d7ccf5d4e16c3500f7ea7ecb73cb122abecf84b11d7f  tools-v3/supervise_campaign.py
3e2db2428dee03b2a75190a29441cc82c7931c5eea01c5bf927cb2563f5f586e  configs-v5/qualification-and-profile-plan.json
3af8a351c66aadd1bcc86a2cc04e303f0cdde27a4aba4b634cede871139901a7  tools/qualification_checks.py
```

The actual YSpec `readTrace` requires 721 rows at 5 s spacing through 3600 s, compares input weights with its Hann window at `1e-12`, and computes `0.5*dt*sum(w*u*u)`. The adapter explicitly validates the full DSpec 1024-row trace on 0–5115 s, finite values, exact 5 s grid, matching Hann1000–1200 weights, and exact zero weights for every removed row. It writes the original header and first 721 data rows without numerical reformatting. Every retained windowed sample is therefore unchanged; every removed contribution has zero weight. Positive finite full/prefix energies must be exactly equal under its common `math.fsum` calculation, in m²·s. That equality establishes adapter energy preservation; it does not claim bitwise equality between Python summation and either solver's summation order.

Input hashes before/after, prefix hash, row/grid/window metadata and energy evidence are saved to a fresh receipt. Exclusive creation refuses existing outputs and path collisions. A failure after writing preserves the partial prefix/receipt for audit rather than deleting evidence. The input remains read-only. The supervisor also hashes declared inputs and outputs, so failure and any detected source change stop progression.

The supervisor change correctly distinguishes `BlockingIOError`/EAGAIN from a zero-length read reporting EOF. A fast process exit no longer allows EAGAIN to masquerade as completed logging. The normal loop continues sampling resource/time limits while draining. After cleanup it drains final output with a deadline checked at the top of every iteration, including continuously readable output; failure to observe EOF is recorded explicitly. The deadline covers the polling/read-loop logic, not a kernel guarantee against arbitrary blocked filesystem I/O. Existing identity-based cleanup, subreaper handling, pidfd signaling, cleanup verification, hash-bound resume, lock, thread/RSS/timeout checks and failed/unclean handling remain unchanged. The reported 100 fast-tail cases, timeout/TERM stderr and continuous-writer synthetic checks support the intended cases but were not rerun in this review.

Config comparison found exactly four added adapter jobs: baseline-l3/l128/l384/l1536. Each hashes the adapter and its matching DSpec v6 native trace, has 60 s/1 GiB limits, and declares both prefix CSV and receipt. Each precedes the corresponding YSpec kernel consumers. YSpec kernel commands/inputs use the prefix; matching radial-node files remain unchanged. YSpec waveform commands are unchanged. The qualifier additionally hashes the L3 adapter CSV/receipt; its command and strict acceptance code are unchanged.

All seven C++ executable paths remain the reviewed existing receiver-v2 plus six solver binaries. Existing DSpec jobs match v4 after only literal profile-02→profile-03 prior-result paths. The plan still runs receiver, complete L3 jobs, strict qualifier, then 64-angle high-degree diagnostics; it has 34 serial jobs and cap 50. Numerical executables retain CPU-0 affinity and numerical library thread limits of one. Adapter overhead is separately timed. No changed baseline file is represented as a recomputed waveform or a different scientific observable.

Parent may launch a fresh profile-03 with the reviewed new logger after verifying these exact transferred script/plan/checker hashes and existing binary receipts, confirming runtime/CPU/disk readiness, and exclusively creating the detached log. Use the new supervisor directly; old task directories remain untouched. Retain the userspace-watchdog caveat: polling is not kernel aggregate accounting and may miss short-lived processes/threads. Numerical and performance conclusions still require the actual gates and controlled saved results; the parser/logging corrections alone prove neither.
