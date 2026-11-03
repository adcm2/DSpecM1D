# Global GSH inner-thread checker status

## Change log

- Added `check_global_threads.py` as a data-only wrapper around `check_global_cache.assess()`, preserving its `2e-13` comparisons and optional exact kernel/trace/spectrum/energy gates.
- Added explicit reference/candidate GSH thread counts and validated the receipt command, result, and manifest fields against the requested counts, GSH execution policy, configured/effective outer workers, requested/effective products, and the producer's 48-thread cap. Keys match the inspected combined-global-v6 producer schema.
- Added an explicit `--reference-legacy-serial` path, valid only for expected reference GSH count 1 when the command flag and all new thread fields are absent and the manifest declares `FFTWpp::Estimate`. The candidate has no legacy fallback.
- Added supervisor-observed thread/RSS/time peaks to reports when available; missing samples are reported as null, without a minimum threshold for short jobs.
- Added synthetic receipt integration for current GSH metadata and legacy serial references, plus wrong CLI/count/product, invalid legacy mode, missing peaks, and exact kernel-regression checks.

## Validation

- `python3 -m py_compile check_global_threads.py test_check_global_threads.py`: passed.
- `python3 -m unittest -v test_check_global_threads.py`: 7 tests passed.
- No build, SSH, or scientific execution was performed.
