# Fresh review: global GSH thread checker

**Result: PASS for the checker and its synthetic receipt tests.** `py_compile` passed and all seven `test_check_global_threads.py` tests passed. No build, SSH, or science was run.

The wrapper delegates numerical comparison to the existing `check_global_cache.assess`. That retains the `2e-13` relative comparison gate and, with `exact=True`, its byte-exact kernel/trace/spectrum and exact energy/qdot requirements. Receipt verification, output hashes, and unchanged-state checks remain in the delegated path. The fixture covers a complete synthetic receipt comparison, including an exact pass and a deliberate exact-kernel failure.

The checked field names and policy strings agree with the combined-global-v6 producer: the result and manifest carry requested GSH threads, execution policy, OpenMP dynamic-team state, requested/effective worker counts and products, and the 48-thread cap. The checker cross-checks these against the receipted `--global-gsh-threads N` argument, validates `workers * N <= 48`, and requires the `FFTWpp::Estimate` planner label. The legacy exception is explicit and serial-only: it permits a reference with no GSH flag or new GSH metadata only when one GSH thread was expected and the manifest says `FFTWpp::Estimate`; the candidate has no legacy exception. The source uses `Execution::Sequential()` at N=1 and `Execution::Parallel(N)` otherwise, matching the checker policy labels.

Observed process peaks are receipt fields reported as observations; missing samples remain null and have no minimum threshold. The checker does not validate `OMP_THREAD_LIMIT` or other inherited environment values. The planned supervisor limit of 1 and per-job inline override therefore need verification in the run plan/receipt path; they are outside this checker’s current assertions. V6 itself also lacks a runtime `omp_get_thread_limit()` check, so this PASS is not a runtime guarantee that the requested OpenMP team actually formed. No scientific equality or performance result is established here.

## SHA-256

- `tools/check_global_threads.py`: `011845c8e2bdf592eb08296306af5227ea519dffd6960854db9c8623faaf3271`
- `tools/test_check_global_threads.py`: `c5cf7037ad8eac0d528af080f206d8cd1d3913a295b9f4f8bdbe6da8b62953b3`
- `tools/check_global_cache.py`: `9890e033bfab6f9bd7a4c3bfdf42cab005fdddb97fa248cbd107d0d82c1a21af`
- `tools/check_global_pair.py`: `340efdf0c60b80eea23fa69d642595cf7b989112d9310ed792b21ebe4b36d51a`
- `tools/check_pair.py`: `6f4d72e7da0489b6b3d0b6fd1d5c15d7ba3b7e6179a14ef453b7e8ccb762e646`
- `candidates/combined-global-v6/global/global_kernel_consumer.hpp`: `1a1600999680edf5f504a37471c80a66eca60db1727af9c4a7512b63ab3c8632`
- `candidates/combined-global-v6/physical_driver.cpp`: `446aec2ba2950bcde89ab8c75f820cef2ec65180f9f761d00258e72c5db2c061`
