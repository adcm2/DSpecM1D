# Final direct native acceleration evidence review

27 October 2026. **PASS for the bounded completed direct-route evidence.** No code edits, SSH, build, or scientific execution were performed. This review does not authorize production integration or promote the experimental global consumer.

## Evidence and provenance

- Before source: `analyses/layout-all-components.json`, backed by `receipts/layout-completed-20261027`.
- Final source: `receipts/final-direct-completed-20261027`, state SHA256 `d512dd41511c4212913cc355fae700a5e42d6984df03045dcb509e2b4b1571f2`; the snapshot verification manifest matches that state hash.
- All 28 jobs completed with return code 0, empty failure reasons, unchanged inputs, and verified cleanup. All 13 comparison gates report PASS. All recorded `byte_exact_files` entries are true, and every reported numeric maximum difference is zero. Reuse gates allow tolerance but the observed result is byte exact, not merely within tolerance.
- Independently rehashed all 155 locally collected files in the 170-entry final artifact manifest: no mismatches. The 15 kernel binary payloads were not collected locally. Their byte equality is established by the hash-verified checker receipts and remote output-hash manifests; this review did not reread those binary payloads locally.
- Final commands use accepted `kernel-acceleration-20261027-final-04/final-profiled/physical-driver` (SHA256 `1be78fc2aae6c1c41b1e14871396d02c236cac2e2a11a2bb4e5302595d0b1075`) and native binary SHA256 `8ccfc45655f73be1a39f0c55e8ca50067c3c570678d1806788ee209136540cbf`. Config `final-v4-build.json` stages `combined-final-v4`; its existing fresh Luna and short Sol source/build reviews provide the source-review chain. This review checks the collected run evidence rather than repeating that source review.
- Every final result requests and uses direct response reuse, reports zero basis-reconstruction error, no capacity/preflight fallback, `response_record_spool_used=false`, and `global_consumer_full_transforms_performed=false`. This is the direct native candidate, not the experimental global route.

## Workload and timing distinctions

The old production baseline used `--forward-only on`. The new matched controls and all final direct candidates used `--forward-only off`; the response capture pass replaces that route. Consequently the old-production comparison answers practical before/after runtime, whereas the matched controls isolate reuse under matched flags. Do not label the old-production timings as matched controls. All are single observations from completed process telemetry, not repeated-run distributions.

| Degree | Old production process s | Matched control process s | Final reuse, symbolic off s | Final reuse, symbolic on s | Old→final speedup | Matched→final speedup |
|---:|---:|---:|---:|---:|---:|---:|
| 3 | 2.203356 | 2.072617 | 1.579783 | 1.601514 | 1.376x | 1.294x |
| 12 | 30.078734 | 29.959697 | 19.778072 | 19.517453 | 1.541x | 1.535x |
| 128 | 41.856418 | 44.221491 | 32.421846 | 27.113679 | 1.544x | 1.631x |
| 384 | 123.951199 | 130.525318 | 85.575150 | 80.840457 | 1.533x | 1.615x |
| 1536 | 532.983888 | 559.869149 | 385.358219 | 370.683470 | 1.438x | 1.510x |

In the requested descending order, old production 1536/384/128 timings are **532.983888 / 123.951199 / 41.856418 s**, versus final symbolic **370.683470 / 80.840457 / 27.113679 s**. The actual full-precision old numbers are retained in the accompanying data JSON; the table is authoritative if a rounded narrative value differs. Matched controls are **559.869149 / 130.525318 / 44.221491 s**.

Native/profiled parity was checked at degrees 3, 12, and 1536. At 1536 native took 369.676990 s versus 370.683470 s profiled (about 0.272% longer profiled). That native/profile gate explicitly required byte equality and passed for trace, spectrum, radial metadata, angular metadata, and the six-component kernel probe, with zero numeric differences.

The physical case has 30 elements, six GLL nodes per element, 180 native element/node sides and 179 available selected kernel sides; 52 kernel frequency bins. Degrees 3 and 12 cover the complete target angular grids (91 and 1225 angles per trace; 16,289 and 219,275 kernel points). Degrees 128, 384, 1536 cover all 179 selected radial sides at 64 diagnostic angles, or 11,456 points each. Those high-degree kernels are angular probes, not complete global outputs or full-Earth volume integrals. Source contacts are disabled in all these runs; no contact-kernel equivalence conclusion follows from them.

## Scientific hash transitivity

For every tested degree 3/12/128/384/1536, old-production reference → matched reference → final symbolic candidate has identical SHA256 for `trace.csv`, `spectrum.csv`, the six-component kernel binary/probe, radial nodes, kernel radial sides, and angular grid. The five non-kernel artifacts were also reread and compared locally. The native 3/12/1536 gates extend the same equality to uninstrumented binaries. Full hashes and local-reread flags are in `native-final-evidence-data.json`.

Result/configuration metadata remains allowed to differ for route flags, cache statistics, timing, and paths. In particular old forward-only=true versus matched/final=false is a documented execution-route difference; byte-identical scientific outputs support the timing comparison without pretending the old flags matched.

## Component attribution

The following entries are **exclusive** stage seconds. They can be compared without adding nested inclusive timers. Rows identify existing timer scopes rather than claiming that a scope contains only the noun in its label.

| Degree | Exclusive component | Old production s | Matched control s | Final symbolic s |
|---:|---|---:|---:|---:|
| 3 | `tensor_coefficient_packing_and_key_reorder` | 0.203912 | 0.203588 | 0.026040 |
| 3 | `sparse_numerical_factorization` | 0.353081 | 0.358686 | 0.166184 |
| 3 | `sparse_symbolic_analysis` | 0.058812 | 0.059861 | 0.000275 |
| 3 | `tensor_angular_gemm_complex_and_packing` | 0.237787 | 0.240490 | 0.235051 |
| 3 | `radial_interpolation_and_strain_numerator_build` | 0.125596 | 0.127636 | 0.131180 |
| 3 | `spheroidal_source_lift_projection_and_solver` | 0.025661 | 0.045554 | 0.060989 |
| 3 | `dynamic_H_assembly` | 0.084406 | 0.079330 | 0.062569 |
| 3 | `sparse_forward_solve_including_rhs_scaling` | 0.061281 | 0.061254 | 0.032752 |
| 3 | `sparse_adjoint_solve_including_rhs_scaling` | 0.020319 | 0.039212 | 0.020958 |
| 12 | `tensor_coefficient_packing_and_key_reorder` | 10.140237 | 10.065260 | 0.588100 |
| 12 | `sparse_numerical_factorization` | 1.135043 | 1.046628 | 0.565420 |
| 12 | `sparse_symbolic_analysis` | 0.211005 | 0.194309 | 0.000277 |
| 12 | `tensor_angular_gemm_complex_and_packing` | 10.137473 | 10.143606 | 9.997484 |
| 12 | `radial_interpolation_and_strain_numerator_build` | 0.492820 | 0.497927 | 0.493984 |
| 12 | `spheroidal_source_lift_projection_and_solver` | 0.102411 | 0.165430 | 0.236737 |
| 12 | `dynamic_H_assembly` | 0.302501 | 0.251276 | 0.192337 |
| 12 | `sparse_forward_solve_including_rhs_scaling` | 0.234871 | 0.234308 | 0.119706 |
| 12 | `sparse_adjoint_solve_including_rhs_scaling` | 0.076507 | 0.152615 | 0.077787 |
| 128 | `tensor_coefficient_packing_and_key_reorder` | 6.445964 | 6.471704 | 0.215656 |
| 128 | `sparse_numerical_factorization` | 10.086872 | 10.087384 | 4.760201 |
| 128 | `sparse_symbolic_analysis` | 1.994701 | 1.989403 | 0.000280 |
| 128 | `tensor_angular_gemm_complex_and_packing` | 5.420350 | 5.434530 | 5.314871 |
| 128 | `radial_interpolation_and_strain_numerator_build` | 5.102765 | 5.157167 | 5.229429 |
| 128 | `spheroidal_source_lift_projection_and_solver` | 1.071299 | 1.783586 | 2.576910 |
| 128 | `dynamic_H_assembly` | 2.729018 | 2.666085 | 1.337878 |
| 128 | `sparse_forward_solve_including_rhs_scaling` | 2.409808 | 2.413757 | 1.208716 |
| 128 | `sparse_adjoint_solve_including_rhs_scaling` | 0.780368 | 1.577901 | 0.792225 |
| 384 | `tensor_coefficient_packing_and_key_reorder` | 21.337348 | 21.444833 | 0.654231 |
| 384 | `sparse_numerical_factorization` | 28.606972 | 27.885943 | 13.299868 |
| 384 | `sparse_symbolic_analysis` | 5.859659 | 5.906983 | 0.000276 |
| 384 | `tensor_angular_gemm_complex_and_packing` | 16.105384 | 16.105166 | 15.979531 |
| 384 | `radial_interpolation_and_strain_numerator_build` | 15.554709 | 15.818464 | 15.744014 |
| 384 | `spheroidal_source_lift_projection_and_solver` | 3.197237 | 5.341091 | 9.432723 |
| 384 | `dynamic_H_assembly` | 8.005265 | 7.915390 | 3.982850 |
| 384 | `sparse_forward_solve_including_rhs_scaling` | 7.037477 | 7.094295 | 3.563593 |
| 384 | `sparse_adjoint_solve_including_rhs_scaling` | 2.276626 | 4.613896 | 2.329357 |
| 1536 | `tensor_coefficient_packing_and_key_reorder` | 134.739959 | 134.936029 | 2.606805 |
| 1536 | `sparse_numerical_factorization` | 98.255247 | 101.432344 | 50.590861 |
| 1536 | `sparse_symbolic_analysis` | 23.257446 | 23.488053 | 0.000279 |
| 1536 | `tensor_angular_gemm_complex_and_packing` | 64.685131 | 64.738395 | 63.746687 |
| 1536 | `radial_interpolation_and_strain_numerator_build` | 64.611900 | 65.683797 | 64.962168 |
| 1536 | `spheroidal_source_lift_projection_and_solver` | 12.756882 | 20.855715 | 75.804807 |
| 1536 | `dynamic_H_assembly` | 31.828148 | 31.349141 | 16.116481 |
| 1536 | `sparse_forward_solve_including_rhs_scaling` | 27.814574 | 27.940565 | 14.157576 |
| 1536 | `sparse_adjoint_solve_including_rhs_scaling` | 8.954861 | 18.087886 | 9.196450 |

At degree 1536 the dominant old exclusive scope was packing/key reorder (134.739957 s), followed by numerical LU (98.255247 s). Final dominant exclusive scopes are spheroidal source/lift/solver remainder (75.804807 s), radial interpolation/strain coefficients (64.962167 s), angular GEMM (63.746687 s), and numerical LU (50.590864 s). Packing/key reorder falls to 2.606803 s; a separate final layout-dedup/permutation-build scope must also be counted. The spheroidal exclusive remainder grows materially; do not describe the overall gain as a uniform acceleration of every component or attribute that remainder entirely to cache copies without finer instrumentation.

| Degree | Inclusive modal workspace old→final s | Inclusive tensor sampling old→final s |
|---:|---:|---:|
| 3 | 0.743852 → 0.452665 | 1.031836 → 0.738665 |
| 12 | 2.497830 → 1.608984 | 26.523591 → 16.919074 |
| 128 | 23.172939 → 14.364137 | 18.187299 → 12.191528 |
| 384 | 67.208395 → 44.189614 | 56.011117 → 35.633436 |
| 1536 | 249.484688 → 220.514462 | 281.367339 → 146.927064 |

Inclusive workspace/tensor rows are overlapping scope totals relative to component rows, not additional runtime to sum into them. All stages/counters are retained in the data JSON.

## Factor, RHS, and cache evidence

| Degree | Factorizations matched→final | Final symbolic analyses / reuses | Response misses / hits | Accounted response cache bytes | Peak RSS GiB, user-space |
|---:|---:|---:|---:|---:|---:|
| 3 | 735 → 423 | 3 / 420 | 318 / 312 | 22,860,384 | 0.0558 |
| 12 | 2,625 → 1,377 | 3 / 1,374 | 1,272 / 1,248 | 91,441,536 | 0.1385 |
| 128 | 26,985 → 13,673 | 3 / 13,670 | 13,568 / 13,312 | 975,376,384 | 0.7892 |
| 384 | 80,745 → 40,809 | 3 / 40,806 | 40,704 / 39,936 | 2,926,129,152 | 2.2789 |
| 1536 | 322,665 → 162,921 | 3 / 162,918 | 162,816 / 159,744 | 11,704,516,608 | 8.7412 |

Every degree has 105 explicit radial-l0 bypass calls. Observed factorization count equals the checker-derived expected count exactly. At degree 1536, response capture has 162,816 nonradial misses and 159,744 kernel hits, plus 105 radial bypass factorizations, yielding 162,921 factors. Three symbolic analyses and 162,918 symbolic reuses are recorded. Accounted cache bytes are not a resident-memory measurement; RSS is separately measured by the user-space watcher. The peak observed thread count is two including the supervising process, consistent with one scientific worker.

| Degree | Actual forward RHS matched→final | Actual unit-adjoint basis RHS matched→final | Weighted-adjoint RHS matched→final |
|---:|---:|---:|---:|
| 3 | 1,995 → 1,059 | 945 → 477 | 105 → 105 |
| 12 | 7,665 → 3,921 | 3,780 → 1,908 | 105 → 105 |
| 128 | 80,745 → 40,809 | 40,320 → 20,352 | 105 → 105 |
| 384 | 242,025 → 122,217 | 120,960 → 61,056 | 105 → 105 |
| 1536 | 967,785 → 488,553 | 483,840 → 244,224 | 105 → 105 |

Use `solver_*_rhs_columns` for actual solve work. Batch/cache-request RHS counters include hits and cannot be presented as columns newly solved. The old forward-only production baseline intentionally solved fewer baseline adjoints than the matched control, so its adjoint-column count is not the matched reference.

Final tensor layout cache records 51 hits/one miss at high degree. At 128/384, angular-plan cache has 102 hits/two misses and complex-basis cache 614 hits/ten misses. At 1536 the same counts remain 51/53 and 103/521 respectively: the fixed angular cache budget still causes many rebuilds. These are unchanged from the old layout profile; the recorded gain is not evidence that high-degree angular cache misses disappeared. Final packing avoids 159,432 / 478,920 / 1,916,616 tree-record lookups at degrees 128/384/1536. The degree-1536 radial source-coefficient cache maximum is 534,637,200 accounted bytes, close to its 512 MiB budget.

## Review boundary

The completed evidence supports direct native response reuse plus layout packing and symbolic reuse as faster for these fixed cases, with byte-identical sampled scientific outputs. It does not establish physical convergence, complete high-degree angular kernel equivalence, source-contact correctness in this campaign, statistical speedup confidence, global-consumer performance, or unrestricted cache scalability. Production sources were not changed by this review.


## Follow-up attribution of the new source-scope cost

A bounded read of the accepted candidate identifies a concrete newly introduced cost outside the response-cache timer. In `candidates/combined-final-v4/engine_headers/Windowed3DFields.h:849` and `:915`, both cache-backed pole-only branches execute `adjoint=Eigen::MatrixXcd::Zero(response.forward.rows(),orders)`, where `orders=2*l+1`. These full-order adjoint matrices are then unused: the global visitor is guarded by `!poleOnly&&!forwardOnly`, and pole blocks retain `forward`, compact `adjointUnit`, lift and reduced coefficients. The ordinary compact pole solver path has empty global weighted RHS and does not need these full-order zero arrays.

For degree 1536, the profile records 105 modal passes, and `105*sum_(l=1)^1536(2*l+1)=248,048,640` allocated complex columns per family. With 456 spheroidal rows and 86 toroidal rows, repeated zero-fill alone entails approximately **1,809,762,877,440 spheroidal bytes plus 341,314,928,640 toroidal bytes = 2.151 TB of writes**. This is cumulative traffic, not peak allocation. At degree 384 the corresponding total is 134.967 GB; at 128 it is 15.152 GB. Therefore it is a strong source-level explanation for growth concentrated at high degree, although the current timers do not measure those allocations separately and cannot assign an exact number of seconds to them.

The measured source-scope exclusive timings at 1536 are:

| Scope | Old production s | Matched control s | Final symbolic s |
|---|---:|---:|---:|
| Spheroidal source/lift/solver remainder | 12.756880 | 20.855711 | 75.804807 |
| Toroidal source/lift/solver remainder | 1.983860 | 2.959266 | 12.744211 |

These allocations are **outside** `batch_response_cache_get_or_compute`, so they accrue to those source-scope exclusive timers. Additional RHS validation, source treatment, receiver matrix construction (including `sem.hS(degree).rows()`/`hTk` materialization), reconstruction copies and pole-block construction also occupy the remainder. Their individual seconds are not separable from current instrumentation. The increase from old production to matched control also includes constructing baseline receiver bases because forward-only was turned off; the larger matched-to-final increase is where the new zero-allocation path is especially relevant.

The proposed alternative attribution to operator/matrix hashing is unsupported by the code: `BatchResponseKey::make` in `solvers/batch_response_cache.hpp:33` contains only SEM pointer identity, family, degree and the two frequency bit patterns, so key construction is constant-size. RHS validation occurs before the cache timer. Full-bit RHS signatures, cache lookup, cached-response copies, insertion/allocation and retained signatures occur inside the cache timer. At degree 1536 its entire measured exclusive time is only **3.362198 s**, so those timed cache operations cannot explain the 75.804807-s spheroidal exclusive remainder. No operator matrix is hashed into that key.

The new layout-dedup/permutation stage is 0.547915 s at 1536. Combined with final packing/key-reorder exclusive 2.606803 s, the two main packing scopes total 3.154718 s, versus old packing exclusive 134.739957 s. The standalone layout stage should be counted in any causal packing comparison.

Recomputed exclusive-time closure: for all five degrees and all old-production, matched-control and final-symbolic profiles, the sum of exclusive stage seconds agrees with `physical_case_total` inclusive seconds to at most 2.85e-14 s. The all-components JSON contains no separate closure field; this review computed the equality directly. This supports consistent nested-timer accounting without proving finer attribution within an exclusive scope.

This finding is a performance limitation of the accepted candidate, not a failed scientific-equality gate. No code was changed. A future scoped change could avoid the unused full-order zeros and benchmark that delta, but no resulting speedup is claimed here.
