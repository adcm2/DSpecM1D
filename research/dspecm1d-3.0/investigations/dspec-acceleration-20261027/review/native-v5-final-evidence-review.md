# Final-v5 direct native evidence review

27 October 2026. **PASS for this completed, bounded evidence campaign.** Read-only review; no SSH, build, scientific execution, or code edits. The prior v4 review remains a separate preserved record. This review does not promote the experimental global consumer or authorize production integration.

## Integrity and acceptance gates

- State and `analyses/native-final-v5-all-components.json` both report completed. All 17 jobs completed with return code zero, empty failure reasons, unchanged input hashes and verified cleanup. State SHA256 `17dbad9b431fc3cd484c46d48f9670d40bf2ac51aad4b8eebbef3bc6b440bbbc` agrees with both the collection verification manifest and all-components analysis.
- All 11 comparison reports PASS: five response-reuse attribution gates, five separate explicitly byte-exact gates, and one explicitly byte-exact native/profiled parity gate at degree 1536. All byte-equality flags are true, and every reported numeric maximum difference is zero. The separate exact gates satisfy a stronger condition than the tolerance-enabled reuse reports alone.
- Independently rehashed all 74 locally collected files in the 80-entry verification manifest; no mismatches. The 6 kernel binary payloads were not collected locally. Their exact equality is supported by hash-verified remote checker reports and output-hash manifests, rather than a new local reread of those binary bytes.
- The five profiled cases use `kernel-acceleration-20261027-final-05/final-profiled/physical-driver`; the native parity case uses `final-05/final-native/physical-driver`. Input executable hashes are preserved in the data JSON. Results use direct response reuse with symbolic reuse enabled, pole-only mode, forward-only off, no capacity/preflight fallback, zero receiver-basis reconstruction error, no response spool, and no global full transforms.
- For every tested degree, scientific output hashes match v4→v5 exactly for trace, spectrum, radial metadata, angular grid and the six-component kernel binary/probe (including probe-angle coordinates where applicable). The previous v4 review established old-production→matched-control→v4 equality, so these completed gates extend that provenance chain to v5. Result files may legitimately differ in timing and route/cache metadata.

## Runtime, including the lower-degree regressions

| Degree | Old production s | Matched control s | v4 symbolic/profiled s | v5 symbolic/profiled s | v5/v4 elapsed ratio |
|---:|---:|---:|---:|---:|---:|
| 3 | 2.203356 | 2.072617 | 1.601514 | 1.495605 | 0.9339 |
| 12 | 30.078734 | 29.959697 | 19.517453 | 19.587391 | 1.0036 |
| 128 | 41.856418 | 44.221491 | 27.113679 | 29.636236 | 1.0930 |
| 384 | 123.951199 | 130.525318 | 80.840457 | 86.567986 | 1.0708 |
| 1536 | 532.983888 | 559.869149 | 370.683470 | 314.662630 | 0.8489 |

At degree 1536, profiled elapsed time decreases from **370.683470 to 314.662630 s** (15.11% less, 1.178x). Native decreases from **369.676990 to 313.360937 s** (15.23% less). The v5 profiled/native difference is about 0.415%. Compared with the earlier production 532.983888-s observation, the final v5 profiled observation is 1.694x faster; the matched 559.869149-s control is a distinct comparison.

**Do not claim improvement at every degree.** Degree 128 increases from 27.113679 to 29.636236 s (about 9.30%); degree 384 increases from 80.840457 to 86.567986 s (about 7.08%). Degree 12 is also slightly slower. These observations occurred during overlapping build/old-serial-global activity in parts of the lower-degree campaign, as documented in `implementation_status.md:75`; they are not isolated repeated benchmarks. In particular old serial global degree 1536 was stopped at Unix 1793117385.8215625: the v5 degree-384 interval 1793117338.9127505–1793117425.5017455 overlaps that activity for approximately 46.91 s. Degree 128 ran entirely before that stop and had overlapping build activity. The main new global 40-thread queue was deliberately deferred until the native queue finished. Neither uniform speedup nor a precise contention penalty can be inferred from these single-run measurements.

The old production timing uses forward-only on; matched control and v4/v5 use forward-only off with response capture. Preserve this distinction when quoting practical before/after versus matched-control gains.

## Component attribution

All entries below are **exclusive** seconds. They are not added to their containing inclusive workspace/tensor totals.

| Degree | Component | v4 s | v5 s |
|---:|---|---:|---:|
| 3 | `spheroidal_source_lift_projection_and_solver` | 0.060989 | 0.061362 |
| 3 | `toroidal_source_lift_projection_and_solver` | 0.006085 | 0.005763 |
| 3 | `sparse_numerical_factorization` | 0.166184 | 0.161765 |
| 3 | `dynamic_H_assembly` | 0.062569 | 0.063890 |
| 3 | `radial_interpolation_and_strain_numerator_build` | 0.131180 | 0.126932 |
| 3 | `tensor_angular_gemm_complex_and_packing` | 0.235051 | 0.233362 |
| 3 | `tensor_coefficient_packing_and_key_reorder` | 0.026040 | 0.025289 |
| 12 | `spheroidal_source_lift_projection_and_solver` | 0.236737 | 0.244715 |
| 12 | `toroidal_source_lift_projection_and_solver` | 0.023428 | 0.021819 |
| 12 | `sparse_numerical_factorization` | 0.565420 | 0.554170 |
| 12 | `dynamic_H_assembly` | 0.192337 | 0.236076 |
| 12 | `radial_interpolation_and_strain_numerator_build` | 0.493984 | 0.498934 |
| 12 | `tensor_angular_gemm_complex_and_packing` | 9.997484 | 10.021242 |
| 12 | `tensor_coefficient_packing_and_key_reorder` | 0.588100 | 0.576746 |
| 128 | `spheroidal_source_lift_projection_and_solver` | 2.576910 | 2.566753 |
| 128 | `toroidal_source_lift_projection_and_solver` | 0.300538 | 0.229709 |
| 128 | `sparse_numerical_factorization` | 4.760201 | 5.276235 |
| 128 | `dynamic_H_assembly` | 1.337878 | 2.458496 |
| 128 | `radial_interpolation_and_strain_numerator_build` | 5.229429 | 5.201062 |
| 128 | `tensor_angular_gemm_complex_and_packing` | 5.314871 | 5.255772 |
| 128 | `tensor_coefficient_packing_and_key_reorder` | 0.215656 | 0.230414 |
| 384 | `spheroidal_source_lift_projection_and_solver` | 9.432723 | 7.417334 |
| 384 | `toroidal_source_lift_projection_and_solver` | 1.313898 | 0.689351 |
| 384 | `sparse_numerical_factorization` | 13.299868 | 14.835895 |
| 384 | `dynamic_H_assembly` | 3.982850 | 7.229768 |
| 384 | `radial_interpolation_and_strain_numerator_build` | 15.744014 | 15.983906 |
| 384 | `tensor_angular_gemm_complex_and_packing` | 15.979531 | 15.983552 |
| 384 | `tensor_coefficient_packing_and_key_reorder` | 0.654231 | 0.684141 |
| 1536 | `spheroidal_source_lift_projection_and_solver` | 75.804807 | 21.301671 |
| 1536 | `toroidal_source_lift_projection_and_solver` | 12.744211 | 2.718416 |
| 1536 | `sparse_numerical_factorization` | 50.590861 | 52.776542 |
| 1536 | `dynamic_H_assembly` | 16.116481 | 16.276080 |
| 1536 | `radial_interpolation_and_strain_numerator_build` | 64.962168 | 66.655361 |
| 1536 | `tensor_angular_gemm_complex_and_packing` | 63.746687 | 63.929936 |
| 1536 | `tensor_coefficient_packing_and_key_reorder` | 2.606805 | 2.715958 |

The v5 source-scope change is visible at 1536: spheroidal exclusive remainder falls from 75.804807 to 21.301671 s, toroidal from 12.744211 to 2.718416 s. Their combined reduction is **64.528931 s**, while total elapsed improves by 56.020840 s. Other scopes become slower in this observation, so the component reduction should not be equated to the net runtime gain. The attribution is consistent with the reviewed v5 removal of the two unused full-order zero matrices identified in the v4 review; current timers still do not measure zero allocation in a separate scope.

The largest remaining exclusive scopes at 1536 are `radial_interpolation_and_strain_numerator_build` (66.655361 s), `tensor_angular_gemm_complex_and_packing` (63.929936 s), `sparse_numerical_factorization` (52.776542 s), `spheroidal_source_lift_projection_and_solver` (21.301671 s), `dynamic_H_assembly` (16.276080 s).

Exclusive-stage closure was recomputed for every v5 profiled case: sum(exclusive seconds) agrees with the physical-case inclusive root within 0 s. No separate closure field was assumed.

## Unchanged solve/cache workload

Every recorded profile counter is identical between v4 and v5 at every tested degree. This strengthens the interpretation that the observed source-scope saving comes from removing unused work rather than dropping numerical solves or sampled output work.

| Degree | Factors | Symbolic analyses / reuses | Nonradial response misses / hits | Accounted response bytes | Peak RSS GiB, user-space |
|---:|---:|---:|---:|---:|---:|
| 3 | 423 | 3 / 420 | 318 / 312 | 22,860,384 | 0.0555 |
| 12 | 1,377 | 3 / 1,374 | 1,272 / 1,248 | 91,441,536 | 0.1376 |
| 128 | 13,673 | 3 / 13,670 | 13,568 / 13,312 | 975,376,384 | 0.7880 |
| 384 | 40,809 | 3 / 40,806 | 40,704 / 39,936 | 2,926,129,152 | 2.2747 |
| 1536 | 162,921 | 3 / 162,918 | 162,816 / 159,744 | 11,704,516,608 | 8.7363 |

Each case retains 105 explicit radial-l0 bypass calls; observed factor count equals expected factor count. At 1536, actual forward/unit-basis/weighted-adjoint solved RHS counts remain 488,553 / 244,224 / 105. Batch request-column counters include cache hits and are not actual substitution counts. Accounted cache bytes are not resident memory; measured peak RSS is reported separately. Peak observed threads remain two including the supervisor, consistent with one scientific worker.

Angular cache counts at 1536 remain 51 hits/53 misses for angular plans and 103 hits/521 misses for complex bases, while tensor layout cache remains 51 hits/one miss. V5 did not remove these remaining cache rebuilds. All complete counter/stage maps are in the data JSON.

## Coverage and limits

The fixed physical/discretization setup remains 30 SEM elements, six GLL nodes per element, 179 selected kernel radial sides and 52 kernel frequency bins. Degrees 3 and 12 contain complete target angular grids; degrees 128/384/1536 contain all 179 radial sides at 64 diagnostic angles (11,456 points), not full high-degree angular outputs. Contacts are off. Exact equality therefore establishes the sampled direct-route outputs and returned waveform/spectrum/energy quantities in these cases; it does not establish physical convergence, whole-Earth integrated kernels, contact-kernel correctness in this campaign, complete high-degree angular equivalence, repeated-run confidence intervals, or global-consumer completion.

No implementation or further run was performed as part of this evidence review.
