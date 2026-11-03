# Phase A final independent evidence review

27 October 2026. **PASS for completed bounded execution, receipt integrity, and internal parity. Scientific promotion remains blocked by the unchanged original-global FAIL gates.** This review performed no SSH, builds, science runs, code edits, figures, or production changes.

## Completed evidence and integrity

- Global v8: 19/19 jobs completed with return code zero; post-analysis: 10/10 jobs completed with return code zero. All jobs record unchanged inputs, empty failure reasons and verified cleanup. The preserved earlier serial queue remains failed because its unfinished degree-1536 run was deliberately stopped; its completed lower-degree records remain preserved.
- Global state SHA256 `17c15c06b5cce5905a96c7c5df9618ad97dbd6233e2cc5e707a7d8451652513a` matches both its verification manifest and the receipt-bound post aggregate. Post state SHA256 `748773b9126b3ed6c57ae371892ea9380ebc67cc8b455c681cede742197b8f6d` matches its verification manifest.
- Independently rehashed 78 collected global artifacts and 53 collected post artifacts with no mismatches. Twenty global binary payloads (ten kernels and ten response spools) and four paired kernel probes were intentionally not copied locally. Their hashes and equality claims are preserved in the verified remote output/checker receipts; this review did not reread those binary bytes locally.
- Nine global internal gates explicitly requested exact comparison and PASS with every trace/spectrum/kernel byte flag true and all reported numeric maximum differences zero. At degrees 3/12/128/384, Generated 32-thread outputs match the preserved one-thread Generated outputs; Stored 32-thread outputs match Generated 32-thread outputs. At 1536 only Generated-versus-Stored 32-thread parity was checked: **there is no completed serial 1536 reference**.
- Four post paired-repeat gates explicitly required byte equality and PASS, including probe-angle metadata and all six kernel channels. The previous native v5 campaign already passed its separate 11 reports and 1536 native/profiled parity. None of these internal equalities repairs the direct-original versus global mismatch.

## Native route: retain the low-degree regressions

The isolated paired sequence began 6.90 s after the global queue finished. Its dependency receipt says terminal completion only, with no scientific acceptance implied.

| Degree | Paired v4 process s | Paired v5 process s | v5 slowdown |
|---:|---:|---:|---:|
| 128 | 26.912080 | 29.867697 | 10.98% |
| 384 | 80.847368 | 86.122401 | 6.52% |

These isolated observations confirm that the earlier v5 regressions at 128 and 384 cannot be dismissed as overlap contamination. All profile counters remain identical between paired v4 and v5. The slower observations localize to assembly/scaling/factorization and solver remainder, partly offset by lower source-scope costs; they do not establish a compiler, allocation, or memory-layout root cause.

| Degree | Exclusive v5−v4 assembly s | Scaling s | Numerical LU s | Shared-solver remainder s | Spheroidal source remainder s | Toroidal source remainder s |
|---:|---:|---:|---:|---:|---:|---:|
| 128 | +1.143823 | +0.532492 | +0.549593 | +0.411604 | +0.011210 | -0.070573 |
| 384 | +3.193279 | +1.249845 | +1.485702 | +1.075520 | -2.061273 | -0.626328 |

At 1536, the already completed native v5 observation is 313.360937 s versus v4 native 369.676990 s (about 15.23% less); profiled v5 is 314.662630 s versus v4 370.683470 s. No new paired 1536 repeat was added. Relative to the earlier production observations, the isolated v5 128/384 runs still improve on 41.856418/123.951199 s, but there is no all-degree v5-over-v4 speedup.

## Global route: full angular output on selected sides

| Degree | Selected native sides | Full angular points per side | Generated process s | Stored process s | Generated / Stored peak RSS GiB | Observed peak threads |
|---:|---|---:|---:|---:|---:|---:|
| 3 | [0, 28, 59, 100, 149, 154, 172, 173] | 91 | 21.938509 | 12.205139 | 0.0369 / 0.0373 | 34 / 34 |
| 12 | [0, 28, 59, 100, 149, 154, 172, 173] | 1,225 | 27.510067 | 21.936852 | 0.0452 / 0.0452 | 34 / 34 |
| 128 | [59, 100] | 133,640 | 31.967635 | 30.423811 | 0.1902 / 0.3416 | 34 / 34 |
| 384 | [59, 100] | 1,184,260 | 104.868045 | 100.432640 | 0.7436 / 4.9272 | 34 / 34 |
| 1536 | [59] | 18,929,680 | 509.233458 | 483.668116 | 6.8598 / 276.7871 | 34 / 34 |

At 1536 Stored is about 5.02% faster end to end than Generated (483.668116 versus 509.233458 s), while observed user-space peak RSS increases from 6.8598 to 276.7871 GiB. Its shared Wigner table is 290,382,151,640 bytes (270.439 GiB), and the response spool is 12,125,305,023 bytes (11.2926 GiB). These are distinct objects; the spool is not the large table allocation.

Both runs request one outer worker and 32 inner GSH threads, with explicit OpenMP policy, dynamic teams off, runtime inner limit 32, successful setup-pool cleanup, and observed peak 34 threads including orchestration. The global queue cap is 40; all observed thread peaks fit it. Stored uses the reviewed 384-GiB working-memory configuration; Generated uses 64 GiB. Preflight bounds and user-space measured RSS are different quantities, and the watcher does not provide kernel-enforced aggregate accounting.

Do not claim a 32x speedup. There is no complete serial 1536 timing. Lower-degree serial Generated observations were 3.501654, 5.242222, 34.400080 and 259.722564 s at 3/12/128/384; corresponding 32-thread observations are 21.938509, 27.510067, 31.967635 and 104.868045 s. Small cases became substantially slower, while 384 improved about 2.48x. These are whole-process observations across reviewed source versions, not an ideal thread-scaling law.

## Component accounting and closure

For all ten global profiles, five original v5 profiled cases and four paired repeats, recomputed sum(exclusive stage seconds) equals the physical-case inclusive root within 3.55e-15 s. This verifies the nested profile ledger. It does not imply that every global-consumer operation has its own component timer: consumer execution in worker threads is separately described in its manifest and appears largely in the root exclusive remainder.

| 1536 consumer component | Generated s | Stored s |
|---|---:|---:|
| Total consumer wall | 235.023090 | 206.644793 |
| Serial grid/context setup | 4.124653 | 40.600084 |
| Worker execution wall | 226.380822 | 161.319530 |
| Worker spool_reads | 0.368226 | 0.356155 |
| Worker coefficient_rotation_layout | 22.185431 | 22.184949 |
| Worker gsh_inverse_fft | 168.187111 | 103.136549 |
| Worker point_kernel_contractions | 28.677233 | 28.689353 |
| Worker frequency_accumulation | 2.781166 | 2.716202 |
| Worker output_writes | 2.864749 | 2.938104 |
| Unsplit consumer wall outside setup+worker execution | 4.517615 | 4.725179 |
| Worker elapsed minus named worker stages | 1.311501 | 1.291658 |

Stored reduces GSH inverse/FFT time from 168.187111 to 103.136549 s, while serial table/context setup rises from 4.124653 to 40.600084 s. Rotation and point contractions remain roughly 22.18 and 28.68 s. The visible residuals above preserve incomplete fine-grained timing coverage; they were not forced to zero. Per-worker duration sums are wall durations and are not parallel CPU-time accounting. The two 1536 cases execute 52 frequency/side transforms, skipping one exactly zero-mixing frequency; the selected side is fluid, so only 104 channels are transformed and 520 structural-zero fluid channels are skipped. The complete component data for every global run are in the accompanying JSON.

## Scientific gate and descriptive resolution

**The original-global fixed 2e−13 acceptance controls still FAIL.** At degree 3, scaled maximum differences in L and kappa are 4.335820e−13 and 3.243901e−13; at 12 they are 3.338206e−13 and 3.098487e−13. Exact thread/storage parity is an internal diagnostic, not independent accuracy evidence. No threshold was relaxed and no global scientific acceptance is claimed.

The resolution helper finished with `DESCRIPTIVE_COMPARISON` and `acceptance_gate=null`. It compares degree 128 with 384 through approximate bilinear interpolation on a 720-point great-circle section, retaining 596 points after excluding 5 degrees at both endpoints. Sources are preserved serial Generated full-angular outputs at sides 59 and 100; source/receiver are (10°,20°) and (45°,75°). It is not a new solve or a convergence certificate.

| Native side / radius | Active channel | Max difference / max amplitude | Difference L2 / larger field L2 |
|---|---|---:|---:|
| 59 / 2915375.000 m | kappa | 8.2234% | 8.6605% |
| 100 / 4755714.286 m | A | 14.8639% | 15.6598% |
| 100 / 4755714.286 m | C | 17.6543% | 16.6250% |
| 100 / 4755714.286 m | F | 16.1752% | 16.3763% |
| 100 / 4755714.286 m | L | 10.8535% | 9.9580% |
| 100 / 4755714.286 m | N | 8.2949% | 8.9255% |

Fluid side 59 has structural-zero A/C/F/L/N; solid side 100 has structural-zero kappa. Those zeros are not convergence evidence. Nonzero differences remain approximately 8–18% under these descriptive normalizations; interpolation error, degree truncation and the unresolved original-global mismatch are not separated.

## Final boundary

The native direct evidence uses 64 angles over all 179 selected radial sides at high degree. The global evidence uses complete angular grids on only eight sides (3/12), two sides (128/384), or one fluid side (1536). These domains and algorithms cannot be compared as equal-output timing workloads. There is no 179-side global volume, no new figures, no contact-kernel validation, no physical convergence claim, and no scientific promotion from bookkeeping completion.
