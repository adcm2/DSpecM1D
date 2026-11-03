# Phase E Stage 1: matched DSpec and YSpec profiling

**Status: complete; Stage 1 only.** The reviewed run matrix and terminal checker passed. No optimization, production integration, discrepancy work, or Stage 2 was performed.

## Workload and provenance

The run used the selected degree-200 NQ6 case (`h=0.02035382429229197`), all 35 frequency bins 7–41, 208 physical one-sided output sides, and 135 along-track plus 19 across-track angles. DSpec used 8 frequency workers with one inner numerical-library thread. YSpec used 8 CLI workers with its internal library threads pinned to one and `DSPEC_FREQUENCY_WORKERS` explicitly unset. The jobs ran sequentially.

Fresh timer-disabled and timer-enabled artifacts came from the same isolated source/compiler inputs. The pinned YSpec source was commit `280eab0200495391539c0f97e62f978083c307bc`; the source-tree hashes matched the recorded clean-source receipt. The earlier Phase D artifacts were immutable output anchors, not the fresh timing controls. Stage 1 v2 changed only the nesting order of two timer stops after the v1 D profile exposed a LIFO scope error; the numerical formulas and scheduling were unchanged.

This is the same requested case and output, but the two backends do different work. DSpec contracts relative `lnVs` directly. YSpec requests `ParameterSet::Velocities`, returns `vp`, `vs`, `qkappa`, and `qmu`, and the driver selects `vs`; it also computes up, north, and east. Its elastic path forms eight intermediate material products before converting to velocity parameters. YSpec's `DepthSolves` radial integration uses its own mesh and settings (relative tolerance `1e-5`, start decay 11 e-folds, step factor 0.25); the DSpec SEM output mesh does not describe that internal discretization.

## Native timings and profile observations

Times below distinguish driver elapsed time from the supervisor's whole-job elapsed time. Three fresh, timer-disabled controls were run for each backend. The single enabled-profile observation is shown separately.

| Backend | Native driver median (min–max) | Native supervisor elapsed (min–max) | Profile driver / supervisor |
|---|---:|---:|---:|
| DSpec | 7.7416 s (7.7229–7.7644) | 8.0308–8.1110 s | 7.7531 / 8.1624 s |
| YSpec | 4.3220 s (4.3010–4.3256) | 4.5558–4.5849 s | 4.3391 / 4.6014 s |

The timer-enabled driver observations were 0.148% (DSpec) and 0.395% (YSpec) above their fresh native medians. Both differences are smaller than the respective three-run native relative ranges (0.536% and 0.569%). One enabled observation per backend does not isolate instrumentation cost more precisely than that.

The profile maps show where each implementation spent worker time; these are elapsed sums across workers and are not wall-clock shares. DSpec's combined worker-scope sum was 51.7512 worker-seconds, with a 7.7527-second master pipeline wall partition. Selected DSpec reconstruction scopes were Eigen angular GEMM (19.2235 worker-seconds; 147,840 calls), radial-mode interpolation/strain construction (4.6090; 32 calls), and coefficient-row lookup/packing (1.9600; 147,840 calls). Angular plan and weight preparation took 0.2200 worker-seconds, tile-weight packing 0.3543, scatter 0.0943, and material evaluation/contraction 0.1303.

YSpec's worker-exclusive sum was 33.5208 worker-seconds: 1.7385 for forward frequency work and 31.7823 for kernel frequency work. The root-scope totals close exactly. Across both forward and kernel frequency phases, spheroidal radial solves used 16.4677 worker-seconds (13,400 calls), toroidal solves 0.9006 (13,400), strain-table construction 1.0412 (6,432), and angular evaluation 13.9817 (2,050,048). Material contraction was 14.2507 seconds inclusive, but only 0.2690 exclusive because angular evaluation is nested beneath it. YSpec covered 35 forward work items and 35 kernel work items on thread IDs 0–7; three kernel bins (7, 8, and 41) were exact-zero skips, leaving 32 active adjoint bins, matching DSpec's 32 active adjoints.

For a coarse functional comparison, the non-overlapping exclusive worker scopes can be grouped as follows. The phase names are only approximate counterparts: DSpec's modal/solver work and YSpec's radial integration do not perform identical internal work.

| Functional grouping of exclusive worker time | DSpec | YSpec |
|---|---:|---:|
| Modal preparation / radial solves | 24.5664 s | 17.3863 s |
| Radial interpolation and strain tables | 4.6090 s | 1.0412 s |
| Angular synthesis / reconstruction | 22.1749 s | 14.7215 s |
| Material evaluation/contraction and fused tile sample emission | 0.1303 s | 0.2690 s |
| Reduction | included in remaining scopes | 0.0126 s |
| Sampling and remaining scopes | 0.2706 s | 0.0902 s |
| **Total worker-summed exclusive time** | **51.7512 s** | **33.5208 s** |

The largest clearly isolated DSpec angular operation was the GEMM itself at 19.2235 worker-seconds; coefficient-row packing was 1.9600 seconds, versus 0.3543 for tile weights, 0.0943 for scatter, and 0.1303 for the fused material-evaluation/contraction and tile-emission scope. That fused scope is not a pure material-only measurement. Other substantial DSpec worker scopes were numerical factorization (8.4778), dynamic H assembly (6.0790), equilibration/scaling (2.1298), forward Eigen solve (2.0221), adjoint Eigen solve (1.6516), and direct measurement receiver-map assembly (1.5254 worker-seconds). These solver costs are within the modal/solver group, not additional to it. This profiling identifies measured work, but Stage 1 made no optimization or speedup claim.

DSpec recorded 35 workspace preparations, 32 active unit-adjoint frequency calls and 3 exact-zero skips at bins 7, 8, and 41. It reported 14,035 sparse factorizations and 14,011 symbolic-analysis reuses. YSpec's separate forward and kernel API scopes recorded 35 forward work items and 35 kernel work items, including the same three exact-zero skips; both backends therefore had 32 active adjoint bins. Worker sums overlap frequency execution and should not be added to, or subtracted from, driver or supervisor wall time. YSpec's outer API scope, master-thread work-item bridge, worker scopes, and small unscoped driver residual are reported as separate accounting domains.

## Parity, resources, and limits

All nine v2 jobs and the checker completed with return code 0, unchanged inputs, and verified cleanup. Geometry and frequency grids matched exactly. Same-backend results met the `1e-11` output gate against their immutable Phase D references: DSpec output relative-L2 values were exactly zero; YSpec trace, spectrum, and components were exact, and its sampled arrays differed by at most `1.39e-16`. Cross-backend raw sampled arrays are diagnostic only (about `0.1211` along and `0.0367` across, with YSpec as the denominator in this checker); they are not a parity gate, and this stage did not investigate the source of those differences. A prior report used DSpec as the denominator, so those percentages are not directly comparable and do not indicate a changed physical discrepancy.

The fresh build used the 8-thread/32-GiB build cap and peaked at 6 observed threads and 2.3101 GiB RSS. Numerical jobs used the 20-thread/16-GiB per-job caps and 8 configured workers; the supervisor observed at most 9 process threads, including its own process accounting. Peak observed RSS was 1.1211 GiB for DSpec and 0.3288 GiB for YSpec. The aggregate ceiling was 50 threads. Resource enforcement used the existing user-space session watchdog; this is not kernel-level aggregate thread or memory accounting.

The v1 attempt is preserved: its three D native controls completed, then its D profile aborted with return code -6 due to the timer-scope LIFO violation. The reviewed v2 stop-order correction resolved that failure; all v2 jobs passed. These timing results are profiling evidence, not an optimization comparison. The separate historical Phase D timings are context only; no cause is inferred from any difference between campaigns.

## Evidence

- [Terminal snapshot with v1 failure, v2 build, all nine v2 jobs, output hashes, and checker summary](../receipts/phase-e-stage1-terminal-snapshots-20261031.json), SHA-256 `1778af89bb3acaca2f95e723502b33fba37fddcbe5475346627bb3f56ece303f`.
- [Reviewed Stage 1 plan](../plans/phase-e-stage1-backend-profiling-plan-20261031.md).
- Run plan SHA-256 `07a69d7996fc57f74f335c23e525ba6dce3fc2c234c2b241290c0c3d50c71d5c`; checker SHA-256 `44095eb7785df315a0443f0c43e24581297aaa57d9fc484f6c60960a39bc6e94`.
- Remote terminal checker summary SHA-256 `b20354e3a8f4e6876e2f4079bbfc43f7526e2e9d6b5d453d1ce1445a9644621c`; its hash-bound remote path is recorded in the terminal snapshot.
