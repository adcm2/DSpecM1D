# DSpecM1D acceleration: Stages 1–2 final evidence

**Status:** Stages 1–2 completed; stop boundary is before Stage 3. Work is isolated to reviewed candidates and L200 evidence. No production checkout or integration was changed.

## Outcome

The measured case is degree 200 with nq=10, refinement 0.025, 35 frequency bins (7–41), attenuation enabled, and 293 native sides across 135 along-path and 19 across-path angles; the recorded kernel parameter is `ln_vs_at_fixed_vp_density_and_effective_Q` in s m^-3, with native-side-major then angle float64 output and no source-shell contacts. Stage 2 frequency scheduling with eight frequency workers reduced median complete-driver elapsed time to **16.817634 s**. This is **3.107×** faster than the contemporaneous degree-eight anchor (52.250099309 s), **3.082×** faster than Stage 1’s native median (51.840044503 s), and **6.359×** faster than the Stage 2 frequency-one median.

The contemporaneous degree-eight anchor took 52.250099309 s overall. Its legacy kernel timer was 40.939090800 s; adding its separately measured 0.519657124 s serial diagnostic time gives an inclusive kernel value of 41.458747924 s. The legacy timer alone is not definitionally comparable to Stage 2 `kernel_s`. The Stage 2 frequency-one route is slower than the Stage 1 native route: its median is 106.946436 s, about 2.063× as long as Stage 1 elapsed time (equivalently, Stage 1 is 0.485× the F1 elapsed time). The speedup therefore comes from parallelizing the whole frequency path, not from the serial frequency path.

## Stage 2 timing

Times below distinguish the driver’s `result.json` elapsed from the external supervisor’s elapsed. `kernel_s` includes diagnostic tensor validation and ordered reduction; `kernel_parallel_pass_wall_s` is parallel pass wall only.

| Route | Driver elapsed min / median / max (s) | Supervisor elapsed min / median / max (s) | Forward median (s) | Kernel inclusive median (s) | Parallel pass wall median (s) | Max RSS (GiB) | Peak threads |
|---|---:|---:|---:|---:|---:|---:|---:|
| Frequency 1 | 106.817846 / 106.946436 / 107.038497 | 107.339535 / 107.471465 / 107.485662 | 62.305996 | 44.625336 | 44.607949 | 2.2698 | 2 |
| Frequency 8 | 16.670132 / 16.817634 / 16.842277 | 17.435563 / 17.589388 / 17.615503 | 9.401452 | 7.356635 | 7.291469 | 4.7398 | 9 |

The reported `forward_s` includes forward work plus unit-adjoint response capture on both routes. Stage 2 native repeats: F1 `106.817845652 / 106.946435628 / 107.038496691 s`; F8 `16.670131539 / 16.817633768 / 16.842276646 s`. Each used eight total frequency workers at F8 and one solve worker; F1 used one frequency worker. Thread peak includes supervisor overhead.

Stage 1 native median was 51.840044503 s (min/max 51.832922741 / 52.140116475 s); its supervisor elapsed median was about 52.480971456 s. Stage 1 and Stage 2 kernel timings have different diagnostic inclusion, so only the complete-driver elapsed is compared across stages.

## Profile interpretation and validation

The profile master wall partition closes at 16.926455471 s with zero residual. Its kernel phase totals 7.438092719 s including diagnostics and ordered reduction; the frequency parallel-pass wall is 7.372566516 s and parent ordered reduction is 0.002826909 s. Summed diagnostic worker elapsed is 0.519680216 s, reported separately and **not subtracted** from wall. Worker exclusive sums below are thread-time, not wall time.

Largest non-diagnostic worker-exclusive sums: modal solver 49.666565 s; tensor/angular synthesis in the kernel phase 30.255005 s; radial sampling/reconstruction 10.825835 s; initial receiver basis projection 8.513978 s; kernel receiver basis projection 5.029144 s; batch response cache 2.183394 s. These are summed across workers and must not be read as elapsed wall time. The raw scope rows preserve names; exact component labels are supplied in the explicit map.

Stage 1 reported 17 symbolic analyses and Stage 2 reported 24; this difference reflects changed TLS/instrumentation history and is not treated as a numerical-work reduction. The logical numerical, RHS and response-cache invariants below remain unchanged. Counter invariants: 14,070 sparse factorizations; 42,070 solver forward RHS columns; 21,000 unit-basis RHS columns; 70 weighted-adjoint RHS; response-cache 14,000 misses and 14,000 hits; 70 radial L0 bypasses; 70 zero-weighted-adjoint skips. The profile run overhead relative to the F8 native median was 0.649% by driver elapsed time; relative to F8a it was 1.538%. Contribution storage peaked at 12,634,160 bytes.

Stage 2 passed all 13 controls and all 16 measurement jobs. The jobs checked the preserved degree-eight anchor, L3 and scientific byte parity, profile accounting, output invariants, forward +/- and finite-difference checks. All task inputs remained unchanged and cleanup was verified. Stage 1 passed all 12 controls and 10 measurement jobs; its finite-difference relative error was `6.440263106e-09`.

## Provenance and evidence

- Stage 1 accepted source and report: `receipts/stages12-stage1-v1-frozen-inputs/candidate/`, `analyses/stages12-stage1-v1-final-report.json`.
- Stage 2 candidate: `candidates/stages12-stage2-frequency-v1/`; exact seed is the Stage 1 frozen candidate. Main source hashes are recorded in the JSON report.
- Stage 2 build validation and binary hashes: `receipts/stages12-stage2-v1-build-completed/build-validation.json`.
- Stage 2 controls state: `receipts/stages12-stage2-v2-controls-completed/raw-task/state.json`.
- Stage 2 terminal measurement tree and state: `receipts/stages12-stage2-v3-measure-completed/raw-task/`; state SHA256 `16d879978fee1682770fc9eda9485a9199d9ad2d8b186d6081e62799a0c5094e`.
- Profile packet: `receipts/stages12-stage2-v3-profile-packet/`, including raw scope rows, raw counters, explicit mapping, component CSV, and the profile accounting analysis.
- Failures were preserved: the initial controls run exposed the supervisor’s OpenMP thread-limit reset; a later degree-anchor checker mismatch was corrected through a fresh strict checker/plan. Neither gate was relaxed, and both failed task roots remain in receipts.
- Historical YSpec elapsed `5.641430309 s` is context only and was not rerun.

**Stop boundary:** no Stage 3 code, testing, integration, or further numerical runs are included.
