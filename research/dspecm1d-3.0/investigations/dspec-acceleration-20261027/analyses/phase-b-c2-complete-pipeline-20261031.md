# Phase B and C2 complete-pipeline validation

**Status:** isolated candidate validated and independently reviewed. No production integration.

**Recommendation:** keep the accepted NQ10 mesh (`h=0.025`, input `relative_error=1e-5`) and full radial domain. The lower meshes and paired principal restriction were measured as controls; there is no agreed numerical accuracy budget for adopting either.

## Phase B: fixed meshes and full-domain restriction control

The NQ5/NQ6 nominal and finer meshes use the standard source-aligned h recipe from input `relative_error=1e-5` and the absolute highest retained complex frequency, including damping. The original B jobs retain their 35-bin waveform and measurement controls, including each mesh's own waveform-derived weights. Those waveform comparisons and timings remain valid. The former common-side tables used an incorrect one-sided mapping at artificial element boundaries; they are superseded by a focused replay at bins 9, 25, and 40. The repaired mapper selects a matching outer/inner element at coincident boundaries and uses a same-material containing-element fallback only at artificial noncoincident boundaries.

| Mesh | Elements / native sides | Corrected q0 rows | Worst solid full-domain vs NQ10 forward / adjoint / product | Worst solid restricted minus full forward / adjoint / product |
|---|---:|---:|---:|---:|
| NQ5 nominal | 70 / 199 | 308 | 1.516e-3 / 1.172e-3 / 3.300e-3 | 5.406e-5 / 1.011e-5 / 6.400e-5 |
| NQ5 finer | 132 / 363 | 2,618 | 1.346e-4 / 6.589e-5 / 1.814e-4 | 8.207e-5 / 1.533e-5 / 1.040e-4 |
| NQ6 nominal | 53 / 189 | 308 | 8.588e-4 / 4.658e-4 / 7.305e-4 | 5.906e-5 / 9.988e-6 / 6.641e-5 |
| NQ6 finer* | 96 / 322 | 0 | 1.642e-4 / 9.998e-5 / 2.900e-4 | 6.689e-5 / 1.297e-5 / 8.285e-5 |

Values are maxima across the 9/25/40 keyed one-sided solid location groups, not a full 35-bin frequency-summed kernel comparison and not an accuracy pass threshold. For each pair, full and restricted runs used identical frozen forcing. Every active bin retains all 45,122 requested locations and 401 actual family/degree cut records; each row is classified against its actual applicable cut radii. The fresh restricted controls executed 1,200 kernel-pass Hᴴ profile calls per mesh across all three bins, with nonzero omitted DOFs and positive principal restrictions. This is a whole-replay count, not a per-bin count. The NQ6 finer row is hash-qualified reuse of the original three-bin full/restricted outputs; its 14,000 Hᴴ calls cover the original 35-bin kernel pass, and it was not rerun.

The uncorrected original B frozen-weight kernel timings were 11.35–17.91 s full and 13.80–22.01 s restricted; full used cached responses while restricted forced actual restricted solves, so these are cache-asymmetric descriptive costs, not a truncation speed comparison. Nominal-to-finer solid self-refinement maxima were NQ5: forward `1.421e-3`, adjoint `1.110e-3`, product `3.119e-3`; NQ6: forward `7.938e-4`, adjoint `4.367e-4`, product `6.828e-4`. The per-bin, grouped-location restricted-minus-full solid kernel-product relative error is at most `1.040e-4` in this three-bin sample; it is not a complete-frequency summed-kernel convergence result. Raw fluid adjoint subgroup ratios remain non-roundoff in parts of the inherited convention: up to `2.51` in the packet-side comparison and `1.43` in the C2 weighting comparison. Absolute errors and reference norms are preserved separately; these subgroup diagnostics do not establish universal subgroup equality or a physical SI-strain discrepancy. The principal control removes cross-cut terms and zero-extends eliminated degrees while retaining the full-gravity formulation; it does not establish a physical boundary condition or an accepted truncation tolerance.

The preserved 35-bin waveform controls were close to the frozen NQ10 reference: NQ5 nominal/finer `9.42e-6` / `3.55e-7`; NQ6 nominal/finer `3.20e-6` / `1.80e-7`. Own-waveform paired-refinement differences were `9.62e-6` (NQ5) and `3.16e-6` (NQ6). Common-field comparisons from the old mapper are superseded; the waveform and timing outputs were not rerun or changed.

## C2: combined preweighted contraction

C2 replaces the validated two-pass C route (22.50 s driver elapsed) with one combined forward, unit-direction Hᴴ adjoint, and reconstruction pass. With `P = adjoint.weight * force`, the actual adjoint weights are `direction * conj(P)`: one unit-direction solve yields actual complex products `P * C_unit`, followed by the ordered real reduction. The pass stores only requested complex contractions (25,268,320 bytes at L200), not solved-response history. Static taper zeros select bins 7, 8, and 41 for unit-adjoint skipping during the combined pass. After the own waveform is available, actual `P` is computed and asserted to match that static zero set before weighted reduction. All 35 forward bins remain; unit-adjoint work is 35 calls with skips off and 32 with skips on.

The independent weighting proof passed: assembled-RHS scaling relative error `2.31e-16`; weighted/unweighted adjoint modes `8.00e-13`; pooled fields `1.56e-13`; complex products `3.29e-13`; stored-product mismatch `0`. A separate common-vs-tensor reconstruction control gives forward and adjoint fields `1.014e-15` and `8.692e-16`, and complex products `1.476e-15`. The complete 35-bin kernel matches the preserved two-pass C route at `4.60e-13` along and `1.99e-13` across, and accepted F8a at `4.22e-13` and `2.29e-13`. Waveform components match A4 exactly; trace, spectrum, components, and retained geometry match the two-pass C output exactly. The L3 finite-difference relative error is `6.44e-9` against the existing `2e-4` gate. L200 gravity error is `2.89e-15`; the L3 gravity error is `2.00e-15`; the time-shift adjoint identity is `0.9999999999999997`.

## Matched native reference and timing limits

The fresh accepted-reference run used the saved F8a DSpec command, case, native executable, libraries, and inputs with only an isolated output path. Its full kernels are byte-identical to saved F8a and match C2 at the errors above. Waveforms, measurement norm, and geometry pass the checker gates.

| Run | Driver elapsed | Supervisor elapsed | Peak observed RSS | Threads |
|---|---:|---:|---:|---:|
| C2 native candidate | 15.891 s | 16.296 s | 1.647 GiB | 9 |
| Fresh accepted reference | 16.935 s | 17.594 s | 4.738 GiB | 9 |

This is one matched pair. The accepted F8a command performs an unconditional packed-v-tensor oracle (about `0.52008 s` summed worker elapsed); C2's timed path has that oracle off. The observed driver and supervisor times are reported as measured. No diagnostic time is subtracted and no distributional speedup is claimed. The older Stage 2 native median was `16.8176 s`.

C2 native components were `15.8591 s` for the combined forward/unit-adjoint/reconstruction wall pass, `0.00284 s` for measurement and actual-adjoint spectrum construction, and `0.01093 s` for weighted complex reduction; driver total was `15.8906 s`. In the separate profile run, pipeline wall was `15.9992 s` and the master combined-pass wall scope was `15.9756 s`. Summed worker scopes include `107.96 s` across 35 combined frequency workspaces and `20.15 s` across 14,035 sparse factorizations. These worker sums overlap the master wall pass and must not be added to it.

## Resources, evidence, and boundaries

L200 C2 and the matched reference used eight frequency workers; the B three-bin replay and L3 audit used one solve worker. All used one inner thread and a 16 GiB per-job memory cap. The campaign cap was 20 threads with an aggregate target of at most 50; builds used an eight-thread and 32 GiB cap. The supervisor enforces a user-space session watchdog, not kernel-level aggregate resource accounting. Peak observed counts were nine threads for eight-worker runs and two for one-worker diagnostics.


- Corrected B terminal state, hashes, summary, and resource receipts: [B replay receipt](../receipts/phase-b-mapper-replay-run-v1-terminal-20261031.json), [B preflight](../receipts/phase-b-mapper-replay-preflight-20261031.json).
- C2 equivalence and skip evidence: [C2 terminal receipt](../receipts/phase-c2-l200-controls-v2-terminal-20261031), [C2 independent review](../reviews/phase-c2-final-independent-review-20261031.txt).
- Matched accepted reference: [fresh-reference terminal receipt](../receipts/phase-c2-accepted-reference-native-v1-terminal-20261031.json), [independent review](../reviews/complete-pipeline-native-reference-independent-review-20261031.txt).
- Candidate implementation log: [campaign status](../implementation_status.md); C2 source/build review: [pre-run review](../reviews/phase-c2-pre-run-independent-review-20261031.txt).
- Large raw replay CSVs remain on Tunya under `/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya/phase-b-fixed-mesh-v1-mapper-replay-v1-controls-v1/continuation/phase-b-mapper-replay-run-v1`. The terminal JSON binds their hashes. The fresh reference raw outputs remain under `/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya/phase-c2-preweighted-v1-reference-match-v1-control/continuation/phase-c2-accepted-reference-native-v1`.
- Failed staging/checker/source attempts remain preserved, including the original B checker-only failure and the first C direct/common validation failure. Corrected retries are separately receipted. No production integration, discrepancy investigation, YSpec rerun, figure generation, deletion, commit, or push occurred.

*`NQ6 finer` common-side and restriction rows are reused hash-bound evidence, not a new numerical replay. The fresh q0 count is zero, so the corrected one-sided mapper has no impact on that mesh.*
