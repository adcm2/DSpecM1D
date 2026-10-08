# Final Sol-high integration review

## OVERALL VERDICT: PASS WITH NON-BLOCKING NOTES

**Final recommendation: READY FOR HUMAN REVIEW WITH NON-BLOCKING NOTES**

Reviewed complete baseline-to-current-feature diff: YES

This was a fresh read-only review of the complete Stages 1–9 implementation, including production code, surrounding solver paths, tests, configuration, saved validation/performance evidence, and untracked validation artifacts. The review did not assume the earlier stage verdicts were correct. This report is generated after the frozen snapshot and is excluded from all statistics below.

## Scope and provenance

- **Baseline:** `889452b6790d2c60afacc156b4c7512990758aee`
- **Feature HEAD:** `8f7e7720bb2201bba035ae4abc8e919857a43e32` on `dev/cowling`
- **Previous approved checkpoint:** `b66cccfc06b922a73e5c0968fff9ea51a0789947`
- **Commits since baseline:** `5bc5ed840b6466b9c7af60bb5cda95ef41027e78`, `b66cccfc06b922a73e5c0968fff9ea51a0789947`, `8f7e7720bb2201bba035ae4abc8e919857a43e32`
- **Frozen Cowling scope, baseline to current:** 39 files: production 13 files, +338/-96; tests 4 files, +1226/-1; validation 20 files, +1150/-0; documentation/status 2 files, +463/-0. Cowling patch: 237,778 bytes and 3,847 patch lines.
- **Frozen raw scope, baseline to current:** 43 files. In addition to the Cowling scope, four historical unsupported benchmark files add 46 textual lines plus two binary PDFs. Raw patch: 1,107,677 bytes and 16,831 patch lines. Those reports are unrelated historical material and are not Cowling validation or performance evidence.
- **Frozen Cowling scope, checkpoint to current:** 26 files: production 3 files, +84/-44; test 1 file, +564/-2; validation 20 files, +1150/-0; documentation/status 2 files, +201/-2. Patch: 163,399 bytes and 2,293 patch lines.
- **Frozen raw scope, checkpoint to current:** 30 files, additionally including the four historical benchmark files (+46 textual lines and two binaries). Patch: 1,033,298 bytes and 15,277 patch lines.

All 97 frozen file hashes matched before this report was written, including `tests/test_solver_api.cpp` at `ff21acaa2921f5c219ad507a31b520fd875bd0fe4a7124597da389e95303ab9b`. The four frozen patch sizes and line counts matched their manifest.

## Findings

### Blocking findings

None.

### Non-blocking findings

1. The scientific validation is intentionally narrow: one PREM/no-attenuation source-receiver geometry, one receiver, three degrees, five nearby sampled bins, and selected two-mesh checks. The 1.5625 mHz bin spacing and 2.6526 mHz imaginary cyclic damping limit resonance interpretation. There is no matched MINEOS full/Cowling comparison and no time-domain ringing test. These limits are accurately stated in `validation/cowling/README.md`.
2. Hard-cutoff evidence supports the sampled 20.3125 mHz transition, not every source, receiver, degree, or cutoff. Maximum measured switch excess is 0.00119%, 0.353%, and 0.748% for l=1, 20, and 100, respectively, versus substantially larger natural adjacent-bin full/full changes.
3. Performance evidence is one single-thread degree-20 case with five warm repetitions. The 8.777 ms construction number builds the shared mesh and both layouts, so it is not a full-only setup comparison. The saved report labels this correctly.
4. The raw PR includes four unrelated historical benchmark artifacts. Their README and provenance correctly mark them unsupported and separate from this feature, but they increase the raw review footprint substantially.

## Physics review

The implementation constructs Cowling matrices explicitly alongside the existing full matrices in the same quadrature loops (`DSpecM1D/src/SEM/SEMConstructor.h:335`). It retains elastic, attenuation, inertia, and background-gravity terms. In particular, the Cowling U-U diagonal contains `-rho*g*r` while omitting the full formulation's `4*pi*G*rho^2*r^2` contribution (`SEMConstructor.h:391`), and U-V retains `rho*g*r` (`SEMConstructor.h:407`). Cowling receives no P indices, U-P/V-P/P-P entries, or exterior-potential boundary term; `hSC` also adds no surface potential condition (`DSpecM1D/src/SEM/SEMMatrices.h:31`). A direct weak-form coefficient oracle independently checks background gravity, mass, attenuation, U-V coupling, and the exact full-versus-Cowling self-gravity difference (`tests/test_sem_component.cpp:201`).

The compact map preserves U at every element boundary and creates the extra V degree of freedom only at fluid-solid transitions (`DSpecM1D/src/SEM/SEMLG.h:34`). Independent quotient-mesh numbering tests cover all-solid, solid-solid, and solid-fluid-solid cases. Source and receiver methods reuse the existing physical interpolation and four source columns, selecting only the layout-specific indices (`SEMForceSpheroidal.h:85`, `SEMReceivers.h:398`).

## Numerical and solver review

`startElementSph` is unchanged; all three `allIndicesSph` overloads map the same selected element through either `ltgS` or `ltgSC` (`DSpecM1D/src/StartElement.h:153`). Single-SEM routing converts the bin frequency back to physical mHz and selects Cowling at `frequency >= cutoff` (`DSpecM1D/src/FullSpecSingleSem.h:155`). Full and Cowling occupy separate contiguous regions, each with its own matrices, truncation cadence, and fresh SparseLU object (`FullSpecSingleSem.h:173`). Multi-SEM applies the same rule within every chunk and handles full-only, Cowling-only, and crossing chunks (`DSpecM1D/src/FullSpecMultiSem.h:248`). Symbolic state cannot cross layouts, and `nskip` restarts within each formulation as required.

The original-baseline regression serialized full `hS`, `hSa`, `pS`, forces, receivers at l=1,2,7, and single/multi spectra at l=1–4 over requested 2–8 mHz with two receivers and attenuation off/on. Baseline and feature files are byte-identical at 9,055,120 bytes with SHA-256 `ec838d8173387a019879565053f52fa7d7dcbe9cc449c8bae81e8370317e9c2b`. This is strong evidence that the disabled full-spheroidal path remains unchanged. Direct diff inspection confirms that the radial and toroidal branches were not modified.

## Explicit answers to the 21 final questions

| # | Answer |
|---:|---|
| 1 | **Yes.** The Cowling weak form is explicitly assembled as U,V rather than extracted from the full principal block. |
| 2 | **Yes.** Background `g(r)` terms remain; perturbed-potential terms are excluded. |
| 3 | **Yes.** The `4*pi*G*rho^2 U U~` contribution appears only in the full U-U entry and is directly tested. |
| 4 | **Yes.** Cowling has no P DOF, P coupling, P-P stiffness, or exterior P boundary contribution. |
| 5 | **Yes.** The two-field map is compact, contiguous, and exhaustively checked on controlled meshes. |
| 6 | **Yes.** U is continuous and V slips at fluid-solid boundaries; V shares elsewhere. |
| 7 | **Yes.** Corresponding physical source rows and receiver columns compare exactly between layouts. |
| 8 | **Yes.** The physical start-element criterion is shared; only final DOF mapping differs. |
| 9 | **Yes.** Default cutoff zero preserves full gravity; the baseline regression is byte-identical. |
| 10 | **Yes.** Above-range cutoff output is exactly equal to disabled output in tests. |
| 11 | **Yes.** Forced/all-Cowling solves are finite and nonzero, including attenuation and truncated l=100. |
| 12 | **Yes.** The lower bin is full and the equal/upper bins are Cowling, exactly as specified. |
| 13 | **Yes.** Single- and multi-SEM agree within `1e-10`, including attenuation and crossing chunks. |
| 14 | **Yes.** Each layout/region owns a fresh SparseLU state. |
| 15 | **Yes.** Cadence reuse remains unchanged within a region and restarts at formulation boundaries. |
| 16 | **Yes, within stated scope.** Errors decrease strongly toward 79.6875 mHz; l=100 is nonmonotone at intermediate bins. |
| 17 | **Yes for the measured case.** Switch excess is small relative to natural adjacent-bin change; broader ringing is untested. |
| 18 | **Yes.** DOFs, bandwidth, nonzeros, factorization, solve, and prepared-spectrum time all fall materially. |
| 19 | **Yes.** Tests combine independent topology/weak-form oracles, exact routing checks, end-to-end solves, and baseline serialization. Scientific breadth remains a documented limitation. |
| 20 | **No clearly removable production code found.** The large disabled diagnostics are test/evidence support and do not affect release execution. |
| 21 | **Yes.** The production change is 434 lines of churn across 13 existing files and keeps the formulation choices explicit. |

## Software-design and code-size review

The design extends the existing SEM object and loops without a new class hierarchy, matrix-builder framework, or solver backend. The `forceCowling` flag is a narrow diagnostic control; `false` follows the configured cutoff, and the public default remains backward compatible. Radial and toroidal production code is untouched. The main cost is test/support volume, not production complexity.

- **Production files changed:** 13
- **Test files changed:** 4
- **Production lines added/deleted:** +338/-96
- **Test lines added/deleted:** +1226/-1
- **Total Cowling diff:** 39 files, 237,778 bytes, 3,847 patch lines
- **Total raw diff:** 43 files, 1,107,677 bytes, 16,831 patch lines, dominated by two unrelated binary PDFs
- **Essential production code:** topology, Cowling matrix storage/assembly/accessors, source/receiver layout selection, cutoff parsing, truncation mapping, and single/multi routing; nearly all production additions fall here.
- **Probably removable production code:** none identified with enough benefit to justify changing reviewed code.
- **Test/support code:** substantial but traceable; it contains the independent oracles, routing regressions, opt-in validation/performance diagnostics, saved logs/tables, and status history.

## Testing and performance review

Fresh review-time checks passed: Release CTest 63/63 active tests (two diagnostics disabled), Debug solver API 9/9, Debug Cowling SEM 8/8, both opt-in diagnostics, and `git diff --check`. The fresh performance rerun reproduced all structural counts and the same direction of timing improvement.

Saved degree-20 measurements are full/Cowling: 2490/1661 total DOFs, 2406/1605 reduced rows, 15/10 half-bandwidth, 35,224/19,213 nonzeros, 2.138/0.609 ms factorization including symbolic analysis, 0.341/0.143 ms for four RHS, and 41.645/16.112 ms for the prepared 13-bin spectrum. The measured gain is large enough to justify the compact path, while remaining a local benchmark rather than a general speedup claim.

No further implementation, cleanup, commit, push, PR, or merge is authorized by this verdict. Stop for human review.
