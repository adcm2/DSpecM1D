# Independent Luna High review — DSpecM1D Phase 1

**Decision: PASS.** I found no required changes within the approved Phase 1 scope. Review was limited to the four authorized files and compared against the recorded pre-Phase-1 state; existing dirty `SensitivityKernels.h` and `validation/phase0/` were left untouched.

## Findings

No blocking or high-priority issues.

The implementation follows the shared strain path in `Windowed3DFields.h`: radial interpolation is keyed by the resolved `(element, radius)`, source lift is applied only on its owning element, spheroidal/toroidal components are combined before strain construction, degree is the inner accumulation, phase is applied after each degree sum, and frame rotation follows the completed sums (lines 911–1018). Source synthesis and each receiver's adjoint synthesis carry distinct set identities; the coefficient cache is owned by one frequency workspace and is byte-budgeted with recomputation fallback (lines 288–302, 926–984, 1079–1085, 1210–1216). The cross-bin sampling cache uses exact frame/point/radius keys and one shared admission counter for geometry, radial weights, Wigner data, and real angular coefficients (lines 134–153, 781–805, 840–918). Cache accounting is an admission budget, not an aggregate process-memory cap; the implementation status correctly says so.

The kernel route retains the public full-field sampler and skips per-volume-point displacement and r-scaled strain synthesis; the high-level evaluator processes 4096-point blocks and samples receiver displacement once per bin (fields lines 767–770, 991–1016, 1110–1126; kernels lines 99–126). Source contacts remain on the existing separate route. The temporary internal result still uses `Field3DSample` before copying required components into `KernelField3DSample`, so the compact API does not yet remove all unused per-sample storage; bounded blocks limit this overhead. This is nonblocking for Phase 1's stated gate, but sampled RSS remains diagnostic rather than a hard total-memory guarantee.

Test changes add compact/full parity, forced cache fallback, permutations and repeated points, receiver/source frame separation for colocated receivers with distinct weights, and scalar composition behavior. Existing tolerances were not relaxed in the reviewed delta. Frozen field/contact, kernel, full-plane, bounded Born, and timing evidence independently inspected in `parent_validation.json` and `attempt-02/`; the recorded numerical comparisons are well inside their frozen limits and the lmax-12 one-objective median improvement is 4.93x. Source contacts are reported identical. These validate the scoped regression and performance checkpoint, not arbitrary receiver-neighborhood accuracy or convergence.

## Independent verification

- Current field/kernel/test source hashes match `after-implementation-v2.json`.
- `phase1-v2-only.patch` contains only `Windowed3DFields.h`, `Windowed3DKernels.h`, `tests/test_sensitivity_kernels.cpp`, and `implementation_status.md`.
- `git diff --check`: passed.
- Ran the existing focused binary with all common numerical thread environment variables set to 1: 31/31 tests passed.
- Reviewed parent evidence for 177 enabled Release tests (2 existing disabled), 336 field comparisons, 240 contact comparisons, 30 old/new kernel groups, five-modulus full-plane equality, and bounded Born checks; all recorded gates pass.

## Nonblocking limitation

The compact kernel sampler temporarily allocates full `Field3DSample` objects (including fields it omits from the returned compact records), then copies the retained fields. This increases bounded block scratch memory. The measured RSS data and explicit cache admission budgets are useful diagnostics but do not prove a hard aggregate-memory ceiling. No Phase 1 acceptance requirement is failed.
