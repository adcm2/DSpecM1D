# OVERALL VERDICT: PASS WITH NON-BLOCKING NOTES

## Scope and reviewed snapshot

Fresh independent Stage 9 review of `dev/cowling` at HEAD
`8f7e7720bb2201bba035ae4abc8e919857a43e32`, against immutable baseline
`889452b6790d2c60afacc156b4c7512990758aee` and the Stage 9 protocol in the
original prompt. The intended dirty-tree scope is present: the Stage 9 test
addition, its implementation-status entry, and validation evidence. There are
no Stage 9 production-source changes. The current test source hash is
`ff21acaa2921f5c219ad507a31b520fd875bd0fe4a7124597da389e95303ab9b`, matching
the README provenance; the PREM input and Release test executable hashes also
match the recorded values.

## Blocking findings

None.

## Non-blocking notes

- The performance figures are local to one degree-20 PREM/no-attenuation case,
  one mesh (`nq=5`, `maxstep=0.005`), one representative frequency, one thread,
  and five warm timed repetitions. They establish the observed benefit for
  this case, not a general speedup guarantee. The README states those limits.
- The 8.777 ms SEM build measurement constructs the mesh and both matrix
  layouts together, so it is not a full-only versus Cowling-only setup-cost
  comparison. The report correctly labels the shared construction scope.
- The opt-in test adds a sizeable but localized diagnostic to the existing
  solver fixture. It is disabled in ordinary CTest, introduces no benchmark
  framework or production API, and uses small local measurement helpers. The
  required diagnostics and existing Stage 8 campaign account for this test
  source growth; I found no Stage 9 simplification that would clearly improve
  readability without hiding the measurement or reducing its evidence.

These notes do not require automatic changes.

## Measurement correctness and interpretation

The source and captured output agree on the structural measurements. Full and
Cowling operators have 2490 and 1661 DOFs (ratio 1.4991); frequency truncation
starts at rows 84 and 56, leaving matrices of size 2406 and 1605. The measured
lower/upper half-bandwidths are 15/15 and 10/10, with 35224 and 19213 stored
nonzeros. The benchmark explicitly checks that the frequency-relative start
indices agree with `startElementSph` and each layout's LTG map, addressing the
earlier global-versus-relative index error recorded in the status history.

Fresh SparseLU computations include COLAMD symbolic analysis and numeric
factorization. Four right-hand sides use the force rows corresponding to the
truncated matrix rows. The test checks finite solutions and relative residuals
below `1e-8`. Its prepared-SEM spectrum measurement uses the same 13-bin band
for each formulation and includes the solver's per-call matrix, source, and
receiver setup while excluding shared model/SEM construction. The saved medians
are 2.138/0.609 ms for factorization, 0.341/0.143 ms for four-RHS solve, and
41.645/16.112 ms for the prepared-SEM spectrum, full/Cowling respectively.
These ratios are consistent with the reported dimension and sparsity changes.

I reran the opt-in performance test with `OMP_NUM_THREADS=1`; it passed and
reproduced all structural counts. Its timing medians varied, as expected for a
short local benchmark, while retaining the same directional result. I also
reran `PreferredSolverApiTests.*`; all 9 active tests passed. The saved parent
logs report 63/63 active Release CTest tests, 9/9 Debug solver API tests, and
8/8 Debug Cowling SEM component tests. `git diff --check` passed.

The saved Stage 8 CSV section after Stage 9 is 89 lines and compares byte for
byte with the original campaign section. This preserves the prior validation
evidence across the Stage 9 addition.

## Simplification and code-size assessment

The accumulated production diff remains limited to the Cowling topology,
operators, input setting, force/receiver layout selection, and single-/multi-
SEM solve routing. The distinct Cowling matrix assembly is necessary because
its radial stiffness term must omit perturbed self-gravity while retaining
background gravity. The solver paths share their frequency traversal and
select the matching maps, matrices, forces, and receivers at the point of use.
I found no clearly beneficial production-code simplification. Stage 9 changed
only the disabled test and its evidence/status records; no code or scientific
results were altered during this review.
