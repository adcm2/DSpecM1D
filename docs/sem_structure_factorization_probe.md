# Structure-aware 1-D SEM factorization probe

This benchmark-only experiment asks whether exact element-interior static
condensation can materially beat the current general-band LAPACK backend.  It
does not change production assembly, solver selection, or public APIs.

## Production topology audit

The audit was made against the current implementation, rather than a generic
SEM description:

- `SEMConstructor.h` creates every stiffness, inertia, attenuation, derivative,
  gravity, and displacement/potential coupling triplet inside one element's
  `idxe` loop. The only extra term is the surface-gravity diagonal on the top
  endpoint. There is no matrix entry connecting interiors of distinct
  elements. Model gravity is integrated during material preprocessing, but it
  does not introduce a nonlocal matrix coupling.
- `SEMLG.h` maps GLL endpoints of adjacent radial elements to the same global
  indices. With `NQ=6` (degree five), radial elements have 12 local DOFs: eight
  private interiors and four endpoints. Toroidal elements have six local DOFs:
  four interiors and two endpoints.
- Toroidal modes use only `sem.el() <= e < sem.eu()`, the contiguous solid
  mantle region selected by `SEMConstructor.h`; neighboring included elements
  share one endpoint DOF.
- Spheroidal elements normally have 18 local DOFs: 12 interiors and six
  endpoints, representing radial displacement `u`, tangential displacement
  `v`, and gravitational potential `p` at each GLL node.
- `ltgS` treats a fluid-solid boundary specially. `u` and `p` remain shared,
  while the tangential `v` on the two sides has distinct global indices. Thus
  each fluid-solid interface adds one retained endpoint DOF. PREM produced
  three such split DOFs in every full spheroidal mesh in this experiment.

The executable independently checks the claimed topology against every nonzero
of each tested dynamic matrix. Each non-boundary global DOF must belong to one
and only one element; interior-interior entries must remain in that element;
and an interior may couple only to the endpoints of its own element. All five
checks passed. The exact decomposition therefore exists.

## Exact prototype

For element-private interiors `i` and the union of all endpoint DOFs `b`, the
prototype factors each `A_ii` with pivoted `LAPACKE_zgetrf`, then uses
`LAPACKE_zgetrs` to form

```text
S = A_bb - sum_e A_bi,e A_ii,e^-1 A_ib,e.
```

No inverse is formed. Starting with the assembled global `A_bb` exactly once
is important: it retains all shared-endpoint contributions without needing to
guess how an already-summed diagonal should be split among elements. This is
algebraically identical to assembling element Schur complements.

After sorting retained global indices, each element couples only its left and
right endpoint blocks. The reduced matrix is block tridiagonal. Its observed
scalar half-bandwidth is one for toroidal, three for radial, and six for
spheroidal (the fluid-solid tangential split accounts for the spheroidal width
being six rather than five). It is factored with the current robust pivoted
`zgbtrf` band solver.

For an arbitrary three-column RHS, the prototype computes

```text
b'_b = b_b - sum_e A_bi,e A_ii,e^-1 b_i,
x_i  = A_ii,e^-1 (b_i - A_ib,e x_b).
```

The same local factors are retained for Schur formation, RHS condensation, and
full interior recovery.

## Systems and reduction

All cases use `data/models/prem.200`, production epsilon, attenuation, and the
exact production `A = H - w*w*P + chi*Ha`. Frequencies are catalog-anchored
near-resonance points already used by the robustness campaign. “Ordinary” and
“fluid-solid” spheroidal label the representative resolution/stress roles;
because PREM contains fluid layers, both exercise the exact fluid-solid map.

| System | n / nnz | Elements | Interior | Interface | Original kl/ku | Reduced n / kl/ku |
|---|---:|---:|---:|---:|---:|---:|
| radial 0S0 | 292 / 4,064 | 29 | 232 | 60 | 11/11 | 60 / 3/3 |
| toroidal 0T12 | 81 / 561 | 16 | 64 | 17 | 5/5 | 17 / 1/1 |
| spheroidal 0S12 | 441 / 7,419 | 29 | 348 | 93 | 18/18 | 93 / 6/6 |
| fluid-solid 0S40 | 876 / 14,814 | 58 | 696 | 180 | 18/18 | 180 / 6/6 |
| large 0S40 | 1,626 / 27,564 | 108 | 1,296 | 330 | 18/18 | 330 / 6/6 |

The estimated complex-real-equivalent operation counts for local LU / Schur
formation / reduced band LU range from `2,731 / 6,144 / 816` for toroidal to
`497,664 / 1,119,744 / 240,240` for the large spheroidal case. The complete
per-system counts are retained in the raw TSV.

## Timing result

The policy was one untimed warm-up followed by one measured batch of three
repeats, GCC 15.2 Release, one outer OpenMP thread, and one BLAS thread. Eigen
uses current NaturalOrdering with threshold 1. LAPACK uses current `zgbtrf`
plus the custom pivoted backsolve.

| System | Condensed intrinsic setup | Condensed full solve/recovery | LAPACK factor | LAPACK solve | LAPACK / condensed total |
|---|---:|---:|---:|---:|---:|
| radial 0S0 | 0.0940 ms | 0.0712 ms | 0.1114 ms | 0.0741 ms | 1.123x |
| toroidal 0T12 | 0.0185 ms | 0.0191 ms | 0.0114 ms | 0.0165 ms | 0.739x |
| spheroidal 0S12 | 0.2209 ms | 0.1341 ms | 0.2374 ms | 0.1367 ms | 1.054x |
| fluid-solid 0S40 | 0.4731 ms | 0.2678 ms | 0.4678 ms | 0.2818 ms | 1.012x |
| large 0S40 | 0.9228 ms | 0.5018 ms | 0.8610 ms | 0.5143 ms | 0.965x |

Across all five systems, intrinsic condensed setup is `0.977x` LAPACK
factorization: 2.3% slower. Including arbitrary-RHS condensation and full
recovery, the ratio is `0.996x`, effectively equal to LAPACK factor-plus-solve.
Radial wins by 12.3%, ordinary spheroidal by 5.4%, fluid-solid by 1.2%, while
the tiny toroidal and largest spheroidal systems lose.

Sparse extraction of the element base H/P/Ha blocks is one-time benchmark
scaffolding and is reported separately (0.0276-3.6443 ms per SEM). The timed
intrinsic setup does include the unavoidable per-frequency local
`H-w*w*P+chi*Ha` assignments. The prototype's reduced-band preparation
(allocation, zeroing, and dense-Schur copy) costs 0.0004-0.0149 ms. It could be
removed by assembling Schur updates directly into a preallocated reusable band
workspace. Even granting that future optimization, intrinsic performance is
only equal to LAPACK; including this current scaffold gives `0.986x`.

## Numerical validation and decision

All local, condensed, Eigen, and LAPACK factorizations/solves returned success
with finite solutions. Maximum normalized backward residual was
`6.07e-17`. Maximum full-solution discrepancy was `4.95e-9` versus Eigen and
`5.35e-9` versus LAPACK, both in the deliberately ill-conditioned radial
near-resonance case, whose condensed residual was `1.99e-21`. The remaining
condensed/LAPACK discrepancies were at most `2.15e-12` for toroidal and
`1.63e-10` for spheroidal systems.

Decision: **C — NO WIN**. Exact condensation is feasible and robust, but its
complete intrinsic kernel is essentially equal to current LAPACK, with the
largest representative system slightly slower. It does not clear the requested
5-10% margin needed to justify bespoke production assembly, storage, and
validation. No production implementation is recommended from this result.

The separate theoretical opportunity remains to condense source and receiver
functionals and avoid full interior recovery when evaluating approximately
`x_R^T A^-1 x_S`. It was not implemented and does not change this stage's
full-solution Decision C.

## Reproduction

Exact compiler, link, environment, run, analysis commands, versions, and file
hashes are in
`results/sem_structure_factorization_probe/sem_structure_factorization_probe_metadata.json`.
The raw TSV retains every topology, operation-count, timing, and numerical
field; the analyzer validates the five-case set, complete partition, narrow
reduced bandwidth, statuses, finiteness, residuals, and solution comparisons.
