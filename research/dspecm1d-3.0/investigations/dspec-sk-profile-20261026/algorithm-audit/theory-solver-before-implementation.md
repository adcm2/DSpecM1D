# Solver landscape addendum, 26 October 2026

Read-only inspection; no compilation or production changes.

**Verdict:** add bounded prototypes for equilibrated general-band LU and element-interior static condensation. Neither is theoretically dominated by shared SparseLU on this one-dimensional SEM. Also benchmark NaturalOrdering SparseLU and verified symbolic-pattern reuse, which are low-effort controls separating ordering/analysis overhead from arithmetic improvement. All candidates must share one factorization between forward and Hermitian-adjoint solves and use the already identified minimal receiver bases.

## Verified matrix structure

SEMLG.h interleaves U,V,gravitational-potential for spheroidal fields (three degrees of freedom per ordinary node), two fields for radial, and one for toroidal. SEMConstructor.h forms all off-diagonal couplings inside individual dense element matrices. SEMMatrices.h adds the outer gravity boundary term on a diagonal. Fluid-solid offsets change local indexing and traces but do not introduce long-range couplings. Therefore the assembled H=K_l-z^2 M_l is narrow banded in its physical indexing. For NN=8, the ordinary spheroidal element index span is at most 3*(8-1)+2=23. Measure actual kl=max(i-j) and ku=max(j-i) over nonzero entries for every family/configuration; do not assume that bound through all boundary-index cases or silently drop out-of-band entries.

H is complex symmetric under the present real symmetric K,M contract, and generally indefinite/non-Hermitian. Keep both triangles and standard pivoting. Cholesky or Hermitian positive-definite band routines are not valid substitutes.

## Equilibrated general-band LU

Preserve the existing scaling exactly: d_i=1/sqrt(max_j |H_ij|), A=D H D. Solve A y=D b, x=D y; for an adjoint solve use A^H y=D b, x=D y. D is real positive.

LAPACK ZGBTRF performs pivoted LU for a general complex band matrix. Column-major zero-based packing is AB[kl+ku+i-j,j]=A[i,j], with ldab>=2*kl+ku+1. The extra rows accommodate pivot fill; allocating only kl+ku+1 rows is wrong. ZGBTRS with TRANS='N' solves forward, and TRANS='C' solves the Hermitian-adjoint system using the same factors. These contracts are verified in the primary [ZGBTRF source](https://www.netlib.org/lapack/complex16/zgbtrf.f) and [ZGBTRS source](https://www.netlib.org/lapack/complex16/zgbtrs.f).

For symmetric bandwidth w and order n, factorization is O(n*w^2), storage O(n*w), and each RHS substitution O(n*w); constants include pivot fill. This avoids COLAMD, sparse symbolic analysis, indirect factor storage, and repeated sparse traversal. Existing sparse structure may have fewer entries than its bounding band, so speedup is plausible rather than certain. At w about 23, the vendor implementation may select largely unblocked kernels; do not promise a Level-3 BLAS gain. Time packing, scaling, factorization, both substitutions, and full end-to-end work separately. Match the linked LAPACK integer/complex ABI through a supported interface or verified declarations. Require INFO success plus original-system forward and adjoint residuals and kernel equivalence, not finite outputs alone.

NaturalOrdering SparseLU is a useful intermediary: physical ordering already exposes the band, so COLAMD may not pay. Reuse symbolic analysis across frequencies only for a verified identical sparsity pattern, retaining structural entries even when numeric coefficients cancel. Numeric factors must still be recomputed for each z.

## Exact element-interior condensation

Partition the *assembled, equilibrated* A into global trace degrees of freedom t and element-exclusive interior degrees of freedom i. Then A_ii is block diagonal by element. Use the actual local-to-global maps and retain all interface-specific degrees of freedom, including duplicated tangential fluid-solid traces and surface/centre constraints. Do not merge degrees of freedom just because their physical radii agree.

For each element factor its interior block A_ii,e with pivoted dense LU. Assemble

S=A_tt - sum_e A_ti,e A_ii,e^{-1} A_it,e,
g=b_t - sum_e A_ti,e A_ii,e^{-1} b_i,e,
x_t=S^{-1}g,
x_i,e=A_ii,e^{-1}(b_i,e-A_it,e x_t).

All b here denote the scaled RHS D*b_original; map the recovered solution back with D. Implement products by solves, never explicit inverses. The condensed graph couples adjacent element-end traces, so it is block tridiagonal after correct trace ordering, with variable interface block sizes. A pivoted general-band LU or SparseLU for S is an appropriate first prototype. An unpivoted block Thomas algorithm is not justified for this indefinite dynamic system.

Reuse interior and Schur factors for all forward/adjoint RHS. The adjoint formulas follow from the same block elimination applied to A^H:

S^H x_t = b_t - sum_e A_it,e^H A_ii,e^{-H} b_i,e,
x_i,e=A_ii,e^{-H}(b_i,e-A_ti,e^H x_t).

If f is fields per ordinary node and q=NN, interior blocks have about m=f*(q-2) unknowns; for spheroidal q=8, m=18. Local dense LU and coupling/recovery costs are O(E*m^3) plus O(E*m^2*f), condensed factorization O(E*f^3), storage O(E*m^2); substitutions/recovery additionally scale with RHS count. Compare this with whole-band n approximately E*f*(q-1), w approximately f*q and O(n*w^2). Both are linear in element count at fixed q and have broadly similar polynomial-order scaling, but condensation can have a better constant and dense-kernel locality. It is especially relevant when all interior fields must be recovered for kernels; include that recovery cost rather than timing just S.

The critical limitation is local invertibility: an interior block can be singular or nearly singular at an element dynamic resonance while the full assembled matrix remains nonsingular. Damping helps but does not provide an unconditional condition guarantee, especially with massless gravitational variables. Check local factor success, condition/residual diagnostics, final original-system residuals, and output equivalence; report/fallback to the whole-system pivoted solve on failure. Do not claim unconditional stability from global solvability. Condensation itself is algebraically exact where its eliminated blocks are invertible, but changed elimination order can amplify rounding.

## Bounded prototype evidence

Use the accepted degrees 3/6/12 and frequencies already selected; fixed mesh and no physical convergence sweep. Record n, nnz, kl/ku, trace/interior sizes, factor memory, pack/scale/analysis/factor/substitution/recovery time, forward and Hermitian-adjoint residuals, and full kernel/objective/receiver/contact agreement. Test fluid-solid boundaries, source-side lift, radial/toroidal families, and the worst-conditioned selected bin. Resource rule remains one runtime thread and at most two compiler jobs across the team. Banded LU is the first solver alternative to implement; condensation is a justified standalone diagnostic prototype before ruling it out or promoting it.
