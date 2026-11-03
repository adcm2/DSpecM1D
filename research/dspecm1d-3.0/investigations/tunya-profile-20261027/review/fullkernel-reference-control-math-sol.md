# Sol mathematical recommendation: full-kernel analytic-endpoint control

27 October 2026. Documentation of the completed source/math assessment only. No implementation edits, rerun, build, scientific execution or SSH performed. Parent-provided full-kernel errors have not been independently recomputed here.

## Recommendation

A separate validation-only v4 control replacing the scalar harmonic value at the exact nonradial pole by its analytic normalized Kronecker delta is a legitimate isolation experiment. Preserve v4's original `receiverAtPole` S/T component aggregation. Do not copy the v6 receiver matrix implementation into the control. Keep the existing fixed `2e-13` full-kernel comparison criterion for control versus v6.

The parent reports original full L3 v4/v6 maximum scaled errors of approximately `3.358578e-13` for C and `4.154004e-13` for F, above `2e-13`, with native/profile output bytes exact and unchanged DSpec trace/energy/q-dot-P evidence. Preserve those original results as **failed under their original contract**. Neither the scalar receiver oracle nor a future control pass retroactively makes that comparison pass. No numerical tolerance is relaxed by this recommendation.

## Exact scope of the control

The frozen scalar `harmonic(int degree,int m,int n,double theta,double phi)` helper is called by the old receiver aggregation at `Windowed3DFields.h:536–538`; it is also called at `:881` to produce radial/source `y00`. Consequently a generic theta=phi=0 substitution would affect more than the intended nonradial receiver path.

Use a branch only for **valid nonradial degrees (`degree>=1`) and valid m,n**, with exact `theta==0.0` and `phi==0.0`. Return `rho*delta(n,m)`, with `rho=sqrt((2*degree+1)/(4*pi))`. Leave degree zero and all other arguments on the original verbatim Wigner fallback. Do not use a near-pole angular tolerance, change the global Wigner library, alter receiver component combinations, or modify source/solver/profiling algorithms. Preserve all original frozen/portable sources and create a separately named, hashed validation-only copy and binary.

This independently locates the analytic substitution at the scalar special-function boundary while retaining v4's S/T assembly and subsequent numerical path. The v6 direct receiver formula and the scalar control share the mathematical normalization identity, already independently checked against a higher-precision oracle; the control is an implementation-path isolation test, not a new proof of that identity or of general Wigner accuracy.

## Interpretation of cancellation

The saved scalar gate shows old L3 pole discrepancies of order `3e-16` relative. A larger relative kernel-output discrepancy is mathematically possible: adjoint basis mixing and kernel contractions can combine terms with cancellation, so small operand perturbations need not produce equally small output-relative errors. A maximum scaled error over a channel does not by itself identify the condition number or the responsible cancellation.

Thus cancellation amplification is presently a hypothesis, not an established explanation. An appropriately isolated control that agrees with v6 at the unchanged criterion would provide causal evidence that replacing the old pole endpoint evaluation accounts for the observed original/control difference. It would not prove broad high-degree kernel accuracy or justify changing any physical threshold.

## Required bounded controls and evidence

1. Freeze and review the exact validation-only source delta. Verify that the old S/T receiver aggregation, scalar fallback, l0 path and all other source files remain unchanged. Complete the specified fresh Luna then short Sol review before execution.
2. Check the control's scalar pole identity and assembled S/T receiver matrices against the independent analytic oracle and v6, with shapes, finite values, signs, support and exact zeros. Retain strict new-code normalization checks. Prefer byte equality for L3 control/v6 receiver matrices when the same floating-point rho expression is intended; otherwise record the actual ulp difference rather than assuming identity. Degree zero remains on the untouched original route.
3. Use identical frozen physical inputs, compiler settings, numerical thread settings and complete L3 output domain. Require control/v6 per-channel/column error `<=2e-13`, complete sizes/layout/radial-side metadata, finite values and existing energy-adjoint checks. Keep native/profile numerical exactness requirements. Do not replace the full-grid gate with selected locations or scalar tests.
4. Preserve and compare all three datasets: original v4, analytic-endpoint v4 control and v6. Record the worst C/F error locations, physical/grid coordinates, both values, absolute error and channel scale for every relevant pair. Retain per-channel metrics for all six parameters, not only C/F. Where available, inspect contributing terms before claiming cancellation amplification. Report control/v6 residuals rather than assuming they vanish.
5. Keep the original v4/v6 failure visible and separately named. If the control passes, describe it as **analytic-endpoint control versus v6 passed** and the original difference as evidence of reference endpoint sensitivity/bias. Do not state “v4/v6 passed.” If the control fails at `2e-13`, stop and investigate further without enlarging the threshold or changing additional source paths.

The existing strict old/v6 checker still stops profile-03 before higher-degree profiles if that comparison fails. This mathematical recommendation does not authorize bypassing it, changing the physical `1e-3` threshold, modifying production/frozen implementations, or claiming a complete high-degree global kernel. Any later launch configuration must accurately identify the new control and its reviewed acceptance contract while preserving all earlier failures and partial evidence.
