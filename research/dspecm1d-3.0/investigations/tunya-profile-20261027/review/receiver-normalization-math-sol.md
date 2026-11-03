# Independent Sol mathematical review: pole receiver reference drift

27 October 2026. Read-only analysis of the exact receiver gate, saved failed Tunya log and transfer-manifest-pinned GSHTrans source. No source/gate edits, builds, solver runs or SSH performed. The existing failed gate remains failed and must be preserved.

## Classification and recommendation

The saved evidence supports **new closed-form receiver correct against its analytic matrix; old recurrence endpoint amplitude drift**. It does not support a sign, support or new-normalization defect. Requiring the new implementation to reproduce the old diagonal errors within `2e-14`, or requiring the old recurrence to meet an eight-ulp analytic test, is an invalid use of that old implementation as a high-degree oracle.

The proposed new gate version is mathematically legitimate: retain the new implementation's eight-ulp, exact-forbidden-zero, shape and finiteness acceptance criteria; independently strengthen its normalization oracle with long-double arithmetic; retain old amplitudes, old/new differences and old strict-test failures explicitly as diagnostics. This changes the role of an imperfect reference rather than weakening new-code accuracy. Do not invent a larger old-reference acceptance tolerance from the observed maximum. Preserve frozen source and failed gate v1, issue a clearly versioned gate revision, and complete the required review before executing it.

## Exact evidence

The transfer manifest identifies GSHTrans `GSHTrans/Wigner.hpp` at
`/home/adcm2/Documents/Research/Projects/DSpecM1D/build/phase-d-build-cache/build-dspec/_deps/gshtrans-src/GSHTrans/Wigner.hpp`.
Its current SHA-256 matches the transferred digest:

```text
1fd2197a9f37ee463b49ecf0c615300cf98ca53612362ec1c3787ad8397b62c4  GSHTrans/Wigner.hpp
04725ac7f8b347c85bfca527593aced212edec2e12e582bdc1c7b9935d6e0720  receipts/receiver-gate-profile01.log
8b2ea3c4512761cd944631cbf8b380f3fd2157f54ed08f72601f75817b30c232  dspec/validation/receiver_gate.cpp
```

The saved log covers degrees 1,2,3,12,384,1536. All payload tests pass. Every new S/T matrix has maximum analytic-oracle error zero and forbidden maximum zero. Old forbidden entries are also zero. The old diagonal errors grow with degree and differ between n values: at L1536 the n=m=+1 error is `1.20587984042686e-10` absolute, or `7.7113015529885342e-12` after division by rho. The old n=0 and n=±2 diagonal entries also drift. Failures already occur at L12 under eight ulps. Consequently the existing aggregate failure is real under its stated contract, but its cause is an over-strong old-reference condition rather than failure of the new analytic receiver.

The original oracle uses the same double rho expression as the new implementation, so zero new-oracle error proves matrix assembly/sign/support agreement and shared-expression agreement. It is not independent evidence of correctly rounded rho. A separately rounded higher-precision rho closes that narrow shared-arithmetic weakness.

## Why endpoint drift is expected from this active source

The active path does **not** use log-gamma cancellation for the tested n=0,±1,±2 seed rows. `ComputeBlock` seeds at l=|n| using the direct small binomial formula, then advances in degree with square-root/reciprocal tables and floating-point recurrences. The retained `WignerMinOrder` log-gamma reference helper is not used by this path. Therefore attributing the observed failure to lgamma is unsupported.

At theta=0, the exact Wigner identity is

```text
d^l_(n,m)(0) = delta(n,m).
Y^l_(n,m)(0,0) = sqrt((2l+1)/(4*pi)) * delta(n,m).
```

GSHTrans recurses in double for `Wigner<double>`. In the interior loop (`Wigner.hpp:538–552`), for diagonal m=n and theta=0, the exact unnormalized recurrence is

```text
p_l = a_l*p_(l-1) - b_l*p_(l-2),
a_l = (2l-1)*(l*(l-1)-n^2) / ((l-1)*(l^2-n^2)),
b_l = l*((l-1)^2-n^2) / ((l-1)*(l^2-n^2)).
```

Algebra gives `a_l-b_l=1`, so exact seeds `p=1` remain exactly one in exact arithmetic. The implementation obtains these coefficients through rounded products of square roots and reciprocals, then a multiply/subtract; their machine values need not preserve `a_l-b_l=1`. Its final normalization pass (`:582–593`) merely multiplies the already drifted p by rho in a different arithmetic factorization. That final factorization alone could cause a few ulps; it does not explain degree-amplified errors of tens of thousands of ulps. Drift is in the recurred endpoint amplitude, even though it appears observationally as a normalized diagonal mismatch.

For the particularly transparent n=m=0 case, write e_l=p_l-1 and let r_l collect local coefficient/arithmetic error. The exact-coefficient error model is

```text
e_l = ((2l-1)/l)*e_(l-1) - ((l-1)/l)*e_(l-2) + r_l.
q_l = e_l-e_(l-1) = ((l-1)/l)*q_(l-1) + r_l.
```

For zero initial error this yields `q_l=(1/l)*sum_(j=2..l)(j*r_j)`. If `abs(r_j)<=rmax`, the resulting bound is

```text
abs(e_L) <= rmax * (L*(L+3)/4 - H_L),
```

where H_L is the harmonic number. Thus coherent local errors can accumulate quadratically in degree. For fixed small nonzero n the coefficient b_l tends to `(l-1)/l`, giving the same endpoint amplification mechanism. This analytically explains why an L-independent eight-ulp requirement on the old recurrence is unjustified. The observed `7.7e-12` relative drift is compatible with this mechanism. It is not a derived rigorous forward-error bound for the actual compiler/libm sequence: obtaining one would require bounding local errors and their feedback. Therefore `O(u*L^2)` must not be turned into an arbitrary numerical acceptance envelope or advertised as a proven error estimate.

## Required controls for the revised gate

1. Compute oracle rho independently as long-double `sqrt((2*L+1)/(4*acos(-1.L)))`, then round to double. Ensure every arithmetic operand is long double before division, and require `numeric_limits<long double>::digits > numeric_limits<double>::digits` on the test platform. Record the precision/finite result. Use the direct analytic S/T matrices with this rho; do not normalize the new matrix using its own output. Keep eight ulps on every present real/imaginary component, exact new forbidden zeros and zero phase components, finite values and exact shapes. No tolerance weakening of the new implementation is justified or needed by the observed log.
2. Retain all old harmonic-indicator rows, old receiver errors, eight-ulp diagnostic outcomes and old/new differences. Explicitly label amplitude-only departures as reference endpoint drift; never rewrite the preserved v1 FAIL as PASS. Retain old finiteness, shapes, forbidden-support checks at the existing bound, expected pure-real/pure-imaginary phase structure and signs. This guards indexing/convention mistakes without forcing old amplitudes into an invented envelope. A sign/support/non-finite failure is a distinct investigation, not automatically this drift classification.
3. Keep exact capped/full payload checks and the same degree/family coverage. The old ±2 harmonic tests remain valuable structural diagnostics even though the receiver itself uses only n=0,±1. Radial l=0 continues through the actual physical path rather than an unsupported old pole-helper call.
4. Keep the complete L3 physical v4/v6 per-channel/column `2e-13` parity and native/profile exact instrumentation gates unchanged. The mathematical endpoint diagnosis does not authorize relaxing physical-output tolerances, skipping full-grid checks, modifying the global Wigner library or claiming full high-degree kernel accuracy. Higher-degree old/new differences remain measured diagnostics with their provenance, not numerical identity claims.

With those controls, a separately versioned receiver gate may accept a correctly normalized new closed form while honestly reporting old-reference drift. It is a justified validation-contract correction, not proof that subsequent physical gates will pass. No implementation was changed by this review.
