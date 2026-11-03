# Standalone receiver gate v2

This directory preserves a v2 generator/template/output without changing any
v1 files. The generator verifies the frozen v4 and v6 header hashes and extracts
the v4 harmonic and receiver methods plus the v6 receiver and Wigner estimator
methods verbatim. The independent expected receiver entries are based on
`Y^L_nm(0,0)=rho*delta(n,m)` for all `m=-1,0,+1` S/T component signs. `rho` is
computed with `long double` `acos(-1.L)` and `sqrt`, then rounded to `double`;
the harness requires long double to have more significand bits than double.

The v6 oracle gates finiteness and shape, present real/imaginary components
within eight ULP of the independent entries, and exact zero for forbidden
components. The v4 path gates finiteness, shape, forbidden support, and present
component phase/sign. Its amplitude error and the v4-v6 norm difference are
reported as diagnostics only. `v4_status` and `old_new_status` mark scaled
amplitude differences above `2e-14` as `REFERENCE_DRIFT`; these labels are not
gates. The old harmonic
indicator checks cover every valid `n,m` through `+/-2`; present entries gate
positive phase/sign, zeros gate forbidden support, and amplitude drift is
reported without an amplitude tolerance.

The payload gate compares the estimator with exact `GSHIndices<All>` storage
counts for both capped `mMax=min(2,L)` and full `mMax=L` bounds, without
allocating large Wigner tables. Timer lambdas consume output through a volatile
checksum so the calls remain observable; timings are diagnostic. The receiver
degrees are 1, 2, 3, 12, 384, and 1536. Degree zero is excluded from this
nonradial receiver comparison; the radial `L=0` case belongs to the separate
full-kernel gate.

`receiver_gate-v2.cpp` is generated source. `receiver_gate-v1-to-v2.diff`,
`receiver_gate-v2.sha256`, and `provenance.json` bind the wrapper changes and
the exact frozen-source hashes. This harness is source-checked only here; it
has not been compiled or executed in this workspace.
