# Standalone receiver gate

`generate_receiver_gate.py` verifies frozen v4/v6 `Windowed3DFields.h` entries against their SHA-256 manifests, extracts exact function bodies, fills `receiver_gate_template.cpp.in`, and checks each generated `#include` line is complete before writing `receiver_gate.cpp`. The extracted v4 `harmonic` and `receiverAtPole`, v6 `receiverAtPole`, and v6 Wigner estimator bodies remain verbatim in independent namespaces. `receiver_gate-v1.cpp` preserves the first reviewed draft, including its malformed include and v4/v6-only check.

The active harness compares spheroidal and toroidal matrices at degrees 1, 2, 3, 12, 384, and 1536. Its independent analytic matrices use the pole identity `Y^L_nm(0,0)=rho*delta(n,m)`, where `rho=sqrt((2L+1)/(4*pi))`, and explicitly encode all six S/T component channels for columns `m=(-1,0,+1)`. Present real/imaginary components must match within eight ulps; v4 forbidden entries must be at most `5e-15*rho`; v6 forbidden entries must be exact zero. It also checks the old Wigner helper against every valid `n,m` pair through ±2 and reports maximum absolute and scaled old/new and analytic errors.

Payload rows compare the v6 estimator with exact `GSHIndices<All>(L,mMax,n).Size()` counts for both capped `mMax=min(2,L)` and full `mMax=L`, without allocating a Wigner table. Receiver timings are separate three-call medians. Degree zero is intentionally excluded because the old nonradial receiver asks for invalid `n=±1` there; the radial l=0 route remains part of the separate full-kernel gate.

After independent review, build and run this gate on Tunya with `OMP_NUM_THREADS=1` and the Eigen/GSHTrans include paths from the frozen source. The generator's include-line check is a source-level syntax check only; the C++ harness has not been compiled or executed locally.

The Tunya portability copies (`source-v4-tunya/` and `source-v6-tunya/`) replace only five absolute includes in `SensitivityKernels.h`; their `Windowed3DFields.h` files are byte-identical to the frozen headers. The gate's v4/v6 source hash bindings therefore remain anchored to the frozen files and still match the portable receiver bodies.
