# Analytic-pole control receiver gate

This isolated validation control copies the portable v4 source tree and adds
one early return to its scalar `harmonic(int degree,int m,int n,...)` helper:
for exact `(theta,phi)=(0,0)` and valid degree/orders with `degree>=1`, return
`rho*delta(m,n)`. Its original Wigner fallback remains byte-for-byte in place,
including the `degree==0` radial helper behavior. The existing `receiverAtPole`
body and all other files are unchanged. This is a one-point numerical control,
not a production change or a new optimization.

`generate_control_receiver_gate.py` verifies the old portable-v4, control-v4,
and frozen-v6 header hashes; extracts the old harmonic helper, control helper,
control receiver, and v6 receiver verbatim; and checks that the receiver body,
radial `harmonic(0,0,0,0.0,0.0)` callsite, and nonradial scalar-helper call graph
are unchanged. `control_receiver_gate.cpp` compares the control-v4 and v6
receiver matrices to an independent long-double `rho*delta` S/T matrix at
degrees 1, 2, 3, 12, 384, and 1536. Present components must agree within eight
ULP; forbidden components must be exact zero. It does not call a nonradial
receiver at degree zero and does not allocate a Wigner payload.

The control source tree's twelve-file hash manifest and exact source diff live
beside the copied tree; harness provenance is in this directory. This workspace
only generated and statically inspected the harness; no compile, scientific run, remote action,
or production-source edit was performed. The owner must review the exact native
L=3 control run and its unchanged fixed per-channel comparison gate before
execution.
