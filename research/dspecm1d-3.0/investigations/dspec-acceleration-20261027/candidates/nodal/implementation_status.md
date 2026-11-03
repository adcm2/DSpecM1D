# DSpec acceleration candidate: batched native recovery

Base: exact copy of the completed `candidates/packing` snapshot, including its 12-file source tree.

## Change log

- Created this isolated candidate from packing; the packing snapshot remains unchanged.
- Replaced planned interpolation's dynamically allocated 6x1 matrix with a six-value stack result. The original quadrature-node loop and per-component accumulation sequence are retained.
- Cached radial, spheroidal-u, spheroidal-v, and solid-element toroidal local-to-global indices in each radial sampling plan. Fluid and out-of-elastic-range elements do not query toroidal indices. Cache byte accounting now includes the index arrays. Added counters for cached layout nodes/lookups, planned recovery calls, and avoided six-by-one allocations.
- In the optional exact nodal r-strain bridge, build node plans once per coefficient build and reuse them across modal rows, keeping the requested element fixed at interfaces.

## Validation

- Static caller inspection confirms both `interpolatePlanned` consumers use the new fixed-size result through `mode`; no other caller expects a matrix.
- Static diff review confirms the coefficient and nodal arithmetic loops remain in their original order. No build, remote execution, or scientific run performed. Parent owns compile/parity/resource checks.
