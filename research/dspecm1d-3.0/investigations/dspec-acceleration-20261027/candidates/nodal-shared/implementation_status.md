# DSpec acceleration candidate: shared packing plus batched native recovery

Base: `candidates/packing-shared` with the nodal recovery change applied. The original `packing` and `nodal` candidates remain unchanged.

## Change log

- Includes the shared-layout packing refinement documented for `candidates/packing-shared`.
- Includes the fixed-size planned recovery, cached family-specific local-to-global layouts, fluid-safe toroidal handling limited to the elastic element range, cache byte accounting, structural counters, and per-coefficient-build reuse of one-sided nodal r-strain plans documented for `candidates/nodal`.
- Patch tooling left `engine_headers/Windowed3DFields.h.orig` in the shared and combined candidate directories as pre-application snapshots. A cleanup attempt was rejected by the automatic approval reviewer because it would delete evidence; both files are preserved.

## Validation

- Static diff inspection confirms the intended combined source is limited to the shared packing and nodal recovery changes in `engine_headers/Windowed3DFields.h`, plus the shared-layout include.
- No build, remote execution, or scientific run performed. Parent owns compile/parity/timing and resource checks.
