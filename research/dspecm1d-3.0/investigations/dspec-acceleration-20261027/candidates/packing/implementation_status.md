# DSpec acceleration candidate: packing

Base: exact 12-file copy of `investigations/tunya-profile-20261027/dspec/source-v6-tunya`.

## Change log

- Created an isolated candidate and preserved the frozen source tree and existing evidence.
- Replaced the per-component ordered key set with a contiguous vector sorted and deduplicated into the same lexicographic `(m,l)` order. Replaced repeated per-record tree maps with one cached record-to-layout permutation per component, reused across radial tiles. The complex and real packing paths retain the prior key/order traversal and write the same source value into the same coefficient slot; contraction and summation loops are unchanged.
- Retained existing profiler scopes and counters, and added structural counts for union label inputs, unique keys, permutation entries, and avoided tree lookups.

## Validation

- Static diff inspection confirms only the packing block in `engine_headers/Windowed3DFields.h` differs from v6-tunya.
- No build, remote execution, or scientific run performed. Parent owns those checks; byte parity and packing timing remain to be measured.
