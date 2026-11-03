# DSpec acceleration candidate: shared packing layouts

Base: exact copy of `candidates/packing` before this refinement; the original packing candidate remains unchanged.

## Change log

- Replaced per-component label flattening/sorting with one exact label-layout deduplication pass per synthesis call. Label hashes only select collision buckets; full label vectors are compared before layouts are shared.
- Built sorted key vectors and record-to-key permutations once for each absolute spin threshold (0, 1, and 2), then reused them across the matching six components. Records with identical `(m,l)` sequences share a permutation. Exact label order and coefficient-slot mappings are preserved.
- Retained the existing coefficient-packing profiler scope and added a separate scope for layout deduplication/permutation construction. Added counts for distinct layouts, threshold union inputs/unique keys, shared permutation entries, and avoided tree lookups.

## Validation

- Static diff inspection confirms only the packing block and `<unordered_map>` include differ from the `packing` base; the matrix contractions and accumulation loops are unchanged.
- No build, remote execution, or scientific run performed. Parent owns compile/parity/timing and resource checks. This refinement shares layouts within each synthesis call; it does not add persistent cross-call layout caching.
