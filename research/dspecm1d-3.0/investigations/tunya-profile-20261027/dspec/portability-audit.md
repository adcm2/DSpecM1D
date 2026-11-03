# Tunya source portability audit

Read-only source audit dated 2026-10-27. The new `source-v4-tunya/` and `source-v6-tunya/` trees are full copies of their frozen source trees, with only the five absolute `DSpecM1D/src/...` directives in `engine_headers/SensitivityKernels.h` rewritten to angle-bracket include-root-relative paths. Each exact per-version diff and complete per-file old/new hash manifest is recorded alongside this note and in `portability-provenance.json`.

An audit over every C/C++ header and source file in both portable copies found zero absolute include directives. The portable `Windowed3DFields.h` files remain byte-identical to their corresponding frozen originals, so the receiver-gate's extracted-body source hash bindings are still valid. Original frozen manifests verify; no build, scientific execution, solver/Earth run, or remote activity occurred.
