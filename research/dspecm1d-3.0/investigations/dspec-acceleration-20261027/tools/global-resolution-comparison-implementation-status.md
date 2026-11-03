# Global resolution comparison helper

Implementation note: added `check_global_resolution.py` as a read-only JSON
reporter for the receipted L128/L384 serial Generated/Estimate global jobs. It
uses the exact receipt-bound theta grids, explicitly checks the receipted full-
grid command and selected native fluid/solid sides, then reports bilinear
great-circle trace differences with the fixed 720-point path and 5-degree
endpoint mask. The output is descriptive and has no acceptance gate.

Tests added: `test_check_global_resolution.py` exercises periodic-seam and
theta interpolation, path masking, structural-zero handling, and input-layout
rejection. Validation is limited to Python syntax and these in-process math and
schema tests; no remote execution, build, scientific run, or figure creation
was performed.

Follow-up implementation note: bound the global job's receipted argv to the
required degree, serial worker count, selected native sides, Generated storage,
contacts-off, pole-only and full-forward settings. Corrected the structural-zero
test to distinguish paired exact zeros from a valid zero-versus-nonzero channel;
the final focused test run passes all four cases.

Review follow-up: exact pre-review checker and test bytes are preserved under
`preserved-global-resolution-v1/`. Theta-grid references are now validated as
receipt-bound 64-angle `kernel-probe` jobs with a full theta CSV and matching
probe sidecar. Receipt hashes are checked before reading result/probe metadata;
the global jobs also pass the legacy serial/Estimate policy validator. The
great-circle arrays are converted to JSON lists at report construction. Added
an end-to-end synthetic receipt fixture that exercises `compare()` and the JSON
writer, including rejection of a changed theta-reference result hash.
The serial-generated unit fixture uses the production three-part verified-job
receipt shape. Final syntax compilation and all five focused tests pass.
Fresh-review follow-up: the theta reference's command receipt must contain
exactly one `kernel-probe` mode token and exactly one `--lmax` matching the
expected degree. The end-to-end fixture uses the actual mode token and rejects
both a `kernel` invocation and mismatched degree.
