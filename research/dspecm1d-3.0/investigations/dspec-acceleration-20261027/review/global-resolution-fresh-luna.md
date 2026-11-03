# Fresh static review: global L128/L384 resolution description

**Result: PASS for the read-only descriptive helper and its fixture tests.** `py_compile` passed; all five `test_check_global_resolution.py` tests passed, including receipt corruption, wrong mode/degree, actual-schema comparison, and JSON writing. No scientific comparison, build, or SSH was run.

The helper binds all inputs through completed-job receipts and output hashes. For the angular nodes it now requires receipt-bound `result.json`, `angular_grid.csv`, and `kernel_probe.json`, plus exactly one `kernel-probe` argv token, matching `--lmax`, and the 64-angle probe flag. For the two global inputs it uses the existing global manifest/result validators, then requires the serial Generated/Estimate route, the explicit selected native side order `(59, 100)`, and the expected fluid/solid identities. The serial check also delegates to the thread checker’s legacy serial validation, requiring no inner-thread override and the Estimate planner. Inputs are rechecked after reading.

Geometry is a 720-sample source-to-receiver great-circle arc from (10°N, 20°E) to (45°N, 75°E), with the first and last 5° masked. The result labels it as an approximate bilinear angular interpolation and explicitly carries no acceptance gate or convergence claim. Interpolation uses the receipt-bound theta nodes and periodic longitude indexing; it does not interpolate radius. The binary is read in radial-side-major, theta-major, phi-major, six-channel order. Sampled traces and metrics reject nonfinite values; exact paired zero channels are marked structural zeros, while an all-zero pair across an entire sampled side is rejected. The JSON writer converts NumPy path arrays to JSON lists.

The comparison remains descriptive: relative max and L2 scales are reported without a pass/fail tolerance. No figure or radial interpolation is produced.

## SHA-256

- `tools/check_global_resolution.py`: `025e52a38b66bf4aec365d6e30ac4bc35260545f322997c5f964c3c8ce684020`
- `tools/test_check_global_resolution.py`: `e1986d7c5e9a373fd2aae813d39afc60dd148d7c3780d8bba21499543fbdca2d`
- `tools/check_global_pair.py`: `340efdf0c60b80eea23fa69d642595cf7b989112d9310ed792b21ebe4b36d51a`
- `tools/check_global_threads.py`: `011845c8e2bdf592eb08296306af5227ea519dffd6960854db9c8623faaf3271`
- `tools/check_pair.py`: `6f4d72e7da0489b6b3d0b6fd1d5c15d7ba3b7e6179a14ef453b7e8ccb762e646`
