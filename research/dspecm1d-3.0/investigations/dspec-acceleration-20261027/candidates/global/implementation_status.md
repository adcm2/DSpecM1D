# Global GSH and bounded accumulation candidate

Base: `investigations/tunya-profile-20261027/dspec/source-v6-tunya` (copied as the candidate starting point).

## Change log

- Candidate created from the approved frozen source tree. Prior evidence and other candidates remain untouched. No build or physical run has been performed by this implementation agent.

## Current status

Implementation is in progress. Correctness, timing, memory and disk claims are pending independent review and parent-owned Tunya validation.
- Added `global_gsh_tile.hpp`: complete-band inverse GSH synthesis with the accepted local-to-global Wigner rotation convention, complex source/unit-adjoint mixing, separate radial L0 addition, and the existing fluid projection. This follows the isolated global-GSH prototype's `SphericalGrid::InverseTransformation` path. It is source-only and has not been compiled or numerically validated.
- Added `pole_response_spool.hpp`: append-only binary storage for two pole records and separate radial-L0 records per frequency/native side, exact original bin and complex `z`, adjoint mixing and receiver-basis residual. Disk index tables are kept on disk; in-memory indexing is limited to per-key completion flags, side keys and frequency audits. A worst-case spool-size/free-space preflight and byte cap are enforced.
- The consumer now transforms complete angular grids per frequency/native side, sums matching-frequency point kernels in a full-side accumulator, and writes each output chunk once after all frequencies. A first-64-angle prefix remains available for reference comparison; its manifest states that full transforms were performed. Working-memory preflight bounds the requested number of independent side workers to 48, with Eigen's internal thread count fixed to one. The estimate explicitly excludes library, FFT, SEM, and runtime overhead.
- Added explicit native-side selection with an eight-side default cap, refusal when a larger spool lacks an explicit index list, and manifest mapping from output-side index to source element/node/radius/width/fluid metadata. Expanded working-memory preflight to include the twelve field arrays, twelve coefficient arrays, side accumulator, two frame Wigner tables, and the estimated complete angular grid table; FFT/GSH/SEM implementation overhead remains excluded.
- The earlier direct per-point/tile evaluator and its phase cache were removed after static complexity review showed they would scale as O(tile-longitudes × L²) per frequency and should not be used as the physical execution path.
- The separate `candidates/responses/physical_driver.cpp` integration invokes the spool capture and consumer; it exposes a 1–48 outer worker count and output chunk size. That integration and its header copies are owned by the response candidate and must be synchronized to the final global headers.
- All additions remain source-only. No build, physical run, or numerical check has been performed. Independent mathematical/code review and parent-owned Tunya validation remain pending.
