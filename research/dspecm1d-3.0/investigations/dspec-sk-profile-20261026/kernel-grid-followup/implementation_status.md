# Native kernel-grid generator status

Scope approved: standalone CSV generator in this directory only. It instantiates `GSHTrans::GaussLegendreGrid::ForBand(L, 2, 2.0, FFTWpp::Estimate)`, validates the legacy angular blocks and per-trace radial metadata, and writes geometry-only rows on that native grid. It does not compute kernel values or change a transform/backend.

The implementation preserves first-seen `(element, radius IEEE-754 bits)` trace order, including same-radius traces from different elements. It strips a trailing CR from the header and rows and filters to `J_Z,A` before validating the remaining columns. Output retains objective, modulus, radius spelling, element, fluid, and `dr_m`; `theta_rad`, `phi_rad`, and `dOmega` come from `GaussLegendreGrid<double, All, All>`. The kernel column is empty and `deltaX_Pa` is zero.

## Change log

- Added `kernel_grid_generator.cpp` with the CLI `INPUT.csv OUTPUT.csv L`, CRLF handling, pre-parse `J_Z,A` filtering, 11-column schema checks for retained rows, finite/range checks, contiguous unique angular-block checks, cross-trace angular and weight checks, per-trace `dr_m` consistency, existing-output refusal, and JSON grid metadata.
- Copied the complete frozen v6 tree to `eigen-candidate/`. The driver `backendId` expression maps `tensor` and `tensor-real` to backend 2; scalar, eigen, and other labels keep their prior mappings and CLI. The dispatch correction below changes only the candidate's `fast_engine.hpp` enum declaration in addition to that driver line.
- Added `build_commands.json` with exact compile-gated commands. The Eigen candidate omits `-DDSFAST_HAVE_CBLAS`, `-lopenblas`, and `EIGEN_USE_BLAS`; it uses no other backend or solver change.
- Generator compile-only succeeded via the shared two-slot gate. Candidate compile-only initially exposed the frozen solver's unresolved `zgbtrf_`/`zgbtrs_` LAPACK symbols; adding `-llapack` resolved them. The candidate still omits `-DDSFAST_HAVE_CBLAS` and direct `-lopenblas`, and defines no `EIGEN_USE_BLAS`. The installed `liblapack.so` has a runtime dependency on `libopenblas.so.0`; that supplies the frozen solver's LAPACK routines, while Eigen handles the selected tensor synthesis multiplication.
- A fresh runtime review found the copied v6 `DotBackend` enum's implicit values were `Scalar=0`, `Eigen=1`, `OpenBlas=2`, while the unchanged driver contract passes backend IDs 1, 2, and 3. Thus the prior tensor/tensor-real ID 2 route selected OpenBlas and would throw with CBLAS disabled. An independent fresh review found the same issue.
- Preserved the previously built binary as `eigen-candidate/eigen_candidate_driver-rejected-dispatch` (SHA-256 recorded in `build_commands.json` and `dispatch_erratum.json`). The isolated candidate now has two intentional source diffs against v6: the driver mapping routes tensor/tensor-real to ID 2, and `fast_engine.hpp` explicitly assigns `Scalar=1`, `Eigen=2`, `OpenBlas=3` to match that contract. No frozen v6 file changed.
- `dispatch_erratum.json` records that earlier comparisons labelled as packed Eigen backend comparisons are not reliable as Eigen measurements; the headline frozen tensor3 route called CBLAS and remains unchanged.
- Corrected candidate compile-only succeeded via the shared two-slot gate, using the same flags and `-llapack`; binary size is 1,556,480 bytes. The rejected binary remains unchanged under its separate name. No numerical run was performed.

Compile and artifact evidence:

- Generator compile: exit 0; binary size 113,928 bytes.
- Prior dispatch candidate binary size: 1,556,480 bytes; preserved as rejected-dispatch evidence.
- Corrected Eigen candidate compile: exit 0; binary size 1,556,480 bytes.

SHA-256:

- `kernel_grid_generator.cpp`: `211961b8a31007071c813c33129836041dbd282b97c1bd54b47e2ae45308e76e`
- `kernel_grid_generator`: `f39d1c246ed16ce9ae7d44106b33192f297c2ff34786fcdf907b40593051a741`
- `eigen-candidate/driver.cpp`: `2d79a7e313d3a9da599a64c240e4feda430b0081389501bd0b89087783907afc`
- `eigen-candidate/fast_engine.hpp`: `e26c438d7e4760e610a10dec5e834f58159eecadea7c5cd87678b38f78629be6`
- rejected `eigen_candidate_driver-rejected-dispatch`: `1eea0b58d3e1429c7fae706540792f36b6c718b6113c09d26a857d0eb4f1bf5b`
- corrected `eigen_candidate_driver`: `ab60c650e58824b94ea0d3cb269c5cb56e1c714ec45052dddd43e6c6c6006320`
- `build_commands.json`: `07b14be17848c8daf1436cf1f879471619bd011a3b57860baa68c92a912c6a5f`
- `dispatch_erratum.json`: `fd20f8e2b2a6a1dc14ad76bdf2af1465ea36b62f0c17643f65f0f9e023be48ef`

No CSV generation, physical calculation, or numerical profiling was performed; runtime validation belongs to the parent phase.

Parent validation: actual GSHTrans native grids generated for L3/L6/L12 (91/325/1225 angles) with exact radial-side metadata preservation,4pi weights and zonal/Fourier quadrature residuals below2e-14. Corrected Eigen candidate passes unchanged strict production-reference and frozen-CBLAS array gates on identical native L3/L12 grids; these are regression-only runs. Physical FFT measurement probe v2 passes true-relative Parseval/adjoint and energy-gradient checks; allthreeFDstep residuals below2.1e-13. Parent TauP PREM ray check selects proposed45N75E receiver to put S1070.84s inside1000–1200s window; this is ray evidence, not a computed DSpec/YSpec wavefield. No physical Earth-kernel timing or convergence claim has been made. Existing production and oldprototype/profilehashes remain unchanged.


Final bounded-phase review: fresh independent Luna and short independent Sol reviews both give qualified PASS; see review/fresh-luna.md and review/short-sol.md. The saved physical configuration is a proposal supported by synthetic measurement identities and TauP ray arrivals. No physical PREM receiver waveform, converged Earth kernel, or physical DSpec–YSpec runtime comparison was computed. DC and upper-bin policy, waveform/energy convergence, and output-memory scheduling remain prerequisites for that next phase. Production files and historical evidence remain unchanged.

Final documentation change: added parent runtime validation to dispatch_erratum.json without altering its historical compile-only record. Added delivery-manifest.json covering current source, binaries, configuration, evidence and reviews.
