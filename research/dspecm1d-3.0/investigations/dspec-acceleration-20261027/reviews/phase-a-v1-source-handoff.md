# Phase A v1 source handoff

Source review only. No build or numerical job has been launched.

## Frozen seed and source scope

- Candidate seed is the exact accepted Stage 2 v1 build archive, extracted to receipts/phase-a-v1-frozen-seed/candidate.
- Candidate source delta from that seed: modified engine_headers/Windowed3DFields.h; new phase_a_contract.hpp, phase_a_contract_test.cpp, phase_a_driver.cpp, phase_a_reconstruction.hpp, phase_a_standard_forward.hpp; and pinned formula provenance files under provenance/pinned_formula_sources/.
- implementation_status.md records each source change in both candidate and campaign roots.
- Formula provenance hashes for YSpec SampledKernels.cpp and DSpec FullSpecSingleSem.h, FullSpecMultiSem.h, StartElement.h are recorded in provenance/pinned_formula_sources/source-manifest.json.

## Mathematical and data contract

- The copied workspace hook is diagnostic-only. It exports already projected radial slots after accepted source/lift or receiver projection. It preserves the Dref/J derivative, canonical radius, one-sided element/node identity, and existing solve/cache code.
- The packet stores complex float64 [U,U',V,V',W,W'], identity slot mapping, exact side/material/frame metadata, per-family slot masks, and source-lift tags. Source and weighted-adjoint frames are distinct; l=0 dual isotropic handling is explicit.
- Common six-spin strain and the separate DS radial replay preserve accepted spin signs, exp(+im*phi), local-to-global rotation, normalized radius, and Hermitian kernelVs contractions. There is no YSpec bilinear-kernel substitution or added blanket zeta normalization.
- Fluid reference, DS replay, and common replay share the inherited rr+h-to-rr projection, with unavailable shear zeroed; raw common channels remain recorded.
- Tiny contract fixtures cover l=0 dilation, l=1 translation and rigid rotation, merged spheroidal/toroidal slots, fluid projection, DS/common radial identity, and common-spin Hermitian contraction invariance.

## Fixed job source paths

- A1 uses the frozen accepted Stage 2 L200 trace adjoint_source, solves only bins 9/25/40 sequentially with per-bin cache lifetime, and saves radial plus same-workspace packed-reference packets.
- A2 consumes those same packet bytes; it replays old DS radial formulas with prebuilt complex Eigen angular bases and the common contiguous degree-inner formula. Both use immutable geometry reused across all three bins, and compare all 293 native sides, both sections/roles, and all active-bin groups.
- A3 uses the existing L3/NQ8/.05 complete-sphere audit grid and radial quadrature, frozen baseline adjoint source and frozen plus/minus traces. It accumulates frequencies in the same order and applies the existing smooth bulk profile/integral. Its gate remains the existing central-FD relative error <=2e-4.
- A4 is the 35-bin forwardOnly=true,poleOnly=true route with 8 frequency workers and receiver spectrum/counter gates against the accepted components receipt.
- A5 runs the pinned standard single-SEM forward code twice sequentially on the same SEM/source/receiver and exact common complex-z grid. The first submode starts every radial mode at zero; the second uses existing source-aware truncation. The standard routine's pinned degree-parallel OpenMP decomposition is retained with 8 workers and serial Eigen/inner solves; reports state this decomposition differs from A4's frequency workers. Both submodes require 14,035 factorization calls, 42,035 forward RHS columns, and zero adjoint RHS columns. Untruncated spectrum must match A4 at relative L2 <=1e-11. Timing remains exploratory.

## Contract checklist

- [x] Pinned YSpec/DSpec formula and standard-forward provenance present.
- [x] Export location and metadata preserve accepted projection, radial derivative, frame, lift, and one-sided native side.
- [x] Tiny independent contract fixtures are present.
- [x] A1/A2/A3/A4/A5 implementations are isolated in the Phase A candidate.
- [x] Frequency order and bounded worker scopes are explicit; A5 standard decomposition difference is disclosed.
- [ ] Sol-high mathematical/source review.
- [ ] Parent final source/config review.
- [ ] Native/profile/tiny builds and binary/source hashes.
- [ ] Tiny contract executable and all planned control gates.
- [ ] Five detached numerical processes, completed receipts, and Phase A report.
- [ ] Stop before Phase B; no production integration.


## 2026-10-31 source review amendments

- Tiny translation and rigid-rotation gates compare complex magnitudes and throw runtime errors, so they remain active under `NDEBUG`; rigid rotation uses W=Omega*r and Wprime=Omega.
- A5 passes `referencePeriodND=2*pi/model.ReferenceFrequency()` to the pinned `SpectraRunContext` and `SpecConstants` APIs. The common YSpec complex z grid remains authoritative for every frequency; this corrects the API's reference-period parameter without regenerating or remapping z.
- Sol-high source/math review and parent final bundle review remain pending. No build or numerical job has been launched.
- A5 uses `nskip=10`, matching the pinned `InputParametersNew` prepared-single-SEM default for both sequential submodes. The summary records this cadence and identifies the implementation as prepared-single-SEM; its eight-worker degree-parallel decomposition differs from the ordinary MultiSEM wrapper and symbolic grouping is reported rather than asserted equal.
- A3 now deduplicates immutable angular tables across native radii by exact frame axes, isotropic tag, and order count; the expected unique key set is checked across all 35 bins and its count is reported.
- A5 writes actual chosen per-family/degree/bin/formulation cut start and retained/full row counts for each declared submode, plus nskip and policy; audit counts are summarized and its small capture overhead is disclosed in elapsed timings.
- A3 summary explicitly labels the reused time-shift identity as coming from the frozen L3 baseline `result.json`; the data checker requires this field.
