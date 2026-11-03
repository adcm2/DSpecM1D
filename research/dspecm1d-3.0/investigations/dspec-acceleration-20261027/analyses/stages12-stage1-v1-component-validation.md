# DSpecM1D Stage 1 component and validation packet

Status: P1 profile and parity gates passed; N2/N3 native repeats are still running at the captured state snapshot.

- Build verification: [`receipts/stages12-stage1-v1-build-validation.json`](investigations/dspec-acceleration-20261027/receipts/stages12-stage1-v1-build-validation.json)
- Controls snapshot: [`receipts/stages12-stage1-v3-controls-completed/`](investigations/dspec-acceleration-20261027/receipts/stages12-stage1-v3-controls-completed) (all 12 jobs passed; state SHA256 `6a15ec61056c4511d2e9dc3bbd66c690df3f2a62c8f9cd7dbbcd211a919dbdfa`)
- Partial P1 packet: [`receipts/stages12-stage1-v3-profile-packet/`](investigations/dspec-acceleration-20261027/receipts/stages12-stage1-v3-profile-packet) (state SHA256 `820abb3d43fffe30d996f26f45152b96323c9a9ade778afc89e741676623a496`)
- Raw profile: [`investigations/dspec-acceleration-20261027/receipts/stages12-stage1-v3-profile-packet/P1-profile-L200/data/profile.json`](investigations/dspec-acceleration-20261027/receipts/stages12-stage1-v3-profile-packet/P1-profile-L200/data/profile.json) (SHA256 `f84f54999ba9f4383fa5b942e65a7ee63b27f370780cb72bffc1e2ae8c23bb81`)
- Scope-name mapping: [`investigations/dspec-acceleration-20261027/analyses/stages12-stage1-v1-scope-component-map.json`](investigations/dspec-acceleration-20261027/analyses/stages12-stage1-v1-scope-component-map.json) (SHA256 `1639270792d11a9c1cbf47cb489737c3da825d7fec8bff2424fa62d505c594ca`)


- Raw profile rows CSV: [`investigations/dspec-acceleration-20261027/analyses/stages12-stage1-v1-profile-raw-rows.csv`](investigations/dspec-acceleration-20261027/analyses/stages12-stage1-v1-profile-raw-rows.csv) (SHA256 `ccba721ffe6362e3bac0b9ea967f61530677fa93beac39db94c6830e502dced9`)
- Raw profile counters CSV: [`investigations/dspec-acceleration-20261027/analyses/stages12-stage1-v1-profile-raw-counters.csv`](investigations/dspec-acceleration-20261027/analyses/stages12-stage1-v1-profile-raw-counters.csv) (SHA256 `fef36dbd527849ffbaa3e37e03a4f81adf41c0207d60852e1c6809b62038424b`)

## Build and controls

Both build variants returned 0. Native binary SHA256: `8f5b47b7e85920b41add04047d268676fed9dfb3cf1d3674d816b2a25ee4e46e`. Profile binary SHA256: `39c60c2cdb3815d5944088d2b0bd42879d0051d62b65df72ab1b88fb5bda403f`. Build supervisor time: `118.670666 s`.

Controls passed candidate scope, scientific parity (C1/C2/C3), profile accounting, and the L3 finite-difference gate. Finite-difference result was `{fd["finite_difference_seconds"]:.12f} s`; integrated kernel was `{fd["integrated_kernel_seconds"]:.12f} s`; relative error `{fd["relative_error"]:.3e}` (limit `{fd["limit"]:.1e}`).

## L200 timing snapshot

| Run | Result elapsed (s) | Supervisor job elapsed (s) | Forward (s) | Kernel (s) | Tensor validation (s) | Workers |
|---|---:|---:|---:|---:|---:|---:|
| N1-native-L200 | 51.840044503 | 52.480971456 | 10.626435160 | 40.678408821 | 0.515909075 | 8 |
| P1-profile-L200 | 52.164914599 | 52.831545025 | 10.706378212 | 40.920608015 | 0.518276328 | 8 |

The P1 raw profile records a `52.164586220 s` pipeline and the same master-exclusive partition with residual `0.000000000 s`; it contains 150 rows and 100 counters. Snapshot overhead is `0.001496908 s` and report serialization is separately named outside the pipeline.

## Largest master-wall exclusive scopes in kernel production

| Scope | Analysis component | Exclusive (s) | Inclusive (s) | Calls |
|---|---|---:|---:|---:|
| `tensor_angular_gemm_complex_and_packing` | tensor_angular_synthesis | 21.861582 | 21.861582 | 231980 |
| `radial_interpolation_and_strain_numerator_build` | radial_sampling_reconstruction | 9.926802 | 9.926802 | 20510 |
| `source_frame_geometry_phase_spin_rotation_and_scatter` | source_receiver_frame_geometry_rotation_scatter | 2.359845 | 8.821923 | 70 |
| `modal_degree_parallel_region_wall` | parallel_degree_region_wall | 2.143165 | 2.143165 | 35 |
| `receiver_frame_geometry_phase_spin_rotation_and_scatter` | source_receiver_frame_geometry_rotation_scatter | 1.509045 | 5.302881 | 70 |
| `workspace_synthesis_vector_construction` | tensor_angular_synthesis | 0.892882 | 1.147600 | 35 |
| `tensor_coefficient_packing_and_key_reorder` | tensor_angular_synthesis | 0.879114 | 22.455780 | 840 |
| `tensor_phase_spin_rotation_and_scatter` | tensor_angular_synthesis | 0.325506 | 0.325506 | 140 |
| `receiver_basis_mix_and_unit_adjoint_projection` | receiver_basis_projection | 0.252446 | 0.252446 | 14000 |
| `tensor_coefficient_layout_dedup_and_permutation_build` | tensor_angular_synthesis | 0.161275 | 0.161275 | 140 |
| `driver_kernel_pipeline_including_diagnostics` | diagnostic_validation | 0.091538 | 41.438884 | 1 |
| `tensor_radial_l0_scatter` | radial_sampling_reconstruction | 0.087368 | 0.087368 | 70 |

The two degree-region wall scopes account for 8.087119 s during initial frequency/capture and 2.143165 s during kernel production. Worker elapsed sums remain a separate overlapping measure; see the JSON for the leading rows by phase/domain.

## Runtime label caveat

The profile JSON retains every raw `scope_name`, `logical_component`, and timing value. Its runtime component classifier uses broad substring matching, so names containing `io` inside words such as “rotation” or “projection” can be labeled `output_io`. The separate mapping file groups by explicit scope names. Do not sum or interpret the raw `output_io` label as physical output work. No profile rows or binaries were edited for this analysis.

## Validation limits

P1 vs N1 and frozen v9 byte-parity checks passed for parameters, trace, spectrum, components, nodes, along kernel, and along angles. The L3 control has no across output by its result metadata. One N1 and one P1 observation are not repeat means. N2 and N3 remain required before Stage 1 is complete.
