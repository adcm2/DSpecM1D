# Fresh independent Luna review: qualification checker

27 October 2026. Read-only static audit of `tools/qualification_checks.py` against its README, status entry, and actual frozen DSpec v4/v6 and YSpec v6 physical drivers. No implementation or frozen source was edited. No scientific calculation, build, or SSH/remote run was performed. I did not treat the earlier synthetic-fixture result as source-schema evidence.

## Decision: HOLD current checker revision

The checker covers the main result families and its expected field names match the actual frozen drivers reviewed. Three parity/gating issues prevent accepting this revision:

1. **DSpec native/profile kernel parity is not exact.** `run()` passes `exact=False` for both DSpec native/profile kernel comparisons, and `compare_job_outputs()` also applies the same tolerance to their trace and spectrum tables. The native/profile outputs are intended to come from the same frozen source and physical inputs; this gate should require exact bytes/values for the kernel and exact numeric values for the tables. The current tolerance can accept unintended instrumentation-build drift.
2. **Cross-version DSpec kernel scaling pools all six channels.** `compare_arrays()` flattens the kernel file and divides the maximum difference by one maximum across all channels. `kernel.bin` is point-major with channels `[A,C,F,L,N,kappa]`; a large channel can therefore hide an error in a smaller channel. For the v4/v6 comparisons, compute and enforce scaled error per channel (retaining byte-exact reporting as a separate fact).
3. **Trace/spectrum comparisons pool unlike columns.** `compare_numeric()` scales a whole table by its largest value. The time/bin/frequency columns can dominate the physical amplitude columns. Compare each physical value column with its own scale and check the coordinate column separately (exactly where the same implementation should reproduce it; explicitly defined tolerance for cross-version coordinates if needed).

## Schema and gate-source audit

The six expected job names agree with the plan. DSpec v4 and v6 `physical_driver.cpp` write `result.json` with `degree`, `mode`, `energy_m2_s`, and `energy_adjoint_qdotp`; `kernel_domain.json` writes the fields the checker reads: `selected_radial_trace_count`, `available_radial_trace_sides`, `full_target_point_count`, `actual_sampled_point_count`, `actual_kernel_binary_bytes`, `n_theta`, and `n_phi`. Their radial CSV header is `radial_index,element,node,radius_m,element_width_m,fluid`, and the six-channel kernel file layout is point-major `[A,C,F,L,N,kappa]`. The DSpec q dot P(B raw)/(2E) calculation uses the actual `energy_adjoint_qdotp` and `energy_m2_s` fields and enforces `1e-10`.

The frozen YSpec v6 `physical_yspec.cpp` writes `result.json` with `solver`=`YSpecpp`, `lmax`, `mode`, `kernel_point_count`, `kernel_radial_blocks`, `kernel_n_theta`, `kernel_n_phi`, `kernel_output_file`, `kernel_baseline_trace`, and `kernel_sampling_input`; the checker uses those exact keys. It emits `kernels.bin` in angle-major, radial-node, `[A,C,F,L,N]` order. In kernel mode it reads the supplied baseline trace and radial sampling file, and does not emit its own trace or spectrum CSV. The checker treats the baseline as a finite input only and does not claim it is a YSpec trace result. Its output `profile_stages` JSON and its single `MATHPROFILE` record use the `YSpecpp.physical_kernel` root with objective 0 and `physical_case_total`, matching the frozen adapter. YSpec native/profile stable result metadata is compared after excluding timing/RSS/profile-stage fields.

Other static checks look directionally sound: finite JSON/CSV/kernel values are enforced; the complete DSpec angular target count must equal the sampled count; binary byte counts derive from points × channel count × 8; all 179 radial trace rows and exact headers are required; and DSpec native/profile and v4/v6 radial metadata plus YSpec shared radial inputs are compared. The cross-code parameter and fluid-shell kappa comparisons remain explicitly descriptive.

The README accurately states that YSpec's shared baseline is only finite-checked, but it does not document the three parity deficiencies above. I did not run the end-to-end synthetic fixture because the current audit already exposes these source-level gate failures; no numerical evidence should be inferred from the status note's prior fixture.

## Disposition

Keep this revision on HOLD until DSpec native/profile parity is exact, v4/v6 kernel parity uses per-channel scales, and trace/spectrum amplitude comparisons use per-column scales. Then review the changed checker revision and its small negative/positive safety probes before treating it as qualification evidence. These checks qualify only the specified L=3 outputs; they do not establish global cross-code accuracy.
