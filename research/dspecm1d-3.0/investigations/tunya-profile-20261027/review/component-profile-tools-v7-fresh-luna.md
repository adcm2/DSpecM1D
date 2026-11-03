# Fresh Luna review: component-profile tools v7

**Verdict: PASS for the narrowly scoped material-aware descriptive comparison.** This review does not certify campaign completion or establish global cross-code accuracy.

## Scope and identity

- Reviewed `tools-v7/collect_component_profile.py` against the preserved v6 collector. The delta is limited to material/radius-aware descriptive kernel comparison and its input reshaping/metadata checks.
- `tools-v7/finalize_profile.py` is byte-identical to the reviewed v6 finalizer (SHA-256 `0e347c21adb38034b9f697da1794db6c9cf44d7e43f7dff242196551f24ca195`).
- Frozen v7 collector SHA-256: `894a9185530c34ff67c2863673884eb5d82b55ebbde6237e92f7531359a1d13b`.
- Synthetic test SHA-256: `ddb2d02b3f99d853f93bbb7aa1af25051ab80580f0b152235ded1722b6363a64`.

## Findings

The collector checks the expected DSpec and YSpec binary layout strings, requires a state-listed and hash-verified DSpec radial-node table with 179 entries, and checks the 64-angle × 179-side point count before reshaping both binaries as angle-major, then radial-side, then channel. For solid sides (`fluid == 0`), it compares A/C/F/L/N channel-by-channel. For fluid sides (`fluid == 1`), it compares DSpec kappa with the explicitly labelled YSpec A+C+F sum. The radial subsets use inclusive radius bounds in metres corresponding to 4,200–6,200 km for solid and 2,200–3,200 km for fluid. Maximum-difference records map the selected local index back to the original radial index and report angle index, element, node, radius, material flag, and both values.

The revised synthetic fixture exercises the fluid channel convention and confirms a deliberately mismatched solid A channel remains visible in both all-solid and regular-solid results. The full synthetic suite passed when rerun (`component profile synthetic checks passed`).

The actual immutable profile-06 snapshot also passed collection and returned **INCOMPLETE**, as expected: four YSpec L1536 jobs remain pending. L128 and L384 produced descriptive outputs with 64 × 179 = 11,456 points; region counts were 125 solid sides, 54 fluid sides, 55 regular-solid sides, and 22 regular outer-core sides. L1536 remains pending, not inferred from other runs. The all-solid scaled maxima can exceed one, which is consistent with the implementation’s explicit descriptive-only status and absence of a tolerance or pass claim.

The v6 acceptance checks remain in place; this delta adds no scientific acceptance threshold. The output labels the kernel comparison as a contacts-off, volume-only, first-64-angle diagnostic and explicitly disclaims a tolerance/pass claim. No source, build, SSH, or scientific computation was performed; checks were static review, local synthetic tests, and data-only collection from the existing immutable snapshot.
