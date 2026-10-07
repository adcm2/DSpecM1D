# OVERALL VERDICT: PASS WITH NON-BLOCKING NOTES

## Scope and reviewed snapshot

Reviewed Stage8 additions at `dev/cowling`, HEAD `8f7e7720bb2201bba035ae4abc8e919857a43e32`, against the Stage8 protocol and immutable baseline `889452b6790d2c60afacc156b4c7512990758aee`. The working diff at HEAD adds the disabled validation campaign to `tests/test_solver_api.cpp` and its entry to `implementation_status.md`; the accompanying `validation/cowling/` evidence contains the README, raw campaign output, separated tables, inputs, timing, and focused-test log. No production source file changes are part of this Stage8 diff.

## Blocking findings

None.

## Non-blocking notes

- The opt-in campaign adds 243 net lines to `tests/test_solver_api.cpp` (246 added, 3 removed), with another 60-line status entry. That is a substantial diagnostic test body, but it remains local to the existing solver fixture, creates no general validation framework, and is explicitly disabled during ordinary GTest runs. Its calculations and output columns are readable. A future code-size pass could consider extracting the repeated nearest-bin scan, though this is not needed for Stage8 correctness.
- The mesh check samples only l=1 and l=100, two frequency bins, and two meshes. The README correctly describes it as selected resolution evidence rather than general convergence.
- The cutoff check covers one 20.3125 mHz cutoff, one receiver/source setup, and three degrees. It supports adequacy for those cases only; the documented absence of a waveform ringing test and matched MINEOS full/Cowling reference is accurate.

## Physics and convergence interpretation

The measured response comparison uses matched full and Cowling spectra on the same PREM model, receiver, frequency window, and mesh. It reports each output component's complex response, relative complex difference, signed amplitude difference, and wrapped phase difference. The relative-error denominator includes a per-component floor scaled to the largest magnitude in the compared spectra, and phase is omitted below that floor. The requested frequencies are correctly identified as nearest sampled bins (4.6875, 9.375, 20.3125, 40.625, and 79.6875 mHz); the README also records the 1.5625 mHz bin spacing and 2.6526 mHz imaginary cyclic damping.

The interpretation is appropriately bounded. The maximum component-wise relative complex errors fall from 2.47% to 0.0359% at l=100 between the lowest and highest sampled bins, while l=100 is nonmonotone over the middle bins. The l=1 and l=20 samples show progressively smaller errors at higher sampled frequencies. This supports a high-frequency trend in these cases; it does not claim universal convergence or a universal acceptance threshold.

For the hard cutoff, the test asserts exact mixed/full routing at the lower adjacent bin and exact mixed/Cowling routing at the cutoff and upper adjacent bin for l=1, 20, and 100. The cutoff calculation selects 20.3125 mHz, with adjacent bins 18.75 and 21.875 mHz. Reported switch excesses are evaluated against the corresponding full-response scale and are materially smaller than the natural full/full changes for these rows. This is a reasonable local check of the hard switch, with the stated limits on other setups and waveform ringing.

The evidence does not independently re-derive Cowling weak-form physics; Stage8 is a solver-output validation. The protocol-specific distinction between retained background gravity and removed perturbed self-gravity remains a production-code review concern from earlier stages. The Stage8 report makes no new claim about that distinction.

## Provenance and verification

The README records the source baseline, immutable baseline, model hash, test-source hash, three input-template hashes, and hashes of both captured logs. I independently checked these listed source/model/input hashes. The saved full-versus-Cowling, resolution, and cutoff tables correspond to their respective sections of the captured campaign output; the saved test source hash is `afcf98eaefc81bb6f7d5fb80a3c0c01106d322c65cee3d95696dc7003cabd919`.

I reran the recorded campaign and ordinary focused filter with `OMP_NUM_THREADS=1` using the current `build/cowling_final_release/bin/dspecm1d_solver_api_tests` binary. The campaign passed, and its complete CSV sections matched the saved campaign CSV data byte-for-byte. The focused filter passed 9/9 tests. Timing lines differ between runs, as expected. `git diff --check` also passed. The independently reported parent full-gravity regression hash is recorded in the README; I did not rerun that separate baseline regression.

The campaign command is documented in `validation/cowling/README.md`. MINEOS comparison was unavailable through a straightforward matched workflow and no infrastructure was added. Stage9 has not started.
