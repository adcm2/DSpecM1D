# Matched 1D source-treatment and elastic-kernel comparison

Run from the DSpecM1D checkout with an absent or empty output directory:

```sh
python3 benchmarks/yspecpp/kernels/run_comparison.py \
  --output build/sensitivity-source-treatment-20261024/za-replacement \
  --build-root /tmp/dspec-source-treatment-20261024-matched --jobs 1
```

The runner verifies the two source checkouts and pinned dependencies, builds
small solver adapters, and writes a self-contained HTML report plus full-
precision CSV/JSON evidence. It performs the coarse/baseline and fine/tight
remote-shell pairs, then fresh source-aware same-element, identical-material
source-node, nine original source-local directions, and one revised Z/A case.
A failed required gate retains the
artifacts and exits with status 1. Matplotlib is required for the plots.

The shared scalar is `dt * sum(j[n] * u[n])`; YSpecpp adjoint samples use the
same spectral weights and Hann window. YSpecpp radial kernels are integrated
with native `r²dr` weights. Sampled curves show `r²K` in the DSpec `dr`
convention. Derivatives are `dJ/dh` for dimensionless fractional modulus
change `h`, with units m s; plotted modulus kernels have units s Pa⁻¹.

## Source-treatment cases

The remote-shell case preserves the earlier homogeneous 2000–4000 km modulus
perturbation and uses the repaired source-aware DSpec mesh. Its source modulus
change is zero, so it remains a useful cross-code control but does not validate
source-local terms by itself.

The same-element fixture is the original unpartitioned homogeneous model where
source and receiver formerly shared one radial element. The fresh run uses the
source-aware mesh and records the new source and receiver cells. The earlier
unaligned result is retained as an unchanged historical reference under the
new output's `inputs/historical_legacy_unaligned/` directory.

A separate fixture adds an identical-material partition at 6321 km, exactly
the 50-km-deep source radius. It checks the ordinary-node path without a
physical material contrast and is not an Earth-interface validation.

The nine unchanged source-local directions apply A/C/F/L/N over the shared
upper-solid support 6200–6371 km (excluding Z/A). The source radius is 6321 km,
interior to that support. The revised Z/A case uses a separate homogeneous deck
with an identical-material partition at 6250 km and applies fractional A only
over 6250–6371 km; the source remains interior and the receiver remains at
6368 km. Density and other moduli stay fixed. Reports label each row by setup.
The prior noise-dominated 6200–6371 km Z/A finite difference is not rerun; its
frozen failure evidence is superseded by this user-requested setup change, not
by a changed gate. Both adapters report positive and negative signed kernel
volume integrals to expose cancellation alongside the net.

The pinned YSpecpp mesh builder adds the source break, then compares later
breaks against the newly growing `edges.back()` value; a receiver above an
interior source can therefore be omitted. The source-local adapter uses a
benchmark-only 1 km maximum element, placing the receiver 47 intervals above
the source for this fixture. This is a quadrature workaround for the pinned
build, not a YSpecpp source change or a production fix. Its exact mesh setting
and rationale are recorded in adapter metadata and provenance.

Both codes use centered fractional steps `0.08, 0.04, 0.02, 0.01, 0.005,
0.0025` for each required source-local direction. The report and
`source-local-fd-sequence.csv` show every D(h), all five adjacent R4 estimates,
and all four triplet R6 estimates. The existing own-kernel gate remains based
on R4(0.005,0.0025), with the same `1e-6` limit; the full sequence is retained
to expose truncation and a small-step noise plateau, not to select a favorable
step.

## Gates and limits

Existing limits remain: metadata `1e-12`, angular frequency `5e-11`, spectrum
L2 and peak-scaled maximum `1e-5` per component, objective `1e-5`, informative
cross-kernel error `1e-3`, and remote-shell kernel-to-Richardson `1e-4`. The
source-local own-kernel-to-Richardson limit is `1e-6` for informative rows.
The informative floor is `1e-6` of the largest absolute derivative for each
observable; weaker rows use scaled absolute error and do not automatically
pass. The source-node diagnostic gates waveforms and objectives. Density,
attenuation, 3D structure, density validation, and broad convergence are
outside scope.

## Bounded 3D kernel benchmark

The separate Phase D 3D fixture, baseline-only command, output guard,
dependency overrides, report, provenance, and YSpec source-contact limits are
documented in [3d/README.md](3d/README.md). Its operational configuration is
`3d/fixture.json`; the reviewed planning protocol remains in `3d/protocol.json`.
