# Stage 1 diffracted-wave window pilot

This is a deliberately small forward-solver smoke only. It runs the pinned
normal `generate_seismogram` multi-SEM path for elastic, no-attenuation PREM
with full gravity, once for toroidal modes and once for spheroidal modes. Both
cases use a 300 km equatorial source and three 3 km-depth equatorial receivers
at longitudes 95, 130, and 165 degrees. The target wave period is approximately
10.5 s (the 80--113 mHz half-octave band); the absolute record is 10 minutes
with 2 s sampling, toroidal degree range 284--292, spheroidal degree range
145--153, and an 80--85--108--113 mHz filter. The model is
`data/models/prem.200.isotropic`.

The frequency/mesh choices are a bounded proxy: they are selected to make the
short window computationally usable and are explicitly unconverged. Finite,
nonzero outputs are the only acceptance checks. This does not add taper code,
change production Cowling behaviour, or claim a scientific result.

After the Release build, run one case at a time with one OpenMP worker:

```sh
python3 validation/diffracted_wave/run_pilot.py \
  /space/adcm2/research/builds/dspecm1d-diffracted-wave-window/20261118-stage1/bin/generate_seismogram \
  /space/adcm2/research/runs/dspecm1d-diffracted-wave-window/20261118-stage1
```

The driver requires a fresh, nonexistent run directory and materialises inputs
and stores outputs, logs, hashes, and timing receipts under it, leaving raw
outputs outside the source tree. Redirect the launcher stdout to a sibling log
so the shell does not create the run directory before the driver starts.
It validates ten semicolon-delimited columns (time plus three components for
each receiver), finite values, a 2 s interval, a maximum time of 602 s (the
writer emits one extra sample before stopping), and family-appropriate
excitation: each family has a nonzero intended component at every receiver;
toroidal vertical is negligible, while spheroidal vertical is nonzero. One
unused horizontal component may be exactly zero on this equatorial geometry.
`summary.json` records the source HEAD, model/input/output/executable SHA-256
values, per-case verbose BusyBox-time RSS receipts (Tunya has no `/usr/bin/time`),
and the unconverged-smoke interpretation. Each solver command is wrapped by
GNU `timeout --kill-after=5s 300s`; OpenMP and BLAS worker limits are pinned to
one.
