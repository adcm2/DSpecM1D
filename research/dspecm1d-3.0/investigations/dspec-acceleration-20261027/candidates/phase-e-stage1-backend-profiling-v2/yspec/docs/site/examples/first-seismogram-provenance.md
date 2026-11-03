# First-seismogram illustration provenance

The current illustration `../assets/first-seismogram.svg` is the Phase 2C
Bolivia render from the 18,001-sample, one-second trace. Its current fixture,
raw-output and SVG hashes are recorded in the “Current-source Bolivia
candidate” and “Phase 2C rendered assets” sections below. It is an explanatory
trace, not a scientific reference result.

## Historical Phase 1 illustration (retired 91-sample render)

The records in this section describe the earlier 91-row, 20-second illustration;
they are historical and are not the provenance for the current SVG.

- Code revision: `9552d2bde73b11551027fb829b17d7b4813253d9` (`docs/phase1-user-guide`; documentation changes were uncommitted during the run).
- Platform/compiler: Ubuntu 24.04; GCC 13.3.0 (`g++-13`).
- Build: CMake 3.28.3, Release, `YSPEC_BUILD_TESTS=OFF`, OpenMP enabled; two build jobs. Five dependencies were selected from the existing read-only cache via explicit source-directory overrides; planetmodel was fetched at its pinned SHA. This was not a cold-network configure.
- FFTW: 3.3.10, with single-, double- and long-double components found.
- Dependency revisions: NumericConcepts `e03a2889d43a03660d0879f46dda1406d90c9517` (v1.0.0); Interpolation `d5a1711ca6a89fc8ad8148fe13212398fe569bc0` (v1.1.2); FFTWpp `cb152c0e85f3822493c2d9e85a56307ac7b4f1c0` (v1.0.1); GaussQuad `7d4786d192f68077cd99a34dcdf1379ed1eeac9a`; GSHTrans `8dd116df8df4222fe9b10a40c7213d34f2a274fb`; planetmodel `fa6e230fde25152d40f27d9839de44c4f3182a24`.
- Runtime command: `OMP_NUM_THREADS=2 ../guide/bin/yspec input.json`, from `build/first-seismogram/`.
- Measured solver runtime: 0.50 s. Progress reported 26 frequencies, one source and one receiver.
- Input: `first-seismogram.json`; SHA-256 `72a30885258026dc1662f7713ac3674a37985bf045165965d41dd308c1f219a4`.
- Output: `first-seismogram.out.1`; 91 rows, four columns, times 0–1800 s in 20 s steps; SHA-256 `606e7c0eece9336b83e3247c97587d5e560c871fef84c90de9d43fe9dc4b5ece`.
- Column meanings: time (s), up displacement (m), north displacement (m), east displacement (m). The vertical peak absolute amplitude is `4.344679e-6 m`; all values were finite and the trace was nonzero.
- Plot library: Matplotlib 3.6.3. The plot script sets a fixed SVG hash salt for reproducible element IDs and removes only trailing horizontal whitespace while preserving line endings. Plot SHA-256: `925f2609ad95bab8d498edcdca73a3ba1add470780bb8c267ff7f7e5878abe61`.
- Plot command from repository root: `python3 docs/site/examples/plot_first_seismogram.py build/first-seismogram/first-seismogram.out.1 build/first-seismogram/first-seismogram.svg`.

## Frozen-main reproduction

The following run was made during the 21 October 2026 documentation campaign
from the exact current manual source tree. It preserves the earlier provenance
above rather than relabeling that historical run. The saved trace and SVG
bytes match the earlier files; byte identity alone does not establish that the
earlier provenance applies to this run.

- Source: `0bc9417fd7a7b38f0987e69969a77ac226522faa` (frozen main selected for this campaign).
- Platform/compiler: Ubuntu 24.04; GCC 13.3.0 (`g++-13`); CMake 3.28.3; FFTW 3.3.10; OpenMP 4.5.
- Build: Release, `YSPEC_BUILD_TESTS=OFF`, OpenMP enabled, FFTW, NetCDF disabled; at most two build jobs. Pinned source dependencies were reused from this campaign checkout's build cache. This was not a cold configure.
- Dependency revisions: NumericConcepts `e03a2889d43a03660d0879f46dda1406d90c9517`; Interpolation `d5a1711ca6a89fc8ad8148fe13212398fe569bc0`; FFTWpp `cb152c0e85f3822493c2d9e85a56307ac7b4f1c0`; GaussQuad `7d4786d192f68077cd99a34dcdf1379ed1eeac9a`; GSHTrans `8dd116df8df4222fe9b10a40c7213d34f2a274fb`; planetmodel `35be4ee54b3d84df1483e647878c5446a886c90b`.
- Input and output: the fixture hash remains `72a30885258026dc1662f7713ac3674a37985bf045165965d41dd308c1f219a4`; the 91-row trace SHA-256 is `606e7c0eece9336b83e3247c97587d5e560c871fef84c90de9d43fe9dc4b5ece`. Runtime was 0.17 s with two OpenMP threads; peak resident memory was about 8.1 MiB.
- Plot: generated from the saved trace with the command above using Matplotlib 3.6.3. SHA-256 remains `925f2609ad95bab8d498edcdca73a3ba1add470780bb8c267ff7f7e5878abe61`.

## Current-source Bolivia candidate (Phase 2A)

This section records the new candidate without replacing the historical runs
above. The canonical fixture follows the user-directed DSpecM1D Figure 2 case:
PREM TI with the ocean, no attenuation, full gravity, displacement; source at
647.1 km with the documented full moment tensor and geographic coordinates
(-13.909919806613164°, -67.25°); receiver at (80.06588632614357°, 0°), depth
3 km; degrees 0–300; 1.5–21 mHz; 300 minutes at 1 s; fixed four-corner filter
[1.9, 2.0, 20.0, 20.1] mHz. The model deck's elastic columns match the reference
model; its SHA-256 is recorded in the machine-readable evidence manifest. The
YSpec++ source coordinate is already geographic. It maps to geocentric latitude
-13.82° at the source and 80° at the receiver; do not apply a second latitude
conversion.

The current-source baseline was run with the source-verified `yspec` binary
identified in [first-seismogram-evidence.json](first-seismogram-evidence.json).
It produced 18,001 four-column samples at every integer second from 0 through
18,000. The unmodified trace is retained at
`/tmp/yspecpp-phase2a-bolivia-20261023/base60/bolivia.out.1`; its SHA-256 is
`196c82b1f2ef926146fe5b32e219891374a4af632e75888e52353214f91d7021`. Runtime
was 52.73 s with at most two OpenMP/BLAS threads. The canonical JSON, model,
binary, inputs and outputs for bounded variants are hashed in the linked
manifest; raw runs are machine-local and are not added to the documentation
site.

The declared trace gates require at most 1% peak-normalized maximum and relative
L2 differences per significant component for two successive degree refinements,
plus a 1% duration check. Halving the time step and tightening solver tolerance
by tenfold are each targeted at 0.2%, and that joint refinement is checked too.
The L300, L450 and L600 serialized traces have identical hashes at the output's
seven-significant-digit precision, so both degree-refinement comparisons pass
for the saved text values. This equality does not bound differences below the
serialization precision. Against L300, tightening tolerance to 1e-6 changes
component maxima by 0.08349% (Z), 0.04051% (N), and 0.05111% (E); halving the
time step gives exactly the same every-other sample values on the 300-minute
record. The combined tighter-tolerance/half-step check has the same maximum
changes as the tolerance-only check and passes the declared 0.2% target.

The duration gate does not pass. On common samples over the full original
record, the 600-minute run differs by 1.43017% (Z), 2.19238% (N), and 2.14357%
(E) in peak-normalized maximum; relative L2 differences are 2.34589%, 3.60058%,
and 3.31844%, respectively. The 300-minute user-selected setting is retained;
its duration/frequency-grid convergence remains unresolved. The baseline
300-minute FFT uses N=32,768, a 32,768 s period, df=0.0305176 mHz and 642 bins
from 1.46484 to 21.02661 mHz. The 600-minute run uses N=65,536, a 65,536 s
period, df=0.0152588 mHz and 1,281 bins from 1.48010 to 21.01135 mHz. At dt=0.5
s over 300 minutes, the frequency grid remains the same; that check is a
sampling check, not a frequency-grid refinement. The nonzero taper support
[1.9, 20.1] mHz lies inside the computed support, but a separate widened-compute-
band comparison was not run, so band-support sensitivity remains untested.

For ray context, PREM TauP at the actual geocentric source/receiver coordinates,
depths and 99.7898974° distance predicts SKS at 1331.108 s, SKKS at 1364.750 s,
and nearly coincident S/ScS at 1393.963/1393.967 s; other predicted phases are
listed in the campaign audit evidence. The inspected bandlimited trace has a
substantial shear packet in 1260–1460 s, with peak 1.8188e-4 m. Since multiple
phases overlap, figures must call it a mixed shear packet (SKS/SKKS/S/ScS), not
an isolated S arrival. The ray predictions and visual inspection do not prove
phase isolation.

A bounded radial-kernel pilot with this deep Bolivia source, L300 and a
1260–1460 s transverse packet completed the forward seismogram but reached only
about 10% of 642 radial-kernel frequencies within its 150 s wall limit. No
kernel artifact was accepted. This does not validate or converge the kernel
examples; their science remains unresolved. No solver, API or output format was
changed. Numerical settings were not modified to conceal the duration result.


## Phase 2C rendered assets

The four existing tracked site image paths were regenerated from current-source saved outputs with the corrected plotting code. These are renderer/output reproducibility records; they do not certify numerical convergence. The first SVG and trace-strip share the Bolivia 300-minute baseline; the radial/section figures use the separate bounded `shear-packet` fixture.

| Asset | Source input and numerical artifact | Plotter | SHA-256 |
|---|---|---|---|
| `docs/site/assets/first-seismogram.svg` | `docs/site/examples/first-seismogram.json`; `/tmp/yspecpp-phase2a-bolivia-20261023/base60/bolivia.out.1` | `docs/site/examples/plot_first_seismogram.py` | `4b9dfc1c8d3f88d002b3f5444de6eb12f41e521532e389bed79aa1720da5e35b` |
| `docs/site/assets/traces-small.png` | Same Bolivia input and trace | `python/yspec/plot/seismograms.py` | `a311892cc0670b71ade5642d23e7b8fad1c21ba724f2445cbaf3738b339d0592` |
| `docs/site/assets/kernels-radial-small.png` | `docs/site/examples/kernels-small.json`; `/tmp/yspecpp-phase2c-runs-20261023/kernels-web/kernel-small.kernels` | `python/yspec/plot/kernels.py` | `768f7ff113e395505e78de7ebcde3b6c5488941b0a2ff8aa943021d3f89d745d` |
| `docs/site/assets/kernels-section-small.png` | Same kernel input; `/tmp/yspecpp-phase2c-runs-20261023/kernels-web/kernel-small.sections` | `python/yspec/plot/sections.py` (CLI: `banana --measurement shear-packet --parameter vs`) | `10d809540bfaa56a9911291adad682a713911e191dd85c739c74eb335269aceb` |

The SVG renderer read the input quantity and receiver sidecar and drew all 18,001 samples. For the trace-strip render, a temporary copy of the canonical JSON changed only the output prefix and source name to match the already-saved `bolivia.out.*` filenames; the physical settings match the canonical JSON. The `shear-packet` forward/adjoint run used two OpenMP threads and reported a 13.08-second runtime; its saved measurement, radial-kernel, section and map hashes are recorded in the Phase 2A evidence handoff.


The section asset was rendered with the documented `--measurement shear-packet --parameter vs` selectors. Its temporary input copy changes filesystem paths only, so it reads the saved kernel output while preserving the fixture physics; the copy SHA-256 is `571793caca81cf85821e331c01937cd33db2d2e1b4a25902347e4e7f0b50477b`. The four assets were rendered with Python 3.12.3, NumPy 2.5.4, Matplotlib 3.11.2, and PlanetModel 1.3.0. Cartopy 0.26.0, SciPy 1.18.1 and ObsPy 1.5.1 were also present in the isolated plotting environment; these four outputs do not use a Cartopy basemap. The plotting sources are the current working-tree files based on source revision `a294fb321cbcf9b6736710eb7ebdcec2c694dc5c`; their hashes are recorded here so later edits cannot be mistaken for the renderers used.

| Input/model/artifact | SHA-256 |
|---|---|
| `docs/site/examples/first-seismogram.json` | `ce09532d22437c8e3266ff855f83dd05e039328bd2fe116522aa8ccf87ac7ab3` |
| `examples/prem.200` | `4c3c53ca04e785ff4fd853cf211fe94e353ef70d6052a526122d0f279a3519ea` |
| Bolivia `bolivia.out.1` | `196c82b1f2ef926146fe5b32e219891374a4af632e75886e52353214f91d7021` |
| `docs/site/examples/plot_first_seismogram.py` | `663ed44a12176d9d20dfa3c6e7536dd7c4a3cbb0615f8b5f7609f2e16d9f14a7` |
| `python/yspec/plot/seismograms.py` | `aac9a45d7ffbc166d97dc42c1536870248bca817c88bcd110fecf98899fe0b5e` |
| `docs/site/examples/kernels-small.json` | `7b0efae842bff4b414b89eaecf209a986d6e918e0641e147cea0881d46009cc0` |
| `python/yspec/plot/kernels.py` | `ebf0b79a4a2787aef2d2450e5ba9e50b153643a8c70de1873d19a92a1801d6db` |
| `python/yspec/plot/_labels.py` | `39ecbede9900941fe68faa4b28effcd55d8f34f3c8ab9771c34de98f30d15e92` |
| `python/yspec/plot/sections.py` | `bb9e1cbfa1483377fc648cd03f48c0c22453da2d8fa0f6c315b81b98ebc5e762` |
| `kernel-small.kernels` | `b9929b9085373b6965c20c58a39831dd493862e2be1438fc6917fb782694dea5` |
| `kernel-small.sections` | `a3dfb4943dc599d45600a1d3190beac03336937e95b6a33641a027e5cca036dd` |

The first SVG and trace-strip were generated from `/tmp/yspecpp-phase2a-bolivia-20261023/base60/`; the renderer-only temporary JSON changed the output prefix and source name to match those saved filenames. The kernel outputs came from `/tmp/yspecpp-phase2c-runs-20261023/kernels-web/`. These absolute paths are provenance pointers for the local audit; the raw outputs are not distributed with the site.

The final section-asset correction selects `vs` explicitly for the transverse `shear-packet`. The renderer now defaults transverse measurements to `vs`; the prior render had selected `vp` because the name lacked an uppercase `S`. The Phase 2B caption correction left the corrected-parameter kernel values unchanged. The along-path and cross-path panels now label the 3000 km radius crop; deeper values are omitted from the display, not shown as zero. Limit, slice-maximum (the evaluated slice maximum before clipping) and saturated-fraction annotations now occupy three short lines in each spatial panel; a single-parameter radial render was checked to keep its four-line scale and source-term footer within the figure canvas. The radial caption says that each parameter panel has an independent horizontal SI scale and that small plotted terms do not establish significance.


## PR5 P2 trace-specific renderer refresh

The tutorial renderer was run from the saved 18,001-sample Bolivia trace and its receiver table using a temporary input copy that changes only the output prefix and source name to locate those saved files. The renderer now resolves the matching receiver row before choosing distance and units. The regenerated SVG is byte-identical to the tracked illustration (SHA-256 `4b9dfc1c8d3f88d002b3f5444de6eb12f41e521532e389bed79aa1720da5e35b`); no saved waveform was changed. The current renderer hash is `663ed44a12176d9d20dfa3c6e7536dd7c4a3cbb0615f8b5f7609f2e16d9f14a7`. The pre-existing 300-to-600-minute duration comparison still fails the 1% gate; the illustration remains an unresolved comparison setup, not a converged reference.
