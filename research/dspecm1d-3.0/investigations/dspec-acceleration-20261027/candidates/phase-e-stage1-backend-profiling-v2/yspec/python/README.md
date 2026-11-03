# yspec-seismo

Synthetic seismograms, snapshots of the wavefield and adjoint sensitivity
kernels for spherically symmetric earth models, from Python. It is the
Python package of [YSpec](https://github.com/da380/YSpecpp), which integrates
the equations of motion of a self-gravitating, anelastic, transversely
isotropic earth directly in radius, in the frequency domain, for every degree
and frequency, and sums the responses: complete seismograms, with every body
wave, surface wave and normal mode the band holds.

```sh
pip install yspec-seismo            # the library and the programs
pip install "yspec-seismo[plot]"    # and yspec-plot, to draw what they write
```

It is imported as `yspec`:

```python
import yspec

model = yspec.Model.prem()  # includes the ocean
source = yspec.Source(
    latitude=-13.909919806613164, longitude=-67.25, depth_km=647.1,
    moment_tensor=dict(rr=-7.59e20, rt=-2.5e21, rp=4.2e19,
                       tt=7.75e20, tp=-2.48e20, pp=-1.6e19),
    name="bolivia")
receivers = [yspec.Receiver(latitude=80.06588632614357, longitude=0,
                            depth_km=3)]
calculation = yspec.Calculation(
    degrees=(0, 300), frequencies_mhz=(1.5, 21), length_min=300,
    time_step_s=1, filter_mhz=(1.9, 2.0, 20.0, 20.1),
    gravity="full", attenuation=False, quantity="displacement")

seismograms = yspec.seismograms(model, calculation, source, receivers,
                                progress=True)
seismograms.time_s           # 18,001 samples, from 0 to 18,000 s
seismograms.data              # [source, receiver, component, sample]: Z, N, E
```

Lengths are in km, angles in degrees, frequencies in mHz and times in s, as
in the programs' input files; the results are numpy arrays in SI units.
The source and receiver latitudes are geographic coordinates converted from the geocentric locations in the published Bolivia comparison. Its 300-to-600-minute waveform difference exceeds the current 1% duration target, so this reproduces the comparison setup without claiming a converged reference.

For plots, `yspec-plot traces --components ZNE` selects the physical Z, N and E
components; the CLI default is ZRT, with R radial away from the source and T
90° anticlockwise from it. Components with compatible units at each receiver
share a physical amplitude scale. `--scale independent` labels a separate
scale for every component. The comparison plot checks exact
time samples, physical quantity, units and receiver geometry before
subtracting runs. Instrument traces retain the output units declared by their
response metadata; plots never share a physical scale across incompatible
units. R and T remain available as radial-away and 90° anticlockwise
transverse components.

- `yspec.seismograms` and `yspec.strains`: at receivers anywhere but the
  centre, in the ocean, the crust or the core, from moment tensors and forces,
  through seismometers given by their poles and zeros.
- `yspec.snapshots`: the wavefield at chosen times on a grid of depths,
  distances and azimuths around each source: cross-sections through the
  earth, or maps.
- `yspec.traveltimes`, `yspec.radial_kernels` and `yspec.sampled_kernels`:
  cross-correlation traveltimes and their sensitivity kernels by the adjoint
  method, to spherically symmetric perturbations of density, the elastic
  parameters and the quality factors, and to laterally varying ones of the
  elastic parameters and quality factors, through the earth, on sections
  about the path or on maps.
- The programs `yspec`, `yspec-snapshot` and `yspec-kernels`, which read
  JSON input files and write their results as text, and `yspec-plot`, which
  draws them (with the `plot` extra).

The calculations run on as many threads as OpenMP gives them
(`OMP_NUM_THREADS`, or `yspec.set_threads`), except on macOS, where the
wheels have one. Ctrl-C stops a calculation.

## References

- Al-Attar, D. & Woodhouse, J. H., 2008. Calculation of seismic displacement
  fields in self-gravitating earth models -- applications of minors vectors
  and symplectic structure. Geophys. J. Int., 175, 1176-1208.
- Friederich, W. & Dalkolmo, J., 1995. Complete synthetic seismograms for a
  spherically symmetric earth by a numerical computation of the Green's
  function in the frequency domain. Geophys. J. Int., 122, 537-550.

BSD 3-Clause licence.
