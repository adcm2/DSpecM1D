# YSpec

YSpec computes synthetic seismograms for a spherically symmetric earth model
by direct radial integration of the equations of motion in the frequency
domain (Friederich & Dalkolmo 1995), including self-gravitation (Al-Attar &
Woodhouse 2008). This is a C++23 re-implementation of the original Fortran
YSpec: the same physics, inputs and outputs, with better numerics and a
structure that is easier to extend.

Start with the [YSpecpp user guide](docs/site/index.html), including a [first seismogram](docs/site/first-seismogram.html) and a [C++ library example](docs/site/library-example.html).

The earth model may be transversely isotropic and anelastic, with any
sequence of solid and fluid regions (a fluid or solid inner core, an ocean or
none). Attenuation is the constant-Q logarithmic dispersion of PREM,
implemented exactly at complex frequency rather than to first order in 1/Q.
Gravity can be left out, included in the Cowling approximation, or included
in full.

It is a C++ library with three programs, and a Python package,
`yspec-seismo` (see the Python section below).

## Building

YSpec needs a C++23 compiler (GCC 13 or Clang 18 and later), CMake 3.24 or
later, FFTW3 (in double, single and long double precision, as FFTWpp
requires; or none, with pocketfft, below), netCDF-C for writing
three-dimensional kernels (`-DYSPEC_WITH_NETCDF=OFF` builds without it), and
optionally OpenMP. Everything else is fetched by CMake:
[NumericConcepts](https://github.com/da380/NumericConcepts),
[Interpolation](https://github.com/da380/Interpolation),
[FFTWpp](https://github.com/da380/FFTWpp),
[GaussQuad](https://github.com/da380/GaussQuad),
[GSHTrans](https://github.com/da380/GSHTrans),
[planetmodel](https://github.com/da380/planetmodel) and, for the tests,
GoogleTest.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build
```

To build against local checkouts of the dependencies instead of fetching
them, pass `-DFETCHCONTENT_SOURCE_DIR_<NAME>=/path/to/checkout` for each
(`GSHTRANS`, `INTERPOLATION`, ...). `-DYSPEC_WITH_OPENMP=OFF` builds a serial
version, which gives the same results. `-DYSPEC_FFT=pocketfft` (the preset
`pocketfft`) makes YSpec's one Fourier transform, of each spectrum to its
time series, with [pocketfft](https://github.com/mreineck/pocketfft) instead
of FFTW, and needs no FFTW at all; the results agree to rounding. The Python
wheels are built so: FFTW is GPL, and a wheel carries the libraries it uses.

## Running

There is a program for each kind of output: `yspec` computes seismograms at
receivers, `yspec-snapshot` the wavefield at chosen times on a grid around
each source, and `yspec-kernels` the sensitivity of traveltimes to the
model. They read JSON files that share most of their members.

```sh
cd build
bin/yspec examples/input.json
bin/yspec-snapshot examples/snapshot.json
bin/yspec-kernels examples/kernels/traveltimes.json
```

The build copies the example inputs and models into `build/examples`, so the
outputs are written there, beside their inputs, and not in the sources.

The number of threads is set with `OMP_NUM_THREADS`. For `yspec`, one file is
written per receiver, `<prefix>.1`, `<prefix>.2`, ..., in the order the
receivers are listed; with several sources, one per source and receiver,
`<prefix>.<source>.1`, .... Each line holds the time after the origin time (s)
and the vertical (up), north and east components, in SI units of the quantity
asked for, or in the units of the instrument's output where there is one. A
table `<prefix>.receivers` lists the files, one line for each source and
receiver: the receiver's latitude, longitude and depth (km), and its distance,
its azimuth at the source and the source's azimuth at it (degrees).

### Examples and plotting

Each directory of `examples/` holds an example, with its input and a script
that runs it and plots the result:

| example | program | what it shows |
|---------|---------|---------------|
| `record-section` | `yspec` | displacement every 5 degrees along a great circle, vertical and transverse |
| `instrument` | `yspec` | one place recorded twice: ground displacement and a seismometer response in its declared output units |
| `depths` | `yspec` | one place at depths from the sea floor to the inner core |
| `gravity` | `yspec` | a long-period seismogram without gravity, in the Cowling approximation and in full |
| `force` | `yspec` | the waves from a point force, of the kind a landslide exerts |
| `cross-section` | `yspec-snapshot` | the wavefield on a great circle through the whole earth at six times |
| `map` | `yspec-snapshot` | the wavefield at the sea floor around the source, over the whole globe |
| `kernels` | `yspec-kernels` | radial kernels for an S-window measurement and a long-period vertical wave packet; neither is claimed to isolate one phase |
| `kernels-3d` | `yspec-kernels` | a three-dimensional kernel for a mixed transverse shear packet in an S-arrival window; its banana and Fresnel-zone sections |
| `kernels-gallery` | `yspec-kernels` | three-dimensional kernels for cross-correlation windows targeting P, PP, S, ScS, SS, Pdiff and Sdiff packets, on sections along and across the path |

```sh
cd build
examples/record-section/run.sh     # or examples/run-all.sh for all of them
```

Each script writes its solver output and figures in its own directory. Runtime
depends on the input and build; kernel galleries are more expensive because
they compute several adjoint fields. The programs are found in the build, on
the `PATH`, or in `YSPEC_BIN`.

The figures are drawn by `yspec-plot`, which the Python package installs
with its `plot` extra (see the Python section below), built on
[planetmodel](https://github.com/da380/planetmodel), whose sections and
colour scales it shares, with maps through cartopy. The scripts find it on
the `PATH`:

```sh
python -m venv .venv && . .venv/bin/activate
pip install -e ".[plot]"
yspec-plot section build/examples/record-section/input.json
```

The package installs the programs too, so the examples also run with it
alone, without the C++ build.

`yspec-plot` takes the input file of a run and finds the output from it:
`section` draws a record section, `traces` the components at a few receivers,
`compare` one receiver in several runs, `snapshot` a section through the
earth or a map around the source, `kernels` the kernels of a measurement
against radius, below its seismogram and window, `volume` its
three-dimensional kernels, from the model file, on a section along the great
circle through the source and the receiver and on maps at some depths,
`banana` its kernel on the sections along the path and across it, below its
seismogram and window, and `gallery` those of every measurement of a run, a
row for each. With several measurements, `--measurement` chooses one by
name. Trace panels share an amplitude scale across compatible components
within each receiver row by default; `traces --scale independent` labels the
deliberately independent alternative. `compare` checks that runs have matching
physical quantities, units, receiver geometry and time samples before it
subtracts them. Instrument output units come from the response metadata when
present; they are not assumed to be counts. The components are Z, N and E, or
R and T, radial away from the source and transverse 90° anticlockwise from it.
Without `-o` the figure is shown.

Snapshot and spatial-kernel colour bars report the displayed units, limits,
per-slice maxima before clipping and saturated fractions. Snapshot panels use independent limits
per time unless `--common` is requested. Banana and gallery panels use their
own symmetric 99th-percentile limits: values above 150 km depth are omitted
from the limit calculation, and along-path limits also omit points within 5°
of the source and receiver. Volume slices use their own percentile limits.
Non-K banana and gallery sections are cropped at a radius of 3000 km; their
titles show the displayed radial extent, and omitted deeper values are not
shown as zero.
These are display scales only: saturation changes colours, not saved numerical
values. Gouraud shading between saved section samples is visual interpolation,
not additional spatial resolution; snapshot maps assign saved values to
nearest-sample cells. Record sections intentionally normalize each trace and
may clip amplified traces to prevent overlap, as disclosed in the plot title.

The example figures demonstrate workflows and output interpretation; running
an example or passing a plotting check does not certify its waveform or kernel
as converged. The <a href="docs/site/examples/all-plots-provenance.md">plot status and provenance record</a>
tracks the source inputs, saved numerical evidence and remaining scientific
limits for each documented figure.

### Input

The calculation is read from a JSON file:

This root input uses the published DSpecM1D Bolivia comparison settings: PREM with ocean, no attenuation, full gravity, displacement, degrees 0–300, 1.5–21 mHz, a 300-minute record at 1 s, and the stated four-corner taper. Its source and receiver latitudes are geographic conversions of the comparison’s geocentric coordinates. Its duration check remains unresolved; the setup is not a converged reference.

```json
{
  "model": "prem.200",
  "attenuation": false,
  "gravity": "full",
  "output": {
    "prefix": "yspec.out",
    "quantity": "displacement",
    "potential_and_tilt_corrections": false
  },
  "degrees": {
    "min": 0,
    "max": 300
  },
  "frequencies_mHz": {
    "min": 1.5,
    "max": 21.0
  },
  "time_series": {
    "length_min": 300.0,
    "time_step_s": 1.0
  },
  "filter_mHz": [
    1.9,
    2.0,
    20.0,
    20.1
  ],
  "sources": [
    {
      "depth_km": 647.1,
      "latitude": -13.909919806613164,
      "longitude": -67.25,
      "moment_tensor_Nm": {
        "rr": -7.59e+20,
        "rt": -2.5e+21,
        "rp": 4.2e+19,
        "tt": 7.75e+20,
        "tp": -2.48e+20,
        "pp": -1.6e+19
      }
    }
  ],
  "receivers": {
    "depth_km": 3.0,
    "locations": [
      [
        80.06588632614357,
        0.0
      ]
    ]
  },
  "solver": {
    "tolerance": 1e-05
  }
}
```

- `model`: an earth model in the deck format of Minos (six or nine columns).
  Relative paths, here and in `output.prefix`, are relative to the input
  file.
- `attenuation`: whether to include attenuation (default `true`).
- `gravity`: `"none"`, `"cowling"` or `"full"` (default `"full"`).
- `output.quantity`: `"displacement"`, `"velocity"` or `"acceleration"`
  (default displacement). `potential_and_tilt_corrections` applies to
  acceleration: it adds the change in gravity seen by an instrument that moves
  with the ground, and the tilt of its axis (default `false`).
- `degrees`: the spherical harmonic degrees summed (`min` defaults to 0).
  What the last degree cuts off is cut off sharply, and the tail of the cut
  runs through the whole record, so enough are needed for two things, the
  larger setting `max`:
  - the fundamental surface waves at the highest frequency, l ≈ 2πfa/c,
    about 10 for each mHz, even where only body waves are wanted;
  - the near field of a shallow source, whose terms fall off as
    e^(−l Δd/a), Δd the difference in depth between the source and the
    receiver: n a/Δd, with n about 6 in bands that reach down to 1 mHz, 8
    at the longest periods, and 2 at higher frequencies if the first few
    minutes do not matter (a/Δd ≈ 240 for the examples' source, 29.5 km
    deep, and receivers at the sea floor). What is cut off of it is a pulse
    at the origin time, at any distance, which a narrow low ramp of the
    filter spreads through the record. Doubling `max` and seeing whether the
    start of the record moves is the check.

  Receivers belong at the sea floor or below. In the water, with gravity,
  the ocean's surface-gravity waves (tsunamis), about 250 degrees for each
  mHz and far more above 5 mHz, would be needed as well; cut off, they swamp
  the horizontal motion there.
- `frequencies_mHz`: the band computed; `time_series` the length of the
  seismograms and their sampling. The Nyquist frequency must lie above the
  highest frequency; without `time_step_s` the step is a quarter of the
  shortest period. The damping ε is 5 e-folds over the length, and the
  frequencies are spaced so that the transform's period, a power of two of
  samples, is at least 2π/ε: what arrives after the period wraps round,
  damped by e^(−εT). The time step sets where in that range the period falls,
  and with it the number of frequencies: a longer one can hold late surface
  waves that would otherwise wrap round onto an early phase.
- `filter_mHz`: an optional raised-cosine taper, rising from zero to one
  between its first two corners and falling back between its last two. The
  taper has zero phase, so its response begins before the origin time, by
  about the inverse of the width of its lowest ramp; with the damping, that
  part wraps round and reappears, amplified, at the end of the time series.
  A low ramp of 1 mHz or more, or a longer time series, keeps it out of the
  record.
- `sources`: one or more, each with its hypocentre (geographic latitude) and
  its moment tensor `moment_tensor_Nm`, a point force `force_N` (members
  `r`, `t`, `p`, each defaulting to zero), or both, in spherical coordinates
  at the source (r up, theta south, phi east), as a step at the origin time.
  Sources at one depth share the radial solves, so extra ones cost little. The optional `name`
  labels the output files when there are several sources (by default the
  sources are numbered from 1).
- `receivers`: a group of receivers, or an array of groups (see below). They
  may lie anywhere but at the centre: in the solid region of the sources, in
  an ocean or in the core. The sources must lie in one solid region of the
  model: in `prem.200`, below its 3 km ocean and above the core.
- `solver.tolerance`: the local error allowed per integration step, relative
  to the size of the solution (default `1e-5`; see below), and
  `solver.start_decay`, the decay in e-folds below a turning point at which an
  integration starts (default 11).

A member not listed here is an error, so that a misspelt one is not ignored.

Each group of receivers has its depth, or several, and its places: either a
list of `locations`, or a `profile` along a great circle, from one point `to`
another or from a point at an azimuth (clockwise from north) between two
distances, with `count` points equally spaced, both ends included:

```json
"receivers": [
  { "depth_km": 3, "locations": [[34.95, -106.46]], "instrument": "ANMO.pz" },
  { "depth_km": [3, 100, 200],
    "profile": { "from": [36.87, 35.58], "azimuth_deg": 60,
                 "distances_deg": [10, 90], "count": 41 } }
]
```

Receivers are numbered group by group and, within a group, depth by depth,
with every place at each depth. Profiles are great circles of the spherical
model, so a profile that starts at a source has points at exactly the
epicentral distances asked for. The optional `instrument` is a SAC pole-zero
file, applied to every component at the group's receivers; its output (counts,
for the files of the IRIS DMC) is then recorded. A comment `* INPUT UNIT : M`
(or `M/S`, or `M/S**2`) in the file says what motion the instrument responds
to, which is otherwise displacement, as in SAC; the conversion from the
quantity computed is exact. The response is applied at the same complex
frequencies as the rest of the calculation, which gives the exact causal
response.

### Snapshots

`yspec-snapshot` reads the same members as `yspec`, except that its `output`
has only `prefix` and `quantity`, and that `snapshots` replaces `receivers`:

```json
"snapshots": {
  "times_s": [300, 600, 900, 1200],
  "depths_km": { "min": 25, "max": 2875, "count": 58 },
  "distances_deg": { "min": 0, "max": 180, "count": 181 },
  "azimuths_deg": [30, 210]
}
```

The grid is in the frame of each source: depths, distances from it and
azimuths at it (clockwise from north). Each of the four is one number, an
array of them, or a range of `count` values from `min` to `max`, both
included. Two azimuths 180 degrees apart give a cross-section through the
source along a great circle; one depth and many azimuths, a map around it.
The points may lie anywhere but at the centre, through the oceans and the
core as well as the solid region of the sources.

One file is written for each time, `<prefix>.1`, `<prefix>.2`, ... (or
`<prefix>.<source>.1`, ... with several sources). After two comment lines,
each line holds a point: its depth (km), distance and azimuth (degrees), its
latitude and longitude, and the up, away and transverse components, away
being along the great circle from the source and transverse 90 degrees
anticlockwise from it, seen from above.

The cost grows with the numbers of depths and distances, while further
azimuths cost almost nothing: around the source only orders |m| <= 2 are
excited, so each depth and distance needs five sums over degree, which any
number of azimuths then share. The values at each time are summed directly
from the spectra, and agree with the samples of seismograms at the same
points to rounding error. The example above, at degrees to 300 and up to 15
mHz, takes about 3 s on 4 threads, and the same section through the whole
earth, from the sea surface to near the centre, about 4 s.

With gravity, points in an ocean need far more degrees than points in the
solid, for the ocean's surface-gravity waves: see `degrees` above.

### Kernels

`yspec-kernels` reads the same members as `yspec`, with one source, except
that its `output` has only `prefix` and `quantity`, and that `receiver`,
`measurements` and `kernels` replace `receivers`:

```json
"receiver": { "depth_km": 3, "distance_deg": 60, "azimuth_deg": 30 },
"measurements": [
  { "name": "S", "component": "T",
    "traveltime": { "window_s": [1000, 1200], "taper_s": 40 } },
  { "name": "long-period-Z-packet", "component": "Z",
    "traveltime": { "window_s": [1550, 2050], "taper_s": 100 } }
],
"kernels": {
  "parameters": "velocities",
  "scaling": "relative",
  "mesh": { "order": 5, "max_element_km": 50 },
  "volume": { "degree": 48, "max_element_km": 100 },
  "sections": { "margin_deg": 15, "half_width_deg": 15, "spacing_deg": 0.25,
                "radial_spacing_km": 30 },
  "maps": { "depths_km": [200, 600, 1500, 2700], "margin_deg": 10,
            "half_width_deg": 15, "spacing_deg": 0.5 }
}
```

The receiver is at a `location`, or at a distance and azimuth from the
source, and may have an `instrument`, as in `yspec`. Each measurement is the
traveltime of the seismogram of its `component` (`Z`, `N`, `E`, or `R` and
`T`, radial and transverse to the source), measured by cross-correlation with
the synthetic itself in its window, with cosine tapers of the given length
inside it; a window of a period or two about one phase keeps its kernel
simple. Each has its own kernels. They are named by `name`, or numbered, and
share the receiver, so that they share the forward field, and the adjoint
fields of at most three directions, of which each measurement's is a
combination. Their sensitivity to spherically symmetric perturbations of the
model is computed by the adjoint method (see
`docs/theory/adjoint-kernels.tex`): the adjoint field is the response, at the
same complex frequencies, to a force at the receiver, so the cost is about
that of the seismograms of a source at a receiver at every node of the mesh
for the forward field and for each direction.

`parameters` chooses the kernels: `love` (density, A, C, F, L, N and their
inverse quality factors, the model's own), `moduli` (density, the bulk and
shear moduli, and their inverse quality factors), `velocities` (the default:
density at fixed velocities, the isotropic speeds and the inverse quality
factors) or `anisotropic` (density at fixed velocities, the five speeds and
eta of transverse isotropy, and the inverse quality factors). `scaling` is
`absolute` (the default), for perturbations of the parameters, or `relative`,
for perturbations of their logarithms. The mesh has Gauss-Legendre nodes in
elements between the discontinuities of the model, the source and the
receiver, no larger than asked.

For each measurement, files `<prefix>.<name>.*` are written, or
`<prefix>.*` if there is only one, as they are called below.
`<prefix>.kernels` has a header that describes the
measurement and the units and gives the kernels' point parts at the source
(on the line beginning `# at_source`), and a line for each node: its radius
and depth (km), its weight in the integral over radius (m^3), and each
kernel, so that the change in the traveltime is the sum over the nodes of
weight times kernel times perturbation, plus the point parts times the
perturbations at the source. The point parts come from the jump in the
displacement across the radius of a moment tensor, which makes the kernels
of the moduli that set it concentrated there. `<prefix>.measurement` has a
line for each sample: the time, the seismogram measured, the window and the
adjoint source. Both kernels and point parts predict finite differences of
the measurement to about 1e-6, in every treatment of gravity; the tests
check them. The current two-window radial run (L ≤ 300, 0.2–15 mHz, 40 minutes at 2 s)
completed in 24.05 s on two threads; that timing is environment-specific. The
long-period Z window is named `long-period-Z-packet`: no surface-wave or
normal-mode comparison has validated it as Rayleigh, and convergence remains
unresolved.

With `volume`, the three-dimensional kernels are computed as well, for the
same parameters at fixed density: laterally varying perturbations of density
would take the model out of hydrostatic equilibrium, and are left out (see
the derivation). They are written as a model file of
[planetmodel](https://github.com/da380/planetmodel), `<prefix>.nc`, on the
skeleton of the model's discontinuities, with nodes in each layer no further
apart than `max_element_km`, both sides of each discontinuity included, and at
each on a Gauss-Legendre grid resolving degree `degree` and at the poles; a
field `K_vs`, say, for each parameter, read linearly between the nodes.
planetmodel draws sections and maps of them, and writes them for ParaView.
Each field is evaluated in its own frame, the forward field with the source at
the pole and the adjoint field with the receiver there, where they need only
the few orders a point source excites, and the two strains are turned into one
frame at each point; the integral of each kernel over every sphere is the
radial kernel, and the kernels at single points agree with those computed in
the source's frame alone, with the adjoint field's every order, which the
tests check. The cost is about that of the seismograms of a source at each
point of the grid at each node, for the forward field and for each
direction, so modest grids are wise: the example, to 10 mHz and degree 200 on
a grid of degree 48, takes about 45 s on four threads. Writing model files
needs netCDF-C (`libnetcdf-dev`); YSpec builds without it, with
`-DYSPEC_WITH_NETCDF=OFF`, and does everything else.

With `sections`, with or without `volume`, the three-dimensional kernels are
computed on two sections: along the great circle through the source and the
receiver, from `margin_deg` before the source to as far beyond the receiver,
and across it at its middle, `half_width_deg` to each side, at points
`spacing_deg` apart, at nodes in each layer no further apart than
`radial_spacing_km`, both sides of each discontinuity included. They are
written to `<prefix>.sections`, a line for each point: the section (`path`
or `across`), the layer, the depth (km), the distances along the path and
across it (degrees), the latitude and longitude, and each kernel. With
`maps`, they are computed on maps at the depths given, over the same stretch
of the path and `half_width_deg` to each side of it, and written to
`<prefix>.maps` in the same way, without the section; a map at a
discontinuity is on its upper side. Sections and maps need no netCDF. The
gallery's run at 60 degrees, five measurements on two components to 40 mHz
and degree 450 on its two sections, takes about two minutes on four
threads.

## Accuracy and speed

The integration controls its step size, so its accuracy is set by
`solver.tolerance`. On the historical benchmark input (degrees to 600, up to 40 mHz, 60
receivers), against seismograms from the Fortran with sixteen times its usual
step density, the largest error relative to the peak of each component is:

| tolerance       | error   | time, 16 threads |
|-----------------|---------|------------------|
| 1e-4            | 1.1e-3  | 4.9 s            |
| 1e-5 (default)  | 2.2e-4  | 8.0 s            |
| 1e-6            | 2.9e-5  | 12.3 s           |
| Fortran YSpec   | 2.3e-3  | 5.9 s            |

(on an 8-core laptop).

## How it works

For each frequency (complex, ω − iε, so that the response is damped and its
time series does not wrap around) and each degree, the radial
equations are solved for unit jumps in displacement and traction at the source
radius. The source is rotated to the north pole, where only orders |m| ≤ 2
are excited; the moment tensor and the receivers enter only
afterwards, through the weights of the jumps and the angular functions.

- **Model** (`YSpec/Model.hpp`, `src/Model.cpp`). The deck is read with
  planetmodel and made non-dimensional. Each property is a natural cubic
  spline within each layer, and gravity is the exact integral of the density
  spline.
- **Equations** (`src/Equations.hpp`). Toroidal and radial systems are
  2x2. Spheroidal systems are 4x4 without self-gravitation and 6x6 with it,
  integrated as minor vectors (second-order minors of the 4x4 systems, the
  fourteen independent third-order minors of the 6x6), with the source carried
  to the receiver by the b-vector system. Solid and fluid regions have their
  own forms, joined at every boundary between them.
- **Solver** (`YSpec/RadialSolver.hpp`, `src/RadialSolver.cpp`). One
  two-point solver serves the toroidal, radial and spheroidal systems, all
  three treatments of gravity and any sequence of regions. The integrations
  from below and from above are shared by every source and receiver depth:
  each is kept as it passes them, and each source then needs only its own
  determinant and the b-vectors carried out to the receivers, which change
  form with the minors at each boundary between solid and fluid, so that a
  receiver may be in an ocean or the core. Integration starts where a JWKB
  estimate says the solution has decayed by e^-11 below its turning point,
  beneath the deepest source, from starting values for a homogeneous sphere;
  below that start the response is negligible and is taken to be zero.
- **Integration** (`src/RungeKutta.hpp`). A Dormand-Prince 5(4) method
  with step-size control, restarted at every discontinuity. Solutions are
  carried with a running logarithmic scale, since they grow and decay by far
  more than a double can hold.
- **Receivers** (`src/Geometry.hpp`, `src/AngularFunctions.cpp`). The angular
  functions come from the generalised Legendre functions of GSHTrans, which
  are regular at both poles.
- **Time series** (`src/TimeSeries.cpp`, `src/Synthetics.cpp`). Spectra
  are rotated to north and east, passed through any instrument, filtered,
  transformed with FFTW and undamped.
- **Snapshots** (`src/Snapshots.cpp`). The coefficients of each degree and
  order are gathered at each frequency, summed over degree at each depth and
  distance a block of distances at a time, and summed over frequency directly
  at each time.
- **Kernels** (`YSpec/Kernels.hpp`, `src/Kernels.cpp`). The forward and
  adjoint fields come from one solver, with the source and the receiver as its
  two source radii and the nodes of the radial mesh as its receivers. The
  adjoint field is expanded in the source's frame, where only orders |m| ≤ 2
  meet the forward field, so the integrals over each sphere are sums over
  degree and order of products of radial functions. The density's effect on
  gravity is a radial integral on the same nodes.

The calculation is threaded over frequency, with a loop over degree inside
each thread. The library's public headers, in `YSpec/`, use only the standard
library; the numerics and the libraries they build on are in `src/`.

## Using it as a library

Link the CMake target `YSpec::YSpec`. The public headers in `YSpec/` use only
the standard library, and the library works in SI units.

```cpp
#include <YSpec/YSpec.hpp>

const auto model = YSpec::Model::ReadDeck("prem.200");

auto calculation = YSpec::Calculation{};
calculation.degrees = {.min = 0, .max = 200};
calculation.band = {.min = 0.2e-3, .max = 15e-3};           // Hz
calculation.sampling = {.duration = 3600, .timeStep = 20};  // s

const auto sources = std::vector<YSpec::Source>{
    {.location = {0, 0}, .depth = 29.5e3, .moment = {.rr = 1e18}}};
const auto receivers = std::vector<YSpec::Receiver>{
    {.location = {0, 30}, .depth = 3e3}, {.location = {0, 60}, .depth = 3e3}};

const auto seismograms =
    YSpec::ComputeSeismograms(model, calculation, sources, receivers);
// The vertical component from the first source at the first receiver.
const auto z = seismograms.Trace(0, 0, YSpec::Vertical);
```

Sources may be at any depths in one solid region of the model (below an
ocean, for instance, and above the core), and receivers anywhere but at the
centre. The radial integrations are shared by all of them: several sources,
or receivers at several depths, cost little more than one. A source may be a
moment tensor, a point force (`YSpec::Force`) or both, and
`YSpec::ComputeStrains` gives seismograms of the strain tensor instead of the
displacement, as strainmeters and fibre-optic cables record; strain and
displacement satisfy reciprocity with each other to 1e-8, which the tests
check. A receiver may carry an `Instrument`, read from a
SAC pole-zero file with `YSpec::ReadSacPoleZero`, and `YSpec::GreatCircle`
gives lines of receivers along great circles.

`YSpec::ComputeRadialKernels` gives the sensitivity of a measurement on a
seismogram to spherically symmetric perturbations of the model, by the adjoint
method: of the density (with the gravity that follows from it), of the five
Love moduli and of their quality factors. A measurement is given by its
adjoint source, and `YSpec::Traveltime` forms that of a cross-correlation
traveltime in a window, on the seismogram along any direction
(`YSpec::Along` a component, `YSpec::RadialDirection`,
`YSpec::TransverseDirection`), and `YSpec::Traveltimes` those of several from
one synthetic; given several measurements at one receiver,
`YSpec::ComputeRadialKernels` gives the kernels of each, sharing the fields.
`YSpec::Reparameterise` gives the kernels in other parameters, absolute or
relative. The adjoint field is the forward response, at the same complex
frequencies, to a force at the receiver, so the kernels cost about what two
sets of seismograms at every node of a radial mesh would; they predict
finite differences of the measurement to about 1e-6, which the tests check.
The derivation is in `docs/theory/adjoint-kernels.tex`.

`YSpec::ComputeSampledKernels` gives the three-dimensional kernels of
measurements at one receiver, at fixed density, on samplings of the model: a
`YSpec::VolumeSampling`, nodes in radius on a Gauss-Legendre grid of each
sphere, `YSpec::PathSection` and `YSpec::CrossSection`, sections along the
path and across it, and `YSpec::PathMaps`, maps about the path at chosen
depths.

`YSpec::ComputeSnapshots` gives the wavefield on a `YSpec::SnapshotGrid` at
chosen times. `YSpec::RunSeismograms("input.json")`,
`YSpec::RunSnapshots("snapshot.json")` and `YSpec::RunKernels("s-wave.json")`
do what the programs do: read an input file, compute what it describes and
write it. `YSpec::RadialSolver` gives the
responses to unit jumps at a set of source radii, at a set of receiver radii,
for one degree and frequency.

## Python

The Python package, `yspec-seismo`, is imported as `yspec`. It gives the
library from Python, with numpy arrays for its results, and installs the
programs, the PREM decks of the examples and, with the `plot` extra,
`yspec-plot`:

```python
import yspec

model = yspec.Model.prem()
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

This Python example matches the Bolivia input above. Its 300-to-600-minute
waveform check exceeds the 1% duration target, so it does not claim a
converged reference; see the <a href="docs/site/examples/first-seismogram-provenance.md">recorded comparison</a>.

It has seismograms and strains, snapshots, traveltimes and their radial and
three-dimensional kernels, in the units of the input files (km, degrees,
mHz, s); `python/README.md` is its page on PyPI. It is built by
[nanobind](https://github.com/wjakob/nanobind) and
[scikit-build-core](https://github.com/scikit-build/scikit-build-core) from
`pyproject.toml`, with the module in `python/src/bindings.cpp`, the Python
around it in `python/yspec` and its tests in `python/tests`:

```sh
pip install -e ".[plot,obspy,test]"   # from the sources, rebuilt on reinstall
pytest
```

The calculations release Python's lock and run on the library's threads;
Ctrl-C, or an exception from a progress callback, stops one at the next
frequency to finish. The wheels, for Linux (x86_64 and aarch64) and macOS
(arm64, on one thread, as Apple's compiler has no OpenMP), are built by
`.github/workflows/wheels.yml` with cibuildwheel and published to PyPI when
a tag `v*` is pushed, as `docs/publishing.md` describes. They are built with pocketfft and without netCDF, so
`volume` kernels are not written as model files from them; sections and
maps are.

## Differences from the Fortran YSpec

- Input is JSON instead of a parameter file read by counting lines.
- The integration controls its error, and is about ten times more accurate at
  the default tolerance.
- Several faults of the original are not reproduced: with gravity switched
  off the degree-zero equations no longer keep a gravity term; an ocean is no
  longer ignored at degree zero when the receiver is below the source; the
  horizontal components at the epicentre no longer have the wrong sign; and
  gravity near the centre is exact.
- Any sequence of solid and fluid regions is handled, including a solid
  region that reaches the centre, and a source or receiver on a solid-fluid
  boundary is placed on its solid side.
- The free-air and tilt corrections use gravity at the receiver rather than at
  the model knot below it, and with full self-gravitation the potential
  correction includes the gravity-flux term that the original omits; it
  vanishes at a free surface but not at a buried receiver (about 1e-3 of the
  peak acceleration on the sea floor).
- Receivers, and the points of snapshots, may be anywhere: in an ocean, in the
  outer core or in the inner core, as well as in the solid region of the
  sources. Their responses agree with direct shooting to 1e-8.
- Without self-gravitation and in the Cowling approximation, the starting
  solution of the S type lacked a factor of zeta on the Bessel-function term
  of its tangential displacement. Its error decays upwards from the start, so
  it only mattered for receivers near the start, which the original could not
  have.
- Degrees above 1827 are allowed: the angular functions are then tabulated in
  long double.

## Tests

`ctest --test-dir build` runs unit tests; comparisons of the solvers with the
Fortran's for single degrees and frequencies over a range of models and
options; checks independent of the Fortran (the normal modes of PREM, the
toroidal modes of a homogeneous sphere against their closed form, responses
in the oceans and the core against direct shooting of the solutions,
continuity at the epicentre, snapshots against seismograms, strain against
displacement by differences and reciprocity, kernels against finite
differences of their measurement in every treatment of gravity); and four
end-to-end cases compared with the Fortran's seismograms.
The scripts in `tests/fortran` regenerate the reference data from a checkout
of the Fortran YSpec. `pytest` runs the Python package's tests, which check
it against the programs' input files, through the same library.

## References

- Friederich, W. & Dalkolmo, J., 1995. Complete synthetic seismograms for a
  spherically symmetric earth by a numerical computation of the Green's
  function in the frequency domain. Geophys. J. Int., 122, 537-550.
- Al-Attar, D. & Woodhouse, J. H., 2008. Calculation of seismic displacement
  fields in self-gravitating earth models -- applications of minors vectors
  and symplectic structure. Geophys. J. Int., 175, 1176-1208.

## Licence

BSD 3-Clause; see `LICENSE`.
