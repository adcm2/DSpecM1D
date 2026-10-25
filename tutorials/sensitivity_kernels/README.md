# Windowed material sensitivity kernels

This small installed-package consumer evaluates the real scalar measurement
\(J=\Re\sum_{r,k}\overline{q_{rk}}y_{rk}\), then compares its discrete kernel
prediction with centered reruns for one solid-density bump and one fluid-bulk
bump. Here \(y\) is the SI displacement spectrum after the configured raw
bins are reconstructed to time, scaled to metres, multiplied by one Hann
window, and transformed back. Density is perturbed independently of the
absolute moduli; geometry, interfaces, source, receiver, frequency grid and
full-gravity formulation stay fixed.

The example uses degrees 0–2, four GLL nodes per element, a maximum element
width of 0.25 times the model radius, and raw FFT bins `[2,4)`. Its nonzero complex weights are on the
vertical component at bins 2 and 3. The source is at 33° N, 10° E and 20 km
depth with tensor order `(Mrr,Mrt,Mrp,Mtt,Mtp,Mpp)` and components
`(1,2,3,4,5,6)×10²⁰ N m`; the receiver is at 30° N, 65° E and 3 km depth.
A single Hann window is
`0.5*(1-cos(2*pi*n/(nt-1)))`. Processing first applies the inverse real FFT
and frequency spacing, compensates the exponential damping once, scales by
the SI radius, applies the Hann window, then uses the forward transform with
`dt` in seconds. The objective is a q-contracted processed displacement
spectrum (m s under the library's transform convention). The report records
the model's actual GLL element count, SI scales, the SI gravity unit constant,
the dimensionless gravity coefficient, and physical gravitational constant. It also records
`epsilon = FreqFull::ep()/TimeNorm` in s⁻¹, `z=2*pi*f-i*epsilon`, and the
Hann support `[0,(nt-1)*dtSI]`. The `f11..f22` corner values are configuration
metadata only: this example applies no frequency taper or filter and solves
only bins `[2,4)`.

Kernel units are objective/(Pa m) for absolute solid moduli and fluid bulk
modulus, and objective/((kg m⁻³) m) for density. Integrate each kernel against
its physical perturbation and the returned GLL `drWeightsMetres`. The density
representation is independent one-sided element-local GLL values with a
polynomial in each SEM element. The directions use `(1-x²)` bumps over one
selected element: +100 kg m⁻³ in a mantle element near 0.7R and a fluid bulk
amplitude study of +10 and +100 GPa in the element near 0.3R. Centered differences use
`[J(+a)-J(-a)]/(2a)` for `a = 0.02, 0.01, 0.005`. The report compares the
kernel result with Richardson `(4*FD(0.005)-FD(0.01))/3`, using relative error
`abs(Richardson-prediction)/max(abs(Richardson),abs(prediction))` (zero if both
are zero). The diagnostic gate is 1e-6. Both fluid amplitudes remain visible;
the 100 GPa direction is selected for the gate to improve the finite-difference
signal against the objective baseline. A failed selected gate remains visible
and makes the runner exit nonzero.

Build against an installed DSpecM1D package, then run the one report entry
point (Python 3 with NumPy and Matplotlib):

```sh
cmake -S tutorials/sensitivity_kernels -B /tmp/dspec-k2-consumer \
  -DCMAKE_PREFIX_PATH=/path/to/dspecm1d/install
cmake --build /tmp/dspec-k2-consumer --target sensitivity_kernels_consumer
python3 tutorials/sensitivity_kernels/run.py \
  --consumer /tmp/dspec-k2-consumer/sensitivity_kernels_consumer \
  --model data/models/prem.200.noatten.txt \
  --output /tmp/dspec-k2-20261023/example
```

The runner and C++ consumer both refuse a nonempty output directory before
writing. It produces a readable `report.html`, three PNG plots, full-precision
`samples.csv` and `directions.csv`, a readable `run.txt`, and a numeric JSON
companion with model and consumer paths and SHA-256 identities. The plots show reference density, solid C and fluid bulk profiles,
the corresponding density, solid-C and fluid-bulk kernels, plus kernel
predictions beside measured localized directional changes. Damping is
compensated once before applying the fixed window; the return transform uses
the physical `dt` in seconds. This is a
fixed-mesh discrete-gradient check for the selected two-bin measurement. It
does not establish mesh or frequency-band convergence, nor should its
selected-band windowed record be read as a broad-band research seismogram.
