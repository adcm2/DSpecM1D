# DSpecM1D baseline numerical conventions

These conventions describe the Phase 0 source tree. They define what later
dependency and documentation changes must preserve unless a separate change
is approved.

## Model and radial discretization

- `prem_norm` uses length 1000 m, density 5515 kg/m3 and
  `T = 1/sqrt(pi * 6.67230e-11 * 5515)` seconds. Mass is density times length cubed.
  Source: `DSpecM1D/src/NormClass.h`.
- `ModelInput` splits duplicated deck radii into separate layers and builds
  natural cubic splines within each layer. It derives `A,C,F,L,N` from density,
  velocities and eta. Sampling must retain each side of a discontinuity.
  Source: `DSpecM1D/src/model_info/ModelInput.h`.
- Elements use GLL nodes and weights. Every layer boundary is retained in the
  radial mesh. Cardinal basis derivatives come from Interpolation.
- `MeshModel` samples material fields by element layer and integrates enclosed
  mass for gravity. Gravity integration splits intervals at available density
  knots. Preserve this algorithm and its constants during dependency migration.
  Source: `DSpecM1D/src/model_info/MeshModel.h`.
- Full spheroidal indexing holds three fields; Cowling holds U/V. The fluid-solid
  interface maps allow tangential slip. Toroidal assembly selects the solid
  shell above the deepest fluid region. These are existing scopes rather than
  a claim of support for every planetary topology.
  Sources: `SEM/SEMLG.h` and `SEM/SEMConstructor.h`.

## Angular representation and observation

- Wigner tables use orthonormal normalization, all signed orders, and upper
  indices including -1, 0 and 1. Receiver colatitude is `90 degrees - latitude`;
  longitude enters through `exp(i*m*phi)`.
- In `ParamInfo::rvFullSph`, the three angular rows are `Y0`, `Y+ - Y-` and
  `-i*(Y- + Y+)`. Toroidal horizontal rows are `i*(Y- + Y+)` and `Y+ - Y-`.
  Preserve these signs, ordering and normalizations when replacing GSHTrans
  types. Do not infer geographic component labels solely from their array indices.
- Spectra have three rows per receiver in parameter-file receiver order.
  Source/receiver radial interpolation and the CMT coefficient routines are
  part of the observable, rather than a postprocessing detail.
  Sources: `ParamInfo.h`, `SEMForce*.h`, `SEMReceivers.h` and `SourceInfo.h`.

## Frequency and output

- The solver uses complex angular frequency `w_real - i*ep`, with the real
  frequency axis converted through the model time scale.
- Dynamic operators have the form `H - w*w*P`, with the attenuation contribution
  added when enabled. Its scalar factor is `i + (2/pi)*log(w_real/w0)`, where
  `w0 = 2*pi/tref`. Source: `SpecHelpers.h` and the spectrum solver headers.
- Output factors are `-i/w` for displacement, 1 for velocity and `i*w` for
  acceleration. The public CLI also applies the corresponding dimensional
  length/time scaling from `InputParametersNew::normFactor()`.
- The default filtering options are pre-taper 0.05, final taper 0.5, one pass
  and `enforceRealSignal=false`. Frequency-to-time conversion includes the
  damping correction and frequency-step normalization. Preserve the full chain.
  Sources: `SignalFiltering.h` and `frequency_info/PostprocessFunctions.h`.
- The public writer emits a semicolon-delimited table: time in the first
  column, then solver rows. Parameter files are ordered value lines rather
  than keyed fields. Sources: `OutputWriters.h` and `InputParser.h`.

## Solver choices

- The baseline uses Eigen SparseLU. A positive Cowling cutoff uses full gravity
  below the actual physical-frequency threshold and Cowling at or above it;
  zero disables the cutoff.
- Single-SEM and multi-SEM calculations retain their existing radial truncation
  and symbolic-factorization cadence. Full/Cowling regions have separate states.
- Multi-SEM accumulates by degree/chunk. Parallel summation may differ at
  roundoff scale, so bitwise equality across thread counts is not a contract.
- The baseline observations exercise selected damped PREM calculations. Kernel
  implementation, lateral structure and derivative conventions are later work.
