# Direct gia3D validation

The validation build uses:

- `da380/gia3D` commit
  `e6a8230fc792416087e530726e3f9f9b41ee0506`;
- `da380/core` commit
  `eb1d36e418ab0ac9e0b34abb717e6d9b2e863ee5`.

The sources are fetched from their HTTPS GitHub repositories unless
`DSPECM1D_GIA3D_SOURCE_DIR` and `DSPECM1D_GIA3D_CORE_SOURCE_DIR` name local
checkouts at those exact revisions. No upstream source is stored here.

`patches/initialize_fluid_radius.patch` is applied only to the copied
`module_matrix.f90` in the build tree. It initializes `rt` from the current
quadrature radius before gia3D evaluates `rt*rt*drho/g` in a fluid element.
The upstream checkout remains unchanged.

Only the modules needed for analytic PREM, the radial mesh, the static
spheroidal matrix, and its loading and tidal forces are compiled. The direct
driver retains gia3D's original degree-dependent centre truncation.

Each paper-validation build writes `gia3d-build-provenance.txt` beside its
validation build files. That generated record contains the resolved source
paths and revisions, Fortran compiler and version, build flags, runtime-check
flags, BLAS and LAPACK libraries, and applied patch path.

The unpatched diagnostic target is excluded from the normal build graph. With
GNU Fortran checks `-O0 -g -fcheck=all -finit-real=snan -fbacktrace`, its
degree-two run reaches the solve and is rejected by the driver's finite-value
check. Repeating those checks with the build-tree patch produces finite
results. The passing paper-validation test does not depend on the unpatched
program failing.

## Controlled homogeneous sphere

The full-domain comparison uses two text representations of the same
homogeneous solid sphere:

- radius: `6371000 m`;
- density: `5500 kg m^-3`;
- P-wave speed: `10000 m s^-1`;
- S-wave speed: `5000 m s^-1`;
- knots: `0`, `6371000/3`, `2*6371000/3`, and `6371000 m`.

`validation/data/homogeneous_sphere.dspec` is the nine-column DSpecM1D
deck. It sets `VPH = VPV`, `VSH = VSV`, `eta = 1`, `Qkappa = 1000`, and
`Qmu = 600`. Its SHA-256 checksum is
`bf5073b8185723bce8ac2326385d5729d7f75ca21a2d86e6a8f0484c9c58731a`.

`validation/data/homogeneous_sphere.gia3d` is gia3D's four-column elastic
DECK format: radius, density, VP, and VS in SI units. Its SHA-256 checksum
is
`214fee40808c4fbec970265759999f02cab2ee27ac030aff6e08f631781fb22a`.
The constant values are deliberate so that both cubic interpolators
reproduce the model exactly.

`full_domain_matrix.f90` is a validation-only adapter around the pinned
gia3D matrix routines. It accepts an explicit starting radius, calls
`build_boolean_spheroidal`, delegates element contributions to gia3D's
existing solid and fluid layer assembly routines, and retains gia3D's
banded factorisation. The controlled driver passes a starting radius of
zero, retaining the complete centre space and gia3D's degree-one
surface-potential elimination. The original analytic-PREM driver and its
degree-dependent centre truncation are unchanged.

The comparison uses five GLL nodes in each code, with DSpecM1D polynomial
degree `p = 4`. It runs `N = 4, 8, 16, 32` uniform elements using the exact
binary fractions `maximum_radial_step = drmax = 1/N`, and compares degrees
1, 2, and 10. The validation script reports model, mesh, dimension,
residual, component, and successive-mesh differences without imposing an
inter-code numerical tolerance.

## Controlled internal-fluid spheres

Stage 12 adds two paired SI decks using the same radius
`a = 6371000 m` and gravitational constant `G = 6.6723e-11 SI`.
The DSpecM1D files use nine columns with `VPH = VPV`, `VSH = VSV`,
`eta = 1`, `Qkappa = 1000`, and `Qmu = 600`. The gia3D files use its
four-column radius, density, VP, VS format.

The solid-fluid-solid model has:

- `0 <= r <= a/4`: `rho = 7000 kg m^-3`, `VP = 11000 m s^-1`,
  `VS = 5500 m s^-1`;
- `a/4 <= r <= a/2`: `rho` decreases linearly from `6000` to
  `5000 kg m^-3`, `VP = 9000 m s^-1`, `VS = 0`;
- `a/2 <= r <= a`: `rho = 4500 kg m^-3`, `VP = 9000 m s^-1`,
  `VS = 4500 m s^-1`.

Duplicated entries at `a/4` and `a/2` preserve both sides of each
material interface. The fluid-side densities are respectively `6000`
and `5000 kg m^-3`; their interface orientations are `-1` and `+1`.
The exact fluid density derivative is
`-4000/a = -6.27844922304190826e-4 kg m^-4`.

The central-fluid model has constant `rho = 6000 kg m^-3`,
`VP = 9000 m s^-1`, and `VS = 0` from the centre to `a/2`, followed by
the same constant outer solid as above. Its fluid density derivative is
zero and its fluid-side interface density is `6000 kg m^-3`.

The paired deck checksums are:

- `solid_fluid_solid.dspec`:
  `9c4ea842a94c1bb18018415a1a3cdc83eb01101511ab3366aacf83108ebca453`;
- `solid_fluid_solid.gia3d`:
  `d184210dc28adf9afc1606d1befffbbc0f097f96f316edc7199e6da1c58a225d`;
- `central_fluid.dspec`:
  `3db06b0148ede2ba121184b57e00bd56a0827e9a560e174c2be890d30bc99574`;
- `central_fluid.gia3d`:
  `0b3360d68add5914fb62f9e94f78c72f60e18b10f145f06454dc438e8f57464f`.

Both comparisons use five GLL nodes, DSpecM1D polynomial degree `p = 4`,
total element counts `N = 8, 16, 32`, and degrees 1, 2, and 10. They
reuse the Stage 11 full-domain adapter and comparison program.

The first runtime-checked central-fluid run independently found a second
pinned gia3D centre defect: `force_for_unit_harmonic_tide_fluid_layer`
formed `drho*r^(l+2)/g` directly at `r = g = 0`. The right-hand side
became non-finite and the IEEE invalid flag was set. The validation build
therefore applies `patches/safe_central_fluid_tide.patch` only to a copied
`module_force.f90`. It sets the combined coefficient to zero when
`g = 0`, before division. The upstream checkout remains unchanged.

## Controlled isotropized PREM

Stage 14 uses the same pinned `gia3D` and `core` revisions above and calls
`elastic_PREM(.false.)`, retaining its 12 solid-surface layers and
`6368000 m` outer radius. The validation-only
`export_isotropic_prem.f90` program accepts:

```text
export_isotropic_prem maximum_knot_spacing_m dspec_file gia3d_file
```

For every solid-layer knot it evaluates the analytic PREM coefficients and
forms

```text
kappa = (C + 4 A - 4 N + 4 F) / 9
mu    = (C + A + 6 L + 5 N - 2 F) / 15
VP    = sqrt((kappa + 4 mu / 3) / rho)
VS    = sqrt(mu / rho).
```

The outer core retains its analytic fluid `kappa` and uses `VS = 0`.
Each layer is divided into the smallest number of equal intervals that does
not exceed the requested maximum spacing, with at least three intervals
(four knots). Every material radius is written once for the layer below and
once for the layer above. The paired files therefore have identical radius,
density, VP, and VS values. The DSpecM1D file additionally uses
`VPH = VPV`, `VSH = VSV`, `eta = 1`, `Qkappa = 1000`, and `Qmu = 600`.

The paper-validation build generates paired models for maximum knot spacings
`100 km`, `50 km`, `25 km`, and `12.5 km` under
`isotropic-prem-models/`. The comparison program reports the SHA-256 checksum
of every generated file, along with the exact per-layer and total knot counts.
The validated GNU Fortran 13.3.0 generation produced:

- `100km`, 95 knots:
  DSpecM1D `688d59063a7ce3a4b5c66cf23089fb9444520aa4491ffd22828628f1c9864337`,
  gia3D `0dc055ef8109627240540fe9ce775eaba29bc22beb36e8e0f6b9c1bf56dd1df4`;
- `50km`, 154 knots:
  DSpecM1D `b9b092aefe80be327f503387d9d1818af34a92d1c462d137062bf4605e48d69e`,
  gia3D `46f7e1afe58287918674cd8de95af033e2b1ee77bdea268769c4fda2a36f3070`;
- `25km`, 278 knots:
  DSpecM1D `b3fba4406d4f4794e7fa93498b265bfe9642552ab7f3c24f1d26acf111aa22a4`,
  gia3D `2c3f80b199158ac727f87a3bf3c1f712cf43e450f6484fb491ea5e79600be769`;
- `12.5km`, 532 knots:
  DSpecM1D `9fcb12eb8434ecb22a4ae90fa8ddefa2a904c3915e3ec72490363e65109556e0`,
  gia3D `a6343309b6cf9381535c96dd4ad1ebee0dc170a8e891d49489bd9766fc4e5342`.

The full-domain three-way comparison uses five GLL nodes, DSpecM1D polynomial
degree `p = 4`, relative radial steps `0.01`, `0.005`, and `0.0025`, and
degrees 1, 2, 10, and 50. It reports separately:

- DSpecM1D sampled deck minus gia3D sampled deck;
- gia3D sampled deck minus gia3D analytic PREM;
- DSpecM1D sampled deck minus gia3D analytic PREM.

Component values, surface radius and gravity, layer element counts, radial
nodes, total degrees of freedom, three-column solve residuals, reciprocity,
mesh convergence, and knot-spacing convergence are diagnostic only. No
inter-code numerical tolerance is imposed. This comparison is isotropic and
does not validate DSpecM1D's full transversely isotropic operator.
