# Love-number implementation

## Problem and outputs

The module solves the static, elastic, self-gravitating equations of a
spherical, radially transversely isotropic model independently for each
spherical-harmonic degree. Solid material is described by density, gravity,
and the five coefficients `A`, `C`, `F`, `L`, and `N`. Internal fluids are
static and carry gravitational potential only.

Three right-hand sides are solved together:

1. unit radial displacement/traction loading;
2. unit surface gravitational-potential loading;
3. unit external tidal potential.

They produce the generalized responses `h_u,k_u`, `h_phi,k_phi`, and the
tidal responses `h_t,k_t`. The physical loading sums are:

```text
h_load = h_u + h_phi
k_load = k_u + k_phi.
```

The public SI units are:

```text
h_u, h_phi, h_load    m^3 kg^-1
k_u, k_phi, k_load    m^4 kg^-1 s^-2
h_t                   s^2 m^-1
k_t                   dimensionless
```

## Radial state and gravity

`RadialState` owns the `EarthMesh::RadialMesh`, sampled `MeshModel`, and
density derivative at every element node. A configuration polynomial order
is the degree `p`; the mesh constructor therefore receives `p + 1` GLL
nodes. The module requires `p >= 3`, so the GLL rule, exact through degree
`2p - 1`, exactly integrates the degree-five product of a cubic density
spline and `r^2` on each spline segment. `maximum_radial_step` is the
existing DSpecM1D relative radial-step argument.

Density derivatives are evaluated as
`model.Density(layer).Derivative(radius)` using each element's own layer.
Duplicated-radius interfaces consequently retain the two one-sided spline
derivatives. A density jump is never inserted into a volume derivative.

The main DSpecM1D `MeshModel` background-gravity integration was corrected
in commit `207b6462ab169a86b7cd4ea183424e2a65d7c34b`. It calculates:

```text
I(r) = integral_0^r rho(s) s^2 ds
g(r) = 4 pi G I(r) / r^2,     g(0) = 0.
```

The integration proceeds outwards, splits at every density-spline knot, and
evaluates density at the quadrature radii. The accumulated density moment is
shared across element and material boundaries, so enclosed mass and gravity
are continuous even when density jumps. `RadialState::gravity` and
`surfaceGravity` delegate directly to the owned `MeshModel`; the former
duplicate private Love-number integration has been removed. The operator and
forcing continue to use the `RadialState` accessors, leaving one authoritative
production gravity calculation. The defect was identified during Love-number
validation; its obsolete output curve is no longer retained as an active
benchmark.

## Degrees of freedom

`LoveDofMap` is built by an outward traversal.

- At `l = 0`, continuous `U,P` coefficients exist at every radial node,
  including fluids, in interleaved order.
- At `l >= 1`, continuous `P` exists everywhere. `U,V` exist only on solid
  elements and are shared across solid--solid boundaries. At a solid--fluid
  boundary the solid has `U,V` and both sides share one `P`.
- Every applicable centre coefficient is retained. A solid centre has
  `U(0),V(0),P(0)` for spheroidal degrees; the radial problem has
  `U(0),P(0)`. No strong centre condition is imposed.
- At `l = 1`, only surface `P(a)` is omitted. Surface `U,V` remain.

For `NE` elements with `NQ` nodes, let
`N_r = NE (NQ - 1) + 1`. If `N_s` counts the unique nodes in each contiguous
solid run independently, the dimensions are:

```text
l = 0:   2 N_r
l = 1:   N_r + 2 N_s - 1
l >= 2:  N_r + 2 N_s
```

## Static operator

The radial `l = 0` matrix is assembled directly from the regular reduced
`U,P` weak form in both solids and fluids. It includes the exterior vacuum
term `a P_hat(a) P(a)/(4 pi G)`.

For `l >= 1`, with `q = l(l+1)`, the solid weak form uses:

```text
D = r U'
S = r V' - V + U
T = 2 U - q V.
```

The implementation assembles the complete TI elastic, self-gravity,
displacement--potential, and potential terms directly at GLL points. It does
not divide `U`, `V`, or `P` by radius, so the reduced form remains regular at
the centre.

Every element contains the potential contribution. A fluid element adds:

```text
+ r^2 rho'(r) / g(r) P_hat P.
```

At `r = 0`, the combined coefficient is set to zero before division. Away
from the centre a non-finite coefficient is rejected.

At each solid--fluid transition, the interface contribution uses the
fluid-side density, the solid-side `U`, and the shared `P`:

```text
s rho_f r^2 [g U_hat U + U_hat P + P_hat U],
s = -1  solid below, fluid above,
s = +1  fluid below, solid above.
```

The exterior vacuum term is present for `l >= 2`. At `l = 1`, surface `P(a)`
is fixed to zero by omission, so no surface-potential index or vacuum term is
requested.

## Forcing and solve

Forcing is accumulated as a load functional and negated once:

```text
K X = -L.
```

The displacement and gravitational columns place `a^2 g(a)` on surface `U`
and `a^2` on surface `P`, respectively. Degree one has no surface `P` test
coefficient, so its gravitational column is exactly zero.

For `l >= 1`, the tidal potential is `psi=(r/a)^l`. Solid volume terms use
`rho r psi (l U_hat + q V_hat)`. Fluid volume terms use
`r^2 rho'/g psi P_hat`, with the same safe centre coefficient as the
operator. Interface tidal terms use `s rho_f r^2 psi U_hat` and the same
explicit orientation convention. At `l = 0`, `psi=1` has zero gradient, so
the tidal column and response are exactly zero.

One `Eigen::SparseLU` factorization solves all three columns. Surface `U` and
`P` are converted from the model normalization to the SI units above.
Degree-one `P(a)` is constrained, so `h_phi`, `k_u`, `k_phi`, and `k_t` are
returned as zero under this convention.

## Public API and command line

`<DSpecM1D/LoveNumbers>` defines:

```cpp
struct Config {
  int maximum_degree;
  int polynomial_order = 6;
  double maximum_radial_step = 0.01;
};

std::vector<DegreeResult>
calculate(const EarthModels::ModelInput<double> &, const Config &);
```

`calculate` validates the configuration and solid surface, creates one
`RadialState`, then solves degrees zero through `maximum_degree` in ascending
order.

The private command-line writer is:

```text
dspecm1d-love MODEL_FILE LMAX OUTPUT_FILE \
              [POLYNOMIAL_ORDER] [MAXIMUM_RADIAL_STEP]
```

`-` writes to standard output. Metadata records the model, configuration,
degree-one convention, columns, and units. Each data row is
`l h_u k_u h_phi k_phi h_load k_load h_t k_t`.

## Source structure

- `include/DSpecM1D/LoveNumbers`: public configuration, results, and
  `calculate`.
- `src/RadialState.h`: mesh, sampled model, density derivative, and
  `MeshModel` gravity access.
- `src/LoveDofMap.h`: degree-dependent private field topology.
- `src/StaticOperator.h`: radial, spheroidal, fluid, and interface matrices.
- `src/StaticForcing.h`: the three load columns.
- `src/StaticSolve.h`: factorization, solve, extraction, and units.
- `src/LoveNumbers.cpp`: public validation and degree loop.
- `app/dspecm1d-love.cpp`: argument parsing and text output.
- `tests/`: focused unit, matrix-oracle, solve, API, and CLI tests.
- `validation/`: optional external-reference programs, data, and provenance.

## Validation and provenance

Normal tests cover topology, centre and degree-one treatment, one-sided
density derivatives, analytic gravity, fluid volume and interface signs,
matrix symmetry, independent forcing assemblies, residuals, linearity,
reciprocity, dimensional conversion, public results, and CLI parsing.
Constant-density matrices are compared with `hR()` and `hS(l)` so those
exact oracles isolate the operator algebra; both paths now use the corrected
`MeshModel::Gravity`.

Paper validation is gated by `DSPECM1D_ENABLE_PAPER_VALIDATION`; the normal
build requires no Fortran, MATLAB, Octave, BLAS, LAPACK, or network access.
It pins:

```text
da380/gia3D  e6a8230fc792416087e530726e3f9f9b41ee0506
da380/core   eb1d36e418ab0ac9e0b34abb717e6d9b2e863ee5
```

Pinned gia3D comparisons cover a common homogeneous sphere, both
solid--fluid interface orientations, a central fluid, and isotropized PREM.
The build-tree gia3D copy receives two documented validation-only patches:
one initializes the fluid quadrature radius and one avoids `0/0` in the
central-fluid tidal force. The original-truncation PREM driver is retained.

The external official Chen--Pan--Bevis ELLN package is identified by archive
SHA-256
`f2904364f9be96239fc85140a032b01b732cba6ea1e9bc323f2f3325834e7da7`.
It is not vendored because no redistribution licence was found. A
validation-only non-GUI patch disables waitbar operations. A separate
build-tree copy changes the sole physical `G` literal to ELLN's value for a
bounded sensitivity diagnostic. ELLN provides conventional loading
`h_l,l_l,k_l`, not the separate generalized or tidal responses returned
here. Its degree-one deformed-centre convention is not equivalent to
`P(a)=0`, so degree-one numerical values are excluded.

The ELLN degree-10 TI-PREM comparison is independent evidence, not an assumed
source of truth. ELLN uses `G=6.67384e-11` SI while production DSpecM1D uses
`6.67230e-11` SI. The constant difference explains most, but not all, of the
published-table discrepancy; the remaining representation, tabulation, or
inter-code difference is recorded as unresolved rather than tuned away.

Detailed checksums, model-generation rules, patches, and diagnostic results
are in `tests/data/PROVENANCE.md`, `validation/gia3d/PROVENANCE.md`, and
`validation/elln/PROVENANCE.md`.

## Technical report

The completed colleague-facing report is:

```text
love_numbers/report/love_numbers_report.tex
love_numbers/report/love_numbers_report.pdf
```

`love_numbers/TESTS.md` records every normal and gated paper-validation
Love-number test and supplies the generated report appendix. Benchmark CSV
files are reproduced from the existing labelled validation outputs with:

```text
python3 love_numbers/report/collect_benchmarks.py \
  --source-root . --build-dir PAPER_BUILD \
  --output-dir love_numbers/report/data \
  --elln-archive /path/to/ggx444_supp.zip
```

The vector figures and PDF are then rebuilt with:

```text
MPLCONFIGDIR=/tmp/dspecm1d-love-matplotlib \
  python3 love_numbers/report/generate_figures.py
love_numbers/report/build_report.sh
```

`love_numbers/report/README.md` gives the full fresh-build commands and maps
each stored dataset to its generating validation command.

## Limits

The module does not support a surface ocean, horizontal output, toroidal
response, rotation, dynamics, or viscoelasticity. The Love-number target and
header are build-tree interfaces only; installation and CMake package export
have not been restructured.

The main-library gravity correction and Love-number delegation are complete;
there is no separate private Love-number gravity implementation.
