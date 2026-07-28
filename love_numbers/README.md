# Static Love numbers

This module will calculate elastic, static, self-gravitating Love numbers for
radially transversely isotropic Earth models using the reference moduli
`A`, `C`, `F`, `L`, and `N`.

For unit displacement/traction, gravitational-potential, and tidal-potential
forcing, respectively, it will return
`h_u, k_u, h_phi, k_phi, h_t, k_t`, with the physical load sums
`h_load = h_u + h_phi` and `k_load = k_u + k_phi`.

The public result units are:

- `h_u`, `h_phi`, `h_load`: m^3 kg^-1;
- `k_u`, `k_phi`, `k_load`: m^4 kg^-1 s^-2;
- `h_t`: s^2 m^-1;
- `k_t`: dimensionless.

## Command-line calculator

Build with `DSPECM1D_BUILD_LOVE_NUMBERS=ON`, then run:

```text
dspecm1d-love MODEL_FILE LMAX OUTPUT_FILE \
              [POLYNOMIAL_ORDER] [MAXIMUM_RADIAL_STEP]
```

The defaults are polynomial degree 6 and maximum relative radial step 0.01.
Use `-` as `OUTPUT_FILE` to write to standard output. For example:

```text
dspecm1d-love data/models/prem.200.no.noatten.txt 10 love.txt
```

The text output has metadata lines beginning with `#`, followed by:

```text
l h_u k_u h_phi k_phi h_load k_load h_t k_t
```

`h_u`, `h_phi`, and `h_load` use m^3 kg^-1; `k_u`, `k_phi`, and
`k_load` use m^4 kg^-1 s^-2; `h_t` uses s^2 m^-1; and `k_t` is
dimensionless. Rows are ordered from degree zero through `LMAX`.
Surface fluids and oceans are unsupported. Degree one uses the exact
surface condition `P(a) = 0`.

For `l = 0`, continuous `U` and `P` fields exist at every radial node,
including in fluids, with interleaved `U, P` ordering. For `l >= 1`, `P` is
continuous everywhere, while `U` and `V` exist only on solid elements and are
shared only across solid-solid boundaries. At a fluid-solid interface the
solid side has `U` and `V`, and both sides share the single `P` degree of
freedom.

The full SEM trial space extends to `r = 0`: a solid centre retains `U(0)`,
`V(0)`, and `P(0)` for `l >= 1`, and `U(0)` and `P(0)` for `l = 0`, without
strong-form centre conditions. For `l = 1`, the surface `P(a)` degree of
freedom is omitted exactly while surface `U` and `V` are retained.

`polynomial_order` is the polynomial degree `p`, so each radial element uses
`p + 1` GLL nodes. `maximum_radial_step` is the existing DSpecM1D relative
radial step. Outermost fluid layers are rejected, while internal fluid regions
are allowed. Density derivatives are evaluated within each model layer, giving
the appropriate one-sided values at duplicated-radius material interfaces.
The Love-number radial state computes background gravity by integrating each
layer density spline exactly between its knots with three-point
Gauss--Legendre quadrature. Enclosed mass and gravity remain continuous across
element and material boundaries, without treating density jumps as volume
terms.
Exact comparisons with the legacy `hS(l)` and `hR()` matrices use
constant-density isotropic and TI test fixtures, where the legacy and
Love-number gravity calculations coincide. Varying-density gravity is tested
against analytic models and the independent controlled gia3D comparison.

The radial state, degree-of-freedom map, and degree-one-and-above TI static
spheroidal operator are implemented. The operator is assembled directly from
the weak form and retains all centre degrees of freedom. Internal fluids use
`P` only and include the `rho'/g` volume term. Both solid-fluid and fluid-solid
interface orientations are supported. Surface fluids remain unsupported, so a
PREM ocean must be removed rather than replaced by artificial solid material.
Forcing vectors, the solver, and dimensional surface extraction are implemented
privately for the three forced problems. Degree zero uses the static radial
`U, P` equations, retains `U(0)` and `P(0)`, and has zero tidal response because
its constant tidal potential has zero gradient. Degree one uses the static
spheroidal equations with the full centre space and the exact surface condition
`P(a) = 0`. Its constrained gravitational right-hand side and all surface
potential responses are zero. The public in-memory `calculate()` function
returns degrees zero through `maximum_degree` in ascending order. The
`dspecm1d-love` executable writes those results to a plain text file or
standard output.

An optional paper-validation diagnostic compares all six public components
with the PREM values in `da380/SLReciprocityGJI` for degrees 1, 2, 3, 10, 20,
50, and 100. The units, load signs, orthonormal-harmonic convention, 6368 km
solid surface, ocean removal, and degree-one `P(a) = 0` frame convention
match. The available repository deck is sampled and rounded, however, while
the reference used analytic PREM, and the exact generalized-data generation
revision is not recorded upstream. Moreover, gia3D's solid matrix is
isotropic: it derives equivalent `kappa` and `mu` from the analytic PREM
`A`, `C`, `F`, `L`, and `N` coefficients rather than assembling the full TI
operator. The existing comparison against a sampled local deck therefore does
not independently validate the full TI implementation. At polynomial order 6
and maximum radial step 0.005, nonzero relative differences range from about
3.2e-6 to 3.3e-2 and mostly plateau under refinement. The stored
`prem_love_reference.dat` values are consequently kept behind
`DSPECM1D_ENABLE_PAPER_VALIDATION` as a published-data diagnostic, not an
exact regression oracle.

The same option also builds a direct Fortran reference from pinned
`da380/gia3D` and `da380/core` revisions without adding Fortran to the normal
build. Its validation-only driver uses analytic, ocean-free PREM and gia3D's
original isotropic solid matrix and centre truncation. It reports independent
displacement, potential, and tidal solutions for degrees 1, 2, and 10, plus
stock combined-load agreement and residuals. A build-tree-only patch
initializes the fluid quadrature radius before evaluating `r^2 rho'/g`.
Published-data differences remain diagnostic and have no permanent numerical
tolerance.

Configure with `-DDSPECM1D_BUILD_LOVE_NUMBERS=ON` to build this module.
