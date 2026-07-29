# Static Love numbers

This optional module calculates elastic, static, self-gravitating Love
numbers for spherical, radially transversely isotropic models. Solids use the
five elastic coefficients `A`, `C`, `F`, `L`, and `N`; internal static fluids
are supported with gravitational potential as their only field. A fluid or
ocean at the surface is rejected.

## C++ API

Enable the module with:

```text
-DDSPECM1D_BUILD_LOVE_NUMBERS=ON
```

The public header is `<DSpecM1D/LoveNumbers>`. A calculation is:

```cpp
EarthModels::ModelInput<double> model(model_path);
DSpecM1D::LoveNumbers::Config config{
    .maximum_degree = 10,
    .polynomial_order = 6,
    .maximum_radial_step = 0.01,
};
const auto results = DSpecM1D::LoveNumbers::calculate(model, config);
```

Results are returned in ascending order from degree zero through
`maximum_degree`. `polynomial_order` is the polynomial degree `p`, so every
element has `p + 1` GLL nodes. `maximum_radial_step` is passed directly to the
existing DSpecM1D radial mesh as its relative radial step; it is not a length
in metres.

Each result contains:

| Value | Units |
| --- | --- |
| `h_u`, `h_phi`, `h_load()` | m^3 kg^-1 |
| `k_u`, `k_phi`, `k_load()` | m^4 kg^-1 s^-2 |
| `h_t` | s^2 m^-1 |
| `k_t` | dimensionless |

The generalized physical-loading response is
`h_load = h_u + h_phi`, `k_load = k_u + k_phi`.

## Command line

The same build provides:

```text
dspecm1d-love MODEL_FILE LMAX OUTPUT_FILE \
              [POLYNOMIAL_ORDER] [MAXIMUM_RADIAL_STEP]
```

The optional arguments default to `6` and `0.01`. Use `-` as
`OUTPUT_FILE` for standard output. For example:

```text
dspecm1d-love data/models/prem.200.no.noatten.txt 10 love.txt
```

Metadata lines begin with `#`. Data rows have nine columns:

```text
l h_u k_u h_phi k_phi h_load k_load h_t k_t
```

They use the units listed above and are written from `l = 0` through `LMAX`.

## Conventions and limits

- Degree zero uses continuous radial displacement `U` and potential `P`,
  retains both centre coefficients, and has exactly zero tidal response.
- Degrees one and above use `U,V,P` in solids and `P` only in fluids. `P` is
  continuous across every interface; solid displacement is not continued
  through a fluid.
- The complete applicable centre SEM space is retained without a strong
  centre condition.
- Degree one removes only the surface `P(a)` coefficient, imposing
  `P(a) = 0` exactly. Its generalized gravitational load and surface
  potential responses are zero in this convention.
- Density derivatives are evaluated inside each density-spline layer, so
  duplicated-radius interfaces retain their one-sided derivatives and
  density jumps are not volume derivatives.
- Surface fluids and oceans are unsupported. Internal fluid layers,
  including a central fluid, are supported.
- Horizontal, toroidal, dynamic, viscoelastic, and rotational responses are
  not provided.

The module owns a corrected background-gravity calculation: it integrates
the model density splines between their knots and preserves enclosed mass
across material boundaries. It does not use `MeshModel::Gravity` in the
operator or forcing.

## Validation summary

Constant-density radial and spheroidal matrices are checked against the
existing `hR()` and `hS(l)` implementations. Gated paper validation uses
pinned `gia3D` and `core` revisions for matched all-solid and internal-fluid
models, and for a three-way isotropized-PREM comparison. The corrected
module-owned gravity brings the varying-density controlled comparisons to
the solver-discretization scale.

The official Chen--Pan--Bevis ELLN TI example is also audited. Its degree-10
loading result is an independent diagnostic, not a ground-truth oracle:
ELLN and DSpecM1D use different physical gravitational constants, have
different model representations, and retain a small unresolved inter-code or
tabulation difference. Detailed provenance and diagnostics live under
`validation/`; they are enabled only by
`DSPECM1D_ENABLE_PAPER_VALIDATION`.
