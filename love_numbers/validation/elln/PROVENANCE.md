# ELLN validation provenance

## Official package

This validation audits the official supporting information for:

J. Y. Chen, E. Pan, and M. Bevis, “Accurate computation of the elastic
load Love numbers to high spectral degree for a finely layered,
transversely isotropic and self-gravitating Earth,” *Geophysical Journal
International* 212(2), 827–838 (2018),
<https://doi.org/10.1093/gji/ggx444>.

The publisher archive is `ggx444_supp.zip`, obtained from the Oxford
Academic supporting-information link. Its stable publisher path is:

```text
https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/gji/212/2/10.1093_gji_ggx444/3/ggx444_supp.zip
```

The CDN may require a temporary query-string signature. The validation
therefore accepts an explicit archive or extracted source directory and
does not download or vendor the package.

Archive SHA-256:

```text
f2904364f9be96239fc85140a032b01b732cba6ea1e9bc323f2f3325834e7da7
```

The outer archive is dated 2017-10-16 and contains:

| File | Purpose | SHA-256 |
| --- | --- | --- |
| `2017 1010 GJI ELLN Supplement FINAL Clean.doc` | Manual and published tables | `cf8174bbb56d06eb53f99f06a2f538da52654785815bc70d70d05924f0ee83a7` |
| `ELLN Direct Run.zip` | Non-GUI MATLAB source and README | `f4c4b10b84d1dfa08d5f4187fffebf9c50d5a860d3d22eeebd801f5befb08f7a` |
| `ELLN Example Input Files.zip` | Mantle, core, and benchmark models | `c5b138d3806c4d6590cd8d417d21b1195315072ff6fd0f02d68cd6ffe9cb0798` |
| `ELLN GUI Run.zip` | GUI installer, source, and GUI README | `56d00b6ed524d0896768bd73aaabe2d5bdb95f280e7e192568f835c9b8799fb4` |

The direct archive contains `ELLN.m`, `ELLN100M.m`, and its README. The
model archive contains `EarthCore26.txt`, `EarthMantle56.txt`,
`EarthMantleTI56.txt`, two benchmark decks, two hard/soft decks, and a
one-metre surface-layer deck. The GUI archive contains the MATLAB app
installer, `ELLN.fig`, `Control.m`, and its README. It contains no
standalone reference-result file; the official numerical values are
tables in the manual.

The files used by the automated audit have these checksums:

| File | SHA-256 |
| --- | --- |
| `ELLN.m` | `6341a62a0c595edfea15a81e0d8ed7d79f601b1ebd0c10448926b66bd0648864` |
| `EarthMantleTI56.txt` | `43d455cf329b7d9db6684b434447b648e3a5311f84873bb283e2a8f22f8f17fe` |
| `EarthCore26.txt` | `b6ddf75db09d15b4d7e065f109134d54cdd68434e7714b75839d1c12c49d9f12` |

The supplied GUI README requires MATLAB 2014 or newer and says that it
was tested on 32-bit and 64-bit Windows. No licence, redistribution
grant, copyright notice, or warranty statement was found in the
manual, READMEs, MATLAB source, models, GUI source, or installer. The
package is consequently kept external.

## Model and mathematical conventions

An isotropic mantle row is

```text
ID  radius(m)  density(kg/m^3)  mu(Pa)  lambda(Pa)
```

and a transversely isotropic mantle row is

```text
ID  radius(m)  density(kg/m^3)  A=c11  C=c33  L=c44  N=c66  F=c13
```

Thus the ELLN-to-DSpecM1D mapping is exactly
`A -> A`, `C -> C`, `F -> F`, `L -> L`, and `N -> N`. For an isotropic
row, `A=C=lambda+2*mu`, `F=lambda`, and `L=N=mu`. A layered core row is

```text
ID  radius(m)  density(kg/m^3)  lambda(Pa)
```

and represents a fluid with zero shear modulus.

For the shell between rows `i` and `i+1`, ELLN uses the properties in
row `i+1`. Within an analytic mantle shell the source uses
`rho(r)=rho[i+1]*r_bar/r`, where `r_bar` is the mass-equivalent radius
formed by `2(r2^3-r1^3)/(3(r2^2-r1^2))`. Elastic coefficients and the
shell gravity are constant at their upper-interface values. Layered
core density is constant per shell and core gravity has the form
`g=k*r`. Gravity is calculated internally from the supplied density by
Newton’s law; the code uses `G=6.67384e-11 m^3 kg^-1 s^-2`.

The supplied `EarthMantleTI56` and `EarthCore26` case has 56 mantle
shells, 26 fluid-core shells, a core-mantle radius of 3,480,000 m, and
a surface radius of 6,371,000 m.

ELLN uses radial displacement positive outward and tangential
displacement positive away from the load centre. Its surface-loading
boundary conditions use a downward traction and the direct potential
of the surface mass. The paper defines all responses using the same
scalar spherical-harmonic load coefficient, so rescaling that harmonic
does not change the reported Love numbers. The software is
degree-based and its examples are axisymmetric.

Degree zero has `l_0=k_0=0`; `h_0` is calculated for a compressible
mantle and set to zero for an incompressible mantle. At degree one the
source subtracts its initially calculated `k_1` from `h_1`, `l_1`, and
`k_1`, giving the Farrell deformed-centre convention and exact
`k_1=0`.

ELLN returns conventional, dimensionless physical loading Love numbers
`h_l`, `l_l`, and `k_l`. It does **not** return the separate generalized
DSpecM1D responses `h_u`, `k_u`, `h_phi`, and `k_phi`, and it does not
calculate tidal Love numbers. Its text output columns are

```text
l  -h_l  l_l  -k_l  l*l_l  -l*k_l
```

The signs in the text file are display conventions; the MATLAB return
values are `h_l`, `l_l`, and `k_l`.

## Dimensional conversion audit

For unconstrained degrees `l >= 2`, the paper’s load boundary
conditions and surface definitions, together with the `+1`
direct-load-potential term visible in `ELLN.m`, give

```text
h_l = (2l+1) g(a) h_load / (4 pi G a)
k_l = -[1 + (2l+1) k_load / (4 pi G a)].
```

Here DSpecM1D `h_load=h_u+h_phi` has units `m^3 kg^-1`,
`k_load=k_u+k_phi` has units `m^4 kg^-1 s^-2`, `g(a)` has units
`m s^-2`, and `G a` has units `m^4 kg^-1 s^-2`; both converted values
are dimensionless. The outer minus sign and direct-potential `+1` are
therefore intentional. DSpecM1D currently has no public horizontal
surface response, so no `l_l` conversion is possible.

Degree one is separate. ELLN first calculates raw `h_1`, `l_1`, and
`k_1` with the generic definition. It then executes

```text
h_1 <- h_1 - k_1
l_1 <- l_1 - k_1
k_1 <- k_1 - k_1
```

to select its deformed-centre frame. DSpecM1D instead omits the surface
`P(a)` degree of freedom, so its dimensional `k_load` is zero. Applying
the generic formula to that value gives `k_1=-1`, not ELLN’s
frame-adjusted `k_1=0`. Exact frame equivalence is therefore not
established: both `h_1` and `k_1` are excluded from the Stage 16
numerical comparison, and the production degree-one calculation is
unchanged.

## Stage 16 TI PREM deck generation

`export_ti_prem.py` reads the audited external `EarthCore26` and
`EarthMantleTI56` files. The 27 core density regions, including the
central sphere, are constant-density fluids. The 56 mantle regions use
the ELLN analytic law

```text
rho(r) = rho_upper * r_bar / r
r_bar = 2 (r2^3-r1^3) / [3 (r2^2-r1^2)].
```

For every sampled mantle knot it writes

```text
VPH = sqrt(A/rho)       VPV = sqrt(C/rho)
VSH = sqrt(N/rho)       VSV = sqrt(L/rho)
eta = F/(A-2L).
```

The core has `VSH=VSV=0` and `VPH=VPV=sqrt(lambda/rho)`. Every physical
boundary is duplicated. Each region has at least four uniformly spaced
knots, with maximum requested spacings 100, 50, 25, 12.5, and 6.25 km.
Rows use scientific notation with 17 digits after the decimal point.
The generated decks remain in the build tree.

Deterministic model checksums are:

| Maximum spacing | Knots | SHA-256 |
| --- | ---: | --- |
| 100 km | 333 | `28145d7ed65b0fe7441af15a1e8995fe60648d03776ecc7d689d8c506593e719` |
| 50 km | 345 | `21e02cf9098b12502edaf0e08a2f9dc86c983471942f3c03ec3e65a5aa2df214` |
| 25 km | 401 | `3ba05031e9960fa340a2c48d044f18cc4e632e12991c98ee9100f0a93b1cff14` |
| 12.5 km | 612 | `7d09b746bfa222981207fd47d75cea79b91ad473a4af9936b58f95bc6b8f6a97` |
| 6.25 km | 1119 | `e182e01b6957bbf9eb06291d63e6527a109d62f22e7b42ddebc4a3f651cfb88f` |

## Stage 16 comparison

The constants are deliberately not reconciled:

```text
G_ELLN     = 6.67384e-11 m^3 kg^-1 s^-2
G_DSpecM1D = 6.67230e-11 m^3 kg^-1 s^-2
(G_DSpecM1D-G_ELLN)/G_ELLN = -2.307517111587586e-4.
```

For the same density moment, the corresponding surface gravities are
`9.825602805567190` and `9.823335530906640 m s^-2`. DSpecM1D’s sampled
gravity approaches the latter as the density knots are refined:

| Spacing | DSpecM1D gravity | Error from DSpecM1D analytic gravity |
| --- | ---: | ---: |
| 100 km | 9.823340437718660 | 4.9068120e-6 |
| 50 km | 9.823340437718660 | 4.9068120e-6 |
| 25 km | 9.823338499990850 | 2.9690842e-6 |
| 12.5 km | 9.823336370545258 | 8.3963862e-7 |
| 6.25 km | 9.823335641111617 | 1.1020498e-7 |

At radial step `0.00125`, the degree-10 comparison with official
Table S7 (`h_10=-1.451586`, `k_10=-0.07066685`) is:

| Spacing | DSpecM1D h | Signed difference | DSpecM1D k | Signed difference |
| --- | ---: | ---: | ---: | ---: |
| 100 km | -1.451095436174441 | 4.9056383e-4 | -0.070643105081278 | 2.3744919e-5 |
| 50 km | -1.451095436174441 | 4.9056383e-4 | -0.070643105081278 | 2.3744919e-5 |
| 25 km | -1.451095052261498 | 4.9094774e-4 | -0.070643100264793 | 2.3749735e-5 |
| 12.5 km | -1.451094603499689 | 4.9139650e-4 | -0.070643090254721 | 2.3759745e-5 |
| 6.25 km | -1.451094469448399 | 4.9153055e-4 | -0.070643084436717 | 2.3765563e-5 |

With the 6.25 km deck, refining the radial step from `0.005` through
`0.0025` to `0.00125` changes `h_10` by `2.45e-10` then `6.18e-12`,
and changes `k_10` by `4.25e-12` then `1.65e-12`. The final relative
differences from Table S7 are `3.38616e-4` for `h_10` and `3.36304e-4`
for `k_10`. Thus the SEM mesh is converged at the reported scale; the
remaining difference includes the known gravitational-constant
mismatch and differences between ELLN’s layerwise analytic treatment
and DSpecM1D’s spline-sampled model. No acceptance tolerance is attached
to this diagnostic.

For the three radial steps, actual element counts are 254, 458, and
845, with 1525, 2749, and 5071 unique radial nodes. Degree-10 dimensions
are 3003, 5415, and 9825; degree-one dimensions are one smaller. The
largest reported solve residual was `2.74e-10`, and reciprocity errors
were at most `7.26e-14`.

## Stage 17 gravitational-constant isolation

The paper-validation build copies the six private Love-number
implementation files, `controlled_love_numbers.cpp`, and the main-library
`MeshModel.h` into `elln-gravity-copy` in the build tree. It places that
build-tree root first on the ELLN-`G` target's include path, so
`RadialState` resolves the copied `MeshModel.h`; normal targets continue
to resolve the production header. It then applies
`patches/elln_gravitational_constant.patch`, whose source changes are

```text
6.67230e-11 / model.GravitationalConstant()
    -> 6.67384e-11 / model.GravitationalConstant()

6.67230 * 10^-11 / model.GravitationalConstant()
    -> 6.67384 * 10^-11 / model.GravitationalConstant().
```

The first value defines
`RadialState::dimensionlessGravitationalConstant()` for the static
operator. The second makes the delegated `MeshModel::Gravity` use the
same validation constant; forcing terms obtain that gravity through
`RadialState`. The Stage 17 comparison separately uses the matching
physical constant in the dimensional-to-conventional `h_l` and `k_l`
conversion. All model values and normalization scales are unchanged.
No source-tree implementation is patched.

For the 6.25 km deck, polynomial order 6, degree 10, and maximum radial
step 0.00125, the two calculations are:

| Quantity | DSpecM1D `G` | ELLN `G` |
| --- | ---: | ---: |
| `h_load` (m3 kg-1) | -3.757599097949508e-5 | -3.758466108007903e-5 |
| `k_load` (m4 kg-1 s-2) | -2.364047912268505e-4 | -2.364554266511174e-4 |
| surface gravity (m s-2) | 9.823335641111617 | 9.825602915797605 |
| conventional `h_10` | -1.451094469448399 | -1.451429287897068 |
| conventional `k_10` | -0.0706430844367174 | -0.0706585225743535 |
| maximum solve residual | 3.43550e-11 | 3.46760e-11 |
| reciprocity error | 3.78312e-14 | 4.43162e-14 |

Against Table S7 (`h_10=-1.451586`, `k_10=-0.07066685`), the default-`G`
signed differences are `+4.91531e-4` and `+2.37656e-5`; the ELLN-`G`
differences are `+1.56712e-4` and `+8.32743e-6`. Thus changing `G`
explains 68.12% of the signed `h_10` discrepancy and 64.96% of the
signed `k_10` discrepancy. The sensitivities, expressed as response
change divided by relative change in `G`, are `-1.45066` for `h_10`
and `-0.0668882` for `k_10`.

The exporter audit reconstructs `rho,A,C,F,L,N` directly from the
official core and mantle tables and compares them with values recovered
from the generated DSpecM1D deck. It checks a fluid-core point, lower-
and upper-mantle points, both sides of the core-mantle boundary, and
both sides of a mantle boundary at 5,871 km. It confirms constant
properties in each core region, upper-row mantle property ownership,
`rho=rho_upper*r_bar/r`, duplicated boundary ownership, and SI units.
Density agrees exactly at the audited points; the largest relative
elastic-coefficient reconstruction error is `6.61e-16`. No exporter
transcription error is demonstrated.

The different gravitational constant therefore explains most, but not
all, of the original discrepancy. With the bounded audit finding no
model-export mismatch, the residual is classified as an ELLN
representation or tabulation difference and/or an unresolved inter-code
difference. This diagnostic does not identify either implementation as
the source of the residual and deliberately introduces no acceptance
tolerance or production change.

## Reproduction status

The official TI PREM example is manual Case II:
`EarthMantleTI56`, `EarthCore26`, `MantleType=3`, and `CoreType=3`.
The validation driver requests degrees 1, 2, 10, and 50. The publisher
Table S7 contains rows 1 and 10 but not rows 2 and 50; the latter two
are reported without manufacturing reference values. Table S7 also
prints inconsistent `l_l` and `l*l_l` exponents at degrees 30,000,
40,000, and 50,000. The transcription preserves those published
values and the parser reports the inconsistency.

The only compatibility patch disables the three graphical waitbar
operations for a non-GUI run. It is stored in `patches/non_gui.patch`
and is applied only to a build-directory copy of `ELLN.m`.

On the Stage 15 validation host, neither MATLAB nor GNU Octave was
installed. The official archive, source, models, and Table S7 were
therefore parsed and audited, while live reproduction was reported as
unavailable. A configured future build prefers MATLAB, then GNU
Octave; Octave compatibility remains conditional on a successful run
of this exact example.
