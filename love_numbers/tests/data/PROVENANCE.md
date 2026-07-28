# PREM Love-number reference

The selected values in `prem_love_reference.dat` are copied from
`data/love.dat` in
<https://github.com/da380/SLReciprocityGJI> at commit
`09330992b8fc47e2e26d81f0551cf53d9b60369a`. The data file was last changed
by commit `9b66a5f1511d91579777d27a0150c075b60a6040`.

The original 4097-row file has SHA-256
`0123d90fb4081598fea88158203a3b64f9274606f3976873fdfd3063b4740af7`.
Only degrees 1, 2, 3, 10, 20, 50, and 100 are copied here.
The selected copy has SHA-256
`66e3d01286cba1300f69c153b3ada6889f70f0c9c5d40d0caa5500824b2c0181`.

`SLmod.py` assigns the columns, and the `gia3D` output conversions establish
their units:

1. spherical-harmonic degree;
2. `h_u`, in m^3 kg^-1;
3. `k_u`, in m^4 kg^-1 s^-2;
4. `h_phi`, in m^3 kg^-1;
5. `k_phi`, in m^4 kg^-1 s^-2;
6. `h_t`, in s^2 m^-1;
7. `k_t`, dimensionless.

The reference README states that the values were calculated for PREM using
<https://github.com/da380/gia3D>, inspected at commit
`e6a8230fc792416087e530726e3f9f9b41ee0506`. Its `love_numbers.f90` selects
`elastic_PREM(.false.)`, which removes the 6368--6371 km ocean and places the
solid surface at 6368 km. That `gia3D` commit pins its `da380/core` submodule
to commit `eb1d36e418ab0ac9e0b34abb717e6d9b2e863ee5`.

The upstream model has twelve retained analytic PREM layers: a solid inner
core, a fluid outer core, and ten solid mantle/crust layers. Density and the
five TI moduli are obtained from analytic PREM density, velocity, and
anisotropy polynomials. The gia3D solid matrix is nevertheless isotropic: it
uses equivalent `kappa` and `mu` derived from `A`, `C`, `F`, `L`, and `N`.
Its normalisation constants are 6371000 m, 5.972e24 kg, and 3600 s, and it
uses `G = 6.6723e-11` SI. `SLmod.py` records a 6368000 m solid surface and
surface gravity 9.825652323 m s^-2.

The comparison uses this repository's existing
`data/models/prem.200.no.txt`, which also ends at 6368 km. It is a sampled
deck rather than the analytic upstream PREM implementation, so residual model
sampling and rounding differences must be distinguished from discretisation
error. That deck has SHA-256
`4e64ef450a0304e6846a453d5d02ed45355fac853ed1f7a4cbad03cf6716dd9f`.
Because the model representations differ and gia3D does not assemble the full
TI solid operator, this comparison is not an independent validation of the
full TI implementation.

`SLmod.py` uses orthonormal spherical harmonics. The Love numbers are ratios
of response and forcing coefficients of the same harmonic, so no additional
normalisation factor is applied. Upstream loads use negative surface
right-hand sides, matching `K X = -L` here. For degree one, upstream matrix
assembly omits the surface potential coefficient and its load/tide extraction
sets surface potential to zero, matching the `P(a) = 0` convention used here.

The inspected `gia3D` program writes the conventional three-component loading
or tidal result, while `love.dat` contains the six generalized components.
The exact program and `da380/core` revision used to generate that generalized
file are not recorded in `SLReciprocityGJI`. Combined with the
analytic-versus-sampled PREM difference, this prevents an exact regression
comparison. The paper-validation test is therefore diagnostic and reports
every signed and relative difference without accepting them against a
tolerance. `prem_love_reference.dat` is a published diagnostic, not an exact
regression oracle.
