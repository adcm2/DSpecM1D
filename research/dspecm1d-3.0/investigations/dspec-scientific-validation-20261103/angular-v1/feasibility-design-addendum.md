# Angular convergence feasibility and run design

Status: revised proposal for independent Sol review before source changes. The full-angle memory policy below supersedes the earlier chunked-visitor proposal in this file. No source, checker, or numerical-run files have been changed for this study.

## Frozen case and scope

Retain the approved Study 2 attenuated anisotropic PREM transverse-velocity traveltime fixture: source depth 29.5 km; solid-side receiver depth 3 km; 60 degree distance and 30 degree azimuth; full radial/gravity solve; FFT bins 7 through 41; the existing 2, 4, 8, and 10 mHz filter; and the existing 1000–1200 s measurement window with its 40 s taper. Each degree case uses its own waveform-derived adjoint weights. The result is specific to this filter and window. The illustrative [10,50] and [10,100] bands are not run controls.

Keep the 168 physical fixed probes and the existing 68 along-section plus 10 cross-section samples. Do not treat source-shell or contact samples as ordinary bounded pointwise convergence samples. Retain and report those raw values and trace classes separately. The current angular section spacing is set by sampling degree 100, independently of solve degree.

## Full-angle memory policy for the new candidate

The current real-Legendre route retains all angle geometry and full-angle output fields for one visitor call. Under the prior reservation, L800 and 78 angles used a declared 27,564,576 B geometry and 39,154,152 B temporary estimate for NQ6, exceeding the prior 16 MiB geometry and 32 MiB temporary reservations per worker. The user has now authorized a new per-worker policy for this isolated candidate: at most 32 MiB geometry and 64 MiB combined angular cache plus temporary work. The shared-system active-memory ceiling is 900 GB, subject to separate admission by the CCDSM/SpecM3D coordination; request a 40 GB reservation. Keep the numerical job cap at 16 GiB and the build cap at 32 GiB. No remote launch may occur before shared admission is agreed.

Retain all 78 section angles (68 along and 10 across) and use the existing unchanged full-angle real-Legendre visitor. Do not add the chunk adapter. This avoids changing visitor-call ordering or repeating radial projection/table construction. In the new candidate only, update the declared per-worker reservation guards and report both formula-derived estimates and observed RSS. At L800/NQ10 with the conservative maximum 350 radial sides, the existing formula gives 27,564,576 B geometry and 44,920,488 B angular temporary/cache use (about 42.84 MiB). The pinned double-precision Wigner constructor also has two local recursion tables of 12,848 B total at L800, in addition to the existing 160,200 B one-angle Wigner-object estimate; report this transient separately. The combined conservative peak remains below the authorized 32/64 MiB limits. NQ6 and NQ8 have no greater side count than this NQ10 bound and must also report their actual calculated values. Recompute from actual degree, angle, and node dimensions with overflow-checked arithmetic and fail closed above either cap. Include simultaneous full-angle output buffers, coefficient/Legendre tiles, Wigner construction scratch, recursion tables, and Eigen packing reserve exactly as the existing accounting does. Do not retain multiple plans or add an L-dependent history.

The all-bin section-product buffer remains separately checked: 78 × 350 × 35 × 16 = 15,288,000 B. The 168 fixed-probe products and selected-bin strain diagnostics remain under their existing independent bounds. The shared resource reservation is not a substitute for the per-worker code cap or per-job 16 GiB cap. Record resource admission, job-level RSS, and the user-space watchdog separately; do not describe the watchdog as kernel aggregate enforcement.

The fixed probes are independent of section sampling and remain at all 168 specified physical points. Keep the two existing selected forward/unit-adjoint strain bins as diagnostics; do not describe these two bins as complete strain convergence evidence.

## High-degree angular stability and action limitation

The existing `LegendreTable` constructs pinned `CachedAngularWigner` values at each requested angle and stores the degree rows through the requested cutoff. This study does not change that recurrence or its conventions. The high-degree result gate must reject non-finite modal rows, fields, spectra, weights, or kernels, and compare DSpec and YSpec on the same fixed points and section grid. The high-degree radial comparison below is required so angular changes are not attributed to an unverified radial discretization.

The current independent `WindowedSensitivityEngine::evaluate` path rejects attenuation and non-displacement output, so it does not support an independent weak derivative action for this fixture. Full three-dimensional angular quadrature at L800 would also exceed the bounded scope. No high-L tangent/action solver or source/receiver-frame pole pairing will be introduced. The action comparison is therefore explicitly unsupported in this campaign; do not imply fixed-point sums or section samples are a full 3D perturbation action.

## Initial physical queue (8 sequential jobs)

Use the accepted fixed absolute mesh-step parameter and NQ6 for the primary ladder. DSpec uses the frozen batched-radial and real-Legendre routes with 8 outer frequency workers and one inner thread. YSpec keeps its accepted radial ODE tolerance and native 8-worker route; its SEM-like output mesh does not refine the YSpec radial ODE.

1. DSpec L400, NQ6, fixed h; full-angle output.
2. YSpec L400, NQ6, fixed h.
3. DSpec L600, NQ6, fixed h; full-angle output.
4. YSpec L600, NQ6, fixed h.
5. DSpec L800, NQ6, fixed h; full-angle output under the new candidate-only memory policy.
6. YSpec L800, NQ6, fixed h.
7. DSpec L800, NQ8, fixed h; radial-control output.
8. DSpec L800, NQ10, fixed h; highest-order radial-control output.

Reuse the sealed L200 and L300 Study 2 outputs only where the frequency grid, case, physical points, sampling degree, and relevant radial controls match exactly. Require exact grid and angle coordinates, exact physical fixed-probe IDs/coordinates/trace classes, and verify actual mesh boundaries for every radial control. NQ8/NQ10 outputs are radial controls, not additional angular ladder levels. Do not add filler jobs or extend above L800.

## Cutoff and tail criteria

All target fractions are 0.001 (0.1%) and apply only to the named frozen measurement filter/window. Preserve absolute differences, reference norms, exact-zero classes, and active-bin counts without arbitrary denominator floors.

- Waveform: compare the 1000–1200 s tapered windowed transverse-velocity samples using `||w*(u_L-u_ref)||2 / ||w*u_ref||2`, where `ref` is the higher cutoff. Also compare the 99 unweighted window samples using `||u_L-u_ref||2 / ||u_ref||2` as a separate diagnostic. Retain whole-trace comparisons separately.
- Measurement: compare the waveform-derived measurement-rate normalization and the band-limited adjoint-weight vector in relative L2 over the frozen 35-bin grid. Preserve zero bins exactly.
- Kernels: the required pointwise-vector acceptance cohort is the 89 far-solid-interior probes and each predeclared depth stratum within that cohort. The depth strata are 10, 24, 29, 30, 50, 100, 300, and 500 km, with 10, 10, 10, 10, 10, 11, 14, and 14 probes respectively. Compare own-weight integrated fixed-probe vectors. For each required cohort and each candidate-versus-higher-cutoff or adjacent increment, report `||K_a-K_b||2 / ||K_b||2` and `max_i |K_a-K_b| / max_i |K_b|`, with `b` the higher cutoff; print the raw numerator and reference denominator too. If a reference denominator is exactly zero, pass only when the corresponding difference is exactly zero and report the zero class. No denominator floor is allowed. Keep material interfaces, source shell, exact coincidences, and positive-distance near-field cohorts separate as diagnostics; they are not required pointwise acceptance cohorts. Do not use a source-shell maximum as a matched one-sided trace claim.
- Weight attribution: retain each case's own P and own-weight kernel as the primary result. From the saved per-bin P and unit complex products, also compute the explicitly labelled L800-weight counterfactual using `P_L800 * C_L` at each lower degree. This separates weight changes from unit-field/product changes without new solves. Do not call the counterfactual a native kernel result.
- Per-bin tails: compare receiver spectrum values per bin across degree cutoffs. Their differences are additive receiver-output increments over the added degree ranges. Compare saved fixed-probe unit products per bin as complete product differences; they include cross terms and are not isolated additive degree-block kernel contributions. Preserve selected-bin strain records as diagnostics.
- Radial control: compare L800/NQ6 against NQ8 and NQ10 using the same fixed h, 35 frequencies, probes, and output locations. The radial error must be shown below 0.1% for the named supported outputs before using L800 as the angular reference.

A lower cutoff is certified only if its agreement with L800 is at most 0.1% for every named required waveform, normalization/adjoint-weight, and kernel-cohort measure, and the next two successive increments at or above that candidate cutoff are each at most 0.1%. Thus L400 requires both 400→600 and 600→800 to pass, and L600 cannot be certified from the single remaining increment. Report 300→400, 400→600, and 600→800, and each candidate against L800. If only the final increment passes, or any required cohort/measure fails, report that the ladder has not established a cutoff. Report DSpec and YSpec separately; do not infer a universal L cutoff from this one case.

## Gate sequence and resources

Sol reviews this addendum and the exact source delta before implementation/build. Source changes are limited to the explicitly bounded degree controls, new-candidate-only 32/64 MiB angular memory checks/metadata, and checker evidence. Parent owns the build package and orchestration. Numerical jobs remain sequential on Tunya under the existing 8-worker/one-inner-thread policy; build caps remain 8 threads/32 GiB and numerical caps 20 threads/16 GiB. The total active Tunya system ceiling is 900 GB under the separately agreed reservation; preserve the accepted at-most-50-thread aggregate policy. Preserve all prior artifacts and failures. No L>800 follow-on, production change, speed optimization, or broader frequency/window claim is in scope.
