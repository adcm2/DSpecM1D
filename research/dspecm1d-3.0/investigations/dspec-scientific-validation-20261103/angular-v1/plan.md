# Angular convergence at L400, L600 and L800

Authorization: 3 November 2026. The user requests Luna high implementation, independent Sol 6.1 high review, isolated Tunya execution, and determination of the angular cutoff that meets the 0.1% target. The requested L400/L600/L800 ladder is authorized through final results review; no production integration, further optimization or automatic extension beyond L800.

## Filter- and window-specific result

The user explicitly notes that a cutoff depends on frequency filtering and the measurement window: passing one frequency band does not certify another. Their [10,50] and [10,100] examples do not replace the existing fixture or specify units for a new case.

For this campaign retain the accepted attenuated anisotropic PREM transverse-velocity traveltime case, source depth29.5km, solid-side receiver at3km, distance60degrees/azimuth30degrees, source tensor, full gravity/full radial operator, frequency grid35bins7..41, 2/4/8/10mHz spectral filter, 1000..1200s measurement window with40s taper, damping and normalization definition. Each cutoff uses its own unperturbed waveform-derived adjoint weights. Save the actual frequency/filter/time-window metadata and weights in the acceptance record. Any accepted cutoff is specific to that case. No broader frequency band, different time window or original displacement-energy observable is implicitly certified.

## Fresh implementation and review gates

Seed a new isolated candidate from the frozen Study2 candidate and record the complete seed manifest before any edit. Preserve all earlier candidates, raw evidence, reports and failed attempts. Keep the accepted LU/modal route and batched radial/real-Legendre reconstruction; source changes must be limited to the explicitly bounded high-L study controls and necessary bounded diagnostic evaluation/output.

1. Luna audits supported APIs and high-L cache/temporary bounds, then proposes an exact design and initial numerical queue. Sol reviews the plan before implementation. Do not promise a directional action or tail diagnostic that the actual APIs do not support.
2. Luna implements only the reviewed scope and logs each source/checker change immediately. Sol independently reviews source, guards, high-degree recurrence behavior, workspace bounds and exact build inputs before build.
3. Parent stages the exact minimal reviewed source/build package into fresh Tunya roots and launches a detached supervised build. Sol reviews actual command/binary/dependency/cleanup evidence before numerical staging.
4. Review the exact run/checker plans and minimal input closure before launch. Sequential physical jobs run under the existing supervisor. A standalone checker reads an immutable finalized producer snapshot. Sol reviews actual raw data, calculations, output hashes, control validity and acceptance claims.
5. Produce tables and readable great-circle DSpec/YSpec/difference wedges using linear percentage difference scales with declared denominators/display saturation. Publish only evidence-supported conclusions, then stop.

The user subsequently authorized a new shared memory policy: total active Tunya RAM must remain below900GB, coordinated with the CCDSM and DSpecM3D Implementation chats. In the fresh isolated candidate only, angular geometry may use32MiB/worker and the combined angular cache+temporary workspace64MiB/worker, with explicit runtime bounds and metadata. This modest reviewed increase permits the unchanged full78-angle reconstruction route at L800 and avoids a new chunk adapter. It does not raise the whole-job limits. No silent fallback, giant Wigner table, high-L full-volume retained grid or production path change.

## Controls and bounded queue

Reuse the hash-qualified L200/L300 results where identical. The primary new cases are DSpec and YSpec at L400, L600, L800. Keep radial settings and physical fixed probes comparable across degree; independently verify radial accuracy at the highest L with higher-order controls, rather than assuming the L200 radial result transfers. YSpec radial ODE tolerance remains1e-7 unless an explicitly reviewed tighter reference control is needed. YSpec's SEM-like output mesh is not its radial-solver discretization.

Retain the168 fixed probes, layer/material/trace identity and source-shell/contact limitations. Report the89 far-solid-interior probes and depth strata, solid interfaces, source-shell and positive-distance near-field cohorts separately. A few shallow probes dominate the previous aggregate metric: preserve that information rather than allowing a bulk norm to hide it.

The initial physical-job ceiling is12, including radial/diagnostic controls. Exact counts/routes are finalized after the API and resource audit and independently reviewed before results. No mandatory filler jobs. No L>800 or unbounded sweep. A failed control may receive a narrowly scoped reviewed correction with its failed artifacts preserved.

## Pre-result acceptance principles

The target fraction is0.001 (0.1%). This is an operational observed-convergence target for defined quantities, not a rigorous universal continuum error bound.

- Compare measurement-window waveform, measurement-rate normalization and adjoint weights, as well as the whole trace as a diagnostic. Preserve exact window coordinates/tapers and distinguish weighted from unweighted sample metrics.
- Compare own-weight integrated fixed-probe kernels, per-frequency contributions and saved forward/unit-adjoint strain diagnostics. Where supported, hold the reference weights fixed in a separately labelled counterfactual to distinguish weighting from field changes. Inspect degree-block tails before attributing common behavior to a particular cause.
- Compare representative smooth-perturbation derivative actions using only a supported, independently reviewed route. Establish the action/quadrature/finite-difference definition and numerical budget before execution; do not present fixed-point or section sums as a full three-dimensional action. If a high-L action is not feasible in the existing API/resource scope, state the gap before results and do not claim it was validated.
- Require two successive cutoff increments below0.1% for the named acceptance quantities, and require the chosen lower cutoff to agree with the highest available reference within0.1%. Retain comparisons300->400,400->600,600->800 and each candidate against800. Cancellation or one small increment alone is insufficient. If only the last increment passes, report that the ladder lacks the additional level needed to certify a cutoff.
- Declare denominators and zero-reference classes explicitly; save absolute differences and reference norms without arbitrary denominator floors. Separate supported near-zero absolute tolerances from relative targets before results.
- Verify the radial error at the largest L is below the stated budget. Keep solve cutoff distinct from output sampling; adequate sampling or exact reconstruction at fixed probes is needed for each diagnostic. Coarse interpolated wedges do not certify angular convergence or global extrema.
- Do not demand a bounded pointwise value at ideal source/receiver support or identify a singular limit solely from growth. No new H1 formulation, smoothing, physical source regularization, receiver averaging or contact correction is authorized.

If the evidence does not establish a cutoff within the ladder, report which observables/cohorts remain unresolved and their measured changes. Do not label L800 converged simply because it is the last level.

## Execution and provenance

All solver builds/numerical runs are detached on Tunya, in fresh isolated build/input/control/results roots. Eight outer workers, one inner library thread; build caps8threads/32GiB, numerical caps20threads/16GiB, campaign aggregate<=50threads. The accepted user-space identity/ancestry watchdog is not kernel aggregate enforcement. Perform quiet-process/resource/dependency preflight; never interrupt other user work or run numerical jobs concurrently with a build.

Reuse the established SSH transport unchanged, with explicit external-access execution approval; do not edit home/authentication/security configuration. Pin external dependencies, exact commands, binaries, inputs and output/log hashes. Retrieve and verify evidence locally. No deletion, commits, pushes or production edits. Record every change in the appropriate implementation_status.md. Parent owns build orchestration and reports; Luna owns scientific source/checkers; Sol reviews independently.

## User-authorized shared memory amendment

The user explicitly permits revised internal reservations provided combined active usage on Tunya does not exceed900GB and asks this chat to liaise with CCDSM and SpecM3D implementation. GB here means decimal10^9bytes, notGiB. Parent contacted the actual chats CCDSM (01a18251-e70e-73a1-9212-9ab80eb289a1) and DSpecM3D Implementation (01a18879-2a7b-7810-bb21-9d8528e356a2). There is already a shared CCDSM global_ram_guard_v4 with a valid active-RAM metric,900GB ceiling and earlier825/875GB protective thresholds. Preserve that distinction from occupied RAM that includes reclaimable Linux cache/evictable ARC.

Parent requests a conservative40GB active-RAM reservation for this DSpec campaign, subject to agreement with the other owners and verified shared admission. This covers at most one32GiB build or one16GiB numerical job, never both. Existing eight-worker/one-inner-thread controls remain. Do not launch remote builds or solves until the shared reservation/admission mechanism is agreed, current valid guard readings are checked and independent review approves the exact wrapper/guard. Local implementation and review may continue. Do not pause, terminate, thaw, relaunch or otherwise modify another campaign's processes. Parent owns cross-chat resource coordination; implementation and numerical plans must consume the agreed mechanism.

The full78-angle L800/NQ10 conservative bound is27,564,576B angular geometry and44,920,488B combined angular cache+temporary work, plus a separately bounded15,288,000B output-product buffer. Fresh-candidate32/64MiB limits leave margin and are far below the whole-job memory reservation. These changes are explicitly authorized budget changes, not scientific/tolerance changes.
