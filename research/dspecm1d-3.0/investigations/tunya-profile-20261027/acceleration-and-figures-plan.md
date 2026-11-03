# DSpec acceleration and kernel figures

Updated 2026-10-27. Phase A is complete. Alex has approved Phase B on Tunya for source–receiver great-circle sections and three selected native-depth slices only. The full global volume is excluded. DSpec only, Generated Wigner values, detached nohup execution and maximum50aggregate campaign threads. Phase B replay controls completed; production/automatic plotting queue launched under cap40. All existing scientific evidence and historical numeric gate results are preserved.

## Scientific objective and fixed reference

Produce interpretable sensitivity-kernel figures promptly for the existing PREM vertical-displacement energy observable: Hann window 1000–1200 s, nominal 0–10 mHz band with the existing taper, damping, DC treatment, FFT convention, source and receiver. Preserve the completed profile-06 evidence and use its L128/L384/L1536 cases for timing comparisons. Existing native/profiler output parity and fixed-L waveform checks passed; degree convergence and pointwise kernel agreement are not established by that result.

Receiver lifts/enrichment, receiver-point convergence and a Sobolev/Riesz gradient formulation are deferred. Endpoint singularities must not become prerequisites for figures. Retain existing source-lift and one-sided interface semantics during acceleration; no new source/contact or receiver formulation is in scope. Existing contacts-off volume outputs must remain explicitly labelled. Source-shell contact contributions remain separate.

## Phase A: remaining accelerations and bounded verification

Implement and measure these changes in sequence, retaining an individually attributable result for each:

1. **Fixed coefficient layout.** Replace repeated degree/order sets and per-record maps with cached structural permutations, contiguous coefficient arrays and reusable buffers. Preserve channel/material mapping. Primary measured target: 134.82 s of exclusive packing at L1536.
2. **Shared frequency responses.** Reuse matrices, factors and forward/unit-adjoint responses across baseline and kernel stages. Retain or spool compact modal responses, determine the energy-adjoint weights from the waveform, then apply them without a second forward factorisation pass. Record memory/disk cost and actual factorisation/RHS counts. The 115.51 s baseline phase indicates the scale of opportunity, not a promised saving.
3. **Batched native-SEM recovery.** Gather element solution blocks and apply cached derivative matrices in batches; replace repeated dynamic six-component allocations and local-to-global lookups. Preserve fluid handling, both sides of interfaces and existing source lifts. Measured radial-recovery target: 64.55 s at L1536.
4. **Reusable angular plans and global GSH execution.** Reuse geometry, phases and angular bases across radial sides and frequencies. Use Eigen internal multiplication. Choose bounded caches/tiles from measured Tunya memory requirements; retain the existing angular-library degree ceiling. Integrate the previously bounded global-GSH prototype into the physical multi-frequency/native-radial pipeline, with explicit frame/spin/material parity checks. The output angular grid must resolve the complete kernel bandwidth, rather than use arbitrary angle counts.
5. **Frequency accumulation with bounded output traffic.** Accumulate matching-frequency field products in memory over an output tile and write the completed tile once. Store/reuse compact frequency-dependent modal data so tile processing does not repeat radial solves. Compare bounded tiling with a full in-memory accumulator on Tunya. At L1536 six float64 volume channels occupy about 151.47 GiB; avoid the current approximately 15.38 TiB of repeated logical read/write traffic. Independently summing the two fields over frequency before forming their product is not an equivalent algorithm.
6. **Symbolic reuse and controlled parallelism.** Reuse symbolic analysis only for verified-compatible sparsity patterns, retaining the accepted numerical factorisation route. After removing duplicate work, add bounded frequency/tile parallelism with aggregate process/thread accounting, deterministic ownership and no nested thread oversubscription. Compare one-thread improvements separately from parallel speedups. Do not reopen previously rejected band/ordering/condensation algorithms in this phase.

Use the established Luna High implementation, fresh independent Luna High review and short Sol High review sequence; add mathematical review for changes to field representation or contraction order. Record each logical implementation change and its result in implementation_status.md. Preserve failed evidence and use fresh exclusive output paths. No commits, pushes, deletion of existing Tunya files, or unrelated production changes are authorised by this plan update.

### Verification and stopping boundary

- Preserve fixed physical settings, DC/FFT/adjoint normalisation, material channels and native radial-side identities. Require byte parity where arithmetic order is unchanged; otherwise assess against frozen tolerances and saved references, with every difference reported. Do not relax a gate to accept an optimisation.
- Retain the existing qdot/(2E), waveform/energy, finiteness, field/frame and material-mapping checks. Compare a bounded set of useful-amplitude points away from source/receiver neighbourhoods. Check smooth interior behaviour and clearly distinguish inherited discrepancies from new regressions; receiver peaks are not a pointwise acceptance target.
- Perform a bounded resolution comparison of the intended figure sections away from the endpoints. Do not require an expensive new global convergence sweep before exploratory figures, and do not label them converged without evidence.
- After each component, measure affected time, allocation/cache behaviour and memory as appropriate. Run the final complete component profile after the accepted changes; no repeated full campaigns without a new reason.
- Retain a component only if its correctness gates pass and its measured benefit or required global scalability justifies its complexity. Preserve and report unsuccessful candidates; do not turn them into an open-ended prerequisite for plotting.
- All long Tunya jobs must be detached under nohup with persistent logs and restart/provenance receipts. Keep all campaign-owned threads within the maximum of 50. Existing supervision is userspace observation/enforcement, not kernel aggregate accounting.

**Phase A deliverable:** reviewed accelerated path, before/after timing and resource comparison, saved verification evidence, and a concrete figure-generation configuration. Stop for Alex's approval before Phase B.

## Phase B: great-circle sections and selected depth slices

### Execution policy clarified by Alex

- Run DSpecM1D only. Do not launch new YSpec calculations in the next stage; retain existing YSpec results as read-only historical context.
- Use Generated Wigner values. Do not run the full Stored-table mode in this stage; preserve its existing benchmark evidence. The measured small runtime saving does not justify the observed memory footprint for these outputs.
- The unresolved original-method comparison is entirely DSpec: original pointwise/native angular reconstruction versus new global GSH/FFT reconstruction on matched points. The largest channel-scaled differences are approximately3.1–4.3e-13 against the frozen2e-13 limit. This is separate from the older DSpec/YSpec comparison. Internal global parity does not clear that gate; no tolerance is changed here.
- Sections use native pole-frame tensor/Eigen replay; selected depth slices use Generated global-GSH replay. Alex regards the tiny original/global differences as roundoff and has approved proceeding. Keep validation bounded and DSpec-only; no new full cross-code profiling campaign. Existing contacts/interface/degree-limit qualifications remain.
- Keep all long Tunya jobs detached under nohup and the aggregate maximum at50threads. Phase B is approved and running under a40-thread campaign watchdog.

1. Produce the source–receiver great-circle section and selected depth slices, together with the waveform, filter and measurement window. Use the same physical definition and accelerated path, with angular sampling derived from the kernel bandwidth. Early sections can be produced before full-volume output finishes.
2. Mark source/receiver positions and documented endpoint display masks. Preserve unmasked numerical data. Use consistent signed scales for matched views and disclose clipping/normalisation. Keep native one-sided interfaces separate rather than interpolating across them.
3. Label initial figures exploratory, giving degree, radial resolution, band/window, parameterisation and omitted contact terms. Report the bounded resolution comparison; inherited off-endpoint discrepancies remain disclosed and are not explained away by the receiver singularity.
4. Save only the great-circle section and three selected depth slices, their raw numerical samples and provenance. Do not produce a full179-side global angular volume. Reuse saved compact responses so output selection does not repeat radial solves.

**Phase B deliverable:** reproducible section and selected-depth figures, their raw sampled data/configurations, actual measured costs and stated scientific limits. No full global volume. Further receiver/Sobolev work remains deferred to a separately approved phase.

## Phase B approval and scope replacement — 2026-10-27

Alex approved Phase B for great-circle sections and the three selected depth slices only. The complete global volume is explicitly excluded and supersedes the earlier Phase B volume deliverable. DSpec only; Generated Wigner values; no YSpec or Stored runs. Reuse immutable Phase A compact responses where possible, retain DC/guard-bin mixing and native radial/interface identities, and compare replay against saved samples before production figures. The measured original/global 3–4e-13 discrepancy remains recorded under its original gate; Alex considers it floating-point error and it is not a blocker for exploratory figures. All numerical/build/render jobs nohup with a maximum of50aggregate campaign threads and preserved provenance. No new physics/source-contact/receiver formulation or production repository modifications.

## Current scope replacement: great-circle sections only

Alex considers great-circle sections more informative and no longer needs constant-radius slices. Preserve completed sections and all partial radius-slice data; stop the unfinished constant-radius numerical queue and dependent wait/plot queue. No new latitude-tiled batching implementation is approved. Great-circle angular bases are already generated once and reused across radii/frequencies (7.729s setup versus478.051s process wall atL1536), so latitude-basis reuse is not the main section target. Any new section contraction/frequency batching is a separately approved optimisation phase.
