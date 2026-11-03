# Independent short Sol review: DSpec Wigner v5

27 October 2026. Reviewed after the fresh Luna review in `wigner-v5-luna.md`.
Scope: the frozen v4/v5 source diff, receiver conventions, actual Wigner allocation arguments, and proposed numerical gates. No source edits, builds, or scientific executable runs were made. The parent owns source/binary-manifest verification and transfer provenance.

## Decision: HOLD general acceptance; pole identity is source-consistent

The pole substitution preserves the existing normalization and all spheroidal/toroidal component signs. However, one changed memory-estimate call does **not** match its allocation. This corrects the fresh Luna review's conclusion that every modified site uses capped order. Do not accept the v5 memory-estimator change globally until this is resolved or explicitly excluded from the authorized execution path with evidence. No numerical v5 equivalence or speedup is established by this review.

### P1: full-order sampler now uses capped-order estimate

In `dspec/source-v5-frozen/engine_headers/Windowed3DFields.h:291`, `sampleModalBlocks` uses `angularWignerBytesPerAngle(maxAngularDegree)`, whose implementation at lines 436–441 fixes `mMax=min(2,L)`. The allocation at line 306 still constructs `AngularWigner(L,L,min(2,L),theta)`. Its order bound is **L**, not 2. The v4 estimator matched that full-order allocation; v5 introduces an underestimate for L>2 and can increase chunk size beyond the stated 32 MiB payload target.

The inspected GSHTrans `Wigner.hpp:694-707` allocates one `GSHIndices<All>(L,mMax,n).Size()` block per stored `(n,theta)` plus one `size_t` offset. On the present 64-bit platform the exact payload counts are:

| L | v5 capped estimate per angle | allocated full-order payload per angle | v5 chunk limit | full-order payload at that limit |
|---:|---:|---:|---:|---:|
| 3 | 520 B | 600 B | 64 | 38,400 B |
| 12 | 2,320 B | 6,720 B | 64 | 430,080 B |
| 384 | 76,720 B | 5,928,960 B | 64 | 379,453,440 B |
| 1536 | 307,120 B | 94,494,720 B | 64 | 6,047,662,080 B |

These are integer counts derived from source, not measured RSS or scientific results. For L>=2 the full-order double count is `5*(L+1)^2-10`; the capped-order count is `25*L-15`. Both tables add five offset entries. Scratch, allocator and container overhead are additional.

The principal pole-frame workspace sampler at lines 1940–1948 deliberately passes only radial blocks to `sampleModalBlocks`, so that route avoids its nonradial Wigner allocation. The actual capped-order allocations at lines 1061–1064 and 1656–1658 agree with the new estimator. Thus this finding need not block transfer or a tightly enforced pole-route numerical micro-gate, but the general frozen v5 source cannot be called fully corrected. If a fix is approved, preserve v5 and issue a new frozen revision; do not silently alter frozen sources. Restoring a full-order estimator at line 291 is sufficient for the observed mismatch; reducing the allocation to mMax2 would require proving all sampled modal orders are bounded and is a different change.

## Pole normalization and six receiver channels

The old accessor at lines 454–457 returns the GSHTrans normalized value times `exp(i*m*phi)`. GSHTrans documents and implements `sqrt((2l+1)/(4*pi))*d^l_(n,m)(theta)`; its normalization pass is in `Wigner.hpp:582-593`. Therefore at theta=phi=0 the mathematical value is `rho*delta(n,m)`, with `rho=sqrt((2l+1)/(4*pi))`. No extra `(-1)^m`, conjugation, or interchange of m/n belongs in the pole replacement.

For l>=1, columns ordered m=(-1,0,+1), the new code produces:

```text
S / rho = [  0   1   0 ]      T / rho = [  0   0   0 ]
          [ -1   0   1 ]                [  i   0   i ]
          [ -i   0  -i ]                [ -1   0   1 ]
```

These are exactly the prior component combinations: S radial=y0, S theta=yp-ym, S phi=-i*(ym+yp); T radial=0, T theta=+i*(ym+yp), T phi=yp-ym. All six S/T component channels agree algebraically. At l=0 only m=0 remains and the tangential indicators vanish. Radial-family l=0 is synthesized separately; `receiverAtPole` has an S-versus-else branch rather than a general three-family contract, so a test should not silently treat a radial argument as a supported pole receiver. Its actual callers use nonradial pole blocks.

The new weight/mixing scopes preserve the algebra: the receiver conjugate, `conj(fz)`, adjoint weights, adjoint basis transform, and unit-adjoint projection remain in the same multiplication order. Explicit intermediate matrices are now materialized; compilation and numerical gates, rather than source inspection alone, must establish evaluation parity. The three scope names should be added to run metadata/stage interpretation because the frozen stage map did not acquire them.

## Required numeric and profiling gates

1. **Old versus new pole receiver:** use the exact pinned GSHTrans dependency and separately compiled or otherwise faithfully isolated old/new implementations. Test S and T at L=1,2,3,12,384,1536, plus the l=0 boundary. Compare every complex entry of all three components and every m column, not just norms. Require identical shape and finite values. Against the analytic matrices above require diagonal/present-entry normalization within 8 ulps and phase/sign consistency, and old forbidden entries <=`5e-15*rho`; v5 forbidden entries must be exact zero. Also report maximum old/new scaled error. A failed high-L old-path test is evidence to investigate, not permission to relax the gate or infer v5 end-to-end correctness. Check radial l=0 through its actual separate route.
2. **Full L=3 physical parity:** identical pinned inputs, frequency grid including damped DC, radial one-sided nodes, all complete native-grid angles and layout. Compare v4/v5 spectra, traces, energy, all six kernel families A/C/F/L/N/kappa, and enabled contact/reference/pole controls. Require complete output sizes, matching metadata/layout, and every value finite. For each nonzero array require `max_abs(new-old)/max_abs(old)<=5e-13`; for a zero/near-zero array also require `max_abs(new-old)<=5e-14*max(1,max_abs(old))`. Report relative L2 for each kernel family as supporting evidence. Preserve the existing accepted energy-adjoint identity tolerance and compare `q dot P(B raw)/(2E)` with its expected identity; do not replace this with a forward-waveform-only gate. Keep exact byte equality where obtained, but a source-only assertion is not a gate result.
3. **Native/profile numerical equality:** for each v5 build pair on identical inputs require numerical output bytes to match exactly, excluding timing/instrumentation metadata. Confirm new receiver-scope calls equal the actual receiver-by-pole-block constructor count summed over workspaces. Use exclusive scope closure against the physical root; do not add nested inclusive scopes together. Save mapping/interpretation of the three new names.
4. **Bounded physical probes after the preceding gates:** proceed to the authorized L128/384/1536 probes only after resource preflight and an explicit route check for the estimator finding. Measure wall time, CPU time, peak/observed memory, chunk dimensions and complete versus partial output status. For any speed claim, run serial same-CPU v4/v5 controls, with identical data, compiler settings, numerical thread settings and physical scope; separately measure native/profile overhead. Record CPU identity/affinity and concurrent load and repeat a pair if noise prevents a stable conclusion. Preserve commands, hashes, logs and failed/partial evidence.

The user-authorized Tunya continuation, nohup persistence and campaign-wide 50-thread maximum remain the execution boundaries; they do not substitute for numerical gates. A watchdog is userspace supervision, not kernel aggregate accounting. No source edit, high-degree library repair or broader mathematical claim is authorized by this report. Full global kernel materialization and broad angular convergence remain separate outstanding scientific requirements.
