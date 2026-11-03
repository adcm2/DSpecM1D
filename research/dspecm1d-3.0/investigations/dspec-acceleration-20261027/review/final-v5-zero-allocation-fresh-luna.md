# Fresh Luna review: combined final v5 pole-only allocation

**Verdict: PASS.** The v4→v5 code diff is exactly two removed zero-matrix constructions in the spheroidal and toroidal cached-response `poleOnly` branches. The only other difference is the implementation-status entry. No build, SSH, or science run was performed.

In each branch, the removed value was `Zero(response.forward.rows(), 2*degree+1)`. The pole-only cached route instead assigns `response.unitAdjoint` to `adjointUnit`, then places that unit basis in the pole block. `Windowed3DFrequencyWorkspace::prepare` filters nonradial modal blocks from its full-block capture when `poleOnly` is true, and `solveFrequencyModes` calls the full modal visitor only under `!poleOnly && !forwardOnly`. The pole-frame path consumes `adjointUnit`; no downstream pole-only consumer reads the default-empty full-order `adjoint`. This closes the stale-adjoint concern for both families.

The change is confined to cache-hit/miss branches under `poleOnly`. The non-pole branch still derives and assigns `weightedAdjoint`; the non-cache solver branch is untouched. The radial L0 branch and its bypass remain untouched. Because the complete source diff contains only these two removals, cache-off, non-pole, forward-only, and R0 behavior are unchanged by this patch. The parent-supplied profiling estimate motivating the removed zero fills was not independently reproduced in this static review.

## Hash binding

- `candidates/combined-final-v4/engine_headers/Windowed3DFields.h`: `79a79bd14d0a8efe7a120ac9637bd5c3169da603d5b9cccd5e03cae721247dbd`
- `candidates/combined-final-v5/engine_headers/Windowed3DFields.h`: `f8fb772dc97c8d19513bbc0b6ff6bfcaac1b0abb95e771af50b489403429700b`
- `candidates/combined-final-v5/implementation_status.md`: `a46e4c68bb7d7b9412155f8a7bd4b1bbf35f5a30cb20a5bbe2ac604b5e351bc1`

