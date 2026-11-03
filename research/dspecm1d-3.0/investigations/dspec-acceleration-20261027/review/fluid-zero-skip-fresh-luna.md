# Fresh review: fluid structural-zero transform skip

**Decision: PASS for this bounded source delta.** Compared combined-final-v2 with reviewed combined-final-v1. No build or science run was performed; this does not validate measured numerical equivalence or performance.

Reviewed SHA-256:

| File | SHA-256 |
|---|---|
| `candidates/combined-final-v2/global/global_gsh_tile.hpp` | `0a884177fc688b6a78d95c04209571b90329124f7d7437af5a063887484671e0` |
| `candidates/combined-final-v2/global/global_kernel_consumer.hpp` | `59315ef04718c1bd8d793b6909ba78f9e9142c7e78edd4604c6d2ee30c8aa4d3` |
| `candidates/combined-final-v2/engine_headers/Windowed3DKernels.h` | `a852dead24ba03c1c4139ef926b807a88551ce2c2d77f443ddf149cf150d5506` |
| `candidates/combined-final-v2/engine_headers/Windowed3DFields.h` | `fcf8c2ec7d428bd3593b0030e75cfcf524382e2d4d5f6496f9aa01824149cb35` |

The candidate diff is limited to `global_gsh_tile.hpp`, `global_kernel_consumer.hpp`, and `implementation_status.md`.

The fluid coefficient helper already defines components 1–5 as exact complex zero and component 0 as the `rr + horizontalTrace` combination. V2 now omits coefficient-array construction/rotation for fluid components 1–5, explicitly fills those ten forward/adjoint field slots with exact zeros, and retains component 0 for both forward and adjoint transforms. The unchanged L0 additions to slots 0, 1, 6, and 7 remain after synthesis. The unchanged fluid `pointKernel` uses only `rr + horizontalTrace`, so slots 2–5 are structurally unused there as well. This is an exact-zero structural skip; no tolerance or threshold was added.

The consumer adds total and per-worker counters for transformed channels and fluid exact-zero channels skipped. They are accumulated from worker-owned plan timing state and written to the manifest. For a fluid frequency-side evaluation the expected counts are 2 channel fields transformed and 10 skipped; for a solid evaluation they are 12 transformed and none skipped. The existing inverse-GSH stage timing includes only calls actually made. The point-kernel loop and all numerical arithmetic on retained fields are unchanged.

The independent radial-cache gate failure and its planned final-v3 response fix remain outside this delta and still require their own review. No code was changed by this review.
