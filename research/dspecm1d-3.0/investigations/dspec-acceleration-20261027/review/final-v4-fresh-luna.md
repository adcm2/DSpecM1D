# Combined-final-v4 driver fix: fresh Luna

**Result: PASS for the bounded source-static fix.** The v3→v4 source diff is confined to physical_driver.cpp; all engine and global source files are byte-identical to the reviewed v3 candidate. No compile or science run was performed.

The result-JSON stream expression now uses valid C++ string literals around the "true"/"false" values. Its bypass metadata is based on the completed responseReuseEnabled state, not on the optional math-profile counter, so native profiling being disabled no longer makes a successful reuse run report zero. For kernel mode, DC is required; the baseline covers bins [0,endBin) and the kernel pass covers [0,last], so the reported radial-L0 bypass count endBin + (last + 1) is the expected total when reuse completes, and is zero when it does not. In the fixed case this is 53 + 52 = 105. The independently recorded profile counter remains available for checker comparison.

## Hash bindings

- combined-final-v4/physical_driver.cpp: 8a6b764e47b2f86389969239e058617fc6bbb7947a34fd0770387ba50edc8267
- combined-final-v4/implementation_status.md: c9bd003dde311fc00dc5855d11275b0aec0536ad4cd2dc5c4534f5604b6bf8f2
- Unchanged combined-final-v4/engine_headers/Windowed3DFields.h: 79a79bd14d0a8efe7a120ac9637bd5c3169da603d5b9cccd5e03cae721247dbd
- Unchanged combined-final-v4/global/global_kernel_consumer.hpp: cfe88cf6b02800152bd3fe6f500e718392b66cc4e9533ee27c4157b2027741f6
- Unchanged combined-final-v4/global/pole_response_spool.hpp: 5e8868439771cb502cb4e3005039f6986e9d0d6ad23a58b168a7a61c2c9152dc

## Gate boundary

This PASS covers the C++ literal repair, source-diff scope, and bypass-metadata derivation. It does not establish compilation or runtime output. No build, SSH, or science run was performed.
