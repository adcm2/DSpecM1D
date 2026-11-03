# Short final Sol review: combined-v4 delta

2026-10-27. **PASS for the bounded source/build delta and updated direct-reuse run plan.** Static review only; no build, scientific execution, SSH, or implementation edits. This carries forward `combined-build-short-sol.md` and `reuse-run-short-sol.md` within their existing boundaries.

Verified the fresh Luna compile-fix hash bindings. Relative to combined-v3, exactly three source headers change (plus the implementation status note); the driver remains `357e0038c5e56188d24b47de99697c79fe4a1b3bd7e3b2b340cfb274b928d5cf`.

- Removing `const` from the spool writer byte accessor permits the existing open-stream `tellp()` query without changing closed-stream file-size behavior or serialization.
- Making only the footprint result/helper public resolves the intended driver preflight access; `private:` resumes immediately before the solver implementation. No footprint or response arithmetic changes.
- Consumer-local angular helpers reproduce the prior private engine formulas for Wigner counts, index bytes, and coefficient-table dimensions. Checked addition/multiplication adds overflow rejection. The existing degree ceiling bounds the remaining small integer expressions. This changes access/accounting implementation, with no rotation, transform, mixing, or kernel-contraction change.

All quoted candidate includes resolve. Independent full-document normalization confirmed build configuration changes are limited to combined-v4 source/staging paths and fresh build root `kernel-acceleration-20261027-combined-02`; compiler flags, dependencies, sequential variants, receipt wrapper, resource limits, and independent control root remain unchanged. The staged manifest name remains `configs/combined-build.json`, matching the build-plan reference; parent staging must preserve that reviewed mapping.

The updated reuse plan is exactly the previously reviewed 25-job plan with candidate binary paths changed from build 01 to build 02 in both command and declared input. Explicitly verified candidate binary at command index 5 and matching input identity. CPU 1, symbolic env 0/1, forward-only off, cache-on/spool-default-off route, original-v6 references, immediate failure gates, and the fixed `2e-13` comparator are unchanged. Run eligibility still requires successful receipted build and parent launch checks; this review is not compilation or numerical evidence.

SHA-256:

- `combined-v4-build.json`: `7f613c638b823defcb28e336fcfbd3cf18df2dd62fa145574645ccd789e0a6e6`
- `combined-v4-build-plan.json`: `1a11789424e4ea99dd076d4c23f342354c8adec5172a7f1ac30c8ad768ea60fd`
- `reuse-run-v2-plan.json`: `46efa86d9be0c7f93f4c01ed15b08e2d373584d3844890de89907bb014fbd347`

Global numerical checker HOLD remains outside scope. No global science run or figure is approved.
