# Luna review: portable DSpec copies and configs v3

**Verdict: PASS for the reviewed portability and config delta.** This records the earlier static review; no implementation was changed and no new build, scientific run, or remote operation was performed for this note.

- The v4 and v6 Tunya source copies differ from their respective frozen trees only in five absolute include directives in `engine_headers/SensitivityKernels.h`, changed to repository-relative `#include <DSpecM1D/src/...>` paths. No math-body changes were present in that portability delta.
- Both portable source trees’ 12-entry SHA-256 manifests verified.
- The v3 `build-plan.json`, `build.json`, and `qualification-and-profile-plan.json` match v2 apart from the intended path substitutions: input root v2 to v3, build root `...-02` to `...-03`, and DSpec source paths to `source-v4-tunya` and `source-v6-tunya`.
- All four DSpec compiler commands preserve `-I/space/adcm2/research/projects/DSpecM1D`, so the new `<DSpecM1D/src/...>` include paths resolve against the reviewed shared DSpec repository root.

This PASS covers only those source-copy and configuration deltas; it is not a build or scientific validation result.
