# Sol delta review: build-03 portable DSpec includes

27 October 2026. Read-only local review after parent-reported fresh Luna PASS (formal report pending). No SSH, compilation, scientific execution or source edits performed.

**PASS for the build-only portability delta.** This supplements the earlier build-only acceptance and does not release the separately held numerical checker/plan.

Independently checked both `source-v4-tunya.sha256` and `source-v6-tunya.sha256` from their intended `dspec/` working directory: 12 of 12 entries passed in each manifest. Each portable tree has exactly the same 12-file set as its frozen original; only `engine_headers/SensitivityKernels.h` differs. Each diff consists of exactly five absolute workstation includes replaced by these include-root-relative names:

```text
DSpecM1D/src/SEM/SEM.h
DSpecM1D/src/ParamInfo.h
DSpecM1D/src/SpecHelpers.h
DSpecM1D/src/frequency_info/FreqFull.h
DSpecM1D/src/frequency_info/PostprocessFunctions.h
```

Both `Windowed3DFields.h` files are byte-identical to their respective reviewed frozen originals. No receiver formula, memory estimate, numerical statement, solver setting, profiling scope or other source content changed. The existing `-I/space/adcm2/research/projects/DSpecM1D` compiler argument resolves the five new `<DSpecM1D/src/...>` includes. The original frozen trees and earlier failed build evidence remain separate.

Parsed all three `configs-v3` JSON files and compared them structurally to `configs-v2`. They are exactly equal after only continuation input v2→v3, fresh build root `kernel-profile-20261027-02`→`kernel-profile-20261027-03`, and the DSpec source paths changing to the two portable copies. Source tracking paths follow the same copies. Compiler flags, existing DSpec include root, seven variants, sequential build model, environment and 50-thread supervision remain unchanged. The run plan's same path substitutions do not constitute acceptance of its unfinished checker.

Reviewed config SHA-256 values:

```text
490c5835a187f4bb34f489d74f3344642f7eb86b30956dab2855d5dedae858b4  configs-v3/build.json
809058d61a5b2e494d4bc33bc3515ff628dc73fa44eea16df7071ecfdaa25ba8  configs-v3/build-plan.json
e9a0e8eb947895c435c87697f9bcafca159b789b8b4614286bae4bb8080b9d7a  configs-v3/qualification-and-profile-plan.json
```

Parent-owned preflight must verify exact transferred hashes, stage local `configs-v3` contents at remote `continuation-inputs/v3/configs/`, retain the earlier failed tasks/build roots, and launch a new build-only supervisor task against a fresh build-03 root. The parent reports the actual build-02 compiler failure was the absolute include; this review verifies the narrow correction, not remote build success. Preserve watchdog limitations and the campaign-wide 50-thread maximum. Numerical execution remains HOLD until the final checker completes its required reviews.
