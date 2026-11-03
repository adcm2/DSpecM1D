# Short Sol review: combined-v4b packaging paths

2026-10-27. **PASS for the corrected packaging/build paths and binary-path-only plan updates.** Read-only static/data checks; no SSH, build, or science execution. Source remains the reviewed combined-v4 candidate at the fresh Luna compile-fix hashes.

Reconstructed the exact file-selection and archive-relative path mapping used by `tools/stage_reviewed_build.py` for candidate `combined-v4`, prefix `combined-v4b`, and `--stage-dir combined-v4b`. The helper includes `configs/combined-v4b-build.json` and `configs/combined-v4b-build-plan.json` by those exact names; the final build plan references that included build manifest in both command and declared inputs. The parent corrected the remaining job cwd to the new `acceleration-inputs/combined-v4b` staging root; this final value was re-read and verified.

Each variant cwd/source directory, candidate-first include, and driver argument maps to the included `candidates/combined-v4` tree underneath the new staging root. All referenced candidate driver files are included. Build config changes are only intended staging substitutions and fresh output root `kernel-acceleration-20261027-combined-03`; flags, dependency paths, sequential variants, control root/cap, limits, and reviewed wrapper/tool paths are unchanged. Existing failed build/input evidence is not overwritten. The helper uses exclusive archive/staging-file writes and verifies local/remote manifest hashes and reviewed tool hashes before launch; this review did not execute it.

Full-document normalization proves `reuse-run-v3-plan.json` and `global-run-v2-plan.json` differ from their immediate predecessors solely by candidate build-root 02→03 in commands and inputs. This preserves the earlier reuse-plan PASS. The global-plan finding is a path-delta check only and does not grant an independent global numerical/run approval or supersede any separate HOLD.

Final reviewed SHA-256:

- `combined-v4b-build.json`: `70027c3a5ea5c2ddec99b30f0464dacd16c5b836874a42ff4658c0e651a0a58d`
- `combined-v4b-build-plan.json`: `a22484a4146705b226017c5b8faa44a88d50bbd28d866ffeb0d61b71098eca2f`
- `reuse-run-v3-plan.json`: `ed7e730158284a86f6b65bd3a7e7a9643d11c649e894cfe4f18e1d5e1900b9e4`
- `global-run-v2-plan.json`: `1e242190dbdf35a05f6c210f2afd8bd949db76af03bcbb0f67e576289908b361`

Actual remote freshness, successful compilation, receipted input identity, and numerical acceptance remain execution gates. No scientific result or figure is established here.
