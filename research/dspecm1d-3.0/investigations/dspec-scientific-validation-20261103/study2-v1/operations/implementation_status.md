# Study 2 execution configuration status

2026-11-03: Created `configs/study2-build-v1.json` from accepted Study1 with only fresh Study2 roots/config/task names. Compiler flags and twelve external dependencies unchanged. SHA-256 `15fcf298aba916fda3944c49605f1ca1be4e66dfaf0fbc2d2674474ca6e5a2f9`. Pending complete independent source/build gate.

2026-11-03: Created `configs/study2-build-v1-plan.json` from accepted Study1 with only fresh Study2 roots/config/task names. Compiler flags and twelve external dependencies unchanged. SHA-256 `20c9cfa6a99d92ca9b679d56be66775a939a20126d74e5690c9df109006487a1`. Pending complete independent source/build gate.

2026-11-03: Copied reviewed Study1 minimal build operator to fresh Study2 context, adjusted investigation-root depth, versioned roots/names, and require the changed controls header hash in the source gate. No remote action. Operator SHA-256 `fa84737ae9066e46b4553fcbcd787dfe69b8834f651d307925af90a0d97f51db`.

2026-11-03: Prepared exact reviewed build allowlist, 23 files; archive SHA-256 b9f8358591309489cb719ed5a8b38041a7ca30828f038641e5389029450cd167. No remote action.

2026-11-03: Read-only preflight validated external pins/fresh roots but correctly stopped at the quiet-process guard: another user-owned DSpecPilot3D compiler was active. No Stage2 remote file/build changes and no interruption of other work. Guard context and failure retained; report preparation continues locally before a bounded quiet recheck.

2026-11-03: Build operation `preflight` completed; receipt saved. No solver job launched by this operation.

2026-11-03: Build operation `stage` completed; receipt saved. No solver job launched by this operation.

2026-11-03: Build operation `launch` completed; receipt saved. Detached build launched.

2026-11-03: Retrieved terminal Study2 build evidence and compiler/binary/dependency/cleanup maps. Status completed; evidence SHA-256 `469fb1ad7bea58074638213f0be305db7eab5b715058b219f92e73d13531ba9b`. Independent build-evidence review pending.

2026-11-03: Added root-owned read-only producer retrieval. It binds the exact finalized 9cd09f3e producer state, verifies all 11 job completions and output/log hashes remotely, streams a manifest-bearing archive without remote writes, and extracts only validated ordinary relative files into fresh local evidence. No numerical source, criterion, result or remote state changes. Local Python compilation passed; retrieval pending.

2026-11-03: Read-only producer retrieval completed: 198 files (197 declared outputs/logs plus exact finalized state), all SHA/size checked, archive ac350ba3b51369f78c559c9e4122f239dd0dfde5d97cef8be6846c700f481590 (4,910,938 bytes). Extracted into evidence/producer-raw-v1; receipt retained. No remote filesystem mutation. Independent AFTER2 scientific review can now inspect raw results.
