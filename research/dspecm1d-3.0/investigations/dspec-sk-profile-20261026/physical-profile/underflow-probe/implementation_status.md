# Underflow probe implementation status

- 2026-10-26: Created this isolated diagnostic scope only. The preload is environment-gated (`DSPEC_PROFILE_FTZ=1`), enables x86 MXCSR FTZ+DAZ and logs the before/after control word; otherwise it returns without changing process state. Unsupported architectures/CPUs fail explicitly if enabled.
- Parent owns all runtime tests and interpretation; no executable or preload was run here.
- 2026-10-26: Frozen source compiled successfully through the shared `run_compile.py` gate with `-Wall -Wextra -Werror` (exit 0). Source and binary hashes are recorded in the v1 manifests and exact build vector in `freeze-v1-build.json`. No preload or scientific executable was run by this agent.
