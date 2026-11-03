# YSpec coarse-timer control

This isolated tree copies the active YSpec v2 profiling sources and frontend. Its local `profile.hpp` keeps `Root`, `Scope`, and `count` unchanged while making `stamp`, `elapsed`, and `addChild` no-ops. The original fused angular block therefore retains its whole-block Scope timer, while the fine per-phase clock reads and child commits are suppressed.

`prepare_commands.py` writes the exact compiler/linker command arrays to `build-commands.json`; it does not compile or execute them. Both actions are held until the parent closes the currently running numerical campaign. The existing YSpec v2 binaries and frozen files are untouched.

`angular_phase_timer_reads` in any future coarse output is the source-side count of attempted fine-timer call sites, not hardware clock reads in this no-op variant.
