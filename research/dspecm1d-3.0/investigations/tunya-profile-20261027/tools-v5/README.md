# v5 qualification snapshot

This snapshot preserves the v4 analytic-reference adapter and its provenance gates. The adapter loads the sibling `qualification_checks.py`; the frozen checker under `tools/` remains unchanged.

The only numerical-reporting change is to the descriptive fluid-kappa comparison shell: it selects regular fluid nodes between 2200 and 3200 km, the outer core, and names the reported metric accordingly. All exact parity gates, per-channel/per-column 2e-13 comparisons, energy checks, full-grid/radial coverage checks, and evidence hashes are unchanged. The cross-code quantity remains descriptive and is not a global scientific-accuracy gate.
