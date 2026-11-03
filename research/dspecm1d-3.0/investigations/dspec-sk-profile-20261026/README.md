# DSpecM1D SK profiling — 26 October 2026

Investigation and proposed phases only. No production solver change is included.

Read `report/output/pdf/dspec-sk-performance-plan.pdf` for findings and the approval-gated plan. Editable sources are in `report/`.

- `timing_summary.json`: first-call native timings; medians/ranges of three independent processes per configuration.
- `records.json`: full run metadata, exact commands, output hashes and scoped counters. DSpec two-objective timings sum two scalar API calls; YSpec batches them. YSpec scoped counters include two calls (first and warm); normalize by two.
- `analysis.json`: disjoint DSpec attribution and bounded chunk diagnostic. Nested counters in records are not additive.
- `validation.json`: all four existing accuracy gates, plus byte-identical instrumented/native outputs for both solvers in each configuration.
- `provenance.json` and `source_snapshot/`: checkout identities, source/dependency/executable hashes, and accepted dirty-source snapshots.
- `runs/`: raw logs and losslessly compressed actual timed kernel output. Failed one-angle YSpec adapter attempt is retained; its API requires at least three angles.

Reproduction uses the standard source checkouts and the pre-existing pinned dependency cache recorded in the build commands. Build diagnostic drivers with `python3 prepare_probes.py`, `python3 prepare_yspec_probe.py`, and `python3 prepare_chunk_probe.py`. The scripts write only to this evidence directory. `run_profiles.py` resumes recorded configurations without overwriting successful results. To obtain independent fresh timings, copy only these diagnostic scripts into a new empty evidence directory; run the builders there, then run `run_profiles.py`. Run `validate_results.py` afterward. The original native YSpec executable is referenced from its existing cache; the timed YSpec translation unit is recompiled from the live clean checkout.

For an individual recorded command, choose a fresh existing empty output directory before replaying it. Do not reuse populated benchmark directories. All profiling runs set numerical thread environment variables to one, and the compile commands disable GSHTrans BLAS by omission of its enable definition. Reproduction assumes the recorded source/dependency revisions; verify provenance if a checkout has since changed.

The chunk driver takes MODEL SAMPLES OUTPUT, uses one prepared degree-12 frequency, and reports API sampling times and old/new kernel discrepancies. Its output is `chunk_results.json`.

The system LaTeX engine successfully compiled the final four-page PDF with no overflow warnings. Every final page was rendered and visually checked. The desktop compiler lacked its cached bundle, so the system engine supplied the PDF.
