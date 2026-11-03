# Scientific report preparation

2026-11-03: Copied the approved standard LaTeX style byte-for-byte for the two scientific study reports. Style SHA-256: de076063c9b60fa57445b12848fdcff5803837e87e22b78255c647c504a90fcc. Numerical results and report conclusions remain pending actual study execution and independent review. Report compilation/rendering will be local artifact preparation; all physical solver builds/runs remain on Tunya. No scientific source changed.

2026-11-03: Added reproducible study-1 report preparation from the frozen actual raw dataset and independently reviewed checker report. It produces a source-trace figure, hash-bound report data and LaTeX source; no physical solves. Script SHA-256 `b90593bfb8b4f4a2f640565636873964756f22106b9ac230863039bc56e95050`. PDF compilation/visual review and report-claim review remain pending.

2026-11-03: Report-only formatting correction: generated tables use proper LaTeX powers of ten for tiny kernel values, epsilon and FD errors. No numerical input/value changed. Script SHA-256 `aec0b05c8c3a476550faa4656b693b1ed8434def1bf7bcdc8e606cac584dcdee`.

2026-11-03: Generated study1.tex, source-trace figure and bound report-data JSON from actual reviewed evidence; compiled LaTeX twice successfully. Current PDF SHA-256 `daafaa79983ec3ebe7022af9ddfae0c6ddb1c0cbc53affeea721d5d6cbefcd65`. Visual layout and final report-claim review pending.

2026-11-03: Study1 report finalized after independent scientific claim review and visual inspection of all three rendered PDF pages. No LaTeX warnings or clipping. Final `output/pdf/dspec-study1-traces-20261103.pdf` SHA-256 `daafaa79983ec3ebe7022af9ddfae0c6ddb1c0cbc53affeea721d5d6cbefcd65`; editable source, figure and source-data hashes retained in study1-delivery-manifest.json. No numerical result changed.

2026-11-03: Added reproducible Study2 report generator bound to the actual checker report aedda47c and sealed raw producer state 9cd09f3e. It rehashes all 198 producer files, uses Sol's independent raw metrics for discrete-vector tables, and independently computes the nonzero-window velocity change. Explicitly distinguishes radial stability, fixed-cutoff code agreement and unresolved degree convergence; no full H1 derivation or singular-limit proof. Report generation, compilation and review pending. No numerical source/result changed.

2026-11-03: Generated Study2 LaTeX, degree-ladder vector figure and bound summary JSON; compiled twice to a three-page PDF. The window-only velocity diagnostic uses 99 positive-window samples and gives 0.8725001394% for L200 to L300. All 198 producer files rehashed. Independent AFTER2/report-claim reviews and visual page QA are pending; no solver source/result changed.

2026-11-03: Report2 formatting-only correction removed an extra blank table row before each bottom rule. Values, source data, scope and scientific claims unchanged. Recompilation/rendering follows.

2026-11-03: Independent report review correction: distinguished DSpec eight frequency workers from YSpec eight outer workers, each with one inner library thread; clarified that changing NQ changes SEM polynomial and quadrature order. No data or numerical conclusions changed. Recompilation and final claim seal follow.

2026-11-03: Both authorized scientific studies complete, with sequential independent Sol AFTER1/AFTER2 and report-claim reviews. Final Study2 PDF aa3b633f6435e7807df005206defc7332042c0edbd5aca48b576c47544c29a01, 3 pages, all pages rendered/visually verified with no LaTeX diagnostics; full delivery hashes in reports/study2-delivery-manifest.json. Radial controls and fixed-cutoff cross-code agreement are small; angular L200-to-L300 stabilization is NOT met (whole velocity 4.162965%, far-interior kernel vector 165.89856%). Point-source/receiver singular-limit questions remain unresolved; no H1 derivation/formulation change implemented. All work remains isolated and stopped before further numerics, optimization or production integration.
