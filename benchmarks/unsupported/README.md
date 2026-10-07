# Unsupported benchmark reports and future QSSP comparison

These historical reports are retained as exploratory research material. They are not supported benchmark targets, automated regression tests, or evidence validating the current DSpecM1D revision. No CMake/CTest or public website integration is provided. The original PDFs are preserved byte-for-byte; hashes and original locations are in `provenance.json`.

## Reports

### [DSpecM1D–QSSP waveform and performance comparison](reports/benchmark_report.pdf)

Five-page report dated 18 August 2026, describing DSpecM1D v1.0.0 (reported commit `9b13dc56`) and QSSP 2024 (reported commit `907706b7`) on earth-tunya. It includes shallow-source body/surface waves, a deep-source PKP-family case, and long-period three-component propagation.

The report explicitly describes a practical consistency comparison, not formal numerical validation. Traces were independently normalized, so absolute source/amplitude normalization was not tested. Source-time functions and Earth-model representation were not exactly identical. Timing rows mix DSpecM1D thread counts with serial QSSP and distinguish QSSP Green-function preparation from synthesis/reuse; they must not be quoted as a controlled general speedup.

### [PKP diffraction exploration](reports/pkp_diffraction_report.pdf)

Three-page report of DSpecM1D analogues of Ivan & He figures 1–3: PKIKP, PKP-Bdiff and PKP-Cdiff. The reported setup uses isotropic PREM, a 300-km explosive source, vertical velocity, and independently normalized filtered traces. It is a qualitative analogue, not an exact reproduction of the paper's ak135/GCMT setup. This report is not itself a QSSP comparison. Its original PKP filename is retained.

The PDF records DSpecM1D v1.0.0 / `9b13dc5`, attenuation off, self-gravity on and rotation off. These are claims recorded in the historical report, not independently rerun or verified during preservation.

## Goal: reproducible QSSP benchmarks

1. Recover the original input files, model conversions, scripts, raw waveforms and execution logs if available. Only the two PDFs are preserved here; these are not yet complete reproducibility packages.
2. Pin exact DSpecM1D and QSSP revisions, dependencies, compilers and run settings. Start with the existing low-frequency comparison cases before extending to the more expensive PKP-diffraction example.
3. Define matched physical and numerical settings: model/layers, source depth and moment-tensor convention, source-time function, receiver components and signs, units, gravity, attenuation/dispersion, frequency band, sampling and filtering. Document unavoidable differences.
4. Add quantitative comparisons of unnormalized amplitudes, phase/arrival timing and waveform errors, alongside qualitative plots. Use convergence studies and investigate discrepancies before choosing acceptance thresholds; neither implementation is automatically ground truth.
5. Measure performance under declared resources, including matched-thread comparisons where possible. Report cold preparation, synthesis/reuse, wall time and memory separately.
6. Provide an opt-in runner with isolated output directories and explicit resource limits. Promote selected cases into supported benchmarks only after review of provenance, reproducibility and numerical evidence.

This plan does not authorize or imply that any QSSP runs, expensive diffraction simulations, or integration work have been performed.
