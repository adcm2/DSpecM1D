# DSpecM1D Test Inventory

This file groups the current tests by library area so it is easy to see which part of `DSpecM1D` each test is protecting.

## Build And Smoke Checks

- `smoke_generate_seismogram`: Runs the public `generate_seismogram` CLI on a small bundled parameter file and checks that a non-empty time-series output file is created.
- `smoke_build_benchmark_ex1`: Verifies that the protected `ex1` paper example still builds in the configured tree.
- `smoke_tutorial_t1`: Runs the self-contained `t1` tutorial/integration workflow to confirm the main synthesis pipeline still executes end to end.

## Input Parsing

- `InputParserTests.GetNextValueLineSkipsCommentsAndQuotes`: Checks that low-level value-line parsing skips blank/comment lines and removes surrounding quotes.
- `InputParserTests.ReadRequiredScalarRejectsTrailingTokens`: Checks that scalar parsing rejects lines with extra trailing tokens instead of accepting malformed input.
- `InputParserTests.ReadRequiredLatLonParsesPair`: Checks that a latitude/longitude pair is parsed correctly from a single value line.
- `InputParserTests.InputParametersParsesKnownGoodFile`: Checks that the legacy ordered parameter-file parser reads a valid synthetic parameter file into the expected fields.
- `InputParserTests.InputParametersRejectsOutOfRangeLatitude`: Checks that the legacy parameter parser rejects a receiver latitude outside the valid `[-90, 90]` range.

## InputParametersNew Workflow Context

- `InputParametersNewTests.ResolvesAbsoluteModelPathAndBuildsContext`: Checks that `InputParametersNew` resolves the Earth-model path correctly and builds the expected bundled workflow objects.
- `InputParametersNewTests.NormFactorMatchesOutputTypeSelection`: Checks that `InputParametersNew::normFactor()` matches the expected SI scaling for displacement, velocity, and acceleration output modes.
- `InputParametersNewTests.SetterClampsToMinimumPositiveIntegers`: Checks that `InputParametersNew` clamps `nq` and `nskip` to valid positive minimum values.

## Reference Data Loading

- `ReferenceSeriesIOTests.LoadYSpecTimeSeriesPadsAndTruncates`: Checks that YSpec time-series loading fills the requested output shape and pads missing samples with zeros.
- `ReferenceSeriesIOTests.LoadSpecnmTimeSeriesParsesSemicolonFormat`: Checks that SpecNM semicolon-delimited time-series files are parsed into the expected three-component layout.
- `ReferenceSeriesIOTests.LoadSpecnmTimeSeriesPadsShortInputs`: Checks that short SpecNM inputs are padded with zeros when more columns are requested than the file provides.
- `ReferenceSeriesIOTests.LoadMineosTimeSeriesAppliesScale`: Checks that Mineos component files are combined into one matrix and scaled by the requested amplitude factor.
- `ReferenceSeriesIOTests.LoadMineosTimeSeriesTruncatesToShortestComponent`: Checks that Mineos loading truncates to the shortest available component before padding the rest with zeros.
- `ReferenceSeriesIOTests.LoadReferenceTimeSeriesBundlesYspecAndMineos`: Checks that the bundled reference loader returns both YSpec and Mineos matrices with the expected contents.
- `ReferenceSeriesIOTests.LoadersReturnEmptyMatricesForZeroRequestedColumns`: Checks that the reference-data loaders return correctly shaped empty matrices when zero output columns are requested.

## Spectral Helper Utilities

- `SpecHelpersTests.ResolveModeFlagsHonorsModeTypeAndAngularRange`: Checks that the mode-selection helper enables and suppresses radial, toroidal, and spheroidal families consistently with the mode flag and angular-degree limits.
- `SpecHelpersTests.OutputFactorMatchesRequestedQuantity`: Checks that the output scaling helper returns the correct displacement, velocity, and acceleration multipliers.
- `SpecHelpersTests.SpecConstantsUsesReferencePeriod`: Checks that the precomputed spectral constants use the supplied reference period and damping parameter consistently.
- `SpecHelpersTests.FactorizeOrComputeChoosesExpectedPath`: Checks that the solver helper calls `compute()` on the configured cadence and `factorize()` on the intervening steps.

## Preferred Solver API

- `PreferredSolverApiTests.SpectraRunContextExposesWorkflowObjects`: Checks that `SpectraRunContext` preserves references to the prepared frequency, source, and parameter objects and clamps `nskip` to a valid value.
- `PreferredSolverApiTests.SemConstructorFromInputParametersNewUsesSettings`: Checks that `SEM(const InputParametersNew&)` reuses the `InputParametersNew` settings consistently when building the mesh.
- `PreferredSolverApiTests.PreferredSolverOverloadsReturnStableShapes`: Checks that the three preferred `SparseFSpec::spectra(...)` overloads produce finite outputs with the expected shape and that the reused-SEM overloads agree exactly.
- `PreferredSolverApiTests.LegacyMultiSemOverloadReturnsFiniteOutput`: Checks that the legacy multi-SEM `SparseFSpec` overload still produces finite output with the expected shape on a tiny setup.
- `PreferredSolverApiTests.LegacySingleSemOverloadReturnsFiniteOutput`: Checks that the legacy single-SEM `SparseFSpec` overload still produces finite output with the expected shape on a tiny setup.

## Output Writing

- `OutputWriterTests.FrequencyWriterTruncatesToSmallestInput`: Checks that the two-way frequency comparison writer truncates to the shortest available input length and writes the expected frequency axis.
- `OutputWriterTests.ThreeWayFrequencyWriterIncludesAllReferenceBlocks`: Checks that the three-way frequency writer emits all expected comparison blocks in one output row.
- `OutputWriterTests.FrequencyWriterThrowsWhenParentDirectoryIsMissing`: Checks that the frequency writer fails loudly instead of silently creating output in a missing parent directory.
- `OutputWriterTests.FrequencyWriterConvenienceOverloadUsesInputParametersNew`: Checks that the `InputParametersNew` two-way frequency writer convenience overload emits the expected output layout.
- `OutputWriterTests.ThreeWayFrequencyConvenienceOverloadUsesInputParametersNew`: Checks that the `InputParametersNew` three-way frequency writer convenience overload emits the expected output layout.
- `OutputWriterTests.TimeWriterStopsAtRequestedOutputWindow`: Checks that the time-domain comparison writer stops at the requested output window rather than writing the full available series.
- `OutputWriterTests.TimeWriterTruncatesToShortestInputAndKeepsLayout`: Checks that the time writer truncates mismatched inputs to the shortest series while preserving the expected column layout.
- `OutputWriterTests.TimeWriterThrowsWhenParentDirectoryIsMissing`: Checks that the time writer reports an error when its parent output directory does not exist.
- `OutputWriterTests.TimeWriterConvenienceOverloadUsesInputParametersNew`: Checks that the `InputParametersNew` time writer convenience overload emits the expected output layout.
- `OutputWriterTests.TimeSeriesWriterWritesAllRows`: Checks that the standalone time-series writer emits one time column plus all signal rows, making it suitable for the public CLI output.

## Filtering

- `SignalFilteringTests.RejectsInvalidPassCount`: Checks that the filter wrapper rejects a pass count less than one.
- `SignalFilteringTests.FilterOptionsExposeExpectedDefaults`: Checks that `FilterOptions` exposes the expected default taper widths, pass count, and real-signal flag.
- `SignalFilteringTests.ZeroFrequencyInputStaysZeroAcrossMultiplePasses`: Checks that an all-zero valid frequency input remains zero after repeated filtering passes.
- `SignalFilteringTests.RejectsFrequencyInputWithUnexpectedColumnCount`: Checks that the frequency-domain filter overload rejects input matrices that do not match the active `FreqFull` width.
- `SignalFilteringTests.RejectsTimeInputWithUnexpectedColumnCount`: Checks that the time-domain filter overload rejects input matrices that do not match the active `FreqFull` time length.
- `SignalFilteringTests.TimeInputPreservesDimensionsAndProducesFiniteOutput`: Checks that filtering a valid synthetic time series preserves the expected output dimensions and produces finite values.
- `SignalFilteringTests.FrequencyInputPreservesDimensionsAndSupportsEnforceRealSignal`: Checks that filtering a valid synthetic frequency series preserves the expected output dimensions and produces finite values when real-signal enforcement is enabled.

## SEM Low-Level Contracts

- `SEMComponentTests.LocalToGlobalMapsAreMonotonicOnMinimalModel`: Checks basic monotonic ordering in the spheroidal and radial local-to-global maps on a minimal SEM setup.
- `SEMComponentTests.ReceiverAndSourceElementsStayWithinMeshBounds`: Checks that receiver and source lookup return element indices that lie within the SEM mesh bounds.
- `SEMComponentTests.SystemMatricesExposeConsistentShapes`: Checks that the main SEM system matrices are square, shape-compatible within each mode family, and structurally non-empty.
- `SEMComponentTests.ReceiverVectorsHaveExpectedDimensions`: Checks that the main SEM receiver-vector helpers return non-zero outputs with dimensions consistent with the active local-to-global indexing.
