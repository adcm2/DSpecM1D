#include <filesystem>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <DSpecM1D/src/FullSpec.h>
#include <DSpecM1D/src/InputParametersNew.h>
#include <DSpecM1D/src/OutputWriters.h>
#include <DSpecM1D/src/SignalFiltering.h>

namespace {

std::filesystem::path
resolveOutputPath(const std::string &paramPath, const InputParameters &params) {
  namespace fs = std::filesystem;

  fs::path outputPath(params.output_prefix());
  if (outputPath.is_absolute()) {
    return outputPath;
  }

  const fs::path inputPath(paramPath);
  const fs::path paramDir = inputPath.parent_path();
  if (paramDir.filename() == "params" && paramDir.parent_path().filename() == "data") {
    return paramDir.parent_path().parent_path() / outputPath;
  }

  return paramDir / outputPath;
}

int
run(const std::string &paramPath) {
  InputParametersNew paramsNew(paramPath);

  SPARSESPEC::SparseFSpec solver;
  auto rawSpectra = solver.spectra(paramsNew);
  rawSpectra *= paramsNew.normFactor();

  DSpecM::FilterOptions filterOptions;
  auto filtered = DSpecM::applyFilter(rawSpectra, paramsNew.freqFull(),
                                      filterOptions);

  std::filesystem::path outputPath =
      resolveOutputPath(paramPath, paramsNew.inputParameters());
  outputPath += "_t.out";
  if (!outputPath.parent_path().empty()) {
    std::filesystem::create_directories(outputPath.parent_path());
  }
  DSpecM::writeTimeSeries(outputPath.string(), paramsNew, filtered.timeSeries);

  std::cout << "Wrote filtered time-domain seismogram to "
            << outputPath.string() << '\n';
  return 0;
}

int
runRaw(const std::string &paramPath, bool writeTime, bool writeSpectrum) {
  InputParametersNew paramsNew(paramPath);
  SPARSESPEC::SparseFSpec solver;
  auto rawSpectra = solver.spectra(paramsNew);

  auto outputPath = resolveOutputPath(paramPath, paramsNew.inputParameters());
  std::cout << std::setprecision(17)
            << "raw_time_undamping_epsilon_per_s="
            << paramsNew.freqFull().ep() /
                   paramsNew.earthModel().TimeNorm()
            << '\n';
  if (writeTime) {
    auto rawTimeSpectrum = rawSpectra * paramsNew.normFactor();
    auto rawTime =
        processfunctions::freq2time(rawTimeSpectrum, paramsNew.freqFull());
    auto timePath = outputPath;
    timePath += "_raw_t.out";
    if (!timePath.parent_path().empty())
      std::filesystem::create_directories(timePath.parent_path());
    DSpecM::writeRawTimeSeries(
        timePath.string(),
        paramsNew.inputParameters().time_step_sec(),
        paramsNew.inputParameters().t_out() * 60.0, rawTime);
    std::cout << "Wrote unfiltered raw SI time series to "
              << timePath.string() << '\n'
              << "raw_time_format=time_s; receiver-major columns Z,N,E; "
                 "semicolon-delimited; scientific max_digits10 precision\n";
  }
  if (writeSpectrum) {
    auto frequencyPath = outputPath;
    frequencyPath += "_raw_f.out";
    if (!frequencyPath.parent_path().empty())
      std::filesystem::create_directories(frequencyPath.parent_path());
    DSpecM::writeRawFrequencySeries(
        frequencyPath.string(), paramsNew.freqFull(),
        paramsNew.earthModel().TimeNorm(),
        paramsNew.normFactor() * paramsNew.earthModel().TimeNorm(), rawSpectra);
    std::cout << "Wrote damped raw SI complex spectrum to "
              << frequencyPath.string() << '\n'
              << "raw_spectrum_format=frequency_Hz; receiver-major columns "
                 "Z_Re,Z_Im,N_Re,N_Im,E_Re,E_Im; semicolon-delimited; full one-sided "
                 "DC-to-Nyquist grid; scientific max_digits10 precision\n"
              << "raw_spectrum_scale=normFactor*TimeNorm (observable SI "
                 "units times seconds)\n";
  }
  return 0;
}

int
runExperimentalSH(const std::string &paramPath) {
  // Use the frequency grid of the validated SH experiment.
  InputParametersNew paramsNew(paramPath, 5, 1, 0.01, 1.0, 0.05, 0.0, 4);
  SPARSESPEC::SparseFSpec solver;
  auto result = solver.postUndampingSH(paramsNew);
  if (!result.timeSeries.allFinite())
    throw std::runtime_error("Post-undamping SH produced a non-finite result.");

  const auto &params = paramsNew.inputParameters();
  auto outputPath = resolveOutputPath(paramPath, params);
  outputPath += "_sh_window_t.out";
  if (!outputPath.parent_path().empty())
    std::filesystem::create_directories(outputPath.parent_path());
  DSpecM::writeTimeSeriesOnGrid(outputPath.string(), result.sampleIntervalSeconds,
                                1.0, params.t_out(), result.timeSeries, 24);
  std::cout << "Wrote experimental post-undamping SH seismogram to "
            << outputPath.string() << '\n'
            << "sh_solved_pairs=" << result.solvedPairs
            << " sh_full_grid_pairs=" << result.fullGridPairs
            << " sh_radial_seconds=" << result.radialSeconds
            << " sh_synthesis_seconds=" << result.synthesisSeconds << '\n';
  return 0;
}

double
parseFiniteDouble(const char *value, const char *description) {
  std::size_t consumed = 0;
  const std::string text(value);
  const double parsed = std::stod(text, &consumed);
  if (consumed != text.size() || !std::isfinite(parsed))
    throw std::invalid_argument(std::string("Invalid ") + description + ".");
  return parsed;
}

int
runExperimentalPhaseWindow(const std::string &paramPath, bool transverse,
                           const std::array<double, 4> &corners,
                           std::optional<double> guard) {
  InputParametersNew paramsNew(paramPath, 5, 1, 0.01, 1.0, 0.05, 0.0, 4);
  SPARSESPEC::PostUndampingWindowOptions options;
  options.component = transverse
                          ? SPARSESPEC::PostUndampingWindowComponent::Transverse
                          : SPARSESPEC::PostUndampingWindowComponent::Vertical;
  options.slownessCorners = corners;
  options.preSolveGuardSPerDegree = guard;
  SPARSESPEC::SparseFSpec solver;
  auto result = solver.postUndampingWindow(paramsNew, options);
  if (!result.timeSeries.allFinite())
    throw std::runtime_error("Post-undamping phase window produced a non-finite result.");

  const auto &params = paramsNew.inputParameters();
  auto outputPath = resolveOutputPath(paramPath, params);
  outputPath += transverse ? "_sdiff_window_t.out" : "_pdiff_window_t.out";
  if (!outputPath.parent_path().empty())
    std::filesystem::create_directories(outputPath.parent_path());
  DSpecM::writeTimeSeriesOnGrid(outputPath.string(), result.sampleIntervalSeconds,
                                1.0, params.t_out(), result.timeSeries, 24);
  std::cout << "Wrote experimental post-undamping "
            << (transverse ? "Sdiff T" : "Pdiff Z") << " seismogram to "
            << outputPath.string() << '\n'
            << "window_solved_pairs=" << result.solvedPairs
            << " window_full_grid_pairs=" << result.fullGridPairs
            << " window_radial_seconds=" << result.radialSeconds
            << " window_synthesis_seconds=" << result.synthesisSeconds << '\n';
  return 0;
}

}   // namespace

int
main(int argc, char **argv) {
  bool rawTime = false;
  bool rawSpectrum = false;
  bool rawOptionsOnly = argc >= 3 && argc <= 4;
  if (rawOptionsOnly) {
    for (int i = 2; i < argc; ++i) {
      const std::string_view option(argv[i]);
      if (option == "--raw-time") {
        if (rawTime)
          rawOptionsOnly = false;
        rawTime = true;
      } else if (option == "--raw-spectrum") {
        if (rawSpectrum)
          rawOptionsOnly = false;
        rawSpectrum = true;
      } else {
        rawOptionsOnly = false;
      }
    }
    rawOptionsOnly = rawOptionsOnly && (rawTime || rawSpectrum);
  }
  const bool experimentalSH =
      argc == 3 && std::string_view(argv[2]) == "--experimental-sh-window";
  const bool transverseWindow =
      argc >= 7 && std::string_view(argv[2]) == "--experimental-sdiff-window";
  const bool verticalWindow =
      argc >= 7 && std::string_view(argv[2]) == "--experimental-pdiff-window";
  if (argc != 2 && !rawOptionsOnly && !experimentalSH && !transverseWindow && !verticalWindow) {
    std::cerr << "Usage: generate_seismogram <parameter-file>"
                 " [--raw-time] [--raw-spectrum]"
                 " [--experimental-sh-window | --experimental-sdiff-window"
                 " p0 p1 p2 p3 [--pre-solve-guard-s-per-deg X] |"
                 " --experimental-pdiff-window p0 p1 p2 p3"
                 " [--pre-solve-guard-s-per-deg X]]\n";
    return 1;
  }

  try {
    if (rawOptionsOnly)
      return runRaw(argv[1], rawTime, rawSpectrum);
    if (experimentalSH)
      return runExperimentalSH(argv[1]);
    if (transverseWindow || verticalWindow) {
      if (argc != 7 && argc != 9)
        throw std::invalid_argument(
            "Phase-window options require four corners and an optional guard.");
      std::array<double, 4> corners;
      for (int i = 0; i < 4; ++i)
        corners[i] = parseFiniteDouble(argv[3 + i], "slowness corner");
      std::optional<double> guard;
      if (argc == 9) {
        if (std::string_view(argv[7]) != "--pre-solve-guard-s-per-deg")
          throw std::invalid_argument("Unknown phase-window option.");
        guard = parseFiniteDouble(argv[8], "pre-solve guard");
      }
      return runExperimentalPhaseWindow(argv[1], transverseWindow, corners,
                                        guard);
    }
    return run(argv[1]);
  } catch (const std::exception &ex) {
    std::cerr << "generate_seismogram failed: " << ex.what() << '\n';
    return 2;
  }
}
