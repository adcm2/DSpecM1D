#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
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

}   // namespace

int
main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "Usage: generate_seismogram <parameter-file>\n";
    return 1;
  }

  try {
    return run(argv[1]);
  } catch (const std::exception &ex) {
    std::cerr << "generate_seismogram failed: " << ex.what() << '\n';
    return 2;
  }
}
