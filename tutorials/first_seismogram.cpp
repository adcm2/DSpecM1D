#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <DSpecM1D/src/FullSpec.h>
#include <DSpecM1D/src/InputParametersNew.h>
#include <DSpecM1D/src/OutputWriters.h>
#include <DSpecM1D/src/SignalFiltering.h>

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "Usage: first_seismogram <parameter-file> <output-file>\n";
    return 1;
  }
  try {
    InputParametersNew params(argv[1]);
    SPARSESPEC::SparseFSpec solver;
    auto spectra = solver.spectra(params);
    spectra *= params.normFactor();
    DSpecM::FilterOptions filterOptions;
    const auto filtered =
        DSpecM::applyFilter(spectra, params.freqFull(), filterOptions);
    const std::filesystem::path output(argv[2]);
    if (!output.parent_path().empty())
      std::filesystem::create_directories(output.parent_path());
    DSpecM::writeTimeSeries(output.string(), params, filtered.timeSeries);
    std::cout << "Wrote time-domain seismogram to " << output.string() << '\n';
  } catch (const std::exception &error) {
    std::cerr << "first_seismogram failed: " << error.what() << '\n';
    return 2;
  }
}
