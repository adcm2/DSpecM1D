#include <DSpecM1D/src/FullSpec.h>
#include <DSpecM1D/src/InputParametersNew.h>

#include <chrono>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

struct Result {
  double seconds;
  Eigen::Index rows;
  Eigen::Index cols;
  double norm;
};

std::vector<Result> runMeasured(const std::string &paramPath, int repetitions) {
  InputParametersNew params(paramPath, 5, 10, 0.05, 1.0, 0.05, 0.0, 1);
  SPARSESPEC::SparseFSpec solver;

  // The warm-up and all timed calls use the same complete preferred adaptive
  // multi-SEM spectra(params) workflow.  Setup, model loading, and this
  // warm-up are deliberately outside the measured interval.
  std::ofstream sink("/dev/null");
  auto *saved = std::cout.rdbuf(sink.rdbuf());
  const auto warmup = solver.spectra(params);
  (void)warmup;
  std::vector<Result> results;
  results.reserve(repetitions);
  for (int rep = 0; rep < repetitions; ++rep) {
    const auto start = Clock::now();
    const auto result = solver.spectra(params);
    const auto stop = Clock::now();
    if (result.size() == 0 || !result.real().array().isFinite().all() ||
        !result.imag().array().isFinite().all()) {
      std::cout.rdbuf(saved);
      throw std::runtime_error(
          "spectra() returned a non-finite or empty result");
    }
    results.push_back({std::chrono::duration<double>(stop - start).count(),
                       result.rows(), result.cols(), result.cwiseAbs().sum()});
  }
  std::cout.rdbuf(saved);
  return results;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 5) {
    std::cerr << "usage: stage_i_solver_benchmark PARAM_FILE BACKEND "
                 "FIRST_REP REPETITIONS\n";
    return EXIT_FAILURE;
  }

  try {
    const int firstRep = std::stoi(argv[3]);
    const int repetitions = std::stoi(argv[4]);
    if (firstRep < 1 || repetitions < 1)
      throw std::runtime_error("rep numbers must be positive");
    const auto results = runMeasured(argv[1], repetitions);
    for (std::size_t index = 0; index < results.size(); ++index) {
      const auto &result = results[index];
      std::cout << std::setprecision(17) << "RESULT\t" << argv[2] << "\t"
                << firstRep + static_cast<int>(index) << "\t"
                << result.seconds << "\t" << result.rows << "\t"
                << result.cols << "\t" << result.norm << "\n";
    }
  } catch (const std::exception &error) {
    std::cerr << "stage_i_solver_benchmark: " << error.what() << "\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
