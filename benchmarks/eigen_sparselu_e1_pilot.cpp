#include "config.h"

#include <DSpecM1D/src/FullSpec.h>
#include <DSpecM1D/src/InputParametersNew.h>

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "usage: eigen_sparselu_e1_pilot PARAM_FILE\n";
    return EXIT_FAILURE;
  }
  try {
    InputParametersNew params(argv[1], 5, 10, 0.05, 1.0, 0.05, 0.0, 1);
    params.setSolverBackend(SPARSESPEC::SolverBackend::EigenSparseLU);
    SPARSESPEC::SparseFSpec solver;
    std::ofstream sink("/dev/null");
    auto *saved = std::cout.rdbuf(sink.rdbuf());
    auto warmup = solver.spectra(params);
    if (warmup.size() == 0 || !warmup.array().isFinite().all())
      throw std::runtime_error("warm-up returned invalid spectra");
    const auto begin = std::chrono::steady_clock::now();
    auto result = solver.spectra(params);
    const auto end = std::chrono::steady_clock::now();
    std::cout.rdbuf(saved);
    if (result.size() == 0 || !result.array().isFinite().all())
      throw std::runtime_error("measured run returned invalid spectra");
    std::cout << std::setprecision(17) << "RESULT\t"
              << std::chrono::duration<double>(end - begin).count() << '\t'
              << result.rows() << '\t' << result.cols() << '\t'
              << result.norm() << '\n';
  } catch (const std::exception &error) {
    std::cerr << "eigen_sparselu_e1_pilot: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
