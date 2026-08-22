#include <DSpecM1D/src/FullSpec.h>
#include <DSpecM1D/src/InputParametersNew.h>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

SPARSESPEC::SolverBackend parseBackend(const std::string &label) {
  if (label == "eigen")
    return SPARSESPEC::SolverBackend::EigenSparseLU;
  if (label == "lapack")
    return SPARSESPEC::SolverBackend::LapackBandLU;
  throw std::invalid_argument("backend must be 'eigen' or 'lapack'");
}

void validateResult(const Eigen::MatrixXcd &result) {
  if (result.size() == 0 || !result.real().array().isFinite().all() ||
      !result.imag().array().isFinite().all())
    throw std::runtime_error("spectra() returned an empty or non-finite result");
}

struct Measurement {
  double seconds;
  Eigen::Index rows;
  Eigen::Index cols;
  double norm;
};

struct BackendRun {
  std::vector<Measurement> measurements;
  Eigen::MatrixXcd finalResult;
};

BackendRun runBackend(InputParametersNew &params, const std::string &backend,
                      int repetitions, SPARSESPEC::SparseFSpec &solver) {
  params.setSolverBackend(parseBackend(backend));

  const auto warmup = solver.spectra(params);
  validateResult(warmup);

  BackendRun run;
  run.measurements.reserve(static_cast<std::size_t>(repetitions));
  for (int repetition = 0; repetition < repetitions; ++repetition) {
    const auto start = Clock::now();
    auto result = solver.spectra(params);
    const auto stop = Clock::now();
    validateResult(result);
    run.measurements.push_back(
        {std::chrono::duration<double>(stop - start).count(), result.rows(),
         result.cols(), result.norm()});
    run.finalResult = std::move(result);
  }
  return run;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 5) {
    std::cerr << "usage: final_solver_benchmark PARAM_FILE FIRST_BACKEND "
                 "REPETITIONS TOLERANCE\n";
    return EXIT_FAILURE;
  }

  try {
#ifndef DSPECM1D_ENABLE_LAPACK_BAND_SOLVER
    throw std::runtime_error(
        "final benchmark requires DSPECM1D_ENABLE_LAPACK_BAND_SOLVER=ON");
#endif
    const std::string firstBackend = argv[2];
    if (firstBackend != "eigen" && firstBackend != "lapack")
      throw std::invalid_argument("FIRST_BACKEND must be 'eigen' or 'lapack'");
    const int repetitions = std::stoi(argv[3]);
    const double tolerance = std::stod(argv[4]);
    if (repetitions < 1 || !std::isfinite(tolerance) || tolerance <= 0.0)
      throw std::invalid_argument(
          "REPETITIONS and TOLERANCE must be positive");

    InputParametersNew params(argv[1], 5, 10, 0.05, 1.0, 0.05, 0.0, 1);
    SPARSESPEC::SparseFSpec solver;
    std::ofstream sink("/dev/null");
    auto *saved = std::cout.rdbuf(sink.rdbuf());

    const std::string secondBackend =
        firstBackend == "eigen" ? "lapack" : "eigen";
    BackendRun first;
    BackendRun second;
    try {
      first = runBackend(params, firstBackend, repetitions, solver);
      second = runBackend(params, secondBackend, repetitions, solver);
    } catch (...) {
      std::cout.rdbuf(saved);
      throw;
    }
    std::cout.rdbuf(saved);

    const BackendRun &eigen = firstBackend == "eigen" ? first : second;
    const BackendRun &lapack = firstBackend == "lapack" ? first : second;
    if (eigen.finalResult.rows() != lapack.finalResult.rows() ||
        eigen.finalResult.cols() != lapack.finalResult.cols())
      throw std::runtime_error("Eigen/LAPACK output shapes differ");

    const double maxDifference =
        (eigen.finalResult - lapack.finalResult).cwiseAbs().maxCoeff();
    if (!std::isfinite(maxDifference) || maxDifference > tolerance)
      throw std::runtime_error("Eigen/LAPACK discrepancy exceeds tolerance");

    for (const auto &backend : {firstBackend, secondBackend}) {
      const BackendRun &run = backend == firstBackend ? first : second;
      for (std::size_t index = 0; index < run.measurements.size(); ++index) {
        const auto &measurement = run.measurements[index];
        std::cout << std::setprecision(17) << "RESULT\t" << backend << "\t"
                  << index + 1 << "\t" << measurement.seconds << "\t"
                  << measurement.rows << "\t" << measurement.cols << "\t"
                  << measurement.norm << '\n';
      }
    }
    std::cout << std::setprecision(17) << "CHECK\t" << eigen.finalResult.rows()
              << "\t" << eigen.finalResult.cols() << "\t"
              << eigen.finalResult.norm() << "\t" << lapack.finalResult.norm()
              << "\t" << maxDifference << "\t" << tolerance << '\n';
    std::cout.flush();
  } catch (const std::exception &error) {
    std::cerr << "final_solver_benchmark: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
