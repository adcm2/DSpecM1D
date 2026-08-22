#include <DSpecM1D/src/FullSpec.h>
#include <DSpecM1D/src/InputParametersNew.h>
#include <DSpecM1D/src/Profiling.h>

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

using SPARSESPEC::detail::profiling::Category;
using SPARSESPEC::detail::profiling::Data;
using SPARSESPEC::detail::profiling::Mode;

SPARSESPEC::SolverBackend parseBackend(const std::string &label) {
  if (label == "eigen")
    return SPARSESPEC::SolverBackend::EigenSparseLU;
  if (label == "lapack")
    return SPARSESPEC::SolverBackend::LapackBandLU;
  throw std::invalid_argument("backend must be 'eigen' or 'lapack'");
}

void writeCounts(std::ostream &out, const Data &data) {
  const auto &c = data.counts;
  out << "\"counts\":{";
  out << "\"sems\":" << c.sems << ",\"degrees\":" << c.degrees
      << ",\"frequency_systems\":" << c.frequencySystems
      << ",\"eigen_compute\":" << c.eigenCompute
      << ",\"eigen_factorize\":" << c.eigenFactorize
      << ",\"lapack_factorize\":" << c.lapackFactorize
      << ",\"solves\":" << c.solves << ",\"band_packs\":"
      << c.bandPacks << ",\"rhs\":" << c.rhs;
  const auto minOrZero = [](long value) { return value == std::numeric_limits<long>::max() ? 0L : value; };
  out << ",\"dimension_min\":" << minOrZero(c.dimensionMin)
      << ",\"dimension_max\":" << c.dimensionMax
      << ",\"nnz_min\":" << minOrZero(c.nonzeroMin)
      << ",\"nnz_max\":" << c.nonzeroMax
      << ",\"kl_min\":" << minOrZero(c.klMin)
      << ",\"kl_max\":" << c.klMax
      << ",\"ku_min\":" << minOrZero(c.kuMin)
      << ",\"ku_max\":" << c.kuMax << "}";
}

void writeData(std::ostream &out, const Data &data) {
  out << std::setprecision(17)
      << "{\"wall_seconds\":"
      << data.seconds[static_cast<std::size_t>(Mode::all)]
                     [static_cast<std::size_t>(Category::total)]
      << ",\"accounting\":{"
      << "\"basis\":\"category timings in OpenMP regions are summed "
         "thread-work seconds; they do not add to wall_seconds\",";
  out << "\"workers\":{";
  for (std::size_t mode = 0; mode < static_cast<std::size_t>(Mode::count); ++mode) {
    if (mode) out << ',';
    out << '\"' << SPARSESPEC::detail::profiling::modeName(static_cast<Mode>(mode))
        << "\":{\"total\":" << data.workerTotal[mode]
        << ",\"categorized\":" << data.workerCategorized[mode]
        << ",\"unclassified\":" << data.workerUnclassified[mode] << '}';
  }
  out << "}},\"timings_seconds\":{";
  for (std::size_t mode = 0; mode < static_cast<std::size_t>(Mode::count); ++mode) {
    if (mode) out << ',';
    out << '\"' << SPARSESPEC::detail::profiling::modeName(static_cast<Mode>(mode))
        << "\":{";
    for (std::size_t category = 0;
         category < static_cast<std::size_t>(Category::count); ++category) {
      if (category) out << ',';
      out << '\"'
          << SPARSESPEC::detail::profiling::categoryName(static_cast<Category>(category))
          << "\":" << data.seconds[mode][category];
    }
    out << '}';
  }
  out << "},";
  writeCounts(out, data);
  out << '}';
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 5) {
    std::cerr << "usage: stage_profile_solver PARAM_FILE BACKEND FIRST_REP REPETITIONS\n";
    return EXIT_FAILURE;
  }
  try {
    const int firstRep = std::stoi(argv[3]);
    const int repetitions = std::stoi(argv[4]);
    if (firstRep < 1 || repetitions < 1)
      throw std::runtime_error("rep numbers must be positive");

    InputParametersNew params(argv[1], 5, 10, 0.05, 1.0, 0.05, 0.0, 1);
    params.setSolverBackend(parseBackend(argv[2]));
    SPARSESPEC::SparseFSpec solver;
    std::ofstream sink("/dev/null");
    auto *saved = std::cout.rdbuf(sink.rdbuf());
    (void)solver.spectra(params); // one untimed warm-up
    for (int repetition = 0; repetition < repetitions; ++repetition) {
      const auto result = solver.spectra(params);
      if (result.size() == 0 || !result.real().array().isFinite().all() ||
          !result.imag().array().isFinite().all())
        throw std::runtime_error("spectra() returned non-finite output");
      const auto &profile = SPARSESPEC::detail::profiling::last();
      std::cout.rdbuf(saved);
      std::cout << "PROFILE\t" << argv[2] << '\t'
                << firstRep + repetition << '\t' << result.rows() << '\t'
                << result.cols() << '\t' << std::setprecision(17)
                << result.cwiseAbs().sum() << '\t';
      writeData(std::cout, profile);
      std::cout << '\n' << std::flush;
      saved = std::cout.rdbuf(sink.rdbuf());
    }
    std::cout.rdbuf(saved);
  } catch (const std::exception &error) {
    std::cerr << "stage_profile_solver: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
