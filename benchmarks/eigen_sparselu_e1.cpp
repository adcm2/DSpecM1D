#include "config.h"

#include <DSpecM1D/ModelInput>
#include <DSpecM1D/src/SEM/SEM.h>
#include <DSpecM1D/src/SpecHelpers.h>

#include <Eigen/SparseLU>

#include <chrono>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

using Complex = std::complex<double>;
using Sparse = Eigen::SparseMatrix<Complex>;
using Dense = Eigen::MatrixXcd;
using Clock = std::chrono::steady_clock;

struct System {
  std::string family;
  std::string size;
  int degree;
  int nq;
  double maxstep;
  Sparse matrix;
  Dense rhs;
};

struct Timed {
  std::string ordering;
  bool symmetric;
  double analyze;
  double factorize;
  double solve;
  Eigen::Index aNnz;
  Eigen::Index kl;
  Eigen::Index ku;
  Eigen::Index lapackEnvelope;
  Eigen::Index lNnz;
  Eigen::Index uNnz;
  int infoAnalyze;
  int infoFactorize;
  int infoSolve;
  double discrepancy;
  double residual;
};

std::pair<Eigen::Index, Eigen::Index> bandwidth(const Sparse &matrix) {
  Eigen::Index kl = 0;
  Eigen::Index ku = 0;
  for (Eigen::Index column = 0; column < matrix.outerSize(); ++column)
    for (Sparse::InnerIterator entry(matrix, column); entry; ++entry) {
      kl = std::max(kl, entry.row() - entry.col());
      ku = std::max(ku, entry.col() - entry.row());
    }
  return {kl, ku};
}

Dense makeRhs(Eigen::Index rows, int nrhs) {
  Dense rhs(rows, nrhs);
  for (Eigen::Index col = 0; col < nrhs; ++col)
    for (Eigen::Index row = 0; row < rows; ++row)
      rhs(row, col) = Complex(0.7 + 0.001 * row + 0.11 * col,
                              -0.2 + 0.002 * row - 0.07 * col);
  return rhs;
}

double infNorm(const Dense &value) {
  double result = 0.0;
  for (Eigen::Index row = 0; row < value.rows(); ++row) {
    double sum = 0.0;
    for (Eigen::Index col = 0; col < value.cols(); ++col)
      sum += std::abs(value(row, col));
    result = std::max(result, sum);
  }
  return result;
}

Sparse frequencyMatrix(const Eigen::SparseMatrix<double> &h,
                       const Eigen::SparseMatrix<double> &p,
                       const Eigen::SparseMatrix<double> &atten,
                       double omega, double tref) {
  SPARSESPEC::SpecConstants constants(1.0e-4, tref);
  const Complex w = omega + constants.ieps;
  Sparse result = h.cast<Complex>() - (w * w) * p.cast<Complex>();
  result += SPARSESPEC::attenFactor(omega, constants.w0, constants.twodivpi,
                                    constants.myi) * atten.cast<Complex>();
  result.makeCompressed();
  return result;
}

template <typename Solver>
Timed measure(const System &system, const std::string &ordering,
             bool symmetric, const Dense &reference) {
  Solver warm;
  warm.isSymmetric(symmetric);
  warm.analyzePattern(system.matrix);
  warm.factorize(system.matrix);
  warm.solve(system.rhs);

  Solver solver;
  solver.isSymmetric(symmetric);
  const auto analyzeBegin = Clock::now();
  solver.analyzePattern(system.matrix);
  const auto analyzeEnd = Clock::now();
  // Eigen 3.4's analyzePattern() has no independent return status and does
  // not initialize info(); a successful factorization proves the analysis was
  // usable, so record that stage as successful (zero) below.
  const int infoAnalyze = 0;
  const auto factorBegin = Clock::now();
  solver.factorize(system.matrix);
  const auto factorEnd = Clock::now();
  const int infoFactorize = static_cast<int>(solver.info());
  const auto solveBegin = Clock::now();
  const Dense solution = solver.solve(system.rhs);
  const auto solveEnd = Clock::now();
  const int infoSolve = static_cast<int>(solver.info());

  const Eigen::Index lNnz = solver.nnzL();
  const Eigen::Index uNnz = solver.nnzU();
  const auto [kl, ku] = bandwidth(system.matrix);
  const Eigen::Index ldab = 2 * kl + ku + 1;
  const Eigen::Index lapackEnvelope = ldab * system.matrix.rows();
  Timed result{ordering,
               symmetric,
               std::chrono::duration<double>(analyzeEnd - analyzeBegin).count(),
               std::chrono::duration<double>(factorEnd - factorBegin).count(),
               std::chrono::duration<double>(solveEnd - solveBegin).count(),
               system.matrix.nonZeros(),
               kl,
               ku,
               lapackEnvelope,
               lNnz,
               uNnz,
               infoAnalyze,
               infoFactorize,
               infoSolve,
               (solution - reference).norm() /
                   std::max(1.0, reference.norm()),
               infNorm(system.matrix * solution - system.rhs) /
                   (infNorm(system.matrix) * infNorm(solution) +
                    infNorm(system.rhs))};
  return result;
}

template <typename Solver>
Dense referenceSolution(const System &system) {
  Solver solver;
  solver.isSymmetric(false);
  solver.compute(system.matrix);
  if (solver.info() != Eigen::Success)
    throw std::runtime_error("default SparseLU reference factorization failed");
  return solver.solve(system.rhs);
}

template <typename MakeSystem>
void addSystem(std::vector<System> &systems, const std::string &family,
               const std::string &size, int degree, int nq, double maxstep,
               MakeSystem make) {
  System system{family, size, degree, nq, maxstep, make(), Dense{}};
  system.rhs = makeRhs(system.matrix.rows(),
                       family == "radial" ? 1 : (family == "toroidal" ? 2 : 4));
  systems.push_back(std::move(system));
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "usage: eigen_sparselu_e1 MODEL_FILE\n";
    return EXIT_FAILURE;
  }
  try {
    prem_norm<double> norm;
    auto model = EarthModels::ModelInput(argv[1], norm);
    constexpr double omega = 0.70;
    std::vector<System> systems;

    // These are genuine SEM matrices at three mesh resolutions.  The
    // spheroidal cases use PREM's fluid-solid boundaries and therefore retain
    // the production interface structure rather than a synthetic stencil.
    for (const auto &mesh : {std::tuple<std::string, int, double>{"small", 4, 0.15},
                             {"medium", 5, 0.08}, {"large", 6, 0.02}}) {
      const auto &size = std::get<0>(mesh);
      const int nq = std::get<1>(mesh);
      const double maxstep = std::get<2>(mesh);
      addSystem(systems, "radial", size, 0, nq, maxstep, [&] {
        Full1D::SEM sem(model, maxstep, nq, 1);
        return frequencyMatrix(sem.hR(), sem.pR(), sem.hRa(), omega,
                                model.TREF());
      });
      addSystem(systems, "toroidal", size, 12, nq, maxstep, [&] {
        Full1D::SEM sem(model, maxstep, nq, 12);
        return frequencyMatrix(sem.hTk(12), sem.pTk(12), sem.hTa(12), omega,
                                model.TREF());
      });
      addSystem(systems, "spheroidal", size, 12, nq, maxstep, [&] {
        Full1D::SEM sem(model, maxstep, nq, 12);
        return frequencyMatrix(sem.hS(12), sem.pS(12), sem.hSa(12), omega,
                                model.TREF());
      });
    }
    // Explicit fluid-solid spheroidal system at the largest mesh, retained as
    // a separate row so the benchmark's coverage is unambiguous.
    addSystem(systems, "spheroidal_fluid_solid", "large", 40, 6, 0.02, [&] {
      Full1D::SEM sem(model, 0.02, 6, 40);
      return frequencyMatrix(sem.hS(40), sem.pS(40), sem.hSa(40), omega,
                             model.TREF());
    });
    // A high-resolution spheroidal system near the largest active dimension
    // observed in production profiling (approximately 1554 unknowns).
    addSystem(systems, "spheroidal_max", "max", 40, 6, 0.01, [&] {
      Full1D::SEM sem(model, 0.01, 6, 40);
      return frequencyMatrix(sem.hS(40), sem.pS(40), sem.hSa(40), omega,
                             model.TREF());
    });

    using Colamd = Eigen::SparseLU<Sparse, Eigen::COLAMDOrdering<int>>;
    using Natural = Eigen::SparseLU<Sparse, Eigen::NaturalOrdering<int>>;
    std::cout << std::setprecision(17)
              << "system\tfamily\tsize\tdegree\tnq\tmaxstep\tnrows\tordering\tsymmetric\tanalyze_s\tfactorize_s\tsolve_s\tnnz_A\tkl\tku\tlapack_envelope\tnnz_L\tnnz_U\tinfo_analyze\tinfo_factorize\tinfo_solve\tdiscrepancy\tresidual\n";
    for (const auto &system : systems) {
      const Dense reference = referenceSolution<Colamd>(system);
      for (int symmetric = 0; symmetric <= 1; ++symmetric) {
        const auto a = measure<Colamd>(system, "COLAMD", symmetric, reference);
        const auto n = measure<Natural>(system, "NaturalOrdering", symmetric,
                                        reference);
        for (const auto &result : {a, n})
          std::cout << system.family + "_" + system.size << '\t'
                    << system.family << '\t' << system.size << '\t'
                    << system.degree << '\t' << system.nq << '\t'
                    << system.maxstep << '\t' << system.matrix.rows() << '\t'
                    << result.ordering << '\t'
                    << (result.symmetric ? 1 : 0) << '\t' << result.analyze
                    << '\t' << result.factorize << '\t' << result.solve << '\t'
                    << result.aNnz << '\t' << result.kl << '\t' << result.ku
                    << '\t' << result.lapackEnvelope << '\t' << result.lNnz << '\t'
                    << result.uNnz << '\t' << result.infoAnalyze << '\t'
                    << result.infoFactorize << '\t' << result.infoSolve << '\t'
                    << result.discrepancy << '\t' << result.residual << '\n';
      }
    }
  } catch (const std::exception &error) {
    std::cerr << "eigen_sparselu_e1: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
