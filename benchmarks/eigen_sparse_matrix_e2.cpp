#include "config.h"

#include <DSpecM1D/ModelInput>
#include <DSpecM1D/src/SEM/SEM.h>
#include <DSpecM1D/src/SpecHelpers.h>

#include <Eigen/SparseLU>

#include <algorithm>
#include <chrono>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

using Complex = std::complex<double>;
using Sparse = Eigen::SparseMatrix<Complex>;
using Dense = Eigen::MatrixXcd;
using Clock = std::chrono::steady_clock;

struct Components {
  std::string name;
  std::string family;
  int degree;
  int nq;
  double maxstep;
  Sparse h;
  Sparse p;
  Sparse ha;
};

Sparse
currentMatrix(const Components &s, Complex w, bool attenuation, double wval,
              double tref) {
  SPARSESPEC::SpecConstants constants(1.0e-4, tref);
  Sparse result = s.h - (w * w) * s.p;
  if (attenuation)
    result += SPARSESPEC::attenFactor(wval, constants.w0, constants.twodivpi,
                                      constants.myi) *
              s.ha;
  result.makeCompressed();
  return result;
}

class FixedPattern {
public:
  FixedPattern(const Components &s, bool attenuation)
      : matrix(s.h.rows(), s.h.cols()), attenuation_(attenuation) {
    std::map<std::pair<Eigen::Index, Eigen::Index>, Eigen::Index> positions;
    auto collect = [&](const Sparse &part) {
      for (Eigen::Index col = 0; col < part.outerSize(); ++col)
        for (Sparse::InnerIterator entry(part, col); entry; ++entry)
          positions.emplace(std::make_pair(entry.row(), entry.col()), 0);
    };
    collect(s.h);
    collect(s.p);
    if (attenuation)
      collect(s.ha);

    std::vector<Eigen::Triplet<Complex>> triplets;
    triplets.reserve(positions.size());
    for (const auto &entry : positions)
      triplets.emplace_back(entry.first.first, entry.first.second, Complex{});
    matrix.setFromTriplets(triplets.begin(), triplets.end());
    matrix.makeCompressed();
    for (Eigen::Index col = 0; col < matrix.outerSize(); ++col)
      for (Eigen::Index k = matrix.outerIndexPtr()[col];
           k < matrix.outerIndexPtr()[col + 1]; ++k)
        positions.at({matrix.innerIndexPtr()[k], col}) = k;

    hValues_.assign(matrix.nonZeros(), Complex{});
    pValues_.assign(matrix.nonZeros(), Complex{});
    haValues_.assign(matrix.nonZeros(), Complex{});
    copyValues(s.h, hValues_, positions);
    copyValues(s.p, pValues_, positions);
    if (attenuation)
      copyValues(s.ha, haValues_, positions);
    outer_ = std::vector<Eigen::Index>(matrix.outerIndexPtr(),
                                       matrix.outerIndexPtr() +
                                           matrix.outerSize() + 1);
    inner_ = std::vector<Eigen::Index>(
        matrix.innerIndexPtr(), matrix.innerIndexPtr() + matrix.nonZeros());
  }

  void update(Complex w, Complex chi) {
    Complex *values = matrix.valuePtr();
    if (attenuation_) {
      for (Eigen::Index k = 0; k < matrix.nonZeros(); ++k)
        values[k] = hValues_[k] - (w * w) * pValues_[k] + chi * haValues_[k];
    } else {
      for (Eigen::Index k = 0; k < matrix.nonZeros(); ++k)
        values[k] = hValues_[k] - (w * w) * pValues_[k];
    }
  }

  const Sparse &get() const { return matrix; }

  bool structureUnchanged() const {
    return std::equal(outer_.begin(), outer_.end(), matrix.outerIndexPtr()) &&
           std::equal(inner_.begin(), inner_.end(), matrix.innerIndexPtr());
  }

private:
  template <typename Positions>
  void copyValues(const Sparse &part, std::vector<Complex> &values,
                  const Positions &positions) {
    for (Eigen::Index col = 0; col < part.outerSize(); ++col)
      for (Sparse::InnerIterator entry(part, col); entry; ++entry)
        values.at(positions.at({entry.row(), col})) = entry.value();
  }

  Sparse matrix;
  std::vector<Complex> hValues_, pValues_, haValues_;
  std::vector<Eigen::Index> outer_, inner_;
  bool attenuation_;
};

Dense
rhs(Eigen::Index n, int nrhs) {
  Dense result(n, nrhs);
  for (Eigen::Index col = 0; col < nrhs; ++col)
    for (Eigen::Index row = 0; row < n; ++row)
      result(row, col) = Complex(0.4 + 0.003 * row + 0.07 * col,
                                 -0.1 + 0.002 * row - 0.03 * col);
  return result;
}

double
infNorm(const Dense &value) {
  double result = 0.0;
  for (Eigen::Index row = 0; row < value.rows(); ++row) {
    double sum = 0.0;
    for (Eigen::Index col = 0; col < value.cols(); ++col)
      sum += std::abs(value(row, col));
    result = std::max(result, sum);
  }
  return result;
}

std::pair<Eigen::Index, Eigen::Index>
bandwidth(const Sparse &matrix) {
  Eigen::Index kl = 0, ku = 0;
  for (Eigen::Index col = 0; col < matrix.outerSize(); ++col)
    for (Sparse::InnerIterator entry(matrix, col); entry; ++entry) {
      kl = std::max(kl, entry.row() - entry.col());
      ku = std::max(ku, entry.col() - entry.row());
    }
  return {kl, ku};
}

template <typename Make>
void
add(std::vector<Components> &systems, const std::string &family,
    const std::string &name, int degree, int nq, double maxstep, Make make) {
  Components s{name, family, degree, nq, maxstep, {}, {}, {}};
  std::tie(s.h, s.p, s.ha) = make();
  s.h.makeCompressed();
  s.p.makeCompressed();
  s.ha.makeCompressed();
  systems.push_back(std::move(s));
}

}   // namespace

int
main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "usage: eigen_sparse_matrix_e2 MODEL_FILE\n";
    return EXIT_FAILURE;
  }
  try {
    prem_norm<double> norm;
    auto model = EarthModels::ModelInput(argv[1], norm);
    std::vector<Components> systems;
    for (const auto &mesh :
         {std::tuple<std::string, int, double>{"small", 4, 0.15},
          {"large", 6, 0.02}}) {
      const auto &name = std::get<0>(mesh);
      const int nq = std::get<1>(mesh);
      const double maxstep = std::get<2>(mesh);
      add(systems, "radial", "radial_" + name, 0, nq, maxstep, [&] {
        Full1D::SEM sem(model, maxstep, nq, 1);
        return std::make_tuple(sem.hR().cast<Complex>().eval(),
                               sem.pR().cast<Complex>().eval(),
                               sem.hRa().cast<Complex>().eval());
      });
      add(systems, "toroidal", "toroidal_" + name, 12, nq, maxstep, [&] {
        Full1D::SEM sem(model, maxstep, nq, 12);
        return std::make_tuple(sem.hTk(12).cast<Complex>().eval(),
                               sem.pTk(12).cast<Complex>().eval(),
                               sem.hTa(12).cast<Complex>().eval());
      });
      add(systems, "spheroidal", "spheroidal_" + name, 12, nq, maxstep, [&] {
        Full1D::SEM sem(model, maxstep, nq, 12);
        return std::make_tuple(sem.hS(12).cast<Complex>().eval(),
                               sem.pS(12).cast<Complex>().eval(),
                               sem.hSa(12).cast<Complex>().eval());
      });
    }
    add(systems, "spheroidal_fluid_solid", "spheroidal_fluid_solid_large", 40,
        6, 0.02, [&] {
          Full1D::SEM sem(model, 0.02, 6, 40);
          return std::make_tuple(sem.hS(40).cast<Complex>().eval(),
                                 sem.pS(40).cast<Complex>().eval(),
                                 sem.hSa(40).cast<Complex>().eval());
        });

    constexpr int batch = 64;
    constexpr double omega0 = 0.4;
    constexpr double domega = 0.003;
    std::cout
        << std::setprecision(17)
        << "system\tfamily\tn\tnnz_union_off\tnnz_union_"
           "on\tkl\tku\tattenuation\tcurrent_s\tdirect_"
           "s\tspeedup\tstructure_ok\tmax_coeff_diff\tmax_relative_coeff_"
           "diff\tcurrent_analyze_completed\tcurrent_factorize_info\tcurrent_"
           "solve_info\tdirect_analyze_completed\tdirect_factorize_"
           "info\tdirect_"
           "solve_info\tsolve_discrepancy\tresidual\n";
    volatile double sink = 0.0;
    for (const auto &system : systems) {
      for (const bool attenuation : {false, true}) {
        FixedPattern direct(system, attenuation);
        FixedPattern directOff(system, false);
        const auto [kl, ku] = bandwidth(system.h);
        for (int i = 0; i < batch; ++i) {
          const double wval = omega0 + domega * i;
          const Complex w =
              wval + SPARSESPEC::SpecConstants(1.0e-4, model.TREF()).ieps;
          const auto chi =
              attenuation
                  ? SPARSESPEC::attenFactor(
                        wval,
                        SPARSESPEC::SpecConstants(1.0e-4, model.TREF()).w0,
                        SPARSESPEC::SpecConstants(1.0e-4, model.TREF())
                            .twodivpi,
                        SPARSESPEC::SpecConstants(1.0e-4, model.TREF()).myi)
                  : Complex{};
          sink += currentMatrix(system, w, attenuation, wval, model.TREF())
                      .coeff(0, 0)
                      .real();
          direct.update(w, chi);
          sink += direct.get().coeff(0, 0).real();
        }
        const Complex warmW =
            omega0 + SPARSESPEC::SpecConstants(1.0e-4, model.TREF()).ieps;
        auto currentWarm =
            currentMatrix(system, warmW, attenuation, omega0, model.TREF());
        direct.update(
            warmW,
            attenuation
                ? SPARSESPEC::attenFactor(
                      omega0,
                      SPARSESPEC::SpecConstants(1.0e-4, model.TREF()).w0,
                      SPARSESPEC::SpecConstants(1.0e-4, model.TREF()).twodivpi,
                      SPARSESPEC::SpecConstants(1.0e-4, model.TREF()).myi)
                : Complex{});
        double currentSeconds = 0.0;
        double directSeconds = 0.0;
        bool structureOk = direct.structureUnchanged();
        double coefficientDifference = 0.0;
        double relativeDifference = 0.0;
        for (int i = 0; i < batch; ++i) {
          const double wval = omega0 + domega * i;
          const Complex w =
              wval + SPARSESPEC::SpecConstants(1.0e-4, model.TREF()).ieps;
          const auto chi =
              attenuation
                  ? SPARSESPEC::attenFactor(
                        wval,
                        SPARSESPEC::SpecConstants(1.0e-4, model.TREF()).w0,
                        SPARSESPEC::SpecConstants(1.0e-4, model.TREF())
                            .twodivpi,
                        SPARSESPEC::SpecConstants(1.0e-4, model.TREF()).myi)
                  : Complex{};
          const auto begin = Clock::now();
          auto current =
              currentMatrix(system, w, attenuation, wval, model.TREF());
          currentSeconds +=
              std::chrono::duration<double>(Clock::now() - begin).count();
          const auto directBegin = Clock::now();
          direct.update(w, chi);
          directSeconds +=
              std::chrono::duration<double>(Clock::now() - directBegin).count();
          structureOk = structureOk &&
                        current.outerSize() == direct.get().outerSize() &&
                        current.nonZeros() == direct.get().nonZeros();
          if (structureOk)
            structureOk =
                std::equal(current.outerIndexPtr(),
                           current.outerIndexPtr() + current.outerSize() + 1,
                           direct.get().outerIndexPtr()) &&
                std::equal(current.innerIndexPtr(),
                           current.innerIndexPtr() + current.nonZeros(),
                           direct.get().innerIndexPtr());
          if (current.rows() == direct.get().rows() &&
              current.nonZeros() == direct.get().nonZeros())
            for (Eigen::Index k = 0; k < current.nonZeros(); ++k) {
              coefficientDifference = std::max(
                  coefficientDifference,
                  std::abs(current.valuePtr()[k] - direct.get().valuePtr()[k]));
              relativeDifference = std::max(
                  relativeDifference,
                  std::abs(current.valuePtr()[k] - direct.get().valuePtr()[k]) /
                      std::max(1.0, std::abs(current.valuePtr()[k])));
            }
        }
        using Solver = Eigen::SparseLU<Sparse, Eigen::NaturalOrdering<int>>;
        Dense b = rhs(system.h.rows(), system.family == "radial" ? 1 : 3);
        Solver oldSolver, newSolver;
        oldSolver.analyzePattern(currentWarm);
        const int currentAnalyzeCompleted = 1;
        oldSolver.factorize(currentWarm);
        const int currentFactorizeInfo = static_cast<int>(oldSolver.info());
        const Dense oldSolution = oldSolver.solve(b);
        const int currentSolveInfo = static_cast<int>(oldSolver.info());
        direct.update(
            warmW,
            attenuation
                ? SPARSESPEC::attenFactor(
                      omega0,
                      SPARSESPEC::SpecConstants(1.0e-4, model.TREF()).w0,
                      SPARSESPEC::SpecConstants(1.0e-4, model.TREF()).twodivpi,
                      SPARSESPEC::SpecConstants(1.0e-4, model.TREF()).myi)
                : Complex{});
        newSolver.analyzePattern(direct.get());
        const int directAnalyzeCompleted = 1;
        newSolver.factorize(direct.get());
        const int directFactorizeInfo = static_cast<int>(newSolver.info());
        const Dense newSolution = newSolver.solve(b);
        const int directSolveInfo = static_cast<int>(newSolver.info());
        const double discrepancy = (newSolution - oldSolution).norm() /
                                   std::max(1.0, oldSolution.norm());
        const double residual =
            infNorm(direct.get() * newSolution - b) /
            (infNorm(direct.get()) * infNorm(newSolution) + infNorm(b));
        std::cout << system.name << '\t' << system.family << '\t'
                  << system.h.rows() << '\t' << directOff.get().nonZeros()
                  << '\t' << direct.get().nonZeros() << '\t' << kl << '\t' << ku
                  << '\t' << (attenuation ? 1 : 0) << '\t' << currentSeconds
                  << '\t' << directSeconds << '\t'
                  << currentSeconds / directSeconds << '\t'
                  << (structureOk ? 1 : 0) << '\t' << coefficientDifference
                  << '\t' << relativeDifference << '\t'
                  << currentAnalyzeCompleted << '\t' << currentFactorizeInfo
                  << '\t' << currentSolveInfo << '\t' << directAnalyzeCompleted
                  << '\t' << directFactorizeInfo << '\t' << directSolveInfo
                  << '\t' << discrepancy << '\t' << residual << '\n';
      }
    }
    std::cerr << "sink=" << sink << '\n';
  } catch (const std::exception &error) {
    std::cerr << "eigen_sparse_matrix_e2: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
