#include "config.h"

#include <algorithm>
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

#include <Eigen/SparseLU>

extern "C" {
#include <slu_zdefs.h>
}

#include <DSpecM1D/ModelInput>
#include <DSpecM1D/src/FullSpecMultiSem.h>
#include <DSpecM1D/src/LapackBandSolver.h>
#include <DSpecM1D/src/SEM/SEM.h>
#include <DSpecM1D/src/SpecHelpers.h>
#include <DSpecM1D/FrequencyTools>

namespace {
using Complex = std::complex<double>;
using Sparse = Eigen::SparseMatrix<Complex>;
using Dense = Eigen::MatrixXcd;
using Clock = std::chrono::steady_clock;

struct NaturalLU : Eigen::SparseLU<Sparse, Eigen::NaturalOrdering<int>> {
  NaturalLU() { setPivotThreshold(1.0); }
};

struct Components {
  std::string family, name;
  int degree;
  double peakMhz;
  Sparse h, p, ha;
};

void add(std::vector<Components> &out, std::string family, std::string name,
         int degree, double peakMhz, Sparse h, Sparse p, Sparse ha) {
  h.makeCompressed(); p.makeCompressed(); ha.makeCompressed();
  out.push_back({std::move(family), std::move(name), degree, peakMhz,
                 std::move(h), std::move(p), std::move(ha)});
}

std::vector<Components> systems(const EarthModels::ModelInput<double, int> &model) {
  std::vector<Components> out;
  auto radial = [&](const char *name, int nq, double step) {
    Full1D::SEM sem(model, step, nq, 1);
    add(out, "radial", name, 0, .823974609375,
        sem.hR().cast<Complex>().eval(), sem.pR().cast<Complex>().eval(),
        sem.hRa().cast<Complex>().eval());
  };
  auto toroidal = [&](const char *name, int degree, int nq, double step) {
    Full1D::SEM sem(model, step, nq, degree);
    add(out, "toroidal", name, degree, 1.8310546875,
        sem.hTk(degree).cast<Complex>().eval(),
        sem.pTk(degree).cast<Complex>().eval(),
        sem.hTa(degree).cast<Complex>().eval());
  };
  auto spheroidal = [&](const char *name, int degree, int nq, double step,
                        double modeMhz) {
    Full1D::SEM sem(model, step, nq, degree);
    add(out, "spheroidal", name, degree, modeMhz,
        sem.hS(degree).cast<Complex>().eval(),
        sem.pS(degree).cast<Complex>().eval(),
        sem.hSa(degree).cast<Complex>().eval());
  };
  radial("radial_small", 4, .15);
  toroidal("toroidal_small", 12, 4, .15);
  spheroidal("spheroidal_small", 2, 4, .15, .30517578125);
  spheroidal("spheroidal_medium", 12, 5, .08, 1.983642578125);
  spheroidal("spheroidal_large", 12, 6, .02, 1.983642578125);
  spheroidal("spheroidal_fluid_solid", 40, 6, .02, 4.730224609375);
  spheroidal("spheroidal_max", 40, 6, .01, 4.730224609375);
  return out;
}

std::vector<Sparse> matrices(const Components &s, double tref,
                             double timeNorm, double gridMhz) {
  SPARSESPEC::detail::FixedPatternFrequencyMatrix fixed(s.h, s.p, s.ha, true);
  const SPARSESPEC::SpecConstants constants(1.e-4, tref);
  std::vector<Sparse> out;
  for (double fMhz : {s.peakMhz, s.peakMhz + gridMhz}) {
    const double wr = 2.0 * EIGEN_PI * fMhz / 1000.0 * timeNorm;
    const Complex w = wr + constants.ieps;
    const Complex chi = SPARSESPEC::attenFactor(
        wr, constants.w0, constants.twodivpi, constants.myi);
    fixed.update(w, chi);
    out.push_back(fixed.matrix());
  }
  return out;
}

Dense rhs(Eigen::Index n) {
  Dense b(n, 3);
  for (Eigen::Index j = 0; j < b.cols(); ++j)
    for (Eigen::Index i = 0; i < b.rows(); ++i)
      b(i, j) = Complex(.4 + .003 * i + .07 * j,
                        -.1 + .002 * i - .03 * j);
  return b;
}

double infNorm(const Dense &m) {
  double result = 0;
  for (Eigen::Index i = 0; i < m.rows(); ++i) {
    double row = 0;
    for (Eigen::Index j = 0; j < m.cols(); ++j) row += std::abs(m(i, j));
    result = std::max(result, row);
  }
  return result;
}

double infNorm(const Sparse &m) {
  Eigen::VectorXd rows = Eigen::VectorXd::Zero(m.rows());
  for (Eigen::Index col = 0; col < m.outerSize(); ++col)
    for (Sparse::InnerIterator it(m, col); it; ++it)
      rows(it.row()) += std::abs(it.value());
  return rows.maxCoeff();
}

double residual(const Sparse &a, const Dense &x, const Dense &b) {
  return infNorm(a * x - b) /
         (infNorm(a) * infNorm(x) + infNorm(b));
}

double relative(const Dense &a, const Dense &b) {
  return infNorm(a - b) / std::max(infNorm(b), 1.e-300);
}

class SuperLUNatural {
public:
  explicit SuperLUNatural(const Sparse &matrix)
      : n_(checkedInt(matrix.rows())), nnz_(checkedInt(matrix.nonZeros())),
        values_(nnz_), rows_(nnz_), columns_(n_ + 1), permC_(n_), permR_(n_),
        etree_(n_) {
    if (matrix.rows() != matrix.cols()) throw std::runtime_error("non-square matrix");
    for (int i = 0; i <= n_; ++i) columns_[i] = checkedInt(matrix.outerIndexPtr()[i]);
    for (int i = 0; i < nnz_; ++i) rows_[i] = checkedInt(matrix.innerIndexPtr()[i]);
    copyValues(matrix);
    zCreate_CompCol_Matrix(&a_, n_, n_, nnz_, values_.data(), rows_.data(),
                           columns_.data(), SLU_NC, SLU_Z, SLU_GE);
    set_default_options(&options_);
    options_.Fact = DOFACT;
    options_.Equil = NO;
    options_.ColPerm = NATURAL;
    options_.DiagPivotThresh = 1.0;
    options_.IterRefine = NOREFINE;
    options_.PivotGrowth = NO;
    options_.ConditionNumber = NO;
    options_.PrintStat = NO;
    get_perm_c(NATURAL, &a_, permC_.data());
    StatInit(&stat_);
  }

  ~SuperLUNatural() {
    clearFactors();
    if (ac_.Store) Destroy_CompCol_Permuted(&ac_);
    if (a_.Store) Destroy_SuperMatrix_Store(&a_);
    StatFree(&stat_);
  }

  SuperLUNatural(const SuperLUNatural &) = delete;
  SuperLUNatural &operator=(const SuperLUNatural &) = delete;

  double copyValues(const Sparse &matrix) {
    const auto begin = Clock::now();
    if (matrix.rows() != n_ || matrix.nonZeros() != nnz_)
      throw std::runtime_error("SuperLU structure changed");
    for (int i = 0; i < nnz_; ++i) {
      values_[i].r = matrix.valuePtr()[i].real();
      values_[i].i = matrix.valuePtr()[i].imag();
    }
    return seconds(begin);
  }

  double setup() {
    if (ac_.Store) Destroy_CompCol_Permuted(&ac_);
    const auto begin = Clock::now();
    options_.Fact = DOFACT;
    sp_preorder(&options_, &a_, permC_.data(), etree_.data(), &ac_);
    return seconds(begin);
  }

  double factor(fact_t mode) {
    clearFactors();
    options_.Fact = mode;
    const int panel = sp_ienv(1);
    const int relax = sp_ienv(2);
    const auto begin = Clock::now();
    zgstrf(&options_, &ac_, relax, panel, etree_.data(), nullptr, 0,
           permC_.data(), permR_.data(), &l_, &u_, &glu_, &stat_, &info_);
    return seconds(begin);
  }

  std::pair<double, Dense> solve(const Dense &b) {
    Dense x(b.rows(), b.cols());
    std::vector<doublecomplex> values(static_cast<std::size_t>(b.size()));
    for (Eigen::Index j = 0; j < b.cols(); ++j)
      for (Eigen::Index i = 0; i < b.rows(); ++i) {
        const auto k = static_cast<std::size_t>(i + j * b.rows());
        values[k].r = b(i, j).real(); values[k].i = b(i, j).imag();
      }
    SuperMatrix dense{};
    zCreate_Dense_Matrix(&dense, n_, checkedInt(b.cols()), values.data(), n_,
                         SLU_DN, SLU_Z, SLU_GE);
    const auto begin = Clock::now();
    zgstrs(NOTRANS, &l_, &u_, permC_.data(), permR_.data(), &dense, &stat_,
           &solveInfo_);
    const double elapsed = seconds(begin);
    for (Eigen::Index j = 0; j < x.cols(); ++j)
      for (Eigen::Index i = 0; i < x.rows(); ++i) {
        const auto k = static_cast<std::size_t>(i + j * x.rows());
        x(i, j) = Complex(values[k].r, values[k].i);
      }
    Destroy_SuperMatrix_Store(&dense);
    return {elapsed, std::move(x)};
  }

  int info() const { return info_; }
  int solveInfo() const { return solveInfo_; }
  int_t nnzL() const { return l_.Store ? static_cast<SCformat *>(l_.Store)->nnz : 0; }
  int_t nnzU() const { return u_.Store ? static_cast<NCformat *>(u_.Store)->nnz : 0; }

private:
  static int checkedInt(Eigen::Index value) {
    if (value < 0 || value > std::numeric_limits<int>::max())
      throw std::runtime_error("matrix exceeds SuperLU int interface");
    return static_cast<int>(value);
  }
  static double seconds(Clock::time_point begin) {
    return std::chrono::duration<double>(Clock::now() - begin).count();
  }
  void clearFactors() {
    if (l_.Store) { Destroy_SuperNode_Matrix(&l_); l_.Store = nullptr; }
    if (u_.Store) { Destroy_CompCol_Matrix(&u_); u_.Store = nullptr; }
  }

  int n_, nnz_, info_ = 0, solveInfo_ = 0;
  std::vector<doublecomplex> values_;
  std::vector<int> rows_, columns_, permC_, permR_, etree_;
  SuperMatrix a_{}, ac_{}, l_{}, u_{};
  superlu_options_t options_{};
  SuperLUStat_t stat_{};
  GlobalLU_t glu_{};
};

} // namespace

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "usage: superlu_factorization_probe MODEL\n";
    return 2;
  }
  try {
    prem_norm<double> norm;
    auto model = EarthModels::ModelInput(argv[1], norm);
    const auto ss = systems(model);
    SpectraSolver::FreqFull grid(.1, 55.0, 1.0, 4.0, 1.0, .05, 0.0, 4.0,
                                  1, model.TimeNorm());
    const double gridMhz = grid.df() * 1000.0 / model.TimeNorm();
    constexpr int repeats = 16;
    volatile double sink = 0;
    std::cout << std::setprecision(17)
      << "system\tfamily\tfreq_index\tfreq_mhz\tn\tnnz_a\trepeats\teigen_setup_s"
      << "\teigen_factor_s\teigen_solve_s\teigen_nnz_l\teigen_nnz_u"
      << "\teigen_factor_info\teigen_solve_info"
      << "\tsuperlu_conversion_s\tsuperlu_setup_s\tsuperlu_dofact_s"
      << "\tsuperlu_samepattern_s\tsuperlu_solve_s\tsuperlu_nnz_l"
      << "\tsuperlu_nnz_u\tsuperlu_dofact_info\tsuperlu_factor_info\tsuperlu_solve_info"
      << "\tlapack_factor_info\tlapack_solve_info\teigen_residual\tsuperlu_residual"
      << "\tlapack_residual\tsolution_discrepancy\teigen_vs_lapack"
      << "\tsuperlu_vs_lapack\tfinite\n";

    for (const auto &system : ss) {
      auto aa = matrices(system, model.TREF(), model.TimeNorm(), gridMhz);
      auto b = rhs(aa.front().rows());

      NaturalLU eigen;
      auto t = Clock::now(); eigen.analyzePattern(aa.front());
      const double eigenSetup = std::chrono::duration<double>(Clock::now() - t).count();
      eigen.factorize(aa.front()); auto eigenWarm = eigen.solve(b); sink += eigenWarm(0, 0).real();

      SuperLUNatural firstWarm(aa.front()); firstWarm.setup();
      firstWarm.factor(DOFACT); auto superWarm = firstWarm.solve(b); sink += superWarm.second(0, 0).real();

      SuperLUNatural first(aa.front());
      const double firstConversion = first.copyValues(aa.front());
      const double superSetup = first.setup();
      const double dofact = first.factor(DOFACT);
      const int dofactInfo = first.info();

      SuperLUNatural repeated(aa.front()); repeated.setup(); repeated.factor(DOFACT);
      repeated.copyValues(aa.back()); repeated.factor(SamePattern);

      for (std::size_t fi = 0; fi < aa.size(); ++fi) {
        double eigenFactor = 0, eigenSolve = 0, conversion = 0, samePattern = 0;
        Dense xe, xs;
        for (int repeat = 0; repeat < repeats; ++repeat) {
          t = Clock::now(); eigen.factorize(aa[fi]);
          eigenFactor += std::chrono::duration<double>(Clock::now() - t).count();
        }
        const int eigenFactorInfo = static_cast<int>(eigen.info());
        for (int repeat = 0; repeat < repeats; ++repeat) {
          t = Clock::now(); xe = eigen.solve(b);
          eigenSolve += std::chrono::duration<double>(Clock::now() - t).count();
        }
        const int eigenSolveInfo = static_cast<int>(eigen.info());
        for (int repeat = 0; repeat < repeats; ++repeat) {
          conversion += repeated.copyValues(aa[fi]);
          samePattern += repeated.factor(SamePattern);
        }
        double superSolve = 0;
        for (int repeat = 0; repeat < repeats; ++repeat) {
          auto result = repeated.solve(b); superSolve += result.first; xs = std::move(result.second);
        }
        const double er = residual(aa[fi], xe, b);
        const double sr = residual(aa[fi], xs, b);
        const double discrepancy = relative(xs, xe);
        SPARSESPEC::detail::LapackBandSolver lapack;
        lapack.compute(aa[fi]); const int lapackFactorInfo = lapack.info();
        Dense xl = lapack.solve(b); const int lapackSolveInfo = lapack.info();
        const double lr = residual(aa[fi], xl, b);
        const double eigenVsLapack = relative(xe, xl);
        const double superVsLapack = relative(xs, xl);
        const bool finite = xe.allFinite() && xs.allFinite() && xl.allFinite() && std::isfinite(er) &&
                            std::isfinite(sr) && std::isfinite(discrepancy) &&
                            eigen.info() == Eigen::Success && repeated.info() == 0 &&
                            repeated.solveInfo() == 0 && lapackFactorInfo == 0 &&
                            lapackSolveInfo == 0;
        std::cout << system.name << '\t' << system.family << '\t' << fi << '\t'
          << system.peakMhz + static_cast<double>(fi) * gridMhz << '\t'
          << aa[fi].rows() << '\t' << aa[fi].nonZeros() << '\t' << repeats << '\t'
          << eigenSetup << '\t' << eigenFactor / repeats << '\t' << eigenSolve / repeats << '\t'
          << eigen.nnzL() << '\t' << eigen.nnzU() << '\t'
          << eigenFactorInfo << '\t' << eigenSolveInfo << '\t'
          << (fi == 0 ? firstConversion : conversion / repeats) << '\t'
          << superSetup << '\t' << dofact << '\t' << samePattern / repeats << '\t'
          << superSolve / repeats << '\t' << repeated.nnzL() << '\t' << repeated.nnzU() << '\t'
          << dofactInfo << '\t' << repeated.info() << '\t' << repeated.solveInfo() << '\t'
          << lapackFactorInfo << '\t' << lapackSolveInfo << '\t'
          << er << '\t' << sr << '\t' << lr << '\t' << discrepancy << '\t'
          << eigenVsLapack << '\t' << superVsLapack << '\t' << (finite ? 1 : 0) << '\n';
      }
    }
    if (sink == std::numeric_limits<double>::infinity()) std::cerr << sink;
  } catch (const std::exception &e) {
    std::cerr << "superlu_factorization_probe: " << e.what() << '\n';
    return 1;
  }
}
