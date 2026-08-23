#include "config.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <Eigen/SparseLU>
#include <umfpack.h>

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

static_assert(std::is_same_v<Sparse::StorageIndex, int>);
static_assert(sizeof(Complex) == 2 * sizeof(double));

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
  h.makeCompressed();
  p.makeCompressed();
  ha.makeCompressed();
  out.push_back({std::move(family), std::move(name), degree, peakMhz,
                 std::move(h), std::move(p), std::move(ha)});
}

std::vector<Components> systems(
    const EarthModels::ModelInput<double, int> &model) {
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
    for (Eigen::Index j = 0; j < m.cols(); ++j)
      row += std::abs(m(i, j));
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

double transferDiscrepancy(const Dense &a, const Dense &b) {
  const std::array<Eigen::Index, 3> rows{0, a.rows() / 2, a.rows() - 1};
  double numerator = 0, denominator = 0;
  for (const auto row : rows)
    for (Eigen::Index col = 0; col < a.cols(); ++col) {
      numerator = std::max(numerator, std::abs(a(row, col) - b(row, col)));
      denominator = std::max(denominator, std::abs(b(row, col)));
    }
  return numerator / std::max(denominator, 1.e-300);
}

double seconds(Clock::time_point begin) {
  return std::chrono::duration<double>(Clock::now() - begin).count();
}

const int32_t *indices(const Sparse::StorageIndex *p) {
  return reinterpret_cast<const int32_t *>(p);
}

const double *packed(const Complex *p) {
  // C++ guarantees std::complex<T> array elements can be accessed as an
  // interleaved array of T. UMFPACK's packed-complex form is identical.
  return reinterpret_cast<const double *>(p);
}

double *packed(Complex *p) { return reinterpret_cast<double *>(p); }

struct UmfResult {
  double elapsed = 0;
  Dense x;
  int status = UMFPACK_OK;
  double refinementSteps = 0;
};

class UmfpackLU {
public:
  UmfpackLU(const Sparse &matrix, int ordering, bool robust)
      : n_(checkedInt(matrix.rows())), ordering_(ordering), robust_(robust) {
    if (!matrix.isCompressed() || matrix.rows() != matrix.cols())
      throw std::runtime_error("UMFPACK requires square compressed CSC");
    umfpack_zi_defaults(control_.data());
    control_[UMFPACK_ORDERING] = ordering;
    if (robust) {
      control_[UMFPACK_PIVOT_TOLERANCE] = 1.0;
      control_[UMFPACK_SYM_PIVOT_TOLERANCE] = 1.0;
    }
    const auto begin = Clock::now();
    symbolicStatus_ = umfpack_zi_symbolic(
        n_, n_, indices(matrix.outerIndexPtr()),
        indices(matrix.innerIndexPtr()), packed(matrix.valuePtr()), nullptr,
        &symbolic_, control_.data(), symbolicInfo_.data());
    symbolicSeconds_ = seconds(begin);
  }

  ~UmfpackLU() {
    umfpack_zi_free_numeric(&numeric_);
    umfpack_zi_free_symbolic(&symbolic_);
  }

  UmfpackLU(const UmfpackLU &) = delete;
  UmfpackLU &operator=(const UmfpackLU &) = delete;

  double factorize(const Sparse &matrix) {
    validatePattern(matrix);
    umfpack_zi_free_numeric(&numeric_);
    const auto begin = Clock::now();
    numericStatus_ = umfpack_zi_numeric(
        indices(matrix.outerIndexPtr()), indices(matrix.innerIndexPtr()),
        packed(matrix.valuePtr()), nullptr, symbolic_, &numeric_,
        control_.data(), numericInfo_.data());
    return seconds(begin);
  }

  UmfResult solve(const Sparse &matrix, const Dense &b, bool refined) const {
    UmfResult result;
    result.x.resize(b.rows(), b.cols());
    auto solveControl = control_;
    solveControl[UMFPACK_IRSTEP] = refined ? UMFPACK_DEFAULT_IRSTEP : 0;
    const auto begin = Clock::now();
    for (Eigen::Index col = 0; col < b.cols(); ++col) {
      std::array<double, UMFPACK_INFO> info{};
      const int status = umfpack_zi_solve(
          UMFPACK_A, indices(matrix.outerIndexPtr()),
          indices(matrix.innerIndexPtr()), packed(matrix.valuePtr()), nullptr,
          packed(result.x.col(col).data()), nullptr, packed(b.col(col).data()),
          nullptr, numeric_, solveControl.data(), info.data());
      if (status != UMFPACK_OK && result.status == UMFPACK_OK)
        result.status = status;
      result.refinementSteps += info[UMFPACK_IR_TAKEN];
    }
    result.elapsed = seconds(begin);
    return result;
  }

  int symbolicStatus() const { return symbolicStatus_; }
  int numericStatus() const { return numericStatus_; }
  double symbolicSeconds() const { return symbolicSeconds_; }
  double info(int index) const { return numericInfo_[index]; }
  double symbolicInfo(int index) const { return symbolicInfo_[index]; }
  int ordering() const { return ordering_; }
  bool robust() const { return robust_; }

  void getLunz(int32_t &lnz, int32_t &unz) const {
    int32_t nr = 0, nc = 0, nzdiag = 0;
    const int status = umfpack_zi_get_lunz(&lnz, &unz, &nr, &nc, &nzdiag,
                                           numeric_);
    if (status != UMFPACK_OK || nr != n_ || nc != n_)
      throw std::runtime_error("umfpack_zi_get_lunz failed");
  }

private:
  static int32_t checkedInt(Eigen::Index n) {
    if (n < 0 || n > std::numeric_limits<int32_t>::max())
      throw std::runtime_error("matrix exceeds umfpack_zi int32 interface");
    return static_cast<int32_t>(n);
  }

  void validatePattern(const Sparse &matrix) const {
    if (!matrix.isCompressed() || matrix.rows() != n_ ||
        matrix.cols() != n_)
      throw std::runtime_error("UMFPACK matrix dimensions changed");
  }

  int32_t n_;
  int ordering_;
  bool robust_;
  void *symbolic_ = nullptr;
  void *numeric_ = nullptr;
  int symbolicStatus_ = UMFPACK_ERROR_internal_error;
  int numericStatus_ = UMFPACK_ERROR_internal_error;
  double symbolicSeconds_ = 0;
  std::array<double, UMFPACK_CONTROL> control_{};
  std::array<double, UMFPACK_INFO> symbolicInfo_{};
  std::array<double, UMFPACK_INFO> numericInfo_{};
};

const char *orderingName(int ordering) {
  return ordering == UMFPACK_ORDERING_NONE ? "natural" : "default_amd";
}

const char *policyName(bool robust) { return robust ? "robust" : "default"; }

} // namespace

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "usage: umfpack_factorization_probe MODEL\n";
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
      << "system\tfamily\tfreq_index\tfreq_mhz\tn\tnnz_a\trepeats\tpolicy\tordering"
      << "\teigen_setup_s\teigen_factor_s\teigen_solve_s\teigen_factor_info\teigen_solve_info"
      << "\tumf_symbolic_s\tumf_numeric_s\tumf_raw_solve_s\tumf_refined_solve_s"
      << "\tumf_symbolic_status\tumf_numeric_status\tumf_raw_status\tumf_refined_status"
      << "\tumf_strategy_used\tumf_ordering_used\tumf_was_scaled\tumf_noff_diag"
      << "\tumf_lnz\tumf_unz\tumf_flops\tumf_numeric_size\tumf_peak_memory\tumf_rcond"
      << "\tumf_raw_ir_steps\tumf_refined_ir_steps\tlapack_factor_info\tlapack_solve_info"
      << "\teigen_residual\tumf_raw_residual\tumf_refined_residual\tlapack_residual"
      << "\tumf_raw_vs_eigen\tumf_refined_vs_eigen\teigen_vs_lapack"
      << "\tumf_raw_vs_lapack\tumf_refined_vs_lapack\ttransfer_raw_vs_eigen"
      << "\ttransfer_refined_vs_eigen\tfinite\n";

    for (const auto &system : ss) {
      auto aa = matrices(system, model.TREF(), model.TimeNorm(), gridMhz);
      const Dense b = rhs(aa.front().rows());
      if (!std::equal(aa[0].outerIndexPtr(), aa[0].outerIndexPtr() + aa[0].outerSize() + 1,
                      aa[1].outerIndexPtr()) ||
          !std::equal(aa[0].innerIndexPtr(), aa[0].innerIndexPtr() + aa[0].nonZeros(),
                      aa[1].innerIndexPtr()))
        throw std::runtime_error("fixed pattern changed between frequencies");

      NaturalLU eigenWarm;
      eigenWarm.analyzePattern(aa.front());
      eigenWarm.factorize(aa.front());
      sink += eigenWarm.solve(b)(0, 0).real();

      NaturalLU eigen;
      auto begin = Clock::now();
      eigen.analyzePattern(aa.front());
      const double eigenSetup = seconds(begin);
      eigen.factorize(aa.front());
      sink += eigen.solve(b)(0, 0).real();

      for (const int ordering : {UMFPACK_ORDERING_NONE,
                                 UMFPACK_DEFAULT_ORDERING}) {
        for (const bool robust : {true, false}) {
          {
            UmfpackLU warm(aa.front(), ordering, robust);
            warm.factorize(aa.front());
            sink += warm.solve(aa.front(), b, false).x(0, 0).real();
          }
          UmfpackLU umf(aa.front(), ordering, robust);
          if (umf.symbolicStatus() != UMFPACK_OK)
            throw std::runtime_error("UMFPACK symbolic failed");

          for (std::size_t fi = 0; fi < aa.size(); ++fi) {
            double eigenFactor = 0, eigenSolve = 0, umfNumeric = 0;
            Dense xe;
            for (int repeat = 0; repeat < repeats; ++repeat) {
              begin = Clock::now();
              eigen.factorize(aa[fi]);
              eigenFactor += seconds(begin);
            }
            const int eigenFactorInfo = static_cast<int>(eigen.info());
            for (int repeat = 0; repeat < repeats; ++repeat) {
              begin = Clock::now();
              xe = eigen.solve(b);
              eigenSolve += seconds(begin);
            }
            const int eigenSolveInfo = static_cast<int>(eigen.info());
            for (int repeat = 0; repeat < repeats; ++repeat)
              umfNumeric += umf.factorize(aa[fi]);

            UmfResult raw, refined;
            for (int repeat = 0; repeat < repeats; ++repeat) {
              auto r = umf.solve(aa[fi], b, false);
              raw.elapsed += r.elapsed;
              raw.x = std::move(r.x);
              raw.status = r.status;
              raw.refinementSteps += r.refinementSteps;
            }
            for (int repeat = 0; repeat < repeats; ++repeat) {
              auto r = umf.solve(aa[fi], b, true);
              refined.elapsed += r.elapsed;
              refined.x = std::move(r.x);
              refined.status = r.status;
              refined.refinementSteps += r.refinementSteps;
            }

            int32_t lnz = 0, unz = 0;
            umf.getLunz(lnz, unz);
            SPARSESPEC::detail::LapackBandSolver lapack;
            lapack.compute(aa[fi]);
            const int lapackFactorInfo = lapack.info();
            const Dense xl = lapack.solve(b);
            const int lapackSolveInfo = lapack.info();
            const double er = residual(aa[fi], xe, b);
            const double urr = residual(aa[fi], raw.x, b);
            const double ufr = residual(aa[fi], refined.x, b);
            const double lr = residual(aa[fi], xl, b);
            const bool finite = xe.allFinite() && raw.x.allFinite() &&
                refined.x.allFinite() && xl.allFinite() &&
                std::isfinite(er) && std::isfinite(urr) &&
                std::isfinite(ufr) && std::isfinite(lr) &&
                eigenFactorInfo == Eigen::Success &&
                eigenSolveInfo == Eigen::Success &&
                umf.numericStatus() == UMFPACK_OK &&
                raw.status == UMFPACK_OK && refined.status == UMFPACK_OK &&
                lapackFactorInfo == 0 && lapackSolveInfo == 0;

            std::cout << system.name << '\t' << system.family << '\t' << fi << '\t'
              << system.peakMhz + static_cast<double>(fi) * gridMhz << '\t'
              << aa[fi].rows() << '\t' << aa[fi].nonZeros() << '\t' << repeats << '\t'
              << policyName(robust) << '\t' << orderingName(ordering) << '\t'
              << eigenSetup << '\t' << eigenFactor / repeats << '\t'
              << eigenSolve / repeats << '\t' << eigenFactorInfo << '\t'
              << eigenSolveInfo << '\t' << umf.symbolicSeconds() << '\t'
              << umfNumeric / repeats << '\t' << raw.elapsed / repeats << '\t'
              << refined.elapsed / repeats << '\t' << umf.symbolicStatus() << '\t'
              << umf.numericStatus() << '\t' << raw.status << '\t' << refined.status << '\t'
              << umf.symbolicInfo(UMFPACK_STRATEGY_USED) << '\t'
              << umf.symbolicInfo(UMFPACK_ORDERING_USED) << '\t'
              << umf.info(UMFPACK_WAS_SCALED) << '\t' << umf.info(UMFPACK_NOFF_DIAG) << '\t'
              << lnz << '\t' << unz << '\t' << umf.info(UMFPACK_FLOPS) << '\t'
              << umf.info(UMFPACK_NUMERIC_SIZE) << '\t' << umf.info(UMFPACK_PEAK_MEMORY) << '\t'
              << umf.info(UMFPACK_RCOND) << '\t' << raw.refinementSteps / repeats << '\t'
              << refined.refinementSteps / repeats << '\t' << lapackFactorInfo << '\t'
              << lapackSolveInfo << '\t' << er << '\t' << urr << '\t' << ufr << '\t'
              << lr << '\t' << relative(raw.x, xe) << '\t'
              << relative(refined.x, xe) << '\t' << relative(xe, xl) << '\t'
              << relative(raw.x, xl) << '\t' << relative(refined.x, xl) << '\t'
              << transferDiscrepancy(raw.x, xe) << '\t'
              << transferDiscrepancy(refined.x, xe) << '\t' << (finite ? 1 : 0) << '\n';
          }
        }
      }
    }
    if (sink == std::numeric_limits<double>::infinity())
      std::cerr << sink;
  } catch (const std::exception &e) {
    std::cerr << "umfpack_factorization_probe: " << e.what() << '\n';
    return 1;
  }
}
