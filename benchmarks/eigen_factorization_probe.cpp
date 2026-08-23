#include "config.h"

#include <complex>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include <Eigen/SparseLU>

#include <DSpecM1D/ModelInput>
#include <DSpecM1D/src/FullSpecMultiSem.h>
#include <DSpecM1D/src/SEM/SEM.h>
#include <DSpecM1D/src/SpecHelpers.h>

namespace {
using Complex = std::complex<double>;
using Sparse = Eigen::SparseMatrix<Complex>;
using Dense = Eigen::MatrixXcd;
using Clock = std::chrono::steady_clock;

template <class Matrix>
class TunableNaturalLU
    : public Eigen::SparseLU<Matrix, Eigen::NaturalOrdering<int>> {
  using Base = Eigen::SparseLU<Matrix, Eigen::NaturalOrdering<int>>;
public:
  void setPanel(int v) { this->m_perfv.panel_size = v; }
  void setRelax(int v) { this->m_perfv.relax = v; }
  void setMaxsuper(int v) { this->m_perfv.maxsuper = v; }
};

struct System {
  std::string family, name;
  Sparse h, p, ha;
};

struct Config {
  std::string name, kind;
  double threshold = 1.0;
  int panel = 16, relax = 1, maxsuper = 128;
};

struct Baseline {
  Dense x;
  double residual = std::numeric_limits<double>::infinity();
};

double infNorm(const Dense &m) {
  double r = 0;
  for (Eigen::Index i = 0; i < m.rows(); ++i) {
    double s = 0;
    for (Eigen::Index j = 0; j < m.cols(); ++j) s += std::abs(m(i, j));
    r = std::max(r, s);
  }
  return r;
}

double matrixInfNorm(const Sparse &m) {
  double r = 0;
  for (Eigen::Index i = 0; i < m.rows(); ++i) {
    double s = 0;
    for (Eigen::Index col = 0; col < m.outerSize(); ++col)
      for (Sparse::InnerIterator it(m, col); it; ++it)
        if (it.row() == i) s += std::abs(it.value());
    r = std::max(r, s);
  }
  return r;
}

double backwardResidual(const Sparse &a, const Dense &x, const Dense &b) {
  return infNorm(a * x - b) /
         (matrixInfNorm(a) * infNorm(x) + infNorm(b));
}

Dense rhs(Eigen::Index n) {
  Dense b(n, 2);
  for (Eigen::Index j = 0; j < b.cols(); ++j)
    for (Eigen::Index i = 0; i < b.rows(); ++i)
      b(i, j) = Complex(.4 + .003 * i + .07 * j, -.1 + .002 * i - .03 * j);
  return b;
}

std::vector<Sparse> fixedPatternMatrices(const System &s,
                                         const std::vector<double> &freqs,
                                         double tref) {
  // Reuse the exact production E2 helper: union pattern and component values
  // are prepared once, then only valuePtr() is updated for each frequency.
  SPARSESPEC::detail::FixedPatternFrequencyMatrix fixed(s.h, s.p, s.ha, true);
  const auto c = SPARSESPEC::SpecConstants(1.e-4, tref);
  std::vector<Sparse> result;
  result.reserve(freqs.size());
  for (double f : freqs) {
    const Complex w = f + c.ieps;
    const Complex chi = SPARSESPEC::attenFactor(f, c.w0, c.twodivpi, c.myi);
    fixed.update(w, chi);
    result.push_back(fixed.matrix());
  }
  return result;
}

void add(std::vector<System> &out, const std::string &family,
         const std::string &name, const Sparse &h, const Sparse &p,
         const Sparse &ha) {
  System s{family, name, h, p, ha};
  s.h.makeCompressed(); s.p.makeCompressed(); s.ha.makeCompressed();
  out.push_back(std::move(s));
}

std::vector<System> systems(const EarthModels::ModelInput<double, int> &model) {
  std::vector<System> out;
  auto radial = [&](const std::string &name, int nq, double step) {
    Full1D::SEM sem(model, step, nq, 1);
    add(out, "radial", name, sem.hR().cast<Complex>().eval(),
        sem.pR().cast<Complex>().eval(), sem.hRa().cast<Complex>().eval());
  };
  auto toroidal = [&](const std::string &name, int degree, int nq,
                      double step) {
    Full1D::SEM sem(model, step, nq, degree);
    add(out, "toroidal", name, sem.hTk(degree).cast<Complex>().eval(),
        sem.pTk(degree).cast<Complex>().eval(),
        sem.hTa(degree).cast<Complex>().eval());
  };
  auto spheroidal = [&](const std::string &name, int degree, int nq,
                        double step) {
    Full1D::SEM sem(model, step, nq, degree);
    add(out, "spheroidal", name, sem.hS(degree).cast<Complex>().eval(),
        sem.pS(degree).cast<Complex>().eval(),
        sem.hSa(degree).cast<Complex>().eval());
  };
  radial("radial_small", 4, .15);
  toroidal("toroidal_small", 12, 4, .15);
  spheroidal("spheroidal_small", 12, 4, .15);
  spheroidal("spheroidal_medium", 12, 5, .08);
  spheroidal("spheroidal_large", 12, 6, .02);
  spheroidal("spheroidal_fluid_solid", 40, 6, .02);
  spheroidal("spheroidal_max", 40, 6, .01);
  return out;
}

std::vector<Config> configs() {
  return { {"pivot_1", "pivot", 1}, {"pivot_0p5", "pivot", .5},
           {"pivot_0p1", "pivot", .1}, {"pivot_0p01", "pivot", .01},
           {"pivot_0", "pivot", 0}, {"natural_default", "perf"},
           {"panel_1", "perf", 1, 1, 1}, {"panel_4", "perf", 1, 4, 1},
           {"panel_8", "perf", 1, 8, 1}, {"panel_32", "perf", 1, 32, 1},
           {"relax_4", "perf", 1, 16, 4},
           {"maxsuper_64", "perf", 1, 16, 1, 64},
           {"maxsuper_256", "perf", 1, 16, 1, 256},
           {"panel_4_relax_4", "perf", 1, 4, 4},
           {"panel_4_maxsuper_64", "perf", 1, 4, 1, 64} };
}

std::string key(std::size_t s, std::size_t f) {
  return std::to_string(s) + ":" + std::to_string(f);
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 2) { std::cerr << "usage: eigen_factorization_probe MODEL\n"; return 2; }
  try {
    prem_norm<double> norm;
    auto model = EarthModels::ModelInput(argv[1], norm);
    const auto ss = systems(model);
    std::vector<std::vector<double>> frequencies(
        ss.size(), std::vector<double>{.18, .41, .67, 1.05});
    // A cheap NaturalOrdering scan selects two additional stress frequencies
    // per real SEM system.  This makes the difficult cases data-driven while
    // retaining the same four representative frequencies for every system.
    const std::vector<double> scanFrequencies{
        .10, .15, .20, .25, .30, .40, .50, .60, .70, .80, .90, 1.00, 1.10, 1.20};
    for (std::size_t si = 0; si < ss.size(); ++si) {
      TunableNaturalLU<Sparse> scan;
      auto first = fixedPatternMatrices(ss[si], {scanFrequencies.front()},
                                        model.TREF()).front();
      scan.analyzePattern(first);
      double worstResidual = -1, worstDisplacement = -1;
      double worstResidualF = scanFrequencies.front();
      double worstDisplacementF = scanFrequencies.front();
      for (double f : scanFrequencies) {
        auto a = fixedPatternMatrices(ss[si], {f}, model.TREF()).front();
        auto b = rhs(a.rows());
        scan.factorize(a);
        auto x = scan.solve(b);
        const double r = backwardResidual(a, x, b);
        int displacement = 0;
        const auto p = scan.rowsPermutation().indices();
        for (Eigen::Index i = 0; i < p.size(); ++i)
          displacement = std::max(displacement,
                                   std::abs(static_cast<int>(p(i) - i)));
        if (r > worstResidual) { worstResidual = r; worstResidualF = f; }
        if (displacement > worstDisplacement) {
          worstDisplacement = displacement; worstDisplacementF = f;
        }
      }
      auto addIfNew = [&](double f) {
        for (double old : frequencies[si])
          if (std::abs(old - f) < 1.e-12) return;
        frequencies[si].push_back(f);
      };
      addIfNew(worstResidualF); addIfNew(worstDisplacementF);
      std::cerr << "difficulty_scan\t" << ss[si].name
                << "\tworst_backward_residual\t" << worstResidualF << "\t"
                << worstResidual << "\tworst_pivot_displacement\t"
                << worstDisplacementF << "\t" << worstDisplacement << '\n';
    }
    std::vector<std::vector<Baseline>> reference(ss.size());
    for (std::size_t si = 0; si < ss.size(); ++si) {
      reference[si].resize(frequencies[si].size());
      for (std::size_t fi = 0; fi < frequencies[si].size(); ++fi) {
        auto a = fixedPatternMatrices(ss[si], {frequencies[si][fi]},
                                      model.TREF()).front();
        auto b = rhs(a.rows());
        TunableNaturalLU<Sparse> solver;
        solver.analyzePattern(a); solver.setPivotThreshold(1.0); solver.factorize(a);
        auto x = solver.solve(b);
        reference[si][fi] = {x, backwardResidual(a, x, b)};
      }
    }

    std::cout << std::setprecision(17)
      << "config\tkind\tthreshold\tpanel\trelax\tmaxsuper\tsystem\tfamily"
      << "\tfreq\tn\tnnz_a\tfactor_s\tsolve_s\tnnz_l\tnnz_u\tnonidentity_pivots"
      << "\tpivot_fraction\tmax_pivot_displacement\tfactor_info\tsolve_info"
      << "\tdiscrepancy\tresidual\tfinite\trepeats\n";
    volatile double sink = 0;
    for (const auto &cfg : configs()) {
      for (std::size_t si = 0; si < ss.size(); ++si) {
        std::vector<Sparse> aa;
        std::vector<Dense> bb;
        aa = fixedPatternMatrices(ss[si], frequencies[si], model.TREF());
        for (const auto &a : aa) bb.push_back(rhs(a.rows()));
        TunableNaturalLU<Sparse> solver;
        solver.setPivotThreshold(cfg.threshold);
        if (cfg.kind == "perf") {
          solver.setPanel(cfg.panel); solver.setRelax(cfg.relax);
          solver.setMaxsuper(cfg.maxsuper);
        }
        solver.analyzePattern(aa.front());
        // One untimed warm-up factorization/solve for this configuration/system.
        solver.factorize(aa.front());
        auto warm = solver.solve(bb.front());
        sink += warm(0, 0).real();
        constexpr int repeats = 4;
        for (std::size_t fi = 0; fi < frequencies[si].size(); ++fi) {
          int factorInfo = 0, solveInfo = 0;
          Dense x;
          const auto tf0 = Clock::now();
          for (int repeat = 0; repeat < repeats; ++repeat) {
            solver.factorize(aa[fi]);
            factorInfo = solver.info();
          }
          const auto tf1 = Clock::now();
          const auto ts0 = Clock::now();
          for (int repeat = 0; repeat < repeats; ++repeat) {
            x = solver.solve(bb[fi]);
            solveInfo = solver.info();
          }
          const auto ts1 = Clock::now();
          const double factorS = std::chrono::duration<double>(tf1 - tf0).count() / repeats;
          const double solveS = std::chrono::duration<double>(ts1 - ts0).count() / repeats;
          const auto p = solver.rowsPermutation().indices();
          int moved = 0, maxdisp = 0;
          for (Eigen::Index i = 0; i < p.size(); ++i) {
            if (p(i) != i) ++moved;
            maxdisp = std::max(maxdisp, std::abs(static_cast<int>(p(i) - i)));
          }
          const double discrepancy =
              infNorm(x - reference[si][fi].x) /
              std::max(1.0, infNorm(reference[si][fi].x));
          const double residual = backwardResidual(aa[fi], x, bb[fi]);
          const bool finite = std::isfinite(factorS) && std::isfinite(solveS) &&
                              std::isfinite(discrepancy) && std::isfinite(residual) &&
                              x.allFinite();
          std::cout << cfg.name << '\t' << cfg.kind << '\t' << cfg.threshold << '\t'
                    << cfg.panel << '\t' << cfg.relax << '\t' << cfg.maxsuper << '\t'
                    << ss[si].name << '\t' << ss[si].family << '\t' << frequencies[si][fi]
                    << '\t' << aa[fi].rows() << '\t' << aa[fi].nonZeros() << '\t'
                    << factorS << '\t' << solveS << '\t' << solver.nnzL() << '\t'
                    << solver.nnzU() << '\t' << moved << '\t'
                    << static_cast<double>(moved) / aa[fi].rows() << '\t' << maxdisp << '\t'
                    << factorInfo << '\t' << solveInfo << '\t' << discrepancy << '\t'
                    << residual << '\t' << (finite ? 1 : 0) << '\t' << repeats << '\n';
        }
      }
    }
    if (sink == std::numeric_limits<double>::infinity()) std::cerr << sink;
  } catch (const std::exception &e) {
    std::cerr << "eigen_factorization_probe: " << e.what() << '\n';
    return 1;
  }
}
