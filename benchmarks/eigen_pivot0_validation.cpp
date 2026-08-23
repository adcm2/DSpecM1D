#include "config.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include <Eigen/SparseLU>

#include <DSpecM1D/ModelInput>
#include <DSpecM1D/src/FullSpecMultiSem.h>
#include <DSpecM1D/src/LapackBandSolver.h>
#include <DSpecM1D/src/SEM/SEM.h>
#include <DSpecM1D/src/SpecHelpers.h>
#include <DSpecM1D/FrequencyTools>

// Focused numerical validation only.  This executable intentionally does not
// alter the production SparseLU type or its default threshold.
namespace {
using C = std::complex<double>;
using S = Eigen::SparseMatrix<C>;
using D = Eigen::MatrixXcd;

struct LU : Eigen::SparseLU<S, Eigen::NaturalOrdering<int>> {
  using Base = Eigen::SparseLU<S, Eigen::NaturalOrdering<int>>;
  void threshold(double x) { setPivotThreshold(x); }
};

struct Case {
  std::string family, name;
  int degree = 0;
  double mode_mhz = 0;
  Full1D::SEM *sem = nullptr;
  S h, p, ha;
};

double inf(const D &x) {
  double v = 0;
  for (Eigen::Index i = 0; i < x.rows(); ++i) {
    double row = 0;
    for (Eigen::Index j = 0; j < x.cols(); ++j) row += std::abs(x(i, j));
    v = std::max(v, row);
  }
  return v;
}
double inf(const S &a) {
  Eigen::VectorXd rows = Eigen::VectorXd::Zero(a.rows());
  for (Eigen::Index j = 0; j < a.outerSize(); ++j)
    for (S::InnerIterator it(a, j); it; ++it) rows(it.row()) += std::abs(it.value());
  return rows.size() ? rows.maxCoeff() : 0;
}
double residual(const S &a, const D &x, const D &b) {
  return inf(a * x - b) / (inf(a) * inf(x) + inf(b));
}
double rel(const D &a, const D &b) {
  return inf(a - b) / std::max(inf(b), 1.e-300);
}
double maxAbs(const S &a) {
  double v = 0;
  for (Eigen::Index j = 0; j < a.outerSize(); ++j)
    for (S::InnerIterator it(a, j); it; ++it) v = std::max(v, std::abs(it.value()));
  return v;
}

struct Growth {
  double min_diag_u = std::numeric_limits<double>::infinity();
  double max_diag_u = 0, max_u = 0, max_l = 0;
};
Growth growth(const LU &s) {
  Growth g;
  auto l = s.matrixL();
  auto u = s.matrixU();
  for (Eigen::Index j = 0; j < u.cols(); ++j)
    for (typename std::decay_t<decltype(u.m_mapU)>::InnerIterator it(u.m_mapU, j); it; ++it) {
      g.max_u = std::max(g.max_u, std::abs(it.value()));
    }
  // L's mapped supernodal storage contains both the dense U diagonal/upper
  // block and the rectangular L rows. Classify by dense block indices, not by
  // the mapped row number (which is why a plain InnerIterator over L is wrong).
  for (Eigen::Index sn = 0; sn <= l.m_mapL.nsuper(); ++sn) {
    const Eigen::Index first = l.m_mapL.supToCol()[sn];
    const Eigen::Index width = l.m_mapL.supToCol()[sn + 1] - first;
    const Eigen::Index rowBegin = l.m_mapL.rowIndexPtr()[first];
    for (Eigen::Index jc = 0; jc < width; ++jc) {
      const Eigen::Index col = first + jc;
      const Eigen::Index begin = l.m_mapL.colIndexPtr()[col];
      const Eigen::Index end = l.m_mapL.colIndexPtr()[col + 1];
      for (Eigen::Index ir = 0; ir < end - begin; ++ir) {
        const double value = std::abs(l.m_mapL.valuePtr()[begin + ir]);
        if (ir < width && ir <= jc) {
          g.max_u = std::max(g.max_u, value);
          if (ir == jc) {
            g.min_diag_u = std::min(g.min_diag_u, value);
            g.max_diag_u = std::max(g.max_diag_u, value);
          }
        } else {
          g.max_l = std::max(g.max_l, value);
        }
      }
      (void)rowBegin;
    }
  }
  if (!std::isfinite(g.min_diag_u)) g.min_diag_u = 0;
  return g;
}

S matrix(const Case &c, double f_mhz, double ep, bool atten, double tref,
         double timeNorm) {
  SPARSESPEC::SpecConstants sc(ep, tref);
  const double wreal = 2.0 * EIGEN_PI * (f_mhz / 1000.0) * timeNorm;
  C w = wreal + sc.ieps;
  S a = c.h - w * w * c.p;
  if (atten) a += SPARSESPEC::attenFactor(wreal, sc.w0, sc.twodivpi, sc.myi) * c.ha;
  a.makeCompressed();
  return a;
}

D receiverProjection(const Case &c, InputParameters &params) {
  D out;
  if (c.family == "radial") {
    out.resize(params.num_receivers(), c.h.rows()); out.setZero();
    for (int idx = 0; idx < params.num_receivers(); ++idx)
      out.row(idx) = c.sem->rvZR(params, idx).col(0).transpose();
  } else if (c.family == "toroidal") {
    auto elems = c.sem->receiverElements(params);
    auto low = c.sem->ltgT(elems[0], 0);
    auto v = c.sem->rvBaseFullT(params, c.degree);
    out.resize(v.rows(), c.h.rows()); out.setZero();
    out.block(0, low, v.rows(), v.cols()) = v;
  } else {
    auto elems = c.sem->receiverElements(params);
    auto low = c.sem->ltgS(0, elems[0], 0);
    auto v = c.sem->rvBaseFull(params, c.degree);
    out.resize(v.rows(), c.h.rows()); out.setZero();
    out.block(0, low, v.rows(), v.cols()) = v;
  }
  return out;
}

void add(std::vector<Case> &out, const std::string &family, const std::string &name,
         int degree, double mode, Full1D::SEM &sem, const S &h, const S &p, const S &ha) {
  Case c{family, name, degree, mode, &sem, h, p, ha};
  c.h.makeCompressed(); c.p.makeCompressed(); c.ha.makeCompressed();
  out.push_back(std::move(c));
}

int movedRows(const LU &s) {
  int n = 0;
  auto p = s.rowsPermutation().indices();
  for (Eigen::Index i = 0; i < p.size(); ++i) if (p(i) != i) ++n;
  return n;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "usage: eigen_pivot0_validation MODEL PARAMS\n";
    return 2;
  }
  try {
    prem_norm<double> norm;
    auto model = EarthModels::ModelInput(argv[1], norm);
    InputParameters params(argv[2]);
    SourceInfo::EarthquakeCMT cmt(params);

    // The frequencies are the n=0 MINEOS catalogue modes for PREM: radial
    // 0S0, toroidal 0T2/0T12, spheroidal 0S2/0S12/0S40.  Three points around
    // each root are the nearest ex3 production-grid frequency and its two
    // adjacent grid points. FreqFull stores angular normalized frequencies;
    // catalogue values are physical mHz and are converted below.
    SpectraSolver::FreqFull productionGrid(
        .1, 55.0, 1.0, 4.0, 1.0, .05, 0.0, 4.0, 1, model.TimeNorm());
    const double gridMhz = productionGrid.df() * 1000.0 / model.TimeNorm();
    Full1D::SEM radial(model, .05, 6, 1);
    Full1D::SEM tor(model, .05, 6, 12);
    Full1D::SEM sph(model, .05, 6, 12);
    Full1D::SEM fluid(model, .01, 6, 40);
    std::vector<Case> cs;
    add(cs, "radial", "radial_0S0", 0, .8146635, radial,
        radial.hR().cast<C>().eval(), radial.pR().cast<C>().eval(), radial.hRa().cast<C>().eval());
    add(cs, "toroidal", "toroidal_0T2", 2, .3830516, tor,
        tor.hTk(2).cast<C>().eval(), tor.pTk(2).cast<C>().eval(), tor.hTa(2).cast<C>().eval());
    add(cs, "toroidal", "toroidal_0T12", 12, 1.881496, tor,
        tor.hTk(12).cast<C>().eval(), tor.pTk(12).cast<C>().eval(), tor.hTa(12).cast<C>().eval());
    add(cs, "spheroidal", "spheroidal_0S2", 2, .3108622, sph,
        sph.hS(2).cast<C>().eval(), sph.pS(2).cast<C>().eval(), sph.hSa(2).cast<C>().eval());
    add(cs, "spheroidal", "spheroidal_0S12", 12, 2.003070, sph,
        sph.hS(12).cast<C>().eval(), sph.pS(12).cast<C>().eval(), sph.hSa(12).cast<C>().eval());
    add(cs, "fluid-solid", "fluid_solid_0S40", 40, 4.761725, fluid,
        fluid.hS(40).cast<C>().eval(), fluid.pS(40).cast<C>().eval(), fluid.hSa(40).cast<C>().eval());

    // A small threshold-1, production-epsilon scan checks that the three
    // tested grid points are anchored to the physical catalogue peak. It is
    // deliberately limited to nine grid points around each catalogue root.
    auto sourceRhs = [&](const Case &c) {
      if (c.family == "radial") return c.sem->calculateForceR(cmt);
      if (c.family == "toroidal") return c.sem->calculateForceAllT(cmt, c.degree);
      return c.sem->calculateForceAll(cmt, c.degree);
    };
    std::vector<double> selectedCenters;
    selectedCenters.reserve(cs.size());
    const double productionEpsilon = productionGrid.ep();
    for (const auto &c : cs) {
      const double catalogGrid = std::round(c.mode_mhz / gridMhz) * gridMhz;
      const D receiver = receiverProjection(c, params);
      const D b = sourceRhs(c);
      double peak = catalogGrid, peakValue = -1;
      for (int off = -4; off <= 4; ++off) {
        const double f = catalogGrid + off * gridMhz;
        S a = matrix(c, f, productionEpsilon, true, model.TREF(), model.TimeNorm());
        LU scan; scan.threshold(1.0); scan.analyzePattern(a); scan.factorize(a);
        if (scan.info() != Eigen::Success) continue;
        const double value = inf(receiver * scan.solve(b));
        if (value > peakValue) { peak = f; peakValue = value; }
      }
      selectedCenters.push_back(peak);
      std::cerr << "physical_peak\t" << c.name << "\t" << c.mode_mhz << "\t"
                << catalogGrid << "\t" << peak << "\t" << peakValue << '\n';
    }

    std::cout << std::setprecision(17)
      << "family\tname\tdegree\tmode_mhz\tfreq_mhz\tepsilon\tattenuation\tthreshold\tn\tnnz_a\tmax_abs_a\tfactor_info\tsolve_info\tfinite\tmin_abs_diag_u\tmax_abs_diag_u\tmax_abs_u\tmax_abs_l\tgrowth_u_over_a\tnonidentity_pivots\tresidual\tsolution_vs_threshold1\tsolution_abs_vs_threshold1\ttransfer\ttransfer_vs_threshold1\ttransfer_relative_vs_threshold1\tlapack_factor_info\tlapack_solve_info\tlapack_finite\tlapack_residual\tlapack_solution_vs_threshold1\tlapack_solution_abs_vs_threshold1\tlapack_transfer_vs_threshold1\tlapack_transfer_relative_vs_threshold1\tsolution_vs_lapack\ttransfer_vs_lapack\ttransfer_relative_vs_lapack\n";

    const std::vector<double> eps{productionEpsilon, productionEpsilon / 10.0,
                                  productionEpsilon / 100.0};
    const std::vector<bool> attenuation{true, false};
    const std::vector<int> gridOffsets{-1, 0, 1};
    volatile double sink = 0;
    for (std::size_t ci = 0; ci < cs.size(); ++ci) for (double ep : eps) for (bool att : attenuation) {
      const auto &c = cs[ci];
      D receiver = receiverProjection(c, params);
      for (int off : gridOffsets) {
        const double f = selectedCenters[ci] + off * gridMhz;
        S a = matrix(c, f, ep, att, model.TREF(), model.TimeNorm());
        D b;
        if (c.family == "radial") b = c.sem->calculateForceR(cmt);
        else if (c.family == "toroidal") b = c.sem->calculateForceAllT(cmt, c.degree);
        else b = c.sem->calculateForceAll(cmt, c.degree);
        // Both threshold settings use exactly the same source-like RHS.
        LU ref; ref.threshold(1.0); ref.analyzePattern(a); ref.factorize(a);
        D xref = ref.solve(b);
        D yref = receiver * xref;
        const double refResidual = residual(a, xref, b);
        (void)refResidual;
        LU zero; zero.threshold(0.0); zero.analyzePattern(a);
        zero.factorize(a); const int zeroFactorInfo = static_cast<int>(zero.info());
        D x = zero.solve(b); const int zeroSolveInfo = static_cast<int>(zero.info());
        D y = receiver * x;
        Growth g = growth(zero);
        const auto bw = SPARSESPEC::detail::lapackBandBandwidth(a);
        SPARSESPEC::detail::LapackBandSolver lapack;
        lapack.compute(a); const int lapackFactorInfo = static_cast<int>(lapack.info());
        D xl = lapack.solve(b); const int lapackSolveInfo = static_cast<int>(lapack.info());
        D yl = receiver * xl;
        const bool lapackFinite = xl.allFinite() && yl.allFinite() &&
                                  lapack.band().coefficients().array().isFinite().all();
        const bool finite = x.allFinite() && y.allFinite() && std::isfinite(maxAbs(a)) &&
                            std::isfinite(residual(a, x, b)) &&
                            std::isfinite(g.min_diag_u) && std::isfinite(g.max_u) &&
                            std::isfinite(g.max_l) && std::isfinite(g.max_diag_u) && lapackFinite &&
                            lapackFactorInfo == 0 && lapackSolveInfo == 0;
        const double transfer = inf(y);
        const double transferRef = inf(yref);
        std::cout << c.family << '\t' << c.name << '\t' << c.degree << '\t' << c.mode_mhz << '\t'
                  << f << '\t' << ep << '\t' << (att ? 1 : 0) << "\t0\t" << a.rows() << '\t'
                  << a.nonZeros() << '\t' << maxAbs(a) << '\t' << zeroFactorInfo << '\t'
                  << zeroSolveInfo << '\t' << (finite ? 1 : 0) << '\t'
                  << g.min_diag_u << '\t' << g.max_diag_u << '\t' << g.max_u << '\t' << g.max_l << '\t'
                  << g.max_u / std::max(maxAbs(a), 1.e-300) << '\t' << movedRows(zero) << '\t'
                  << residual(a, x, b) << '\t' << rel(x, xref) << '\t' << inf(x - xref) << '\t' << transfer << '\t'
                  << inf(y - yref) << '\t'
                  << inf(y - yref) / std::max(transferRef, 1.e-300) << '\t'
                  << lapackFactorInfo << '\t' << lapackSolveInfo << '\t' << (lapackFinite ? 1 : 0) << '\t' << residual(a, xl, b) << '\t'
                  << rel(xl, xref) << '\t' << inf(xl - xref) << '\t' << inf(yl - yref) << '\t'
                  << inf(yl - yref) / std::max(transferRef, 1.e-300) << '\t'
                  << rel(x, xl) << '\t' << inf(y - yl) << '\t'
                  << inf(y - yl) / std::max(inf(yl), 1.e-300) << '\n';

        // Emit a second record for threshold 1.  Keeping the same schema makes
        // direct threshold-0 versus threshold-1 factor diagnostics auditable.
        LU one; one.threshold(1.0); one.analyzePattern(a); one.factorize(a);
        const int oneFactorInfo = static_cast<int>(one.info());
        D xo = one.solve(b); const int oneSolveInfo = static_cast<int>(one.info());
        Growth go = growth(one);
        D yo = receiver * xo;
        const bool oneFinite = xo.allFinite() && yo.allFinite() && std::isfinite(maxAbs(a)) &&
                               std::isfinite(residual(a, xo, b)) && std::isfinite(go.min_diag_u) &&
                               std::isfinite(go.max_diag_u) && std::isfinite(go.max_u) &&
                               std::isfinite(go.max_l);
        std::cout << c.family << '\t' << c.name << '\t' << c.degree << '\t' << c.mode_mhz << '\t'
                  << f << '\t' << ep << '\t' << (att ? 1 : 0) << "\t1\t" << a.rows() << '\t'
                  << a.nonZeros() << '\t' << maxAbs(a) << '\t' << oneFactorInfo << '\t'
                  << oneSolveInfo << '\t' << (oneFinite ? 1 : 0) << '\t'
                  << go.min_diag_u << '\t' << go.max_diag_u << '\t' << go.max_u << '\t' << go.max_l << '\t'
                  << go.max_u / std::max(maxAbs(a), 1.e-300) << '\t' << movedRows(one) << '\t'
                  << residual(a, xo, b) << "\t0\t0\t" << inf(yo) << "\t0\t0\t"
                  << lapackFactorInfo << '\t' << lapackSolveInfo << '\t' << (lapackFinite ? 1 : 0) << '\t' << residual(a, xl, b) << '\t'
                  << rel(xl, xo) << '\t' << inf(xl - xo) << '\t'
                  << inf(yl - yo) << '\t'
                  << inf(yl - yo) / std::max(inf(yo), 1.e-300) << '\t'
                  << rel(xo, xl) << '\t'
                  << inf(yo - yl) << '\t'
                  << inf(yo - yl) / std::max(inf(yl), 1.e-300) << '\n';
        sink += transfer + transferRef + static_cast<double>(bw.first + bw.second);
      }
    }
    if (sink == std::numeric_limits<double>::infinity()) std::cerr << sink;
  } catch (const std::exception &e) {
    std::cerr << "eigen_pivot0_validation: " << e.what() << '\n';
    return 1;
  }
}
