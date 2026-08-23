// Focused numerical-error characterization.  The accepted physical stress
// machinery is included verbatim so this benchmark uses the same SEM source,
// receiver projections, matrices, and LAPACK/custom backsolve.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#define main eigen_pivot001_validation_reference_main
#include "eigen_pivot001_validation.cpp"
#undef main
#pragma GCC diagnostic pop

namespace {
struct TransferEntry {
  Eigen::Index row = -1, col = -1;
  C value = C(0, 0);
};

TransferEntry largestEntry(const D &y) {
  TransferEntry out;
  double best = -1;
  for (Eigen::Index i = 0; i < y.rows(); ++i)
    for (Eigen::Index j = 0; j < y.cols(); ++j)
      if (std::abs(y(i, j)) > best) {
        best = std::abs(y(i, j)); out.row = i; out.col = j; out.value = y(i, j);
      }
  return out;
}

std::string responseClass(int off, double value, double left, double right) {
  if (off == 0) return "resonance";
  if (value <= left && value <= right) return "antiresonance";
  if (std::abs(off) <= 2) return "flank";
  return "elsewhere";
}
}

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "usage: eigen_pivot_error_characterization MODEL PARAMS\n";
    return 2;
  }
  try {
    prem_norm<double> norm;
    auto model = EarthModels::ModelInput(argv[1], norm);
    InputParameters params(argv[2]);
    SourceInfo::EarthquakeCMT cmt(params);
    SpectraSolver::FreqFull productionGrid(
        .1, 55.0, 1.0, 4.0, 1.0, .05, 0.0, 4.0, 1, model.TimeNorm());
    const double gridMhz = productionGrid.df() * 1000.0 / model.TimeNorm();
    const double productionEpsilon = productionGrid.ep();
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

    auto rhs = [&](const Case &c) {
      if (c.family == "radial") return c.sem->calculateForceR(cmt);
      if (c.family == "toroidal") return c.sem->calculateForceAllT(cmt, c.degree);
      return c.sem->calculateForceAll(cmt, c.degree);
    };
    const std::vector<double> eps{productionEpsilon, productionEpsilon / 10.0,
                                  productionEpsilon / 100.0};
    const std::vector<bool> attenuation{true, false};
    // The two extra +/-2 stencils target the largest-output radial resonance
    // and the prior worst threshold-.01 transfer-error spheroidal resonance.
    const std::vector<bool> sharp{true, false, false, true, false, false};
    std::vector<double> centers;
    std::cout << std::setprecision(17)
      << "family\tname\tdegree\tmode_mhz\tfreq_mhz\tgrid_offset\tstencil_member\tresponse_class\tepsilon\tattenuation\tthreshold\tn\tnnz_a\tmax_abs_a\tfactor_info\tsolve_info\tfinite\tmin_abs_diag_u\tmax_abs_diag_u\tmax_abs_u\tmax_abs_l\tgrowth_u_over_a\tnonidentity_pivots\tpivot_fraction\tmax_pivot_displacement\tresidual\tsolution_vs_threshold1\tsolution_abs_vs_threshold1\ttransfer_norm\ttransfer_abs_error_vs_threshold1\ttransfer_relative_error_vs_threshold1\tlocal_response_scale\tlocal_response_normalized_error\ttransfer_max_row\ttransfer_max_col\ttransfer_re\ttransfer_im\tthreshold1_transfer_re\tthreshold1_transfer_im\tlapack_transfer_re\tlapack_transfer_im\ttransfer_component_abs_error_vs_threshold1\tlapack_factor_info\tlapack_solve_info\tlapack_finite\tlapack_residual\tlapack_solution_vs_threshold1\tlapack_solution_abs_vs_threshold1\tlapack_transfer_abs_error_vs_threshold1\tlapack_transfer_relative_error_vs_threshold1\tsolution_vs_lapack\ttransfer_abs_error_vs_lapack\ttransfer_relative_error_vs_lapack\n";

    // Select each physical response peak with the accepted +/-4 production
    // grid scan.  This scan is reported on stderr for the machine record.
    for (const auto &c : cs) {
      const double catalogGrid = std::round(c.mode_mhz / gridMhz) * gridMhz;
      const D receiver = receiverProjection(c, params);
      const D b = rhs(c);
      double peak = catalogGrid, peakValue = -1;
      for (int off = -4; off <= 4; ++off) {
        S a = matrix(c, catalogGrid + off * gridMhz, productionEpsilon, true,
                     model.TREF(), model.TimeNorm());
        LU scan; scan.threshold(1.0); scan.analyzePattern(a); scan.factorize(a);
        if (scan.info() != Eigen::Success) continue;
        const double value = inf(receiver * scan.solve(b));
        if (value > peakValue) { peak = catalogGrid + off * gridMhz; peakValue = value; }
      }
      centers.push_back(peak);
      std::cerr << "physical_peak\t" << c.name << "\t" << c.mode_mhz << "\t"
                << catalogGrid << "\t" << peak << "\t" << peakValue << '\n';
    }

    volatile double sink = 0;
    for (std::size_t ci = 0; ci < cs.size(); ++ci)
      for (double ep : eps) for (bool att : attenuation) {
        const auto &c = cs[ci];
        const D receiver = receiverProjection(c, params);
        const D b = rhs(c);
        // Baseline responses at +/-2 are computed for the local scale for all
        // modes, while only the two marked sharp systems emit those rows.
        std::vector<D> localY(5);
        std::vector<double> localNorm(5, 0);
        for (int oi = -2; oi <= 2; ++oi) {
          S a = matrix(c, centers[ci] + oi * gridMhz, ep, att,
                       model.TREF(), model.TimeNorm());
          LU one; one.threshold(1.0); one.analyzePattern(a); one.factorize(a);
          localY[oi + 2] = receiver * one.solve(b);
          localNorm[oi + 2] = inf(localY[oi + 2]);
        }
        const double localScale = *std::max_element(localNorm.begin(), localNorm.end());
        const std::vector<int> offsets = sharp[ci] ? std::vector<int>{-2,-1,0,1,2}
                                                    : std::vector<int>{-1,0,1};
        for (int off : offsets) {
          const int oi = off + 2;
          const double f = centers[ci] + off * gridMhz;
          const double value = localNorm[oi];
          const double left = localNorm[std::max(0, oi - 1)];
          const double right = localNorm[std::min(4, oi + 1)];
          const std::string klass = responseClass(off, value, left, right);
          S a = matrix(c, f, ep, att, model.TREF(), model.TimeNorm());
          const D xref = [&] { LU one; one.threshold(1.0); one.analyzePattern(a); one.factorize(a); D out = one.solve(b); return out; }();
          const D yref = receiver * xref;
          SPARSESPEC::detail::LapackBandSolver lapack;
          lapack.compute(a); const int lfi = static_cast<int>(lapack.info());
          const D xl = lapack.solve(b); const int lsi = static_cast<int>(lapack.info());
          const D yl = receiver * xl;
          const bool lf = xl.allFinite() && yl.allFinite() &&
                          lapack.band().coefficients().array().isFinite().all();
          // Use one physical receiver/output quantity for every backend: the
          // threshold-1 largest entry index.  Selecting a separate maximum
          // in each solve would compare different receiver/RHS components.
          const TransferEntry refEntry = largestEntry(yref);
          for (double threshold : {1.0, 0.1, 0.01}) {
            LU solver; solver.threshold(threshold); solver.analyzePattern(a); solver.factorize(a);
            const int fi = static_cast<int>(solver.info());
            const D x = solver.solve(b); const int si = static_cast<int>(solver.info());
            const D y = receiver * x;
            const Growth g = growth(solver);
            const PivotStats ps = pivotStats(solver);
            const C candidateTransfer = y(refEntry.row, refEntry.col);
            const C lapackTransfer = yl(refEntry.row, refEntry.col);
            const double transferErr = inf(y - yref);
            const bool finite = x.allFinite() && y.allFinite() &&
              std::isfinite(residual(a, x, b)) && std::isfinite(g.min_diag_u) &&
              std::isfinite(g.max_diag_u) && std::isfinite(g.max_u) &&
              std::isfinite(g.max_l) && lf && lfi == 0 && lsi == 0;
            const double componentErr = std::abs(candidateTransfer - refEntry.value);
            std::cout << c.family << '\t' << c.name << '\t' << c.degree << '\t' << c.mode_mhz << '\t'
              << f << '\t' << off << '\t' << (sharp[ci] ? 1 : 0) << '\t' << klass << '\t' << ep << '\t'
              << (att ? 1 : 0) << '\t' << threshold << '\t' << a.rows() << '\t' << a.nonZeros() << '\t'
              << maxAbs(a) << '\t' << fi << '\t' << si << '\t' << (finite ? 1 : 0) << '\t'
              << g.min_diag_u << '\t' << g.max_diag_u << '\t' << g.max_u << '\t' << g.max_l << '\t'
              << g.max_u / std::max(maxAbs(a), 1.e-300) << '\t' << ps.nonidentity << '\t'
              << static_cast<double>(ps.nonidentity) / a.rows() << '\t' << ps.max_displacement << '\t'
              << residual(a, x, b) << '\t' << rel(x, xref) << '\t' << inf(x - xref) << '\t'
              << inf(y) << '\t' << transferErr << '\t' << transferErr / std::max(inf(yref), 1.e-300) << '\t'
              << localScale << '\t' << transferErr / std::max(localScale, 1.e-300) << '\t'
              << refEntry.row << '\t' << refEntry.col << '\t' << candidateTransfer.real() << '\t' << candidateTransfer.imag() << '\t'
              << refEntry.value.real() << '\t' << refEntry.value.imag() << '\t'
              << lapackTransfer.real() << '\t' << lapackTransfer.imag() << '\t' << componentErr << '\t'
              << lfi << '\t' << lsi << '\t' << (lf ? 1 : 0) << '\t' << residual(a, xl, b) << '\t'
              << rel(xl, xref) << '\t' << inf(xl - xref) << '\t' << inf(yl - yref) << '\t'
              << inf(yl - yref) / std::max(inf(yref), 1.e-300) << '\t' << rel(x, xl) << '\t'
              << inf(y - yl) << '\t' << inf(y - yl) / std::max(inf(yl), 1.e-300) << '\n';
            sink += inf(y) + localScale;
          }
        }
      }
    if (sink == std::numeric_limits<double>::infinity()) std::cerr << sink;
  } catch (const std::exception &e) {
    std::cerr << "eigen_pivot_error_characterization: " << e.what() << '\n';
    return 1;
  }
}
