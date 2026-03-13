#ifndef FULL_SPEC_SINGLE_SEM_H
#define FULL_SPEC_SINGLE_SEM_H

#include "FullSpec.h"

namespace SPARSESPEC {

template <class model1d>
auto
Sparse_F_Spec::Spectra(SpectraSolver::FreqFull &myff, Full1D::specsem &sem,
                       model1d &inp_model, SourceInfo::EarthquakeCMT &cmt,
                       InputParameters &params, int nskip) {
  using Complex = std::complex<double>;
  using MATRIX = Eigen::MatrixXcd;
  using SMATRIX = Eigen::SparseMatrix<Complex>;
  using SLU = Eigen::SparseLU<SMATRIX, Eigen::COLAMDOrdering<int>>;

  auto vec_w = myff.w();
  SpecConstants sc(myff.ep(), inp_model.TREF());
  auto myi = sc.myi;
  auto ieps = sc.ieps;
  auto w0 = sc.w0;
  auto twodivpi = sc.twodivpi;

  MATRIX vec_raw = MATRIX::Zero(3 * params.num_receivers(), vec_w.size());
  SLU solver, solver1;

  int lmin = params.lmin();
  int lmax = params.lmax();
  ParamInfo param_info(params, lmax);
  auto NQ = sem.mesh().NN();

  auto mtype = params.type();
  ModeFlags flags = resolveModeFlags(mtype, lmin, lmax);
  bool inc_rad = flags.inc_rad;
  bool inc_tor = flags.inc_tor;
  bool inc_sph = flags.inc_sph;
  lmin = std::max(lmin, 1);

  auto num_rec = params.num_receivers();
  auto idx_source = sem.Source_Element(cmt);

  // radials
  if (inc_rad) {
    auto rec_elems = sem.Receiver_Elements(params);
    auto lowidx = sem.LtG_R(0, rec_elems[0], 0);
    auto upidx = sem.LtG_R(1, rec_elems.back(), NQ - 1);
    int lenidx = upidx - lowidx + 1;
    SMATRIX ke_r = sem.H_R().cast<Complex>();
    SMATRIX in_r = sem.P_R().cast<Complex>();
    SMATRIX ke_r_atten = sem.H_RA().cast<Complex>();
    ke_r.makeCompressed();
    in_r.makeCompressed();
    ke_r_atten.makeCompressed();
    MATRIX f_r = sem.CalculateForce_R(cmt);
    std::vector<MATRIX> vec_RV_Z;
    for (int idxr = 0; idxr < num_rec; ++idxr)
      vec_RV_Z.push_back(sem.RV_Z_R(params, idxr).block(lowidx, 0, lenidx, 1));

#pragma omp parallel default(shared) private(solver)
    {
#pragma omp for schedule(dynamic, 10)
      for (int idx = myff.i1(); idx < myff.i2(); ++idx) {
        Complex w = vec_w[idx] + ieps;
        SMATRIX w_r = -w * w * in_r + ke_r;
        if (params.attenuation())
          w_r += attenFactor(vec_w[idx], w0, twodivpi, myi) * ke_r_atten;
        w_r.makeCompressed();
        solver.compute(w_r);
        MATRIX vec_x = solver.solve(f_r);
        for (int idxr = 0; idxr < num_rec; ++idxr) {
#pragma omp critical(torvecadd)
          {
            vec_raw(3 * idxr, idx) +=
                vec_RV_Z[idxr]
                    .cwiseProduct(vec_x.block(lowidx, 0, lenidx, 1))
                    .sum();
          }
        }
      }
    }
  }

  // toroidals
  if (inc_tor) {
    auto rec_elems = sem.Receiver_Elements(params);
    auto lowidx = sem.LtG_T(rec_elems[0], 0);
    auto upidx = sem.LtG_T(rec_elems.back(), NQ - 1) + 1;
    int lenidx = upidx - lowidx;
    auto lentor = sem.LtG_T(sem.mesh().NE() - 1, NQ - 1) + 1;
#pragma omp parallel default(shared) private(solver1)
    {
#pragma omp for schedule(dynamic)
      for (int idxl = lmin; idxl < lmax + 1; ++idxl) {
        SMATRIX H_tor = sem.H_TK(idxl).cast<Complex>();
        SMATRIX P_tor = sem.P_TK(idxl).cast<Complex>();
        SMATRIX H_tor_atten = sem.H_TA(idxl).cast<Complex>();
        H_tor.makeCompressed();
        P_tor.makeCompressed();
        H_tor_atten.makeCompressed();
        MATRIX RV_VALS = param_info.RV_FULL_TOR(idxl);
        MATRIX F_VALS = sem.CalculateForce_Coefficients_T(cmt, idxl);
        MATRIX RED_C = RV_VALS * F_VALS;
        MATRIX F_BASE = sem.CalculateForce_All_T(cmt, idxl);
        MATRIX RV_BASE = sem.RV_BASE_FULL_T(params, idxl);
        auto ridxtor =
            SpectralTools::AllIndices_TOR(sem, idxl, myff, idx_source, nskip);
        MATRIX vec_raw_l =
            MATRIX::Zero(3 * params.num_receivers(), vec_w.size());
        for (int idx = myff.i2() - 1; idx > myff.i1() - 1; --idx) {
          std::size_t ridx = ridxtor[idx - myff.i1()];
          std::size_t len_ms = lentor - ridx;
          Complex w = vec_w[idx] + ieps;
          SMATRIX mat_w_tor_red =
              H_tor.block(ridx, ridx, len_ms, len_ms) -
              w * w * P_tor.block(ridx, ridx, len_ms, len_ms);
          if (params.attenuation())
            mat_w_tor_red += attenFactor(vec_w[idx], w0, twodivpi, myi) *
                             H_tor_atten.block(ridx, ridx, len_ms, len_ms);
          mat_w_tor_red.makeCompressed();
          auto f_red = F_BASE.block(ridx, 0, len_ms, F_BASE.cols());
          factorizeOrCompute(solver1, mat_w_tor_red, myff.i2() - idx - 1,
                             nskip);
          MATRIX vec_sol = solver1.solve(f_red);
          auto lidx = lowidx - ridx;
          vec_raw_l.col(idx) +=
              RED_C.cwiseProduct(RV_BASE * vec_sol.block(lidx, 0, lenidx, 2))
                  .rowwise()
                  .sum();
        }
#pragma omp critical(torvecadd)
        {
          vec_raw += vec_raw_l;
        }
      }
    }
  }

  // spheroidals
  if (inc_sph) {
    auto rec_elems = sem.Receiver_Elements(params);
    auto lowidx = sem.LtG_S(0, rec_elems[0], 0);
    auto upidx = sem.LtG_S(1, rec_elems.back(), NQ - 1);
    int lenidx = upidx - lowidx + 1;
    auto lensph = sem.LtG_S(2, sem.mesh().NE() - 1, NQ - 1) + 1;
#pragma omp parallel default(shared) private(solver1)
    {
#pragma omp for schedule(dynamic)
      for (int idxl = lmin; idxl < lmax + 1; ++idxl) {
        SMATRIX h_s = sem.H_S(idxl).cast<Complex>();
        SMATRIX p_s = sem.P_S(idxl).cast<Complex>();
        SMATRIX h_sa = sem.H_SA(idxl).cast<Complex>();
        h_s.makeCompressed();
        p_s.makeCompressed();
        h_sa.makeCompressed();
        MATRIX RV_VALS = param_info.RV_FULL_SPH(idxl);
        MATRIX F_VALS = sem.CalculateForce_Coefficients(cmt, idxl);
        MATRIX RED_C = RV_VALS * F_VALS;
        MATRIX F_BASE = sem.CalculateForce_All(cmt, idxl);
        MATRIX RV_BASE = sem.RV_BASE_FULL(params, idxl);
        auto vec_ridx =
            SpectralTools::AllIndices_SPH(sem, idxl, myff, idx_source, nskip);
        MATRIX vec_raw_l =
            MATRIX::Zero(3 * params.num_receivers(), vec_w.size());
        for (int idx = myff.i2() - 1; idx > myff.i1() - 1; --idx) {
          std::size_t idx_rs = vec_ridx[idx - myff.i1()];
          std::size_t len_ms = lensph - idx_rs;
          Complex w = vec_w[idx] + ieps;
          SMATRIX mat_sph = h_s.block(idx_rs, idx_rs, len_ms, len_ms) -
                            w * w * p_s.block(idx_rs, idx_rs, len_ms, len_ms);
          if (params.attenuation())
            mat_sph += attenFactor(vec_w[idx], w0, twodivpi, myi) *
                       h_sa.block(idx_rs, idx_rs, len_ms, len_ms);
          mat_sph.makeCompressed();
          auto f_red = F_BASE.block(idx_rs, 0, len_ms, F_BASE.cols());
          factorizeOrCompute(solver1, mat_sph, myff.i2() - idx - 1, nskip);
          MATRIX vec_sol = solver1.solve(f_red);
          auto lidx = lowidx - idx_rs;
          vec_raw_l.col(idx) +=
              RED_C.cwiseProduct(RV_BASE * vec_sol.block(lidx, 0, lenidx, 4))
                  .rowwise()
                  .sum();
        }
#pragma omp critical(sphvecadd)
        {
          vec_raw += vec_raw_l;
        }
      }
    }
  }

  for (int idx = 0; idx < vec_raw.cols(); ++idx) {
    Complex w = vec_w[idx] + ieps;
    vec_raw.col(idx) *= outputFactor(params.output_type(), w, myi);
  }
  return vec_raw;
}

}   // namespace SPARSESPEC
#endif   // FULL_SPEC_SINGLE_SEM_H
