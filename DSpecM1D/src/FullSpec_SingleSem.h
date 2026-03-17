#ifndef FULL_SPEC_SINGLE_SEM_H
#define FULL_SPEC_SINGLE_SEM_H

#include "FullSpec.h"

namespace SPARSESPEC {

template <class model1d>
auto
SparseFSpec::spectra(SpectraSolver::FreqFull &myff, Full1D::SEM &sem,
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
  ParamInfo paramInfo(params, lmax);
  auto NQ = sem.mesh().NN();

  auto mtype = params.type();
  ModeFlags flags = resolveModeFlags(mtype, lmin, lmax);
  bool inc_rad = flags.inc_rad;
  bool inc_tor = flags.inc_tor;
  bool inc_sph = flags.inc_sph;
  lmin = std::max(lmin, 1);

  auto numRec = params.num_receivers();
  auto idxSource = sem.sourceElement(cmt);

  // radials
  if (inc_rad) {
    auto recElems = sem.receiverElements(params);
    auto lowidx = sem.ltgR(0, recElems[0], 0);
    auto upidx = sem.ltgR(1, recElems.back(), NQ - 1);
    int lenidx = upidx - lowidx + 1;
    SMATRIX ke_r = sem.hR().cast<Complex>();
    SMATRIX in_r = sem.pR().cast<Complex>();
    SMATRIX ke_r_atten = sem.hRa().cast<Complex>();
    ke_r.makeCompressed();
    in_r.makeCompressed();
    ke_r_atten.makeCompressed();
    MATRIX f_r = sem.calculateForceR(cmt);
    std::vector<MATRIX> vecRvZ;
    for (int idxr = 0; idxr < numRec; ++idxr)
      vecRvZ.push_back(sem.rvZR(params, idxr).block(lowidx, 0, lenidx, 1));

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
        for (int idxr = 0; idxr < numRec; ++idxr) {
#pragma omp critical(torvecadd)
          {
            vec_raw(3 * idxr, idx) +=
                vecRvZ[idxr]
                    .cwiseProduct(vec_x.block(lowidx, 0, lenidx, 1))
                    .sum();
          }
        }
      }
    }
  }

  // toroidals
  if (inc_tor) {
    auto recElems = sem.receiverElements(params);
    auto lowidx = sem.ltgT(recElems[0], 0);
    auto upidx = sem.ltgT(recElems.back(), NQ - 1) + 1;
    int lenidx = upidx - lowidx;
    auto lentor = sem.ltgT(sem.mesh().NE() - 1, NQ - 1) + 1;
#pragma omp parallel default(shared) private(solver1)
    {
#pragma omp for schedule(dynamic)
      for (int idxl = lmin; idxl < lmax + 1; ++idxl) {
        SMATRIX H_tor = sem.hTk(idxl).cast<Complex>();
        SMATRIX P_tor = sem.pTk(idxl).cast<Complex>();
        SMATRIX H_tor_atten = sem.hTa(idxl).cast<Complex>();
        H_tor.makeCompressed();
        P_tor.makeCompressed();
        H_tor_atten.makeCompressed();
        MATRIX rvVals = paramInfo.rvFullTor(idxl);
        MATRIX fVals = sem.calculateForceCoefficientsT(cmt, idxl);
        MATRIX redC = rvVals * fVals;
        MATRIX fBase = sem.calculateForceAllT(cmt, idxl);
        MATRIX rvBase = sem.rvBaseFullT(params, idxl);
        auto ridxtor =
            SpectralTools::allIndicesTor(sem, idxl, myff, idxSource, nskip);
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
          auto f_red = fBase.block(ridx, 0, len_ms, fBase.cols());
          factorizeOrCompute(solver1, mat_w_tor_red, myff.i2() - idx - 1,
                             nskip);
          MATRIX vec_sol = solver1.solve(f_red);
          auto lidx = lowidx - ridx;
#pragma omp critical(torvecadd)
          {
            vec_raw.col(idx) +=
                redC.cwiseProduct(rvBase * vec_sol.block(lidx, 0, lenidx, 2))
                    .rowwise()
                    .sum();
          }
        }
      }
    }
  }

  // spheroidals
  if (inc_sph) {
    auto recElems = sem.receiverElements(params);
    auto lowidx = sem.ltgS(0, recElems[0], 0);
    auto upidx = sem.ltgS(1, recElems.back(), NQ - 1);
    int lenidx = upidx - lowidx + 1;
    auto lensph = sem.ltgS(2, sem.mesh().NE() - 1, NQ - 1) + 1;
#pragma omp parallel default(shared) private(solver1)
    {
#pragma omp for schedule(dynamic)
      for (int idxl = lmin; idxl < lmax + 1; ++idxl) {
        SMATRIX h_s = sem.hS(idxl).cast<Complex>();
        SMATRIX p_s = sem.pS(idxl).cast<Complex>();
        SMATRIX h_sa = sem.hSa(idxl).cast<Complex>();
        h_s.makeCompressed();
        p_s.makeCompressed();
        h_sa.makeCompressed();
        MATRIX rvVals = paramInfo.rvFullSph(idxl);
        MATRIX fVals = sem.calculateForceCoefficients(cmt, idxl);
        MATRIX redC = rvVals * fVals;
        MATRIX fBase = sem.calculateForceAll(cmt, idxl);
        MATRIX rvBase = sem.rvBaseFull(params, idxl);
        auto vecRidx =
            SpectralTools::allIndicesSph(sem, idxl, myff, idxSource, nskip);
        // MATRIX vec_raw_l =
        //     MATRIX::Zero(3 * params.num_receivers(), vec_w.size());
        for (int idx = myff.i2() - 1; idx > myff.i1() - 1; --idx) {
          std::size_t idxRs = vecRidx[idx - myff.i1()];
          std::size_t len_ms = lensph - idxRs;
          Complex w = vec_w[idx] + ieps;
          SMATRIX mat_sph = h_s.block(idxRs, idxRs, len_ms, len_ms) -
                            w * w * p_s.block(idxRs, idxRs, len_ms, len_ms);
          if (params.attenuation())
            mat_sph += attenFactor(vec_w[idx], w0, twodivpi, myi) *
                       h_sa.block(idxRs, idxRs, len_ms, len_ms);
          mat_sph.makeCompressed();
          auto f_red = fBase.block(idxRs, 0, len_ms, fBase.cols());
          factorizeOrCompute(solver1, mat_sph, myff.i2() - idx - 1, nskip);
          MATRIX vec_sol = solver1.solve(f_red);
          auto lidx = lowidx - idxRs;

#pragma omp critical(sphvecadd)
          {
            vec_raw.col(idx) +=
                redC.cwiseProduct(rvBase * vec_sol.block(lidx, 0, lenidx, 4))
                    .rowwise()
                    .sum();
          }
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
