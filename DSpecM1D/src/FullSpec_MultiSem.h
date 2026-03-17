#ifndef FULL_SPEC_MULTI_SEM_H
#define FULL_SPEC_MULTI_SEM_H

#include "FullSpec.h"

namespace SPARSESPEC {

template <class model1d>
auto
Sparse_F_Spec::Spectra(SpectraSolver::FreqFull &myff, model1d &inp_model,
                       SourceInfo::EarthquakeCMT &cmt, InputParameters &params,
                       int NQ, SRInfo &srinfo, double relerr) {
  using Complex = std::complex<double>;
  using MATRIX = Eigen::MatrixXcd;
  using SMATRIX = Eigen::SparseMatrix<Complex>;
  using SLU = Eigen::SparseLU<SMATRIX, Eigen::COLAMDOrdering<int>>;
  Timer timer1;

  auto vec_w = myff.w();
  SpecConstants sc(myff.ep(), inp_model.TREF());
  auto myi = sc.myi;
  auto ieps = sc.ieps;
  auto w0 = sc.w0;
  auto twodivpi = sc.twodivpi;
  auto twopid = sc.twopid;

  auto num_rec = params.num_receivers();
  int nskip = std::max(1, (myff.i2() - myff.i1()) / 20);
  auto i_begin = myff.i1();
  auto i_end = myff.i2();   // one past last valid index
  if (i_begin < 0 || i_end > static_cast<int>(vec_w.size()) ||
      i_begin >= i_end) {
    throw std::runtime_error(
        "Invalid frequency index range in FullSpec_MultiSem.");
  }

  int num_chunks = std::max(
      1, static_cast<int>(std::floor((myff.f22() - myff.f11()) / 10.0) + 1.0));

  std::vector<std::vector<double>> freq_chunks(num_chunks);
  std::vector<std::vector<int>> idx_chunks(num_chunks);

  double wl = vec_w[i_begin];
  double wh = vec_w[i_end - 1];
  double wdiff = (wh - wl) / static_cast<double>(num_chunks);

  int idxw = i_begin;
  for (int c = 0; c < num_chunks; ++c) {
    double upper = (c == num_chunks - 1) ? wh : (wl + wdiff * (c + 1));
    while (idxw < i_end && vec_w[idxw] <= upper) {
      freq_chunks[c].push_back(vec_w[idxw]);
      idx_chunks[c].push_back(idxw);
      ++idxw;
    }
  }

  // Any stragglers (numeric edge cases) go to last chunk
  while (idxw < i_end) {
    freq_chunks.back().push_back(vec_w[idxw]);
    idx_chunks.back().push_back(idxw);
    ++idxw;
  }

  // Build max_steps safely
  std::vector<double> max_steps;
  max_steps.reserve(num_chunks);
  double base_len = 0.78 * 1000 / inp_model.TimeNorm() * twopid;
  double newlen = std::pow(100.0 * relerr, 1.0 / (NQ - 1.0)) * base_len;

  for (int c = 0; c < num_chunks; ++c) {
    double fref = freq_chunks[c].empty() ? wh : freq_chunks[c].back();
    if (fref <= 0.0) {
      throw std::runtime_error("Non-positive reference frequency in chunking.");
    }
    double step = std::min(newlen / fref, 0.05);
    max_steps.push_back(step);
  }

  MATRIX vec_raw = MATRIX::Zero(3 * params.num_receivers(), vec_w.size());

  for (auto s : max_steps)
    std::cout << "Max step: " << s << "\n";

  // build one specsem per chunk
  std::vector<Full1D::specsem> sems;
  for (int idx = 0; idx < num_chunks; ++idx)
    sems.emplace_back(inp_model, max_steps[idx], NQ, params.lmax());

  timer1.start();
  SLU solver, solver1;

  int lmin = params.lmin();
  int lmax = params.lmax();
  ParamRedInfo param_info(params, lmax);

  ModeFlags flags = resolveModeFlags(params.type(), lmin, lmax);
  bool inc_rad = flags.inc_rad;
  bool inc_tor = flags.inc_tor;
  bool inc_sph = flags.inc_sph;
  lmin = std::max(lmin, 1);

  // radials
  if (inc_rad) {
    timer1.start();
    for (int idx_chunk = 0; idx_chunk < num_chunks; ++idx_chunk) {
      Full1D::specsem &sem = sems[idx_chunk];
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
      MATRIX f_r = sem.CalculateForce_Red_R(cmt);
      MATRIX vec_red_z = sem.RV_RED_Z_R(params);
#pragma omp parallel default(shared) private(solver)
      {
#pragma omp for schedule(dynamic)
        for (int idx = 0; idx < (int) idx_chunks[idx_chunk].size(); ++idx) {
          double wval = freq_chunks[idx_chunk][idx];
          int idxw = idx_chunks[idx_chunk][idx];
          Complex w = wval + ieps;
          SMATRIX w_r = -w * w * in_r + ke_r;
          if (params.attenuation())
            w_r += attenFactor(wval, w0, twodivpi, myi) * ke_r_atten;
          w_r.makeCompressed();
          solver.compute(w_r);
          MATRIX vec_x = solver.solve(f_r);
          auto tmp = vec_x.block(lowidx, 0, lenidx, 1);
          auto resval = (vec_red_z.transpose() * tmp).sum();
          for (int idxr = 0; idxr < num_rec; ++idxr) {
#pragma omp critical(torvecadd)
            {
              vec_raw(3 * idxr, idxw) += resval;
            }
          }
        }
      }
    }
    timer1.stop("Time for Radial Modes");
    std::cout << "\n";
  }

  // toroidals
  if (inc_tor) {
    std::cout << "Doing Toroidal Modes\n";
#pragma omp parallel default(shared) private(solver1)
    {
#pragma omp for schedule(dynamic)
      for (int idxl = lmin; idxl < lmax + 1; ++idxl) {
        MATRIX RV_VALS = srinfo.RV_RED_TOR(idxl);
        for (int idx_chunk = 0; idx_chunk < num_chunks; ++idx_chunk) {
          Full1D::specsem &sem = sems[idx_chunk];
          auto idx_source = sem.Source_Element(cmt);
          auto rec_elems = sem.Receiver_Elements(params);
          auto lowidx = sem.LtG_T(rec_elems[0], 0);
          auto upidx = sem.LtG_T(rec_elems.back(), NQ - 1) + 1;
          int lenidx = upidx - lowidx;
          auto lentor = sem.LtG_T(sem.mesh().NE() - 1, NQ - 1) + 1;
          SMATRIX H_tor = sem.H_TK(idxl).cast<Complex>();
          SMATRIX P_tor = sem.P_TK(idxl).cast<Complex>();
          SMATRIX H_tor_atten = sem.H_TA(idxl).cast<Complex>();
          H_tor.makeCompressed();
          P_tor.makeCompressed();
          H_tor_atten.makeCompressed();
          MATRIX F_VALS = sem.CalculateForce_RED_Coefficients_T(cmt, idxl, 0.0);
          MATRIX RED_C = RV_VALS * F_VALS;
          MATRIX F_BASE = sem.CalculateForce_All_T(cmt, idxl);
          MATRIX RV_BASE = sem.RV_BASE_FULL_T(params, idxl);
          auto ridxtor = SpectralTools::AllIndices_TOR(
              sem, idxl, freq_chunks[idx_chunk], idx_source, nskip);
          auto i1 = idx_chunks[idx_chunk][0];
          auto len_chunk = freq_chunks[idx_chunk].size();
          MATRIX vec_raw_l =
              MATRIX::Zero(3 * params.num_receivers(), len_chunk);
          for (int idx = (int) len_chunk - 1; idx > -1; --idx) {
            auto wval = freq_chunks[idx_chunk][idx];
            Complex w = wval + ieps;
            std::size_t ridx = ridxtor[idx];
            std::size_t len_ms = lentor - ridx;
            SMATRIX mat = H_tor.block(ridx, ridx, len_ms, len_ms) -
                          w * w * P_tor.block(ridx, ridx, len_ms, len_ms);
            if (params.attenuation())
              mat += attenFactor(wval, w0, twodivpi, myi) *
                     H_tor_atten.block(ridx, ridx, len_ms, len_ms);
            mat.makeCompressed();
            auto f_red = F_BASE.block(ridx, 0, len_ms, F_BASE.cols());
            factorizeOrCompute(solver1, mat, (int) len_chunk - idx - 1, nskip);
            MATRIX vec_sol = solver1.solve(f_red);
            auto lidx = lowidx - ridx;
            vec_raw_l.col(idx) +=
                RED_C.cwiseProduct(RV_BASE * vec_sol.block(lidx, 0, lenidx, 2))
                    .rowwise()
                    .sum();
          }
#pragma omp critical(torvecadd)
          {
            vec_raw.block(0, i1, 3 * params.num_receivers(), len_chunk) +=
                vec_raw_l;
          }
        }
      }
    }
  }

  // spheroidals
  if (inc_sph) {
    std::cout << "\nDoing Spheroidal Modes\n";
#pragma omp parallel default(shared) private(solver1)
    {
#pragma omp for schedule(dynamic)
      for (int idxl = lmin; idxl < lmax + 1; ++idxl) {
        MATRIX RV_VALS = srinfo.RV_RED_SPH(idxl);
        for (int idx_chunk = 0; idx_chunk < num_chunks; ++idx_chunk) {
          Full1D::specsem &sem = sems[idx_chunk];
          auto idx_source = sem.Source_Element(cmt);
          auto rec_elems = sem.Receiver_Elements(params);
          auto lowidx = sem.LtG_S(0, rec_elems[0], 0);
          auto upidx = sem.LtG_S(1, rec_elems.back(), NQ - 1);
          int lenidx = upidx - lowidx + 1;
          auto lensph = sem.LtG_S(2, sem.mesh().NE() - 1, NQ - 1) + 1;
          SMATRIX p_s = sem.P_S(idxl).cast<Complex>();
          SMATRIX h_s = sem.H_S(idxl).cast<Complex>();
          SMATRIX h_sa = sem.H_SA(idxl).cast<Complex>();
          p_s.makeCompressed();
          h_s.makeCompressed();
          h_sa.makeCompressed();
          MATRIX F_VALS = sem.CalculateForce_RED_Coefficients(cmt, idxl, 0.0);
          MATRIX RED_C = RV_VALS * F_VALS;
          MATRIX F_BASE = sem.CalculateForce_All(cmt, idxl);
          MATRIX RV_BASE = sem.RV_BASE_FULL(params, idxl);
          auto ridxsph = SpectralTools::AllIndices_SPH(
              sem, idxl, freq_chunks[idx_chunk], idx_source, nskip);
          auto i1 = idx_chunks[idx_chunk][0];
          auto len_chunk = freq_chunks[idx_chunk].size();
          MATRIX vec_raw_l =
              MATRIX::Zero(3 * params.num_receivers(), len_chunk);
          for (int idx = (int) len_chunk - 1; idx > -1; --idx) {
            auto wval = freq_chunks[idx_chunk][idx];
            std::size_t ridx = ridxsph[idx];
            std::size_t len_ms = lensph - ridx;
            Complex w = wval + ieps;
            SMATRIX w_s = h_s.block(ridx, ridx, len_ms, len_ms) -
                          w * w * p_s.block(ridx, ridx, len_ms, len_ms);
            if (params.attenuation())
              w_s += attenFactor(wval, w0, twodivpi, myi) *
                     h_sa.block(ridx, ridx, len_ms, len_ms);
            w_s.makeCompressed();
            MATRIX f_red = F_BASE.block(ridx, 0, len_ms, F_BASE.cols());
            factorizeOrCompute(solver1, w_s, (int) len_chunk - idx - 1, nskip);
            MATRIX vec_sol = solver1.solve(f_red);
            auto lidx = lowidx - ridx;
            vec_raw_l.col(idx) +=
                RED_C.cwiseProduct(RV_BASE * vec_sol.block(lidx, 0, lenidx, 4))
                    .rowwise()
                    .sum();
          }
#pragma omp critical(torvecadd)
          {
            vec_raw.block(0, i1, 3 * params.num_receivers(), len_chunk) +=
                vec_raw_l;
          }
        }
      }
    }
  }

  // backazimuth rotation
  for (int idxr = 0; idxr < params.num_receivers(); ++idxr) {
    auto baz = srinfo.backazimuths(idxr);
    auto sbaz = std::sin(baz), cbaz = std::cos(baz);
    for (int j = 0; j < vec_raw.cols(); ++j) {
      auto tmp1 = vec_raw(3 * idxr + 1, j);
      auto tmp2 = vec_raw(3 * idxr + 2, j);
      vec_raw(3 * idxr + 1, j) = -tmp1 * cbaz + tmp2 * sbaz;
      vec_raw(3 * idxr + 2, j) = tmp1 * sbaz + tmp2 * cbaz;
    }
  }

  // output factor
  for (int idx_chunk = 0; idx_chunk < num_chunks; ++idx_chunk) {
    for (int idx = 0; idx < (int) idx_chunks[idx_chunk].size(); ++idx) {
      auto wval = freq_chunks[idx_chunk][idx];
      auto idxw = idx_chunks[idx_chunk][idx];
      Complex w = wval + ieps;
      vec_raw.col(idxw) *= outputFactor(params.output_type(), w, myi);
    }
  }
  return vec_raw;
}

}   // namespace SPARSESPEC
#endif   // FULL_SPEC_MULTI_SEM_H
