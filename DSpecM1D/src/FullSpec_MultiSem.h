#ifndef FULL_SPEC_MULTI_SEM_H
#define FULL_SPEC_MULTI_SEM_H

#include "FullSpec.h"

namespace SPARSESPEC {

template <class model1d>
auto
SparseFSpec::spectra(SpectraSolver::FreqFull &myff, model1d &inp_model,
                     SourceInfo::EarthquakeCMT &cmt, InputParameters &params,
                     int NQ, SRInfo &srInfo, double relerr) {
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

  auto numRec = params.num_receivers();
  int nskip = std::max(1, (myff.i2() - myff.i1()) / 20);
  auto iBegin = myff.i1();
  auto iEnd = myff.i2();   // one past last valid index
  if (iBegin < 0 || iEnd > static_cast<int>(vec_w.size()) || iBegin >= iEnd) {
    throw std::runtime_error(
        "Invalid frequency index range in FullSpec_MultiSem.");
  }

  int num_chunks = std::max(
      1, static_cast<int>(std::floor((myff.f22() - myff.f11()) / 10.0) + 1.0));

  std::vector<std::vector<double>> freq_chunks(num_chunks);
  std::vector<std::vector<int>> idx_chunks(num_chunks);

  double wl = vec_w[iBegin];
  double wh = vec_w[iEnd - 1];
  double wdiff = (wh - wl) / static_cast<double>(num_chunks);

  int idxw = iBegin;
  for (int c = 0; c < num_chunks; ++c) {
    double upper = (c == num_chunks - 1) ? wh : (wl + wdiff * (c + 1));
    while (idxw < iEnd && vec_w[idxw] <= upper) {
      freq_chunks[c].push_back(vec_w[idxw]);
      idx_chunks[c].push_back(idxw);
      ++idxw;
    }
  }

  // Any stragglers (numeric edge cases) go to last chunk
  while (idxw < iEnd) {
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

  // build one SEM per chunk
  std::vector<Full1D::SEM> sems;
  for (int idx = 0; idx < num_chunks; ++idx)
    sems.emplace_back(inp_model, max_steps[idx], NQ, params.lmax());

  timer1.start();
  SLU solver, solver1;

  int lmin = params.lmin();
  int lmax = params.lmax();
  ParamRedInfo paramInfo(params, lmax);

  ModeFlags flags = resolveModeFlags(params.type(), lmin, lmax);
  bool inc_rad = flags.inc_rad;
  bool inc_tor = flags.inc_tor;
  bool inc_sph = flags.inc_sph;
  lmin = std::max(lmin, 1);

  // radials
  if (inc_rad) {
    timer1.start();
    for (int idx_chunk = 0; idx_chunk < num_chunks; ++idx_chunk) {
      Full1D::SEM &sem = sems[idx_chunk];
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
      MATRIX fR = sem.calculateForceRedR(cmt);
      MATRIX vecRedZ = sem.rvRedZR(params);
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
          MATRIX vecX = solver.solve(fR);
          auto tmp = vecX.block(lowidx, 0, lenidx, 1);
          auto resval = (vecRedZ.transpose() * tmp).sum();
          for (int idxr = 0; idxr < numRec; ++idxr) {
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
        MATRIX rvVals = srInfo.rvRedTor(idxl);
        for (int idx_chunk = 0; idx_chunk < num_chunks; ++idx_chunk) {
          Full1D::SEM &sem = sems[idx_chunk];
          auto idxSource = sem.sourceElement(cmt);
          auto recElems = sem.receiverElements(params);
          auto lowidx = sem.ltgT(recElems[0], 0);
          auto upidx = sem.ltgT(recElems.back(), NQ - 1) + 1;
          int lenidx = upidx - lowidx;
          auto lentor = sem.ltgT(sem.mesh().NE() - 1, NQ - 1) + 1;
          SMATRIX H_tor = sem.hTk(idxl).cast<Complex>();
          SMATRIX P_tor = sem.pTk(idxl).cast<Complex>();
          SMATRIX H_tor_atten = sem.hTa(idxl).cast<Complex>();
          H_tor.makeCompressed();
          P_tor.makeCompressed();
          H_tor_atten.makeCompressed();
          MATRIX fVals = sem.calculateForceRedCoefficientsT(cmt, idxl, 0.0);
          MATRIX redC = rvVals * fVals;
          MATRIX fBase = sem.calculateForceAllT(cmt, idxl);
          MATRIX rvBase = sem.rvBaseFullT(params, idxl);
          auto ridxtor = SpectralTools::allIndicesTor(
              sem, idxl, freq_chunks[idx_chunk], idxSource, nskip);
          auto i1 = idx_chunks[idx_chunk][0];
          auto lenChunk = freq_chunks[idx_chunk].size();
          MATRIX vecRawL = MATRIX::Zero(3 * params.num_receivers(), lenChunk);
          for (int idx = (int) lenChunk - 1; idx > -1; --idx) {
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
            auto fRed = fBase.block(ridx, 0, len_ms, fBase.cols());
            factorizeOrCompute(solver1, mat, (int) lenChunk - idx - 1, nskip);
            MATRIX vecSol = solver1.solve(fRed);
            auto lidx = lowidx - ridx;
            vecRawL.col(idx) +=
                redC.cwiseProduct(rvBase * vecSol.block(lidx, 0, lenidx, 2))
                    .rowwise()
                    .sum();
          }
#pragma omp critical(torvecadd)
          {
            vec_raw.block(0, i1, 3 * params.num_receivers(), lenChunk) +=
                vecRawL;
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
        MATRIX rvVals = srInfo.rvRedSph(idxl);
        for (int idx_chunk = 0; idx_chunk < num_chunks; ++idx_chunk) {
          Full1D::SEM &sem = sems[idx_chunk];
          auto idxSource = sem.sourceElement(cmt);
          auto recElems = sem.receiverElements(params);
          auto lowidx = sem.ltgS(0, recElems[0], 0);
          auto upidx = sem.ltgS(1, recElems.back(), NQ - 1);
          int lenidx = upidx - lowidx + 1;
          auto lensph = sem.ltgS(2, sem.mesh().NE() - 1, NQ - 1) + 1;
          SMATRIX p_s = sem.pS(idxl).cast<Complex>();
          SMATRIX h_s = sem.hS(idxl).cast<Complex>();
          SMATRIX h_sa = sem.hSa(idxl).cast<Complex>();
          p_s.makeCompressed();
          h_s.makeCompressed();
          h_sa.makeCompressed();
          MATRIX fVals = sem.calculateForceRedCoefficients(cmt, idxl, 0.0);
          MATRIX redC = rvVals * fVals;
          MATRIX fBase = sem.calculateForceAll(cmt, idxl);
          MATRIX rvBase = sem.rvBaseFull(params, idxl);
          auto ridxsph = SpectralTools::allIndicesSph(
              sem, idxl, freq_chunks[idx_chunk], idxSource, nskip);
          auto i1 = idx_chunks[idx_chunk][0];
          auto lenChunk = freq_chunks[idx_chunk].size();
          MATRIX vecRawL = MATRIX::Zero(3 * params.num_receivers(), lenChunk);
          for (int idx = (int) lenChunk - 1; idx > -1; --idx) {
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
            MATRIX fRed = fBase.block(ridx, 0, len_ms, fBase.cols());
            factorizeOrCompute(solver1, w_s, (int) lenChunk - idx - 1, nskip);
            MATRIX vecSol = solver1.solve(fRed);
            auto lidx = lowidx - ridx;
            vecRawL.col(idx) +=
                redC.cwiseProduct(rvBase * vecSol.block(lidx, 0, lenidx, 4))
                    .rowwise()
                    .sum();
          }
#pragma omp critical(torvecadd)
          {
            vec_raw.block(0, i1, 3 * params.num_receivers(), lenChunk) +=
                vecRawL;
          }
        }
      }
    }
  }

  // backazimuth rotation
  for (int idxr = 0; idxr < params.num_receivers(); ++idxr) {
    auto baz = srInfo.backazimuths(idxr);
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
