#ifndef DSPECM1D_FULL_SPEC_MULTI_SEM_H
#define DSPECM1D_FULL_SPEC_MULTI_SEM_H

#include "FullSpec.h"

namespace SPARSESPEC {

template <class model1d>
auto
SparseFSpec::spectra(SpectraSolver::FreqFull &myff, model1d &inp_model,
                     SourceInfo::EarthquakeCMT &cmt, InputParameters &params,
                     int NQ, SRInfo &srInfo, double relerr) {
  using Complex = std::complex<double>;
  using MatrixC = Eigen::MatrixXcd;
  using SparseMatrixC = Eigen::SparseMatrix<Complex>;
  using SparseLUType =
      Eigen::SparseLU<SparseMatrixC, Eigen::COLAMDOrdering<int>>;
  Timer timer1;

  auto vecW = myff.w();
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
  if (iBegin < 0 || iEnd > static_cast<int>(vecW.size()) || iBegin >= iEnd) {
    throw std::runtime_error(
        "Invalid frequency index range in FullSpec_MultiSem.");
  }

  int numChunks = std::max(
      1, static_cast<int>(std::floor((myff.f22() - myff.f11()) / 10.0) + 1.0));

  std::vector<std::vector<double>> freqChunks(numChunks);
  std::vector<std::vector<int>> idxChunks(numChunks);

  double wl = vecW[iBegin];
  double wh = vecW[iEnd - 1];
  double wdiff = (wh - wl) / static_cast<double>(numChunks);

  int idxw = iBegin;
  for (int c = 0; c < numChunks; ++c) {
    double upper = (c == numChunks - 1) ? wh : (wl + wdiff * (c + 1));
    while (idxw < iEnd && vecW[idxw] <= upper) {
      freqChunks[c].push_back(vecW[idxw]);
      idxChunks[c].push_back(idxw);
      ++idxw;
    }
  }

  // Any stragglers (numeric edge cases) go to last chunk
  while (idxw < iEnd) {
    freqChunks.back().push_back(vecW[idxw]);
    idxChunks.back().push_back(idxw);
    ++idxw;
  }

  // Build maxSteps safely
  std::vector<double> maxSteps;
  maxSteps.reserve(numChunks);
  double baseLen = 0.78 * 1000 / inp_model.TimeNorm() * twopid;
  double newLen = std::pow(100.0 * relerr, 1.0 / (NQ - 1.0)) * baseLen;

  for (int c = 0; c < numChunks; ++c) {
    double fref = freqChunks[c].empty() ? wh : freqChunks[c].back();
    if (fref <= 0.0) {
      throw std::runtime_error("Non-positive reference frequency in chunking.");
    }
    double step = std::min(newLen / fref, 0.05);
    maxSteps.push_back(step);
  }

  MatrixC vecRaw = MatrixC::Zero(3 * params.num_receivers(), vecW.size());

  for (auto s : maxSteps)
    std::cout << "Max step: " << s << "\n";

  // build one SEM per chunk
  std::vector<Full1D::SEM> sems;
  for (int idx = 0; idx < numChunks; ++idx)
    sems.emplace_back(inp_model, maxSteps[idx], NQ, params.lmax());

  timer1.start();
  SparseLUType solver, solver1;

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
    for (int idxChunk = 0; idxChunk < numChunks; ++idxChunk) {
      Full1D::SEM &sem = sems[idxChunk];
      auto recElems = sem.receiverElements(params);
      auto lowidx = sem.ltgR(0, recElems[0], 0);
      auto upidx = sem.ltgR(1, recElems.back(), NQ - 1);
      int lenidx = upidx - lowidx + 1;
      SparseMatrixC keR = sem.hR().cast<Complex>();
      SparseMatrixC inR = sem.pR().cast<Complex>();
      SparseMatrixC keRAtten = sem.hRa().cast<Complex>();
      keR.makeCompressed();
      inR.makeCompressed();
      keRAtten.makeCompressed();
      MatrixC fR = sem.calculateForceRedR(cmt);
      MatrixC vecRedZ = sem.rvRedZR(params);
#pragma omp parallel default(shared) private(solver)
      {
#pragma omp for schedule(dynamic)
        for (int idx = 0; idx < (int) idxChunks[idxChunk].size(); ++idx) {
          double wval = freqChunks[idxChunk][idx];
          int idxw = idxChunks[idxChunk][idx];
          Complex w = wval + ieps;
          SparseMatrixC wR = -w * w * inR + keR;
          if (params.attenuation())
            wR += attenFactor(wval, w0, twodivpi, myi) * keRAtten;
          wR.makeCompressed();
          solver.compute(wR);
          MatrixC vecX = solver.solve(fR);
          auto tmp = vecX.block(lowidx, 0, lenidx, 1);
          auto resval = (vecRedZ.transpose() * tmp).sum();
          for (int idxr = 0; idxr < numRec; ++idxr) {
#pragma omp critical(torvecadd)
            {
              vecRaw(3 * idxr, idxw) += resval;
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
        MatrixC rvVals = srInfo.rvRedTor(idxl);
        for (int idxChunk = 0; idxChunk < numChunks; ++idxChunk) {
          Full1D::SEM &sem = sems[idxChunk];
          auto idxSource = sem.sourceElement(cmt);
          auto recElems = sem.receiverElements(params);
          auto lowidx = sem.ltgT(recElems[0], 0);
          auto upidx = sem.ltgT(recElems.back(), NQ - 1) + 1;
          int lenidx = upidx - lowidx;
          auto lentor = sem.ltgT(sem.mesh().NE() - 1, NQ - 1) + 1;
          SparseMatrixC hTor = sem.hTk(idxl).cast<Complex>();
          SparseMatrixC pTor = sem.pTk(idxl).cast<Complex>();
          SparseMatrixC hTorAtten = sem.hTa(idxl).cast<Complex>();
          hTor.makeCompressed();
          pTor.makeCompressed();
          hTorAtten.makeCompressed();
          MatrixC fVals = sem.calculateForceRedCoefficientsT(cmt, idxl, 0.0);
          MatrixC redC = rvVals * fVals;
          MatrixC fBase = sem.calculateForceAllT(cmt, idxl);
          MatrixC rvBase = sem.rvBaseFullT(params, idxl);
          auto ridxtor = SpectralTools::allIndicesTor(
              sem, idxl, freqChunks[idxChunk], idxSource, nskip);
          auto i1 = idxChunks[idxChunk][0];
          auto lenChunk = freqChunks[idxChunk].size();
          MatrixC vecRawL = MatrixC::Zero(3 * params.num_receivers(), lenChunk);
          for (int idx = (int) lenChunk - 1; idx > -1; --idx) {
            auto wval = freqChunks[idxChunk][idx];
            Complex w = wval + ieps;
            std::size_t ridx = ridxtor[idx];
            std::size_t len_ms = lentor - ridx;
            SparseMatrixC mat = hTor.block(ridx, ridx, len_ms, len_ms) -
                                w * w * pTor.block(ridx, ridx, len_ms, len_ms);
            if (params.attenuation())
              mat += attenFactor(wval, w0, twodivpi, myi) *
                     hTorAtten.block(ridx, ridx, len_ms, len_ms);
            mat.makeCompressed();
            auto fRed = fBase.block(ridx, 0, len_ms, fBase.cols());
            factorizeOrCompute(solver1, mat, (int) lenChunk - idx - 1, nskip);
            MatrixC vecSol = solver1.solve(fRed);
            auto lidx = lowidx - ridx;
            vecRawL.col(idx) +=
                redC.cwiseProduct(rvBase * vecSol.block(lidx, 0, lenidx, 2))
                    .rowwise()
                    .sum();
          }
#pragma omp critical(torvecadd)
          {
            vecRaw.block(0, i1, 3 * params.num_receivers(), lenChunk) +=
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
        MatrixC rvVals = srInfo.rvRedSph(idxl);
        for (int idxChunk = 0; idxChunk < numChunks; ++idxChunk) {
          Full1D::SEM &sem = sems[idxChunk];
          auto idxSource = sem.sourceElement(cmt);
          auto recElems = sem.receiverElements(params);
          auto lowidx = sem.ltgS(0, recElems[0], 0);
          auto upidx = sem.ltgS(1, recElems.back(), NQ - 1);
          int lenidx = upidx - lowidx + 1;
          auto lensph = sem.ltgS(2, sem.mesh().NE() - 1, NQ - 1) + 1;
          SparseMatrixC pS = sem.pS(idxl).cast<Complex>();
          SparseMatrixC hS = sem.hS(idxl).cast<Complex>();
          SparseMatrixC hSa = sem.hSa(idxl).cast<Complex>();
          pS.makeCompressed();
          hS.makeCompressed();
          hSa.makeCompressed();
          MatrixC fVals = sem.calculateForceRedCoefficients(cmt, idxl, 0.0);
          MatrixC redC = rvVals * fVals;
          MatrixC fBase = sem.calculateForceAll(cmt, idxl);
          MatrixC rvBase = sem.rvBaseFull(params, idxl);
          auto ridxsph = SpectralTools::allIndicesSph(
              sem, idxl, freqChunks[idxChunk], idxSource, nskip);
          auto i1 = idxChunks[idxChunk][0];
          auto lenChunk = freqChunks[idxChunk].size();
          MatrixC vecRawL = MatrixC::Zero(3 * params.num_receivers(), lenChunk);
          for (int idx = (int) lenChunk - 1; idx > -1; --idx) {
            auto wval = freqChunks[idxChunk][idx];
            std::size_t ridx = ridxsph[idx];
            std::size_t len_ms = lensph - ridx;
            Complex w = wval + ieps;
            SparseMatrixC wS = hS.block(ridx, ridx, len_ms, len_ms) -
                               w * w * pS.block(ridx, ridx, len_ms, len_ms);
            if (params.attenuation())
              wS += attenFactor(wval, w0, twodivpi, myi) *
                    hSa.block(ridx, ridx, len_ms, len_ms);
            wS.makeCompressed();
            MatrixC fRed = fBase.block(ridx, 0, len_ms, fBase.cols());
            factorizeOrCompute(solver1, wS, (int) lenChunk - idx - 1, nskip);
            MatrixC vecSol = solver1.solve(fRed);
            auto lidx = lowidx - ridx;
            vecRawL.col(idx) +=
                redC.cwiseProduct(rvBase * vecSol.block(lidx, 0, lenidx, 4))
                    .rowwise()
                    .sum();
          }
#pragma omp critical(torvecadd)
          {
            vecRaw.block(0, i1, 3 * params.num_receivers(), lenChunk) +=
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
    for (int j = 0; j < vecRaw.cols(); ++j) {
      auto tmp1 = vecRaw(3 * idxr + 1, j);
      auto tmp2 = vecRaw(3 * idxr + 2, j);
      vecRaw(3 * idxr + 1, j) = -tmp1 * cbaz + tmp2 * sbaz;
      vecRaw(3 * idxr + 2, j) = tmp1 * sbaz + tmp2 * cbaz;
    }
  }

  // output factor
  for (int idxChunk = 0; idxChunk < numChunks; ++idxChunk) {
    for (int idx = 0; idx < (int) idxChunks[idxChunk].size(); ++idx) {
      auto wval = freqChunks[idxChunk][idx];
      auto idxw = idxChunks[idxChunk][idx];
      Complex w = wval + ieps;
      vecRaw.col(idxw) *= outputFactor(params.output_type(), w, myi);
    }
  }
  return vecRaw;
}

}   // namespace SPARSESPEC
#endif   // DSPECM1D_FULL_SPEC_MULTI_SEM_H
