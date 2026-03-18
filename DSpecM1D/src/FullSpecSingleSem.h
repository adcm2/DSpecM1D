#ifndef DSPECM1D_FULL_SPEC_SINGLE_SEM_H
#define DSPECM1D_FULL_SPEC_SINGLE_SEM_H

#include "FullSpec.h"

namespace SPARSESPEC {

inline Eigen::MatrixXcd
SparseFSpec::spectra(InputParametersNew &paramsNew, Full1D::SEM &sem) {
  SpectraRunContext request(paramsNew.freqFull(), paramsNew.cmt(),
                            paramsNew.inputParameters(), paramsNew.tref(),
                            paramsNew.nskip());
  return spectra(request, sem);
}

inline Eigen::MatrixXcd
SparseFSpec::spectra(const SpectraRunContext &request, Full1D::SEM &sem) {
  using Complex = std::complex<double>;
  using MatrixC = Eigen::MatrixXcd;
  using SparseMatrixC = Eigen::SparseMatrix<Complex>;
  using SparseLUType =
      Eigen::SparseLU<SparseMatrixC, Eigen::COLAMDOrdering<int>>;

  auto &myff = request.freqFull();
  auto &cmt = request.cmt();
  auto &params = request.params();
  const int nskip = request.nskip();

  auto vecW = myff.w();
  SpecConstants sc(myff.ep(), request.tref());
  auto myi = sc.myi;
  auto ieps = sc.ieps;
  auto w0 = sc.w0;
  auto twodivpi = sc.twodivpi;

  MatrixC vecRaw = MatrixC::Zero(3 * params.num_receivers(), vecW.size());
  SparseLUType solver, solver1;

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
    SparseMatrixC keR = sem.hR().cast<Complex>();
    SparseMatrixC inR = sem.pR().cast<Complex>();
    SparseMatrixC keRAtten = sem.hRa().cast<Complex>();
    keR.makeCompressed();
    inR.makeCompressed();
    keRAtten.makeCompressed();
    MatrixC fR = sem.calculateForceR(cmt);
    std::vector<MatrixC> vecRvZ;
    for (int idxr = 0; idxr < numRec; ++idxr)
      vecRvZ.push_back(sem.rvZR(params, idxr).block(lowidx, 0, lenidx, 1));

#pragma omp parallel default(shared) private(solver)
    {
#pragma omp for schedule(dynamic, 10)
      for (int idx = myff.i1(); idx < myff.i2(); ++idx) {
        Complex w = vecW[idx] + ieps;
        SparseMatrixC wR = -w * w * inR + keR;
        if (params.attenuation())
          wR += attenFactor(vecW[idx], w0, twodivpi, myi) * keRAtten;
        wR.makeCompressed();
        solver.compute(wR);
        MatrixC vecX = solver.solve(fR);
        for (int idxr = 0; idxr < numRec; ++idxr) {
#pragma omp critical(torvecadd)
          {
            vecRaw(3 * idxr, idx) +=
                vecRvZ[idxr]
                    .cwiseProduct(vecX.block(lowidx, 0, lenidx, 1))
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
        SparseMatrixC hTor = sem.hTk(idxl).cast<Complex>();
        SparseMatrixC pTor = sem.pTk(idxl).cast<Complex>();
        SparseMatrixC hTorAtten = sem.hTa(idxl).cast<Complex>();
        hTor.makeCompressed();
        pTor.makeCompressed();
        hTorAtten.makeCompressed();
        MatrixC rvVals = paramInfo.rvFullTor(idxl);
        MatrixC fVals = sem.calculateForceCoefficientsT(cmt, idxl);
        MatrixC redC = rvVals * fVals;
        MatrixC fBase = sem.calculateForceAllT(cmt, idxl);
        MatrixC rvBase = sem.rvBaseFullT(params, idxl);
        auto ridxtor =
            SpectralTools::allIndicesTor(sem, idxl, myff, idxSource, nskip);
        for (int idx = myff.i2() - 1; idx > myff.i1() - 1; --idx) {
          std::size_t ridx = ridxtor[idx - myff.i1()];
          std::size_t len_ms = lentor - ridx;
          Complex w = vecW[idx] + ieps;
          SparseMatrixC matWTorRed =
              hTor.block(ridx, ridx, len_ms, len_ms) -
              w * w * pTor.block(ridx, ridx, len_ms, len_ms);
          if (params.attenuation())
            matWTorRed += attenFactor(vecW[idx], w0, twodivpi, myi) *
                          hTorAtten.block(ridx, ridx, len_ms, len_ms);
          matWTorRed.makeCompressed();
          auto fRed = fBase.block(ridx, 0, len_ms, fBase.cols());
          factorizeOrCompute(solver1, matWTorRed, myff.i2() - idx - 1, nskip);
          MatrixC vecSol = solver1.solve(fRed);
          auto lidx = lowidx - ridx;
#pragma omp critical(torvecadd)
          {
            vecRaw.col(idx) +=
                redC.cwiseProduct(rvBase * vecSol.block(lidx, 0, lenidx, 2))
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
        SparseMatrixC hS = sem.hS(idxl).cast<Complex>();
        SparseMatrixC pS = sem.pS(idxl).cast<Complex>();
        SparseMatrixC hSa = sem.hSa(idxl).cast<Complex>();
        hS.makeCompressed();
        pS.makeCompressed();
        hSa.makeCompressed();
        MatrixC rvVals = paramInfo.rvFullSph(idxl);
        MatrixC fVals = sem.calculateForceCoefficients(cmt, idxl);
        MatrixC redC = rvVals * fVals;
        MatrixC fBase = sem.calculateForceAll(cmt, idxl);
        MatrixC rvBase = sem.rvBaseFull(params, idxl);
        auto vecRidx =
            SpectralTools::allIndicesSph(sem, idxl, myff, idxSource, nskip);
        // MatrixC vecRawL =
        //     MatrixC::Zero(3 * params.num_receivers(), vecW.size());
        for (int idx = myff.i2() - 1; idx > myff.i1() - 1; --idx) {
          std::size_t idxRs = vecRidx[idx - myff.i1()];
          std::size_t len_ms = lensph - idxRs;
          Complex w = vecW[idx] + ieps;
          SparseMatrixC matSph = hS.block(idxRs, idxRs, len_ms, len_ms) -
                                 w * w * pS.block(idxRs, idxRs, len_ms, len_ms);
          if (params.attenuation())
            matSph += attenFactor(vecW[idx], w0, twodivpi, myi) *
                      hSa.block(idxRs, idxRs, len_ms, len_ms);
          matSph.makeCompressed();
          auto fRed = fBase.block(idxRs, 0, len_ms, fBase.cols());
          factorizeOrCompute(solver1, matSph, myff.i2() - idx - 1, nskip);
          MatrixC vecSol = solver1.solve(fRed);
          auto lidx = lowidx - idxRs;

#pragma omp critical(sphvecadd)
          {
            vecRaw.col(idx) +=
                redC.cwiseProduct(rvBase * vecSol.block(lidx, 0, lenidx, 4))
                    .rowwise()
                    .sum();
          }
        }
      }
    }
  }

  for (int idx = 0; idx < vecRaw.cols(); ++idx) {
    Complex w = vecW[idx] + ieps;
    vecRaw.col(idx) *= outputFactor(params.output_type(), w, myi);
  }
  return vecRaw;
}

template <class model1d>
auto
SparseFSpec::spectra(SpectraSolver::FreqFull &myff, Full1D::SEM &sem,
                     model1d &inp_model, SourceInfo::EarthquakeCMT &cmt,
                     InputParameters &params, int nskip) {
  SpectraRunContext request(myff, cmt, params, inp_model, nskip);
  return spectra(request, sem);
}

}   // namespace SPARSESPEC
#endif   // DSPECM1D_FULL_SPEC_SINGLE_SEM_H
