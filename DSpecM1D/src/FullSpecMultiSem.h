#ifndef DSPECM1D_FULL_SPEC_MULTI_SEM_H
#define DSPECM1D_FULL_SPEC_MULTI_SEM_H

#include "FullSpec.h"

namespace SPARSESPEC {

inline Eigen::MatrixXcd
SparseFSpec::spectra(InputParametersNew &paramsNew) {
  auto &params = paramsNew.inputParameters();
  return spectra(paramsNew.freqFull(), paramsNew.earthModel(), paramsNew.cmt(),
                 params, paramsNew.nq(), paramsNew.srInfo(),
                 params.relative_error());
}

inline Eigen::MatrixXcd SparseFSpec::spectra(
    const InputParametersNew &prepared,
    const EarthModels::IsotropicModel &model) {
  auto trial = prepared.withEarthModel(model);
  return spectra(trial);
}

template <class model1d>
Eigen::MatrixXcd
SparseFSpec::spectra(SpectraSolver::FreqFull &myff, model1d &inp_model,
                     SourceInfo::EarthquakeCMT &cmt, InputParameters &params,
                     int NQ, SRInfo &srInfo, double relerr) {
  return multiSemSpectra(myff, inp_model, cmt, params, NQ, srInfo, relerr,
                         nullptr);
}

inline PostUndampingSHResult
SparseFSpec::postUndampingSH(InputParametersNew &paramsNew) {
  auto &params = paramsNew.inputParameters();
  if (params.type() != 2 || params.lmax() < 1)
    throw std::invalid_argument(
        "Post-undamping SH requires explicit toroidal mode (type 2).");

  Timer timer;
  detail::SHDegreeSpectra spectra;
  auto &freq = paramsNew.freqFull();
  timer.start();
  multiSemSpectra(freq, paramsNew.earthModel(), paramsNew.cmt(), params,
                  paramsNew.nq(), paramsNew.srInfo(), params.relative_error(),
                  &spectra);
  timer.stop();

  PostUndampingSHResult result;
  result.radialSeconds = timer.duration_seconds();
  result.sampleIntervalSeconds = 0.5 * freq.dt() * freq.timeNorm();
  result.solvedPairs = spectra.solvedPairs;
  result.fullGridPairs = spectra.fullGridPairs;
  timer.start();
  result.timeSeries = detail::synthesizePostUndampingSH(
      freq, paramsNew.srInfo(), params.num_receivers(), paramsNew.normFactor(),
      spectra);
  timer.stop();
  result.synthesisSeconds = timer.duration_seconds();
  return result;
}

inline PostUndampingWindowResult
SparseFSpec::postUndampingWindow(
    InputParametersNew &paramsNew, const PostUndampingWindowOptions &options) {
  auto &params = paramsNew.inputParameters();
  const auto expectedType = options.component ==
                                    PostUndampingWindowComponent::Transverse
                                ? 4
                                : 3;
  if (params.type() != expectedType || params.lmax() < 1)
    throw std::invalid_argument(
        "Post-undamping phase window requires explicit type 4 for T or type 3 for Z.");
  if (options.component == PostUndampingWindowComponent::Transverse) {
    for (int receiver = 0; receiver < params.num_receivers(); ++receiver) {
      const double angle = paramsNew.srInfo().epicentralAngles(receiver);
      if (!std::isfinite(angle) || std::abs(std::sin(angle)) < 1e-12)
        throw std::invalid_argument(
            "Post-undamping T is undefined for coincident or antipodal geometry.");
    }
  }
  detail::validatePostUndampingWindowOptions(options);
  Timer timer;
  detail::WindowDegreeSpectra spectra;
  spectra.options = options;
  auto &freq = paramsNew.freqFull();
  timer.start();
  multiSemSpectra(freq, paramsNew.earthModel(), paramsNew.cmt(), params,
                  paramsNew.nq(), paramsNew.srInfo(), params.relative_error(),
                  nullptr, &spectra);
  timer.stop();

  PostUndampingWindowResult result;
  result.radialSeconds = timer.duration_seconds();
  result.sampleIntervalSeconds = 0.5 * freq.dt() * freq.timeNorm();
  result.solvedPairs = spectra.solvedPairs;
  result.fullGridPairs = spectra.fullGridPairs;
  timer.start();
  result.timeSeries = detail::synthesizePostUndampingWindow(
      freq, params.num_receivers(), paramsNew.normFactor(), spectra);
  timer.stop();
  result.synthesisSeconds = timer.duration_seconds();
  return result;
}

inline PostUndampingWindowResult SparseFSpec::postUndampingWindow(
    const InputParametersNew &prepared,
    const EarthModels::IsotropicModel &model,
    const PostUndampingWindowOptions &options) {
  auto trial = prepared.withEarthModel(model);
  return postUndampingWindow(trial, options);
}

template <class model1d>
Eigen::MatrixXcd
SparseFSpec::multiSemSpectra(SpectraSolver::FreqFull &myff, model1d &inp_model,
                             SourceInfo::EarthquakeCMT &cmt,
                             InputParameters &params, int NQ, SRInfo &srInfo,
                             double relerr, detail::SHDegreeSpectra *sh,
                             detail::WindowDegreeSpectra *window) {
  using Complex = std::complex<double>;
  using MatrixC = Eigen::MatrixXcd;
  using SparseMatrixC = Eigen::SparseMatrix<Complex>;
  using SparseLUType =
      Eigen::SparseLU<SparseMatrixC, Eigen::COLAMDOrdering<int>>;
  Timer timer1;

  if (params.attenuation())
    validateAttenuationReferencePeriod(inp_model.TREF());

  auto vecW = myff.w();
  SpecConstants sc(myff.ep(), inp_model.TREF());
  auto myi = sc.myi;
  auto ieps = sc.ieps;
  auto w0 = sc.w0;
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

  MatrixC vecRaw;
  if (!sh)
    vecRaw = MatrixC::Zero(3 * numRec, vecW.size());

  // std::cout << "Just before max step printout\n";
  for (auto s : maxSteps)
    std::cout << "Max step: " << s << "\n";
  // std::cout << "Just after max step printout\n\n";

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

  if (window) {
    inc_rad = false;
    inc_tor = window->options.component ==
              PostUndampingWindowComponent::Transverse;
    inc_sph = true;
  }

  if (sh) {
    sh->fullGridPairs = static_cast<std::size_t>(lmax - lmin + 1) *
                        (iEnd - iBegin);
    sh->degrees.reserve(lmax - lmin + 1);
    int firstBin = iBegin;
    for (int degree = lmin; degree <= lmax; ++degree) {
      const double cutoff = detail::shFirstSolveFrequencyHz(degree);
      while (firstBin < iEnd && myff.f(firstBin) / myff.timeNorm() < cutoff)
        ++firstBin;
      sh->degrees.push_back(
          {degree, firstBin, MatrixC(2 * numRec, iEnd - firstBin)});
    }
  }
  if (window) {
    const auto bounds =
        detail::validatePostUndampingWindowOptions(window->options);
    const std::size_t familyCount =
        window->options.component == PostUndampingWindowComponent::Transverse
            ? 2
            : 1;
    window->fullGridPairs = static_cast<std::size_t>(lmax - lmin + 1) *
                            static_cast<std::size_t>(iEnd - iBegin) *
                            familyCount;
    window->degrees.reserve(lmax - lmin + 1);
    for (int degree = lmin; degree <= lmax; ++degree) {
      int firstBin = iBegin;
      while (firstBin < iEnd &&
             detail::degreeSlownessSPerDegree(
                 degree, myff.f(firstBin) / myff.timeNorm()) > bounds.second)
        ++firstBin;
      int lastBin = firstBin;
      while (lastBin < iEnd &&
             detail::degreeSlownessSPerDegree(
                 degree, myff.f(lastBin) / myff.timeNorm()) >= bounds.first)
        ++lastBin;
      window->degrees.push_back(
          {degree, firstBin, lastBin,
           MatrixC::Zero(numRec, lastBin - firstBin)});
    }
  }

  // Validate each chunk's receiver support before entering any OpenMP region.
  std::vector<std::vector<int>> receiverElementsByChunk;
  receiverElementsByChunk.reserve(sems.size());
  for (const auto &sem : sems) {
    auto receiverElems = sem.receiverElements(params, inc_tor);
    receiverElementsByChunk.push_back(std::move(receiverElems));
  }

  // radials
  if (inc_rad) {
    timer1.start();
    for (int idxChunk = 0; idxChunk < numChunks; ++idxChunk) {
      Full1D::SEM &sem = sems[idxChunk];
      const auto &recElems = receiverElementsByChunk[idxChunk];
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
            wR += attenFactor(w, w0) * keRAtten;
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
    std::size_t solvedPairs = 0;
#pragma omp parallel default(shared) private(solver1) reduction(+:solvedPairs)
    {
#pragma omp for schedule(dynamic)
      for (int idxl = lmin; idxl < lmax + 1; ++idxl) {
        auto *degree = sh ? &sh->degrees[idxl - lmin] : nullptr;
        auto *windowDegree = window ? &window->degrees[idxl - lmin] : nullptr;
        if (degree && degree->firstBin == iEnd)
          continue;

        MatrixC rvVals = srInfo.rvRedTor(idxl);
        for (int idxChunk = 0; idxChunk < numChunks; ++idxChunk) {
          const auto &indices = idxChunks[idxChunk];
          if (indices.empty())
            continue;
          const int i1 = indices.front();
          const int lenChunk = static_cast<int>(indices.size());
          int first = degree ? std::max(0, degree->firstBin - i1) : 0;
          if (window) {
            while (first < lenChunk) {
              const int globalIdx = indices[first];
              if (globalIdx >= windowDegree->firstBin)
                break;
              ++first;
            }
          }
          if (first >= lenChunk ||
              (window && indices[first] >= windowDegree->lastBin))
            continue;
          Full1D::SEM &sem = sems[idxChunk];
          auto idxSource = sem.sourceElement(cmt);
          const auto &recElems = receiverElementsByChunk[idxChunk];
          // std::cout << "Debug 0\n";
          auto lowidx = sem.ltgT(recElems[0], 0);
          // std::cout << "Debug 1\n";
          auto upidx = sem.ltgT(recElems.back(), NQ - 1) + 1;
          // std::cout << "Debug 2\n";
          int lenidx = upidx - lowidx;
          auto lentor = sem.ltgT(sem.eu() - 1, NQ - 1) + 1;
          // std::cout << "Debug 3\n";

          SparseMatrixC hTor = sem.hTk(idxl).cast<Complex>();
          SparseMatrixC pTor = sem.pTk(idxl).cast<Complex>();
          SparseMatrixC hTorAtten = sem.hTa(idxl).cast<Complex>();
          hTor.makeCompressed();
          pTor.makeCompressed();
          hTorAtten.makeCompressed();

          // std::cout << "Debug 4\n";
          MatrixC fVals = sem.calculateForceRedCoefficientsT(cmt, idxl, 0.0);
          // std::cout << "Debug 5\n";
          MatrixC redC = rvVals * fVals;
          MatrixC fBase = sem.calculateForceAllT(cmt, idxl);
          // std::cout << "Debug 6\n";
          MatrixC rvBase = sem.rvBaseFullT(params, idxl);
          // std::cout << "Debug 7\n";
          auto ridxtor = SpectralTools::allIndicesTor(
              sem, idxl, freqChunks[idxChunk], idxSource, nskip);
          const bool capture = degree || windowDegree;
          MatrixC vecRawL = MatrixC::Zero(3 * numRec, capture ? 1 : lenChunk);
          // Retain the high-frequency suffix of the original chunk. The
          // index and factorization cadence still start at its highest bin.
          for (int idx = lenChunk - 1; idx >= first; --idx) {
            auto wval = freqChunks[idxChunk][idx];
            const int globalIdx = indices[idx];
            if (window) {
              if (globalIdx < windowDegree->firstBin ||
                  globalIdx >= windowDegree->lastBin)
                continue;
            }
            Complex w = wval + ieps;
            std::size_t ridx = ridxtor[idx];
            std::size_t len_ms = lentor - ridx;
            SparseMatrixC mat = hTor.block(ridx, ridx, len_ms, len_ms) -
                                w * w * pTor.block(ridx, ridx, len_ms, len_ms);
            if (params.attenuation())
              mat += attenFactor(w, w0) *
                     hTorAtten.block(ridx, ridx, len_ms, len_ms);
            mat.makeCompressed();
            auto fRed = fBase.block(ridx, 0, len_ms, fBase.cols());
            const int cadence = window &&
                                        globalIdx == windowDegree->lastBin - 1
                                    ? 0
                                    : lenChunk - idx - 1;
            factorizeOrCompute(solver1, mat, cadence, nskip);
            MatrixC vecSol = solver1.solve(fRed);
            auto lidx = lowidx - ridx;
            auto contribution = vecRawL.col(capture ? 0 : idx);
            if (capture)
              contribution.setZero();
            contribution +=
                redC.cwiseProduct(rvBase * vecSol.block(lidx, 0, lenidx, 2))
                    .rowwise()
                    .sum();
            if (degree) {
              const int out = indices[idx] - degree->firstBin;
              const auto factor = outputFactor(params.output_type(), w, myi);
              for (int receiver = 0; receiver < numRec; ++receiver) {
                degree->horizontal(2 * receiver, out) =
                    factor * contribution(3 * receiver + 1);
                degree->horizontal(2 * receiver + 1, out) =
                    factor * contribution(3 * receiver + 2);
              }
              ++solvedPairs;
            } else if (windowDegree) {
              const auto factor = outputFactor(params.output_type(), w, myi);
              for (int receiver = 0; receiver < numRec; ++receiver) {
                // True transverse is clockwise from outward R, which is the
                // negative of the stored source-centred transverse row.
                windowDegree->components(receiver,
                                         globalIdx - windowDegree->firstBin) +=
                    -factor * contribution(3 * receiver + 2);
              }
              ++solvedPairs;
            }
          }
          if (!capture) {
#pragma omp critical(torvecadd)
            {
              vecRaw.block(0, i1, 3 * numRec, lenChunk) += vecRawL;
            }
          }
        }
      }
    }
    if (sh) {
      sh->solvedPairs = solvedPairs;
      return vecRaw;
    }
    if (window)
      window->solvedPairs = solvedPairs;
  }

  // spheroidals
  if (inc_sph) {
    std::cout << "\nDoing Spheroidal Modes\n";
    std::size_t windowSolvedPairs = window ? window->solvedPairs : 0;
#pragma omp parallel default(shared) private(solver1) reduction(+:windowSolvedPairs)
    {
#pragma omp for schedule(dynamic)
      for (int idxl = lmin; idxl < lmax + 1; ++idxl) {
        auto *windowDegree = window ? &window->degrees[idxl - lmin] : nullptr;
        MatrixC rvVals = srInfo.rvRedSph(idxl);
        for (int idxChunk = 0; idxChunk < numChunks; ++idxChunk) {
          Full1D::SEM &sem = sems[idxChunk];
          auto idxSource = sem.sourceElement(cmt);
          const auto &recElems = receiverElementsByChunk[idxChunk];
          const int lenChunk = static_cast<int>(idxChunks[idxChunk].size());
          int selectedBegin = 0;
          int selectedEnd = lenChunk;
          if (window) {
            while (selectedBegin < lenChunk) {
              const int globalIdx = idxChunks[idxChunk][selectedBegin];
              if (globalIdx >= windowDegree->firstBin)
                break;
              ++selectedBegin;
            }
            selectedEnd = selectedBegin;
            while (selectedEnd < lenChunk) {
              const int globalIdx = idxChunks[idxChunk][selectedEnd];
              if (globalIdx >= windowDegree->lastBin)
                break;
              ++selectedEnd;
            }
            if (selectedBegin >= selectedEnd)
              continue;
          }
          MatrixC fVals = sem.calculateForceRedCoefficients(cmt, idxl, 0.0);
          MatrixC redC = rvVals * fVals;
          MatrixC vecRawL = MatrixC::Zero(3 * numRec, lenChunk);
          const double cutoffMhz = params.cowlingFrequencyMhz();
          int cowlingFirst = static_cast<int>(idxChunks[idxChunk].size());
          if (cutoffMhz > 0.0) {
            for (int idx = 0; idx < static_cast<int>(idxChunks[idxChunk].size());
                 ++idx) {
              const int globalIdx = idxChunks[idxChunk][idx];
              const double frequencyMhz =
                  myff.f(globalIdx) * 1000.0 / myff.timeNorm();
              if (frequencyMhz >= cutoffMhz) {
                cowlingFirst = idx;
                break;
              }
            }
          }

          // Each contiguous region owns its matrix layout and restarts the
          // truncation/factorization cadence at the region boundary.
          for (int region = 0; region < 2; ++region) {
            const bool useCowling = region == 1;
            const int regionBegin = useCowling ? cowlingFirst : 0;
            const int regionEnd = useCowling
                                      ? static_cast<int>(idxChunks[idxChunk].size())
                                      : cowlingFirst;
            if (regionBegin >= regionEnd)
              continue;
            if (window &&
                std::max(regionBegin, selectedBegin) >=
                    std::min(regionEnd, selectedEnd))
              continue;

            const auto lowidx = useCowling
                                    ? sem.ltgSC(0, recElems[0], 0)
                                    : sem.ltgS(0, recElems[0], 0);
            const auto upidx = useCowling
                                   ? sem.ltgSC(1, recElems.back(), NQ - 1)
                                   : sem.ltgS(1, recElems.back(), NQ - 1);
            const int lenidx = upidx - lowidx + 1;
            const auto lensph = useCowling
                                    ? sem.ltgSC(1, sem.mesh().NE() - 1, NQ - 1) + 1
                                    : sem.ltgS(2, sem.mesh().NE() - 1, NQ - 1) + 1;
            SparseMatrixC pS =
                (useCowling ? sem.pSC(idxl) : sem.pS(idxl)).cast<Complex>();
            SparseMatrixC hS =
                (useCowling ? sem.hSC(idxl) : sem.hS(idxl)).cast<Complex>();
            SparseMatrixC hSa =
                (useCowling ? sem.hSCa(idxl) : sem.hSa(idxl)).cast<Complex>();
            pS.makeCompressed();
            hS.makeCompressed();
            hSa.makeCompressed();
            MatrixC fBase = sem.calculateForceAll(cmt, idxl, useCowling);
            MatrixC rvBase = sem.rvBaseFull(params, idxl, useCowling);

            std::vector<double> regionW(freqChunks[idxChunk].begin() + regionBegin,
                                        freqChunks[idxChunk].begin() + regionEnd);
            auto ridxsph = SpectralTools::allIndicesSph(
                sem, idxl, regionW, idxSource, nskip, useCowling);
            SparseLUType regionSolver;
            for (int idx = regionEnd - 1; idx >= regionBegin; --idx) {
              if (window && (idx < selectedBegin || idx >= selectedEnd))
                continue;
              const int localRegionIdx = idx - regionBegin;
              const double wval = freqChunks[idxChunk][idx];
              const std::size_t ridx = ridxsph[localRegionIdx];
              const std::size_t len_ms = lensph - ridx;
              const Complex w = wval + ieps;
              SparseMatrixC wS = hS.block(ridx, ridx, len_ms, len_ms) -
                                 w * w * pS.block(ridx, ridx, len_ms, len_ms);
              if (params.attenuation())
                wS += attenFactor(w, w0) *
                      hSa.block(ridx, ridx, len_ms, len_ms);
              wS.makeCompressed();
              MatrixC fRed = fBase.block(ridx, 0, len_ms, fBase.cols());
              const int globalIdx = idxChunks[idxChunk][idx];
              const int cadence = window &&
                                          globalIdx == windowDegree->lastBin - 1
                                      ? 0
                                      : regionEnd - idx - 1;
              factorizeOrCompute(regionSolver, wS, cadence, nskip);
              MatrixC vecSol = regionSolver.solve(fRed);
              const auto lidx = lowidx - ridx;
              const auto contribution =
                  redC.cwiseProduct(rvBase * vecSol.block(lidx, 0, lenidx, 4))
                      .rowwise()
                      .sum();
              if (windowDegree) {
                const auto receiverContribution = contribution.eval();
                const auto factor = outputFactor(params.output_type(), w, myi);
                if (window->options.component ==
                    PostUndampingWindowComponent::Transverse) {
                  for (int receiver = 0; receiver < numRec; ++receiver)
                    windowDegree->components(
                        receiver, globalIdx - windowDegree->firstBin) +=
                        -factor * receiverContribution(3 * receiver + 2);
                } else {
                  for (int receiver = 0; receiver < numRec; ++receiver)
                    windowDegree->components(
                        receiver, globalIdx - windowDegree->firstBin) +=
                        factor * receiverContribution(3 * receiver);
                }
                ++windowSolvedPairs;
              } else {
                vecRawL.col(idx) += contribution;
              }
            }
          }
          if (lenChunk > 0) {
            const int i1 = idxChunks[idxChunk].front();
#pragma omp critical(torvecadd)
            {
              vecRaw.block(0, i1, 3 * numRec, lenChunk) += vecRawL;
            }
          }
        }
      }
    }
    if (window) {
      window->solvedPairs = windowSolvedPairs;
      return vecRaw;
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
