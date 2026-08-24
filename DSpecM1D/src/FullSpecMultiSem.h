#ifndef DSPECM1D_FULL_SPEC_MULTI_SEM_H
#define DSPECM1D_FULL_SPEC_MULTI_SEM_H

#include "FullSpec.h"
#include "Profiling.h"

#include <map>
#include <memory>
#include <limits>
#include <vector>

#ifdef DSPECM1D_ENABLE_LAPACK_BAND_SOLVER
#include "LapackBandSolver.h"
#endif

namespace SPARSESPEC {

namespace detail {

using EigenMultiSemSolver =
    Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>,
                    Eigen::NaturalOrdering<int>>;

/// Owns one compressed union-pattern matrix and updates only its values for
/// each frequency.  The component matrices are copied into aligned value
/// arrays once when the active truncation changes; the hot path then performs
/// no sparse insertion, merging, or compression.
class FixedPatternFrequencyMatrix {
public:
  using Complex = std::complex<double>;
  using SparseMatrix = Eigen::SparseMatrix<Complex>;

  FixedPatternFrequencyMatrix(const SparseMatrix &h, const SparseMatrix &p,
                              const SparseMatrix &ha, bool attenuation)
      : matrix_(h.rows(), h.cols()), attenuation_(attenuation) {
    std::map<std::pair<Eigen::Index, Eigen::Index>, Eigen::Index> positions;
    auto collect = [&](const SparseMatrix &part) {
      for (Eigen::Index col = 0; col < part.outerSize(); ++col)
        for (typename SparseMatrix::InnerIterator entry(part, col); entry;
             ++entry)
          positions.emplace(std::make_pair(entry.row(), entry.col()), 0);
    };
    collect(h);
    collect(p);
    if (attenuation)
      collect(ha);

    std::vector<Eigen::Triplet<Complex>> triplets;
    triplets.reserve(positions.size());
    for (const auto &entry : positions)
      triplets.emplace_back(entry.first.first, entry.first.second, Complex{});
    matrix_.setFromTriplets(triplets.begin(), triplets.end());
    matrix_.makeCompressed();
    for (Eigen::Index col = 0; col < matrix_.outerSize(); ++col)
      for (Eigen::Index k = matrix_.outerIndexPtr()[col];
           k < matrix_.outerIndexPtr()[col + 1]; ++k)
        positions.at({matrix_.innerIndexPtr()[k], col}) = k;

    hValues_.assign(matrix_.nonZeros(), Complex{});
    pValues_.assign(matrix_.nonZeros(), Complex{});
    haValues_.assign(matrix_.nonZeros(), Complex{});
    copyValues(h, hValues_, positions);
    copyValues(p, pValues_, positions);
    if (attenuation)
      copyValues(ha, haValues_, positions);
  }

  void update(Complex w, Complex chi) {
    Complex *values = matrix_.valuePtr();
    if (attenuation_) {
      for (Eigen::Index k = 0; k < matrix_.nonZeros(); ++k)
        values[k] = hValues_[k] - (w * w) * pValues_[k] + chi * haValues_[k];
    } else {
      for (Eigen::Index k = 0; k < matrix_.nonZeros(); ++k)
        values[k] = hValues_[k] - (w * w) * pValues_[k];
    }
  }

  const SparseMatrix &matrix() const { return matrix_; }

private:
  template <typename Positions>
  void copyValues(const SparseMatrix &part, std::vector<Complex> &values,
                  const Positions &positions) {
    for (Eigen::Index col = 0; col < part.outerSize(); ++col)
      for (typename SparseMatrix::InnerIterator entry(part, col); entry;
           ++entry)
        values.at(positions.at({entry.row(), col})) = entry.value();
  }

  SparseMatrix matrix_;
  std::vector<Complex> hValues_, pValues_, haValues_;
  bool attenuation_;
};

#ifdef DSPECM1D_ENABLE_PROFILING
// Eigen::SparseLU::compute() combines symbolic analysis and numerical
// factorization.  The profiling-only route below expands that same operation
// so the two costs can be measured independently; ordinary builds retain the
// established compute()/factorize() calls below.
template <typename Solver, typename Matrix>
void profileEigenFactorizeOrCompute(Solver &solver, const Matrix &matrix,
                                    int idxn, int nskip,
                                    profiling::Context &profile,
                                    profiling::Mode mode,
                                    bool patternRebuilt) {
  if ((idxn % nskip) == 0) {
    {
      profiling::Scope analyzeProfile(
          profile, profiling::Category::analyze_pattern, mode);
      solver.analyzePattern(matrix);
    }
    profile.countAnalyzePattern(patternRebuilt);
    {
      profiling::Scope factorProfile(profile, profiling::Category::factorization,
                                     mode);
      solver.factorize(matrix);
    }
    profile.countNumericalFactorize();
    profile.countCompute();
  } else {
    profiling::Scope factorProfile(profile, profiling::Category::factorization,
                                   mode);
    solver.factorize(matrix);
    profile.countNumericalFactorize();
    profile.countFactorize(false);
  }
}
#endif

#ifdef DSPECM1D_ENABLE_LAPACK_BAND_SOLVER
using LapackMultiSemSolver = LapackBandSolver;
#endif

} // namespace detail

inline Eigen::MatrixXcd
SparseFSpec::spectra(InputParametersNew &paramsNew) {
  auto &params = paramsNew.inputParameters();
  return spectra(paramsNew.freqFull(), paramsNew.earthModel(), paramsNew.cmt(),
                 params, paramsNew.nq(), paramsNew.srInfo(),
                 params.relative_error(), paramsNew.solverBackend());
}

template <class model1d>
Eigen::MatrixXcd
SparseFSpec::spectra(SpectraSolver::FreqFull &myff, model1d &inp_model,
                     SourceInfo::EarthquakeCMT &cmt, InputParameters &params,
                     int NQ, SRInfo &srInfo, double relerr,
                     SolverBackend backend) {
  using Complex = std::complex<double>;
  using MatrixC = Eigen::MatrixXcd;
  using SparseMatrixC = Eigen::SparseMatrix<Complex>;
  const bool useLapack = backend == SolverBackend::LapackBandLU;
#ifndef DSPECM1D_ENABLE_LAPACK_BAND_SOLVER
  if (useLapack)
    throw std::runtime_error(
        "LapackBandLU was requested, but DSpecM1D was built without LAPACK support.");
#endif
  detail::profiling::Context profile;
  profile.activate();
#ifdef DSPECM1D_ENABLE_PROFILING
  const auto profileStart = std::chrono::steady_clock::now();
#endif
  Timer timer1;

  auto vecW = myff.w();
  SpecConstants sc(myff.ep(), inp_model.TREF());
  auto myi = sc.myi;
  auto ieps = sc.ieps;
  auto w0 = sc.w0;
  auto twodivpi = sc.twodivpi;
  auto twopid = sc.twopid;

  auto numRec = params.num_receivers();
  // Multi-SEM keeps its existing internally derived cadence.  Single-SEM
  // requests receive nskip through SpectraRunContext instead.
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
#ifdef DSPECM1D_ENABLE_PROFILING
  profile.addTime(
      detail::profiling::Category::unclassified,
      detail::profiling::Mode::all,
      std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                    profileStart)
          .count());
#endif

  MatrixC vecRaw = MatrixC::Zero(3 * params.num_receivers(), vecW.size());

  // std::cout << "Just before max step printout\n";
  for (auto s : maxSteps)
    std::cout << "Max step: " << s << "\n";
  // std::cout << "Just after max step printout\n\n";

  // build one SEM per chunk
  std::vector<Full1D::SEM> sems;
  {
    detail::profiling::Scope profileMesh(profile, detail::profiling::Category::sem_mesh);
    for (int idx = 0; idx < numChunks; ++idx) {
      sems.emplace_back(inp_model, maxSteps[idx], NQ, params.lmax());
      profile.countSem();
    }
  }

  timer1.start();
  detail::EigenMultiSemSolver eigenSolver, eigenSolver1;
#ifdef DSPECM1D_ENABLE_LAPACK_BAND_SOLVER
  detail::LapackMultiSemSolver lapackSolver, lapackSolver1;
#endif

  int lmin = params.lmin();
  int lmax = params.lmax();
  ParamRedInfo paramInfo(params, lmax);

#ifdef DSPECM1D_ENABLE_PROFILING
  auto recordSystem = [&](const SparseMatrixC &matrix) {
    Eigen::Index kl = 0, ku = 0;
    for (Eigen::Index column = 0; column < matrix.outerSize(); ++column)
      for (SparseMatrixC::InnerIterator entry(matrix, column); entry;
           ++entry) {
        kl = std::max(kl, entry.row() - entry.col());
        ku = std::max(ku, entry.col() - entry.row());
      }
    profile.countFrequencySystem(static_cast<long>(matrix.rows()),
                                 static_cast<long>(matrix.nonZeros()),
                                 static_cast<long>(kl), static_cast<long>(ku));
  };
#endif

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
      SparseMatrixC keR;
      SparseMatrixC inR;
      SparseMatrixC keRAtten;
      {
        detail::profiling::Scope profileBase(profile, detail::profiling::Category::base_operator,
                                     detail::profiling::Mode::radial);
        keR = sem.hR().cast<Complex>();
        inR = sem.pR().cast<Complex>();
        keRAtten = sem.hRa().cast<Complex>();
      }
      {
        detail::profiling::Scope compressionProfile(profile, detail::profiling::Category::compression,
                                            detail::profiling::Mode::radial);
        keR.makeCompressed();
        inR.makeCompressed();
        keRAtten.makeCompressed();
      }
#ifdef DSPECM1D_ENABLE_LAPACK_BAND_SOLVER
      std::pair<Eigen::Index, Eigen::Index> radialBandwidth{0, 0};
      detail::LapackBandMatrix keRBand, inRBand, keRAttenBand;
      if (useLapack) {
        radialBandwidth = [&] {
          auto result = detail::lapackBandBandwidth(keR);
          const auto pBandwidth = detail::lapackBandBandwidth(inR);
          result.first = std::max(result.first, pBandwidth.first);
          result.second = std::max(result.second, pBandwidth.second);
          if (params.attenuation()) {
            const auto aBandwidth = detail::lapackBandBandwidth(keRAtten);
            result.first = std::max(result.first, aBandwidth.first);
            result.second = std::max(result.second, aBandwidth.second);
          }
          return result;
        }();
        detail::profiling::Scope bandBaseProfile(
            profile, detail::profiling::Category::base_operator,
            detail::profiling::Mode::radial);
        detail::packLapackBandInto(keR, keRBand, radialBandwidth.first,
                                   radialBandwidth.second);
        detail::packLapackBandInto(inR, inRBand, radialBandwidth.first,
                                   radialBandwidth.second);
        if (params.attenuation())
          detail::packLapackBandInto(keRAtten, keRAttenBand,
                                     radialBandwidth.first,
                                     radialBandwidth.second);
      }
#endif
      MatrixC fR;
      MatrixC vecRedZ;
      {
        detail::profiling::Scope projectionProfile(profile, detail::profiling::Category::projection,
                                           detail::profiling::Mode::radial);
        fR = sem.calculateForceRedR(cmt);
        vecRedZ = sem.rvRedZR(params);
      }
      std::vector<std::unique_ptr<detail::FixedPatternFrequencyMatrix>>
          radialFixed(omp_get_max_threads());
#ifdef DSPECM1D_ENABLE_LAPACK_BAND_SOLVER
#pragma omp parallel default(shared) private(eigenSolver, lapackSolver)
#else
#pragma omp parallel default(shared) private(eigenSolver)
#endif
      {
        profile.setMode(detail::profiling::Mode::radial);
        detail::profiling::WorkerScope workerProfile(
            profile, detail::profiling::Mode::radial);
#ifdef DSPECM1D_ENABLE_LAPACK_BAND_SOLVER
        if (useLapack)
          lapackSolver.ensureBandWorkspace(keR.rows(), radialBandwidth.first,
                                           radialBandwidth.second);
#endif
#pragma omp for schedule(dynamic)
        for (int idx = 0; idx < (int) idxChunks[idxChunk].size(); ++idx) {
          double wval = freqChunks[idxChunk][idx];
          int idxw = idxChunks[idxChunk][idx];
          Complex w = wval + ieps;
#ifdef DSPECM1D_ENABLE_PROFILING
          bool patternRebuilt = false;
#endif
#ifdef DSPECM1D_ENABLE_PROFILING
          const auto matrixStart = std::chrono::steady_clock::now();
#endif
          MatrixC vecX;
          const SparseMatrixC *matrixForSolver = nullptr;
#ifdef DSPECM1D_ENABLE_LAPACK_BAND_SOLVER
          if (useLapack) {
            auto &band = lapackSolver.bandWorkspace();
            auto activeA = band.coefficients();
            const auto hActive = keRBand.coefficients();
            const auto pActive = inRBand.coefficients();
            if (params.attenuation()) {
              const auto haActive = keRAttenBand.coefficients();
              activeA = hActive - (w * w) * pActive +
                        attenFactor(wval, w0, twodivpi, myi) * haActive;
            } else {
              activeA = hActive - (w * w) * pActive;
            }
          } else
#endif
          {
            auto &fixed = radialFixed[omp_get_thread_num()];
            if (!fixed) {
              fixed = std::make_unique<detail::FixedPatternFrequencyMatrix>(
                  keR, inR, keRAtten, params.attenuation());
#ifdef DSPECM1D_ENABLE_PROFILING
              patternRebuilt = true;
#endif
              if (auto *active = detail::profiling::Context::active())
                active->countPatternRebuild();
            }
            fixed->update(
                w, params.attenuation()
                       ? attenFactor(wval, w0, twodivpi, myi)
                       : Complex{});
            matrixForSolver = &fixed->matrix();
          }
#ifdef DSPECM1D_ENABLE_PROFILING
          profile.addTime(
              detail::profiling::Category::dynamic_matrix,
              detail::profiling::Mode::radial,
              std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                            matrixStart)
                  .count());
#endif
#ifdef DSPECM1D_ENABLE_PROFILING
          if (!useLapack)
            recordSystem(*matrixForSolver);
#ifdef DSPECM1D_ENABLE_LAPACK_BAND_SOLVER
          else {
            const auto stats = detail::lapackBandActiveStats(
                lapackSolver.bandWorkspace());
            profile.countFrequencySystem(
                static_cast<long>(lapackSolver.bandWorkspace().n),
                std::get<0>(stats), static_cast<long>(std::get<1>(stats)),
                static_cast<long>(std::get<2>(stats)));
          }
#endif
#endif
#ifdef DSPECM1D_ENABLE_LAPACK_BAND_SOLVER
          if (useLapack) {
            lapackSolver.factorize();
            vecX = lapackSolver.solve(fR);
          } else
#endif
          {
#ifdef DSPECM1D_ENABLE_PROFILING
            detail::profileEigenFactorizeOrCompute(
                eigenSolver, *matrixForSolver, 0, 1, profile,
                detail::profiling::Mode::radial, patternRebuilt);
#else
            eigenSolver.compute(*matrixForSolver);
#endif
            {
              detail::profiling::Scope solveProfile(
                  profile, detail::profiling::Category::solve,
                  detail::profiling::Mode::radial);
              vecX = eigenSolver.solve(fR);
              if (auto *active = detail::profiling::Context::active())
                active->countSolve(static_cast<long>(fR.cols()));
            }
          }
          detail::profiling::Scope projectionProfile(profile, detail::profiling::Category::projection,
                                             detail::profiling::Mode::radial);
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
    std::vector<std::unique_ptr<detail::FixedPatternFrequencyMatrix>>
        torFixed(omp_get_max_threads());
    std::vector<std::size_t> torFixedRidx(
        omp_get_max_threads(), std::numeric_limits<std::size_t>::max());
#ifdef DSPECM1D_ENABLE_LAPACK_BAND_SOLVER
#pragma omp parallel default(shared) private(eigenSolver1, lapackSolver1)
#else
#pragma omp parallel default(shared) private(eigenSolver1)
#endif
    {
      profile.setMode(detail::profiling::Mode::toroidal);
      detail::profiling::WorkerScope workerProfile(
          profile, detail::profiling::Mode::toroidal);
#pragma omp for schedule(dynamic)
      for (int idxl = lmin; idxl < lmax + 1; ++idxl) {
        profile.countDegree();
        MatrixC rvVals = srInfo.rvRedTor(idxl);
        for (int idxChunk = 0; idxChunk < numChunks; ++idxChunk) {
          torFixed[omp_get_thread_num()].reset();
          torFixedRidx[omp_get_thread_num()] =
              std::numeric_limits<std::size_t>::max();
          Full1D::SEM &sem = sems[idxChunk];
          auto idxSource = sem.sourceElement(cmt);
          auto recElems = sem.receiverElements(params);
          for (auto idxelem : recElems) {
            if (idxelem < sem.el() || idxelem >= sem.eu()) {
              throw std::runtime_error(
                  "Receiver element out of SEM range in toroidal mode loop.");
            }
          }
          // std::cout << "Debug 0\n";
          auto lowidx = sem.ltgT(recElems[0], 0);
          // std::cout << "Debug 1\n";
          auto upidx = sem.ltgT(recElems.back(), NQ - 1) + 1;
          // std::cout << "Debug 2\n";
          int lenidx = upidx - lowidx;
          auto lentor = sem.ltgT(sem.eu() - 1, NQ - 1) + 1;
          // std::cout << "Debug 3\n";

          SparseMatrixC hTor;
          SparseMatrixC pTor;
          SparseMatrixC hTorAtten;
          {
            detail::profiling::Scope baseProfile(profile, detail::profiling::Category::base_operator,
                                         detail::profiling::Mode::toroidal);
            hTor = sem.hTk(idxl).cast<Complex>();
            pTor = sem.pTk(idxl).cast<Complex>();
            hTorAtten = sem.hTa(idxl).cast<Complex>();
          }
          {
            detail::profiling::Scope compressionProfile(profile, detail::profiling::Category::compression,
                                                detail::profiling::Mode::toroidal);
            hTor.makeCompressed();
            pTor.makeCompressed();
            hTorAtten.makeCompressed();
          }
#ifdef DSPECM1D_ENABLE_LAPACK_BAND_SOLVER
          std::pair<Eigen::Index, Eigen::Index> torBandwidth{0, 0};
          detail::LapackBandMatrix hTorBand, pTorBand, hTorAttenBand;
          if (useLapack) {
            torBandwidth = [&] {
              auto result = detail::lapackBandBandwidth(hTor);
              const auto pBandwidth = detail::lapackBandBandwidth(pTor);
              result.first = std::max(result.first, pBandwidth.first);
              result.second = std::max(result.second, pBandwidth.second);
              if (params.attenuation()) {
                const auto aBandwidth = detail::lapackBandBandwidth(hTorAtten);
                result.first = std::max(result.first, aBandwidth.first);
                result.second = std::max(result.second, aBandwidth.second);
              }
              return result;
            }();
            detail::profiling::Scope bandBaseProfile(
                profile, detail::profiling::Category::base_operator,
                detail::profiling::Mode::toroidal);
            detail::packLapackBandInto(hTor, hTorBand, torBandwidth.first,
                                       torBandwidth.second);
            detail::packLapackBandInto(pTor, pTorBand, torBandwidth.first,
                                       torBandwidth.second);
            if (params.attenuation())
              detail::packLapackBandInto(hTorAtten, hTorAttenBand,
                                         torBandwidth.first,
                                         torBandwidth.second);
          }
#endif

          // std::cout << "Debug 4\n";
          MatrixC fVals;
          // std::cout << "Debug 5\n";
          MatrixC redC;
          MatrixC fBase;
          // std::cout << "Debug 6\n";
          MatrixC rvBase;
          // std::cout << "Debug 7\n";
          {
            detail::profiling::Scope projectionProfile(profile, detail::profiling::Category::projection,
                                               detail::profiling::Mode::toroidal);
            fVals = sem.calculateForceRedCoefficientsT(cmt, idxl, 0.0);
            redC = rvVals * fVals;
            fBase = sem.calculateForceAllT(cmt, idxl);
            rvBase = sem.rvBaseFullT(params, idxl);
          }
          std::vector<int> ridxtor;
          {
            detail::profiling::Scope truncationProfile(profile, detail::profiling::Category::truncation,
                                               detail::profiling::Mode::toroidal);
            ridxtor = SpectralTools::allIndicesTor(
                sem, idxl, freqChunks[idxChunk], idxSource, nskip);
          }
          auto i1 = idxChunks[idxChunk][0];
          auto lenChunk = freqChunks[idxChunk].size();
          MatrixC vecRawL = MatrixC::Zero(3 * params.num_receivers(), lenChunk);
#ifdef DSPECM1D_ENABLE_LAPACK_BAND_SOLVER
          if (useLapack)
            lapackSolver1.ensureBandWorkspace(hTor.rows(), torBandwidth.first,
                                              torBandwidth.second);
#endif
          for (int idx = (int) lenChunk - 1; idx > -1; --idx) {
            auto wval = freqChunks[idxChunk][idx];
            Complex w = wval + ieps;
            std::size_t ridx = ridxtor[idx];
            std::size_t len_ms = lentor - ridx;
#ifdef DSPECM1D_ENABLE_PROFILING
            bool patternRebuilt = false;
            profile.countRidx(static_cast<long>(ridx));
#endif
#ifdef DSPECM1D_ENABLE_PROFILING
            const auto matrixStart = std::chrono::steady_clock::now();
#endif
            const SparseMatrixC *matrixForSolver = nullptr;
#ifdef DSPECM1D_ENABLE_LAPACK_BAND_SOLVER
            Eigen::Index ridxEigen = static_cast<Eigen::Index>(ridx);
            if (useLapack) {
              auto &band = lapackSolver1.bandWorkspace();
              const auto activeSize = static_cast<Eigen::Index>(len_ms);
              auto activeA = band.coefficients().middleCols(
                  ridxEigen, activeSize);
              const auto hActive = hTorBand.coefficients().middleCols(
                  ridxEigen, activeSize);
              const auto pActive = pTorBand.coefficients().middleCols(
                  ridxEigen, activeSize);
              if (params.attenuation()) {
                const auto haActive = hTorAttenBand.coefficients().middleCols(
                    ridxEigen, activeSize);
                activeA = hActive - (w * w) * pActive +
                          attenFactor(wval, w0, twodivpi, myi) * haActive;
              } else {
                activeA = hActive - (w * w) * pActive;
              }
            } else
#endif
            {
              auto &fixed = torFixed[omp_get_thread_num()];
              if (!fixed || torFixedRidx[omp_get_thread_num()] != ridx) {
                SparseMatrixC hActive =
                    hTor.block(ridx, ridx, len_ms, len_ms);
                SparseMatrixC pActive =
                    pTor.block(ridx, ridx, len_ms, len_ms);
                SparseMatrixC haActive =
                    hTorAtten.block(ridx, ridx, len_ms, len_ms);
                fixed = std::make_unique<detail::FixedPatternFrequencyMatrix>(
                    hActive, pActive, haActive, params.attenuation());
                torFixedRidx[omp_get_thread_num()] = ridx;
#ifdef DSPECM1D_ENABLE_PROFILING
                patternRebuilt = true;
#endif
                if (auto *active = detail::profiling::Context::active())
                  active->countPatternRebuild();
              }
              fixed->update(
                  w, params.attenuation()
                         ? attenFactor(wval, w0, twodivpi, myi)
                         : Complex{});
              matrixForSolver = &fixed->matrix();
            }
#ifdef DSPECM1D_ENABLE_PROFILING
            profile.addTime(
                detail::profiling::Category::dynamic_matrix,
                detail::profiling::Mode::toroidal,
                std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                              matrixStart)
                    .count());
#endif
#ifdef DSPECM1D_ENABLE_PROFILING
            if (!useLapack) {
              recordSystem(*matrixForSolver);
#ifdef DSPECM1D_ENABLE_LAPACK_BAND_SOLVER
            } else {
              const auto stats = detail::lapackBandActiveStats(
                  lapackSolver1.bandWorkspace(), ridxEigen);
              profile.countFrequencySystem(
                  static_cast<long>(len_ms), std::get<0>(stats),
                  static_cast<long>(std::get<1>(stats)),
                  static_cast<long>(std::get<2>(stats)));
#endif
            }
#endif
            auto fRed = fBase.block(ridx, 0, len_ms, fBase.cols());
            MatrixC vecSol;
#ifdef DSPECM1D_ENABLE_LAPACK_BAND_SOLVER
            if (useLapack) {
              lapackSolver1.factorize(ridxEigen);
              vecSol = lapackSolver1.solve(fRed);
            } else
#endif
            {
#ifdef DSPECM1D_ENABLE_PROFILING
              detail::profileEigenFactorizeOrCompute(
                  eigenSolver1, *matrixForSolver,
                  (int) lenChunk - idx - 1, nskip, profile,
                  detail::profiling::Mode::toroidal, patternRebuilt);
#else
              factorizeOrCompute(eigenSolver1, *matrixForSolver,
                                 (int) lenChunk - idx - 1, nskip);
#endif
              detail::profiling::Scope solveProfile(profile, detail::profiling::Category::solve,
                                            detail::profiling::Mode::toroidal);
              vecSol = eigenSolver1.solve(fRed);
              profile.countSolve(static_cast<long>(fRed.cols()));
            }
            auto lidx = lowidx - ridx;
            detail::profiling::Scope projectionProfile(profile, detail::profiling::Category::projection,
                                               detail::profiling::Mode::toroidal);
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
    std::vector<std::unique_ptr<detail::FixedPatternFrequencyMatrix>>
        sphFixed(omp_get_max_threads());
    std::vector<std::size_t> sphFixedRidx(
        omp_get_max_threads(), std::numeric_limits<std::size_t>::max());
#ifdef DSPECM1D_ENABLE_LAPACK_BAND_SOLVER
#pragma omp parallel default(shared) private(eigenSolver1, lapackSolver1)
#else
#pragma omp parallel default(shared) private(eigenSolver1)
#endif
    {
      profile.setMode(detail::profiling::Mode::spheroidal);
      detail::profiling::WorkerScope workerProfile(
          profile, detail::profiling::Mode::spheroidal);
#pragma omp for schedule(dynamic)
      for (int idxl = lmin; idxl < lmax + 1; ++idxl) {
        profile.countDegree();
        MatrixC rvVals = srInfo.rvRedSph(idxl);
        for (int idxChunk = 0; idxChunk < numChunks; ++idxChunk) {
          sphFixed[omp_get_thread_num()].reset();
          sphFixedRidx[omp_get_thread_num()] =
              std::numeric_limits<std::size_t>::max();
          Full1D::SEM &sem = sems[idxChunk];
          auto idxSource = sem.sourceElement(cmt);
          auto recElems = sem.receiverElements(params);
          auto lowidx = sem.ltgS(0, recElems[0], 0);
          auto upidx = sem.ltgS(1, recElems.back(), NQ - 1);
          int lenidx = upidx - lowidx + 1;
          auto lensph = sem.ltgS(2, sem.mesh().NE() - 1, NQ - 1) + 1;
          SparseMatrixC pS;
          SparseMatrixC hS;
          SparseMatrixC hSa;
          {
            detail::profiling::Scope baseProfile(profile, detail::profiling::Category::base_operator,
                                         detail::profiling::Mode::spheroidal);
            pS = sem.pS(idxl).cast<Complex>();
            hS = sem.hS(idxl).cast<Complex>();
            hSa = sem.hSa(idxl).cast<Complex>();
          }
          {
            detail::profiling::Scope compressionProfile(profile, detail::profiling::Category::compression,
                                                detail::profiling::Mode::spheroidal);
            pS.makeCompressed();
            hS.makeCompressed();
            hSa.makeCompressed();
          }
#ifdef DSPECM1D_ENABLE_LAPACK_BAND_SOLVER
          std::pair<Eigen::Index, Eigen::Index> sphBandwidth{0, 0};
          detail::LapackBandMatrix hSBand, pSBand, hSaBand;
          if (useLapack) {
            sphBandwidth = [&] {
              auto result = detail::lapackBandBandwidth(hS);
              const auto pBandwidth = detail::lapackBandBandwidth(pS);
              result.first = std::max(result.first, pBandwidth.first);
              result.second = std::max(result.second, pBandwidth.second);
              if (params.attenuation()) {
                const auto aBandwidth = detail::lapackBandBandwidth(hSa);
                result.first = std::max(result.first, aBandwidth.first);
                result.second = std::max(result.second, aBandwidth.second);
              }
              return result;
            }();
            detail::profiling::Scope bandBaseProfile(
                profile, detail::profiling::Category::base_operator,
                detail::profiling::Mode::spheroidal);
            detail::packLapackBandInto(hS, hSBand, sphBandwidth.first,
                                       sphBandwidth.second);
            detail::packLapackBandInto(pS, pSBand, sphBandwidth.first,
                                       sphBandwidth.second);
            if (params.attenuation())
              detail::packLapackBandInto(hSa, hSaBand, sphBandwidth.first,
                                         sphBandwidth.second);
          }
#endif
          MatrixC fVals;
          MatrixC redC;
          MatrixC fBase;
          MatrixC rvBase;
          {
            detail::profiling::Scope projectionProfile(profile, detail::profiling::Category::projection,
                                               detail::profiling::Mode::spheroidal);
            fVals = sem.calculateForceRedCoefficients(cmt, idxl, 0.0);
            redC = rvVals * fVals;
            fBase = sem.calculateForceAll(cmt, idxl);
            rvBase = sem.rvBaseFull(params, idxl);
          }
          std::vector<int> ridxsph;
          {
            detail::profiling::Scope truncationProfile(profile, detail::profiling::Category::truncation,
                                               detail::profiling::Mode::spheroidal);
            ridxsph = SpectralTools::allIndicesSph(
                sem, idxl, freqChunks[idxChunk], idxSource, nskip);
          }
          auto i1 = idxChunks[idxChunk][0];
          auto lenChunk = freqChunks[idxChunk].size();
          MatrixC vecRawL = MatrixC::Zero(3 * params.num_receivers(), lenChunk);
#ifdef DSPECM1D_ENABLE_LAPACK_BAND_SOLVER
          if (useLapack)
            lapackSolver1.ensureBandWorkspace(hS.rows(), sphBandwidth.first,
                                              sphBandwidth.second);
#endif
          for (int idx = (int) lenChunk - 1; idx > -1; --idx) {
            auto wval = freqChunks[idxChunk][idx];
            std::size_t ridx = ridxsph[idx];
            std::size_t len_ms = lensph - ridx;
#ifdef DSPECM1D_ENABLE_PROFILING
            bool patternRebuilt = false;
            profile.countRidx(static_cast<long>(ridx));
#endif
            Complex w = wval + ieps;
#ifdef DSPECM1D_ENABLE_PROFILING
            const auto matrixStart = std::chrono::steady_clock::now();
#endif
            const SparseMatrixC *matrixForSolver = nullptr;
#ifdef DSPECM1D_ENABLE_LAPACK_BAND_SOLVER
            const auto ridxEigen = static_cast<Eigen::Index>(ridx);
            if (useLapack) {
              auto &band = lapackSolver1.bandWorkspace();
              const auto activeSize = static_cast<Eigen::Index>(len_ms);
              auto activeA = band.coefficients().middleCols(
                  ridxEigen, activeSize);
              const auto hActive = hSBand.coefficients().middleCols(
                  ridxEigen, activeSize);
              const auto pActive = pSBand.coefficients().middleCols(
                  ridxEigen, activeSize);
              if (params.attenuation()) {
                const auto haActive = hSaBand.coefficients().middleCols(
                    ridxEigen, activeSize);
                activeA = hActive - (w * w) * pActive +
                          attenFactor(wval, w0, twodivpi, myi) * haActive;
              } else {
                activeA = hActive - (w * w) * pActive;
              }
            } else
#endif
            {
              auto &fixed = sphFixed[omp_get_thread_num()];
              if (!fixed || sphFixedRidx[omp_get_thread_num()] != ridx) {
                SparseMatrixC hActive =
                    hS.block(ridx, ridx, len_ms, len_ms);
                SparseMatrixC pActive =
                    pS.block(ridx, ridx, len_ms, len_ms);
                SparseMatrixC haActive =
                    hSa.block(ridx, ridx, len_ms, len_ms);
                fixed = std::make_unique<detail::FixedPatternFrequencyMatrix>(
                    hActive, pActive, haActive, params.attenuation());
                sphFixedRidx[omp_get_thread_num()] = ridx;
#ifdef DSPECM1D_ENABLE_PROFILING
                patternRebuilt = true;
#endif
                if (auto *active = detail::profiling::Context::active())
                  active->countPatternRebuild();
              }
              fixed->update(
                  w, params.attenuation()
                         ? attenFactor(wval, w0, twodivpi, myi)
                         : Complex{});
              matrixForSolver = &fixed->matrix();
            }
#ifdef DSPECM1D_ENABLE_PROFILING
            profile.addTime(
                detail::profiling::Category::dynamic_matrix,
                detail::profiling::Mode::spheroidal,
                std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                              matrixStart)
                    .count());
#endif
#ifdef DSPECM1D_ENABLE_PROFILING
            if (!useLapack) {
              recordSystem(*matrixForSolver);
#ifdef DSPECM1D_ENABLE_LAPACK_BAND_SOLVER
            } else {
              const auto stats = detail::lapackBandActiveStats(
                  lapackSolver1.bandWorkspace(), ridxEigen);
              profile.countFrequencySystem(
                  static_cast<long>(len_ms), std::get<0>(stats),
                  static_cast<long>(std::get<1>(stats)),
                  static_cast<long>(std::get<2>(stats)));
#endif
            }
#endif
#ifdef DSPECM1D_ENABLE_PROFILING
            const auto rhsStart = std::chrono::steady_clock::now();
#endif
            MatrixC fRed = fBase.block(ridx, 0, len_ms, fBase.cols());
#ifdef DSPECM1D_ENABLE_PROFILING
            profile.addTime(
                detail::profiling::Category::truncation,
                detail::profiling::Mode::spheroidal,
                std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                              rhsStart)
                    .count());
#endif
            MatrixC vecSol;
#ifdef DSPECM1D_ENABLE_LAPACK_BAND_SOLVER
            if (useLapack) {
              lapackSolver1.factorize(ridxEigen);
              vecSol = lapackSolver1.solve(fRed);
            } else
#endif
            {
#ifdef DSPECM1D_ENABLE_PROFILING
              detail::profileEigenFactorizeOrCompute(
                  eigenSolver1, *matrixForSolver,
                  (int) lenChunk - idx - 1, nskip, profile,
                  detail::profiling::Mode::spheroidal, patternRebuilt);
#else
              factorizeOrCompute(eigenSolver1, *matrixForSolver,
                                 (int) lenChunk - idx - 1, nskip);
#endif
              detail::profiling::Scope solveProfile(profile, detail::profiling::Category::solve,
                                            detail::profiling::Mode::spheroidal);
              vecSol = eigenSolver1.solve(fRed);
              profile.countSolve(static_cast<long>(fRed.cols()));
            }
            auto lidx = lowidx - ridx;
            detail::profiling::Scope projectionProfile(profile, detail::profiling::Category::projection,
                                               detail::profiling::Mode::spheroidal);
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

  {
    detail::profiling::Scope outputProfile(
        profile, detail::profiling::Category::projection,
        detail::profiling::Mode::all);
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
  }
#ifdef DSPECM1D_ENABLE_PROFILING
  profile.addTime(
      detail::profiling::Category::total, detail::profiling::Mode::all,
      std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                    profileStart)
          .count());
  detail::profiling::publish(profile.finish());
#endif
  profile.deactivate();
  return vecRaw;
}

}   // namespace SPARSESPEC
#endif   // DSPECM1D_FULL_SPEC_MULTI_SEM_H
