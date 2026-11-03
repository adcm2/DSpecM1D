#ifndef DSPECM1D_PHASE_A_STANDARD_FORWARD_HPP
#define DSPECM1D_PHASE_A_STANDARD_FORWARD_HPP

#include "DSpecM1D/src/SpectraRunContext.h"
#include "DSpecM1D/src/StartElement.h"

// Isolated Phase A copy of the pinned FullSpecSingleSem implementation.
// The only algorithm switch fills source-start indices with zero for the
// explicitly untruncated diagnostic submode; all matrix/force/receiver/output
// expressions remain the pinned standard implementation.
namespace SPARSESPEC {
struct PhaseACutIndexRecord {
  int family=-1; // 0 radial, 1 toroidal, 2 spheroidal
  int degree=-1;
  int bin=-1;
  int formulation=0; // spheroidal: 1 Cowling, 2 full; otherwise 0
  int startDof=-1;
  int retainedRows=0;
  int fullRows=0;
  int nskip=0;
};

inline Eigen::MatrixXcd
phaseAStandardForward(const SpectraRunContext &request, Full1D::SEM &sem,
                     bool forceCowling, bool untruncated,
                     const std::vector<std::complex<double>> &exactZ,
                     std::vector<PhaseACutIndexRecord> *cutIndexAudit=nullptr) {
  using Complex = std::complex<double>;
  using MatrixC = Eigen::MatrixXcd;
  using SparseMatrixC = Eigen::SparseMatrix<Complex>;
  using SparseLUType = EquilibratedSparseLU<Complex>;

  auto &myff = request.freqFull();
  auto &cmt = request.cmt();
  auto &params = request.params();
  const int nskip = request.nskip();

  if (!sem.sourceMatches(cmt))
    throw std::invalid_argument(
        "Prepared SEM is not source-aligned to the current source depth; construct it with InputParametersNew or the source-aware SEM overload.");

  if (params.attenuation())
    validateAttenuationReferencePeriod(request.tref());

  auto vecW = myff.w();
  SpecConstants sc(myff.ep(), request.tref());
  if (exactZ.size()!=static_cast<std::size_t>(myff.i2()-myff.i1()))
    throw std::invalid_argument("Phase A exact frequency-grid size differs from standard request");
  if (exactZ.empty())throw std::invalid_argument("Phase A exact frequency grid is empty");
  for(std::size_t k=0;k<exactZ.size();++k){
    if(!std::isfinite(exactZ[k].real())||!std::isfinite(exactZ[k].imag())||
       exactZ[k].imag()!=sc.ieps.imag())
      throw std::invalid_argument("Phase A exact frequency damping differs from the shared case");
    vecW[static_cast<std::size_t>(myff.i1())+k]=exactZ[k].real();
  }
  std::vector<double> toroidalFrequencies(
      vecW.begin()+myff.i1(),vecW.begin()+myff.i2());
  auto myi = sc.myi;
  auto ieps = sc.ieps;
  auto w0 = sc.w0;

  MatrixC vecRaw = MatrixC::Zero(3 * params.num_receivers(), vecW.size());
  SparseLUType solver, solver1;

  std::vector<mathprofile::WorkerData> profileWorkers(
      static_cast<std::size_t>(std::max(1,omp_get_max_threads())));
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

  const std::size_t frequencyCount=static_cast<std::size_t>(myff.i2()-myff.i1());
  const std::size_t degreeCount=lmax>=lmin?static_cast<std::size_t>(lmax-lmin+1):0;
  const std::size_t radialCutOffset=0;
  const std::size_t toroidalCutOffset=inc_rad?frequencyCount:0;
  const std::size_t spheroidalCutOffset=toroidalCutOffset+
      (inc_tor?degreeCount*frequencyCount:0);
  if(cutIndexAudit)cutIndexAudit->assign(
      spheroidalCutOffset+(inc_sph?degreeCount*2*frequencyCount:0),{});

  // Validate every receiver before entering any OpenMP region.
  const auto receiverElems = sem.receiverElements(params, inc_tor);

  auto numRec = params.num_receivers();
  auto idxSource = sem.sourceElement(cmt);

  // radials
  if (inc_rad) {
    mathprofile::Scope radialSetup("phasea_standard_radial_static_operator_and_receiver_setup");
    const auto &recElems = receiverElems;
    auto lowidx = sem.ltgR(0, recElems[0], 0);
    auto upidx = sem.ltgR(1, recElems.back(), NQ - 1);
    int lenidx = upidx - lowidx + 1;
    SparseMatrixC keR = sem.hR().cast<Complex>();
    SparseMatrixC inR = sem.pR().cast<Complex>();
    SparseMatrixC keRAtten = sem.hRa().cast<Complex>();
    keR.makeCompressed();
    inR.makeCompressed();
    keRAtten.makeCompressed();
    std::vector<MatrixC> vecRvZ;
    for (int idxr = 0; idxr < numRec; ++idxr)
      vecRvZ.push_back(sem.rvZR(params, idxr).block(lowidx, 0, lenidx, 1));
    radialSetup.stop();

    mathprofile::Scope radialRegionWall("phasea_standard_radial_parallel_region_wall");
#pragma omp parallel default(shared) private(solver)
    {
      mathprofile::WorkerCapture workerCapture(true,profileWorkers.at(static_cast<std::size_t>(omp_get_thread_num())));
      mathprofile::Scope workerRegion("phasea_standard_radial_worker_region");
#pragma omp for schedule(dynamic, 10)
      for (int idx = myff.i1(); idx < myff.i2(); ++idx) {
        mathprofile::Scope sourceTreatment("phasea_standard_radial_frequency_source_treatment");
        Complex w = vecW[idx] + ieps;
        const Complex sourceAttenuation = params.attenuation()
                                              ? attenFactor(w, w0)
                                              : Complex(0.0, 0.0);
        auto source = sem.sourceTreatedRadial(cmt, w, sourceAttenuation);
        sourceTreatment.stop();
        mathprofile::Scope dynamicOperator("phasea_standard_radial_dynamic_operator_assembly");
        SparseMatrixC wR = -w * w * inR + keR;
        if (params.attenuation())
          wR += attenFactor(w, w0) * keRAtten;
        wR.makeCompressed();
        dynamicOperator.stop();
        solver.compute(wR);
        MatrixC vecX = solver.solve(source.force);
        if(cutIndexAudit){
          const int fullRows=static_cast<int>(source.force.rows());
          (*cutIndexAudit)[radialCutOffset+static_cast<std::size_t>(idx-myff.i1())]=
              {0,0,idx,0,0,fullRows,fullRows,nskip};
        }
        mathprofile::Scope receiverProjection("phasea_standard_radial_receiver_lift_and_projection");
        MatrixC liftOutput = sem.sourceLiftReceiverRadial(
            params, source.lift, source.sourceElement);
        for (int idxr = 0; idxr < numRec; ++idxr) {
#pragma omp critical(torvecadd)
          {
            vecRaw(3 * idxr, idx) +=
                vecRvZ[idxr]
                    .cwiseProduct(vecX.block(lowidx, 0, lenidx, 1))
                    .sum() + liftOutput(3 * idxr, 0);
          }
        }
        receiverProjection.stop();
      }
      workerRegion.stop();
    }
    radialRegionWall.stop();
    for(const auto &data:profileWorkers)mathprofile::mergeWorker(data,"standard_worker.");
  }

  // toroidals
  if (inc_tor) {
    const auto &recElems = receiverElems;
    auto lowidx = sem.ltgT(recElems[0], 0);
    auto upidx = sem.ltgT(recElems.back(), NQ - 1) + 1;
    int lenidx = upidx - lowidx;
    auto lentor = sem.ltgT(sem.eu() - 1, NQ - 1) + 1;
    mathprofile::Scope toroidalRegionWall("phasea_standard_toroidal_parallel_region_wall");
#pragma omp parallel default(shared) private(solver1)
    {
      mathprofile::WorkerCapture workerCapture(true,profileWorkers.at(static_cast<std::size_t>(omp_get_thread_num())));
      mathprofile::Scope workerRegion("phasea_standard_toroidal_worker_region");
#pragma omp for schedule(dynamic)
      for (int idxl = lmin; idxl < lmax + 1; ++idxl) {
        mathprofile::Scope staticDegree("phasea_standard_toroidal_degree_static_operators_and_coefficients");
        SparseMatrixC hTor = sem.hTk(idxl).cast<Complex>();
        SparseMatrixC pTor = sem.pTk(idxl).cast<Complex>();
        SparseMatrixC hTorAtten = sem.hTa(idxl).cast<Complex>();
        hTor.makeCompressed();
        pTor.makeCompressed();
        hTorAtten.makeCompressed();
        MatrixC rvVals = paramInfo.rvFullTor(idxl);
        MatrixC fVals = sem.calculateForceCoefficientsT(cmt, idxl);
        MatrixC redC = rvVals * fVals;
        MatrixC rvBase = sem.rvBaseFullT(params, idxl);
        staticDegree.stop();
        std::vector<int> ridxtor;
        mathprofile::Scope cutoffSelection("phasea_standard_toroidal_cutoff_policy_selection");
        if (untruncated)
          ridxtor.assign(static_cast<std::size_t>(myff.i2()-myff.i1()),0);
        else
          ridxtor=SpectralTools::allIndicesTor(
              sem, idxl, toroidalFrequencies,
              idxSource, nskip);
        cutoffSelection.stop();
        for (int idx = myff.i2() - 1; idx > myff.i1() - 1; --idx) {
          std::size_t ridx = ridxtor[idx - myff.i1()];
          std::size_t len_ms = lentor - ridx;
          if(cutIndexAudit){
            const std::size_t slot=toroidalCutOffset+
                static_cast<std::size_t>(idxl-lmin)*frequencyCount+
                static_cast<std::size_t>(idx-myff.i1());
            (*cutIndexAudit)[slot]={1,idxl,idx,0,static_cast<int>(ridx),
                static_cast<int>(len_ms),static_cast<int>(lentor),nskip};
          }
          Complex w = vecW[idx] + ieps;
          mathprofile::Scope sourceTreatment("phasea_standard_toroidal_frequency_source_treatment");
          const Complex sourceAttenuation = params.attenuation()
                                                ? attenFactor(w, w0)
                                                : Complex(0.0, 0.0);
          auto source = sem.sourceTreatedToroidal(
              cmt, idxl, w, sourceAttenuation);
          sourceTreatment.stop();
          mathprofile::Scope dynamicOperator("phasea_standard_toroidal_dynamic_operator_assembly");
          SparseMatrixC matWTorRed =
              hTor.block(ridx, ridx, len_ms, len_ms) -
              w * w * pTor.block(ridx, ridx, len_ms, len_ms);
          if (params.attenuation())
            matWTorRed += attenFactor(w, w0) *
                          hTorAtten.block(ridx, ridx, len_ms, len_ms);
          matWTorRed.makeCompressed();
          dynamicOperator.stop();
          auto fRed = source.force.block(ridx, 0, len_ms, source.force.cols());
          factorizeOrCompute(solver1, matWTorRed, myff.i2() - idx - 1, nskip);
          MatrixC vecSol = solver1.solve(fRed);
          mathprofile::Scope receiverProjection("phasea_standard_toroidal_receiver_lift_and_projection");
          MatrixC receiverLift = sem.sourceLiftReceiverToroidal(
              params, idxl, source.lift, source.sourceElement);
          auto lidx = lowidx - ridx;
#pragma omp critical(torvecadd)
          {
            vecRaw.col(idx) +=
                redC.cwiseProduct(rvBase * vecSol.block(lidx, 0, lenidx, 2) +
                                  receiverLift)
                    .rowwise()
                    .sum();
          }
          receiverProjection.stop();
        }
      }
      workerRegion.stop();
    }
    toroidalRegionWall.stop();
    for(const auto &data:profileWorkers)mathprofile::mergeWorker(data,"standard_worker.");
  }

  // spheroidals
  if (inc_sph) {
    const auto &recElems = receiverElems;
    mathprofile::Scope cutoffBoundary("phasea_standard_spheroidal_cowling_cutoff_boundary_selection");
    int cowlingFirstIndex = forceCowling ? myff.i1() : myff.i2();
    const double cutoffMhz = params.cowlingFrequencyMhz();
    if (!forceCowling && cutoffMhz > 0.0) {
      for (int idx = myff.i1(); idx < myff.i2(); ++idx) {
        const double frequencyMhz = myff.f(idx) * 1000.0 / myff.timeNorm();
        if (frequencyMhz >= cutoffMhz) {
          cowlingFirstIndex = idx;
          break;
        }
      }
    }
    cutoffBoundary.stop();
    mathprofile::Scope spheroidalRegionWall("phasea_standard_spheroidal_parallel_region_wall");
#pragma omp parallel default(shared) private(solver1)
    {
      mathprofile::WorkerCapture workerCapture(true,profileWorkers.at(static_cast<std::size_t>(omp_get_thread_num())));
      mathprofile::Scope workerRegion("phasea_standard_spheroidal_worker_region");
#pragma omp for schedule(dynamic)
      for (int idxl = lmin; idxl < lmax + 1; ++idxl) {
        mathprofile::Scope degreeCoefficients("phasea_standard_spheroidal_degree_angular_and_source_coefficients");
        MatrixC rvVals = paramInfo.rvFullSph(idxl);
        MatrixC fVals = sem.calculateForceCoefficients(cmt, idxl);
        MatrixC redC = rvVals * fVals;
        degreeCoefficients.stop();
        // Separate slices restart each formulation's truncation and LU cadence.
        for (int region = 1; region >= 0; --region) {
          const bool useCowling = region == 1;
          const int regionBegin = useCowling ? cowlingFirstIndex : myff.i1();
          const int regionEnd = useCowling ? myff.i2() : cowlingFirstIndex;
          if (regionBegin >= regionEnd)
            continue;
          mathprofile::Scope regionSetup("phasea_standard_spheroidal_formulation_static_operators_and_receiver_basis");
          const auto lowidx = useCowling ? sem.ltgSC(0, recElems[0], 0)
                                         : sem.ltgS(0, recElems[0], 0);
          const auto upidx = useCowling
                                 ? sem.ltgSC(1, recElems.back(), NQ - 1)
                                 : sem.ltgS(1, recElems.back(), NQ - 1);
          const int lenidx = upidx - lowidx + 1;
          const auto lensph = useCowling
                                  ? sem.ltgSC(1, sem.mesh().NE() - 1, NQ - 1) + 1
                                  : sem.ltgS(2, sem.mesh().NE() - 1, NQ - 1) + 1;
          SparseMatrixC hS =
              (useCowling ? sem.hSC(idxl) : sem.hS(idxl)).cast<Complex>();
          SparseMatrixC pS =
              (useCowling ? sem.pSC(idxl) : sem.pS(idxl)).cast<Complex>();
          SparseMatrixC hSa =
              (useCowling ? sem.hSCa(idxl) : sem.hSa(idxl)).cast<Complex>();
          hS.makeCompressed();
          pS.makeCompressed();
          hSa.makeCompressed();
          MatrixC rvBase = sem.rvBaseFull(params, idxl, useCowling);
          regionSetup.stop();
          std::vector<double> regionW(vecW.begin() + regionBegin,
                                      vecW.begin() + regionEnd);
          std::vector<int> regionRidx;
          mathprofile::Scope cutoffSelection("phasea_standard_spheroidal_frequency_start_index_policy_selection");
          if (untruncated)
            regionRidx.assign(regionW.size(),0);
          else
            regionRidx=SpectralTools::allIndicesSph(
                sem, idxl, regionW, idxSource, nskip, useCowling);
          cutoffSelection.stop();
          SparseLUType solver;
          for (int idx = regionEnd - 1; idx >= regionBegin; --idx) {
            const std::size_t idxRs = regionRidx[idx - regionBegin];
            std::size_t len_ms = lensph - idxRs;
            if(cutIndexAudit){
              const std::size_t slot=spheroidalCutOffset+
                  (static_cast<std::size_t>(idxl-lmin)*2+static_cast<std::size_t>(region))*frequencyCount+
                  static_cast<std::size_t>(idx-regionBegin);
              (*cutIndexAudit)[slot]={2,idxl,idx,region==1?1:2,
                  static_cast<int>(idxRs),static_cast<int>(len_ms),lensph,nskip};
            }
            Complex w = vecW[idx] + ieps;
            mathprofile::Scope sourceTreatment("phasea_standard_spheroidal_frequency_source_treatment");
            const Complex sourceAttenuation = params.attenuation()
                                                  ? attenFactor(w, w0)
                                                  : Complex(0.0, 0.0);
            auto source = sem.sourceTreatedSpheroidal(
                cmt, idxl, w, sourceAttenuation, useCowling);
            sourceTreatment.stop();
            mathprofile::Scope dynamicOperator("phasea_standard_spheroidal_dynamic_operator_assembly");
            SparseMatrixC matSph = hS.block(idxRs, idxRs, len_ms, len_ms) -
                                   w * w *
                                       pS.block(idxRs, idxRs, len_ms, len_ms);
            if (params.attenuation())
              matSph += attenFactor(w, w0) *
                        hSa.block(idxRs, idxRs, len_ms, len_ms);
            matSph.makeCompressed();
            dynamicOperator.stop();
            auto fRed = source.force.block(idxRs, 0, len_ms,
                                           source.force.cols());
            factorizeOrCompute(solver, matSph, regionEnd - idx - 1, nskip);
            MatrixC vecSol = solver.solve(fRed);
            auto lidx = lowidx - idxRs;
            mathprofile::Scope receiverProjection("phasea_standard_spheroidal_receiver_lift_and_projection");
            MatrixC receiverLift = sem.sourceLiftReceiver(
                params, idxl, source.lift, source.sourceElement, useCowling);

#pragma omp critical(sphvecadd)
            {
              vecRaw.col(idx) +=
                  redC.cwiseProduct(
                      rvBase * vecSol.block(lidx, 0, lenidx, 4) +
                      receiverLift)
                      .rowwise()
                      .sum();
            }
            receiverProjection.stop();
          }
        }
      }
      workerRegion.stop();
    }
    spheroidalRegionWall.stop();
    for(const auto &data:profileWorkers)mathprofile::mergeWorker(data,"standard_worker.");
  }

  for (int idx = 0; idx < vecRaw.cols(); ++idx) {
    Complex w = vecW[idx] + ieps;
    vecRaw.col(idx) *= outputFactor(params.output_type(), w, myi);
  }
  return vecRaw;
}
} // namespace SPARSESPEC

#endif
