#ifndef FULL_SPEC_GUARD_H
#define FULL_SPEC_GUARD_H

#include <iostream>
#include <PlanetaryModel/All>
#include <DSpecM1D/Timer>
#include "ReadStation.h"
#include "SourceInfo.h"
#include "StartElement.h"
#include "SpectraMaster.h"
#include "ParamInfo.h"
#include <omp.h>
#include "FEMPreconditioner.h"
#include "BiCGSTABT.h"
#include "ParamRedInfo.h"
#include "SRInfo.h"
#include "SpecHelpers.h"
#include "SpectraRunContext.h"
#include "InputParametersNew.h"

namespace SPARSESPEC {

class SparseFSpec {
public:
  SparseFSpec() {};
  ~SparseFSpec() {};

  Eigen::MatrixXcd spectra(InputParametersNew &, Full1D::SEM &);
  Eigen::MatrixXcd spectra(const SpectraRunContext &, Full1D::SEM &);

  template <class model1d>
  auto spectra(SpectraSolver::FreqFull &, Full1D::SEM &, model1d &,
               SourceInfo::EarthquakeCMT &, InputParameters &, int = 10);

  template <class model1d>
  auto spectra(SpectraSolver::FreqFull &, model1d &,
               SourceInfo::EarthquakeCMT &, InputParameters &, int, SRInfo &,
               double = 1e-4);

private:
};

}   // namespace SPARSESPEC

#include "FullSpecSingleSem.h"
#include "FullSpecMultiSem.h"

#endif   // FULL_SPEC_GUARD_H