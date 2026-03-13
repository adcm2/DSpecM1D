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
#include "FEM_Preconditioner.h"
#include "BiCGSTABT.h"
#include "ParamRedInfo.h"
#include "SR_Info.h"
#include "SpecHelpers.h"

namespace SPARSESPEC {

class Sparse_F_Spec {
public:
  Sparse_F_Spec() {};
  ~Sparse_F_Spec() {};

  template <class model1d>
  auto Spectra(SpectraSolver::FreqFull &, Full1D::specsem &, model1d &,
               SourceInfo::EarthquakeCMT &, InputParameters &, int = 10);

  template <class model1d>
  auto Spectra(SpectraSolver::FreqFull &, model1d &,
               SourceInfo::EarthquakeCMT &, InputParameters &, int, SRInfo &,
               double = 1e-4);

private:
};

}   // namespace SPARSESPEC

#include "FullSpec_SingleSem.h"
#include "FullSpec_MultiSem.h"

#endif   // FULL_SPEC_GUARD_H