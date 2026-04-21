#ifndef DSPECM1D_FULL_SPEC_H
#define DSPECM1D_FULL_SPEC_H

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

/**
 * @brief Public orchestration entry point for sparse frequency-domain
 * seismogram synthesis.
 *
 * `SparseFSpec` is the preferred high-level solver interface for release
 * users. The `InputParametersNew` overloads are the recommended path; the
 * legacy overloads remain available for backward compatibility with older
 * workflows and paper examples.
 */
class SparseFSpec {
public:
  SparseFSpec() {};
  ~SparseFSpec() {};

  /// Preferred release-facing overload taking a fully prepared workflow
  /// context.
  Eigen::MatrixXcd spectra(InputParametersNew &);
  /// Preferred release-facing overload that reuses a prepared SEM object.
  Eigen::MatrixXcd spectra(InputParametersNew &, Full1D::SEM &);
  /// Preferred release-facing overload with reversed convenience arguments.
  Eigen::MatrixXcd spectra(Full1D::SEM &, InputParametersNew &);
  /// Preferred internal/reuse overload using an explicit run context.
  Eigen::MatrixXcd spectra(const SpectraRunContext &, Full1D::SEM &);

  /// Legacy overload retained for compatibility with existing workflows.
  template <class model1d>
  Eigen::MatrixXcd spectra(SpectraSolver::FreqFull &, Full1D::SEM &, model1d &,
                           SourceInfo::EarthquakeCMT &, InputParameters &,
                           int = 10);

  /// Legacy overload retained for compatibility with existing workflows.
  template <class model1d>
  Eigen::MatrixXcd spectra(SpectraSolver::FreqFull &, model1d &,
                           SourceInfo::EarthquakeCMT &, InputParameters &, int,
                           SRInfo &, double = 1e-4);

private:
};

}   // namespace SPARSESPEC

#include "FullSpecSingleSem.h"
#include "FullSpecMultiSem.h"

#endif   // DSPECM1D_FULL_SPEC_H
