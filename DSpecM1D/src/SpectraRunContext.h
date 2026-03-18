#ifndef DSPECM1D_SPECTRA_RUN_CONTEXT_H
#define DSPECM1D_SPECTRA_RUN_CONTEXT_H

#include <algorithm>
#include "SourceInfo.h"
#include "InputParser.h"

namespace SpectraSolver {
class FreqFull;
}

namespace SPARSESPEC {

// Transitional request object for spectra calls.
// This keeps legacy InputParameters-based workflows working while enabling
// a single-argument API surface for spectrum requests.
class SpectraRunContext {
public:
  SpectraRunContext(SpectraSolver::FreqFull &freqFull,
                    SourceInfo::EarthquakeCMT &cmt, InputParameters &params,
                    double tref, int nskip = 10)
      : m_freqFull(&freqFull), m_cmt(&cmt), m_params(&params), m_tref(tref),
        m_nskip(std::max(1, nskip)) {}

  template <class model1d>
  SpectraRunContext(SpectraSolver::FreqFull &freqFull,
                    SourceInfo::EarthquakeCMT &cmt, InputParameters &params,
                    const model1d &model, int nskip = 10)
      : SpectraRunContext(freqFull, cmt, params, model.TREF(), nskip) {}

  SpectraSolver::FreqFull &freqFull() const { return *m_freqFull; }
  SourceInfo::EarthquakeCMT &cmt() const { return *m_cmt; }
  InputParameters &params() const { return *m_params; }
  double tref() const { return m_tref; }
  int nskip() const { return m_nskip; }

private:
  SpectraSolver::FreqFull *m_freqFull;
  SourceInfo::EarthquakeCMT *m_cmt;
  InputParameters *m_params;
  double m_tref;
  int m_nskip;
};

}   // namespace SPARSESPEC

#endif   // DSPECM1D_SPECTRA_RUN_CONTEXT_H
