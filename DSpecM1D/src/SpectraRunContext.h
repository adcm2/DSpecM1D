#ifndef DSPECM1D_SPECTRA_RUN_CONTEXT_H
#define DSPECM1D_SPECTRA_RUN_CONTEXT_H

#include <algorithm>
#include "SourceInfo.h"
#include "InputParser.h"

namespace SpectraSolver {
class FreqFull;
}

namespace SPARSESPEC {

/**
 * @brief Transitional request object for spectra calls.
 *
 * This keeps legacy `InputParameters`-based workflows working while enabling
 * a single-argument API surface for spectrum requests.
 */
class SpectraRunContext {
public:
  /**
   * @brief Constructs a context from explicit frequency, source, and parameter
   * objects.
   *
   * @param freqFull Frequency helper used by the solve.
   * @param cmt Earthquake source object.
   * @param params Legacy parameter object.
   * @param tref Reference period of the active Earth model.
   * @param nskip Re-factorization cadence for repeated solves.
   */
  SpectraRunContext(SpectraSolver::FreqFull &freqFull,
                    SourceInfo::EarthquakeCMT &cmt, InputParameters &params,
                    double tref, int nskip = 10)
      : m_freqFull(&freqFull), m_cmt(&cmt), m_params(&params), m_tref(tref),
        m_nskip(std::max(1, nskip)) {}

  /**
   * @brief Convenience constructor that derives the reference period from a
   * model object exposing `TREF()`.
   */
  template <class model1d>
  SpectraRunContext(SpectraSolver::FreqFull &freqFull,
                    SourceInfo::EarthquakeCMT &cmt, InputParameters &params,
                    const model1d &model, int nskip = 10)
      : SpectraRunContext(freqFull, cmt, params, model.TREF(), nskip) {}

  /// Returns the frequency helper used during the solve.
  SpectraSolver::FreqFull &freqFull() const { return *m_freqFull; }
  /// Returns the source object.
  SourceInfo::EarthquakeCMT &cmt() const { return *m_cmt; }
  /// Returns the legacy parameter object.
  InputParameters &params() const { return *m_params; }
  /// Returns the Earth-model reference period used in attenuation terms.
  double tref() const { return m_tref; }
  /// Returns the solve re-factorization cadence.
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
