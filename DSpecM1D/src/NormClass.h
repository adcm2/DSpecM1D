#ifndef NORM_CLASS_GUARD_H
#define NORM_CLASS_GUARD_H

#include <complex>
#include <cmath>

/// Mathematical constant π used throughout the library.
constexpr double PI = 3.1415926535897932;
/// Convenience constant equal to 2π.
constexpr double TWO_PI = 2.0 * PI;
/// Imaginary unit used in spectral-domain expressions.
constexpr std::complex<double> I(0.0, 1.0);

/**
 * @brief Non-dimensionalisation scales for PREM-style Earth models.
 *
 * @tparam FLOAT Floating-point type used for the stored scales.
 *
 * This helper provides the length, mass, and time scales used by DSpecM1D
 * when mapping between dimensional and non-dimensional quantities.
 */
template <typename FLOAT> class prem_norm {
public:
  prem_norm() = default;

  /// Returns the length normalisation scale in metres.
  FLOAT LengthNorm() const { return m_lengthNorm; }
  /// Returns the mass normalisation scale in kilograms.
  FLOAT MassNorm() const { return m_massNorm; }
  /// Returns the time normalisation scale in seconds.
  FLOAT TimeNorm() const { return m_timeNorm; }

private:
  FLOAT m_lengthNorm = 1000.0;
  FLOAT m_massNorm = 5515.0 * std::pow(m_lengthNorm, 3.0);
  FLOAT m_timeNorm = 1.0 / std::sqrt(PI * 6.67230e-11 * 5515.0);
};

#endif   // NORM_CLASS_GUARD_H
