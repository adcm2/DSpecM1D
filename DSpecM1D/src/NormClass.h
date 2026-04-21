#ifndef NORM_CLASS_GUARD_H
#define NORM_CLASS_GUARD_H

#include <cmath>

// Use constexpr for mathematical constants
constexpr double PI = 3.1415926535897932;
constexpr double TWO_PI = 2.0 * PI;
constexpr std::complex<double> I(0.0, 1.0);

template <typename FLOAT> class prem_norm {
public:
  prem_norm() = default;

  FLOAT LengthNorm() const { return m_lengthNorm; }
  FLOAT MassNorm() const { return m_massNorm; }
  FLOAT TimeNorm() const { return m_timeNorm; }

private:
  FLOAT m_lengthNorm = 1000.0;
  FLOAT m_massNorm = 5515.0 * std::pow(m_lengthNorm, 3.0);
  FLOAT m_timeNorm = 1.0 / std::sqrt(PI * 6.67230e-11 * 5515.0);
};

#endif   // NORM_CLASS_GUARD_H
