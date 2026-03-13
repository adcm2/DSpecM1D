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

  FLOAT LengthNorm() const { return _length_norm; }
  FLOAT MassNorm() const { return _mass_norm; }
  FLOAT TimeNorm() const { return _time_norm; }

private:
  FLOAT _length_norm = 1000.0;
  FLOAT _mass_norm = 5515.0 * std::pow(_length_norm, 3.0);
  FLOAT _time_norm = 1.0 / std::sqrt(PI * 6.67230e-11 * 5515.0);
};

#endif   // NORM_CLASS_GUARD_H