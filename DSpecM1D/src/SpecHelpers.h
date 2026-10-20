#ifndef SPEC_HELPERS_GUARD_H
#define SPEC_HELPERS_GUARD_H

#include <complex>
#include <cmath>
#include <stdexcept>

namespace SPARSESPEC {

/**
 * @brief Flags indicating which mode families should be included in a solve.
 */
struct ModeFlags {
  bool inc_rad = false, inc_tor = false, inc_sph = false;
};

/// Resolves which mode types (radial, toroidal, spheroidal) to include
/// based on the mode type flag and the angular degree range [lmin, lmax].
inline ModeFlags
resolveModeFlags(int mtype, int lmin, int lmax) {
  ModeFlags f;
  if (mtype == 4) {
    f.inc_rad = f.inc_tor = f.inc_sph = true;
  } else if (mtype == 1) {
    f.inc_rad = true;
  } else if (mtype == 2) {
    f.inc_tor = true;
  } else if (mtype == 3) {
    f.inc_sph = true;
  }
  if (lmin > 0)
    f.inc_rad = false;
  if (lmax < 1) {
    f.inc_tor = false;
    f.inc_sph = false;
  }
  return f;
}

/// Returns the complex attenuation correction factor: i + (2/pi)*log(w/w0).
inline std::complex<double>
attenFactor(double w, double w0, double twodivpi, std::complex<double> myi) {
  return myi + twodivpi * std::log(w / w0);
}

/// Validates a reference period only for attenuation-enabled solve paths.
inline void
validateAttenuationReferencePeriod(double tref) {
  if (!std::isfinite(tref) || tref <= 0.0)
    throw std::invalid_argument(
        "Attenuation requires a finite positive reference period.");
  const double w0 = 2.0 * std::acos(-1.0) / tref;
  if (!std::isfinite(w0) || w0 <= 0.0)
    throw std::invalid_argument(
        "Attenuation reference period produces an invalid reference frequency.");
}

/// Returns c(z)=(2/pi)Log(i*z/w0) on the principal branch, for the
/// lower-half-plane solver frequency z=omega-i*epsilon. The existing
/// material correction remains first-order in 1/Q.
inline std::complex<double>
attenFactor(std::complex<double> z, double w0) {
  if (!std::isfinite(z.real()) || !std::isfinite(z.imag()) ||
      z == std::complex<double>{})
    throw std::invalid_argument(
        "Complex attenuation frequency must be finite and nonzero.");
  if (!std::isfinite(w0) || w0 <= 0.0)
    throw std::invalid_argument(
        "Attenuation reference frequency must be finite and positive.");
  const std::complex<double> imaginary{0.0, 1.0};
  return (2.0 / std::acos(-1.0)) * std::log(imaginary * z / w0);
}

/// Returns the output-type multiplier applied to the final spectrum.
/// output_type == 0: displacement  =>  -i/w
/// output_type == 1: velocity      =>   1
/// output_type == 2: acceleration  =>   iw
inline std::complex<double>
outputFactor(int output_type, std::complex<double> w,
             std::complex<double> myi) {
  if (output_type == 0)
    return -myi / w;
  if (output_type == 2)
    return myi * w;
  return 1.0;
}

/// Pre-computed spectral constants derived from the damping parameter
/// and the Earth model reference period.
struct SpecConstants {
  std::complex<double> myi{0.0, 1.0};
  std::complex<double> ieps;
  double twopid, w0,
      twodivpi;   // twopid declared before w0 so init order is correct

  /**
   * @brief Pre-computes frequency-domain constants from damping and the
   * reference period.
   *
   * @param ep Exponential damping parameter.
   * @param tref Reference period from the Earth model.
   */
  SpecConstants(double ep, double tref)
      : ieps{-ep * myi}, twopid{2.0 * M_PI}, w0{twopid / tref},
        twodivpi{2.0 / M_PI} {}
};

/// Calls solver.compute() every nskip steps (full symbolic + numeric
/// factorization) and solver.factorize() otherwise (numeric only).
/// This amortises the cost of symbolic analysis over many solves.
template <class Solver, class Matrix>
void
factorizeOrCompute(Solver &solver, const Matrix &mat, int idxn, int nskip) {
  if ((idxn % nskip) == 0) {
    solver.compute(mat);
  } else {
    solver.factorize(mat);
  }
}

}   // namespace SPARSESPEC

#endif   // SPEC_HELPERS_GUARD_H
