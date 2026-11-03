#pragma once

/**
 * @file Source.hpp
 * @brief A point source: a moment tensor, a force, or both.
 */

#include "Location.hpp"

namespace YSpec {

/**
 * @brief The components of a moment tensor (N m) in spherical coordinates
 * @f$(r, \theta, \phi)@f$ at the source: r up, @f$\theta@f$ south and
 * @f$\phi@f$ east.
 */
struct MomentTensor {
  double rr = 0;  ///< @f$M_{rr}@f$.
  double rt = 0;  ///< @f$M_{r\theta}@f$.
  double rp = 0;  ///< @f$M_{r\phi}@f$.
  double tt = 0;  ///< @f$M_{\theta\theta}@f$.
  double tp = 0;  ///< @f$M_{\theta\phi}@f$.
  double pp = 0;  ///< @f$M_{\phi\phi}@f$.

  /** @brief The tensor divided by a scale. */
  MomentTensor operator/(double scale) const {
    return {rr / scale, rt / scale, rp / scale,
            tt / scale, tp / scale, pp / scale};
  }
};

/**
 * @brief The components of a force (N) in spherical coordinates at the
 * source: r up, @f$\theta@f$ south and @f$\phi@f$ east.
 */
struct Force {
  double r = 0;  ///< @f$f_r@f$.
  double t = 0;  ///< @f$f_\theta@f$.
  double p = 0;  ///< @f$f_\phi@f$.

  /** @brief The force divided by a scale. */
  Force operator/(double scale) const {
    return {r / scale, t / scale, p / scale};
  }
};

/**
 * @brief A point source: a step at the origin time in a moment tensor, a
 * force, or both.
 */
struct Source {
  Location location;    ///< The epicentre.
  double depth = 0;     ///< The depth (m).
  MomentTensor moment;  ///< The moment tensor (N m).
  Force force = {};     ///< The force (N).
};

}  // namespace YSpec
