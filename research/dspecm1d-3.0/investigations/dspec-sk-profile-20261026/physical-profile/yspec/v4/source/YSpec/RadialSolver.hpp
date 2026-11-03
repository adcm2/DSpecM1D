#pragma once

/**
 * @file RadialSolver.hpp
 * @brief The responses at a set of receiver radii to unit jumps at a set of
 * source radii, for one degree and one complex frequency.
 *
 * @details The sources lie in one solid region; the receivers may lie
 * anywhere. For each system a solution regular at the bottom (integrated
 * upwards from a starting level) and one satisfying the conditions at the top
 * (integrated downwards from the surface) are combined with the jump at each
 * source, in the manner of Friederich & Dalkolmo (1995) and Al-Attar &
 * Woodhouse (2008).
 *
 * The integration starts at the deepest point at which the solution can
 * matter: below the turning point of the slowest wave at the given degree and
 * frequency, beneath the deepest source, the solution decays, and the WKB
 * estimate of that decay is accumulated downwards until it reaches a given
 * number of e-folds. The response at a receiver below the start is taken to
 * be zero, to the same accuracy.
 *
 * Integrations restart at every discontinuity of the model, and the systems
 * change from solid to fluid form and back wherever the phase does, so a fluid
 * inner core, an ocean, or any other sequence of solid and fluid regions needs
 * no special treatment. The b-vectors that carry the source out to the
 * receivers change form with them, so a receiver may be in another region
 * from the source: in an ocean, or in the core. The toroidal response is zero
 * outside the solid region of the sources, and in a fluid the consoidal
 * displacement follows from the vanishing of the tangential traction.
 *
 * The solver works in the model's non-dimensional units. The jumps at the
 * source are unit jumps in the solution vector of the original YSpec, whose
 * displacements are scaled by the radius (@f$rU@f$, @f$\zeta rV@f$,
 * @f$\zeta rW@f$, with @f$\zeta = \sqrt{l(l+1)}@f$) and whose tractions are
 * the matching ones; the responses are the physical displacement
 * coefficients at the receiver.
 *
 * A solver holds no state that a solve changes, so one solver can be shared
 * by any number of threads.
 */

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

#include "Calculation.hpp"
#include "Model.hpp"
#include "Types.hpp"

namespace YSpec {

/** @brief A radius in a specified layer. */
struct Point {
  Int layer;  ///< The layer, which settles which side of a discontinuity.
  double r;   ///< The radius.
};

/**
 * @brief The toroidal response at the receiver for the two unit jumps: W, and
 * its radial derivative, which is zero in a fluid.
 */
struct ToroidalResponse {
  std::array<Complex, 2> W;   ///< The toroidal displacement coefficient.
  std::array<Complex, 2> DW;  ///< Its radial derivative.
};

/** @brief The radial (l = 0) response: U and its radial derivative. */
struct RadialResponse {
  std::array<Complex, 2> U;   ///< The radial displacement.
  std::array<Complex, 2> DU;  ///< Its radial derivative.
};

/**
 * @brief The spheroidal response at the receiver for the four unit jumps,
 * in the displacement and the traction.
 *
 * @details The potential and its gradient are zero unless the solver is
 * self-gravitating. The radial derivatives of the displacement follow from
 * the tractions; in a fluid, where there is no shear, that of V is not given
 * and is zero.
 */
struct SpheroidalResponse {
  std::array<Complex, 4> U;     ///< The radial displacement coefficient.
  std::array<Complex, 4> V;     ///< The consoidal displacement coefficient.
  std::array<Complex, 4> Phi;   ///< The perturbation of the potential.
  std::array<Complex, 4> DPhi;  ///< Its radial derivative.
  std::array<Complex, 4> DU;    ///< The radial derivative of U.
  std::array<Complex, 4> DV;    ///< The radial derivative of V.
};

/**
 * @brief A value for every pair of a source radius and a receiver radius.
 */
template <class T>
class ResponseTable {
 public:
  /** @brief A table for the given numbers of sources and receivers. */
  ResponseTable(std::size_t sources, std::size_t receivers)
      : receivers_{receivers}, values_(sources * receivers) {}

  /** @brief The number of source radii. */
  std::size_t Sources() const {
    return receivers_ == 0 ? 0 : values_.size() / receivers_;
  }
  /** @brief The number of receiver radii. */
  std::size_t Receivers() const { return receivers_; }

  /** @brief The value for a source and a receiver. */
  T& operator()(std::size_t source, std::size_t receiver) & {
    return values_[source * receivers_ + receiver];
  }
  /** @brief The value for a source and a receiver. */
  const T& operator()(std::size_t source, std::size_t receiver) const& {
    return values_[source * receivers_ + receiver];
  }
  /**
   * @brief The value for a source and a receiver, from a temporary table:
   * a copy, so that it cannot outlive the table it came from.
   */
  T operator()(std::size_t source, std::size_t receiver) && {
    return values_[source * receivers_ + receiver];
  }

 private:
  std::size_t receivers_;
  std::vector<T> values_;
};

/**
 * @brief Solves the radial equations for a set of source radii and a set of
 * receiver radii.
 *
 * @details The integrations upwards from the starting level and downwards
 * from the surface do not depend on the source, so one of each serves every
 * source and every receiver: the solutions are kept at each radius as the
 * integrations pass it. Each source then needs only its own determinant and
 * the integration of its b-vectors out to the receivers. Many receiver
 * radii, for a snapshot of the wavefield, or many source depths, cost little
 * more than one.
 */
class RadialSolver {
 public:
  /**
   * @param model The earth model, which must outlive the solver.
   * @param sourceRadii The radii of the sources, non-dimensional.
   * @param receiverRadii The radii of the receivers, non-dimensional.
   * @param gravity How gravity enters.
   * @param attenuation Whether to include attenuation.
   * @param options Integration settings.
   * @throws std::invalid_argument unless every source lies in one solid
   * region, if a source or receiver is at the centre, or if the options are
   * out of range. A source or receiver on a boundary between solid and fluid
   * is taken to be on the solid side.
   */
  RadialSolver(const Model& model, std::vector<double> sourceRadii,
               std::vector<double> receiverRadii, Gravity gravity,
               bool attenuation = true, SolverOptions options = {});

  /**
   * @brief A solver with its receivers given as points, each in its layer,
   * so that one at a discontinuity is on the side named.
   * @throws std::invalid_argument as the other constructor does, or if a
   * receiver is outside its layer.
   */
  RadialSolver(const Model& model, std::vector<double> sourceRadii,
               std::vector<Point> receivers, Gravity gravity,
               bool attenuation = true, SolverOptions options = {});

  /** @brief A solver for one source radius and one receiver radius. */
  RadialSolver(const Model& model, double sourceRadius, double receiverRadius,
               Gravity gravity, bool attenuation = true,
               SolverOptions options = {})
      : RadialSolver(model, std::vector{sourceRadius},
                     std::vector{receiverRadius}, gravity, attenuation,
                     options) {}

  /** @brief The model. */
  const Model& GetModel() const;
  /** @brief The number of source radii. */
  std::size_t Sources() const;
  /** @brief The number of receiver radii. */
  std::size_t Receivers() const;
  /** @brief A source radius, in its layer. */
  Point Source(std::size_t i = 0) const;
  /** @brief A receiver radius, in its layer. */
  Point Receiver(std::size_t i = 0) const;
  /** @brief How gravity enters. */
  Gravity GetGravity() const;
  /** @brief Whether attenuation is included. */
  bool Attenuation() const;

  /** @brief The toroidal responses for degree @p l >= 1. */
  ResponseTable<ToroidalResponse> Toroidal(int l, Complex omega) const;

  /** @brief The radial (l = 0) responses. */
  ResponseTable<RadialResponse> Radial(Complex omega) const;

  /** @brief The spheroidal responses for degree @p l >= 1. */
  ResponseTable<SpheroidalResponse> Spheroidal(int l, Complex omega) const;

  /**
   * @brief Where the upward integration starts for a degree and frequency.
   * @param l The degree.
   * @param omega The real part of the frequency.
   * @param top The point the integration goes to.
   * @param lowestRegion The deepest region to look in.
   * @return A knot below @p top, or the bottom of @p lowestRegion if the
   * solution does not decay enough above it.
   */
  Point StartPoint(int l, double omega, Point top, Int lowestRegion) const;

 private:
  class Impl;
  std::shared_ptr<const Impl> impl_;
};

}  // namespace YSpec
