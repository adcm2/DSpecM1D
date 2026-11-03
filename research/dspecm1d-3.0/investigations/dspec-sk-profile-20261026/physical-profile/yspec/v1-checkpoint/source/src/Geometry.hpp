#pragma once

/**
 * @file Geometry.hpp
 * @brief Where a receiver is seen from a source.
 *
 * @details In a frame with the source at the north pole a receiver is at
 * colatitude @f$\Delta@f$, the epicentral distance, and longitude
 * @f$\phi = \pi - \zeta_s@f$ with @f$\zeta_s@f$ the azimuth of the receiver
 * from the source; the horizontal components there are rotated to north and
 * east with the azimuth of the source from the receiver. Geographic latitudes
 * are converted to geocentric with a flattening of about 1/297, and a point
 * exactly at a pole is moved 0.02 degrees off it, where the azimuths are
 * undefined.
 */

#include <YSpec/Location.hpp>
#include <array>
#include <cmath>
#include <numbers>

#include "Sphere.hpp"

namespace YSpec {

/**
 * @brief A receiver as seen from the source.
 */
struct ReceiverGeometry {
  double distance;  ///< The epicentral distance (radians).
  double phi;       ///< The longitude in the frame with the source at the pole.
  double rotation;  ///< The angle that rotates the horizontal components.
};

/**
 * @brief The geometry of a point relative to a source, both given by their
 * geocentric colatitudes and longitudes (radians).
 */
inline ReceiverGeometry GeometryAt(double thetaR, double phiR, double thetaS,
                                   double phiS) {
  constexpr auto pi = std::numbers::pi;

  // The distance from the cross and dot products of the unit vectors, which
  // keeps its precision near 0 and pi where an arc cosine does not.
  const auto xR = Sphere::Unit(thetaR, phiR);
  const auto xS = Sphere::Unit(thetaS, phiS);
  const auto cross =
      std::hypot(xR[1] * xS[2] - xR[2] * xS[1], xR[2] * xS[0] - xR[0] * xS[2],
                 xR[0] * xS[1] - xR[1] * xS[0]);
  const auto dot = xR[0] * xS[0] + xR[1] * xS[1] + xR[2] * xS[2];
  const auto delta = std::atan2(cross, dot);

  // Where the two points coincide or are antipodal the azimuths are
  // undefined. Those chosen are the limits seen from a receiver approaching
  // along any great circle: at the source the receiver's azimuth plus pi is
  // the source's azimuth seen from the receiver, and at the antipode the two
  // azimuths are equal.
  if (cross == 0)
    return dot > 0 ? ReceiverGeometry{0, pi, pi}
                   : ReceiverGeometry{pi, pi, 2 * pi};

  // The azimuth at point a of point b, clockwise from north, in [0, 2 pi).
  auto azimuth = [](double thetaA, double phiA, double thetaB, double phiB) {
    const auto dPhi = phiB - phiA;
    const auto angle =
        std::atan2(std::sin(dPhi) * std::sin(thetaB),
                   std::sin(thetaA) * std::cos(thetaB) -
                       std::cos(thetaA) * std::sin(thetaB) * std::cos(dPhi));
    return angle < 0 ? angle + 2 * pi : angle;
  };
  const auto atReceiver = azimuth(thetaR, phiR, thetaS, phiS);
  const auto atSource = azimuth(thetaS, phiS, thetaR, phiR);
  return {delta, pi - atSource, 2 * pi - atReceiver};
}

/**
 * @brief The geometry of a receiver relative to a source.
 */
inline ReceiverGeometry Geometry(const Location& receiver,
                                 const Location& source) {
  using Sphere::Degree;
  auto offPole = [](double latitude) {
    if (latitude == 90) return 89.98;
    if (latitude == -90) return -89.98;
    return latitude;
  };
  return GeometryAt(Sphere::Colatitude(offPole(receiver.latitude)),
                    receiver.longitude * Degree,
                    Sphere::Colatitude(offPole(source.latitude)),
                    source.longitude * Degree);
}

}  // namespace YSpec
