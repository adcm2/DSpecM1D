#pragma once

/**
 * @file Location.hpp
 * @brief A point on the surface of the earth, and lines of points along great
 * circles.
 *
 * @details Locations are in geographic latitude. Since the model is
 * spherical, the great circles here, and their distances and azimuths, are
 * those of geocentric latitude, which is what the epicentral distances of the
 * seismograms are measured in.
 */

#include <cstddef>
#include <vector>

namespace YSpec {

/** @brief A point on the surface, in geographic degrees. */
struct Location {
  double latitude = 0;   ///< Geographic latitude, degrees north.
  double longitude = 0;  ///< Longitude, degrees east.
};

/**
 * @brief Points equally spaced along the shorter great circle between two
 * points, both included.
 * @param from The first point.
 * @param to The last point.
 * @param count The number of points; one gives @p from alone.
 * @throws std::invalid_argument if the count is zero, or if the points are
 * antipodal, so that the great circle is undefined.
 */
std::vector<Location> GreatCircle(const Location& from, const Location& to,
                                  std::size_t count);

/**
 * @brief Points along the great circle leaving a point at an azimuth, at
 * distances equally spaced between two limits, both included.
 * @param from The point the distances are measured from.
 * @param azimuth The azimuth at @p from, degrees clockwise from north.
 * @param minDistance The distance of the first point (degrees).
 * @param maxDistance The distance of the last point (degrees).
 * @param count The number of points; one gives the point at @p minDistance.
 * @throws std::invalid_argument if the count is zero.
 */
std::vector<Location> GreatCircle(const Location& from, double azimuth,
                                  double minDistance, double maxDistance,
                                  std::size_t count);

}  // namespace YSpec
