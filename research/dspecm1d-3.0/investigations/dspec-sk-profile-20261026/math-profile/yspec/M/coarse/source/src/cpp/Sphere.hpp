#pragma once

// Geographic and geocentric latitudes. The model is spherical, so distances
// and great circles are those of geocentric latitude; locations are given in
// geographic latitude, converted with a flattening of about 1/297.

#include <array>
#include <cmath>
#include <numbers>

namespace YSpec::Sphere {

inline constexpr auto Degree = std::numbers::pi / 180;

// The ratio of the tangents of geocentric and geographic latitude.
inline constexpr auto GeocentricRatio = 0.993277;

// The geocentric colatitude (radians) of a geographic latitude (degrees).
inline double Colatitude(double latitude) {
  return std::numbers::pi / 2 -
         std::atan(GeocentricRatio * std::tan(latitude * Degree));
}

// The geographic latitude (degrees) of a geocentric colatitude (radians).
inline double Latitude(double colatitude) {
  return std::atan(std::tan(std::numbers::pi / 2 - colatitude) /
                   GeocentricRatio) /
         Degree;
}

using Vector = std::array<double, 3>;

// The unit vector at a colatitude and longitude (radians).
inline Vector Unit(double theta, double phi) {
  return {std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi),
          std::cos(theta)};
}

}  // namespace YSpec::Sphere
