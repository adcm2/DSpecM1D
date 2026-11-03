#include <YSpec/Location.hpp>
#include <cmath>
#include <numbers>
#include <ranges>
#include <stdexcept>

#include "Sphere.hpp"

namespace YSpec {

namespace {

using Sphere::Degree;
using Sphere::Vector;

Vector UnitVector(const Location& location) {
  return Sphere::Unit(Sphere::Colatitude(location.latitude),
                      location.longitude * Degree);
}

Location ToLocation(const Vector& x) {
  const auto theta = std::atan2(std::hypot(x[0], x[1]), x[2]);
  return {Sphere::Latitude(theta), std::atan2(x[1], x[0]) / Degree};
}

Vector Combine(double a, const Vector& x, double b, const Vector& y) {
  return {a * x[0] + b * y[0], a * x[1] + b * y[1], a * x[2] + b * y[2]};
}

// The fraction of the way along of point i of count.
double Fraction(std::size_t i, std::size_t count) {
  return count == 1 ? 0.0
                    : static_cast<double>(i) / static_cast<double>(count - 1);
}

}  // namespace

std::vector<Location> GreatCircle(const Location& from, const Location& to,
                                  std::size_t count) {
  if (count == 0) throw std::invalid_argument("a profile needs a point");
  const auto a = UnitVector(from);
  const auto b = UnitVector(to);
  const auto cross =
      std::hypot(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                 a[0] * b[1] - a[1] * b[0]);
  const auto dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const auto delta = std::atan2(cross, dot);
  if (std::numbers::pi - delta < 1e-9) {
    throw std::invalid_argument(
        "the ends of a profile are antipodal, so the great circle between "
        "them is undefined: give an azimuth instead");
  }
  auto points = std::vector<Location>{};
  for (auto i : std::views::iota(std::size_t{0}, count)) {
    if (delta == 0) {
      points.push_back(from);
      continue;
    }
    // Spherical linear interpolation.
    const auto t = Fraction(i, count);
    points.push_back(ToLocation(
        Combine(std::sin((1 - t) * delta), a, std::sin(t * delta), b)));
  }
  // The ends exactly as given.
  points.front() = from;
  if (count > 1) points.back() = to;
  return points;
}

std::vector<Location> GreatCircle(const Location& from, double azimuth,
                                  double minDistance, double maxDistance,
                                  std::size_t count) {
  if (count == 0) throw std::invalid_argument("a profile needs a point");
  const auto theta = Sphere::Colatitude(from.latitude);
  const auto phi = from.longitude * Degree;
  const auto x = Sphere::Unit(theta, phi);
  // North and east at the starting point.
  const auto north = Vector{-std::cos(theta) * std::cos(phi),
                            -std::cos(theta) * std::sin(phi), std::sin(theta)};
  const auto east = Vector{-std::sin(phi), std::cos(phi), 0};
  const auto heading = Combine(std::cos(azimuth * Degree), north,
                               std::sin(azimuth * Degree), east);
  auto points = std::vector<Location>{};
  for (auto i : std::views::iota(std::size_t{0}, count)) {
    const auto t = Fraction(i, count);
    const auto delta = ((1 - t) * minDistance + t * maxDistance) * Degree;
    if (delta == 0) {
      points.push_back(from);
      continue;
    }
    points.push_back(
        ToLocation(Combine(std::cos(delta), x, std::sin(delta), heading)));
  }
  return points;
}

}  // namespace YSpec
