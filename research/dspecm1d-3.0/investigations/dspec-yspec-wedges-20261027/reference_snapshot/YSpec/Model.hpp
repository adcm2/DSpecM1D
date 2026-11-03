#pragma once

/**
 * @file Model.hpp
 * @brief A spherically symmetric, transversely isotropic, anelastic earth
 * model, read from a deck and made non-dimensional.
 *
 * @details The model is a stack of layers separated by discontinuities. Within
 * a layer each property is a natural cubic spline through the deck's knots,
 * which is what Minos and the original YSpec use. All the properties of a
 * layer share its knots, so a single search finds the segment and every
 * property is then a cubic in the distance from its lower knot.
 *
 * Gravity is computed from the density splines exactly: the mass below a
 * radius is the integral of a cubic times @f$r^2@f$, a polynomial integrated
 * in closed form.
 *
 * Consecutive layers of the same phase form a region. The equations change
 * form only between regions, where solid meets fluid.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <ranges>
#include <span>
#include <string_view>
#include <vector>

#include "Types.hpp"
#include "Units.hpp"

namespace YSpec {

/** @brief Whether a layer can support shear. */
enum class Phase { Solid, Fluid };

/** @brief Which of two layers meeting at a discontinuity is meant. */
enum class Side { Below, Above };

/**
 * @brief The properties of the model at one radius, non-dimensional.
 *
 * @details The anelastic quantities are inverse quality factors, zero where
 * the deck gives no attenuation. Following Minos, @f$A@f$ and @f$C@f$ are
 * attenuated with the effective inverse quality factor @c qA, @f$F@f$ with
 * @c qF, and @f$L@f$ and @f$N@f$ with @c qMu; for an isotropic material this
 * is the same as attenuating the bulk and shear moduli.
 */
struct Material {
  double rho;     ///< Density.
  double A;       ///< Love modulus A.
  double C;       ///< Love modulus C.
  double F;       ///< Love modulus F.
  double L;       ///< Love modulus L.
  double N;       ///< Love modulus N.
  double kappa;   ///< Voigt average bulk modulus.
  double mu;      ///< Voigt average shear modulus.
  double qKappa;  ///< Inverse quality factor in bulk.
  double qMu;     ///< Inverse quality factor in shear.
  double qA;      ///< Effective inverse quality factor of A and C.
  double qF;      ///< Effective inverse quality factor of F.
  double g;       ///< Gravitational acceleration.

  /** @brief The isotropic P-wave speed. */
  auto Alpha() const { return std::sqrt((kappa + 4 * mu / 3) / rho); }
  /** @brief The isotropic S-wave speed. */
  auto Beta() const { return std::sqrt(mu / rho); }
};

/** @brief The properties a layer splines: those of Material less gravity. */
enum class Property : std::size_t {
  Rho,
  A,
  C,
  F,
  L,
  N,
  Kappa,
  Mu,
  QKappa,
  QMu,
  QA,
  QF,
  Count
};

/** @brief A value for each splined property, indexed by Property. */
class Properties {
 public:
  /** @brief The value of a property. */
  double& operator[](Property p) {
    return values_[static_cast<std::size_t>(p)];
  }
  /** @brief The value of a property. */
  double operator[](Property p) const {
    return values_[static_cast<std::size_t>(p)];
  }
  /** @brief The values, in the order of Property. */
  auto& Values() { return values_; }
  /** @brief The values, in the order of Property. */
  const auto& Values() const { return values_; }

 private:
  std::array<double, static_cast<std::size_t>(Property::Count)> values_{};
};

/**
 * @brief One layer of the model: the region between two discontinuities.
 */
class Layer {
 public:
  /**
   * @brief Builds a layer from its knots and the properties at each.
   * @param radii The knots, strictly increasing, at least two.
   * @param values The properties at each knot.
   * @param phase Whether the layer is solid or fluid.
   * @param massBelow The integral of @f$\rho r^2@f$ from the centre to the
   * layer's lower radius.
   * @throws std::invalid_argument for fewer than two knots or a value
   * missing.
   */
  Layer(std::vector<double> radii, const std::vector<Properties>& values,
        Phase phase, double massBelow);

  /** @brief The lower radius. */
  auto Lower() const { return r_.front(); }
  /** @brief The upper radius. */
  auto Upper() const { return r_.back(); }
  /** @brief Whether the layer is solid or fluid. */
  auto GetPhase() const { return phase_; }
  /** @brief Whether the layer is fluid. */
  auto IsFluid() const { return phase_ == Phase::Fluid; }
  /** @brief The knots. */
  std::span<const double> Knots() const { return r_; }
  /** @brief The integral of @f$\rho r^2@f$ from the centre to the top. */
  auto MassBelowTop() const { return massBelow_.back(); }

  /**
   * @brief The properties at a radius.
   * @details A radius outside the layer is given the value of the end cubic
   * continued, which is what an integration step overshooting the layer by
   * rounding needs. Evaluated for every step of every integration, so it is
   * defined here, to be inlined.
   */
  Material operator()(double r) const {
    const auto k = Segment(r);
    const auto t = r - r_[k];
    const auto& c = coefficients_[k];
    auto p = Properties{};
    for (auto i : std::views::iota(std::size_t{0}, p.Values().size())) {
      p.Values()[i] =
          ((c[3].Values()[i] * t + c[2].Values()[i]) * t + c[1].Values()[i]) *
              t +
          c[0].Values()[i];
    }
    return {.rho = p[Property::Rho],
            .A = p[Property::A],
            .C = p[Property::C],
            .F = p[Property::F],
            .L = p[Property::L],
            .N = p[Property::N],
            .kappa = p[Property::Kappa],
            .mu = p[Property::Mu],
            .qKappa = p[Property::QKappa],
            .qMu = p[Property::QMu],
            .qA = p[Property::QA],
            .qF = p[Property::QF],
            .g = Gravity(k, t, r)};
  }

  /** @brief The gravitational acceleration at a radius. */
  double Gravity(double r) const {
    const auto k = Segment(r);
    return Gravity(k, r - r_[k], r);
  }

 private:
  // Per segment, the coefficients of the cubic in the distance from the
  // segment's lower knot, constant term first, for each property.
  using Cubic = std::array<Properties, 4>;

  std::vector<double> r_;
  Phase phase_;
  std::vector<Cubic> coefficients_;
  std::vector<double> massBelow_;

  std::size_t Segment(double r) const {
    const auto it = std::ranges::upper_bound(r_, r);
    const auto k = std::distance(r_.begin(), it) - 1;
    return static_cast<std::size_t>(
        std::clamp<std::ptrdiff_t>(k, 0, std::ssize(r_) - 2));
  }

  // The integral of rho r^2 over [r_k, r_k + t], for the density cubic
  // rho = a0 + a1 t + a2 t^2 + a3 t^3 and r^2 = r_k^2 + 2 r_k t + t^2.
  double MassInSegment(std::size_t k, double t) const {
    const auto& c = coefficients_[k];
    const auto r0 = r_[k];
    auto sum = 0.0;
    auto tPower = t;  // t^(i+1)
    for (auto i : std::views::iota(0, 4)) {
      const auto a = c[i][Property::Rho];
      sum += a * (r0 * r0 * tPower / (i + 1) + 2 * r0 * tPower * t / (i + 2) +
                  tPower * t * t / (i + 3));
      tPower *= t;
    }
    return sum;
  }

  double Gravity(std::size_t k, double t, double r) const {
    if (r <= 0) return 0;
    return 4 * PiG * (massBelow_[k] + MassInSegment(k, t)) / (r * r);
  }
};

/**
 * @brief A run of consecutive layers of the same phase.
 */
struct Region {
  Int first;     ///< The index of the lowest layer.
  Int last;      ///< The index of the highest layer.
  Phase phase;   ///< The phase of every layer in it.
  double lower;  ///< The lower radius.
  double upper;  ///< The upper radius.

  /** @brief Whether the region is fluid. */
  auto IsFluid() const { return phase == Phase::Fluid; }
  /** @brief The layer indices, lowest first. */
  auto Layers() const { return std::views::iota(first, last + 1); }
};

/**
 * @brief The earth model.
 */
class Model {
 public:
  /**
   * @brief Reads a model in the deck format of Minos (six or nine columns,
   * SI units).
   * @param path The deck.
   * @param scales The scales with which to make it non-dimensional.
   * @throws std::runtime_error if the deck cannot be read, and
   * std::invalid_argument if it does not describe a usable model.
   */
  static Model ReadDeck(const std::filesystem::path& path,
                        Scales scales = Scales{});

  /**
   * @brief Reads a model from the text of a deck.
   * @param text The deck.
   * @param scales The scales with which to make it non-dimensional.
   * @throws as ReadDeck.
   */
  static Model ReadDeckText(std::string_view text, Scales scales = Scales{});

  /** @brief The scales with which the model was made non-dimensional. */
  const Scales& GetScales() const { return scales_; }

  /** @brief The frequency at which the moduli are given, non-dimensional. */
  auto ReferenceFrequency() const { return referenceFrequency_; }

  /** @brief The number of layers. */
  auto NumberOfLayers() const { return std::ssize(layers_); }
  /** @brief A layer, numbered from the centre. */
  const Layer& GetLayer(Int i) const { return layers_[i]; }
  /** @brief The layers, from the centre. */
  std::span<const Layer> Layers() const { return layers_; }

  /** @brief The regions of one phase, from the centre. */
  std::span<const Region> Regions() const { return regions_; }

  /** @brief The radius of the centre (zero unless the model is hollow). */
  auto InnerRadius() const { return layers_.front().Lower(); }
  /** @brief The radius of the surface. */
  auto SurfaceRadius() const { return layers_.back().Upper(); }

  /**
   * @brief The layer containing a radius.
   * @param r The radius.
   * @param side At a discontinuity, which of the two layers to return.
   * @throws std::out_of_range if the radius is outside the model.
   */
  Int LayerIndex(double r, Side side = Side::Below) const;

  /** @brief The index of the region containing a layer. */
  Int RegionIndex(Int layer) const;

  /** @brief The properties at a radius, in the layer given. */
  Material operator()(Int layer, double r) const { return layers_[layer](r); }

  /**
   * @brief Builds a model from its layers.
   * @param scales The scales with which it was made non-dimensional.
   * @param referenceFrequency The frequency at which the moduli are given.
   * @param layers The layers, from the centre, each beginning where the
   * last ends.
   * @throws std::invalid_argument if there are no layers.
   */
  Model(Scales scales, double referenceFrequency, std::vector<Layer> layers);

 private:
  Scales scales_;
  double referenceFrequency_;
  std::vector<Layer> layers_;
  std::vector<Region> regions_;
};

}  // namespace YSpec
