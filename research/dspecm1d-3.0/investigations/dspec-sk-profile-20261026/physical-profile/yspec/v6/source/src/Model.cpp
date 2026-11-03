#include <Interpolation/CubicSplineSystem.hpp>
#include <YSpec/Model.hpp>
#include <format>
#include <numbers>
#include <planetmodel/deck.hpp>
#include <stdexcept>
#include <utility>

namespace YSpec {

Layer::Layer(std::vector<double> radii, const std::vector<Properties>& values,
             Phase phase, double massBelow)
    : r_{std::move(radii)}, phase_{phase} {
  if (r_.size() < 2 || values.size() != r_.size()) {
    throw std::invalid_argument(
        "Layer: need at least two knots, with a value at each");
  }
  const auto knots = r_.size();
  const auto system = Interpolation::CubicSplineSystem{std::views::all(r_)};
  auto y = std::vector<double>(knots);
  auto curvature = std::vector<double>(knots);
  coefficients_.assign(knots - 1, {});
  for (auto property :
       std::views::iota(std::size_t{0}, Properties{}.Values().size())) {
    std::ranges::transform(values, y.begin(), [property](const auto& v) {
      return v.Values()[property];
    });
    system.Solve(y, std::span{curvature});
    for (auto k : std::views::iota(std::size_t{0}, knots - 1)) {
      const auto width = r_[k + 1] - r_[k];
      auto& segment = coefficients_[k];
      segment[0].Values()[property] = y[k];
      segment[1].Values()[property] =
          (y[k + 1] - y[k]) / width -
          width * (2 * curvature[k] + curvature[k + 1]) / 6;
      segment[2].Values()[property] = curvature[k] / 2;
      segment[3].Values()[property] =
          (curvature[k + 1] - curvature[k]) / (6 * width);
    }
  }
  // The mass below each knot, from the exact integral of the density.
  massBelow_.reserve(knots);
  massBelow_.push_back(massBelow);
  for (auto k : std::views::iota(std::size_t{0}, knots - 1)) {
    massBelow_.push_back(massBelow_.back() +
                         MassInSegment(k, r_[k + 1] - r_[k]));
  }
}

namespace {

// The layers of a deck, made non-dimensional, with their Love moduli, Voigt
// averages and Minos's effective quality factors.
std::vector<Layer> Layers(const planetmodel::Deck& deck, const Scales& scales) {
  const auto anisotropic = deck.HasColumn("vpv");
  const auto& radius = deck.Radius();
  const auto& rho = deck.Column("rho");
  const auto& vpv = deck.Column(anisotropic ? "vpv" : "vp");
  const auto& vsv = deck.Column(anisotropic ? "vsv" : "vs");
  const auto& vph = deck.Column(anisotropic ? "vph" : "vp");
  const auto& vsh = deck.Column(anisotropic ? "vsh" : "vs");
  const auto& qKappa = deck.Column("qkappa");
  const auto& qMu = deck.Column("qmu");

  const auto velocity = scales.Velocity();
  auto inverse = [](double q) { return q > 0 ? 1 / q : 0.0; };
  auto layers = std::vector<Layer>{};
  auto massBelow = 0.0;
  for (auto [begin, end] : deck.Layers()) {
    auto r = std::vector<double>{};
    auto values = std::vector<Properties>{};
    auto shear = 0.0;
    for (auto i : std::views::iota(begin, end)) {
      const auto density = rho[i] / scales.Density();
      auto modulus = [&](double v) {
        const auto u = v / velocity;
        return density * u * u;
      };
      auto p = Properties{};
      p[Property::Rho] = density;
      p[Property::A] = modulus(vph[i]);
      p[Property::C] = modulus(vpv[i]);
      p[Property::L] = modulus(vsv[i]);
      p[Property::N] = modulus(vsh[i]);
      const auto eta = anisotropic ? deck.Column("eta")[i] : 1.0;
      p[Property::F] = eta * (p[Property::A] - 2 * p[Property::L]);
      const auto kappa = (p[Property::C] + 4 * p[Property::A] -
                          4 * p[Property::N] + 4 * p[Property::F]) /
                         9;
      const auto mu = (p[Property::C] + p[Property::A] + 6 * p[Property::L] +
                       5 * p[Property::N] - 2 * p[Property::F]) /
                      15;
      p[Property::Kappa] = kappa;
      p[Property::Mu] = mu;
      const auto inverseQKappa = inverse(qKappa[i]);
      const auto inverseQMu = inverse(qMu[i]);
      p[Property::QKappa] = inverseQKappa;
      p[Property::QMu] = inverseQMu;
      // Minos's effective inverse quality factors for A, C and F.
      const auto ratio = 4 * mu / (3 * kappa + 4 * mu);
      p[Property::QA] = (1 - ratio) * inverseQKappa + ratio * inverseQMu;
      p[Property::QF] = ((1 - ratio) * inverseQKappa - ratio * inverseQMu / 2) /
                        (1 - 1.5 * ratio);
      r.push_back(radius[i] / scales.Length());
      values.push_back(p);
      shear = std::max({shear, p[Property::L], p[Property::N]});
    }
    const auto phase = shear > 0 ? Phase::Solid : Phase::Fluid;
    layers.emplace_back(std::move(r), values, phase, massBelow);
    massBelow = layers.back().MassBelowTop();
  }
  return layers;
}

Model FromDeck(const planetmodel::Deck& deck, Scales scales) {
  const auto header =
      planetmodel::ReadMineosHeader(deck.Header(), deck.Source());
  const auto referenceFrequency =
      2 * std::numbers::pi / header.tref / scales.Frequency();
  return Model(scales, referenceFrequency, Layers(deck, scales));
}

}  // namespace

Model Model::ReadDeck(const std::filesystem::path& path, Scales scales) {
  return FromDeck(planetmodel::ReadDeck(path, planetmodel::MINEOS()), scales);
}

Model Model::ReadDeckText(std::string_view text, Scales scales) {
  return FromDeck(planetmodel::ReadDeckText(text, planetmodel::MINEOS()),
                  scales);
}

Model::Model(Scales scales, double referenceFrequency,
             std::vector<Layer> layers)
    : scales_{scales},
      referenceFrequency_{referenceFrequency},
      layers_{std::move(layers)} {
  if (layers_.empty()) throw std::invalid_argument("Model: no layers");
  // Group the layers into regions of one phase.
  for (auto i : std::views::iota(Int{0}, std::ssize(layers_))) {
    const auto& layer = layers_[i];
    if (regions_.empty() || regions_.back().phase != layer.GetPhase()) {
      regions_.push_back({.first = i,
                          .last = i,
                          .phase = layer.GetPhase(),
                          .lower = layer.Lower(),
                          .upper = layer.Upper()});
    } else {
      regions_.back().last = i;
      regions_.back().upper = layer.Upper();
    }
  }
}

Int Model::LayerIndex(double r, Side side) const {
  if (r < InnerRadius() || r > SurfaceRadius()) {
    throw std::out_of_range(std::format("radius {} is outside the model", r));
  }
  auto it = side == Side::Below
                ? std::ranges::lower_bound(layers_, r, {}, &Layer::Upper)
                : std::ranges::upper_bound(layers_, r, {}, &Layer::Upper);
  if (it == layers_.end()) --it;
  return std::distance(layers_.begin(), it);
}

Int Model::RegionIndex(Int layer) const {
  const auto it = std::ranges::find_if(
      regions_, [layer](const auto& region) { return layer <= region.last; });
  return std::distance(regions_.begin(), it);
}

}  // namespace YSpec
