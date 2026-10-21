#ifndef DSPECM1D_MODEL_INFO_MODEL_INPUT_H
#define DSPECM1D_MODEL_INFO_MODEL_INPUT_H

#include <cmath>
#include <format>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <planetmodel/deck.hpp>

namespace EarthModels {

template <typename FLOAT = double>
class EarthConstants {
public:
  using value_type = FLOAT;
  FLOAT LengthNorm() const { return m_lengthNorm; }
  FLOAT MassNorm() const { return m_massNorm; }
  FLOAT TimeNorm() const { return m_timeNorm; }
  FLOAT DensityNorm() const { return m_massNorm / std::pow(m_lengthNorm, 3.0); }
  FLOAT InertiaNorm() const { return m_massNorm * std::pow(m_lengthNorm, 2.0); }
  FLOAT VelocityNorm() const { return m_lengthNorm / m_timeNorm; }
  FLOAT AccelerationNorm() const { return m_lengthNorm / std::pow(m_timeNorm, 2.0); }
  FLOAT ForceNorm() const { return m_massNorm * m_lengthNorm / std::pow(m_timeNorm, 2.0); }
  FLOAT StressNorm() const { return m_massNorm / (m_lengthNorm * std::pow(m_timeNorm, 2.0)); }
  FLOAT GravitationalConstant() const {
    return std::pow(m_lengthNorm, 3.0) / (m_massNorm * std::pow(m_timeNorm, 2.0));
  }
private:
  const FLOAT m_lengthNorm = 6.371e6;
  const FLOAT m_massNorm = 5.972e24;
  const FLOAT m_timeNorm = 3600.0;
};

/** Owned, unit-aware view of one planetmodel natural-spline deck field. */
class ModelField {
public:
  ModelField() = default;
  ModelField(planetmodel::RadialFunction field, double lowerRadius, double upperRadius)
      : field_(std::make_shared<planetmodel::RadialFunction>(std::move(field))),
        lowerRadius_(lowerRadius), upperRadius_(upperRadius) {}

  double operator()(double radius) const { return Evaluate(radius); }
  double Evaluate(double radius) const { return evaluate(radius, 0); }
  template <int Derivative>
  double Evaluate(double radius) const {
    static_assert(Derivative == 0 || Derivative == 1);
    return evaluate(radius, Derivative);
  }

private:
  double evaluate(double radius, int derivative) const {
    if (!(radius >= lowerRadius_ && radius <= upperRadius_)) {
      throw std::out_of_range(std::format(
          "Normalized radius {} is outside model layer [{}, {}]", radius, lowerRadius_, upperRadius_));
    }
    return derivative == 0 ? field_->Value(radius) : field_->Derivative(radius);
  }

  std::shared_ptr<const planetmodel::RadialFunction> field_;
  double lowerRadius_ = 0.0, upperRadius_ = 0.0;
};

template <typename FLOAT = double, typename INTEGRAL = int>
class ModelInput {
public:
  using size_type = INTEGRAL;
  using value_type = FLOAT;
  using myvector = std::vector<FLOAT>;
  // Compatibility name: interpolation is now owned by planetmodel::RadialFunction.
  using InterpA = ModelField;

  ModelInput() = default;
  explicit ModelInput(const std::string &pathToFile);
  template <template <typename> class ParameterModel>
  ModelInput(const std::string &pathToFile, const ParameterModel<FLOAT> &modelConstants);

  FLOAT LengthNorm() const { return m_lengthNorm; }
  FLOAT MassNorm() const { return m_massNorm; }
  FLOAT TimeNorm() const { return m_timeNorm; }
  FLOAT DensityNorm() const { return m_densityNorm; }
  FLOAT InertiaNorm() const { return m_inertiaNorm; }
  FLOAT VelocityNorm() const { return m_velocityNorm; }
  FLOAT AccelerationNorm() const { return m_accelerationNorm; }
  FLOAT ForceNorm() const { return m_forceNorm; }
  FLOAT StressNorm() const { return m_stressNorm; }
  FLOAT GravitationalConstant() const { return m_gravitationalConstant; }
  FLOAT TREF() const { return m_tref; }

  int NumberOfLayers() const { return m_numLayers; }
  int LayerLowerIndex(int i) const { return m_layerIndices.at(i).first; }
  int LayerUpperIndex(int i) const { return m_layerIndices.at(i).second; }
  int LayerIndexDifference(int i) const { return LayerUpperIndex(i) - LayerLowerIndex(i); }
  auto LayerRadii() const { return m_layeredRadii; }
  auto LayerRadii(int i) const { return m_layeredRadii.at(i); }
  auto LayerRadiiNumber(int i) const { return m_layeredRadii.at(i).size(); }

  FLOAT LowerRadius(INTEGRAL i) const { return m_layerBounds.at(validLayer(i)); }
  FLOAT UpperRadius(INTEGRAL i) const { return m_layerBounds.at(validLayer(i) + 1); }
  FLOAT OuterRadius() const { return m_layerBounds.back(); }
  bool IsIsotropic() const { return m_isIsotropic; }
  bool IsSolid(INTEGRAL i) const { return m_isSolid.at(validLayer(i)); }
  bool IsFluid(INTEGRAL i) const { return !IsSolid(i); }

  InterpA Density(INTEGRAL i) const { return checkedLayer(m_density, i); }
  InterpA VPV(INTEGRAL i) const { return checkedLayer(m_vpv, i); }
  InterpA VPH(INTEGRAL i) const { return checkedLayer(m_vph, i); }
  InterpA VSV(INTEGRAL i) const { return checkedLayer(m_vsv, i); }
  InterpA VSH(INTEGRAL i) const { return checkedLayer(m_vsh, i); }
  InterpA QKappa(INTEGRAL i) const { return checkedLayer(m_qkappa, i); }
  InterpA QMu(INTEGRAL i) const { return checkedLayer(m_qmu, i); }
  InterpA Eta(INTEGRAL i) const { return checkedLayer(m_eta, i); }

  auto VS(INTEGRAL i) const {
    auto mu = Mu(i);
    auto rho = Density(i);
    return [mu = std::move(mu), rho = std::move(rho)](FLOAT x) { return std::sqrt(mu(x) / rho(x)); };
  }
  auto VP(INTEGRAL i) const {
    auto kappa = Kappa(i);
    auto mu = Mu(i);
    auto rho = Density(i);
    return [kappa = std::move(kappa), mu = std::move(mu), rho = std::move(rho)](FLOAT x) {
      return std::sqrt((kappa(x) + 4.0 / 3.0 * mu(x)) / rho(x));
    };
  }
  auto A(INTEGRAL i) const {
    auto rho = Density(i);
    auto vph = VPH(i);
    return [rho = std::move(rho), vph = std::move(vph)](FLOAT x) { return rho(x) * vph(x) * vph(x); };
  }
  auto C(INTEGRAL i) const {
    auto rho = Density(i);
    auto vpv = VPV(i);
    return [rho = std::move(rho), vpv = std::move(vpv)](FLOAT x) { return rho(x) * vpv(x) * vpv(x); };
  }
  auto N(INTEGRAL i) const {
    auto rho = Density(i);
    auto vsh = VSH(i);
    return [rho = std::move(rho), vsh = std::move(vsh)](FLOAT x) { return rho(x) * vsh(x) * vsh(x); };
  }
  auto L(INTEGRAL i) const {
    auto rho = Density(i);
    auto vsv = VSV(i);
    return [rho = std::move(rho), vsv = std::move(vsv)](FLOAT x) { return rho(x) * vsv(x) * vsv(x); };
  }
  auto F(INTEGRAL i) const {
    auto eta = Eta(i);
    auto a = A(i);
    auto l = L(i);
    return [eta = std::move(eta), a = std::move(a), l = std::move(l)](FLOAT x) {
      return eta(x) * (a(x) - 2 * l(x));
    };
  }
  auto Kappa(INTEGRAL i) const {
    auto c = C(i);
    auto a = A(i);
    auto n = N(i);
    auto f = F(i);
    return [c = std::move(c), a = std::move(a), n = std::move(n), f = std::move(f)](FLOAT x) {
      return (c(x) + 4.0 * (a(x) - n(x) + f(x))) / 9.0;
    };
  }
  auto Mu(INTEGRAL i) const {
    auto c = C(i);
    auto a = A(i);
    auto l = L(i);
    auto n = N(i);
    auto f = F(i);
    return [c = std::move(c), a = std::move(a), l = std::move(l), n = std::move(n), f = std::move(f)](FLOAT x) {
      return (c(x) + a(x) + 6.0 * l(x) + 5.0 * n(x) - 2.0 * f(x)) / 15.0;
    };
  }

private:
  template <template <typename> class ParameterModel>
  void initializeNorms(const ParameterModel<FLOAT> &modelConstants) {
    m_lengthNorm = modelConstants.LengthNorm();
    m_massNorm = modelConstants.MassNorm();
    m_timeNorm = modelConstants.TimeNorm();
    m_densityNorm = m_massNorm / std::pow(m_lengthNorm, 3.0);
    m_velocityNorm = m_lengthNorm / m_timeNorm;
    m_accelerationNorm = m_lengthNorm / std::pow(m_timeNorm, 2.0);
    m_forceNorm = m_massNorm * m_lengthNorm / std::pow(m_timeNorm, 2.0);
    m_stressNorm = m_massNorm / (m_lengthNorm * std::pow(m_timeNorm, 2.0));
    m_inertiaNorm = m_massNorm * std::pow(m_lengthNorm, 2.0);
    m_gravitationalConstant = std::pow(m_lengthNorm, 3.0) / (m_massNorm * std::pow(m_timeNorm, 2.0));
  }
  int validLayer(INTEGRAL i) const {
    if (i < 0 || i >= m_numLayers) throw std::invalid_argument("Invalid layer index");
    return static_cast<int>(i);
  }
  template <typename T> T checkedLayer(const std::vector<T> &values, INTEGRAL i) const {
    return values.at(validLayer(i));
  }
  void loadModel(const std::string &pathToFile);

  std::string m_modelTitle;
  int m_numLayers = 0;
  FLOAT m_tref = 0.0;
  FLOAT m_lengthNorm = 0.0, m_massNorm = 0.0, m_timeNorm = 0.0;
  FLOAT m_densityNorm = 0.0, m_inertiaNorm = 0.0, m_velocityNorm = 0.0;
  FLOAT m_accelerationNorm = 0.0, m_forceNorm = 0.0, m_stressNorm = 0.0;
  FLOAT m_gravitationalConstant = 0.0;
  bool m_isIsotropic = true;
  std::vector<bool> m_isSolid;
  std::vector<std::pair<int, int>> m_layerIndices;
  std::vector<FLOAT> m_layerBounds;
  std::vector<std::vector<FLOAT>> m_layeredRadii;
  std::vector<InterpA> m_density, m_vpv, m_vsv, m_qkappa, m_qmu, m_vph, m_vsh, m_eta;
};

template <typename FLOAT, typename INTEGRAL>
ModelInput<FLOAT, INTEGRAL>::ModelInput(const std::string &pathToFile) {
  initializeNorms(EarthConstants<FLOAT>());
  loadModel(pathToFile);
}

template <typename FLOAT, typename INTEGRAL>
template <template <typename> class ParameterModel>
ModelInput<FLOAT, INTEGRAL>::ModelInput(const std::string &pathToFile,
                                         const ParameterModel<FLOAT> &modelConstants) {
  initializeNorms(modelConstants);
  loadModel(pathToFile);
}

template <typename FLOAT, typename INTEGRAL>
void ModelInput<FLOAT, INTEGRAL>::loadModel(const std::string &pathToFile) {
  const auto deck = planetmodel::ReadDeck(pathToFile, planetmodel::MINEOS());
  const auto header = planetmodel::ReadMineosHeader(deck.Header(), pathToFile);
  m_modelTitle = header.name;
  m_tref = header.tref / m_timeNorm;
  m_isIsotropic = header.ifanis == 0;
  const auto rawRadii = deck.Radius();
  const auto ranges = deck.Layers();
  const auto boundaries = deck.Boundaries();
  m_numLayers = static_cast<int>(ranges.size());
  m_layerBounds.reserve(boundaries.size());
  for (double radius : boundaries) m_layerBounds.push_back(radius / m_lengthNorm);
  m_layeredRadii.resize(ranges.size());
  m_layerIndices.resize(ranges.size());
  m_isSolid.resize(ranges.size());

  const auto& rho = deck.Column("rho");
  const auto isotropic = deck.HasColumn("vp");
  const auto& vpv = deck.Column(isotropic ? "vp" : "vpv");
  const auto& vsv = deck.Column(isotropic ? "vs" : "vsv");
  const auto& qk = deck.Column("qkappa");
  const auto& qm = deck.Column("qmu");
  if ((header.ifanis == 0) != isotropic) {
    throw std::invalid_argument(std::format(
        "{}: ifanis={} is inconsistent with the {}-column MINEOS rows",
        pathToFile, header.ifanis, isotropic ? "isotropic" : "TI"));
  }
  for (std::size_t row = 0; row < rawRadii.size(); ++row) {
    for (const auto [name, value] : {std::pair{"qkappa", qk[row]}, std::pair{"qmu", qm[row]}}) {
      if (!std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument(std::format("{}: row {} has invalid {}={}; expected finite Q >= 0",
                                                pathToFile, row + 1, name, value));
      }
    }
  }

  int start = 0;
  for (std::size_t layer = 0; layer < ranges.size(); ++layer) {
    const auto [begin, end] = ranges[layer];
    auto &radii = m_layeredRadii[layer];
    radii.reserve(end - begin);
    for (auto row = begin; row < end; ++row) radii.push_back(rawRadii[row] / m_lengthNorm);
    m_layerIndices[layer] = {start, start + static_cast<int>(end - begin) - 1};
    start += static_cast<int>(end - begin);
    m_isSolid[layer] = false;
    for (auto row = begin; row < end; ++row) {
      if (vsv[row] != 0.0 || (!isotropic && deck.Column("vsh")[row] != 0.0)) m_isSolid[layer] = true;
    }

    const auto lo = rawRadii[begin], hi = rawRadii[end - 1];
    const auto make = [&](const char* name, double scale) {
      const auto& source = deck.Column(name);
      std::vector<double> knots, values;
      knots.reserve(end - begin);
      values.reserve(end - begin);
      // Normalize before spline construction to preserve legacy floating-point evaluation.
      for (auto row = begin; row < end; ++row) {
        knots.push_back(rawRadii[row] / m_lengthNorm);
        values.push_back(source[row] / scale);
      }
      auto field = planetmodel::ThroughKnots(
          std::move(knots), std::move(values), planetmodel::KnotKind::Natural, 0.0);
      return InterpA(std::move(field), lo / m_lengthNorm, hi / m_lengthNorm);
    };
    m_density.push_back(make("rho", m_densityNorm));
    m_vpv.push_back(make(isotropic ? "vp" : "vpv", m_velocityNorm));
    m_vsv.push_back(make(isotropic ? "vs" : "vsv", m_velocityNorm));
    m_qkappa.push_back(make("qkappa", 1.0));
    m_qmu.push_back(make("qmu", 1.0));
    m_vph.push_back(make(isotropic ? "vp" : "vph", m_velocityNorm));
    m_vsh.push_back(make(isotropic ? "vs" : "vsh", m_velocityNorm));
    if (isotropic) {
      std::vector<double> knots;
      knots.reserve(end - begin);
      for (auto row = begin; row < end; ++row)
        knots.push_back(rawRadii[row] / m_lengthNorm);
      std::vector<double> ones(end - begin, 1.0);
      auto unitEta = planetmodel::ThroughKnots(
          std::move(knots),
          std::move(ones), planetmodel::KnotKind::Natural, 0.0);
      m_eta.emplace_back(std::move(unitEta), lo / m_lengthNorm, hi / m_lengthNorm);
    } else {
      m_eta.push_back(make("eta", 1.0));
    }
  }
}

} // namespace EarthModels

#endif // DSPECM1D_MODEL_INFO_MODEL_INPUT_H
