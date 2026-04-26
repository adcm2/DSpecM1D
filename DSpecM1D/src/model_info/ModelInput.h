#ifndef DSPECM1D_MODEL_INFO_MODEL_INPUT_H
#define DSPECM1D_MODEL_INFO_MODEL_INPUT_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <Interpolation/CubicSpline>

namespace EarthModels {

template <typename FLOAT = double>
class EarthConstants {
public:
  using value_type = FLOAT;

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

private:
  const FLOAT m_lengthNorm = 6.371 * std::pow(10.0, 6.0);
  const FLOAT m_massNorm = 5.972 * std::pow(10.0, 24.0);
  const FLOAT m_timeNorm = 3600.0;

  const FLOAT m_densityNorm = m_massNorm / std::pow(m_lengthNorm, 3.0);
  const FLOAT m_inertiaNorm = m_massNorm * std::pow(m_lengthNorm, 2.0);
  const FLOAT m_velocityNorm = m_lengthNorm / m_timeNorm;
  const FLOAT m_accelerationNorm =
      m_lengthNorm / std::pow(m_timeNorm, 2.0);
  const FLOAT m_forceNorm =
      m_massNorm * m_lengthNorm / std::pow(m_timeNorm, 2.0);
  const FLOAT m_stressNorm =
      m_massNorm / (std::pow(m_timeNorm, 2.0) * m_lengthNorm);
  const FLOAT m_gravitationalConstant =
      std::pow(m_lengthNorm, 3.0) / (m_massNorm * std::pow(m_timeNorm, 2.0));
};

template <typename FLOAT = double, typename INTEGRAL = int>
class ModelInput {
public:
  using size_type = INTEGRAL;
  using value_type = FLOAT;
  using myvector = std::vector<FLOAT>;
  using myiter = typename myvector::iterator;
  using InterpA = Interpolation::CubicSpline<myiter, myiter>;

  ModelInput() = default;
  explicit ModelInput(const std::string &pathToFile);

  template <template <typename> class ParameterModel>
  ModelInput(const std::string &pathToFile,
             const ParameterModel<FLOAT> &modelConstants);

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
  auto TREF() const { return m_tref; }

  int NumberOfLayers() const { return m_numLayers; }
  int LayerLowerIndex(int i) const { return m_layerIndices[i][0]; }
  int LayerUpperIndex(int i) const { return m_layerIndices[i][1]; }
  int LayerIndexDifference(int i) const {
    return m_layerIndices[i][1] - m_layerIndices[i][0];
  }
  auto LayerRadii() const { return m_layeredRadii; }
  auto LayerRadii(int i) const { return m_layeredRadii[i]; }
  auto LayerRadiiNumber(int i) const { return m_layeredRadii[i].size(); }

  FLOAT LowerRadius(INTEGRAL i) const {
    if (i < 0) {
      throw std::invalid_argument("Negative layer index");
    } else if (i > m_numLayers - 1) {
      throw std::invalid_argument("Layer index greater than number of layers");
    }
    return m_layerBounds[i];
  }

  FLOAT UpperRadius(INTEGRAL i) const {
    if (i < 0) {
      throw std::invalid_argument("Negative layer index");
    } else if (i > m_numLayers - 1) {
      throw std::invalid_argument("Layer index greater than number of layers");
    }
    return m_layerBounds[i + 1];
  }

  FLOAT OuterRadius() const { return m_layerBounds[m_numLayers]; }

  bool IsIsotropic() const { return m_isIsotropic; }
  bool IsSolid(INTEGRAL i) const { return m_isSolid[i]; }
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
    auto fn = [i, this](FLOAT x) { return std::sqrt(Mu(i)(x) / Density(i)(x)); };
    return fn;
  }

  auto VP(INTEGRAL i) const {
    auto fn = [i, this](FLOAT x) {
      return std::sqrt((Kappa(i)(x) + 4.0 / 3.0 * Mu(i)(x)) / Density(i)(x));
    };
    return fn;
  }

  auto A(INTEGRAL i) const {
    auto fn = [i, this](FLOAT x) {
      return Density(i)(x) * VPH(i)(x) * VPH(i)(x);
    };
    return fn;
  }

  auto C(INTEGRAL i) const {
    auto fn = [i, this](FLOAT x) {
      return Density(i)(x) * VPV(i)(x) * VPV(i)(x);
    };
    return fn;
  }

  auto N(INTEGRAL i) const {
    auto fn = [i, this](FLOAT x) {
      return Density(i)(x) * VSH(i)(x) * VSH(i)(x);
    };
    return fn;
  }

  auto L(INTEGRAL i) const {
    auto fn = [i, this](FLOAT x) {
      return Density(i)(x) * VSV(i)(x) * VSV(i)(x);
    };
    return fn;
  }

  auto F(INTEGRAL i) const {
    auto fn = [i, this](FLOAT x) { return Eta(i)(x) * (A(i)(x) - 2 * L(i)(x)); };
    return fn;
  }

  auto Kappa(INTEGRAL i) const {
    auto fn = [i, this](FLOAT x) {
      return (C(i)(x) + 4.0 * (A(i)(x) - N(i)(x) + F(i)(x))) / 9.0;
    };
    return fn;
  }

  auto Mu(INTEGRAL i) const {
    auto fn = [i, this](FLOAT x) {
      return (C(i)(x) + A(i)(x) + 6.0 * L(i)(x) + 5.0 * N(i)(x) -
              2.0 * F(i)(x)) /
             15.0;
    };
    return fn;
  }

private:
  using VecDouble = std::vector<double>;
  using LayeredValues = std::vector<VecDouble>;

  template <typename LayerContainer>
  auto checkedLayer(const LayerContainer &values, INTEGRAL i) const {
    if (i < 0) {
      throw std::invalid_argument("Negative layer index");
    } else if (i > m_numLayers - 1) {
      throw std::invalid_argument("Layer index greater than number of layers");
    }
    return values[i];
  }

  template <template <typename> class ParameterModel>
  void initializeNorms(const ParameterModel<FLOAT> &modelConstants) {
    m_lengthNorm = modelConstants.LengthNorm();
    m_massNorm = modelConstants.MassNorm();
    m_timeNorm = modelConstants.TimeNorm();
    m_densityNorm = modelConstants.MassNorm() /
                    std::pow(modelConstants.LengthNorm(), 3.0);
    m_velocityNorm = modelConstants.LengthNorm() / modelConstants.TimeNorm();
    m_accelerationNorm =
        modelConstants.LengthNorm() / std::pow(modelConstants.TimeNorm(), 2.0);
    m_forceNorm = modelConstants.MassNorm() * modelConstants.LengthNorm() /
                  std::pow(modelConstants.TimeNorm(), 2.0);
    m_stressNorm =
        modelConstants.MassNorm() /
        (modelConstants.LengthNorm() * std::pow(modelConstants.TimeNorm(), 2.0));
    m_inertiaNorm =
        modelConstants.MassNorm() * std::pow(modelConstants.LengthNorm(), 2.0);
    m_gravitationalConstant =
        std::pow(modelConstants.LengthNorm(), 3.0) /
        (modelConstants.MassNorm() * std::pow(modelConstants.TimeNorm(), 2.0));
  }

  void loadModel(const std::string &pathToFile,
                 bool respectIfanisColumns,
                 bool printModelKind);

  std::string m_modelTitle;
  int m_ifanis = 0;
  int m_ifdeck = 0;
  int m_numNodes = 0;
  int m_nic = 0;
  int m_noc = 0;
  int m_numLayers = 0;
  double m_tref = 0.0;

  FLOAT m_lengthNorm = 0.0;
  FLOAT m_massNorm = 0.0;
  FLOAT m_timeNorm = 0.0;
  FLOAT m_densityNorm = 0.0;
  FLOAT m_inertiaNorm = 0.0;
  FLOAT m_velocityNorm = 0.0;
  FLOAT m_accelerationNorm = 0.0;
  FLOAT m_forceNorm = 0.0;
  FLOAT m_stressNorm = 0.0;
  FLOAT m_gravitationalConstant = 0.0;

  bool m_isIsotropic = true;
  std::vector<bool> m_isSolid = std::vector<bool>(1, true);
  std::vector<std::vector<int>> m_layerIndices;

  VecDouble m_layerBounds;
  LayeredValues m_layeredRadii = LayeredValues(1, VecDouble());
  LayeredValues m_layeredRho = LayeredValues(1, VecDouble());
  LayeredValues m_layeredVpv = LayeredValues(1, VecDouble());
  LayeredValues m_layeredVsv = LayeredValues(1, VecDouble());
  LayeredValues m_layeredQkappa = LayeredValues(1, VecDouble());
  LayeredValues m_layeredQmu = LayeredValues(1, VecDouble());
  LayeredValues m_layeredVph = LayeredValues(1, VecDouble());
  LayeredValues m_layeredVsh = LayeredValues(1, VecDouble());
  LayeredValues m_layeredEta = LayeredValues(1, VecDouble());

  std::vector<InterpA> m_density;
  std::vector<InterpA> m_vpv;
  std::vector<InterpA> m_vsv;
  std::vector<InterpA> m_qkappa;
  std::vector<InterpA> m_qmu;
  std::vector<InterpA> m_vph;
  std::vector<InterpA> m_vsh;
  std::vector<InterpA> m_eta;
};

template <typename FLOAT, typename INTEGRAL>
ModelInput<FLOAT, INTEGRAL>::ModelInput(const std::string &pathToFile) {
  initializeNorms(EarthConstants<FLOAT>());
  loadModel(pathToFile, false, false);
}

template <typename FLOAT, typename INTEGRAL>
template <template <typename> class ParameterModel>
ModelInput<FLOAT, INTEGRAL>::ModelInput(
    const std::string &pathToFile, const ParameterModel<FLOAT> &modelConstants) {
  initializeNorms(modelConstants);
  loadModel(pathToFile, true, false);
}

template <typename FLOAT, typename INTEGRAL>
void
ModelInput<FLOAT, INTEGRAL>::loadModel(const std::string &pathToFile,
                                       bool respectIfanisColumns,
                                       bool printModelKind) {
  std::fstream modelFile(pathToFile, std::ios::in);
  if (!modelFile.is_open()) {
    throw std::runtime_error("Could not open Earth model file: " + pathToFile);
  }

  std::getline(modelFile, m_modelTitle);

  modelFile >> m_ifanis >> m_tref >> m_ifdeck;
  modelFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  m_tref *= 1.0 / m_timeNorm;

  modelFile >> m_numNodes >> m_nic >> m_noc;
  modelFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  if (printModelKind) {
    if (m_ifanis) {
      std::cout << "Anisotropic model loaded from file: " << pathToFile << "\n";
    } else {
      std::cout << "Isotropic model loaded from file: " << pathToFile << "\n";
    }
  }

  int layerNumber = 0;
  int outerIndex = 0;
  int innerIndex = 0;

  while (outerIndex < m_numNodes) {
    double radius = 0.0;
    double rho = 0.0;
    double vpv = 0.0;
    double vsv = 0.0;
    double qkappa = 0.0;
    double qshear = 0.0;
    double vph = 0.0;
    double vsh = 0.0;
    double eta = 0.0;

    if (!respectIfanisColumns || m_ifanis) {
      modelFile >> radius >> rho >> vpv >> vsv >> qkappa >> qshear >> vph >>
          vsh >> eta;
    } else {
      modelFile >> radius >> rho >> vpv >> vsv >> qkappa >> qshear;
      vph = vpv;
      vsh = vsv;
      eta = 1.0;
    }

    const double normalizedRadius = radius / m_lengthNorm;
    if (innerIndex > 0 &&
        normalizedRadius == m_layeredRadii[layerNumber][innerIndex - 1]) {
      m_layeredRadii.push_back({normalizedRadius});
      m_layeredRho.push_back({rho / m_densityNorm});
      m_layeredVpv.push_back({vpv / m_velocityNorm});
      m_layeredVsv.push_back({vsv / m_velocityNorm});
      m_layeredQkappa.push_back({qkappa});
      m_layeredQmu.push_back({qshear});
      m_layeredVph.push_back({vph / m_velocityNorm});
      m_layeredVsh.push_back({vsh / m_velocityNorm});
      m_layeredEta.push_back({eta});

      if (m_isIsotropic && vpv != vsv) {
        m_isIsotropic = false;
      }

      m_isSolid.push_back(!(vsv == 0.0 && vsh == 0.0));
      innerIndex = 0;
      ++layerNumber;
    } else {
      m_layeredRadii[layerNumber].push_back(normalizedRadius);
      m_layeredRho[layerNumber].push_back(rho / m_densityNorm);
      m_layeredVpv[layerNumber].push_back(vpv / m_velocityNorm);
      m_layeredVsv[layerNumber].push_back(vsv / m_velocityNorm);
      m_layeredQkappa[layerNumber].push_back(qkappa);
      m_layeredQmu[layerNumber].push_back(qshear);
      m_layeredVph[layerNumber].push_back(vph / m_velocityNorm);
      m_layeredVsh[layerNumber].push_back(vsh / m_velocityNorm);
      m_layeredEta[layerNumber].push_back(eta);

      if (outerIndex == 0 && vsv == 0.0 && vsh == 0.0) {
        m_isSolid[layerNumber] = false;
      }
    }

    modelFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    ++innerIndex;
    ++outerIndex;
  }

  m_numLayers = static_cast<int>(m_layeredRadii.size());
  m_layerBounds.reserve(m_numLayers + 1);
  m_layerBounds.push_back(0.0);

  m_layerIndices =
      std::vector<std::vector<int>>(m_numLayers, std::vector<int>(2, 0));
  m_layerIndices[0][0] = 0;

  for (int idx = 0; idx < m_numLayers; ++idx) {
    m_layerBounds.push_back(m_layeredRadii[idx].back());
    if (idx != 0) {
      m_layerIndices[idx][0] = m_layerIndices[idx - 1][1] + 1;
    }
    m_layerIndices[idx][1] =
        m_layerIndices[idx][0] + static_cast<int>(m_layeredRadii[idx].size()) - 1;
  }

  for (int idx = 0; idx < m_numLayers; ++idx) {
    auto it1 = m_layeredRadii[idx].begin();
    auto it2 = m_layeredRadii[idx].end();

    m_density.push_back(InterpA(it1, it2, m_layeredRho[idx].begin()));
    m_vpv.push_back(InterpA(it1, it2, m_layeredVpv[idx].begin()));
    m_vsv.push_back(InterpA(it1, it2, m_layeredVsv[idx].begin()));
    m_qkappa.push_back(InterpA(it1, it2, m_layeredQkappa[idx].begin()));
    m_qmu.push_back(InterpA(it1, it2, m_layeredQmu[idx].begin()));
    m_vph.push_back(InterpA(it1, it2, m_layeredVph[idx].begin()));
    m_vsh.push_back(InterpA(it1, it2, m_layeredVsh[idx].begin()));
    m_eta.push_back(InterpA(it1, it2, m_layeredEta[idx].begin()));
  }
}

} // namespace EarthModels

#endif // DSPECM1D_MODEL_INFO_MODEL_INPUT_H
