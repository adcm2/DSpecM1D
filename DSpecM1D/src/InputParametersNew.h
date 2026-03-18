#ifndef DSPECM1D_INPUT_PARAMETERS_NEW_H
#define DSPECM1D_INPUT_PARAMETERS_NEW_H

#include <algorithm>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <utility>
#include <PlanetaryModel/All>
#include <SpectraSolver/FF>
#include "InputParser.h"
#include "SRInfo.h"
#include "SourceInfo.h"

// Transitional aggregate input object for spectra workflows.
// It intentionally coexists with InputParameters to keep adoption non-breaking.
class InputParametersNew {
private:
  using ModelType = decltype(EarthModels::ModelInput(
      std::declval<std::string>(), std::declval<prem_norm<double>>(), "true"));

  static std::string resolveModelPath(const std::string &paramPath,
                                      const InputParameters &params) {
    namespace fs = std::filesystem;

    fs::path inputPath(paramPath);
    fs::path dataDir = inputPath.parent_path().parent_path();
    fs::path modelPath = dataDir / params.earth_model();

    if (!fs::exists(modelPath)) {
      // Fallback: allow direct/absolute earth_model paths from the param file.
      modelPath = fs::path(params.earth_model());
    }

    if (!fs::exists(modelPath)) {
      throw std::runtime_error(
          "Could not resolve earth model path for InputParametersNew: " +
          modelPath.string());
    }

    return modelPath.string();
  }

  static ModelType loadEarthModel(const std::string &modelPath) {
    prem_norm<double> normClass;
    return EarthModels::ModelInput(modelPath, normClass, "true");
  }

  InputParameters m_params;
  std::string m_modelPath;
  ModelType m_model;
  SourceInfo::EarthquakeCMT m_cmt;
  SRInfo m_srInfo;
  SpectraSolver::FreqFull m_freqFull;
  double m_maxstep;
  int m_nq;
  int m_nskip;

public:
  explicit InputParametersNew(const std::string &paramPath, int nq = 5,
                              int nskip = 10, double maxstep = 0.05,
                              double df0 = 1.0, double wtb = 0.05,
                              double t1 = 0.0, int qex = 1)
      : m_params(paramPath), m_modelPath(resolveModelPath(paramPath, m_params)),
        m_model(loadEarthModel(m_modelPath)), m_cmt(m_params),
        m_srInfo(m_params),
        m_freqFull(m_params.f1(), m_params.f2(), m_params.f11(), m_params.f12(),
                   m_params.f21(), m_params.f22(), m_params.time_step_sec(),
                   m_params.t_out() / 60.0, df0, wtb, t1,
                   m_params.t_out() / 60.0, qex, m_model.TimeNorm()),
        m_maxstep(maxstep), m_nq(std::max(1, nq)), m_nskip(std::max(1, nskip)) {
  }

  InputParameters &inputParameters() { return m_params; }
  const InputParameters &inputParameters() const { return m_params; }

  SourceInfo::EarthquakeCMT &cmt() { return m_cmt; }
  const SourceInfo::EarthquakeCMT &cmt() const { return m_cmt; }

  SRInfo &srInfo() { return m_srInfo; }
  const SRInfo &srInfo() const { return m_srInfo; }

  ModelType &earthModel() { return m_model; }
  const ModelType &earthModel() const { return m_model; }
  const std::string &earthModelPath() const { return m_modelPath; }

  SpectraSolver::FreqFull &freqFull() { return m_freqFull; }
  const SpectraSolver::FreqFull &freqFull() const { return m_freqFull; }

  double timeNorm() const { return m_model.TimeNorm(); }
  double tref() const { return m_model.TREF(); }

  double normFactor() const {
    if (m_params.output_type() == 0)
      return m_model.LengthNorm();
    if (m_params.output_type() == 1)
      return m_model.LengthNorm() / m_model.TimeNorm();
    if (m_params.output_type() == 2)
      return m_model.LengthNorm() / (m_model.TimeNorm() * m_model.TimeNorm());
    return 1.0;
  }

  double maxstep() const { return m_maxstep; }
  void setMaxstep(double maxstep) { m_maxstep = maxstep; }

  int nq() const { return m_nq; }
  void setNq(int nq) { m_nq = std::max(1, nq); }

  int nskip() const { return m_nskip; }
  void setNskip(int nskip) { m_nskip = std::max(1, nskip); }
};

#endif   // DSPECM1D_INPUT_PARAMETERS_NEW_H
