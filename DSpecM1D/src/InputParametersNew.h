#ifndef DSPECM1D_INPUT_PARAMETERS_NEW_H
#define DSPECM1D_INPUT_PARAMETERS_NEW_H

#include <algorithm>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <PlanetaryModel/All>
#include <SpectraSolver/FF>
#include "InputParser.h"
#include "SRInfo.h"
#include "SourceInfo.h"

// Transitional aggregate input object for spectra workflows.
// It intentionally coexists with InputParameters to keep adoption non-breaking.
class InputParametersNew {
private:
  struct ModelScales {
    double timeNorm;
    double tref;
  };

  static ModelScales inferModelScales(const std::string &paramPath,
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

    prem_norm<double> normClass;
    auto model = EarthModels::ModelInput(modelPath.string(), normClass, "true");
    return {model.TimeNorm(), model.TREF()};
  }

  InputParameters m_params;
  ModelScales m_scales;
  SourceInfo::EarthquakeCMT m_cmt;
  SRInfo m_srInfo;
  SpectraSolver::FreqFull m_freqFull;
  int m_nskip;

public:
  explicit InputParametersNew(const std::string &paramPath, int nskip = 10,
                              double df0 = 1.0, double wtb = 0.05,
                              double t1 = 0.0, int qex = 1)
      : m_params(paramPath), m_scales(inferModelScales(paramPath, m_params)),
        m_cmt(m_params), m_srInfo(m_params),
        m_freqFull(m_params.f1(), m_params.f2(), m_params.f11(), m_params.f12(),
                   m_params.f21(), m_params.f22(), m_params.time_step_sec(),
                   m_params.t_out() / 60.0, df0, wtb, t1,
                   m_params.t_out() / 60.0, qex, m_scales.timeNorm),
        m_nskip(std::max(1, nskip)) {}

  InputParameters &inputParameters() { return m_params; }
  const InputParameters &inputParameters() const { return m_params; }

  SourceInfo::EarthquakeCMT &cmt() { return m_cmt; }
  const SourceInfo::EarthquakeCMT &cmt() const { return m_cmt; }

  SRInfo &srInfo() { return m_srInfo; }
  const SRInfo &srInfo() const { return m_srInfo; }

  SpectraSolver::FreqFull &freqFull() { return m_freqFull; }
  const SpectraSolver::FreqFull &freqFull() const { return m_freqFull; }

  double timeNorm() const { return m_scales.timeNorm; }
  double tref() const { return m_scales.tref; }

  int nskip() const { return m_nskip; }
  void setNskip(int nskip) { m_nskip = std::max(1, nskip); }
};

#endif   // DSPECM1D_INPUT_PARAMETERS_NEW_H
