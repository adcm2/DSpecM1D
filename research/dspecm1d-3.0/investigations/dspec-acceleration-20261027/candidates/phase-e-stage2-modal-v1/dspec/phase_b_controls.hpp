#ifndef PHASE_B_CONTROLS_HPP
#define PHASE_B_CONTROLS_HPP

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace phaseb {

struct MeshSelection {
  int nq=0;
  bool finer=false;
  double relativeError=0.0;
  double timeNorm=0.0;
  std::complex<double> highestComplexOmega{};
  double omegaReference=0.0;
  double nominalStep=0.0;
  double selectedStep=0.0;
};

inline MeshSelection selectMesh(int nq,double relativeError,double timeNorm,
                                std::complex<double> highestComplexOmega,
                                bool finer) {
  if(nq<2||!(std::isfinite(relativeError)&&relativeError>0.0)||
     !(std::isfinite(timeNorm)&&timeNorm>0.0)||
     !std::isfinite(highestComplexOmega.real())||
     !std::isfinite(highestComplexOmega.imag()))
    throw std::invalid_argument("Invalid Phase B mesh-selection inputs.");
  const double omegaReference=std::abs(highestComplexOmega);
  if(!(std::isfinite(omegaReference)&&omegaReference>0.0))
    throw std::invalid_argument("Highest retained complex frequency must be finite and nonzero.");
  const double nominal=std::min(
      std::pow(100.0*relativeError,1.0/static_cast<double>(nq-1)) *
      (0.78*1000.0/timeNorm*2.0*std::acos(-1.0))/omegaReference,
      0.05);
  if(!(std::isfinite(nominal)&&nominal>0.0))
    throw std::runtime_error("Phase B nominal mesh step is invalid.");
  return {nq,finer,relativeError,timeNorm,highestComplexOmega,omegaReference,
          nominal,finer?nominal*0.5:nominal};
}

} // namespace phaseb

#endif
