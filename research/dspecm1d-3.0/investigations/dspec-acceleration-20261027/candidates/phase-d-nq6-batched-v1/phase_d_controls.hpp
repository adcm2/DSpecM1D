#ifndef DSPECM1D_PHASE_D_CONTROLS_HPP
#define DSPECM1D_PHASE_D_CONTROLS_HPP

#include "phase_b_controls.hpp"
#include <cstdlib>
#include <string>

namespace phased {
inline double requestedAccuracy() {
  const char *value=std::getenv("DSPEC_PHASE_D_REQUESTED_RELATIVE_ACCURACY");
  const double result=value?std::stod(value):1e-3;
  if(!(std::isfinite(result)&&result>0.0&&result<1.0))
    throw std::invalid_argument("Phase D requested relative accuracy must be in (0,1).");
  return result;
}
inline double meshScale() {
  const char *value=std::getenv("DSPEC_PHASE_D_MESH_SCALE");
  const double result=value?std::stod(value):1.0;
  if(result!=1.0&&result!=0.9&&result!=0.8)
    throw std::invalid_argument("Phase D mesh scale must be exactly 1.0, 0.9, or 0.8.");
  return result;
}
}
#endif
