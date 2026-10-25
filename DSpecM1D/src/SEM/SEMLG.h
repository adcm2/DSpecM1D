#ifndef SEM_LTG_H
#define SEM_LTG_H

#include "SEM.h"

namespace Full1D {

auto
SEM::ltgS(int neig, int idx_e, int idx_n) const {
  assert((neig >= 0) && (neig <= 2) &&
         "Error: neig must be 0, 1 or 2 in ltgS");
  assert((idx_e >= 0) && (idx_e < m_mesh.NE()) &&
         "Error: idx_e out of range in ltgS");
  assert((idx_n >= 0) && (idx_n < m_mesh.NN()) &&
         "Error: idx_n out of range in ltgS");
  std::size_t retval = 0;

  int offset_val = m_vecOffset[idx_e];
  if (idx_n == 0) {
    if ((std::find(m_fsb.begin(), m_fsb.end(), idx_e - 1) != m_fsb.end())) {
      if (neig == 1) {
        offset_val += 1;
      } else {
        offset_val -= 1;
      }
    }
  }
  retval = (3 * idx_e * (m_mesh.NN() - 1) + idx_n * 3 + neig) + offset_val;
  return retval;
};

auto
SEM::ltgSC(int field, int idx_e, int idx_n) const {
  assert((field >= 0) && (field <= 1) &&
         "Error: field must be 0 (U) or 1 (V) in ltgSC");
  assert((idx_e >= 0) && (idx_e < m_mesh.NE()) &&
         "Error: idx_e out of range in ltgSC");
  assert((idx_n >= 0) && (idx_n < m_mesh.NN()) &&
         "Error: idx_n out of range in ltgSC");

  int offset_val = m_vecOffset[idx_e];
  if (field == 0 && idx_n == 0 &&
      std::find(m_fsb.begin(), m_fsb.end(), idx_e - 1) != m_fsb.end())
    offset_val -= 1;

  return 2 * idx_e * (m_mesh.NN() - 1) + 2 * idx_n + field + offset_val;
};

auto
SEM::ltgR(int neig, int idx_e, int idx_n) const {
  assert((neig >= 0) && (neig < 2) && "Error: neig must be 0 or 1 in ltgR");
  assert((idx_e >= 0) && (idx_e < m_mesh.NE()) &&
         "Error: idx_e out of range in ltgR");
  assert((idx_n >= 0) && (idx_n < m_mesh.NN()) &&
         "Error: idx_n out of range in ltgR");

  std::size_t retval = 2 * idx_e * (m_mesh.NN() - 1);
  retval += idx_n * 2;
  retval += neig;
  return retval;
};

auto
SEM::ltgT(int idx_e, int idx_n) const {
  assert((idx_e >= m_el) && (idx_e < m_eu) && "idxe out of range in ltgT");
  assert((idx_n >= 0) && (idx_n < m_mesh.NN()) &&
         "Error: idx_n out of range in ltgT");
  std::size_t retval = (idx_e - m_el) * (m_mesh.NN() - 1) + idx_n;
  return retval;
};

auto
SEM::sourceElement(SourceInfo::EarthquakeCMT &cmt) const {
  if (sourceAligned() && !sourceMatches(cmt))
    throw std::invalid_argument(
        "Source depth differs from the radius used to align this SEM.");
  const double rad_source = sourceAligned()
      ? m_sourceRadius
      : m_mesh.PR() - 1000.0 * cmt.Depth() / m_lengthNorm;

  int idxout = m_mesh.NE() - 1;
  for (int idx = 0; idx < m_mesh.NE(); ++idx) {
    if ((m_mesh.ELR(idx) <= rad_source) && (m_mesh.EUR(idx) > rad_source)) {
      idxout = idx;
      break;
    };
  };
  return idxout;
};

auto
SEM::receiverElements(InputParameters &param, bool toroidal) const {
  if (param.num_receivers() == 0)
    throw std::invalid_argument("At least one receiver is required.");
  std::vector<int> containing_elems;
  double rad_receiver =
      m_mesh.PR() - 1000.0 * param.receiver_depth() / m_lengthNorm;

  for (int idx = 0; idx < m_mesh.NE(); ++idx) {
    if ((m_mesh.ELR(idx) <= rad_receiver) && (m_mesh.EUR(idx) >= rad_receiver)) {
      containing_elems.push_back(idx);
    }
  }

  if (containing_elems.empty())
    throw std::out_of_range("Receiver lies outside the SEM radial mesh.");

  // At a fluid-solid interface the receiver belongs to the physical rock
  // side. Preserve the established element selection at ordinary boundaries.
  bool has_solid = false;
  bool has_fluid = false;
  for (int idx : containing_elems) {
    has_fluid = has_fluid || m_vecFluid[idx] != 0;
    has_solid = has_solid || m_vecFluid[idx] == 0;
  }
  std::vector<int> receiver_elems;
  for (int idx : containing_elems) {
    if (!(has_solid && has_fluid) || m_vecFluid[idx] == 0)
      receiver_elems.push_back(idx);
  }
  if (toroidal) {
    for (int idx : receiver_elems) {
      if (idx < m_el || idx >= m_eu)
        throw std::invalid_argument(
            "Toroidal receiver lies outside the solid SEM range.");
    }
  }
  return receiver_elems;
};

}   // namespace Full1D

#endif   // SEM_LTG_H
