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
  double rad_source = m_mesh.PR() - 1000.0 * cmt.Depth() / m_lengthNorm;

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
SEM::receiverElements(InputParameters &param) const {
  std::vector<int> receiver_elems;
  double rad_receiver =
      m_mesh.PR() - 1000.0 * param.receiver_depth() / m_lengthNorm;

  for (int idx = 0; idx < m_mesh.NE(); ++idx) {
    if ((m_mesh.ELR(idx) <= rad_receiver) && (m_mesh.EUR(idx) >= rad_receiver)) {
      receiver_elems.push_back(idx);
    }
  }
  return receiver_elems;
};

}   // namespace Full1D

#endif   // SEM_LTG_H
