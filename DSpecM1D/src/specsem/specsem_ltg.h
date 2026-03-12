#ifndef SPECSEM_LTG_H
#define SPECSEM_LTG_H

#include "specsem.h"

namespace Full1D {

auto
specsem::LtG_S(int neig, int idx_e, int idx_n) const {
  assert((neig >= 0) && (neig <= 2) &&
         "Error: neig must be 0, 1 or 2 in LtG_S");
  assert((idx_e >= 0) && (idx_e < _mesh.NE()) &&
         "Error: idx_e out of range in LtG_S");
  assert((idx_n >= 0) && (idx_n < _mesh.NN()) &&
         "Error: idx_n out of range in LtG_S");
  std::size_t retval = 0;

  int offset_val = vec_offset[idx_e];
  if (idx_n == 0) {
    if ((std::find(fsb.begin(), fsb.end(), idx_e - 1) != fsb.end())) {
      if (neig == 1) {
        offset_val += 1;
      } else {
        offset_val -= 1;
      }
    }
  }
  retval = (3 * idx_e * (_mesh.NN() - 1) + idx_n * 3 + neig) + offset_val;
  return retval;
};

auto
specsem::LtG_R(int neig, int idx_e, int idx_n) const {
  assert((neig >= 0) && (neig < 2) && "Error: neig must be 0 or 1 in LtG_R");
  assert((idx_e >= 0) && (idx_e < _mesh.NE()) &&
         "Error: idx_e out of range in LtG_R");
  assert((idx_n >= 0) && (idx_n < _mesh.NN()) &&
         "Error: idx_n out of range in LtG_R");

  std::size_t retval = 2 * idx_e * (_mesh.NN() - 1);
  retval += idx_n * 2;
  retval += neig;
  return retval;
};

auto
specsem::LtG_T(int idx_e, int idx_n) const {
  assert((idx_e >= _el) && (idx_e < _eu) && "idxe out of range in LtG_T");
  assert((idx_n >= 0) && (idx_n < _mesh.NN()) &&
         "Error: idx_n out of range in LtG_T");
  std::size_t retval = (idx_e - _el) * (_mesh.NN() - 1) + idx_n;
  return retval;
};

auto
specsem::Source_Element(SourceInfo::EarthquakeCMT &cmt) const {
  double rad_source = _mesh.PR() - 1000.0 * cmt.Depth() / _length_norm;

  int idxout = _mesh.NE() - 1;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) > rad_source)) {
      idxout = idx;
      break;
    };
  };
  return idxout;
};

auto
specsem::Receiver_Elements(InputParameters &param) const {
  std::vector<int> receiver_elems;
  double rad_receiver =
      _mesh.PR() - 1000.0 * param.receiver_depth() / _length_norm;

  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_receiver) && (_mesh.EUR(idx) >= rad_receiver)) {
      receiver_elems.push_back(idx);
    }
  }
  return receiver_elems;
};

}   // namespace Full1D

#endif   // SPECSEM_LTG_H
