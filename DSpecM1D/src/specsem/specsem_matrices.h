#ifndef SPECSEM_MATRICES_H
#define SPECSEM_MATRICES_H

#include "specsem.h"

namespace Full1D {

Eigen::SparseMatrix<double>
specsem::H_S(int idxl) const {
  auto k2 = idxl * (idxl + 1);
  Eigen::SparseMatrix<double> tmp = vec_ke_s_base[0];
  tmp += vec_ke_s_base[1] * k2;
  tmp += vec_ke_s_base[2] * k2 * k2;

  auto idx_pp = this->LtG_S(2, _mesh.NE() - 1, _NQ - 1);
  double rpb = _mesh.NodeRadius(_mesh.NE() - 1, _NQ - 1);
  tmp.coeffRef(idx_pp, idx_pp) += rpb * (idxl + 1) / (4.0 * pi_db * bigg_db);

  return tmp;
};

Eigen::SparseMatrix<double>
specsem::H_SA(int idxl) const {
  auto k2 = idxl * (idxl + 1);
  Eigen::SparseMatrix<double> tmp = vec_ke_s_atten[0];
  tmp += vec_ke_s_atten[1] * k2;
  tmp += vec_ke_s_atten[2] * k2 * k2;
  return tmp;
};

Eigen::SparseMatrix<double>
specsem::P_S(int idxl) const {
  auto k2 = idxl * (idxl + 1);
  Eigen::SparseMatrix<double> tmp = vec_in_s_base[0];
  tmp += vec_in_s_base[1] * k2;
  return tmp;
};

Eigen::SparseMatrix<double>
specsem::H_TK(int idxl) const {
  auto k2 = idxl * (idxl + 1);
  Eigen::SparseMatrix<double> tmp = vec_ke_t_base[0] * k2;
  tmp += vec_ke_t_base[1] * k2 * k2;
  return tmp;
};

Eigen::SparseMatrix<double>
specsem::H_TA(int idxl) const {
  auto k2 = idxl * (idxl + 1);
  Eigen::SparseMatrix<double> tmp = vec_ke_t_atten[0] * k2;
  tmp += vec_ke_t_atten[1] * k2 * k2;
  return tmp;
};

Eigen::SparseMatrix<double>
specsem::P_TK(int idxl) const {
  auto k2 = idxl * (idxl + 1);
  Eigen::SparseMatrix<double> tmp = mat_in_t_base * k2;
  return tmp;
};

}   // namespace Full1D

#endif   // SPECSEM_MATRICES_H
