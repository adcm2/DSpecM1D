#ifndef SEM_MATRICES_H
#define SEM_MATRICES_H

#include "SEM.h"

namespace Full1D {

Eigen::SparseMatrix<double>
SEM::hS(int idxl) const {
  auto k2 = idxl * (idxl + 1);
  Eigen::SparseMatrix<double> tmp = vec_ke_s_base[0];
  tmp += vec_ke_s_base[1] * k2;
  tmp += vec_ke_s_base[2] * k2 * k2;

  auto idx_pp = this->ltgS(2, _mesh.NE() - 1, _NQ - 1);
  double rpb = _mesh.NodeRadius(_mesh.NE() - 1, _NQ - 1);
  tmp.coeffRef(idx_pp, idx_pp) += rpb * (idxl + 1) / (4.0 * pi_db * bigg_db);

  return tmp;
};

Eigen::SparseMatrix<double>
SEM::hSa(int idxl) const {
  auto k2 = idxl * (idxl + 1);
  Eigen::SparseMatrix<double> tmp = vec_ke_s_atten[0];
  tmp += vec_ke_s_atten[1] * k2;
  tmp += vec_ke_s_atten[2] * k2 * k2;
  return tmp;
};

Eigen::SparseMatrix<double>
SEM::pS(int idxl) const {
  auto k2 = idxl * (idxl + 1);
  Eigen::SparseMatrix<double> tmp = vec_in_s_base[0];
  tmp += vec_in_s_base[1] * k2;
  return tmp;
};

Eigen::SparseMatrix<double>
SEM::hTk(int idxl) const {
  auto k2 = idxl * (idxl + 1);
  Eigen::SparseMatrix<double> tmp = vec_ke_t_base[0] * k2;
  tmp += vec_ke_t_base[1] * k2 * k2;
  return tmp;
};

Eigen::SparseMatrix<double>
SEM::hTa(int idxl) const {
  auto k2 = idxl * (idxl + 1);
  Eigen::SparseMatrix<double> tmp = vec_ke_t_atten[0] * k2;
  tmp += vec_ke_t_atten[1] * k2 * k2;
  return tmp;
};

Eigen::SparseMatrix<double>
SEM::pTk(int idxl) const {
  auto k2 = idxl * (idxl + 1);
  Eigen::SparseMatrix<double> tmp = mat_in_t_base * k2;
  return tmp;
};

}   // namespace Full1D

#endif   // SEM_MATRICES_H
