#ifndef SEM_MATRICES_H
#define SEM_MATRICES_H

#include "SEM.h"

namespace Full1D {

Eigen::SparseMatrix<double>
SEM::hS(int idxl) const {
  auto k2 = idxl * (idxl + 1);
  Eigen::SparseMatrix<double> tmp = m_vecKeSBase[0];
  tmp += m_vecKeSBase[1] * k2;
  tmp += m_vecKeSBase[2] * k2 * k2;

  auto idx_pp = this->ltgS(2, m_mesh.NE() - 1, m_nq - 1);
  double rpb = m_mesh.NodeRadius(m_mesh.NE() - 1, m_nq - 1);
  tmp.coeffRef(idx_pp, idx_pp) += rpb * (idxl + 1) / (4.0 * pi_db * bigg_db);

  return tmp;
};

Eigen::SparseMatrix<double>
SEM::hSa(int idxl) const {
  auto k2 = idxl * (idxl + 1);
  Eigen::SparseMatrix<double> tmp = m_vecKeSAtten[0];
  tmp += m_vecKeSAtten[1] * k2;
  tmp += m_vecKeSAtten[2] * k2 * k2;
  return tmp;
};

Eigen::SparseMatrix<double>
SEM::pS(int idxl) const {
  auto k2 = idxl * (idxl + 1);
  Eigen::SparseMatrix<double> tmp = m_vecInSBase[0];
  tmp += m_vecInSBase[1] * k2;
  return tmp;
};

Eigen::SparseMatrix<double>
SEM::hTk(int idxl) const {
  auto k2 = idxl * (idxl + 1);
  Eigen::SparseMatrix<double> tmp = m_vecKeTBase[0] * k2;
  tmp += m_vecKeTBase[1] * k2 * k2;
  return tmp;
};

Eigen::SparseMatrix<double>
SEM::hTa(int idxl) const {
  auto k2 = idxl * (idxl + 1);
  Eigen::SparseMatrix<double> tmp = m_vecKeTAtten[0] * k2;
  tmp += m_vecKeTAtten[1] * k2 * k2;
  return tmp;
};

Eigen::SparseMatrix<double>
SEM::pTk(int idxl) const {
  auto k2 = idxl * (idxl + 1);
  Eigen::SparseMatrix<double> tmp = m_matInTBase * k2;
  return tmp;
};

}   // namespace Full1D

#endif   // SEM_MATRICES_H
