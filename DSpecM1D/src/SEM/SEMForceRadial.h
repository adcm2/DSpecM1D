#ifndef SEM_FORCE_RADIAL_H
#define SEM_FORCE_RADIAL_H

#include "SEM.h"

namespace Full1D {

Eigen::MatrixXcd
SEM::calculateForceR(SourceInfo::EarthquakeCMT &cmt) {
  int nq = m_mesh.NN();
  totlen = this->ltgR(1, m_mesh.NE() - 1, nq - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totlen, 1);

  double rad_source = _SourceRadius(cmt);

  auto wigdmat2 =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(0, 0, 0, 0.0);
  Complex y0c = wigdmat2[0][0, 0];

  double lprefac = 1.0;

  for (int idx = 0; idx < m_mesh.NE(); ++idx) {
    if ((m_mesh.ELR(idx) <= rad_source) && (m_mesh.EUR(idx) > rad_source)) {
      stdvec vec_nodes(m_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < m_mesh.NN(); ++idxn)
        vec_nodes[idxn] = m_mesh.NodeRadius(idx, idxn);
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

      for (int idxq = 0; idxq < nq; ++idxq) {
        auto w_val = pleg(idxq, rad_source) / rad_source;
        auto w_deriv = pleg.Derivative(idxq, rad_source);
        auto idx_u = this->ltgR(0, idx, idxq);

        Complex tmp_pm = 2.0 * cmt.MCmp() * y0c;
        Complex tmp_u = -w_val * tmp_pm + cmt.MC00() * w_deriv * y0c;

        vec_lforce(idx_u, 0) = tmp_u * lprefac;
      };
    };
  };

  vec_lforce *= (1.0 / m_momentNorm);
  return vec_lforce;
};

Eigen::MatrixXcd
SEM::calculateForceRedR(SourceInfo::EarthquakeCMT &cmt) {
  int nq = m_mesh.NN();
  totlen = this->ltgR(1, m_mesh.NE() - 1, nq - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totlen, 1);

  double rad_source = _SourceRadius(cmt);

  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(0, 0, 0, 0.0);
  Complex y0c = wigdmat[0][0, 0];

  for (int idx = 0; idx < m_mesh.NE(); ++idx) {
    if ((m_mesh.ELR(idx) <= rad_source) && (m_mesh.EUR(idx) > rad_source)) {
      stdvec vec_nodes(m_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < m_mesh.NN(); ++idxn)
        vec_nodes[idxn] = m_mesh.NodeRadius(idx, idxn);
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

      for (int idxq = 0; idxq < nq; ++idxq) {
        auto w_val = pleg(idxq, rad_source) / rad_source;
        auto w_deriv = pleg.Derivative(idxq, rad_source);
        auto idx_u = this->ltgR(0, idx, idxq);

        Complex tmp_pm = 2.0 * cmt.MCmp() * y0c;
        Complex tmp_u = -w_val * tmp_pm + cmt.MC00() * w_deriv * y0c;

        vec_lforce(idx_u, 0) = tmp_u;
      };
    };
  };

  vec_lforce *= (1.0 / m_momentNorm);
  return vec_lforce;
};

}   // namespace Full1D

#endif   // SEM_FORCE_RADIAL_H
