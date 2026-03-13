#ifndef SPECSEM_FORCE_RADIAL_H
#define SPECSEM_FORCE_RADIAL_H

#include "specsem.h"

namespace Full1D {

Eigen::MatrixXcd
specsem::CalculateForce_R(SourceInfo::EarthquakeCMT &cmt) {
  int NQ = _mesh.NN();
  totlen = this->LtG_R(1, _mesh.NE() - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totlen, 1);

  double rad_source = _SourceRadius(cmt);

  auto wigdmat2 =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(0, 0, 0, 0.0);
  Complex y0c = wigdmat2[0][0, 0];

  double lprefac = 1.0;

  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) > rad_source)) {
      stdvec vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn)
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

      for (int idxq = 0; idxq < NQ; ++idxq) {
        auto w_val = pleg(idxq, rad_source) / rad_source;
        auto w_deriv = pleg.Derivative(idxq, rad_source);
        auto idx_u = this->LtG_R(0, idx, idxq);

        Complex tmp_pm = 2.0 * cmt.MCmp() * y0c;
        Complex tmp_u = -w_val * tmp_pm + cmt.MC00() * w_deriv * y0c;

        vec_lforce(idx_u, 0) = tmp_u * lprefac;
      };
    };
  };

  vec_lforce *= (1.0 / _moment_norm);
  return vec_lforce;
};

Eigen::MatrixXcd
specsem::CalculateForce_Red_R(SourceInfo::EarthquakeCMT &cmt) {
  int NQ = _mesh.NN();
  totlen = this->LtG_R(1, _mesh.NE() - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totlen, 1);

  double rad_source = _SourceRadius(cmt);

  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(0, 0, 0, 0.0);
  Complex y0c = wigdmat[0][0, 0];

  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) > rad_source)) {
      stdvec vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn)
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

      for (int idxq = 0; idxq < NQ; ++idxq) {
        auto w_val = pleg(idxq, rad_source) / rad_source;
        auto w_deriv = pleg.Derivative(idxq, rad_source);
        auto idx_u = this->LtG_R(0, idx, idxq);

        Complex tmp_pm = 2.0 * cmt.MCmp() * y0c;
        Complex tmp_u = -w_val * tmp_pm + cmt.MC00() * w_deriv * y0c;

        vec_lforce(idx_u, 0) = tmp_u;
      };
    };
  };

  vec_lforce *= (1.0 / _moment_norm);
  return vec_lforce;
};

}   // namespace Full1D

#endif   // SPECSEM_FORCE_RADIAL_H
