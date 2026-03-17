#ifndef SEM_FORCE_TOROIDAL_H
#define SEM_FORCE_TOROIDAL_H

#include "SEM.h"

namespace Full1D {

Eigen::MatrixXcd
SEM::calculateForceT(SourceInfo::EarthquakeCMT &cmt, int idxl) {
  int NQ = _mesh.NN();
  totlen = this->ltgT(_eu - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totlen, 2 * idxl + 1);
  double kval =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));

  double rad_source = _SourceRadius(cmt);
  double theta_s = (90.0 - cmt.Latitude()) * EIGEN_PI / (180.0);
  double phi_s = cmt.Longitude() * EIGEN_PI / (180.0);

  int maxn = 2;
  if (maxn > _lmax)
    maxn = _lmax;
  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(_lmax, _lmax,
                                                                maxn, theta_s);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    return tmp * std::exp(Complex(0.0, m * phi));
  };

  Complex isq2 = Complex(0.0, 1.0 / std::sqrt(2.0));
  double omegal2 = std::sqrt((idxl + 2) * (idxl - 1) / 2.0);

  for (int idx = _el; idx < _eu; ++idx) {
    if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) > rad_source)) {
      stdvec vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn)
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

      for (int idxq = 0; idxq < NQ; ++idxq) {
        auto w_val = pleg(idxq, rad_source) / rad_source;
        auto w_prefactor = pleg.Derivative(idxq, rad_source) - w_val;
        std::size_t ridx = this->ltgT(idx, idxq);

        for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
          Complex ymc = std::conj(ylmn(idxl, idxm, -1, phi_s));
          Complex ypc = std::conj(ylmn(idxl, idxm, 1, phi_s));
          Complex ymmc = 0.0, yppc = 0.0;
          if (idxl > 1) {
            ymmc = std::conj(ylmn(idxl, idxm, -2, phi_s));
            yppc = std::conj(ylmn(idxl, idxm, 2, phi_s));
          }

          Complex tmp = w_prefactor * (cmt.MC0m() * ymc - cmt.MC0p() * ypc);
          tmp += w_val * omegal2 * (cmt.MCmm() * ymmc - cmt.MCpp() * yppc);
          tmp *= isq2;

          vec_lforce(ridx, idxm + idxl) = tmp;
        };
      };
    };
  };
  vec_lforce *= (1.0 / _moment_norm);
  return vec_lforce;
};

Eigen::MatrixXcd
SEM::calculateForceAllT(SourceInfo::EarthquakeCMT &cmt, int idxl) {
  int NQ = _mesh.NN();
  totlen = this->ltgT(_mesh.NE() - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totlen, 2);
  double kval =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  double rad_source = _SourceRadius(cmt);

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
        auto idx_v = this->ltgT(idx, idxq);
        vec_lforce(idx_v, 0) = kval * w_deriv;
        vec_lforce(idx_v, 1) = kval * w_val;
      };
    };
  };
  return vec_lforce;
};

Eigen::MatrixXcd
SEM::calculateForceCoefficientsT(SourceInfo::EarthquakeCMT &cmt,
                                       int idxl) {
  int NQ = _mesh.NN();
  totlen = this->ltgT(_mesh.NE() - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(2 * idxl + 1, 2);

  double theta_s = (90.0 - cmt.Latitude()) * EIGEN_PI / (180.0);
  double phi_s = cmt.Longitude() * EIGEN_PI / (180.0);

  int maxn = 2;
  if (maxn > idxl)
    maxn = idxl;
  auto wigdmat =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Single, GSHTrans::ColumnMajor>(idxl, idxl,
                                                                maxn, theta_s);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    auto dl = wigdmat[N];
    auto tmp = dl[l, m];
    return tmp * std::exp(Complex(0.0, m * phi));
  };

  Complex isq2 = Complex(0.0, 1.0 / std::sqrt(2.0));
  double omegal2 = std::sqrt((idxl + 2) * (idxl - 1) / 2.0);

  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    Complex ymc = std::conj(ylmn(idxl, idxm, -1, phi_s));
    Complex ypc = std::conj(ylmn(idxl, idxm, 1, phi_s));
    Complex ymmc = 0.0, yppc = 0.0;
    if (idxl > 1) {
      ymmc = std::conj(ylmn(idxl, idxm, -2, phi_s));
      yppc = std::conj(ylmn(idxl, idxm, 2, phi_s));
    }
    auto tmp1 = cmt.MC0m() * ymc - cmt.MC0p() * ypc;
    auto tmp2 = cmt.MCmm() * ymmc - cmt.MCpp() * yppc;
    vec_lforce(idxm + idxl, 0) = isq2 * tmp1;
    vec_lforce(idxm + idxl, 1) = isq2 * (omegal2 * tmp2 - tmp1);
  };

  vec_lforce *= (1.0 / _moment_norm);
  return vec_lforce;
};

Eigen::MatrixXcd
SEM::calculateForceRedCoefficientsT(SourceInfo::EarthquakeCMT &cmt,
                                           int idxl, double az) {
  int maxn = 2;
  (maxn > idxl) ? maxn = idxl : maxn = maxn;
  Eigen::MatrixXcd vec_force = Eigen::MatrixXcd::Zero(2 * maxn + 1, 2);

  double lv = static_cast<double>(idxl);
  auto mfact = std::sqrt((2.0 * lv + 1.0) / (4.0 * EIGEN_PI));
  auto omegal2 = std::sqrt((lv + 2) * (lv - 1) / 2.0);
  auto isq2 = Complex(0.0, 1.0 / std::sqrt(2.0));

  Complex expm2 = 1.0, expm1 = 1.0, expp1 = 1.0, expp2 = 1.0;

  Complex tmp_mm = cmt.MCmm() * expm2;
  Complex tmp_0m = cmt.MC0m() * expm1;
  Complex tmp_0p = cmt.MC0p() * expp1;
  Complex tmp_pp = cmt.MCpp() * expp2;

  if (idxl == 1) {
    vec_force(0, 0) = tmp_0m;
    vec_force(2, 0) = -tmp_0p;
    vec_force(0, 1) = -tmp_0m;
    vec_force(2, 1) = tmp_0p;
  } else if (idxl > 1) {
    vec_force(1, 0) = tmp_0m;
    vec_force(3, 0) = -tmp_0p;
    vec_force(1, 1) = -tmp_0m;
    vec_force(3, 1) = tmp_0p;
    vec_force(0, 1) = omegal2 * tmp_mm;
    vec_force(4, 1) = -omegal2 * tmp_pp;
  }
  vec_force *= (double) mfact * isq2 * (1.0 / _moment_norm);
  return vec_force;
};

}   // namespace Full1D

#endif   // SEM_FORCE_TOROIDAL_H
