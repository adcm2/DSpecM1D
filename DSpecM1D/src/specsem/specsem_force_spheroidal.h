#ifndef SPECSEM_FORCE_SPHEROIDAL_H
#define SPECSEM_FORCE_SPHEROIDAL_H

#include "specsem.h"

namespace Full1D {

Eigen::MatrixXcd
specsem::CalculateForce(SourceInfo::EarthquakeCMT &cmt, int idxl) {
  int NQ = _mesh.NN();
  totlen = this->LtG_S(2, _mesh.NE() - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totlen, 2 * idxl + 1);
  double kval =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  double kd2 = kval / std::sqrt(2.0);

  double rad_source = _mesh.PR() - 1000.0 * cmt.Depth() / _length_norm;
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
  double lprefac = 1.0;

  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) > rad_source)) {
      stdvec vec_nodes(NQ, 0.0);
      for (int idxn = 0; idxn < NQ; ++idxn)
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

      for (int idxq = 0; idxq < NQ; ++idxq) {
        auto w_val = pleg(idxq, rad_source) / rad_source;
        auto w_deriv = pleg.Derivative(idxq, rad_source);
        auto idx_u = this->LtG_S(0, idx, idxq);
        auto idx_v = this->LtG_S(1, idx, idxq);

        for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
          Complex y0c = std::conj(ylmn(idxl, idxm, 0, phi_s));
          Complex ymc = std::conj(ylmn(idxl, idxm, -1, phi_s));
          Complex ypc = std::conj(ylmn(idxl, idxm, 1, phi_s));
          Complex ymmc = 0.0, yppc = 0.0;
          if (idxl > 1) {
            ymmc = std::conj(ylmn(idxl, idxm, -2, phi_s));
            yppc = std::conj(ylmn(idxl, idxm, 2, phi_s));
          }

          Complex tmp_0pm = cmt.MC0p() * ypc + cmt.MC0m() * ymc;
          Complex tmp_pm = 2.0 * cmt.MCmp() * y0c;
          Complex tmp_ppmm = cmt.MCpp() * yppc + cmt.MCmm() * ymmc;

          Complex tmp_u =
              w_val * (kd2 * tmp_0pm - tmp_pm) + cmt.MC00() * w_deriv * y0c;
          Complex tmp_v =
              kd2 * (w_deriv * tmp_0pm +
                     w_val * (kd2 * tmp_pm - tmp_0pm + omegal2 * tmp_ppmm));

          vec_lforce(idx_u, idxm + idxl) = tmp_u * lprefac;
          vec_lforce(idx_v, idxm + idxl) = tmp_v * lprefac;
        };
      };
    };
  };

  vec_lforce *= (1.0 / _moment_norm);
  return vec_lforce;
};

Eigen::MatrixXcd
specsem::CalculateForce_All(SourceInfo::EarthquakeCMT &cmt, int idxl) {
  int NQ = _mesh.NN();
  totlen = this->LtG_S(2, _mesh.NE() - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(totlen, 4);
  double kval =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  double kd2 = kval / std::sqrt(2.0);
  double rad_source = _mesh.PR() - 1000.0 * cmt.Depth() / _length_norm;

  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_source) && (_mesh.EUR(idx) > rad_source)) {
      stdvec vec_nodes(NQ, 0.0);
      for (int idxn = 0; idxn < NQ; ++idxn)
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());

      for (int idxq = 0; idxq < NQ; ++idxq) {
        auto w_val = pleg(idxq, rad_source) / rad_source;
        auto w_deriv = pleg.Derivative(idxq, rad_source);
        auto idx_u = this->LtG_S(0, idx, idxq);
        auto idx_v = this->LtG_S(1, idx, idxq);

        vec_lforce(idx_u, 0) = w_deriv;
        vec_lforce(idx_u, 1) = w_val;
        vec_lforce(idx_v, 2) = kd2 * w_deriv;
        vec_lforce(idx_v, 3) = kd2 * w_val;
      };
    };
  };
  return vec_lforce;
};

Eigen::MatrixXcd
specsem::CalculateForce_Coefficients(SourceInfo::EarthquakeCMT &cmt, int idxl) {
  int NQ = _mesh.NN();
  totlen = this->LtG_S(2, _mesh.NE() - 1, NQ - 1) + 1;
  Eigen::MatrixXcd vec_lforce = Eigen::MatrixXcd::Zero(2 * idxl + 1, 4);
  double kval =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  double kd2 = kval / std::sqrt(2.0);

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

  double omegal2 = std::sqrt((idxl + 2) * (idxl - 1) / 2.0);

  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    Complex y0c = std::conj(ylmn(idxl, idxm, 0, phi_s));
    Complex ymc = std::conj(ylmn(idxl, idxm, -1, phi_s));
    Complex ypc = std::conj(ylmn(idxl, idxm, 1, phi_s));
    Complex ymmc = 0.0, yppc = 0.0;
    if (idxl > 1) {
      ymmc = std::conj(ylmn(idxl, idxm, -2, phi_s));
      yppc = std::conj(ylmn(idxl, idxm, 2, phi_s));
    }

    Complex tmp_0pm = cmt.MC0p() * ypc + cmt.MC0m() * ymc;
    Complex tmp_pm = 2.0 * cmt.MCmp() * y0c;
    Complex tmp_ppmm = cmt.MCpp() * yppc + cmt.MCmm() * ymmc;

    vec_lforce(idxm + idxl, 0) = cmt.MC00() * y0c;
    vec_lforce(idxm + idxl, 1) = kd2 * tmp_0pm - tmp_pm;
    vec_lforce(idxm + idxl, 2) = tmp_0pm;
    vec_lforce(idxm + idxl, 3) = kd2 * tmp_pm - tmp_0pm + omegal2 * tmp_ppmm;
  };

  vec_lforce *= (1.0 / _moment_norm);
  return vec_lforce;
};

Eigen::MatrixXcd
specsem::CalculateForce_RED_Coefficients(SourceInfo::EarthquakeCMT &cmt,
                                         int idxl, double az) {
  int maxn = 2;
  (maxn > idxl) ? maxn = idxl : maxn = maxn;
  Eigen::MatrixXcd vec_force = Eigen::MatrixXcd::Zero(2 * maxn + 1, 4);

  double lv = static_cast<double>(idxl);
  auto mfact = std::sqrt((2.0 * lv + 1.0) / (4.0 * EIGEN_PI));
  auto omegal2 = std::sqrt((lv + 2) * (lv - 1) / 2.0);
  auto kval = std::sqrt(lv * (lv + 1.0));
  auto kd2 = kval / std::sqrt(2.0);
  auto isq2 = Complex(0.0, 1.0 / std::sqrt(2.0));

  Complex expm2 = 1.0, expm1 = 1.0, expp1 = 1.0, expp2 = 1.0;

  Complex tmp_mm = cmt.MCmm() * expm2;
  Complex tmp_0m = cmt.MC0m() * expm1;
  Complex tmp_0p = cmt.MC0p() * expp1;
  Complex tmp_pp = cmt.MCpp() * expp2;

  if (idxl == 1) {
    vec_force(1, 0) = cmt.MC00();
    vec_force(0, 1) = kd2 * tmp_0m;
    vec_force(1, 1) = -2.0 * cmt.MCmp();
    vec_force(2, 1) = kd2 * tmp_0p;
    vec_force(0, 2) = tmp_0m;
    vec_force(2, 2) = tmp_0p;
    vec_force(0, 3) = -tmp_0m;
    vec_force(1, 3) = kd2 * 2.0 * cmt.MCmp();
    vec_force(2, 3) = -tmp_0p;
  } else if (idxl > 1) {
    vec_force(2, 0) = cmt.MC00();
    vec_force(1, 1) = kd2 * tmp_0m;
    vec_force(2, 1) = -2.0 * cmt.MCmp();
    vec_force(3, 1) = kd2 * tmp_0p;
    vec_force(1, 2) = tmp_0m;
    vec_force(3, 2) = tmp_0p;
    vec_force(0, 3) = omegal2 * tmp_mm;
    vec_force(1, 3) = -tmp_0m;
    vec_force(2, 3) = kd2 * 2.0 * cmt.MCmp();
    vec_force(3, 3) = -tmp_0p;
    vec_force(4, 3) = omegal2 * tmp_pp;
  }
  vec_force *= (double) mfact * (1.0 / _moment_norm);
  return vec_force;
};

}   // namespace Full1D

#endif   // SPECSEM_FORCE_SPHEROIDAL_H
