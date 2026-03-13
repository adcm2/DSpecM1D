#ifndef PARAM_INFO_H
#define PARAM_INFO_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Core>
#include <GSHTrans/Core>
#include <Interpolation/Lagrange>
#include "SourceInfo.h"
#include <EarthMesh/All>
#include "InputParser.h"
#include "MeshModel.h"

class ParamInfo {
public:
  ParamInfo(InputParameters &param, int lmax) {
    using namespace GSHTrans;
    // vector of receiver angle values:
    nrec = param.num_receivers();
    vec_theta_r = std::vector<double>(nrec, 0.0);
    vec_phi_r = std::vector<double>(nrec, 0.0);
    for (int idx = 0; idx < nrec; ++idx) {
      auto rec = param.receivers()[idx];
      vec_theta_r[idx] = (90.0 - rec.first) * deg2rad;
      vec_phi_r[idx] = rec.second * deg2rad;
    }

    // wigner class
    wigdmat = Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                     GSHTrans::Multiple, GSHTrans::ColumnMajor>(lmax, lmax, 1,
                                                                vec_theta_r);
  };

  Eigen::MatrixXcd RV_FULL_SPH(int idxl) {
    auto ylmn = [&](int l, int m, int N, int idxr) {
      auto dl = wigdmat[N, idxr];
      auto tmp = dl[l, m];
      auto ylm = tmp * std::exp(std::complex<double>(0.0, m * vec_phi_r[idxr]));
      return ylm;
    };
    Eigen::MatrixXcd vec_receiver =
        Eigen::MatrixXcd::Zero(3 * nrec, 2 * idxl + 1);

    // loop through receivers
    for (int idxr = 0; idxr < nrec; ++idxr) {
      for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
        // spherical harmonic values
        auto yl0 = ylmn(idxl, idxm, 0, idxr);
        auto ylm = ylmn(idxl, idxm, -1, idxr);
        auto ylp = ylmn(idxl, idxm, 1, idxr);

        // adding to receiver vector
        vec_receiver(3 * idxr, idxm + idxl) = yl0;
        vec_receiver(3 * idxr + 1, idxm + idxl) = ylp - ylm;
        vec_receiver(3 * idxr + 2, idxm + idxl) = -i1 * (ylm + ylp);
      }
    }

    return vec_receiver;
  };

  Eigen::MatrixXcd RV_FULL_TOR(int idxl) {
    auto ylmn = [&](int l, int m, int N, int idxr) {
      auto dl = wigdmat[N, idxr];
      auto tmp = dl[l, m];
      auto ylm = tmp * std::exp(std::complex<double>(0.0, m * vec_phi_r[idxr]));
      return ylm;
    };
    Eigen::MatrixXcd vec_receiver =
        Eigen::MatrixXcd::Zero(3 * nrec, 2 * idxl + 1);

    // loop through receivers
    for (int idxr = 0; idxr < nrec; ++idxr) {
      for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
        // spherical harmonic values
        auto ylm = ylmn(idxl, idxm, -1, idxr);
        auto ylp = ylmn(idxl, idxm, 1, idxr);

        // adding to receiver vector
        vec_receiver(3 * idxr + 1, idxm + idxl) = i1 * (ylm + ylp);
        vec_receiver(3 * idxr + 2, idxm + idxl) = ylp - ylm;
      }
    }

    return vec_receiver;
  };

private:
  GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                   GSHTrans::Multiple, GSHTrans::ColumnMajor>
      wigdmat;
  std::vector<double> vec_theta_r, vec_phi_r;
  int nrec;
  double deg2rad = EIGEN_PI / 180.0;
  std::complex<double> i1 = std::complex<double>(0.0, 1.0);
};

#endif   // PARAM_INFO_H