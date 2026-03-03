#ifndef PARAMRED_INFO_H
#define PARAMRED_INFO_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <Eigen/Core>
#include <GSHTrans/Core>
#include <Interpolation/Lagrange>
#include "SourceInfo.h"
#include <EarthMesh/All>
#include "input_parser.h"
#include "mesh_model.h"

class ParamRedInfo {

private:
  using COMPLEX = std::complex<double>;
  using MATRIX = Eigen::MatrixXcd;
  using WIGNER =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Multiple, GSHTrans::ColumnMajor>;
  WIGNER wigdmat;
  std::vector<double> vec_theta_r, vec_phi_r;
  int nrec, _mmax;
  double deg2rad = EIGEN_PI / 180.0;
  std::complex<double> i1 = std::complex<double>(0.0, 1.0);

  COMPLEX ylmn_trial(int l, int m, int N, int idxr) {
    auto dl = this->wigdmat[N, idxr];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(COMPLEX(0.0, m * this->vec_phi_r[idxr]));
    return ylm;
  };

public:
  ParamRedInfo(InputParameters &param, int lmax) {
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
    _mmax = std::min(lmax, 2);

    // wigner class
    wigdmat = WIGNER(lmax, _mmax, 1, vec_theta_r);
  };

  MATRIX RV_RED_SPH(int idxl) {
    // ylmn function for given l, m, N and receiver index
    auto ylmn = [&](int l, int m, int N, int idxr) {
      auto dl = wigdmat[N, idxr];
      auto tmp = dl[l, m];
      auto ylm = tmp * std::exp(std::complex<double>(0.0, m * vec_phi_r[idxr]));
      return ylm;
    };

    // maximum m value for given l
    int max_m = std::min(idxl, _mmax);
    int nm = 2 * max_m + 1;

    // receiver vector to be filled
    MATRIX vec_receiver = MATRIX::Zero(3 * nrec, nm);

    // loop through receivers
    for (int idxr = 0; idxr < nrec; ++idxr) {
      int idx_col = 0;
      for (int idxm = -max_m; idxm < max_m + 1; ++idxm) {
        // spherical harmonic values
        auto yl0 = ylmn_trial(idxl, idxm, 0, idxr);
        auto ylm = ylmn(idxl, idxm, -1, idxr);
        auto ylp = ylmn(idxl, idxm, 1, idxr);

        // adding to receiver vector
        vec_receiver(3 * idxr, idx_col) = yl0;
        vec_receiver(3 * idxr + 1, idx_col) = ylp - ylm;
        vec_receiver(3 * idxr + 2, idx_col) = -i1 * (ylm + ylp);

        // incrementing column index
        ++idx_col;
      }
    }

    return vec_receiver;
  };

  MATRIX RV_RED_TOR(int idxl) {
    auto ylmn = [&](int l, int m, int N, int idxr) {
      auto dl = wigdmat[N, idxr];
      auto tmp = dl[l, m];
      auto ylm = tmp * std::exp(std::complex<double>(0.0, m * vec_phi_r[idxr]));
      return ylm;
    };

    // maximum m value for given l
    int max_m = std::min(idxl, _mmax);
    int nm = 2 * max_m + 1;

    // receiver vector to be filled
    MATRIX vec_receiver = MATRIX::Zero(3 * nrec, nm);

    // loop through receivers
    for (int idxr = 0; idxr < nrec; ++idxr) {
      int idx_col = 0;
      for (int idxm = -max_m; idxm < max_m + 1; ++idxm) {
        // spherical harmonic values
        auto ylm = ylmn(idxl, idxm, -1, idxr);
        auto ylp = ylmn(idxl, idxm, 1, idxr);

        // adding to receiver vector
        vec_receiver(3 * idxr + 1, idx_col) = i1 * (ylm + ylp);
        vec_receiver(3 * idxr + 2, idx_col) = ylp - ylm;

        // incrementing column index
        ++idx_col;
      }
    }

    return vec_receiver;
  };
};

#endif   // PARAM_INFO_H