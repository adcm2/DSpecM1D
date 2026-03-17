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
    m_numReceivers = param.num_receivers();
    m_thetaR = std::vector<double>(m_numReceivers, 0.0);
    m_phiR = std::vector<double>(m_numReceivers, 0.0);
    for (int idx = 0; idx < m_numReceivers; ++idx) {
      auto rec = param.receivers()[idx];
      m_thetaR[idx] = (90.0 - rec.first) * m_deg2rad;
      m_phiR[idx] = rec.second * m_deg2rad;
    }

    // wigner class
    m_wigner = Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                      GSHTrans::Multiple, GSHTrans::ColumnMajor>(lmax, lmax, 1,
                                                                 m_thetaR);
  };

  Eigen::MatrixXcd rvFullSph(int idxl) {
    auto ylmn = [&](int l, int m, int N, int idxr) {
      auto dl = m_wigner[N, idxr];
      auto tmp = dl[l, m];
      auto ylm = tmp * std::exp(std::complex<double>(0.0, m * m_phiR[idxr]));
      return ylm;
    };
    Eigen::MatrixXcd vec_receiver =
        Eigen::MatrixXcd::Zero(3 * m_numReceivers, 2 * idxl + 1);

    // loop through receivers
    for (int idxr = 0; idxr < m_numReceivers; ++idxr) {
      for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
        // spherical harmonic values
        auto yl0 = ylmn(idxl, idxm, 0, idxr);
        auto ylm = ylmn(idxl, idxm, -1, idxr);
        auto ylp = ylmn(idxl, idxm, 1, idxr);

        // adding to receiver vector
        vec_receiver(3 * idxr, idxm + idxl) = yl0;
        vec_receiver(3 * idxr + 1, idxm + idxl) = ylp - ylm;
        vec_receiver(3 * idxr + 2, idxm + idxl) = -m_i1 * (ylm + ylp);
      }
    }

    return vec_receiver;
  };

  Eigen::MatrixXcd rvFullTor(int idxl) {
    auto ylmn = [&](int l, int m, int N, int idxr) {
      auto dl = m_wigner[N, idxr];
      auto tmp = dl[l, m];
      auto ylm = tmp * std::exp(std::complex<double>(0.0, m * m_phiR[idxr]));
      return ylm;
    };
    Eigen::MatrixXcd vec_receiver =
        Eigen::MatrixXcd::Zero(3 * m_numReceivers, 2 * idxl + 1);

    // loop through receivers
    for (int idxr = 0; idxr < m_numReceivers; ++idxr) {
      for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
        // spherical harmonic values
        auto ylm = ylmn(idxl, idxm, -1, idxr);
        auto ylp = ylmn(idxl, idxm, 1, idxr);

        // adding to receiver vector
        vec_receiver(3 * idxr + 1, idxm + idxl) = m_i1 * (ylm + ylp);
        vec_receiver(3 * idxr + 2, idxm + idxl) = ylp - ylm;
      }
    }

    return vec_receiver;
  };

private:
  GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                   GSHTrans::Multiple, GSHTrans::ColumnMajor>
      m_wigner;
  std::vector<double> m_thetaR, m_phiR;
  int m_numReceivers;
  double m_deg2rad = EIGEN_PI / 180.0;
  std::complex<double> m_i1 = std::complex<double>(0.0, 1.0);
};

#endif   // PARAM_INFO_H