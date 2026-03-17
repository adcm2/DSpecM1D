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
#include "InputParser.h"
#include "MeshModel.h"

class ParamRedInfo {

private:
  using Complex = std::complex<double>;
  using MatrixC = Eigen::MatrixXcd;
  using WignerType =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Multiple, GSHTrans::ColumnMajor>;
  WignerType m_wigner;
  std::vector<double> m_thetaR, m_phiR;
  int m_numReceivers;
  int m_mmax;
  double m_deg2rad = EIGEN_PI / 180.0;
  Complex m_i1 = Complex(0.0, 1.0);

  Complex ylmnTrial(int l, int m, int N, int idxr) {
    auto dl = this->m_wigner[N, idxr];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(Complex(0.0, m * this->m_phiR[idxr]));
    return ylm;
  };

public:
  ParamRedInfo(InputParameters &param, int lmax) {
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
    m_mmax = std::min(lmax, 2);

    // wigner class
    m_wigner = WignerType(lmax, m_mmax, 1, m_thetaR);
  };

  MatrixC rvRedSph(int idxl) {
    // ylmn function for given l, m, N and receiver index
    auto ylmn = [&](int l, int m, int N, int idxr) {
      auto dl = m_wigner[N, idxr];
      auto tmp = dl[l, m];
      auto ylm = tmp * std::exp(std::complex<double>(0.0, m * m_phiR[idxr]));
      return ylm;
    };

    // maximum m value for given l
    int max_m = std::min(idxl, m_mmax);
    int nm = 2 * max_m + 1;

    // receiver vector to be filled
    MatrixC vec_receiver = MatrixC::Zero(3 * m_numReceivers, nm);

    // loop through receivers
    for (int idxr = 0; idxr < m_numReceivers; ++idxr) {
      int idx_col = 0;
      for (int idxm = -max_m; idxm < max_m + 1; ++idxm) {
        // spherical harmonic values
        auto yl0 = ylmnTrial(idxl, idxm, 0, idxr);
        auto ylm = ylmn(idxl, idxm, -1, idxr);
        auto ylp = ylmn(idxl, idxm, 1, idxr);

        // adding to receiver vector
        vec_receiver(3 * idxr, idx_col) = yl0;
        vec_receiver(3 * idxr + 1, idx_col) = ylp - ylm;
        vec_receiver(3 * idxr + 2, idx_col) = -m_i1 * (ylm + ylp);

        // incrementing column index
        ++idx_col;
      }
    }

    return vec_receiver;
  };

  MatrixC rvRedTor(int idxl) {
    auto ylmn = [&](int l, int m, int N, int idxr) {
      auto dl = m_wigner[N, idxr];
      auto tmp = dl[l, m];
      auto ylm = tmp * std::exp(std::complex<double>(0.0, m * m_phiR[idxr]));
      return ylm;
    };

    // maximum m value for given l
    int max_m = std::min(idxl, m_mmax);
    int nm = 2 * max_m + 1;

    // receiver vector to be filled
    MatrixC vec_receiver = MatrixC::Zero(3 * m_numReceivers, nm);

    // loop through receivers
    for (int idxr = 0; idxr < m_numReceivers; ++idxr) {
      int idx_col = 0;
      for (int idxm = -max_m; idxm < max_m + 1; ++idxm) {
        // spherical harmonic values
        auto ylm = ylmn(idxl, idxm, -1, idxr);
        auto ylp = ylmn(idxl, idxm, 1, idxr);

        // adding to receiver vector
        vec_receiver(3 * idxr + 1, idx_col) = m_i1 * (ylm + ylp);
        vec_receiver(3 * idxr + 2, idx_col) = ylp - ylm;

        // incrementing column index
        ++idx_col;
      }
    }

    return vec_receiver;
  };
};

#endif   // PARAMRED_INFO_H