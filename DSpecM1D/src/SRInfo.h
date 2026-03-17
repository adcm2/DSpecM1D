#ifndef DSPECM1D_SR_INFO_H
#define DSPECM1D_SR_INFO_H

#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <EarthMesh/All>
#include <GSHTrans/Core>
#include <Interpolation/Lagrange>
#include "SourceInfo.h"
#include "InputParser.h"
#include "MeshModel.h"

class SRInfo {

private:
  using ComplexD = std::complex<double>;
  using MatrixC = Eigen::MatrixXcd;
  using WignerType =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Multiple, GSHTrans::ColumnMajor>;
  WignerType m_wigner;
  int m_numReceivers;
  int m_mmax;
  int m_lmax;
  const double m_deg2rad = EIGEN_PI / 180.0;
  const ComplexD m_i1 = ComplexD(0.0, 1.0);
  std::vector<double> m_epicentralAngles;
  std::vector<double> m_azimuths;
  std::vector<double> m_backazimuths;
  std::vector<std::vector<ComplexD>> m_azimuthExp;

  ComplexD ylmn(int l, int m, int idxN, int idxr) {
    auto tmp = this->m_wigner[idxN, idxr][l, m];
    auto ylm = tmp * m_azimuthExp[idxr][m + m_mmax];
    return ylm;
  };

public:
  SRInfo(InputParameters &);

  MatrixC rvRedSph(int);
  MatrixC rvRedTor(int);
  // MATRIX rvRedRad(int);

  auto epicentralAngles() const { return m_epicentralAngles; };
  auto epicentralAngles(int idx) const { return m_epicentralAngles[idx]; };

  auto azimuths() const { return m_azimuths; };
  auto azimuths(int idx) const { return m_azimuths[idx]; };

  auto backazimuths() const { return m_backazimuths; };
  auto backazimuths(int idx) const { return m_backazimuths[idx]; };
};

inline SRInfo::SRInfo(InputParameters &param)
    : m_numReceivers{param.num_receivers()}, m_lmax{param.lmax()} {
  using namespace GSHTrans;

  // finding epicentral angle, azimuth, backazimuth
  auto twopi = 2.0 * EIGEN_PI;
  auto thetaS = (90.0 - param.source_lat_deg()) * m_deg2rad;
  auto phiS = param.source_lon_deg() * m_deg2rad;
  auto cts = std::cos(thetaS);
  auto sts = std::sin(thetaS);

  // get the azimuth and takeoff angles for the receiver
  for (int idx = 0; idx < param.num_receivers(); ++idx) {
    auto rec = param.receivers()[idx];
    auto thetaR = (90.0 - rec.first) * m_deg2rad;
    auto phiR = rec.second * m_deg2rad;
    auto ctr = std::cos(thetaR);
    auto str = std::sin(thetaR);

    // get epicentral angle
    auto cosep = cts * ctr + sts * str * std::cos(phiR - phiS);
    (std::abs(cosep) > 1.0) ? cosep = std::copysign(1.0, cosep) : cosep = cosep;
    auto ep = std::acos(cosep);
    auto sep = std::sin(ep);

    // get backazimuth
    auto cbaz = (str * cts - ctr * sts * std::cos(phiR - phiS)) / sep;
    (std::abs(cbaz) > 1.0) ? cbaz = std::copysign(1.0, cbaz) : cbaz = cbaz;
    auto baz = std::acos(cbaz);

    // get azimuth
    auto caz = (sts * ctr - cts * str * std::cos(phiR - phiS)) / sep;
    (std::abs(caz) > 1.0) ? caz = std::copysign(1.0, caz) : caz = caz;
    auto az = std::acos(caz);

    // correction for quadrant
    auto dif = std::sin(phiS - phiR);
    (dif < 0.0) ? baz = twopi - baz : baz = baz;
    (dif > 0.0) ? az = twopi - az : az = az;

    // final change?
    baz = 2.0 * EIGEN_PI - baz;
    az = EIGEN_PI - az;

    // push backs
    m_epicentralAngles.push_back(ep);
    m_backazimuths.push_back(baz);
    m_azimuths.push_back(az);
  }

  // declare the Wigner class
  m_mmax = std::min(m_lmax, 2);
  m_wigner = WignerType(m_lmax, m_mmax, 1, m_epicentralAngles);

  // define exponential of m of azimuths
  for (int idxr = 0; idxr < m_numReceivers; ++idxr) {
    std::vector<ComplexD> vec_exp_m;
    auto az = m_azimuths[idxr];
    for (int m = -m_mmax; m < m_mmax + 1; ++m) {
      vec_exp_m.push_back(std::exp(ComplexD(0.0, m * az)));
    }
    m_azimuthExp.push_back(vec_exp_m);
  }
};

inline Eigen::MatrixXcd
SRInfo::rvRedSph(int idxl) {
  // maximum m value for given l
  int maxM = std::min(idxl, 2);
  int nm = 2 * maxM + 1;

  // receiver vector to be filled
  MatrixC vec_receiver = MatrixC::Zero(3 * m_numReceivers, nm);

  // loop through receivers
  for (int idxr = 0; idxr < m_numReceivers; ++idxr) {
    int idxCol = 0;
    for (int idxm = -maxM; idxm < maxM + 1; ++idxm) {
      // spherical harmonic values
      auto yl0 = ylmn(idxl, idxm, 0, idxr);
      auto ylm = ylmn(idxl, idxm, -1, idxr);
      auto ylp = ylmn(idxl, idxm, 1, idxr);

      // adding to receiver vector
      vec_receiver(3 * idxr, idxCol) = yl0;
      vec_receiver(3 * idxr + 1, idxCol) = ylm - ylp;
      vec_receiver(3 * idxr + 2, idxCol) = -m_i1 * (ylm + ylp);

      // incrementing column index
      ++idxCol;
    }
  }

  return vec_receiver;
};

inline Eigen::MatrixXcd
SRInfo::rvRedTor(int idxl) {
  // maximum m value for given l
  int maxM = std::min(idxl, 2);
  int nm = 2 * maxM + 1;

  // receiver vector to be filled
  MatrixC vec_receiver = MatrixC::Zero(3 * m_numReceivers, nm);

  // loop through receivers
  for (int idxr = 0; idxr < m_numReceivers; ++idxr) {
    int idxCol = 0;
    for (int idxm = -maxM; idxm < maxM + 1; ++idxm) {
      // spherical harmonic values
      auto ylm = ylmn(idxl, idxm, -1, idxr);
      auto ylp = ylmn(idxl, idxm, 1, idxr);

      // adding to receiver vector
      vec_receiver(3 * idxr + 1, idxCol) = -m_i1 * (ylm + ylp);
      vec_receiver(3 * idxr + 2, idxCol) = ylp - ylm;

      // incrementing column index
      ++idxCol;
    }
  }

  return vec_receiver;
};
#endif   // DSPECM1D_SR_INFO_H
