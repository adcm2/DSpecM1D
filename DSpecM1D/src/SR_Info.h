#ifndef SR_INFO_H
#define SR_INFO_H

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

class SRInfo {

private:
  // using namespace SourceInfo;
  using COMPLEX = std::complex<double>;
  using MATRIX = Eigen::MatrixXcd;
  using WIGNER =
      GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                       GSHTrans::Multiple, GSHTrans::ColumnMajor>;
  WIGNER wigdmat;
  int nrec, _mmax, _lmax;
  double deg2rad = EIGEN_PI / 180.0;
  COMPLEX i1 = std::complex<double>(0.0, 1.0);
  std::vector<double> ep_angles, vec_azimuths, vec_backazimuths;
  std::vector<std::vector<COMPLEX>> vec_az_exp;

  COMPLEX ylmn(int l, int m, int N, int idxr) {
    auto tmp = this->wigdmat[N, idxr][l, m];
    auto ylm = tmp * vec_az_exp[idxr][m + _mmax];
    return ylm;
  };

public:
  SRInfo(InputParameters &);

  MATRIX RV_RED_SPH(int);
  MATRIX RV_RED_TOR(int);
  // MATRIX RV_RED_RAD(int);

  auto epicentralAngles() const { return ep_angles; };
  auto epicentralAngles(int idx) const { return ep_angles[idx]; };

  auto azimuths() const { return vec_azimuths; };
  auto azimuths(int idx) const { return vec_azimuths[idx]; };

  auto backazimuths() const { return vec_backazimuths; };
  auto backazimuths(int idx) const { return vec_backazimuths[idx]; };
};

inline SRInfo::SRInfo(InputParameters &param)
    : _lmax{param.lmax()}, nrec{param.num_receivers()} {
  using namespace GSHTrans;

  // finding epicentral angle, azimuth, backazimuth
  auto twopi = 2.0 * EIGEN_PI;
  auto deg2rad = EIGEN_PI / 180.0;
  auto theta_s = (90.0 - param.source_lat_deg()) * deg2rad;
  auto phi_s = param.source_lon_deg() * deg2rad;
  auto cps = std::cos(phi_s);
  auto cts = std::cos(theta_s);
  auto sps = std::sin(phi_s);
  auto sts = std::sin(theta_s);

  // get the azimuth and takeoff angles for the receiver
  for (int idx = 0; idx < param.num_receivers(); ++idx) {
    auto rec = param.receivers()[idx];
    auto theta_r = (90.0 - rec.first) * deg2rad;
    auto phi_r = rec.second * deg2rad;
    auto cpr = std::cos(phi_r);
    auto ctr = std::cos(theta_r);
    auto spr = std::sin(phi_r);
    auto str = std::sin(theta_r);

    // get epicentral angle
    auto cosep = cts * ctr + sts * str * std::cos(phi_r - phi_s);
    (std::abs(cosep) > 1.0) ? cosep = std::copysign(1.0, cosep) : cosep = cosep;
    auto ep = std::acos(cosep);
    auto sep = std::sin(ep);

    // get backazimuth
    auto cbaz = (str * cts - ctr * sts * std::cos(phi_r - phi_s)) / sep;
    (std::abs(cbaz) > 1.0) ? cbaz = std::copysign(1.0, cbaz) : cbaz = cbaz;
    auto baz = std::acos(cbaz);

    // get azimuth
    auto caz = (sts * ctr - cts * str * std::cos(phi_r - phi_s)) / sep;
    (std::abs(caz) > 1.0) ? caz = std::copysign(1.0, caz) : caz = caz;
    auto az = std::acos(caz);

    // correction for quadrant
    auto dif = std::sin(phi_s - phi_r);
    (dif < 0.0) ? baz = twopi - baz : baz = baz;
    (dif > 0.0) ? az = twopi - az : az = az;

    // final change?
    baz = 2.0 * EIGEN_PI - baz;
    az = EIGEN_PI - az;

    // push backs
    ep_angles.push_back(ep);
    vec_backazimuths.push_back(baz);
    vec_azimuths.push_back(az);
  }

  // declare the Wigner class
  _mmax = std::min(_lmax, 2);
  wigdmat = WIGNER(_lmax, _mmax, 1, ep_angles);

  // define exponential of m of azimuths
  for (int idxr = 0; idxr < nrec; ++idxr) {
    std::vector<COMPLEX> vec_exp_m;
    auto az = vec_azimuths[idxr];
    for (int m = -_mmax; m < _mmax + 1; ++m) {
      vec_exp_m.push_back(std::exp(COMPLEX(0.0, m * az)));
    }
    vec_az_exp.push_back(vec_exp_m);
  }
};

inline Eigen::MatrixXcd
SRInfo::RV_RED_SPH(int idxl) {
  // maximum m value for given l
  int max_m = std::min(idxl, 2);
  int nm = 2 * max_m + 1;

  // receiver vector to be filled
  MATRIX vec_receiver = MATRIX::Zero(3 * nrec, nm);

  // loop through receivers
  for (int idxr = 0; idxr < nrec; ++idxr) {
    int idx_col = 0;
    for (int idxm = -max_m; idxm < max_m + 1; ++idxm) {
      // spherical harmonic values
      auto yl0 = ylmn(idxl, idxm, 0, idxr);
      auto ylm = ylmn(idxl, idxm, -1, idxr);
      auto ylp = ylmn(idxl, idxm, 1, idxr);

      // adding to receiver vector
      vec_receiver(3 * idxr, idx_col) = yl0;
      vec_receiver(3 * idxr + 1, idx_col) = ylm - ylp;
      vec_receiver(3 * idxr + 2, idx_col) = -i1 * (ylm + ylp);

      // incrementing column index
      ++idx_col;
    }
  }

  return vec_receiver;
};

inline Eigen::MatrixXcd
SRInfo::RV_RED_TOR(int idxl) {
  // maximum m value for given l
  int max_m = std::min(idxl, 2);
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
      vec_receiver(3 * idxr + 1, idx_col) = -i1 * (ylm + ylp);
      vec_receiver(3 * idxr + 2, idx_col) = ylp - ylm;

      // incrementing column index
      ++idx_col;
    }
  }

  return vec_receiver;
};
#endif   // SR_INFO_H
