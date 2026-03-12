#ifndef SPECSEM_RECEIVERS_H
#define SPECSEM_RECEIVERS_H

#include "specsem.h"

namespace Full1D {

Eigen::MatrixXcd
specsem::RV_FULL(InputParameters &param, int idxl) {
  auto nrec = param.num_receivers();
  using MATRIX = Eigen::MatrixXcd;
  MATRIX vec_receiver = MATRIX::Zero(3 * nrec, 2 * idxl + 1);
  using namespace GSHTrans;
  int maxn = 1;
  auto i1 = std::complex<double>(0.0, 1.0);
  auto deg2rad = EIGEN_PI / 180.0;
  std::vector<double> vec_theta_r(nrec, 0.0), vec_phi_r(nrec, 0.0);
  for (int idx = 0; idx < nrec; ++idx) {
    auto rec = param.receivers()[idx];
    vec_theta_r[idx] = (90.0 - rec.first) * deg2rad;
    vec_phi_r[idx] = rec.second * deg2rad;
  }
  auto wigdmat_all = Wigner<double, Ortho, All, All, Multiple, ColumnMajor>(
      idxl, idxl, maxn, vec_theta_r);
  auto ylmn_all = [&wigdmat_all, &vec_phi_r](int l, int m, int N, int idxr) {
    auto dl = wigdmat_all[N, idxr];
    auto tmp = dl[l, m];
    return tmp * std::exp(std::complex<double>(0.0, m * vec_phi_r[idxr]));
  };
  for (int idxr = 0; idxr < nrec; ++idxr) {
    for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
      Complex yl0 = ylmn_all(idxl, idxm, 0, idxr);
      Complex ylm = ylmn_all(idxl, idxm, -1, idxr);
      Complex ylp = ylmn_all(idxl, idxm, 1, idxr);
      vec_receiver(3 * idxr, idxm + idxl) = yl0;
      vec_receiver(3 * idxr + 1, idxm + idxl) = ylp - ylm;
      vec_receiver(3 * idxr + 2, idxm + idxl) = -i1 * (ylm + ylp);
    }
  }
  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_FULL_T(InputParameters &param, int idxl) {
  auto nrec = param.num_receivers();
  using MATRIX = Eigen::MatrixXcd;
  MATRIX vec_receiver = MATRIX::Zero(3 * nrec, 2 * idxl + 1);
  using namespace GSHTrans;
  auto i1 = std::complex<double>(0.0, 1.0);
  auto deg2rad = EIGEN_PI / 180.0;
  std::vector<double> vec_theta_r(nrec, 0.0), vec_phi_r(nrec, 0.0);
  for (int idx = 0; idx < nrec; ++idx) {
    auto rec = param.receivers()[idx];
    vec_theta_r[idx] = (90.0 - rec.first) * deg2rad;
    vec_phi_r[idx] = rec.second * deg2rad;
  }
  auto wigdmat_all = Wigner<double, Ortho, All, All, Multiple, ColumnMajor>(
      idxl, idxl, 1, vec_theta_r);
  auto ylmn_all = [&wigdmat_all, &vec_phi_r](int l, int m, int N, int idxr) {
    auto dl = wigdmat_all[N, idxr];
    auto tmp = dl[l, m];
    return tmp * std::exp(std::complex<double>(0.0, m * vec_phi_r[idxr]));
  };
  for (int idxr = 0; idxr < nrec; ++idxr) {
    for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
      Complex ylm = ylmn_all(idxl, idxm, -1, idxr);
      Complex ylp = ylmn_all(idxl, idxm, 1, idxr);
      vec_receiver(3 * idxr + 1, idxm + idxl) = i1 * (ylm + ylp);
      vec_receiver(3 * idxr + 2, idxm + idxl) = ylp - ylm;
    }
  }
  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_BASE_Z(InputParameters &param, int idxl, int idxr) {
  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;
  std::size_t flen = this->LtG_S(2, _mesh.NE() - 1, _mesh.NN() - 1) + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(1, flen);
  bool evaluated = false;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn)
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        auto idx_u = this->LtG_S(0, idx, idxq);
        vec_receiver(0, idx_u) = pleg(idxq, rad_r);
      }
    }
  }
  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_VAL_Z(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(2 * idxl + 1, 1);
  using namespace GSHTrans;
  int maxn = std::min(2, _lmax);
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    return wigdmat[N][l, m] * std::exp(std::complex<double>(0.0, m * phi));
  };
  for (int idxm = -idxl; idxm < idxl + 1; ++idxm)
    vec_receiver(idxm + idxl, 0) = ylmn(idxl, idxm, 0, phi_r);
  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_BASE_THETA(InputParameters &param, int idxl, int idxr) {
  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;
  auto k =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  std::size_t flen = this->LtG_S(2, _mesh.NE() - 1, _mesh.NN() - 1) + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(1, flen);
  bool evaluated = false;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn)
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        auto idx_v = this->LtG_S(1, idx, idxq);
        vec_receiver(0, idx_v) = k / 2.0 * pleg(idxq, rad_r);
      }
    }
  }
  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_VAL_THETA(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(2 * idxl + 1, 1);
  using namespace GSHTrans;
  int maxn = std::min(2, _lmax);
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    return wigdmat[N][l, m] * std::exp(std::complex<double>(0.0, m * phi));
  };
  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    Complex ylm = ylmn(idxl, idxm, -1, phi_r);
    Complex ylp = ylmn(idxl, idxm, 1, phi_r);
    vec_receiver(idxm + idxl, 0) = ylm - ylp;
  }
  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_BASE_THETA_T(InputParameters &param, int idxl, int idxr) {
  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;
  std::size_t flen = this->LtG_T(_mesh.NE() - 1, _mesh.NN() - 1) + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(1, flen);
  bool evaluated = false;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn)
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        auto idx_v = this->LtG_T(idx, idxq);
        vec_receiver(0, idx_v) = 0.5 * pleg(idxq, rad_r);
      }
    }
  }
  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_VAL_THETA_T(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(2 * idxl + 1, 1);
  using namespace GSHTrans;
  int maxn = std::min(2, _lmax);
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    return wigdmat[N][l, m] * std::exp(std::complex<double>(0.0, m * phi));
  };
  auto i1 = std::complex<double>(0.0, 1.0);
  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    Complex ylm = ylmn(idxl, idxm, -1, phi_r);
    Complex ylp = ylmn(idxl, idxm, 1, phi_r);
    vec_receiver(idxm + idxl, 0) = -i1 * (ylm + ylp);
  }
  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_VAL_PHI(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(2 * idxl + 1, 1);
  using namespace GSHTrans;
  int maxn = std::min(2, _lmax);
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    return wigdmat[N][l, m] * std::exp(std::complex<double>(0.0, m * phi));
  };
  Complex i1 = std::complex<double>(0.0, 1.0);
  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    Complex ylm = ylmn(idxl, idxm, -1, phi_r);
    Complex ylp = ylmn(idxl, idxm, 1, phi_r);
    vec_receiver(idxm + idxl, 0) = -i1 * (ylm + ylp);
  }
  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_BASE_PHI_T(InputParameters &param, int idxl, int idxr) {
  return this->RV_BASE_THETA_T(param, idxl, idxr);
};

Eigen::MatrixXcd
specsem::RV_VAL_PHI_T(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(2 * idxl + 1, 1);
  using namespace GSHTrans;
  int maxn = std::min(2, _lmax);
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    return wigdmat[N][l, m] * std::exp(std::complex<double>(0.0, m * phi));
  };
  for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
    Complex ylm = ylmn(idxl, idxm, -1, phi_r);
    Complex ylp = ylmn(idxl, idxm, 1, phi_r);
    vec_receiver(idxm + idxl, 0) = ylp - ylm;
  }
  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_THETA_T(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;
  std::size_t flen = this->LtG_T(_eu - 1, _mesh.NN() - 1) + 1;
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(flen, fcols);
  using namespace GSHTrans;
  int maxn = std::min(2, _lmax);
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    return (std::complex<double>) wigdmat[N][l, m] *
           std::exp(std::complex<double>(0.0, m * phi));
  };
  Complex i1 = std::complex<double>(0.0, 1.0);
  bool evaluated = false;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn)
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        auto idx_v = this->LtG_T(idx, idxq);
        for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
          Complex ylm = ylmn(idxl, idxm, -1, phi_r);
          Complex ylp = ylmn(idxl, idxm, 1, phi_r);
          vec_receiver(idx_v, idxm + idxl) =
              -i1 * pleg(idxq, rad_r) * (ylp + ylm) / 2.0;
        }
      }
    }
  }
  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_PHI_T(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;
  std::size_t flen = this->LtG_T(_eu - 1, _mesh.NN() - 1) + 1;
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(flen, fcols);
  using namespace GSHTrans;
  int maxn = std::min(2, _lmax);
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      _lmax, _lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    return (std::complex<double>) wigdmat[N][l, m] *
           std::exp(std::complex<double>(0.0, m * phi));
  };
  bool evaluated = false;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn)
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        auto idx_v = this->LtG_T(idx, idxq);
        for (int idxm = -idxl; idxm < idxl + 1; ++idxm) {
          Complex ylm = ylmn(idxl, idxm, -1, phi_r);
          Complex ylp = ylmn(idxl, idxm, 1, phi_r);
          vec_receiver(idx_v, idxm + idxl) =
              pleg(idxq, rad_r) * (ylp - ylm) / 2.0;
        }
      }
    }
  }
  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_Z_R(InputParameters &param, int idxr) {
  auto rec = param.receivers()[idxr];
  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;
  std::size_t flen = this->LtG_R(1, _mesh.NE() - 1, _mesh.NN() - 1) + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(flen, 1);
  using namespace GSHTrans;
  auto wigdmat2 =
      Wigner<double, Ortho, All, All, Single, ColumnMajor>(0, 0, 0, 0.0);
  Complex yl0 = wigdmat2[0][0, 0];
  bool evaluated = false;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn)
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        auto idx_u = this->LtG_R(0, idx, idxq);
        vec_receiver(idx_u, 0) = yl0 * pleg(idxq, rad_r);
      }
    }
  }
  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_RED_Z_R(InputParameters &param) {
  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;
  auto rec_elems = this->Receiver_Elements(param);
  auto lowidx = this->LtG_R(0, rec_elems[0], 0);
  auto upidx = this->LtG_R(1, rec_elems.back(), _mesh.NN() - 1);
  int lenidx = upidx - lowidx + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(lenidx, 1);
  using namespace GSHTrans;
  auto wigdmat =
      Wigner<double, Ortho, All, All, Single, ColumnMajor>(0, 0, 0, 0.0);
  Complex yl0 = wigdmat[0][0, 0];
  bool evaluated = false;
  for (int idx = 0; idx < _mesh.NE(); ++idx) {
    if ((_mesh.ELR(idx) <= rad_r) && (_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;
      std::vector<double> vec_nodes(_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn)
        vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
        auto idx_u = this->LtG_R(0, idx, idxq) - lowidx;
        vec_receiver(idx_u, 0) = yl0 * pleg(idxq, rad_r);
      }
    }
  }
  return vec_receiver;
};

Eigen::MatrixXcd
specsem::RV_BASE_FULL(InputParameters &param, int idxl) {
  auto rec_elems = this->Receiver_Elements(param);
  auto lowidx = this->LtG_S(0, rec_elems[0], 0);
  auto upidx = this->LtG_S(1, rec_elems.back(), _NQ - 1);
  int lenidx = upidx - lowidx + 1;
  auto nrec = param.num_receivers();
  Eigen::MatrixXcd mat_base = Eigen::MatrixXcd::Zero(3 * nrec, lenidx);
  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;
  double k = std::sqrt(1.0 * idxl * (idxl + 1.0));
  for (int idx = rec_elems[0]; idx < rec_elems.back() + 1; ++idx) {
    std::vector<double> vec_nodes(_mesh.NN(), 0.0);
    for (int idxn = 0; idxn < _mesh.NN(); ++idxn)
      vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
    auto pleg =
        Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
    for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
      auto idx_u = this->LtG_S(0, idx, idxq) - lowidx;
      auto idx_v = this->LtG_S(1, idx, idxq) - lowidx;
      auto zv = pleg(idxq, rad_r);
      auto tv = k / 2.0 * zv;
      for (int idxr = 0; idxr < nrec; ++idxr) {
        mat_base(3 * idxr, idx_u) = zv;
        mat_base(3 * idxr + 1, idx_v) = tv;
        mat_base(3 * idxr + 2, idx_v) = tv;
      }
    }
  }
  return mat_base;
};

Eigen::MatrixXcd
specsem::RV_BASE_FULL_T(InputParameters &param, int idxl) {
  auto rec_elems = this->Receiver_Elements(param);
  auto lowidx = this->LtG_T(rec_elems[0], 0);
  auto upidx = this->LtG_T(rec_elems.back(), _NQ - 1);
  int lenidx = upidx - lowidx + 1;
  auto nrec = param.num_receivers();
  Eigen::MatrixXcd mat_base = Eigen::MatrixXcd::Zero(3 * nrec, lenidx);
  double rad_r = _mesh.PR() - param.receiver_depth() * 1000.0 / _length_norm;
  double k = std::sqrt(1.0 * idxl * (idxl + 1.0));
  for (int idx = rec_elems[0]; idx < rec_elems.back() + 1; ++idx) {
    std::vector<double> vec_nodes(_mesh.NN(), 0.0);
    for (int idxn = 0; idxn < _mesh.NN(); ++idxn)
      vec_nodes[idxn] = _mesh.NodeRadius(idx, idxn);
    auto pleg =
        Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
    for (int idxq = 0; idxq < _mesh.NN(); ++idxq) {
      auto idx_v = this->LtG_T(idx, idxq) - lowidx;
      auto tv = k / 2.0 * pleg(idxq, rad_r);
      for (int idxr = 0; idxr < nrec; ++idxr) {
        mat_base(3 * idxr + 1, idx_v) = tv;
        mat_base(3 * idxr + 2, idx_v) = tv;
      }
    }
  }
  return mat_base;
};

}   // namespace Full1D

#endif   // SPECSEM_RECEIVERS_H
