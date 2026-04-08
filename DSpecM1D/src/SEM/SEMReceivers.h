#ifndef SEM_RECEIVERS_H
#define SEM_RECEIVERS_H

#include "SEM.h"

namespace Full1D {

Eigen::MatrixXcd
SEM::rvFull(InputParameters &param, int idxl) {
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
SEM::rvFullT(InputParameters &param, int idxl) {
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
SEM::rvBaseZ(InputParameters &param, int idxl, int idxr) {
  double rad_r = _ReceiverRadius(param);
  std::size_t flen = this->ltgS(2, m_mesh.NE() - 1, m_mesh.NN() - 1) + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(1, flen);
  bool evaluated = false;
  for (int idx = 0; idx < m_mesh.NE(); ++idx) {
    if ((m_mesh.ELR(idx) <= rad_r) && (m_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;
      std::vector<double> vec_nodes(m_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < m_mesh.NN(); ++idxn)
        vec_nodes[idxn] = m_mesh.NodeRadius(idx, idxn);
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      for (int idxq = 0; idxq < m_mesh.NN(); ++idxq) {
        auto idx_u = this->ltgS(0, idx, idxq);
        vec_receiver(0, idx_u) = pleg(idxq, rad_r);
      }
    }
  }
  return vec_receiver;
};

Eigen::MatrixXcd
SEM::rvValZ(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(2 * idxl + 1, 1);
  using namespace GSHTrans;
  int maxn = std::min(2, m_lmax);
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      m_lmax, m_lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    return wigdmat[N][l, m] * std::exp(std::complex<double>(0.0, m * phi));
  };
  for (int idxm = -idxl; idxm < idxl + 1; ++idxm)
    vec_receiver(idxm + idxl, 0) = ylmn(idxl, idxm, 0, phi_r);
  return vec_receiver;
};

Eigen::MatrixXcd
SEM::rvBaseTheta(InputParameters &param, int idxl, int idxr) {
  double rad_r = _ReceiverRadius(param);
  auto k =
      std::sqrt(static_cast<double>(idxl) * (static_cast<double>(idxl) + 1.0));
  std::size_t flen = this->ltgS(2, m_mesh.NE() - 1, m_mesh.NN() - 1) + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(1, flen);
  bool evaluated = false;
  for (int idx = 0; idx < m_mesh.NE(); ++idx) {
    if ((m_mesh.ELR(idx) <= rad_r) && (m_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;
      std::vector<double> vec_nodes(m_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < m_mesh.NN(); ++idxn)
        vec_nodes[idxn] = m_mesh.NodeRadius(idx, idxn);
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      for (int idxq = 0; idxq < m_mesh.NN(); ++idxq) {
        auto idx_v = this->ltgS(1, idx, idxq);
        vec_receiver(0, idx_v) = k / 2.0 * pleg(idxq, rad_r);
      }
    }
  }
  return vec_receiver;
};

Eigen::MatrixXcd
SEM::rvValTheta(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(2 * idxl + 1, 1);
  using namespace GSHTrans;
  int maxn = std::min(2, m_lmax);
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      m_lmax, m_lmax, maxn, theta_r);
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
SEM::rvBaseThetaT(InputParameters &param, int idxl, int idxr) {
  double rad_r = _ReceiverRadius(param);
  std::size_t flen = this->ltgT(m_mesh.NE() - 1, m_mesh.NN() - 1) + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(1, flen);
  bool evaluated = false;
  for (int idx = 0; idx < m_mesh.NE(); ++idx) {
    if ((m_mesh.ELR(idx) <= rad_r) && (m_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;
      std::vector<double> vec_nodes(m_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < m_mesh.NN(); ++idxn)
        vec_nodes[idxn] = m_mesh.NodeRadius(idx, idxn);
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      for (int idxq = 0; idxq < m_mesh.NN(); ++idxq) {
        auto idx_v = this->ltgT(idx, idxq);
        vec_receiver(0, idx_v) = 0.5 * pleg(idxq, rad_r);
      }
    }
  }
  return vec_receiver;
};

Eigen::MatrixXcd
SEM::rvValThetaT(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(2 * idxl + 1, 1);
  using namespace GSHTrans;
  int maxn = std::min(2, m_lmax);
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      m_lmax, m_lmax, maxn, theta_r);
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
SEM::rvValPhi(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(2 * idxl + 1, 1);
  using namespace GSHTrans;
  int maxn = std::min(2, m_lmax);
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      m_lmax, m_lmax, maxn, theta_r);
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
SEM::rvBasePhiT(InputParameters &param, int idxl, int idxr) {
  return this->rvBaseThetaT(param, idxl, idxr);
};

Eigen::MatrixXcd
SEM::rvValPhiT(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(2 * idxl + 1, 1);
  using namespace GSHTrans;
  int maxn = std::min(2, m_lmax);
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      m_lmax, m_lmax, maxn, theta_r);
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
SEM::rvThetaT(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  double rad_r = _ReceiverRadius(param);
  std::size_t flen = this->ltgT(m_eu - 1, m_mesh.NN() - 1) + 1;
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(flen, fcols);
  using namespace GSHTrans;
  int maxn = std::min(2, m_lmax);
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      m_lmax, m_lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    return (std::complex<double>) wigdmat[N][l, m] *
           std::exp(std::complex<double>(0.0, m * phi));
  };
  Complex i1 = std::complex<double>(0.0, 1.0);
  bool evaluated = false;
  for (int idx = 0; idx < m_mesh.NE(); ++idx) {
    if ((m_mesh.ELR(idx) <= rad_r) && (m_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;
      std::vector<double> vec_nodes(m_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < m_mesh.NN(); ++idxn)
        vec_nodes[idxn] = m_mesh.NodeRadius(idx, idxn);
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      for (int idxq = 0; idxq < m_mesh.NN(); ++idxq) {
        auto idx_v = this->ltgT(idx, idxq);
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
SEM::rvPhiT(InputParameters &param, int idxl, int idxr) {
  auto rec = param.receivers()[idxr];
  double theta_r = (90.0 - rec.first) * EIGEN_PI / 180.0;
  double phi_r = rec.second * EIGEN_PI / 180.0;
  double rad_r = _ReceiverRadius(param);
  std::size_t flen = this->ltgT(m_eu - 1, m_mesh.NN() - 1) + 1;
  std::size_t fcols = 2 * idxl + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(flen, fcols);
  using namespace GSHTrans;
  int maxn = std::min(2, m_lmax);
  auto wigdmat = Wigner<double, Ortho, All, All, Single, ColumnMajor>(
      m_lmax, m_lmax, maxn, theta_r);
  auto ylmn = [&wigdmat](int l, int m, int N, double phi) {
    return (std::complex<double>) wigdmat[N][l, m] *
           std::exp(std::complex<double>(0.0, m * phi));
  };
  bool evaluated = false;
  for (int idx = 0; idx < m_mesh.NE(); ++idx) {
    if ((m_mesh.ELR(idx) <= rad_r) && (m_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;
      std::vector<double> vec_nodes(m_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < m_mesh.NN(); ++idxn)
        vec_nodes[idxn] = m_mesh.NodeRadius(idx, idxn);
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      for (int idxq = 0; idxq < m_mesh.NN(); ++idxq) {
        auto idx_v = this->ltgT(idx, idxq);
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
SEM::rvZR(InputParameters &param, int idxr) {
  auto rec = param.receivers()[idxr];
  double rad_r = _ReceiverRadius(param);
  std::size_t flen = this->ltgR(1, m_mesh.NE() - 1, m_mesh.NN() - 1) + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(flen, 1);
  using namespace GSHTrans;
  auto wigdmat2 =
      Wigner<double, Ortho, All, All, Single, ColumnMajor>(0, 0, 0, 0.0);
  Complex yl0 = wigdmat2[0][0, 0];
  bool evaluated = false;
  for (int idx = 0; idx < m_mesh.NE(); ++idx) {
    if ((m_mesh.ELR(idx) <= rad_r) && (m_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;
      std::vector<double> vec_nodes(m_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < m_mesh.NN(); ++idxn)
        vec_nodes[idxn] = m_mesh.NodeRadius(idx, idxn);
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      for (int idxq = 0; idxq < m_mesh.NN(); ++idxq) {
        auto idx_u = this->ltgR(0, idx, idxq);
        vec_receiver(idx_u, 0) = yl0 * pleg(idxq, rad_r);
      }
    }
  }
  return vec_receiver;
};

Eigen::MatrixXcd
SEM::rvRedZR(InputParameters &param) {
  double rad_r = _ReceiverRadius(param);
  auto rec_elems = this->receiverElements(param);
  auto lowidx = this->ltgR(0, rec_elems[0], 0);
  auto upidx = this->ltgR(1, rec_elems.back(), m_mesh.NN() - 1);
  int lenidx = upidx - lowidx + 1;
  Eigen::MatrixXcd vec_receiver = Eigen::MatrixXcd::Zero(lenidx, 1);
  using namespace GSHTrans;
  auto wigdmat =
      Wigner<double, Ortho, All, All, Single, ColumnMajor>(0, 0, 0, 0.0);
  Complex yl0 = wigdmat[0][0, 0];
  bool evaluated = false;
  for (int idx = 0; idx < m_mesh.NE(); ++idx) {
    if ((m_mesh.ELR(idx) <= rad_r) && (m_mesh.EUR(idx) >= rad_r) &&
        (!evaluated)) {
      evaluated = true;
      std::vector<double> vec_nodes(m_mesh.NN(), 0.0);
      for (int idxn = 0; idxn < m_mesh.NN(); ++idxn)
        vec_nodes[idxn] = m_mesh.NodeRadius(idx, idxn);
      auto pleg =
          Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
      for (int idxq = 0; idxq < m_mesh.NN(); ++idxq) {
        auto idx_u = this->ltgR(0, idx, idxq) - lowidx;
        vec_receiver(idx_u, 0) = yl0 * pleg(idxq, rad_r);
      }
    }
  }
  return vec_receiver;
};

Eigen::MatrixXcd
SEM::rvBaseFull(InputParameters &param, int idxl) {
  auto rec_elems = this->receiverElements(param);
  auto lowidx = this->ltgS(0, rec_elems[0], 0);
  auto upidx = this->ltgS(1, rec_elems.back(), m_nq - 1);
  int lenidx = upidx - lowidx + 1;
  auto nrec = param.num_receivers();
  Eigen::MatrixXcd mat_base = Eigen::MatrixXcd::Zero(3 * nrec, lenidx);
  double rad_r = _ReceiverRadius(param);
  double k = std::sqrt(1.0 * idxl * (idxl + 1.0));
  for (int idx = rec_elems[0]; idx < rec_elems.back() + 1; ++idx) {
    std::vector<double> vec_nodes(m_mesh.NN(), 0.0);
    for (int idxn = 0; idxn < m_mesh.NN(); ++idxn)
      vec_nodes[idxn] = m_mesh.NodeRadius(idx, idxn);
    auto pleg =
        Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
    for (int idxq = 0; idxq < m_mesh.NN(); ++idxq) {
      auto idx_u = this->ltgS(0, idx, idxq) - lowidx;
      auto idx_v = this->ltgS(1, idx, idxq) - lowidx;
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
SEM::rvBaseFullT(InputParameters &param, int idxl) {
  auto rec_elems = this->receiverElements(param);
  auto lowidx = this->ltgT(rec_elems[0], 0);
  auto upidx = this->ltgT(rec_elems.back(), m_nq - 1);
  int lenidx = upidx - lowidx + 1;
  auto nrec = param.num_receivers();
  Eigen::MatrixXcd mat_base = Eigen::MatrixXcd::Zero(3 * nrec, lenidx);
  double rad_r = _ReceiverRadius(param);
  double k = std::sqrt(1.0 * idxl * (idxl + 1.0));
  for (int idx = rec_elems[0]; idx < rec_elems.back() + 1; ++idx) {
    std::vector<double> vec_nodes(m_mesh.NN(), 0.0);
    for (int idxn = 0; idxn < m_mesh.NN(); ++idxn)
      vec_nodes[idxn] = m_mesh.NodeRadius(idx, idxn);
    auto pleg =
        Interpolation::LagrangePolynomial(vec_nodes.begin(), vec_nodes.end());
    for (int idxq = 0; idxq < m_mesh.NN(); ++idxq) {
      auto idx_v = this->ltgT(idx, idxq) - lowidx;
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

#endif   // SEM_RECEIVERS_H
