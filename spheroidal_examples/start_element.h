#ifndef START_ELEMENT_H
#define START_ELEMENT_H

#include <iostream>
#include <string>
#include "sem_full.h"
#include "../SpectraSolver/SpectraSolver/FF"

namespace SpectralTools {

template <class model1d>
auto
StartElement(Full1D::sem &sem, model1d &inp_model, int l, double omega_in,
             bool isP) {

  auto mesh = sem.mesh();
  int NE = mesh.NE();
  int NQ = mesh.NN();

  // omega_in is dimensional frequency. Divide by frequency norm = 1/T
  //   double omega = omega_in * inp_model.TimeNorm();
  double omega = omega_in;

  double slowval2 =
      static_cast<double>(l) * static_cast<double>(l + 1) / (omega * omega);

  if (isP) {
    // std::cout << "Calculating Start Element for P waves, l=" << l << "\n";
  } else {
    // std::cout << "Calculating Start Element for S waves, l=" << l << "\n";
  }
  auto Slow_VAL = [slowval2, &inp_model, isP](int laynum, double r) {
    double cval;
    if (isP) {
      cval = inp_model.VP(laynum)(r);
    } else {
      cval = inp_model.VS(laynum)(r);
    }
    double tmp = 1.0 / (cval * cval) - slowval2 / (r * r);
    return tmp;
  };

  double sum = 0.0;
  double tolval = 11.0;
  int idx = 0;
  bool found_elem = false;

  // increase number of points compared to the spectral element nodes?
  int npts = 20;
  std::vector<double> vec_vals(npts);
  for (int i = 0; i < npts; ++i) {
    vec_vals[i] = static_cast<double>(i) / static_cast<double>(npts - 1);
  }

  // loop through elements from top to bottom
  for (int idxe = NE - 1; idxe > -1; --idxe) {
    auto laynum = mesh.LayerNumber(idxe);
    auto elr = mesh.ELR(idxe);
    auto eur = mesh.EUR(idxe);
    auto ew = eur - elr;

    // loop through nodes within element from top to bottom
    for (int idxq = npts - 1; idxq > 0; --idxq) {
      // node radii
      auto rq1 = elr + ew * vec_vals[idxq];
      auto rq0 = elr + ew * vec_vals[idxq - 1];
      auto rqd = rq1 - rq0;

      // slowness
      auto rs1 = Slow_VAL(laynum, rq1);
      auto rs0 = Slow_VAL(laynum, rq0);

      // testing
      if ((rs1 < 0) && (rs0 < 0)) {
        auto qs1 = std::sqrt(-rs1);
        auto qs0 = std::sqrt(-rs0);
        auto rs = 0.5 * (qs1 + qs0) * rqd;
        sum += omega * std::abs(rs);   // convert to minutes

      } else {
        sum = 0.0;
      }
    }
    if ((sum > tolval)) {
      //   std::cout << "l: " << l << ", idxe: " << idxe << ", sum: " << sum
      //             << "\n\n";
      idx = idxe;
      break;
    }
  }
  return idx;
};

template <class semtype>
auto
StartElementClean(semtype &sem, int l, double omega_in, bool isP) {

  auto mesh = sem.mesh();
  auto mesh_model = sem.mesh_model();
  auto q = mesh.GLL();
  int NE = mesh.NE();
  int NQ = mesh.NN();
  double omega = omega_in;
  double slowval2 = l * (l + 1.0) / (omega * omega);
  double sum = 0.0;
  double tolval = 12.0;
  int idx = 0;

  // auto rturnmax = std::sqrt(slowval2) * 5800000.0 / sem.VelocityNorm();
  // std::cout << "Max turning radius for l=" << l << ", omega: " << omega_in
  //           << ": " << rturnmax << "\n";
  // loop through elements from top to bottom
  for (int idxe = NE - 1; idxe > -1; --idxe) {

    // auto tmp = 0.0;
    double tmpmult = 0.5 * omega * mesh.EW(idxe);
    for (int idxq = 0; idxq < NQ; ++idxq) {
      double r = mesh.NodeRadius(idxe, idxq);
      // auto sv = Slow_VAL(idxe, idxq);
      auto sv = -slowval2 / (r * r);
      auto vv = 0.0;
      if (isP) {
        // vv = mesh_model.VP(idxe, idxq);
        sv += mesh_model.PSlow(idxe, idxq);
      } else {
        // vv = mesh_model.VS(idxe, idxq);
        sv += mesh_model.SSlow(idxe, idxq);
      }
      // sv += 1.0 / (vv * vv);

      if (sv < 0) {
        sum += tmpmult * std::sqrt(-sv) * q.W(idxq);
      } else {
        sum = 0.0;
        break;
      }
    }
    // sum += omega * tmp * mesh.EW(idxe) * 0.5;

    if ((sum > tolval)) {
      idx = idxe;
      break;
    }
  }
  return idx;
};

template <class semtype>
auto
StartRadiusClean(semtype &sem, int l, double omega_in, bool isP) {

  auto mesh = sem.mesh();
  auto mesh_model = sem.mesh_model();
  int NE = mesh.NE();
  int NQ = mesh.NN();
  double omega = omega_in;
  double slowval2 = l * (l + 1.0) / (omega * omega);
  double sum = 0.0;
  double tolval = 12.0;
  double rtmp = 0.0;
  bool texit = false;

  // loop through elements from top to bottom
  for (int idxe = NE - 1; idxe > -1; --idxe) {
    for (int idxq = NQ - 1; idxq > 0; --idxq) {
      double ru = mesh.NodeRadius(idxe, idxq);
      double rl = mesh.NodeRadius(idxe, idxq - 1);
      // auto sv = Slow_VAL(idxe, idxq);
      auto sv = -slowval2 / (ru * ru);
      auto sl = -slowval2 / (rl * rl);
      auto vv = 0.0;
      auto vl = 0.0;
      if (isP) {
        vv = mesh_model.VP(idxe, idxq);
        vl = mesh_model.VP(idxe, idxq - 1);
        // sv += mesh_model.PSlow(idxe, idxq);
      } else {
        vv = mesh_model.VS(idxe, idxq);
        vl = mesh_model.VS(idxe, idxq - 1);
        // sv += mesh_model.SSlow(idxe, idxq);
      }
      sv += 1.0 / (vv * vv);
      sl += 1.0 / (vl * vl);

      if ((sv < 0) && (sl < 0)) {
        sum += omega * (std::sqrt(-sv) + std::sqrt(-sl)) * 0.5 *
               (ru - rl);   // convert to minutes

        if (sum > tolval) {
          rtmp = rl;
          texit = true;
          break;
        }
      } else {
        sum = 0.0;
        break;
      }
    }

    if (texit) {
      break;
    }
  }
  return rtmp;
};

template <class model1d, class semtype>
auto
StartElement_Tor(semtype &sem, model1d &inp_model, int l, double omega_in) {
  auto idx1 = StartElementClean(sem, l, omega_in, 0);
  auto idx2 = sem.EL();
  return std::max(idx1, idx2);
};

template <class model1d, class semtype>
auto
StartElement_Sph(semtype &sem, model1d &inp_model, int l, double omega_in) {
  auto idx1 = StartElementClean(sem, l, omega_in, 0);
  auto idx2 = StartElementClean(sem, l, omega_in, 1);
  return std::min(idx1, idx2);
};

template <class model1d, class semtype>
auto
AllIndices_TOR(semtype &sem, model1d &inp_model, int l,
               SpectraSolver::FreqFull &myff, int nskip = 10) {
  auto i1 = myff.i1();
  auto i2 = myff.i2();
  std::vector<int> tmp(i2 - i1, 0);
  for (int idx = i2 - 1; idx > i1 - 1; --idx) {
    auto ovidx = idx - i1;
    auto idxn = i2 - 1 - idx;
    if ((idxn % nskip) == 0) {
      auto idxlow_e = StartElement_Tor(sem, inp_model, l, myff.w(idx));
      std::size_t ridx = sem.LtG_T(idxlow_e, 0);
      tmp[ovidx] = ridx;
    } else {
      tmp[ovidx] = tmp[ovidx + 1];
    }
  }
  return tmp;
};

template <class model1d, class semtype>
auto
AllIndices_SPH(semtype &sem, model1d &inp_model, int l,
               SpectraSolver::FreqFull &myff, int nskip = 10) {
  auto i1 = myff.i1();
  auto i2 = myff.i2();
  std::vector<int> tmp(i2 - i1, 0);
  for (int idx = i2 - 1; idx > i1 - 1; --idx) {
    auto ovidx = idx - i1;
    auto idxn = i2 - 1 - idx;
    if ((idxn % nskip) == 0) {
      auto idxlow_e = StartElement_Sph(sem, inp_model, l, myff.w(idx));
      std::size_t ridx = sem.LtG_S(0, idxlow_e, 0);
      tmp[ovidx] = ridx;
    } else {
      tmp[ovidx] = tmp[ovidx + 1];
    }
  }
  return tmp;
};

template <class model1d>
auto
StartElement_S(Full1D::sem &sem, model1d &inp_model, int l, double omega_in) {

  auto mesh = sem.mesh();
  int NE = mesh.NE();
  int NQ = mesh.NN();

  // omega_in is dimensional frequency. Divide by frequency norm = 1/T
  double omega = omega_in * inp_model.TimeNorm();

  double slowval2 =
      static_cast<double>(l) * static_cast<double>(l + 1) / (omega * omega);
  //   std::cout << "Slowval2: " << slowval2 << "\n";
  auto Slow_P = [slowval2, &inp_model](int laynum, double r) {
    double cval = inp_model.VP(laynum)(r);
    double tmp = 1.0 / (cval * cval) - slowval2 / (r * r);
    return tmp;
  };

  auto Slow_S = [slowval2, &inp_model](int laynum, double r) {
    double cval = inp_model.VS(laynum)(r);
    double tmp = 1.0 / (cval * cval) - slowval2 / (r * r);
    return tmp;
  };

  double sum_P = 0.0;
  double sum_S = 0.0;
  double tolval = 12.0;
  int idx_P = 0;
  int idx_S = 0;
  bool found_elem = false;

  // increase number of points compared to the spectral element nodes?
  int npts = 10;
  std::vector<double> vec_vals(npts);
  for (int i = 0; i < npts; ++i) {
    vec_vals[i] = static_cast<double>(i) / static_cast<double>(npts - 1);
  }

  // loop through elements from top to bottom
  for (int idxe = NE - 1; idxe > -1; --idxe) {
    auto laynum = mesh.LayerNumber(idxe);
    auto elr = mesh.ELR(idxe);
    auto eur = mesh.EUR(idxe);
    auto ew = eur - elr;

    // loop through nodes within element from top to bottom
    for (int idxq = npts - 1; idxq > 0; --idxq) {
      // node radii
      auto rq1 = elr + ew * vec_vals[idxq];
      auto rq0 = elr + ew * vec_vals[idxq - 1];
      auto rqd = rq1 - rq0;

      // slowness
      auto rp1 = Slow_P(laynum, rq1);
      auto rp0 = Slow_P(laynum, rq0);
      auto rs1 = Slow_S(laynum, rq1);
      auto rs0 = Slow_S(laynum, rq0);
      //   if (idxe == (NE - 1)) {
      //     std::cout << "Top element slowness values at l=" << l << ": " <<
      //     rs0
      //               << ", " << rs1 << "\n";
      //   }

      // testing
      if ((rp1 < 0) && (rp0 < 0)) {
        auto qs1 = std::sqrt(-rp1);
        auto qs0 = std::sqrt(-rp0);
        auto rs = 0.5 * (qs1 + qs0) * rqd;
        sum_P += omega * std::abs(rs);   // convert to minutes

      } else {
        sum_P = 0.0;
      }

      if ((rs1 < 0) && (rs0 < 0)) {
        auto qs1 = std::sqrt(-rs1);
        auto qs0 = std::sqrt(-rs0);
        auto rs = 0.5 * (qs1 + qs0) * rqd;
        sum_S += omega * std::abs(rs);   // convert to minutes

      } else {
        sum_S = 0.0;
      }

      if ((sum_P > tolval) && (sum_S > tolval)) {
        // std::cout << "l: " << l << ", idxe: " << idxe << ", idxq: " << idxq
        //           << ", sum_P: " << sum_P << ", sum_S: " << sum_S << "\n\n";
        idx_P = idxe;
        idx_S = idxe;
        found_elem = true;
        break;
      }
    }
    if (found_elem) {
      break;
    }
  }
  return idx_P;
};
};   // namespace SpectralTools

#endif   // START_ELEMENT_H