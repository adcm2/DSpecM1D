#ifndef START_ELEMENT_H
#define START_ELEMENT_H

#include <algorithm>   // std::max, std::min
#include <cmath>       // std::sqrt, std::abs
#include <vector>      // std::vector
#include "SourceInfo.h"
#include <SpectraSolver/FF>

namespace SpectralTools {

namespace detail {

/// Shared implementation for allIndicesTor and allIndicesSph.
/// GetStartIdx is a callable: (double omega) -> std::size_t
/// that returns the global DOF index of the lowest element at that frequency.
template <class GetStartIdx>
std::vector<int>
allIndicesImpl(const std::vector<double> &vec_w,   // frequencies to evaluate
               int i_begin,                        // first index into vec_w
               int i_end,                          // one-past-last index
               int nskip, GetStartIdx getStartIdx) {
  int len = i_end - i_begin;
  std::vector<int> tmp(len, 0);

  for (int idx = i_end - 1; idx >= i_begin; --idx) {
    int ovidx = idx - i_begin;
    int idxn = i_end - 1 - idx;   // counts up from 0 as idx decreases

    if ((idxn % nskip) == 0) {
      tmp[ovidx] = static_cast<int>(getStartIdx(vec_w[idx]));
    } else {
      tmp[ovidx] = tmp[ovidx + 1];   // reuse neighbour
    }
  }
  return tmp;
}

}   // namespace detail

template <class semtype>
auto
startElementClean(semtype &sem, int l, double omega_in, bool isP,
                  int idx_source = -1) {
  auto mesh = sem.mesh();
  auto meshModel = sem.meshModel();
  auto q = mesh.GLL();
  int NE = mesh.NE();
  int NQ = mesh.NN();
  double omega = omega_in;
  double slowval2 = l * (l + 1.0) / (omega * omega);
  double sum = 0.0;
  // use a slightly tighter tolerance when starting from source element
  double tolval = (idx_source < 0) ? 12.0 : 14.0;
  int idx = 0;

  // start from source element if provided, otherwise from the top
  int idx_start =
      (idx_source >= 0 && idx_source < NE) ? idx_source - 1 : NE - 1;

  for (int idxe = idx_start; idxe > -1; --idxe) {
    double tmpmult = 0.5 * omega * mesh.EW(idxe);
    for (int idxq = 0; idxq < NQ; ++idxq) {
      double r = mesh.NodeRadius(idxe, idxq);
      auto sv = -slowval2 / (r * r);
      if (isP) {
        sv += meshModel.PSlow(idxe, idxq);
      } else {
        sv += meshModel.SSlow(idxe, idxq);
      }
      if (sv < 0) {
        sum += tmpmult * std::sqrt(-sv) * q.W(idxq);
      } else {
        sum = 0.0;
        break;
      }
    }
    if (sum > tolval) {
      idx = idxe;
      break;
    }
  }
  return idx;
}

template <class semtype>
auto
startElementTor(semtype &sem, int l, double omega_in) {
  auto idx1 = startElementClean(sem, l, omega_in, 0);
  auto idx2 = sem.el();
  return std::max(idx1, idx2);
}

template <class semtype>
auto
startElementTor(semtype &sem, int l, double omega_in, int idx_source) {
  auto idx1 = startElementClean(sem, l, omega_in, 0, idx_source);
  auto idx2 = sem.el();
  return std::max(idx1, idx2);
}

template <class semtype>
auto
startElementSph(semtype &sem, int l, double omega_in) {
  auto idx1 = startElementClean(sem, l, omega_in, 0);
  auto idx2 = startElementClean(sem, l, omega_in, 1);
  return std::min(idx1, idx2);
}

template <class semtype>
auto
startElementSph(semtype &sem, int l, double omega_in, int idx_source) {
  auto idx1 = startElementClean(sem, l, omega_in, 0, idx_source);
  auto idx2 = startElementClean(sem, l, omega_in, 1, idx_source);
  return std::min(idx1, idx2);
}

template <class semtype>
auto
allIndicesTor(semtype &sem, int l, SpectraSolver::FreqFull &myff,
               int nskip = 10) {
  auto vec_w = myff.w();
  return detail::allIndicesImpl(vec_w, myff.i1(), myff.i2(), nskip,
                                [&](double w) {
                                  auto idxlow_e = startElementTor(sem, l, w);
                                  return sem.ltgT(idxlow_e, 0);
                                });
}

template <class semtype>
auto
allIndicesTor(semtype &sem, int l, SpectraSolver::FreqFull &myff,
               int idx_source, int nskip) {
  auto vec_w = myff.w();
  return detail::allIndicesImpl(
      vec_w, myff.i1(), myff.i2(), nskip, [&](double w) {
        auto idxlow_e = startElementTor(sem, l, w, idx_source);
        return sem.ltgT(idxlow_e, 0);
      });
}

template <class semtype>
auto
allIndicesTor(semtype &sem, int l, std::vector<double> &vec_w, int idx_source,
               int nskip) {
  return detail::allIndicesImpl(
      vec_w, 0, static_cast<int>(vec_w.size()), nskip, [&](double w) {
        auto idxlow_e = startElementTor(sem, l, w, idx_source);
        return sem.ltgT(idxlow_e, 0);
      });
}

template <class semtype>
auto
allIndicesSph(semtype &sem, int l, SpectraSolver::FreqFull &myff,
               int nskip = 10) {
  auto vec_w = myff.w();
  return detail::allIndicesImpl(vec_w, myff.i1(), myff.i2(), nskip,
                                [&](double w) {
                                  auto idxlow_e = startElementSph(sem, l, w);
                                  return sem.ltgS(0, idxlow_e, 0);
                                });
}

template <class semtype>
auto
allIndicesSph(semtype &sem, int l, SpectraSolver::FreqFull &myff,
               int idx_source, int nskip) {
  auto vec_w = myff.w();
  return detail::allIndicesImpl(
      vec_w, myff.i1(), myff.i2(), nskip, [&](double w) {
        auto idxlow_e = startElementSph(sem, l, w, idx_source);
        return sem.ltgS(0, idxlow_e, 0);
      });
}

template <class semtype>
auto
allIndicesSph(semtype &sem, int l, std::vector<double> &vec_w, int idx_source,
               int nskip) {
  return detail::allIndicesImpl(
      vec_w, 0, static_cast<int>(vec_w.size()), nskip, [&](double w) {
        auto idxlow_e = startElementSph(sem, l, w, idx_source);
        return sem.ltgS(0, idxlow_e, 0);
      });
};

template <class model1d, class semtype>
auto
startElementS(semtype &sem, model1d &inp_model, int l, double omega_in) {

  auto mesh = sem.mesh();
  int NE = mesh.NE();
  int NQ = mesh.NN();

  double omega = omega_in * inp_model.TimeNorm();

  double slowval2 =
      static_cast<double>(l) * static_cast<double>(l + 1) / (omega * omega);

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
  int idx_S = 0;   // computed but never returned — is this intentional?
  bool found_elem = false;

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
      auto rq1 = elr + ew * vec_vals[idxq];
      auto rq0 = elr + ew * vec_vals[idxq - 1];
      auto rqd = rq1 - rq0;

      auto rp1 = Slow_P(laynum, rq1);
      auto rp0 = Slow_P(laynum, rq0);
      auto rs1 = Slow_S(laynum, rq1);
      auto rs0 = Slow_S(laynum, rq0);

      if ((rp1 < 0) && (rp0 < 0)) {
        auto qs1 = std::sqrt(-rp1);
        auto qs0 = std::sqrt(-rp0);
        sum_P += omega * std::abs(0.5 * (qs1 + qs0) * rqd);
      } else {
        sum_P = 0.0;
      }

      if ((rs1 < 0) && (rs0 < 0)) {
        auto qs1 = std::sqrt(-rs1);
        auto qs0 = std::sqrt(-rs0);
        sum_S += omega * std::abs(0.5 * (qs1 + qs0) * rqd);
      } else {
        sum_S = 0.0;
      }

      if ((sum_P > tolval) && (sum_S > tolval)) {
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
  return idx_P;   // idx_S is silently discarded
};

template <class semtype>
auto
startRadiusClean(semtype &sem, int l, double omega_in, bool isP) {

  auto mesh = sem.mesh();
  auto meshModel = sem.meshModel();
  int NE = mesh.NE();
  int NQ = mesh.NN();
  double omega = omega_in;
  double slowval2 = l * (l + 1.0) / (omega * omega);
  double sum = 0.0;
  double tolval = 13.0;
  double rtmp = 0.0;
  bool texit = false;

  // loop through elements from top to bottom
  for (int idxe = NE - 1; idxe > -1; --idxe) {
    for (int idxq = NQ - 1; idxq > 0; --idxq) {
      double ru = mesh.NodeRadius(idxe, idxq);
      double rl = mesh.NodeRadius(idxe, idxq - 1);
      auto sv = -slowval2 / (ru * ru);
      auto sl = -slowval2 / (rl * rl);
      auto vv = 0.0;
      auto vl = 0.0;
      if (isP) {
        vv = meshModel.VP(idxe, idxq);
        vl = meshModel.VP(idxe, idxq - 1);
      } else {
        vv = meshModel.VS(idxe, idxq);
        vl = meshModel.VS(idxe, idxq - 1);
      }
      sv += 1.0 / (vv * vv);
      sl += 1.0 / (vl * vl);

      if ((sv < 0) && (sl < 0)) {
        sum += omega * (std::sqrt(-sv) + std::sqrt(-sl)) * 0.5 * (ru - rl);
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
}

template <class model1d, class semtype>
auto
StartElement(semtype &sem, model1d &inp_model, int l, double omega_in,
             bool isP) {

  auto mesh = sem.mesh();
  int NE = mesh.NE();
  int NQ = mesh.NN();
  double omega = omega_in;

  double slowval2 =
      static_cast<double>(l) * static_cast<double>(l + 1) / (omega * omega);

  auto Slow_VAL = [slowval2, &inp_model, isP](int laynum, double r) {
    double cval;
    if (isP) {
      cval = inp_model.VP(laynum)(r);
    } else {
      cval = inp_model.VS(laynum)(r);
    }
    return 1.0 / (cval * cval) - slowval2 / (r * r);
  };

  double sum = 0.0;
  double tolval = 11.0;
  int idx = 0;
  bool found_elem = false;

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
      auto rq1 = elr + ew * vec_vals[idxq];
      auto rq0 = elr + ew * vec_vals[idxq - 1];
      auto rqd = rq1 - rq0;

      auto rs1 = Slow_VAL(laynum, rq1);
      auto rs0 = Slow_VAL(laynum, rq0);

      if ((rs1 < 0) && (rs0 < 0)) {
        auto qs1 = std::sqrt(-rs1);
        auto qs0 = std::sqrt(-rs0);
        sum += omega * std::abs(0.5 * (qs1 + qs0) * rqd);
      } else {
        sum = 0.0;
      }
    }
    if ((sum > tolval)) {
      idx = idxe;
      break;
    }
  }
  return idx;
};

}   // namespace SpectralTools

#endif   // START_ELEMENT_H
