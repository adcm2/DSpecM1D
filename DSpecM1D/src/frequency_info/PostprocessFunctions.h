#ifndef DSPECM1D_FREQUENCY_INFO_POSTPROCESSFUNCTIONS_H
#define DSPECM1D_FREQUENCY_INFO_POSTPROCESSFUNCTIONS_H

#include <FFTWpp/Ranges>
#include <Eigen/Core>
#include <cmath>
#include <complex>

#include "FreqFull.h"

namespace filters {

template <typename xtype>
xtype
hannref(const xtype &f, const xtype &f1, const xtype &f2, const xtype &fac) {
  const xtype f11 = f1;
  const xtype f22 = f2;
  const xtype f12 = f11 + fac * (f22 - f11);
  const xtype f21 = f22 - fac * (f22 - f11);
  const bool ltmp = (f11 == 0.0 && f12 == 0.0 && f21 == 0.0 && f22 == 0.0);

  if (!ltmp) {
    if (f < f11) {
      return static_cast<xtype>(0.0);
    } else if (f >= f11 && f < f12) {
      const xtype tmp =
          static_cast<xtype>(3.1415926535) * (f - f11) / (f12 - f11);
      return static_cast<xtype>(0.5 * (1.0 - std::cos(tmp)));
    } else if (f >= f12 && f < f21) {
      return static_cast<xtype>(1.0);
    } else if (f >= f21 && f < f22) {
      const xtype tmp =
          static_cast<xtype>(3.1415926535) * (f22 - f) / (f22 - f21);
      return static_cast<xtype>(0.5 * (1.0 - std::cos(tmp)));
    } else {
      return static_cast<xtype>(0.0);
    }
  }

  return static_cast<xtype>(0.0);
}

template <typename xtype>
xtype
hannref(const xtype &f, const xtype &f11, const xtype &f12, const xtype &f21,
        const xtype &f22) {
  const bool ltmp = (f11 == 0.0 && f12 == 0.0 && f21 == 0.0 && f22 == 0.0);

  if (!ltmp) {
    if (f < f11) {
      return static_cast<xtype>(0.0);
    } else if (f >= f11 && f < f12) {
      const xtype tmp = static_cast<xtype>(3.1415926535897932) * (f - f11) /
                        (f12 - f11);
      return static_cast<xtype>(0.5 * (1.0 - std::cos(tmp)));
    } else if (f >= f12 && f < f21) {
      return static_cast<xtype>(1.0);
    } else if (f >= f21 && f < f22) {
      const xtype tmp = static_cast<xtype>(3.1415926535897932) * (f22 - f) /
                        (f22 - f21);
      return static_cast<xtype>(0.5 * (1.0 - std::cos(tmp)));
    } else {
      return static_cast<xtype>(0.0);
    }
  }

  return static_cast<xtype>(0.0);
}

}  // namespace filters

namespace processfunctions {

using Float = double;
using Complex = std::complex<Float>;
using RealVector = FFTWpp::vector<Float>;
using ComplexVector = FFTWpp::vector<Complex>;

inline Eigen::MatrixXd
rawfreq2time(const Eigen::MatrixXcd &rawspec, const int nt) {
  using namespace FFTWpp;

  Eigen::MatrixXd tmp(rawspec.rows(), nt);
  ComplexVector inFL(nt / 2 + 1);
  RealVector outFL(nt);
  auto planForward =
      Ranges::Plan(Ranges::View(inFL), Ranges::View(outFL), FFTWpp::Estimate);

  for (int idx = 0; idx < rawspec.rows(); ++idx) {
    for (int idx2 = 0; idx2 < nt / 2 + 1; ++idx2) {
      inFL[idx2] = rawspec(idx, idx2);
    }
    planForward.Execute();
    for (int idx2 = 0; idx2 < nt; ++idx2) {
      tmp(idx, idx2) = outFL[idx2];
    }
  }

  return tmp;
}

inline Eigen::MatrixXcd
rawtime2freq(const Eigen::MatrixXd &rawspec, const int nt, double dt) {
  using namespace FFTWpp;

  RealVector inFL(nt);
  ComplexVector outFL(nt / 2 + 1);
  auto planForward =
      Ranges::Plan(Ranges::View(inFL), Ranges::View(outFL), FFTWpp::Estimate);

  Eigen::MatrixXcd tmp(rawspec.rows(), nt / 2 + 1);
  for (int idx = 0; idx < rawspec.rows(); ++idx) {
    for (int idx2 = 0; idx2 < nt; ++idx2) {
      inFL[idx2] = rawspec(idx, idx2);
    }
    planForward.Execute();
    for (int idx2 = 0; idx2 < nt / 2 + 1; ++idx2) {
      tmp(idx, idx2) = outFL[idx2] * dt;
    }
  }

  return tmp;
}

inline Eigen::MatrixXd
freq2time(const Eigen::MatrixXcd &rawspec,
          const SpectraSolver::FreqFull &calcdata, bool undo_exponent = true) {
  auto tmp = rawfreq2time(rawspec, calcdata.nt());

  if (undo_exponent) {
    for (int idx = 0; idx < calcdata.nt(); ++idx) {
      tmp.col(idx) *= std::exp(calcdata.ep() * calcdata.dt() * idx) *
                      calcdata.df();
    }
  } else {
    tmp *= calcdata.df();
  }

  return tmp;
}

inline Eigen::MatrixXd
filtfreq2time(const Eigen::MatrixXcd &rawspec,
              const SpectraSolver::FreqFull &calcdata,
              bool undo_exponent = true) {
  const int nrow = rawspec.rows();
  Eigen::MatrixXd tmp(nrow, calcdata.nt());
  Eigen::MatrixXcd tmpraw = rawspec;

  for (int idx = 0; idx < calcdata.nt() / 2 + 1; ++idx) {
    tmpraw.block(0, idx, nrow, 1) *=
        filters::hannref(calcdata.df() * idx, calcdata.f11(), calcdata.f12(),
                         calcdata.f21(), calcdata.f22());
  }

  tmp = rawfreq2time(tmpraw, calcdata.nt());

  if (undo_exponent) {
    for (int idx = 0; idx < calcdata.nt(); ++idx) {
      tmp.block(0, idx, nrow, 1) *=
          std::exp(calcdata.ep() * calcdata.dt() * idx) * calcdata.df();
    }
  } else {
    tmp *= calcdata.df();
  }

  return tmp;
}

inline Eigen::MatrixXcd
simptime2freq(const Eigen::MatrixXd &rawspec, const double dt, const double t2) {
  const int nrow = rawspec.rows();
  const int nt = rawspec.cols();
  Eigen::MatrixXd tmpraw = rawspec;

  for (int idx = 0; idx < nt; ++idx) {
    tmpraw.block(0, idx, nrow, 1) *= filters::hannref(idx * dt, 0.0, t2, 0.5);
  }

  return rawtime2freq(tmpraw, nt, dt);
}

inline Eigen::MatrixXcd
fulltime2freq(const Eigen::MatrixXd &rawspec,
              const SpectraSolver::FreqFull &calcdata) {
  const int nrow = rawspec.rows();
  Eigen::MatrixXd tmpraw = Eigen::MatrixXd::Zero(nrow, calcdata.nt0());
  tmpraw.block(0, 0, nrow, calcdata.nt()) = rawspec;
  return simptime2freq(tmpraw, calcdata.dt(), calcdata.t2());
}

inline Eigen::MatrixXcd
simptime2freq(const Eigen::MatrixXd &rawspec, const double dt, const double t2,
              double hval) {
  const int nrow = rawspec.rows();
  const int nt = rawspec.cols();
  Eigen::MatrixXd tmpraw = rawspec;

  for (int idx = 0; idx < nt; ++idx) {
    tmpraw.block(0, idx, nrow, 1) *=
        filters::hannref(idx * dt, 0.0, t2, hval);
  }

  return rawtime2freq(tmpraw, nt, dt);
}

inline Eigen::MatrixXcd
fulltime2freq(const Eigen::MatrixXd &rawspec,
              const SpectraSolver::FreqFull &calcdata, double hval) {
  const int nrow = rawspec.rows();
  Eigen::MatrixXd tmpraw = Eigen::MatrixXd::Zero(nrow, calcdata.nt0());
  tmpraw.block(0, 0, nrow, calcdata.nt()) = rawspec;
  return simptime2freq(tmpraw, calcdata.dt(), calcdata.t2(), hval);
}

}  // namespace processfunctions

#endif  // DSPECM1D_FREQUENCY_INFO_POSTPROCESSFUNCTIONS_H
