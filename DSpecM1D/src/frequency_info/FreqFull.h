#ifndef DSPECM1D_FREQUENCY_INFO_FREQFULL_H
#define DSPECM1D_FREQUENCY_INFO_FREQFULL_H

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

namespace SpectraSolver {

class FreqFull {
public:
  using VECTOR = Eigen::VectorXcd;
  using MATRIX = Eigen::MatrixXcd;
  using COMPLEX = std::complex<double>;

  FreqFull() {};
  ~FreqFull() {};

  FreqFull(double f1, double f2, double dt, double tout, double df0,
           double wtb, double t1, double t2, int qex,
           double TimeNorm = 1.0)
      : FreqFull(f1, f2, f1, f2, f2, f2, dt, tout, df0, wtb, t1, t2, qex,
                 TimeNorm) {}

  FreqFull(double f1, double f2, double f11, double f12, double f21,
           double f22, double dt, double tout, double df0, double wtb,
           double t1, double t2, int qex, double TimeNorm = 1.0)
      : m_f1{f1 / 1000.0 * TimeNorm},
        m_f2{f2 / 1000.0 * TimeNorm},
        m_f11{f11 / 1000.0 * TimeNorm},
        m_f12{f12 / 1000.0 * TimeNorm},
        m_f21{f21 / 1000.0 * TimeNorm},
        m_f22{f22 / 1000.0 * TimeNorm},
        m_dt{dt / TimeNorm},
        m_df0{df0 / 1000.0 * TimeNorm},
        m_tout{tout * 3600.0 / TimeNorm},
        m_t1{t1 * 3600.0 / TimeNorm},
        m_t2{std::min(t2, tout) * 3600.0 / TimeNorm},
        m_wtb{wtb * 3.1415926535 / 500.0 * TimeNorm},
        m_timenorm{TimeNorm},
        m_frequencynorm{1.0 / TimeNorm} {
    double fn = 0.5 / m_dt;
    if (fn < m_f2) {
      std::cout << "f2 is greater than the Nyquist frequency for the time "
                   "step. Behaviour may be unexpected"
                << std::endl;
      m_f2 = fn;
    }

    int mex = 5;
    m_ep = mex / m_tout;
    m_df = m_ep / (6.28318530718 * qex);

    m_nt = std::ceil(1.0 / (m_df * m_dt));
    int ne = static_cast<int>(std::log(static_cast<double>(m_nt)) /
                                  std::log(2.0) +
                              1);
    m_nt = std::pow(2, ne);
    m_df = 1.0 / (m_nt * m_dt);

    m_i1 = std::max(static_cast<int>(std::floor(m_f1 / m_df)), 2);
    m_i2 = std::min(static_cast<int>(std::floor(m_f2 / m_df)) + 2, m_nt);
    m_i1 -= 1;

    m_w.reserve(m_nt / 2 + 1);
    for (int idx = 0; idx < m_nt / 2 + 1; ++idx) {
      m_w.push_back(2.0 * 3.1415926535 * m_df * static_cast<double>(idx));
    }

    m_t.reserve(m_nt);
    for (int idx = 0; idx < m_nt; ++idx) {
      m_t.push_back(m_dt * static_cast<double>(idx));
    }

    m_nt0 = std::floor(1.0 / m_df0 * m_dt);
    if (m_nt0 > m_nt) {
      int ne2 = std::log(static_cast<double>(m_nt0)) / std::log(2.0) + 1;
      m_nt0 = std::pow(2, ne2);
    } else {
      m_nt0 = m_nt;
    }

    m_df2 = 1.0 / (m_nt0 * m_dt);
    m_i12 = std::max(static_cast<int>(std::floor(m_f1 / m_df2)) - 1, 0);
    m_i22 = static_cast<int>(std::floor(m_f2 / m_df2)) + 1;
  }

  double f(int idx) const { return static_cast<double>(idx) * m_df; }
  double f0(int idx) const { return static_cast<double>(idx) * m_df0; }
  double f2(int idx) const { return static_cast<double>(idx) * m_df2; }
  double f1() const { return m_f1; }
  double f2() const { return m_f2; }
  double f12() const { return m_f12; }
  double f11() const { return m_f11; }
  double f21() const { return m_f21; }
  double tout() const { return m_tout; }
  double f22() const { return m_f22; }
  double df() const { return m_df; }
  double df0() const { return m_df0; }
  double df2() const { return m_df2; }
  double wtb() const { return m_wtb; }
  double t1() const { return m_t1; }
  double t2() const { return m_t2; }
  double ep() const { return m_ep; }
  double dt() const { return m_dt; }

  int nt() const { return m_nt; }
  int nt0() const { return m_nt0; }
  int i1() const { return m_i1; }
  int i2() const { return m_i2; }
  int i12() const { return m_i12; }
  int i22() const { return m_i22; }

  std::vector<double> w() const { return m_w; }
  std::vector<double> t() const { return m_t; }
  double w(int idx) const { return m_w[idx]; }
  double t(int idx) const { return m_t[idx]; }

private:
  double m_f1 = 0.0;
  double m_f2 = 0.0;
  double m_tout = 0.0;
  double m_df0 = 0.0;
  double m_wtb = 0.0;
  double m_t1 = 0.0;
  double m_t2 = 0.0;
  double m_df = 0.0;
  double m_ep = 0.0;
  double m_df2 = 0.0;
  double m_dt = 0.0;
  double m_timenorm = 1.0;
  double m_frequencynorm = 1.0;
  double m_f12 = 0.0;
  double m_f21 = 0.0;
  double m_f11 = 0.0;
  double m_f22 = 0.0;
  int m_nt = 0;
  int m_nt0 = 0;
  int m_i1 = 0;
  int m_i2 = 0;
  int m_i12 = 0;
  int m_i22 = 0;
  std::vector<double> m_w;
  std::vector<double> m_t;
};

}  // namespace SpectraSolver

#endif  // DSPECM1D_FREQUENCY_INFO_FREQFULL_H
