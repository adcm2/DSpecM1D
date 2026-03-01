#include <iostream>
#include <fstream>
#include <vector>
#include <utility>   // for std::pair
#include <complex>
#include <cmath>

// class sph_basis {
// public:
//   sph_basis() {};

// private:
//   std::vector<std::complex<double>> m_vec(3);
// };

// class thetavec {
// public:
//   thetavec() {};

// private:
//   std::vector<std::complex<double>> m_vec{1 / std::sqrt{2}, 0,
//                                           -1 / std::sqrt{2}};
// };
namespace TensorClasses {
class tensormult {
public:
  tensormult(int v1, int v2, int v3, int v4)
      : _v1(v1), _v2(v2), _v3(v3), _v4(v4) {
    std::vector<std::complex<double>> tmp1{0.0, 1.0, 0.0};
    std::vector<std::complex<double>> tmp2{1.0 / std::sqrt(2.0), 0.0,
                                           -1.0 / std::sqrt(2.0)};
    std::vector<std::complex<double>> tmp3{
        std::complex<double>(0.0, -1.0 / std::sqrt(2.0)), 0.0,
        std::complex<double>(0.0, -1.0 / std::sqrt(2.0))};
    m_tensor.push_back(tmp1);
    m_tensor.push_back(tmp2);
    m_tensor.push_back(tmp3);
  };

  auto compval(int i1, int i2, int i3, int i4) {
    return m_tensor[_v1][i1] * m_tensor[_v2][i2] * m_tensor[_v3][i3] *
           m_tensor[_v4][i4];
  };

private:
  std::vector<std::vector<std::complex<double>>> m_tensor;
  int _v1, _v2, _v3, _v4;
};

class elastic_tensor {
public:
  elastic_tensor(double C, double A, double N, double F, double L)
      : _C(C), _A(A), _N(N), _F(F), _L(L) {};

  auto canonical_component(double i1, double i2, double i3, double i4) {
    std::complex<double> tmp1 =
        static_cast<std::complex<double>>(_C) * m_rrrr.compval(i1, i2, i3, i4);
    auto tmp2 =
        static_cast<std::complex<double>>(_A) *
        (m_tttt.compval(i1, i2, i3, i4) + m_pppp.compval(i1, i2, i3, i4));
    auto tmp3 =
        static_cast<std::complex<double>>(_N) *
        (m_tptp.compval(i1, i2, i3, i4) + m_ptpt.compval(i1, i2, i3, i4) +
         m_pttp.compval(i1, i2, i3, i4) + m_tppt.compval(i1, i2, i3, i4));
    auto tmp4 =
        static_cast<std::complex<double>>(_F) *
        (m_rrtt.compval(i1, i2, i3, i4) + m_ttrr.compval(i1, i2, i3, i4) +
         m_rrpp.compval(i1, i2, i3, i4) + m_pprr.compval(i1, i2, i3, i4));
    auto tmp5 =
        static_cast<std::complex<double>>(_L) *
        (m_rtrt.compval(i1, i2, i3, i4) + m_trtr.compval(i1, i2, i3, i4) +
         m_rprp.compval(i1, i2, i3, i4) + m_prpr.compval(i1, i2, i3, i4) +
         m_rttr.compval(i1, i2, i3, i4) + m_trrt.compval(i1, i2, i3, i4) +
         m_rppr.compval(i1, i2, i3, i4) + m_prrp.compval(i1, i2, i3, i4));
    auto tmp6 =
        static_cast<std::complex<double>>(_A - 2.0 * _N) *
        (m_ttpp.compval(i1, i2, i3, i4) + m_pptt.compval(i1, i2, i3, i4));
    return tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6;
    // return tmp1;
  };

private:
  double _C, _A, _N, _F, _L;
  tensormult m_rrrr{0, 0, 0, 0};

  tensormult m_tttt{1, 1, 1, 1};
  tensormult m_pppp{2, 2, 2, 2};

  tensormult m_ttpp{1, 1, 2, 2};
  tensormult m_pptt{2, 2, 1, 1};

  tensormult m_rrtt{0, 0, 1, 1};
  tensormult m_ttrr{1, 1, 0, 0};
  tensormult m_rrpp{0, 0, 2, 2};
  tensormult m_pprr{2, 2, 0, 0};

  tensormult m_tptp{1, 2, 1, 2};
  tensormult m_ptpt{2, 1, 2, 1};
  tensormult m_pttp{2, 1, 1, 2};
  tensormult m_tppt{1, 2, 2, 1};

  tensormult m_rtrt{0, 1, 0, 1};
  tensormult m_trtr{1, 0, 1, 0};
  tensormult m_rprp{0, 2, 0, 2};
  tensormult m_prpr{2, 0, 2, 0};
  tensormult m_rttr{0, 1, 1, 0};
  tensormult m_trrt{1, 0, 0, 1};
  tensormult m_rppr{0, 2, 2, 0};
  tensormult m_prrp{2, 0, 0, 2};
};
}   // namespace TensorClasses

int
main() {
  using namespace TensorClasses;
  tensormult tm(0, 1, 2, 0);

  elastic_tensor ec(1, 0, 0, 0, 0);
  elastic_tensor ea(0, 1, 0, 0, 0);
  elastic_tensor en(0, 0, 1, 0, 0);
  elastic_tensor ef(0, 0, 0, 1, 0);
  elastic_tensor el(0, 0, 0, 0, 1);

  //   std::cout << "C_0000: " << et.canonical_component(0, 0, 0, 0) <<
  //   std::endl;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        for (int l = 0; l < 3; ++l) {
          auto tmpc = ec.canonical_component(i, j, k, l);
          auto tmpa = ea.canonical_component(i, j, k, l);
          auto tmpn = en.canonical_component(i, j, k, l);
          auto tmpf = ef.canonical_component(i, j, k, l);
          auto tmpl = el.canonical_component(i, j, k, l);
          if (std::abs(tmpc) > 1e-10) {
            std::cout << "C_" << i - 1 << j - 1 << k - 1 << l - 1 << ": "
                      << tmpc << std::endl;
          }
          if (std::abs(tmpa) > 1e-10) {
            std::cout << "A_" << i - 1 << j - 1 << k - 1 << l - 1 << ": "
                      << tmpa << std::endl;
          }
          if (std::abs(tmpn) > 1e-10) {
            std::cout << "N_" << i - 1 << j - 1 << k - 1 << l - 1 << ": "
                      << tmpn << std::endl;
          }
          if (std::abs(tmpf) > 1e-10) {
            std::cout << "F_" << i - 1 << j - 1 << k - 1 << l - 1 << ": "
                      << tmpf << std::endl;
          }
          if (std::abs(tmpl) > 1e-10) {
            std::cout << "L_" << i - 1 << j - 1 << k - 1 << l - 1 << ": "
                      << tmpl << std::endl;
          }
          std::cout << "\n";
        }
      }
    }
  }
  return 0;
}