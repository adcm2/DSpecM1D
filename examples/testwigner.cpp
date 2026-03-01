#include <iostream>
#include <cmath>
#include <functional>
#include <filesystem>
#include <fstream>
#include <Eigen/Core>
// #include <new_coupling/All>
// #include <PlanetaryModel/All>
#include <GSHTrans/All>
#include "SourceInfo.h"

/*
class EarthquakeCMT {
public:
  EarthquakeCMT(const std::string &pathtocmt) {
    std::fstream modelfile;
    modelfile.open(pathtocmt, std::ios::in);
    if (modelfile.is_open()) {
      // extract information
      modelfile >> cmt_id >> cmt_yr >> cmt_day >> cmt_hr >> cmt_min >>
          cmt_sec >> cmt_lat >> cmt_lon >> cmt_depth >> cmt_step >> cmt_halfd >>
          cmt_m0 >> cmt_mrr >> cmt_mtt >> cmt_mpp >> cmt_mrt >> cmt_mrp >>
          cmt_mtp >> crt_mn >> cmt_s1 >> cmt_d1 >> cmt_sl1 >> cmt_s2 >>
          cmt_d2 >> cmt_sl2;
    }
  };

  // Function to print CMT information
  void PrintCMTInfo() const {
    std::cout << "CMT ID: " << cmt_id << "\n";
    std::cout << "CMT Year: " << cmt_yr << "\n";
    std::cout << "CMT Day: " << cmt_day << "\n";
    std::cout << "CMT Hour: " << cmt_hr << "\n";
    std::cout << "CMT Minute: " << cmt_min << "\n";
    std::cout << "CMT Second: " << cmt_sec << "\n";
    std::cout << "CMT Latitude: " << cmt_lat << "\n";
    std::cout << "CMT Longitude: " << cmt_lon << "\n";
    std::cout << "CMT Depth: " << cmt_depth << "\n";
    std::cout << "CMT Step: " << cmt_step << "\n";
    std::cout << "CMT Half Duration: " << cmt_halfd << "\n";
    std::cout << "CMT M0: " << cmt_m0 << "\n";
    std::cout << "CMT Mrr: " << cmt_mrr << "\n";
    std::cout << "CMT Mtt: " << cmt_mtt << "\n";
    std::cout << "CMT Mpp: " << cmt_mpp << "\n";
    std::cout << "CMT Mrt: " << cmt_mrt << "\n";
    std::cout << "CMT Mrp: " << cmt_mrp << "\n";
    std::cout << "CMT Mtp: " << cmt_mtp << "\n";
    std::cout << "CMT Moment Tensor Norm: " << crt_mn << "\n";
    std::cout << "CMT S1: " << cmt_s1 << "\n";
    std::cout << "CMT D1: " << cmt_d1 << "\n";
    std::cout << "CMT Slip 1: " << cmt_sl1 << "\n";
    std::cout << "CMT S2: " << cmt_s2 << "\n";
    std::cout << "CMT D2: " << cmt_d2 << "\n";
    std::cout << "CMT Slip 2: " << cmt_sl2 << "\n";
  }

  // print information
  auto ID() const { return cmt_id; }
  auto Year() const { return cmt_yr; }
  auto Day() const { return cmt_day; }
  auto Hour() const { return cmt_hr; }
  auto Minute() const { return cmt_min; }
  auto Second() const { return cmt_sec; }
  auto Latitude() const { return cmt_lat; }
  auto Longitude() const { return cmt_lon; }
  auto Depth() const { return cmt_depth; }
  auto Step() const { return cmt_step; }
  auto HalfDuration() const { return cmt_halfd; }
  auto M0() const { return cmt_m0; }
  auto Mrr() const { return cmt_mrr; }
  auto Mtt() const { return cmt_mtt; }
  auto Mpp() const { return cmt_mpp; }
  auto Mrt() const { return cmt_mrt; }
  auto Mrp() const { return cmt_mrp; }
  auto Mtp() const { return cmt_mtp; }
  auto MomentTensorNorm() const { return crt_mn; }
  auto Strike1() const { return cmt_s1; }
  auto Dip1() const { return cmt_d1; }
  auto Slip1() const { return cmt_sl1; }
  auto Strike2() const { return cmt_s2; }
  auto Dip2() const { return cmt_d2; }
  auto Slip2() const { return cmt_sl2; }

  // canonical moment tensor
  std::complex<double> MC00() { return cmt_mrr; };
  std::complex<double> MC0p() {
    std::complex<double> tmp(-cmt_mrt, cmt_mrp);
    return tmp / sqrt(2.0);
  };
  std::complex<double> MC0m() {
    std::complex<double> tmp(cmt_mrt, cmt_mrp);
    return tmp / sqrt(2.0);
  };
  std::complex<double> MCpp() {
    std::complex<double> tmp(0.5 * (cmt_mtt - cmt_mpp), -cmt_mtp);
    return tmp;
  };
  std::complex<double> MCmm() {
    std::complex<double> tmp(0.5 * (cmt_mtt - cmt_mpp), cmt_mtp);
    return tmp;
  };
  std::complex<double> MCmp() {
    std::complex<double> tmp(-0.5 * (cmt_mtt + cmt_mpp), 0.0);
    return tmp;
  }

private:
  std::string cmt_id;
  int cmt_yr, cmt_day, cmt_hr, cmt_min;
  double cmt_sec, cmt_lat, cmt_lon, cmt_depth, cmt_step, cmt_halfd, cmt_m0,
      cmt_mrr, cmt_mtt, cmt_mpp, cmt_mrt, cmt_mrp, cmt_mtp, crt_mn;
  double cmt_s1, cmt_d1, cmt_sl1, cmt_s2, cmt_d2, cmt_sl2;
};
*/

int
main() {
  using namespace SourceInfo;
  //   int N = 1;   // order

  //   double theta = EIGEN_PI / 2.0;   // rotation angle
  std::vector<double> vec_theta = {0.0, EIGEN_PI / 6.0, EIGEN_PI / 4.0,
                                   EIGEN_PI / 3.0, EIGEN_PI / 2.0};

  for (int idx = 0; idx < vec_theta.size(); ++idx) {
    int l = 1;   // degree
    double theta = vec_theta[idx];
    auto wigdmat =
        GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                         GSHTrans::Single, GSHTrans::ColumnMajor>(l, l, l,
                                                                  theta);
    Eigen::MatrixXcd mat_val = Eigen::MatrixXcd::Zero(2 * l + 1, 2 * l + 1);
    int cidx = 0;   // row index
    for (int N = -l; N <= l; ++N) {
      // auto wigtemp = GSHTrans::Wigner(l, l, N, theta);   // Wigner D-matrix
      auto dl = wigdmat[N];
      int ridx = 0;   // reset row index
      for (int mp = -l; mp <= l; ++mp) {
        std::complex<double> tmp =
            dl[l, mp] * std::sqrt((4.0 * EIGEN_PI) / (2 * l + 1));
        mat_val(ridx, cidx) = tmp;
        ++ridx;
      }
      ++cidx;
    }
    std::cout << mat_val << "\n\n";
  }

  //   we want to be able to get Y_{lm}^N(theta,phi). Consequently let's test
  //   this out
  auto ylmn = [](int l, int m, int N, double theta, double phi) {
    // auto wigtemp = GSHTrans::Wigner(l, l, N, theta);   // Wigner D-matrix
    auto wigtemp =
        GSHTrans::Wigner<double, GSHTrans::Ortho, GSHTrans::All, GSHTrans::All,
                         GSHTrans::Single, GSHTrans::ColumnMajor>(l, l, l,
                                                                  theta);
    auto dl = wigtemp[N];
    auto tmp = dl[l, m];
    auto ylm = tmp * std::exp(std::complex<double>(0.0, m * phi));
    return ylm;
  };

  // test
  //   for (int idx = 0; idx < vec_theta.size(); ++idx) {
  double theta = EIGEN_PI / 2.0;
  double phi = EIGEN_PI / 4.0;
  //   l = 2;
  for (int l = 0; l < 4; ++l) {
    for (int m = -l; m <= l; ++m) {
      // for (int N = -l; N <= l; ++N) {
      int N = 0;
      auto ylm = ylmn(l, m, N, theta, phi);
      std::cout << "Y_" << l << m << "(" << theta << "," << phi << ") = " << ylm
                << "\n";
    }
  }

  /////////////////////////////////////////////////
  // test out read in of CMT
  std::string filename = "examples/china_cmt_event";
  std::fstream modelfile;
  modelfile.open(filename, std::ios::in);
  std::string cmt_id;
  int cmt_yr, cmt_day, cmt_hr, cmt_min;
  double cmt_sec, cmt_lat, cmt_lon, cmt_depth, cmt_step, cmt_halfd, cmt_m0,
      cmt_mrr, cmt_mtt, cmt_mpp, cmt_mrt, cmt_mrp, cmt_mtp, crt_mn;
  double cmt_s1, cmt_d1, cmt_sl1, cmt_s2, cmt_d2, cmt_sl2;
  // getting information out of file
  if (modelfile.is_open()) {
    // get first line (title)
    // getline(modelfile, modeltitle);

    // extract information from second line and move to next line
    modelfile >> cmt_id >> cmt_yr >> cmt_day >> cmt_hr >> cmt_min >> cmt_sec >>
        cmt_lat >> cmt_lon >> cmt_depth >> cmt_step >> cmt_halfd >> cmt_m0 >>
        cmt_mrr >> cmt_mtt >> cmt_mpp >> cmt_mrt >> cmt_mrp >> cmt_mtp >>
        crt_mn >> cmt_s1 >> cmt_d1 >> cmt_sl1 >> cmt_s2 >> cmt_d2 >> cmt_sl2;
    // modelfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  // test extracting CMT information
  EarthquakeCMT cmt(filename);
  cmt.PrintCMTInfo();

  return 0;
}
