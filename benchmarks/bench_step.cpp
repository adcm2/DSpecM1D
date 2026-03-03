#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#include <iostream>
#include <PlanetaryModel/All>
#include <DSpecM1D/Timer>
#include <DSpecM1D/All>
#include <SpectraSolver/FF>

template <typename FLOAT> class prem_norm {
public:
  FLOAT LengthNorm() const { return _length_norm; };
  FLOAT MassNorm() const { return _mass_norm; };
  FLOAT TimeNorm() const { return _time_norm; };

  prem_norm() {};

private:
  FLOAT _length_norm = 1000.0;
  FLOAT _mass_norm = 5515 * std::pow(_length_norm, 3.0);
  FLOAT _time_norm = 1 / std::sqrt(3.1415926535897932 * 6.67230e-11 * 5515);
};

int
main() {
  using Complex = std::complex<double>;
  using MATRIX = Eigen::MatrixXcd;
  using SMATRIX = Eigen::SparseMatrix<Complex>;
  using SLU = Eigen::SparseLU<SMATRIX, Eigen::COLAMDOrdering<int>>;
  Timer timer1;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Read all parameters from the input file_w
  InputParameters params("./bench_params/input_bench_step.txt");

  // earth model
  std::string cpath = params.earth_model();
  std::string earth_model_path = params.earth_model();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // parameters of sem
  int lval = params.lmax(), NQ = 6;
  bool toaug = false;
  double maxstep = 0.05;
  // std::vector<double> vec_step{0.1, 0.05, 0.01, 0.001};

  // std::cout << "Max step: \n";
  // std::cin >> maxstep;
  const double twopi = 2.0 * 3.1415926535897932;
  std::string pathpert;

  //////////////////////////////
  // frequency solver parameters
  double dt = params.time_step_sec(), tout = params.t_out() / 60.0, df0 = 1.0,
         wtb = 0.05, t1 = 0, t2 = tout;
  int qex = 1;
  Complex myi(0.0, 1.0);
  double droptol = 1e-4;
  auto tend = 0.6 * params.t_out();   // length of time series in seconds
  // std::cout << "tend: " << tend << "\n";
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // prem
  prem_norm<double> norm_class;
  auto prem = EarthModels::ModelInput(earth_model_path, norm_class, "true");

  // frequency class
  SpectraSolver::FreqFull myff(params.f1(), params.f2(), params.f11(),
                               params.f12(), params.f21(), params.f22(), dt,
                               tout, df0, wtb, t1, t2, qex, prem.TimeNorm());
  auto vec_w = myff.w();

  // source information
  auto cmt = SourceInfo::EarthquakeCMT(params);

  int nsteps = 50;
  // auto lambda_min = 0.63;
  auto step_0 = 2.0 * 0.63 / myff.f22();
  std::cout << "Initial step: " << step_0 << "\n";
  std::vector<double> vec_step;
  for (int idx = 0; idx < nsteps; ++idx) {
    vec_step.push_back(
        step_0 / std::pow(10.0, 2.0 * idx / static_cast<double>(nsteps - 1)));
  }

  //////////////////////////////////////////////////////////////////////////////
  // normalise
  double norm_factor;
  auto accel_norm = prem.LengthNorm() / (prem.TimeNorm() * prem.TimeNorm());
  if (params.output_type() == 0) {
    norm_factor = prem.LengthNorm();
  } else if (params.output_type() == 1) {
    norm_factor = prem.LengthNorm() / prem.TimeNorm();
  } else if (params.output_type() == 2) {
    norm_factor = accel_norm;
  }
  double hann_w = 0.2;

  // test new functionality
  SPARSESPEC::Sparse_F_Spec mytest;
  std::vector<Eigen::MatrixXcd> vec_final_w;
  std::vector<Eigen::MatrixXd> vec_final_t;
  // omp_set_num_threads(10);
  for (int idx = 0; idx < nsteps; ++idx) {
    maxstep = vec_step[idx];
    int nskip = 3 * std::floor(maxstep / ((vec_w[1] - vec_w[0]) * 0.003));
    nskip += 1;
    // for this maxstep set up the SEM
    Full1D::specsem sem(prem, maxstep, NQ, lval);
    std::cout << "\nDoing step: " << maxstep << ", nskip: " << nskip << "\n";
    MATRIX vec_raw = mytest.Spectra(myff, sem, prem, cmt, params, nskip);

    vec_raw *= norm_factor;

    // process responses
    auto vec_r2t_b = processfunctions::freq2time(vec_raw, myff);
    auto a_filt0 = processfunctions::fulltime2freq(vec_r2t_b, myff, 0.05);
    auto vec_filt_t = processfunctions::filtfreq2time(a_filt0, myff, false);
    auto a_filt = processfunctions::fulltime2freq(vec_filt_t, myff, hann_w);

    vec_final_w.push_back(a_filt);
    vec_final_t.push_back(vec_filt_t);
  }

  //////////////////////////////////////////////////////////////////////////////
  // find error
  std::vector<std::vector<double>> vec_err_t(vec_step.size() - 1,
                                             std::vector<double>(3, 0.0)),
      vec_err_w(vec_step.size() - 1, std::vector<double>(3, 0.0));
  std::vector<std::vector<double>> vec_l2_err_t(vec_step.size() - 1,
                                                std::vector<double>(3, 0.0)),
      vec_l2_err_w(vec_step.size() - 1, std::vector<double>(3, 0.0));
  int idxout = vec_final_t[0].cols() - 1;
  for (int idx = 0; idx < vec_final_t[0].cols(); ++idx) {
    auto tval = idx * myff.dt() * prem.TimeNorm();
    if (tval > t2 * 3600.0) {
      idxout = idx;
      break;
    }
  }

  for (int idx = 0; idx < nsteps; ++idx) {
    vec_final_t[idx]
        .block(0, idxout, 3, vec_final_t[idx].cols() - idxout)
        .setZero();
  }

  // "exact" solution is the one with smallest step
  auto idxback = nsteps - 1;
  auto vec_t_ex = vec_final_t[idxback];
  auto vec_w_ex = vec_final_w[idxback];

  // infinity and l2 norms
  auto tmp_t = 100.0 / idxout;
  auto mv_t_1 = tmp_t / vec_t_ex.row(0).lpNorm<Eigen::Infinity>();
  auto mv_t_2 = tmp_t / vec_t_ex.row(1).lpNorm<Eigen::Infinity>();
  auto mv_t_3 = tmp_t / vec_t_ex.row(2).lpNorm<Eigen::Infinity>();

  auto mv_t_1_l2 = 100.0 / vec_t_ex.row(0).norm();
  auto mv_t_2_l2 = 100.0 / vec_t_ex.row(1).norm();
  auto mv_t_3_l2 = 100.0 / vec_t_ex.row(2).norm();

  auto wcols = vec_final_w[nsteps - 1].cols();
  auto tmp_w = 100.0 / wcols;
  auto mv_w_1 = tmp_w / vec_w_ex.row(0).lpNorm<Eigen::Infinity>();
  auto mv_w_2 = tmp_w / vec_w_ex.row(1).lpNorm<Eigen::Infinity>();
  auto mv_w_3 = tmp_w / vec_w_ex.row(2).lpNorm<Eigen::Infinity>();
  auto mv_w_1_l2 = 100.0 / vec_w_ex.row(0).norm();
  auto mv_w_2_l2 = 100.0 / vec_w_ex.row(1).norm();
  auto mv_w_3_l2 = 100.0 / vec_w_ex.row(2).norm();

  for (int idx = 0; idx < vec_step.size() - 1; ++idx) {

    // responses at idx
    auto vec_t = vec_final_t[idx];
    auto w_res = vec_final_w[idx];

    // error vectors
    auto vec_t_err = vec_t - vec_t_ex;
    auto vec_w_err = w_res - vec_w_ex;

    // infinity norms
    vec_err_t[idx][0] = vec_t_err.row(0).lpNorm<1>() * mv_t_1;
    vec_err_t[idx][1] = vec_t_err.row(1).lpNorm<1>() * mv_t_2;
    vec_err_t[idx][2] = vec_t_err.row(2).lpNorm<1>() * mv_t_3;
    vec_err_w[idx][0] = vec_w_err.row(0).lpNorm<1>() * mv_w_1;
    vec_err_w[idx][1] = vec_w_err.row(1).lpNorm<1>() * mv_w_2;
    vec_err_w[idx][2] = vec_w_err.row(2).lpNorm<1>() * mv_w_3;

    // l2 norms
    vec_l2_err_t[idx][0] = vec_t_err.row(0).norm() * mv_t_1_l2;
    vec_l2_err_t[idx][1] = vec_t_err.row(1).norm() * mv_t_2_l2;
    vec_l2_err_t[idx][2] = vec_t_err.row(2).norm() * mv_t_3_l2;
    vec_l2_err_w[idx][0] = vec_w_err.row(0).norm() * mv_w_1_l2;
    vec_l2_err_w[idx][1] = vec_w_err.row(1).norm() * mv_w_2_l2;
    vec_l2_err_w[idx][2] = vec_w_err.row(2).norm() * mv_w_3_l2;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // output
  // outputting result
  // std::cout << "Second output\n";
  std::string ptf_w = "./work/benchmarks/bench_w_NQ";
  ptf_w += std::to_string(NQ) + "_step.out";
  std::ofstream file_w(ptf_w);
  if (!file_w) {
    std::cerr << "Error: unable to open output file_w: " << ptf_w << "\n";
    return 1;
  }

  double nval = 1.0 / prem.TimeNorm();
  const auto &vec_w_ref = vec_w;             // alias
  const std::size_t nw = vec_w_ref.size();   // use actual vector size

  // Write header (optional)
  // file_w <<
  // "#freq_mHz;Z_re;Z_im;Z_abs;TH_re;TH_im;TH_abs;PH_re;PH_im;PH_abs\n";
  file_w.setf(std::ios::fixed);
  file_w << std::setprecision(22);

  for (std::size_t idx = 0; idx < myff.i2() + 100; ++idx) {
    // index bounds are safe since nw == vec_w.size() and a_filt has
    // compatible columns
    file_w << (vec_w_ref[idx] * nval * 1000.0 / (twopi));
    for (int idx2 = 0; idx2 < vec_step.size(); ++idx2) {
      file_w << ";" << vec_final_w[idx2](0, idx).real() << ';'
             << vec_final_w[idx2](0, idx).imag() << ';'
             << std::abs(vec_final_w[idx2](0, idx)) << ';'
             << vec_final_w[idx2](1, idx).real() << ';'
             << vec_final_w[idx2](1, idx).imag() << ';'
             << std::abs(vec_final_w[idx2](1, idx)) << ';'
             << vec_final_w[idx2](2, idx).real() << ';'
             << vec_final_w[idx2](2, idx).imag() << ';'
             << std::abs(vec_final_w[idx2](2, idx));
    }
    file_w << '\n';
  }
  file_w.close();

  //////////////////////////////////////////////////////////////////////////////
  // output mineos time series
  std::string ptf_t = "./work/benchmarks/bench_t_NQ";
  ptf_t += std::to_string(NQ) + "_step.out";
  std::ofstream file_t(ptf_t);
  if (!file_t) {
    std::cerr << "Error: unable to open output file_t: " << ptf_t << "\n";
    return 1;
  }
  file_t.setf(std::ios::fixed);
  file_t << std::setprecision(22);
  for (std::size_t idx = 0; idx < idxout; ++idx) {
    // index bounds are safe since nw == vec_w.size() and a_filt has
    // compatible columns
    auto tval = idx * myff.dt() * prem.TimeNorm();
    file_t << (idx * myff.dt()) * prem.TimeNorm();
    for (int idx2 = 0; idx2 < nsteps; ++idx2) {
      file_t << ";" << vec_final_t[idx2](0, idx) << ";"
             << vec_final_t[idx2](1, idx) << ";" << vec_final_t[idx2](2, idx);
    }
    file_t << '\n';
  }
  file_t.close();

  //////////////////////////////////////////////////////////////////////////////
  // output mineos time series
  std::string ptf_err = "./work/benchmarks/bench_NQ";
  ptf_err += std::to_string(NQ) + "_step_error_";
  ptf_err += std::to_string(int(params.f22())) + ".out";
  std::ofstream file_err(ptf_err);
  if (!file_err) {
    std::cerr << "Error: unable to open output file_t: " << ptf_err << "\n";
    return 1;
  }
  file_err.setf(std::ios::fixed);
  file_err << std::setprecision(22);
  for (int idx = 0; idx < vec_step.size() - 1; ++idx) {
    file_err << vec_step[idx] << ";" << vec_err_t[idx][0] << ";"
             << vec_err_t[idx][1] << ";" << vec_err_t[idx][2] << ";"
             << vec_err_w[idx][0] << ";" << vec_err_w[idx][1] << ";"
             << vec_err_w[idx][2] << "\n";
  }
  file_err.close();

  //////////////////////////////////////////////////////////////////////////////
  // output mineos time series
  std::string ptf_l2 = "./work/benchmarks/bench_NQ";
  ptf_l2 += std::to_string(NQ) + "_step_error_l2_";
  ptf_l2 += std::to_string(int(params.f22())) + ".out";
  std::ofstream file_l2(ptf_l2);
  if (!file_l2) {
    std::cerr << "Error: unable to open output file_t: " << ptf_l2 << "\n";
    return 1;
  }
  file_l2.setf(std::ios::fixed);
  file_l2 << std::setprecision(22);
  for (int idx = 0; idx < vec_step.size() - 1; ++idx) {
    file_l2 << vec_step[idx] << ";" << vec_l2_err_t[idx][0] << ";"
            << vec_l2_err_t[idx][1] << ";" << vec_l2_err_t[idx][2] << ";"
            << vec_l2_err_w[idx][0] << ";" << vec_l2_err_w[idx][1] << ";"
            << vec_l2_err_w[idx][2] << "\n";
  }
  file_l2.close();

  return 0;
}
#endif
