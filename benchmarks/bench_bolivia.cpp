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
  InputParameters params("bench_params/input_bench_bolivia_mf.txt");
  SRInfo sr_info(params);
  auto cmt = SourceInfo::EarthquakeCMT(params);

  // earth model
  std::string cpath = params.earth_model();
  std::string earth_model_path = params.earth_model();
  prem_norm<double> norm_class;
  auto prem = EarthModels::ModelInput(earth_model_path, norm_class, "true");

  // frequency solver parameters
  double dt = params.time_step_sec(), tout = params.t_out() / 60.0, df0 = 1.0,
         wtb = 0.05, t1 = 0, t2 = tout;
  int qex = 1, nskip = 20, num_chunks = 5, NQ = 6;
  const double twopi = 2.0 * 3.1415926535897932;
  Complex myi(0.0, 1.0);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // frequency class
  SpectraSolver::FreqFull myff(params.f1(), params.f2(), params.f11(),
                               params.f12(), params.f21(), params.f22(), dt,
                               tout, df0, wtb, t1, t2, qex, prem.TimeNorm());
  auto vec_w = myff.w();

  //////////////////////////////////////////////////////////////////////////////

  // test new functionality
  SPARSESPEC::Sparse_F_Spec spec;

  timer1.start();
  MATRIX vec_raw =
      spec.Spectra(myff, prem, cmt, params, NQ, nskip, num_chunks, sr_info);
  timer1.stop("Total time for sparse frequency spectrum");

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
  vec_raw *= norm_factor;

  double hann_w = 0.2;
  // process responses
  // don't do hann filter until you have time domain response etc
  auto vec_r2t_b = processfunctions::freq2time(vec_raw, myff);
  auto a_filt0 = processfunctions::fulltime2freq(vec_r2t_b, myff, 0.01);

  // now do the convolution for a particular source time function
  double hd = 60.0 / prem.TimeNorm();   // 60 seconds half duration
  hd = 1e-9;                            // set to zero for testing
  double dec_t = 1.628;                 // 1.628 decay don't know units
  auto pi_d = 3.1415926535897932;
  double kap_val = dec_t / (std::sqrt(pi_d) * hd);
  double st_time = 2.0 / kap_val;   // set source time

  std::cout << "Source time function parameters: \n";
  std::cout << "Half duration (s): " << hd << "\n";
  std::cout << "Kappa: " << kap_val << "\n";
  std::cout << "df: " << myff.df() << "\n";
  std::cout << "Source time (s): " << st_time * prem.TimeNorm() << "\n";

  // multiply in frequency domain
  auto a_filt_stf0 = a_filt0;
  for (int idx = 0; idx < a_filt0.cols(); ++idx) {
    auto wval = myff.w(idx);
    Complex stf_factor =
        exp(-myi * wval * st_time) *
        std::exp(-(1.0 / (4.0 * pi_d)) * std::pow(wval / kap_val, 2.0));
    // stf_factor = 1.0;
    a_filt_stf0.col(idx) *= stf_factor;
  }
  auto vec_filt_t = processfunctions::filtfreq2time(a_filt_stf0, myff, false);
  auto a_filt = processfunctions::fulltime2freq(vec_filt_t, myff, hann_w);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // testing read in of yspec
  // std::cout << "First output\n";
  std::string yspec_path = "../YSpec/output/yspec.out.bolivia.mf.1";
  YSPECREADER::DataColumns yspec_data(yspec_path);

  Eigen::MatrixXd yspec_t = Eigen::MatrixXd::Zero(3, vec_r2t_b.cols());
  for (int idx = 0; idx < yspec_data.getColumn1().size(); ++idx) {
    yspec_t(0, idx) += yspec_data.getColumn2()[idx];
    yspec_t(1, idx) += yspec_data.getColumn3()[idx];
    yspec_t(2, idx) += yspec_data.getColumn4()[idx];
  }
  auto a_filt_yspec0 = processfunctions::fulltime2freq(yspec_t, myff, 0.01);
  auto a_yspec_stf0 = a_filt_yspec0;
  for (int idx = 1; idx < a_yspec_stf0.cols(); ++idx) {
    auto wval = myff.w(idx);
    Complex stf_factor =
        exp(-myi * wval * st_time) *
        std::exp(-(1.0 / (4.0 * pi_d)) * std::pow(wval / kap_val, 2.0));
    a_yspec_stf0.col(idx) *= stf_factor;
  }
  auto vec_filt_t_yspec =
      processfunctions::filtfreq2time(a_yspec_stf0, myff, false);
  auto a_filt_yspec =
      processfunctions::fulltime2freq(vec_filt_t_yspec, myff, hann_w);

  // read in mineos output
  std::string mineos_path = "../mineos/DEMO/MYEX/Syndat_ASC_BOLIVIA/"
                            "Syndat2.2000160: 0:33:16.TLY.LHZ.ASC";
  MINEOSREADER::DataColumns mineos_data(mineos_path);
  std::string mineos_path1 = "../mineos/DEMO/MYEX/Syndat_ASC_BOLIVIA/"
                             "Syndat2.2000160: 0:33:16.TLY.LHN.ASC";
  MINEOSREADER::DataColumns mineos_data1(mineos_path1);
  std::string mineos_path2 = "../mineos/DEMO/MYEX/Syndat_ASC_BOLIVIA/"
                             "Syndat2.2000160: 0:33:16.TLY.LHE.ASC";
  MINEOSREADER::DataColumns mineos_data2(mineos_path2);

  Eigen::MatrixXd mineos_t = Eigen::MatrixXd::Zero(3, vec_r2t_b.cols());
  std::size_t maxcol;
  if (mineos_data.getColumn1().size() > vec_r2t_b.cols()) {
    maxcol = vec_r2t_b.cols();
  } else {
    maxcol = mineos_data.getColumn1().size();
  }

  for (int idx = 0; idx < maxcol; ++idx) {
    mineos_t(0, idx) += mineos_data.getColumn2()[idx] * 1e-9;
    mineos_t(1, idx) += mineos_data1.getColumn2()[idx] * 1e-9;
    mineos_t(2, idx) += mineos_data2.getColumn2()[idx] * 1e-9;
  }
  auto a_mineos_0 = processfunctions::fulltime2freq(mineos_t, myff, 0.01);
  auto a_mineos_stf0 = a_mineos_0;
  for (int idx = 1; idx < a_mineos_0.cols(); ++idx) {
    auto wval = myff.w(idx);
    Complex stf_factor =
        -1.0 / (wval * wval) * exp(-myi * wval * st_time) *
        std::exp(-(1.0 / (4.0 * pi_d)) * std::pow(wval / kap_val, 2.0));
    a_mineos_stf0.col(idx) *= stf_factor;
  }
  a_mineos_stf0 *= prem.TimeNorm() * prem.TimeNorm();   // correct units
  auto vec_filt_t_mineos =
      processfunctions::filtfreq2time(a_mineos_stf0, myff, false);
  auto a_filt_mineos =
      processfunctions::fulltime2freq(vec_filt_t_mineos, myff, hann_w);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // output
  // outputting result
  // std::cout << "Second output\n";
  std::string ptf_w = "./work/benchmarks/bench_mf_bolivia_w_50_500s.out";
  std::ofstream file_w(ptf_w);
  if (!file_w) {
    std::cerr << "Error: unable to open output file_w: " << ptf_w << "\n";
    return 1;
  }

  double nval = 1.0 / prem.TimeNorm();
  const auto &vec_w_ref = vec_w;             // alias
  const std::size_t nw = vec_w_ref.size();   // use actual vector size
  file_w.setf(std::ios::fixed);
  file_w << std::setprecision(16);

  for (std::size_t idx = 0; idx < myff.i2() + 100; ++idx) {
    // index bounds are safe since nw == vec_w.size() and a_filt has
    // compatible columns
    file_w << (vec_w_ref[idx] * nval * 1000.0 / (twopi)) << ';'
           << a_filt(0, idx).real() << ';' << a_filt(0, idx).imag() << ';'
           << std::abs(a_filt(0, idx)) << ';' << a_filt(1, idx).real() << ';'
           << a_filt(1, idx).imag() << ';' << std::abs(a_filt(1, idx)) << ';'
           << a_filt(2, idx).real() << ';' << a_filt(2, idx).imag() << ';'
           << std::abs(a_filt(2, idx)) << ";" << a_filt_yspec(0, idx).real()
           << ";" << a_filt_yspec(0, idx).imag() << ";"
           << std::abs(a_filt_yspec(0, idx)) << ";"
           << a_filt_yspec(1, idx).real() << ";" << a_filt_yspec(1, idx).imag()
           << ";" << std::abs(a_filt_yspec(1, idx)) << ";"
           << a_filt_yspec(2, idx).real() << ";" << a_filt_yspec(2, idx).imag()
           << ";" << std::abs(a_filt_yspec(2, idx)) << ";"
           << a_filt_mineos(0, idx).real() << ";"
           << a_filt_mineos(0, idx).imag() << ";"
           << std::abs(a_filt_mineos(0, idx)) << ";"
           << a_filt_mineos(1, idx).real() << ";"
           << a_filt_mineos(1, idx).imag() << ";"
           << std::abs(a_filt_mineos(1, idx)) << ";"
           << a_filt_mineos(2, idx).real() << ";"
           << a_filt_mineos(2, idx).imag() << ";"
           << std::abs(a_filt_mineos(2, idx)) << '\n';
  }
  file_w.close();

  //////////////////////////////////////////////////////////////////////////////
  // output mineos time series
  std::string ptf_t = "./work/benchmarks/bench_mf_bolivia_t_50_500s.out";
  std::ofstream file_t(ptf_t);
  if (!file_t) {
    std::cerr << "Error: unable to open output file_t: " << ptf_t << "\n";
    return 1;
  }
  file_t.setf(std::ios::fixed);
  file_t << std::setprecision(16);
  for (std::size_t idx = 0; idx < vec_filt_t_mineos.cols(); ++idx) {
    // index bounds are safe since nw == vec_w.size() and a_filt has
    // compatible columns
    auto tval = idx * myff.dt() * prem.TimeNorm();
    file_t << (idx * myff.dt() - st_time) * prem.TimeNorm() << ";"
           << vec_filt_t(0, idx) << ";" << vec_filt_t(1, idx) << ";"
           << vec_filt_t(2, idx) << ";" << vec_filt_t_yspec(0, idx) << ";"
           << vec_filt_t_yspec(1, idx) << ";" << vec_filt_t_yspec(2, idx) << ';'
           << vec_filt_t_mineos(0, idx) << ';' << vec_filt_t_mineos(1, idx)
           << ';' << vec_filt_t_mineos(2, idx)

           << '\n';
    if (tval > t2 * 3600.0) {
      break;
    }
  }
  file_t.close();

  //////////////////////////////////////////////////////////////////////////////
  // output mineos time series
  std::string ptf_t_raw =
      "./work/benchmarks/bench_mf_bolivia_t_raw_50_500s.out";
  std::ofstream file_t_raw(ptf_t_raw);
  if (!file_t_raw) {
    std::cerr << "Error: unable to open output file_t_raw: " << ptf_t_raw
              << "\n";
    return 1;
  }
  file_t_raw.setf(std::ios::fixed);
  file_t_raw << std::setprecision(16);
  for (std::size_t idx = 0; idx < vec_filt_t_mineos.cols(); ++idx) {
    // index bounds are safe since nw == vec_w.size() and a_filt has
    // compatible columns
    auto tval = idx * myff.dt() * prem.TimeNorm();
    file_t_raw << (idx * myff.dt() * prem.TimeNorm()) << ";"
               << vec_r2t_b(0, idx) << ";" << vec_r2t_b(1, idx) << ";"
               << vec_r2t_b(2, idx) << ";" << yspec_t(0, idx) << ";"
               << yspec_t(1, idx) << ";" << yspec_t(2, idx) << ';'
               << mineos_t(0, idx) << ';' << mineos_t(1, idx) << ';'
               << mineos_t(2, idx)

               << '\n';
    if (tval > t2 * 3600.0) {
      break;
    }
  }
  file_t_raw.close();

  return 0;
}
#endif
