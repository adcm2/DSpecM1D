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

  //////////////////////////////////////////////////////////////////////////////
  // Read all parameters from the input file
  InputParameters params("bench_params/input_bench_hf.txt");
  SRInfo sr_info(params);
  // Earth model
  std::string earth_model_path = params.earth_model();

  //////////////////////////////////////////////////////////////////////////////
  // Parameters of SEM
  int lval = params.lmax(), NQ = 5;
  double maxstep = 0.05;
  std::cout << "Enter max step size for SEM integration: \n";
  std::cin >> maxstep;
  const double twopi = 2.0 * 3.1415926535897932;
  auto cmt = SourceInfo::EarthquakeCMT(params);
  // Frequency solver parameters
  //   double dt = params.time_step_sec(), tout = params.t_out() / 60.0, df0
  //   = 1.0,
  //          wtb = 0.05, t1 = 0, t2 = tout + 1.0;
  //   int qex = 1;

  // Initialize Earth model and frequency grid
  timer1.start();
  prem_norm<double> norm_class;
  auto prem = EarthModels::ModelInput(earth_model_path, norm_class, "true");

  ////////////////////////////////
  // Initialize SEM
  timer1.start();
  Full1D::specsem sem(prem, maxstep, NQ, lval);
  auto _mesh = sem.mesh();
  auto meshmodel = sem.mesh_model();
  timer1.stop("Total time for setting up SEM class");

  ////////////////////////////////
  // calculate the Brunt-Vaisala frequency
  std::vector<std::vector<double>> vec_N2;

  // indices of fluid-solid boundaries
  auto fsb = _mesh.FS_Boundaries();
  int idxlow = fsb[0] + 1;
  int idxup = fsb[1] + 1;

  double maxn2 = 0.0, minn2 = 0.0;
  // calculate N2 in the fluid region
  for (int idxe = idxlow; idxe < idxup; ++idxe) {
    std::vector<double> tmp;
    for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
      auto crad = _mesh.NodeRadius(idxe, idxn);
      auto rho = meshmodel.Density(idxe, idxn);
      auto drho_dr = prem.Density(1).Derivative(crad);
      auto gval = meshmodel.Gravity(idxe, idxn);
      auto tn2 =
          -drho_dr * gval / rho - rho * gval * gval / prem.Kappa(1)(crad);
      tmp.push_back(tn2);
      if (tn2 > maxn2) {
        maxn2 = tn2;
      }
      if (tn2 < minn2) {
        minn2 = tn2;
      }
    }
    vec_N2.push_back(tmp);
  }
  timer1.stop("Total time for reading PREM and setting up frequency class");

  auto n_hz = std::sqrt(maxn2) / (twopi * prem.TimeNorm());
  auto n_mhz = n_hz * 1e3;
  std::cout << "Max Brunt-Vaisala frequency (mHz): " << n_mhz << "\n";

  // max imaginary
  auto max_imag = std::sqrt(-minn2) / (twopi * prem.TimeNorm());
  auto max_imag_mhz = max_imag * 1e3;
  std::cout << "Max imaginary Brunt-Vaisala frequency (mHz): " << max_imag_mhz
            << "\n";
  // BV frequency in days
  auto N_hours = max_imag * 3600.0;   // convert to cycles per hour
  std::cout << "Max Brunt-Vaisala frequency (cycles per hour): " << N_hours
            << "\n";
  auto n_hperiod = 1.0 / N_hours;

  //////////////////////////////////////////////////////////////////////////////
  // frequency setup
  int num_freqs = 4;
  // std::vector<double> vec_tperiods;
  // for (int idx = 0; idx < num_freqs; ++idx) {
  //   double period = 0.1 * std::exp(static_cast<double>(idx) * 3.0) *
  //   n_hperiod; vec_tperiods.push_back(period);
  // }
  std::vector<double> vec_tperiods{12.0, 48.0, 168.0, 1680.0};   // in hours
  std::vector<Complex> vec_omega;
  for (int idx = 0; idx < vec_tperiods.size(); ++idx) {
    auto period = vec_tperiods[idx];
    auto freq = 1.0 / (period * 3600.0);   // frequency in Hz
    Complex omega_c = twopi * freq * prem.TimeNorm();
    vec_omega.push_back(omega_c);
  }

  ////////////////////////////////
  auto idxl = 3;   // l = 2 or 3 for tidal forcing
  SMATRIX ke_s = sem.H_S(idxl).cast<Complex>();
  SMATRIX in_s = sem.P_S(idxl).cast<Complex>();
  Eigen::MatrixXcd fval = sem.CalculateForce(cmt, idxl);

  ////////////////////////////////
  std::cout << "Setting up tidal forcing term...\n";
  timer1.start();
  // tidal forcing term (maybe with correct coefficients)
  Eigen::VectorXcd rhs = Eigen::VectorXcd::Zero(ke_s.rows());

  {

    for (int idx = 0; idx < _mesh.NE(); ++idx) {
      auto jacval = 2.0 / (_mesh.EUR(idx) - _mesh.ELR(idx));
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
        auto crad = _mesh.NodeRadius(idx, idxn);
        auto idx_u = sem.LtG_S(0, idx, idxn);
        auto idx_v = sem.LtG_S(1, idx, idxn);
        auto tmp = _mesh.GLL().W(idxn) * jacval * std::pow(crad, idxl + 1);
        rhs(idx_u) += idxl * tmp;
        rhs(idx_v) += idxl * (idxl + 1.0) * tmp * std::sqrt(2.0);
      }
    };
  };

  for (int idxe = idxlow; idxe < idxup; ++idxe) {
    for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
      auto idx_u = sem.LtG_S(0, idxe, idxn);
      auto idx_v = sem.LtG_S(1, idxe, idxn);
      //   rhs(idx_u) *= 0.0;
      //   rhs(idx_v) *= 0.0;
    }
  }
  // rhs = fval.col(0);

  timer1.stop("Total time for setting up tidal forcing term");

  ////////////////////////////////
  timer1.start();
  Eigen::SparseLU<SMATRIX, Eigen::COLAMDOrdering<int>> solver;
  std::vector<Eigen::VectorXcd> vec_solutions;
  for (int idx = 0; idx < vec_omega.size(); ++idx) {
    auto omega_c = vec_omega[idx];
    SMATRIX A = -omega_c * omega_c * in_s + ke_s;
    A.makeCompressed();
    solver.compute(A);
    if (solver.info() != Eigen::Success) {
      std::cerr << "Decomposition failed for frequency index: " << idx << "\n";
      return 1;
    }
    Eigen::VectorXcd sol = solver.solve(rhs);
    vec_solutions.push_back(sol);
  }
  timer1.stop("Total time for solving linear system");

  //////////////////////////////////////////////////////////////////////////////
  // calculate the total kinetic energy in the fluid region
  std::vector<double> vec_total_KE, vec_fluid_KE;
  auto energynorm = prem.MassNorm() * prem.LengthNorm() * prem.LengthNorm() /
                    (prem.TimeNorm() * prem.TimeNorm());
  for (int idxf = 0; idxf < vec_omega.size(); ++idxf) {
    auto sol = vec_solutions[idxf];
    double total_ke = 0.0;
    double fluid_ke = 0.0;
    for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
      for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
        auto crad = _mesh.NodeRadius(idxe, idxn);
        auto rho = meshmodel.Density(idxe, idxn);
        auto jacval = 2.0 / (_mesh.EUR(idxe) - _mesh.ELR(idxe));
        auto wval = _mesh.GLL().W(idxn);
        auto uidx = sem.LtG_S(0, idxe, idxn);
        auto vidx = sem.LtG_S(1, idxe, idxn);
        auto uval = sol(uidx);
        auto vval = sol(vidx);
        auto ke_local = jacval * wval * rho *
                        (std::norm(uval) + std::norm(vval)) *
                        std::pow(crad, 2.0);
        total_ke += ke_local;
        // check if in fluid region
        if ((idxe >= idxlow) && (idxe < idxup)) {
          fluid_ke += ke_local;
        }
      }
    }
    vec_total_KE.push_back(total_ke * energynorm);
    vec_fluid_KE.push_back(fluid_ke * energynorm);
  }
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  // output Brunt-Vaisala frequency
  std::string N2_out_path = "./work/benchmarks/tidal_N2_";
  N2_out_path += std::to_string(maxstep) + ".out";
  std::ofstream N2_file(N2_out_path);
  if (!N2_file) {
    std::cerr << "Error: unable to open output file: " << N2_out_path << "\n";
    return 1;
  }
  N2_file.setf(std::ios::fixed);
  N2_file << std::setprecision(22);
  auto fnorm2 = 1.0 / (prem.TimeNorm() * prem.TimeNorm());
  for (int idxe = idxlow; idxe < idxup; ++idxe) {
    for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
      auto crad = _mesh.NodeRadius(idxe, idxn) * prem.LengthNorm();
      N2_file << crad << ";" << vec_N2[idxe - idxlow][idxn] * fnorm2 << "\n";
    }
  }
  N2_file.close();

  //////////////////////////////////////////////////////////////////////////////

  // output frequencies
  std::string freq_out_path = "./work/benchmarks/tidal_frequency_";
  freq_out_path += std::to_string(maxstep) + ".out";
  std::ofstream freq_file(freq_out_path);
  if (!freq_file) {
    std::cerr << "Error: unable to open output file: " << freq_out_path << "\n";
    return 1;
  }
  freq_file.setf(std::ios::fixed);
  freq_file << std::setprecision(22);

  for (int idx = 0; idx < vec_omega.size(); ++idx) {
    freq_file << std::setprecision(22) << vec_tperiods[idx] / n_hperiod << ";"
              << vec_fluid_KE[idx] / vec_total_KE[idx] << "\n";
  }
  freq_file.close();

  //////////////////////////////////////////////////////////////////////////////
  // Output radial solution
  std::string pathtofile = "./work/benchmarks/tidal_radial_response";
  pathtofile += std::to_string(maxstep) + ".out";
  std::ofstream file(pathtofile);
  if (!file) {
    std::cerr << "Error: unable to open output file: " << pathtofile << "\n";
    return 1;
  }

  //   double nval = 1.0 / prem.TimeNorm();
  //   const auto &vec_w_ref = vec_w;
  {
    file.setf(std::ios::fixed);
    file << std::setprecision(22);
    auto _mesh = sem.mesh();
    for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
      int maxq = NQ;
      for (int idxq = 0; idxq < maxq; ++idxq) {
        auto crad = _mesh.NodeRadius(idxe, idxq) * prem.LengthNorm();
        auto uidx = sem.LtG_S(0, idxe, idxq);
        auto vidx = sem.LtG_S(1, idxe, idxq);

        // outputting to file
        file << crad;
        for (int idxf = 0; idxf < vec_omega.size(); ++idxf) {
          auto sol = vec_solutions[idxf];
          file << ";" << sol(uidx).real() << ";" << sol(uidx).imag() << ";"
               << std::abs(sol(uidx)) << ";" << sol(vidx).real() << ";"
               << sol(vidx).imag() << ";" << std::abs(sol(vidx));
        }
        file << '\n';
      }
    }
    file.close();
  }

  return 0;
}
#endif