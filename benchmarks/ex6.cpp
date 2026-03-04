#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

// Standard library includes
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <complex>
#include <cmath>
#include <vector>

// Project-specific includes
#include <PlanetaryModel/All>
#include <DSpecM1D/Timer>
#include <DSpecM1D/All>
#include <SpectraSolver/FF>

constexpr double PI = 3.1415926535897932;
constexpr double TWO_PI = 2.0 * PI;

template <typename FLOAT> class prem_norm {
public:
  prem_norm() = default;

  FLOAT LengthNorm() const { return _length_norm; }
  FLOAT MassNorm() const { return _mass_norm; }
  FLOAT TimeNorm() const { return _time_norm; }

private:
  FLOAT _length_norm = 1000.0;
  FLOAT _mass_norm = 5515.0 * std::pow(_length_norm, 3.0);
  FLOAT _time_norm = 1.0 / std::sqrt(PI * 6.67230e-11 * 5515.0);
};

int
main() {
  using Complex = std::complex<double>;
  using SMATRIX = Eigen::SparseMatrix<Complex>;

  Timer timer1;

  // --- 1. Read Inputs & Earth Model ---
  InputParameters params("bench_params/ex6.txt");
  std::string earth_model_path = params.earth_model();
  auto cmt = SourceInfo::EarthquakeCMT(params);

  // --- 2. Parameters of SEM ---
  int lval = params.lmax();
  int NQ = 5;
  double maxstep = 0.05;

  std::cout << "Enter max step size for SEM integration: \n";
  std::cin >> maxstep;

  // Initialize Earth model
  timer1.start();
  prem_norm<double> norm_class;
  auto prem = EarthModels::ModelInput(earth_model_path, norm_class, "true");

  // Initialize SEM
  Full1D::specsem sem(prem, maxstep, NQ, lval);
  auto _mesh = sem.mesh();
  auto meshmodel = sem.mesh_model();

  // --- 3. Calculate the Brunt-Väisälä frequency ---
  std::vector<std::vector<double>> vec_N2;

  // Indices of fluid-solid boundaries
  auto fsb = _mesh.FS_Boundaries();
  int idxlow = fsb[0] + 1;
  int idxup = fsb[1] + 1;

  double maxn2 = 0.0;
  double minn2 = 0.0;

  // Calculate N^2 in the fluid region
  for (int idxe = idxlow; idxe < idxup; ++idxe) {
    std::vector<double> tmp;
    tmp.reserve(_mesh.NN());
    for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
      auto crad = _mesh.NodeRadius(idxe, idxn);
      auto rho = meshmodel.Density(idxe, idxn);
      auto drho_dr = prem.Density(1).Derivative(crad);
      auto gval = meshmodel.Gravity(idxe, idxn);

      auto tn2 =
          -drho_dr * gval / rho - rho * gval * gval / prem.Kappa(1)(crad);
      tmp.push_back(tn2);

      if (tn2 > maxn2)
        maxn2 = tn2;
      if (tn2 < minn2)
        minn2 = tn2;
    }
    vec_N2.push_back(tmp);
  }

  // BV frequency conversions
  auto max_imag = std::sqrt(-minn2) / (TWO_PI * prem.TimeNorm());
  auto N_hours = max_imag * 3600.0;   // convert to cycles per hour
  auto n_hperiod = 1.0 / N_hours;

  // --- 4. Frequency Setup ---
  std::vector<double> vec_tperiods{12.0, 48.0, 168.0, 1680.0};   // in hours
  std::vector<Complex> vec_omega;
  vec_omega.reserve(vec_tperiods.size());

  for (double period : vec_tperiods) {
    auto freq = 1.0 / (period * 3600.0);   // frequency in Hz
    Complex omega_c = TWO_PI * freq * prem.TimeNorm();
    vec_omega.push_back(omega_c);
  }

  // --- 5. System Setup (Tidal Forcing) ---
  auto idxl = 3;   // l = 2 or 3 for tidal forcing
  SMATRIX ke_s = sem.H_S(idxl).cast<Complex>();
  SMATRIX in_s = sem.P_S(idxl).cast<Complex>();

  Eigen::VectorXcd rhs = Eigen::VectorXcd::Zero(ke_s.rows());

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
  }

  // --- 6. Solve Linear System ---
  timer1.start();
  Eigen::SparseLU<SMATRIX, Eigen::COLAMDOrdering<int>> solver;
  std::vector<Eigen::VectorXcd> vec_solutions;
  vec_solutions.reserve(vec_omega.size());

  for (std::size_t idx = 0; idx < vec_omega.size(); ++idx) {
    Complex omega_c = vec_omega[idx];
    SMATRIX A = -omega_c * omega_c * in_s + ke_s;
    A.makeCompressed();
    solver.compute(A);

    if (solver.info() != Eigen::Success) {
      std::cerr << "Decomposition failed for frequency index: " << idx << "\n";
      return 1;
    }
    vec_solutions.push_back(solver.solve(rhs));
  }
  timer1.stop("Total time for solving linear system");

  // --- 7. Calculate Kinetic Energy ---
  std::vector<double> vec_total_KE, vec_fluid_KE;
  vec_total_KE.reserve(vec_omega.size());
  vec_fluid_KE.reserve(vec_omega.size());

  auto energynorm = prem.MassNorm() * prem.LengthNorm() * prem.LengthNorm() /
                    (prem.TimeNorm() * prem.TimeNorm());

  for (std::size_t idxf = 0; idxf < vec_omega.size(); ++idxf) {
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

        if (idxe >= idxlow && idxe < idxup) {
          fluid_ke += ke_local;
        }
      }
    }
    vec_total_KE.push_back(total_ke * energynorm);
    vec_fluid_KE.push_back(fluid_ke * energynorm);
  }

  // --- 8. Outputs ---

  // 8a. Output Brunt-Väisälä frequency
  std::string N2_out_path =
      "./plotting/outputs/ex6_N2_" + std::to_string(maxstep) + ".out";
  std::ofstream N2_file(N2_out_path);
  if (!N2_file) {
    std::cerr << "Error: unable to open output file: " << N2_out_path << "\n";
    return 1;
  }

  N2_file << std::fixed << std::setprecision(22);
  auto fnorm2 = 1.0 / (prem.TimeNorm() * prem.TimeNorm());

  for (int idxe = idxlow; idxe < idxup; ++idxe) {
    for (int idxn = 0; idxn < _mesh.NN(); ++idxn) {
      auto crad = _mesh.NodeRadius(idxe, idxn) * prem.LengthNorm();
      N2_file << crad << ";" << vec_N2[idxe - idxlow][idxn] * fnorm2 << "\n";
    }
  }
  N2_file.close();

  // 8b. Output Frequencies
  std::string freq_out_path =
      "./plotting/outputs/ex6_w_" + std::to_string(maxstep) + ".out";
  std::ofstream freq_file(freq_out_path);
  if (!freq_file) {
    std::cerr << "Error: unable to open output file: " << freq_out_path << "\n";
    return 1;
  }

  freq_file << std::fixed << std::setprecision(22);
  for (std::size_t idx = 0; idx < vec_omega.size(); ++idx) {
    freq_file << vec_tperiods[idx] / n_hperiod << ";"
              << vec_fluid_KE[idx] / vec_total_KE[idx] << "\n";
  }
  freq_file.close();

  // 8c. Output Radial Solution
  std::string pathtofile = "./plotting/outputs/ex6_radial_response_" +
                           std::to_string(maxstep) + ".out";
  std::ofstream file(pathtofile);
  if (!file) {
    std::cerr << "Error: unable to open output file: " << pathtofile << "\n";
    return 1;
  }

  file << std::fixed << std::setprecision(22);
  for (int idxe = 0; idxe < _mesh.NE(); ++idxe) {
    for (int idxq = 0; idxq < NQ; ++idxq) {
      auto crad = _mesh.NodeRadius(idxe, idxq) * prem.LengthNorm();
      auto uidx = sem.LtG_S(0, idxe, idxq);
      auto vidx = sem.LtG_S(1, idxe, idxq);

      file << crad;
      for (std::size_t idxf = 0; idxf < vec_omega.size(); ++idxf) {
        auto sol = vec_solutions[idxf];
        file << ";" << sol(uidx).real() << ";" << sol(uidx).imag() << ";"
             << std::abs(sol(uidx)) << ";" << sol(vidx).real() << ";"
             << sol(vidx).imag() << ";" << std::abs(sol(vidx));
      }
      file << '\n';
    }
  }
  file.close();

  return 0;
}