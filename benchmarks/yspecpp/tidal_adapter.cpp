#include <DSpecM1D/All>

#include <Eigen/SparseLU>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using Complex = std::complex<double>;
using SparseMatrixC = Eigen::SparseMatrix<Complex>;

constexpr int Degree = 2;
constexpr double ForceSign = 1.0;
constexpr double DimensionlessPotentialCoefficient = 1.0;
constexpr double Pi = 3.141592653589793238462643383279502884;

struct LocalLoad {
  int element;
  int node;
  int u_global;
  int v_global;
  double radius;
  double density;
  double u_contribution;
  double v_contribution;
};

template <class Model>
Eigen::VectorXd assembleTidalLoad(Full1D::SEM &sem, const Model &model,
                                  std::vector<LocalLoad> *localLoads = nullptr) {
  const auto &mesh = sem.mesh();
  const auto abar = model.OuterRadius();
  if (!(abar > 0.0)) throw std::runtime_error("model radius must be positive");

  const auto dofs = sem.hS(Degree).rows();
  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(dofs);
  const auto quadrature =
      GaussQuad::GaussLobattoLegendreQuadrature1D<double>(mesh.NN());
  for (int element = 0; element < mesh.NE(); ++element) {
    const auto jacobian = mesh.EW(element) / 2.0;
    for (int node = 0; node < mesh.NN(); ++node) {
      const auto r = mesh.NodeRadius(element, node);
      const auto rho = sem.meshModel().Density(element, node);
      const auto weightedIntegral = jacobian * quadrature.W(node) * rho *
                                    std::pow(r, 3.0) / std::pow(abar, 2.0);
      const auto uValue = ForceSign * DimensionlessPotentialCoefficient * 2.0 *
                          weightedIntegral;
      const auto vValue = ForceSign * DimensionlessPotentialCoefficient * 6.0 *
                          weightedIntegral;
      const auto uIndex = static_cast<int>(sem.ltgS(0, element, node));
      const auto vIndex = static_cast<int>(sem.ltgS(1, element, node));
      rhs(uIndex) += uValue;
      rhs(vIndex) += vValue;
      if (localLoads) {
        localLoads->push_back(LocalLoad{element, node, uIndex, vIndex,
                                        static_cast<double>(r),
                                        static_cast<double>(rho),
                                        static_cast<double>(uValue),
                                        static_cast<double>(vValue)});
      }
    }
  }
  return rhs;
}

void writeVectorCsv(const std::string &path, const Eigen::VectorXcd &values) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("could not open " + path);
  out << "global_dof,real,imag\n" << std::setprecision(17);
  for (Eigen::Index i = 0; i < values.size(); ++i)
    out << i << ',' << values(i).real() << ',' << values(i).imag() << '\n';
}

void writeLocalLoadsCsv(const std::string &path,
                        const std::vector<LocalLoad> &loads,
                        const Eigen::VectorXd &rhs) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("could not open " + path);
  out << "element,node,u_global,v_global,radius_nondimensional,density_nondimensional,u_local,v_local,u_global_rhs,v_global_rhs\n"
      << std::setprecision(17);
  for (const auto &load : loads) {
    out << load.element << ',' << load.node << ',' << load.u_global << ','
        << load.v_global << ',' << load.radius << ',' << load.density << ','
        << load.u_contribution << ',' << load.v_contribution << ','
        << rhs(load.u_global) << ',' << rhs(load.v_global) << '\n';
  }
}

std::pair<double, double> writeN2Csv(const std::string &path,
                                     const Full1D::SEM &sem,
                                     const auto &model) {
  const auto &mesh = sem.mesh();
  const auto fluidBounds = mesh.FS_Boundaries();
  if (fluidBounds.size() < 2)
    throw std::runtime_error("PREM model does not expose an outer-core fluid shell");
  const int firstFluid = fluidBounds[0] + 1;
  const int afterFluid = fluidBounds[1] + 1;
  std::ofstream out(path);
  if (!out) throw std::runtime_error("could not open " + path);
  out << "radius_m,n_squared_s_minus2\n" << std::setprecision(17);

  double maxN2 = -std::numeric_limits<double>::infinity();
  for (int element = firstFluid; element < afterFluid; ++element) {
    for (int node = 0; node < mesh.NN(); ++node) {
      const auto r = mesh.NodeRadius(element, node);
      const auto rho = sem.meshModel().Density(element, node);
      const auto drho = model.Density(1).template Evaluate<1>(r);
      const auto gravity = sem.meshModel().Gravity(element, node);
      const auto kappa = model.Kappa(1)(r);
      const auto n2 = -drho * gravity / rho - rho * gravity * gravity / kappa;
      if (!std::isfinite(n2))
        throw std::runtime_error("non-finite outer-core N^2 sample");
      maxN2 = std::max(maxN2, n2);
      out << r * model.LengthNorm() << ','
          << n2 / (model.TimeNorm() * model.TimeNorm()) << '\n';
    }
  }
  if (!std::isfinite(maxN2) || !(maxN2 > 0.0))
    throw std::runtime_error("outer-core N^2 has no positive values");
  return {mesh.ELR(firstFluid) * model.LengthNorm(),
          mesh.EUR(afterFluid - 1) * model.LengthNorm()};
}

void runManufactured(const std::string &modelPath, double maxstep,
                     const std::string &prefix) {
  prem_norm<double> norm;
  EarthModels::ModelInput<double> model(modelPath, norm);
  constexpr int nq = 5;
  Full1D::SEM sem(model, maxstep, nq, Degree);
  std::vector<LocalLoad> localLoads;
  const auto rhs = assembleTidalLoad(sem, model, &localLoads);
  Eigen::VectorXcd complexRhs = rhs.cast<Complex>();
  writeVectorCsv(prefix + ".rhs.csv", complexRhs);
  writeLocalLoadsCsv(prefix + ".local-loads.csv", localLoads, rhs);
  std::ofstream meta(prefix + ".meta.json");
  if (!meta) throw std::runtime_error("could not open manufactured metadata");
  meta << std::setprecision(17)
       << "{\n  \"mode\": \"manufactured_load_only\",\n"
       << "  \"degree\": 2,\n  \"nq\": 5,\n"
       << "  \"mesh_elements\": " << sem.mesh().NE() << ",\n"
       << "  \"global_dofs\": " << rhs.size() << ",\n"
       << "  \"a_nondimensional\": " << model.OuterRadius() << ",\n"
       << "  \"force_sign\": 1,\n  \"c20_dimensionless\": 1,\n"
       << "  \"sum_u_rhs\": ";
  double sumU = 0.0;
  double sumV = 0.0;
  for (Eigen::Index index = 0; index < rhs.size(); ++index) {
    // The manufactured SEM map interleaves U, V, and potential DOFs; the
    // potential load is zero for this external-force assembly.
    if (index % 3 == 0) sumU += rhs(index);
    if (index % 3 == 1) sumV += rhs(index);
  }
  meta << sumU << ",\n  \"sum_v_rhs\": " << sumV << "\n}\n";
}

void runPhysical(const std::string &modelPath, double maxstep,
                 double periodHours, const std::string &prefix) {
  using Clock = std::chrono::steady_clock;
  const auto started = Clock::now();
  prem_norm<double> norm;
  EarthModels::ModelInput<double> model(modelPath, norm);
  constexpr int nq = 5;
  const auto meshStarted = Clock::now();
  Full1D::SEM sem(model, maxstep, nq, Degree);
  const auto meshSeconds =
      std::chrono::duration<double>(Clock::now() - meshStarted).count();

  std::vector<LocalLoad> localLoads;
  const auto realRhs = assembleTidalLoad(sem, model, &localLoads);
  const Eigen::VectorXcd rhs = realRhs.cast<Complex>();
  const auto timeNorm = model.TimeNorm();
  const auto omega = Complex(2.0 * Pi / (periodHours * 3600.0) * timeNorm, 0.0);
  SparseMatrixC stiffness = sem.hS(Degree).cast<Complex>();
  SparseMatrixC inertia = sem.pS(Degree).cast<Complex>();
  SparseMatrixC system = stiffness - omega * omega * inertia;
  system.makeCompressed();

  const auto solveStarted = Clock::now();
  Eigen::SparseLU<SparseMatrixC, Eigen::COLAMDOrdering<int>> solver;
  solver.compute(system);
  if (solver.info() != Eigen::Success)
    throw std::runtime_error("SparseLU factorization failed");
  const Eigen::VectorXcd solution = solver.solve(rhs);
  if (solver.info() != Eigen::Success)
    throw std::runtime_error("SparseLU solve failed");
  const auto solveSeconds =
      std::chrono::duration<double>(Clock::now() - solveStarted).count();
  const auto residual = system * solution - rhs;
  const auto relativeResidual = residual.norm() / std::max(rhs.norm(), 1e-300);
  double matrixInfinityNorm = 0.0;
  Eigen::VectorXd matrixRowSums = Eigen::VectorXd::Zero(system.rows());
  for (int outer = 0; outer < system.outerSize(); ++outer) {
    for (SparseMatrixC::InnerIterator entry(system, outer); entry; ++entry)
      matrixRowSums(entry.row()) += std::abs(entry.value());
  }
  matrixInfinityNorm = matrixRowSums.size() == 0 ? 0.0 : matrixRowSums.maxCoeff();
  const auto residualInfinityNorm = residual.size() == 0 ? 0.0 : residual.cwiseAbs().maxCoeff();
  const auto solutionInfinityNorm = solution.size() == 0 ? 0.0 : solution.cwiseAbs().maxCoeff();
  const auto rhsInfinityNorm = rhs.size() == 0 ? 0.0 : rhs.cwiseAbs().maxCoeff();
  const auto backwardDenominator = matrixInfinityNorm * solutionInfinityNorm + rhsInfinityNorm;
  const auto backwardError = backwardDenominator == 0.0
                                 ? (residualInfinityNorm == 0.0
                                        ? 0.0
                                        : std::numeric_limits<double>::infinity())
                                 : residualInfinityNorm / backwardDenominator;
  constexpr double BackwardErrorLimit = 1e-12;
  if (!solution.allFinite() || !std::isfinite(relativeResidual) ||
      !std::isfinite(backwardError))
    throw std::runtime_error("non-finite solution or residual diagnostic");
  if (backwardError > BackwardErrorLimit)
    throw std::runtime_error("normwise backward error exceeds 1e-12");

  std::ofstream radial(prefix + ".radial.csv");
  if (!radial) throw std::runtime_error("could not open radial output");
  radial << "element,node,radius_m,U_scaled_m_real,U_scaled_m_imag,V_scaled_m_real,V_scaled_m_imag\n"
         << std::setprecision(17);
  const auto lengthScale = model.LengthNorm();
  for (int element = 0; element < sem.mesh().NE(); ++element) {
    for (int node = 0; node < sem.mesh().NN(); ++node) {
      const auto r = sem.mesh().NodeRadius(element, node) * lengthScale;
      const auto u = solution(sem.ltgS(0, element, node)) * lengthScale;
      const auto v = solution(sem.ltgS(1, element, node)) * lengthScale;
      radial << element << ',' << node << ',' << r << ',' << u.real() << ','
             << u.imag() << ',' << v.real() << ',' << v.imag() << '\n';
    }
  }
  radial.close();
  writeVectorCsv(prefix + ".rhs.csv", rhs);
  const auto [icbRadiusM, cmbRadiusM] = writeN2Csv(prefix + ".n2.csv", sem, model);

  const auto maxU = [&]() {
    double value = 0.0;
    for (int element = 0; element < sem.mesh().NE(); ++element)
      for (int node = 0; node < sem.mesh().NN(); ++node)
        value = std::max(value, std::abs(solution(sem.ltgS(0, element, node))));
    return value * lengthScale;
  }();
  const auto maxV = [&]() {
    double value = 0.0;
    for (int element = 0; element < sem.mesh().NE(); ++element)
      for (int node = 0; node < sem.mesh().NN(); ++node)
        value = std::max(value, std::abs(solution(sem.ltgS(1, element, node))));
    return value * lengthScale;
  }();
  const auto elapsed = std::chrono::duration<double>(Clock::now() - started).count();
  std::ofstream meta(prefix + ".meta.json");
  if (!meta) throw std::runtime_error("could not open metadata output");
  meta << std::setprecision(17)
       << "{\n  \"mode\": \"physical_fig7_diagnostic\",\n"
       << "  \"status\": \"success\",\n"
       << "  \"degree\": 2,\n  \"nq\": 5,\n"
       << "  \"maxstep_fraction_of_planet_radius\": " << maxstep << ",\n"
       << "  \"period_hours\": " << periodHours << ",\n"
       << "  \"omega_rad_s\": " << 2.0 * Pi / (periodHours * 3600.0) << ",\n"
       << "  \"gravity\": \"full self-gravity\",\n"
       << "  \"attenuation\": false,\n"
       << "  \"model_units\": \"PREM normalized by prem_norm, physical outputs rescaled to SI\",\n"
       << "  \"external_potential\": \"cbar_20*(r/a)^2*Y_20; unit cbar_20=1\",\n"
       << "  \"force_sign\": 1,\n"
       << "  \"force_sign_note\": \"chosen explicitly; paper does not specify the Sun/Moon coefficient or sign\",\n"
       << "  \"displacement_units\": \"m per unit dimensionless cbar_20; shape only\",\n"
       << "  \"mesh_elements\": " << sem.mesh().NE() << ",\n"
       << "  \"icb_radius_m\": " << icbRadiusM << ",\n"
       << "  \"cmb_radius_m\": " << cmbRadiusM << ",\n"
       << "  \"local_gll_nodes\": " << sem.mesh().NE() * sem.mesh().NN() << ",\n"
       << "  \"global_dofs\": " << solution.size() << ",\n"
       << "  \"rhs_nonzeros\": " << (rhs.array() != Complex(0.0, 0.0)).count() << ",\n"
       << "  \"matrix_nonzeros\": " << system.nonZeros() << ",\n"
       << "  \"linear_residual_relative_l2\": " << relativeResidual << ",\n"
       << "  \"normwise_backward_error_inf\": " << backwardError << ",\n"
       << "  \"normwise_backward_error_limit\": " << BackwardErrorLimit << ",\n"
       << "  \"linear_solver\": \"Eigen SparseLU COLAMD\",\n"
       << "  \"condition_estimate\": null,\n"
       << "  \"max_abs_U_scaled_m\": " << maxU << ",\n"
       << "  \"max_abs_V_scaled_m\": " << maxV << ",\n"
       << "  \"mesh_build_s\": " << meshSeconds << ",\n"
       << "  \"factor_and_solve_s\": " << solveSeconds << ",\n"
       << "  \"total_s\": " << elapsed << "\n}\n";
}
}  // namespace

int main(int argc, char **argv) {
  try {
    if (argc == 5 && std::string(argv[1]) == "manufactured") {
      runManufactured(argv[2], std::stod(argv[3]), argv[4]);
      return 0;
    }
    if (argc == 6 && std::string(argv[1]) == "run") {
      runPhysical(argv[2], std::stod(argv[3]), std::stod(argv[4]), argv[5]);
      return 0;
    }
    std::cerr << "usage:\n"
              << "  fig7_tidal_adapter manufactured MODEL MAXSTEP OUTPUT_PREFIX\n"
              << "  fig7_tidal_adapter run MODEL MAXSTEP PERIOD_H OUTPUT_PREFIX\n";
    return 2;
  } catch (const std::exception &error) {
    std::cerr << "fig7 tidal adapter: " << error.what() << '\n';
    return 1;
  }
}
