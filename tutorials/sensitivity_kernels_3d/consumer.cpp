#include <DSpecM1D/All>
#include <DSpecM1D/ModelInput>
#include <DSpecM1D/Windowed3DKernels>
#include <GaussQuad/All>
#include <DSpecM1D/src/InputParser.h>
#include <DSpecM1D/src/NormClass.h>
#include <DSpecM1D/src/SourceInfo.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using Complex = std::complex<double>;
constexpr int nt = 128;
constexpr double r0 = 4200000.0;
constexpr double r1 = 6200000.0;
constexpr double sourceLatitude = 10.0;
constexpr double sourceLongitude = 20.0;
constexpr double receiverLatitude = 45.0;
constexpr double receiverLongitude = 90.0;
constexpr std::array<const char *, 5> moduli{"A", "C", "F", "L", "N"};
constexpr std::array<double, 5> baselineModuliPa{
    352960000000.0, 352960000000.0, 129602500000.0,
    111678750000.0, 111678750000.0};

struct AngularRule {
  GaussQuad::Quadrature1D<double> polar;
  int nphi = 0;
};

AngularRule angularRule(int lmax) {
  const int n = std::max(8, 2 * lmax + 2);
  return {GaussQuad::GaussLegendreQuadrature1D<double>(n),
          std::max(24, 4 * lmax + 4)};
}

std::string paramsText(const fs::path &model, const fs::path &out, int lmax) {
  return "\"" + (out / "unused").string() + "\"\n\"" + model.string() +
      "\"\n4\n0\n2\n0\n0\n1e-5\n0\n" + std::to_string(lmax) +
      "\n0.1\n0.8\n30\n20\n0\n0\n0\n0\n50\n10\n20\n1e19\n2e18\n-3e18\n-5e18\n4e18\n-5e18\n3\n1\n45 90\n";
}

Eigen::VectorXd hannWindow() {
  Eigen::VectorXd window = Eigen::VectorXd::Zero(nt);
  for (int n = 0; n <= 90; ++n)
    window[n] = 0.5 - 0.5 * std::cos(2.0 * std::numbers::pi * n / 90.0);
  return window;
}

Eigen::MatrixXcd objectiveWeights() {
  Eigen::MatrixXcd q = Eigen::MatrixXcd::Zero(3, nt / 2 + 1);
  q(0, 1) = Complex(0.7, 0.2);
  q(0, 2) = Complex(-0.25, 0.45);
  return q;
}

double angularPerturbation(double theta, double phi) {
  const double x = std::cos(theta), s = std::sin(theta);
  const double n21 = std::sqrt(5.0 / (24.0 * std::numbers::pi));
  const double n32 = std::sqrt(7.0 / (480.0 * std::numbers::pi));
  return -3.0 * x * s * n21 * std::cos(phi) +
      7.5 * x * (1.0 - x * x) * n32 * std::cos(2.0 * phi);
}

struct SectionPoint {
  std::string plane;
  double distanceMetres;
  double radiusMetres;
  double theta;
  double phi;
};

struct RadialNode {
  int element;
  int node;
  double radiusMetres;
};

double initialBearing(double lat1, double lon1, double lat2, double lon2) {
  const double dlon = lon2 - lon1;
  return std::atan2(std::sin(dlon) * std::cos(lat2),
      std::cos(lat1) * std::sin(lat2) -
      std::sin(lat1) * std::cos(lat2) * std::cos(dlon));
}

std::array<double, 2> destination(double lat, double lon, double bearing,
                                  double angularDistance) {
  const double lat2 = std::asin(std::sin(lat) * std::cos(angularDistance) +
      std::cos(lat) * std::sin(angularDistance) * std::cos(bearing));
  const double lon2 = lon + std::atan2(
      std::sin(bearing) * std::sin(angularDistance) * std::cos(lat),
      std::cos(angularDistance) - std::sin(lat) * std::sin(lat2));
  return {lat2, lon2};
}

std::vector<SectionPoint> makeSections(const std::vector<double> &radii,
                                       double earthRadius) {
  const double deg = static_cast<double>(EIGEN_PI) / 180.0;
  const double lat1 = sourceLatitude * deg, lon1 = sourceLongitude * deg;
  const double lat2 = receiverLatitude * deg, lon2 = receiverLongitude * deg;
  const double azimuth = initialBearing(lat1, lon1, lat2, lon2);
  std::vector<SectionPoint> sections;
  constexpr int nDistance = 41;
  constexpr double maxAngle = 60.0 * static_cast<double>(EIGEN_PI) / 180.0;
  const std::array<std::pair<const char *, double>, 2> planes{{
      {"great_circle", azimuth}, {"cross_plane", azimuth + static_cast<double>(EIGEN_PI) / 2.0}}};
  for (const auto &[name, bearing] : planes) {
    for (int j = 0; j < nDistance; ++j) {
      const double signedDistance = -maxAngle + 2.0 * maxAngle * j / (nDistance - 1);
      const double direction = signedDistance < 0.0 ? bearing + static_cast<double>(EIGEN_PI) : bearing;
      const auto [lat, lon] = destination(lat1, lon1, direction, std::abs(signedDistance));
      const double theta = static_cast<double>(EIGEN_PI) / 2.0 - lat;
      for (double radius : radii)
        sections.push_back({name, signedDistance * earthRadius, radius, theta, lon});
    }
  }
  return sections;
}

int main(int argc, char **argv) try {
  if (argc != 4)
    throw std::runtime_error("usage: dspec_3d_kernel_consumer MODEL OUTPUT LMAX");
  const fs::path modelPath = fs::absolute(argv[1]);
  const fs::path out = fs::absolute(argv[2]);
  const int lmax = std::stoi(argv[3]);
  if (!fs::exists(out) || !fs::is_directory(out) ||
      fs::directory_iterator(out) != fs::directory_iterator{})
    throw std::runtime_error("output must be a newly created empty directory");
  if (lmax != 3 && lmax != 12)
    throw std::runtime_error("lmax preset must be 3 or 12");

  {
    std::ofstream file(out / "params.txt");
    file << paramsText(modelPath, out, lmax);
  }
  constexpr double maxstep = 0.05;
  constexpr int nq = 8;
  InputParametersNew input((out / "params.txt").string(), nq, 1, maxstep, 1.0, 0.1, 0.0, 1);
  auto &parameters = input.inputParameters();
  auto &frequency = input.freqFull();
  auto &source = input.cmt();
  Full1D::SEM sem(input);
  if (frequency.nt() != nt || frequency.i1() > 1 || frequency.i2() < 4)
    throw std::runtime_error("resolved grid must contain raw bins 1, 2, and 3");

  const auto rule = angularRule(lmax);
  const double dphi = 2.0 * std::numbers::pi / rule.nphi;
  std::vector<std::array<double, 2>> angles;
  for (int it = 0; it < rule.polar.N(); ++it) {
    const double theta = std::acos(rule.polar.X(it));
    for (int ip = 0; ip < rule.nphi; ++ip)
      angles.push_back({theta, dphi * ip});
  }

  std::vector<SPARSESPEC::Field3DPoint> volumePoints;
  std::vector<double> radii;
  std::vector<std::array<int, 2>> meshNodes;
  for (int e = 0; e < static_cast<int>(sem.mesh().NE()); ++e) {
    const double lower = sem.mesh().ELR(e) * sem.lengthScaleMetres();
    const double upper = sem.mesh().EUR(e) * sem.lengthScaleMetres();
    if (lower < r0 - 1e-7 || upper > r1 + 1e-7)
      continue;
    for (int q = 0; q < sem.mesh().NN(); ++q) {
      const double radius = sem.mesh().NodeRadius(e, q) * sem.lengthScaleMetres();
      radii.push_back(radius);
      meshNodes.push_back({e, q});
      for (const auto &angle : angles)
        volumePoints.push_back({radius, angle[0], angle[1], e});
    }
  }
  if (volumePoints.empty())
    throw std::runtime_error("the model has no radial samples inside the regular support");
  std::vector<double> sectionRadii;
  for (int j = 0; j < 21; ++j)
    sectionRadii.push_back(r0 + (r1 - r0) * j / 20.0);
  const auto sectionMeta = makeSections(sectionRadii, sem.mesh().PR() * sem.lengthScaleMetres());
  const double sourceRadius = sem.mesh().SourceCutRadius() * sem.lengthScaleMetres();
  std::vector<RadialNode> sourceLocalNodes;
  std::vector<SPARSESPEC::Field3DPoint> sourceLocalPoints;
  for (int e = 0; e < static_cast<int>(sem.mesh().NE()); ++e) {
    const double lower = sem.mesh().ELR(e) * sem.lengthScaleMetres();
    const double upper = sem.mesh().EUR(e) * sem.lengthScaleMetres();
    if (upper < sourceRadius - 10000.0 || lower > sourceRadius + 10000.0)
      continue;
    for (int q = 0; q < sem.mesh().NN(); ++q) {
      const double radius = sem.mesh().NodeRadius(e, q) * sem.lengthScaleMetres();
      if (radius < sourceRadius - 10000.0 || radius > sourceRadius + 10000.0)
        continue;
      sourceLocalNodes.push_back({e, q, radius});
      for (const auto &angle : angles)
        sourceLocalPoints.push_back({radius, angle[0], angle[1], e});
    }
  }
  if (sourceLocalNodes.empty())
    throw std::runtime_error("the mesh has no GLL nodes in the source-local support");
  std::vector<SPARSESPEC::Field3DPoint> evaluationPoints = volumePoints;
  for (const auto &point : sectionMeta)
    evaluationPoints.push_back({point.radiusMetres, point.theta, point.phi, -1});
  evaluationPoints.insert(evaluationPoints.end(), sourceLocalPoints.begin(),
                          sourceLocalPoints.end());

  const auto q = objectiveWeights();
  const auto result = SPARSESPEC::Windowed3DKernelEngine::evaluate(
      sem, source, parameters, frequency, hannWindow(), q,
      1, 4, evaluationPoints, angles);
  for (const auto &kernel : result.pointKernels)
    for (double value : {kernel.A, kernel.C, kernel.F, kernel.L, kernel.N})
      if (!std::isfinite(value))
        throw std::runtime_error("3D kernel field contains a non-finite value");
  for (const auto &contact : result.sourceContact)
    for (double value : {contact.A, contact.C, contact.F, contact.L, contact.N})
      if (!std::isfinite(value))
        throw std::runtime_error("source-shell contact contains a non-finite value");
  std::ofstream samples(out / "samples.csv");
  samples << "objective,modulus,radius_m,theta_rad,phi_rad,element,fluid,"
             "kernel_objective_per_Pa_m3,deltaX_Pa,dr_m,dOmega\n"
          << std::setprecision(17);
  std::size_t sampleIndex = 0;
  for (const auto &[element, node] : meshNodes) {
    const double radius = sem.mesh().NodeRadius(element, node) * sem.lengthScaleMetres();
    const double dr = 0.5 * sem.mesh().EW(element) *
        const_cast<EarthMesh::RadialMesh &>(sem.mesh()).GLL().W(node) *
        sem.lengthScaleMetres();
    const double radialWeight = std::pow(std::sin(
        std::numbers::pi * (radius - r0) / (r1 - r0)), 2);
    for (std::size_t ia = 0; ia < angles.size(); ++ia, ++sampleIndex) {
      const auto &angle = angles[ia];
      const double perturbation = angularPerturbation(angle[0], angle[1]);
      const auto &kernel = result.pointKernels[sampleIndex];
      const std::array<double, 5> values{kernel.A, kernel.C, kernel.F, kernel.L, kernel.N};
      for (int m = 0; m < 5; ++m)
        samples << "J_Z," << moduli[m] << ',' << radius << ',' << angle[0] << ','
                << angle[1] << ',' << element << ',' << sem.elementIsFluid(element)
                << ',' << values[m] << ',' << baselineModuliPa[m] * radialWeight *
                    perturbation << ',' << dr << ',' << rule.polar.W(ia / rule.nphi) * dphi << '\n';
    }
  }

  std::ofstream sections(out / "sections.csv");
  sections << "plane,distance_m,radius_m,theta_rad,phi_rad,modulus,"
              "kernel_objective_per_Pa_m3\n" << std::setprecision(17);
  for (std::size_t j = 0; j < sectionMeta.size(); ++j) {
    const auto &point = sectionMeta[j];
    const auto &kernel = result.pointKernels[volumePoints.size() + j];
    const std::array<double, 5> values{kernel.A, kernel.C, kernel.F, kernel.L, kernel.N};
    for (int m = 0; m < 5; ++m)
      sections << point.plane << ',' << point.distanceMetres << ',' << point.radiusMetres
               << ',' << point.theta << ',' << point.phi << ',' << moduli[m] << ','
               << values[m] << '\n';
  }

  std::ofstream contacts(out / "source_contacts.csv");
  contacts << "objective,modulus,theta_rad,phi_rad,dOmega,"
              "contact_objective_per_Pa,deltaX_Pa\n" << std::setprecision(17);
  for (std::size_t j = 0; j < angles.size(); ++j) {
    const auto &angle = angles[j];
    const auto &contact = result.sourceContact[j];
    const std::array<double, 5> values{contact.A, contact.C, contact.F, contact.L, contact.N};
    const double dOmega = rule.polar.W(j / rule.nphi) * dphi;
    for (int m = 0; m < 5; ++m)
      contacts << "J_Z," << moduli[m] << ',' << angle[0] << ',' << angle[1] << ','
               << dOmega << ',' << values[m] << ','
               << baselineModuliPa[m] * angularPerturbation(angle[0], angle[1]) << '\n';
  }

  std::array<double, 3> sourceVolume{};
  std::array<double, 3> sourceContact{};
  const std::array<int, 3> sourceModulusIndices{1, 2, 3};
  for (std::size_t ir = 0; ir < sourceLocalNodes.size(); ++ir) {
    const auto &radial = sourceLocalNodes[ir];
    const double dr = 0.5 * sem.mesh().EW(radial.element) *
        const_cast<EarthMesh::RadialMesh &>(sem.mesh()).GLL().W(radial.node) *
        sem.lengthScaleMetres();
    const double offset = (radial.radiusMetres - sourceRadius) / 10000.0;
    const double envelope = std::pow(std::cos(0.5 * std::numbers::pi * offset), 2);
    for (std::size_t ia = 0; ia < angles.size(); ++ia) {
      const double chi = angularPerturbation(angles[ia][0], angles[ia][1]);
      const double dOmega = rule.polar.W(ia / rule.nphi) * dphi;
      const auto &kernel = result.pointKernels[volumePoints.size() + sectionMeta.size() +
                                                ir * angles.size() + ia];
      const std::array<double, 5> values{kernel.A, kernel.C, kernel.F, kernel.L, kernel.N};
      for (std::size_t k = 0; k < sourceModulusIndices.size(); ++k) {
        const int m = sourceModulusIndices[k];
        sourceVolume[k] += values[m] * baselineModuliPa[m] * envelope * chi *
                           radial.radiusMetres * radial.radiusMetres * dr * dOmega;
      }
    }
  }
  for (std::size_t ia = 0; ia < angles.size(); ++ia) {
    const double chi = angularPerturbation(angles[ia][0], angles[ia][1]);
    const double dOmega = rule.polar.W(ia / rule.nphi) * dphi;
    const auto &contact = result.sourceContact[ia];
    const std::array<double, 5> values{contact.A, contact.C, contact.F, contact.L, contact.N};
    for (std::size_t k = 0; k < sourceModulusIndices.size(); ++k) {
      const int m = sourceModulusIndices[k];
      sourceContact[k] += values[m] * baselineModuliPa[m] * chi * dOmega;
    }
  }
  std::ofstream sourceSplit(out / "source_split.csv");
  sourceSplit << "objective,modulus,source_radius_m,half_width_m,volume_integral,"
                 "contact_integral,total_integral,units\n" << std::setprecision(17);
  for (std::size_t k = 0; k < sourceModulusIndices.size(); ++k) {
    const int m = sourceModulusIndices[k];
    sourceSplit << "J_Z," << moduli[m] << ',' << sourceRadius << ",10000,"
                << sourceVolume[k] << ',' << sourceContact[k] << ','
                << sourceVolume[k] + sourceContact[k]
                << ",objective per unit fractional perturbation\n";
  }

  std::ofstream meta(out / "run.json");
  const auto &model = input.earthModel();
  meta << std::setprecision(17)
       << "{\"solver\":\"DSpecM1D\",\"objective\":\"J_Z=Re(sum(conj(q)*P u))\","
       << "\"objective_units\":\"m s\","
       << "\"lmax\":" << lmax << ",\"frequency_bins_zero_based\":[1,2,3],"
       << "\"degrees\":[0," << lmax << "],\"solver_bins\":[1,2,3],"
       << "\"nt\":" << frequency.nt() << ",\"dt_s\":"
       << frequency.dt() * frequency.timeNorm() << ",\"df_hz\":"
       << frequency.df() / frequency.timeNorm() << ",\"epsilon_s_inv\":"
       << frequency.ep() / frequency.timeNorm() << ",\"radius_m\":"
       << sem.mesh().PR() * sem.lengthScaleMetres() << ",\"sample_count\":"
       << volumePoints.size() << ",\"angular_nodes\":" << angles.size()
       << ",\"reference_moduli_Pa\":[352960000000,352960000000,129602500000,111678750000,111678750000]"
       << ",\"unit_norms\":{\"length_m\":" << model.LengthNorm()
       << ",\"mass_kg\":" << model.MassNorm() << ",\"time_s\":"
       << model.TimeNorm() << ",\"density_kg_m3\":" << model.DensityNorm()
       << ",\"gravitational_constant_unit_m3_kg_s2\":" << model.GravitationalConstant()
       << ",\"dimensionless_G\":" << sem.gravityCoefficient()
       << ",\"attenuation_reference_period_model_units\":" << model.TREF()
       << ",\"attenuation_reference_period_s\":" << model.TREF() * model.TimeNorm() << "}"
       << ",\"physical_G_SI\":" << model.GravitationalConstant() * sem.gravityCoefficient()
       << ",\"source\":{\"latitude_deg\":" << parameters.source_lat_deg()
       << ",\"longitude_deg\":" << parameters.source_lon_deg() << ",\"depth_km\":"
       << parameters.source_depth_km() << ",\"radius_m\":" << sourceRadius << "}"
       << ",\"moment_tensor_Nm\":[" << parameters.m_rr() << ',' << parameters.m_rt()
       << ',' << parameters.m_rp() << ',' << parameters.m_tt() << ',' << parameters.m_tp()
       << ',' << parameters.m_pp() << ']'
       << ",\"receiver\":{\"latitude_deg\":45,\"longitude_deg\":90,\"depth_km\":"
       << parameters.receiver_depth() << "}"
       << ",\"gravity\":\"full\",\"attenuation\":"
       << (parameters.attenuation() != 0 ? "true" : "false")
       << ",\"time_window\":{\"type\":\"Hann\",\"sample_start_inclusive\":0,"
          "\"sample_end_inclusive\":90,\"sample_count\":128}"
       << ",\"source_local_profile\":\"deltaX=X0*cos^2(pi*(r-rs)/20000 m)*chi for |r-rs|<=10000 m; C,F,L\""
       << ",\"angular_rule\":\"Gauss-Legendre max(8,2*lmax+2) x periodic max(24,4*lmax+4)\","
       << "\"coordinate_convention\":\"geocentric colatitude theta and east longitude phi, radians\","
       << "\"kernel_units\":\"objective/(Pa m^3); integrate with r^2 dr dOmega\","
       << "\"contact_units\":\"objective/Pa per steradian; integrate with dOmega\","
       << "\"source_shell_contact_saved_separately\":true,"
       << "\"source_split_units\":\"objective per unit fractional perturbation\","
       << "\"q_nonzero\":{\"Z,1\":[0.7,0.2],\"Z,2\":[-0.25,0.45]},"
       << "\"objective_scale_S\":" << q.norm() * result.processedSISpectrum.norm() << ','
       << "\"convergence_claim\":false}\n";
  return 0;
} catch (const std::exception &error) {
  std::cerr << "3D kernel consumer: " << error.what() << '\n';
  return 1;
}
