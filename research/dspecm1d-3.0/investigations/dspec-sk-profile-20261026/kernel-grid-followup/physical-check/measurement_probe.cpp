#include "/home/adcm2/Documents/Research/Projects/DSpecM1D/DSpecM1D/src/SensitivityKernels.h"

#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
using Complex = std::complex<double>;
using MatrixC = Eigen::MatrixXcd;
using MatrixR = Eigen::MatrixXd;
using VectorR = Eigen::VectorXd;

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kLengthScaleMetres = 6371000.0;

void require(bool ok, const char *message) {
  if (!ok) throw std::runtime_error(message);
}

double relativeError(double a, double b) {
  return std::abs(a - b) /
         std::max({std::abs(a), std::abs(b), std::numeric_limits<double>::min()});
}

MatrixC makeRawSpectrum(const VectorR &time, double dt) {
  MatrixR row(1, time.size());
  row.row(0) = time.transpose();
  return processfunctions::rawtime2freq(row, static_cast<int>(time.size()), dt);
}

MatrixC applyBand(const MatrixC &spectrum, const SpectraSolver::FreqFull &frequency) {
  MatrixC filtered = spectrum;
  for (int k = 0; k < filtered.cols(); ++k) {
    const double fMhz = 1000.0 * frequency.f(k);
    double weight = 0.0;
    if (fMhz < 8.0) {
      weight = 1.0;
    } else if (fMhz < 10.0) {
      weight = 0.5 * (1.0 + std::cos(kPi * (fMhz - 8.0) / 2.0));
    }
    filtered(0, k) *= weight;
  }
  return filtered;
}

VectorR makeWindow(const SpectraSolver::FreqFull &frequency) {
  VectorR window(frequency.nt());
  const double dtSI = frequency.dt() * frequency.timeNorm();
  for (int n = 0; n < frequency.nt(); ++n) {
    const double t = n * dtSI;
    if (t < 1000.0 || t > 1200.0) {
      window[n] = 0.0;
    } else {
      const double phase = 2.0 * kPi * (t - 1000.0) / 200.0;
      window[n] = std::sqrt(0.5 * (1.0 - std::cos(phase)));
    }
  }
  return window;
}

VectorR wavePacket(const SpectraSolver::FreqFull &frequency, double centre,
                   double width, double fHz, double phase) {
  VectorR values(frequency.nt());
  const double dtSI = frequency.dt() * frequency.timeNorm();
  for (int n = 0; n < frequency.nt(); ++n) {
    const double t = n * dtSI;
    const double x = (t - centre) / width;
    const double envelope = std::exp(-0.5 * x * x);
    values[n] = 1.0e-8 * envelope * std::cos(2.0 * kPi * fHz * t + phase);
  }
  return values;
}

MatrixC dampedSpectrum(const VectorR &undampedTime,
                       const SpectraSolver::FreqFull &frequency) {
  VectorR damped = undampedTime;
  for (int n = 0; n < frequency.nt(); ++n)
    damped[n] *= std::exp(-frequency.ep() * frequency.dt() * n);
  return makeRawSpectrum(damped, frequency.dt());
}

double directEnergy(const MatrixC &rawSpectrum,
                    const SpectraSolver::FreqFull &frequency,
                    const VectorR &window) {
  MatrixR undamped = processfunctions::freq2time(rawSpectrum, frequency);
  undamped.row(0) *= kLengthScaleMetres;
  undamped.row(0).array() *= window.transpose().array();
  const double dtSI = frequency.dt() * frequency.timeNorm();
  return 0.5 * dtSI * undamped.squaredNorm();
}

double spectralEnergy(const MatrixC &processed,
                      const SpectraSolver::FreqFull &frequency) {
  const int nt = frequency.nt();
  VectorR eta = VectorR::Constant(processed.cols(), 2.0);
  eta[0] = 1.0;
  if (nt % 2 == 0) eta[processed.cols() - 1] = 1.0;
  return 0.5 * frequency.df() *
         (processed.cwiseAbs2().array().rowwise() * eta.transpose().array()).sum();
}

double realInner(const MatrixC &a, const MatrixC &b) {
  return (a.conjugate().cwiseProduct(b)).sum().real();
}

} // namespace

int main() {
  try {
    // f1/f2 are in mHz; dt is seconds, tout/t2 are hours, qex=1.
    const SpectraSolver::FreqFull frequency(0.0, 10.0, 5.0, 1.0,
                                             0.1, 0.0, 0.0, 1.0, 1);
    require(frequency.nt() == 1024, "fixture FFT length changed");
    require(std::abs(frequency.dt() - 5.0) < 1e-14,
            "fixture dt changed");
    require(std::abs(frequency.df() - 1.0 / 5120.0) < 1e-15,
            "fixture df changed");
    require(std::abs(frequency.ep() - 5.0 / 3600.0) < 1e-15,
            "fixture damping changed");

    const VectorR window = makeWindow(frequency);
    SPARSESPEC::WindowedDisplacementMap measurement(
        frequency, kLengthScaleMetres, window);

    const double dtSI = frequency.dt() * frequency.timeNorm();
    const VectorR primary = wavePacket(frequency, 1100.0, 75.0, 0.0034, 0.3);
    const VectorR directionTime = wavePacket(frequency, 1080.0, 82.0, 0.0041, -0.6);
    MatrixC raw = dampedSpectrum(primary, frequency);
    MatrixC direction = dampedSpectrum(directionTime, frequency);
    direction *= raw.cwiseAbs().maxCoeff() / direction.cwiseAbs().maxCoeff();
    require(raw.allFinite() && direction.allFinite() && window.allFinite(),
            "synthetic fixture contains non-finite values");
    require(std::abs(raw(0, 0).imag()) < 1e-13 &&
                std::abs(raw(0, raw.cols() - 1).imag()) < 1e-13,
            "real-time fixture endpoints must be real");
    require(std::abs(direction(0, 0).imag()) < 1e-13 &&
                std::abs(direction(0, direction.cols() - 1).imag()) < 1e-13,
            "real-time direction endpoints must be real");

    const MatrixC filtered = applyBand(raw, frequency);
    const MatrixC processed = measurement.forward(filtered);
    const double eTime = directEnergy(filtered, frequency, window);
    const double eParseval = spectralEnergy(processed, frequency);
    require(std::isfinite(eTime) && std::isfinite(eParseval) && eTime > 1e-12,
            "fixture energy must be finite and meaningfully nonzero");
    const double parsevalError = relativeError(eTime, eParseval);
    require(parsevalError < 2e-11,
            "direct time energy and one-sided Parseval energy disagree");

    MatrixC testWeights(1, measurement.frequencyCount());
    for (int k = 0; k < testWeights.cols(); ++k)
      testWeights(0, k) = Complex(std::sin(0.37 * k), std::cos(0.19 * k));
    testWeights(0, 0) = Complex(testWeights(0, 0).real(), 0.0);
    testWeights(0, testWeights.cols() - 1) =
        Complex(testWeights(0, testWeights.cols() - 1).real(), 0.0);
    const double lhs = realInner(measurement.forward(filtered), testWeights);
    const double rhs = realInner(filtered, measurement.adjoint(testWeights));
    const double adjointError = relativeError(lhs, rhs);
    require(adjointError < 2e-10,
            "WindowedDisplacementMap real-inner-product adjoint check failed");

    VectorR eta = VectorR::Constant(processed.cols(), 2.0);
    eta[0] = 1.0;
    eta[processed.cols() - 1] = 1.0;
    MatrixC q = processed;
    for (int k = 0; k < q.cols(); ++k)
      q(0, k) *= frequency.df() * eta[k];
    const MatrixC gradient = applyBand(measurement.adjoint(q), frequency);
    const double directionalDerivative = realInner(gradient, direction);
    require(std::isfinite(directionalDerivative) &&
                std::abs(directionalDerivative) > 1e-12,
            "fixture directional derivative must be meaningfully nonzero");
    const std::array<double, 3> steps{1e-2, 1e-3, 1e-4};
    std::array<double, 3> fdValues{};
    std::array<double, 3> fdErrors{};
    for (std::size_t i = 0; i < steps.size(); ++i) {
      const double h = steps[i];
      const double ePlus = directEnergy(
          applyBand(raw + h * direction, frequency), frequency, window);
      const double eMinus = directEnergy(
          applyBand(raw - h * direction, frequency), frequency, window);
      fdValues[i] = (ePlus - eMinus) / (2.0 * h);
      fdErrors[i] = relativeError(fdValues[i], directionalDerivative);
      require(std::isfinite(fdValues[i]) && fdErrors[i] < 2e-7,
              "central finite difference disagrees with frozen-baseline gradient");
    }

    std::cout << std::setprecision(17)
              << "PASS synthetic measurement probe\n"
              << "nt=" << frequency.nt() << " dt_s=" << dtSI
              << " df_hz=" << frequency.df() << " ep_per_s=" << frequency.ep()
              << " tfft_s=" << frequency.nt() * dtSI << '\n'
              << "energy_time=" << eTime << " energy_parseval=" << eParseval
              << " parseval_absolute_error=" << std::abs(eTime - eParseval)
              << " parseval_relative_error=" << parsevalError << '\n'
              << "adjoint_lhs=" << lhs << " adjoint_rhs=" << rhs
              << " adjoint_absolute_error=" << std::abs(lhs - rhs)
              << " adjoint_relative_error=" << adjointError << '\n'
              << "frozen_gradient_derivative=" << directionalDerivative << '\n'
              << "fd_h=" << steps[0] << " fd_derivative=" << fdValues[0]
              << " fd_absolute_error=" << std::abs(fdValues[0] - directionalDerivative)
              << " fd_relative_error=" << fdErrors[0] << '\n'
              << "fd_h=" << steps[1] << " fd_derivative=" << fdValues[1]
              << " fd_absolute_error=" << std::abs(fdValues[1] - directionalDerivative)
              << " fd_relative_error=" << fdErrors[1] << '\n'
              << "fd_h=" << steps[2] << " fd_derivative=" << fdValues[2]
              << " fd_absolute_error=" << std::abs(fdValues[2] - directionalDerivative)
              << " fd_relative_error=" << fdErrors[2] << '\n'
              << "qdotp_equals_2E="
              << (frequency.df() *
                  (processed.cwiseAbs2().array().rowwise() * eta.transpose().array()).sum())
              << "\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "FAIL: " << error.what() << '\n';
    return 1;
  }
}
