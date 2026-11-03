#pragma once

// What the radial and the volume kernels share: the checks of their inputs,
// the solver of the forward and adjoint fields, the spectrum of the adjoint
// field's force, and the parameters in which the kernels are expressed.

#include <YSpec/Calculation.hpp>
#include <YSpec/Kernels.hpp>
#include <YSpec/Model.hpp>
#include <YSpec/RadialSolver.hpp>
#include <YSpec/Source.hpp>
#include <array>
#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "SourceWeights.hpp"

namespace YSpec {

// The radii of the source and of the receiver, non-dimensional, as the
// seismograms have them.
struct KernelRadii {
  double source;
  double receiver;
};

// Checks what kernels are asked of: the degrees, the corrections, which have
// no adjoint here, an adjoint source at every sample, and the radii.
KernelRadii CheckKernelInputs(const Model& model,
                              const Calculation& calculation,
                              const Source& source,
                              const Measurement& measurement);

// The same for several measurements, which must share their receiver: its
// place and depth, though not its instrument.
KernelRadii CheckKernelInputs(const Model& model,
                              const Calculation& calculation,
                              const Source& source,
                              std::span<const Measurement> measurements);

// An orthonormal basis of the measurements' directions, of one, two or three
// directions, and the coordinates of each measurement's direction in it: the
// adjoint field is linear in its force, so those of the basis serve all.
struct DirectionBasis {
  std::vector<Direction> directions;
  std::vector<std::vector<double>> coordinates;  // [measurement][direction]
};

DirectionBasis MakeBasis(std::span<const Measurement> measurements);

// The solver of the forward field, from source 0, and of the adjoint field,
// from source AdjointIndex, at the points given, and at the source's own
// radius, last, if asked.
RadialSolver KernelSolver(const Model& model, const Calculation& calculation,
                          const KernelRadii& radii, std::vector<Point> points,
                          std::optional<double> atSource = std::nullopt);

inline std::size_t AdjointIndex(const KernelRadii& radii) {
  return radii.receiver != radii.source ? 1 : 0;
}

// The kernels of the model's own parameters accumulated at points,
// non-dimensional, density's part through gravity left out, and after them
// the two angular integrals from which that part follows, M_rr and the flux;
// and the point parts at the source, at the point given.
using RadialAccumulator =
    std::array<double, static_cast<std::size_t>(RadialParameter::Count) + 2>;

struct RadialAccumulation {
  std::vector<RadialAccumulator> points;
  RadialAccumulator atSource;
  Point source;
};

// For several measurements at one receiver, one accumulation each.
std::vector<RadialAccumulation> AccumulateRadial(
    const Model& model, const Calculation& calculation, const Source& source,
    std::span<const Measurement> measurements, std::span<const Point> points,
    const Progress& progress = {});

// The spectrum of the adjoint field's force at each frequency of the grid,
// A_k, and the weight of each frequency in the sum that replaces the
// integral over frequency.
struct AdjointSpectrum {
  std::vector<Complex> force;
  std::vector<double> weight;
};

AdjointSpectrum MakeAdjointSpectrum(const Model& model,
                                    const Calculation& calculation,
                                    const Measurement& measurement);

// Multiplies every weight of a source.
void Scale(SourceWeights& weights, Complex factor);

// A force along a direction at the receiver, in (r, theta, phi) of the frame
// with the source at the pole, where the angle given turns (theta, phi) to
// (south, east).
Force InSourceFrame(const Direction& direction, double rotation);

// The model at a point, in SI units, as the reparameterisation needs it.
struct ModelValues {
  double rho, A, C, F, L, N, kappa, mu, qA, qF, qMu, qKappa;
};

// At a radius, non-dimensional, in a layer.
ModelValues ValuesAt(const Model& model, Int layer, double x);

// The names of the parameters of a set, density first.
std::vector<std::string> Names(ParameterSet set);

// The kernels of a set of parameters at a point, and the values of the
// parameters, from the kernels of the model's own.
std::pair<std::vector<double>, std::vector<double>> Transform(
    const std::array<double, static_cast<std::size_t>(RadialParameter::Count)>&
        k,
    const ModelValues& v, ParameterSet set);

}  // namespace YSpec
