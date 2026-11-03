#include "Coefficients.hpp"

#include <format>
#include <map>
#include <ranges>
#include <stdexcept>

namespace YSpec {

using namespace std::complex_literals;

FrequencyGrid MakeGrid(const Calculation& calculation, const Scales& scales) {
  const auto time = scales.Time();
  const auto& sampling = calculation.sampling;
  return FrequencyGrid(calculation.band.min * time, calculation.band.max * time,
                       sampling.timeStep / time, sampling.duration / time);
}

std::optional<Taper> MakeTaper(const Calculation& calculation,
                               const Scales& scales) {
  auto taper = calculation.taper;
  if (taper) {
    for (auto* f : {&taper->f1, &taper->f2, &taper->f3, &taper->f4}) {
      *f *= scales.Time();
    }
  }
  return taper;
}

// The SI unit of a quantity, in the model's scales.
double Unit(const Scales& scales, Quantity quantity) {
  switch (quantity) {
    case Quantity::Displacement:
      return scales.Length();
    case Quantity::Velocity:
      return scales.Velocity();
    case Quantity::Acceleration:
      return scales.Acceleration();
  }
  return scales.Length();
}

// The number of time derivatives of displacement a quantity is.
int Derivatives(Quantity quantity) {
  switch (quantity) {
    case Quantity::Displacement:
      return 0;
    case Quantity::Velocity:
      return 1;
    case Quantity::Acceleration:
      return 2;
  }
  return 0;
}

namespace {

// Groups depths (m), as radii.
std::vector<DepthGroup> GroupByDepth(std::span<const double> depths,
                                     const Model& model, const char* what) {
  const auto surface = model.SurfaceRadius();
  const auto length = model.GetScales().Length();
  auto groups = std::map<double, std::vector<std::size_t>>{};
  for (auto i : std::views::iota(std::size_t{0}, depths.size())) {
    groups[depths[i]].push_back(i);
  }
  auto result = std::vector<DepthGroup>{};
  for (auto& [depth, members] : groups) {
    const auto r = surface - depth / length;
    if (depth < 0 || r < model.InnerRadius()) {
      throw std::invalid_argument(std::format(
          "a {} depth of {} km is outside the model, whose radius is {} km",
          what, depth / 1e3, surface * length / 1e3));
    }
    result.push_back({r, std::move(members)});
  }
  return result;
}

std::vector<double> Radii(const std::vector<DepthGroup>& groups) {
  auto radii = std::vector<double>{};
  for (const auto& group : groups) radii.push_back(group.radius);
  return radii;
}

std::vector<double> SourceDepths(std::span<const Source> sources) {
  auto depths = std::vector<double>{};
  for (const auto& source : sources) depths.push_back(source.depth);
  return depths;
}

}  // namespace

DepthSolves::DepthSolves(const Model& model, const Calculation& calculation,
                         std::span<const Source> sources,
                         std::span<const double> depths)
    : sourceDepths{GroupByDepth(SourceDepths(sources), model, "source")},
      receiverDepths{GroupByDepth(depths, model, "receiver")},
      solver{model,
             Radii(sourceDepths),
             Radii(receiverDepths),
             calculation.gravity,
             calculation.attenuation,
             calculation.solver} {
  for (auto i : std::views::iota(std::size_t{0}, solver.Receivers())) {
    const auto receiver = solver.Receiver(i);
    const auto at = model(receiver.layer, receiver.r);
    gravity.push_back({.r = receiver.r,
                       .g = at.g,
                       .dg = 4 * PiG * at.rho - 2 * at.g / receiver.r});
  }
}

Responses Record(const DepthSolves& solves, int l, Complex omega,
                 const Calculation& calculation) {
  const auto& solver = solves.solver;
  const auto& gravity = solves.gravity;
  const auto derivative = [&](Complex y) {
    switch (calculation.quantity) {
      case Quantity::Displacement:
        return y;
      case Quantity::Velocity:
        return 1i * omega * y;
      case Quantity::Acceleration:
        return -omega * omega * y;
    }
    return y;
  };
  const auto corrections =
      calculation.corrections && calculation.quantity == Quantity::Acceleration;
  const auto full = calculation.gravity == Gravity::Full;
  auto result = Responses{};
  if (l == 0) {
    result.radial = solver.Radial(omega);
  } else {
    result.toroidal = solver.Toroidal(l, omega);
    result.spheroidal = solver.Spheroidal(l, omega);
  }
  for (auto source : std::views::iota(std::size_t{0}, solver.Sources())) {
    for (auto receiver : std::views::iota(std::size_t{0}, solver.Receivers())) {
      const auto& at = gravity[receiver];
      if (l == 0) {
        auto& radial = result.radial(source, receiver);
        for (auto& u : radial.U) {
          // For l = 0 the change in gravity reduces to the free-air term.
          u = derivative(u) - (corrections ? 2 * at.g * u / at.r : 0.0);
        }
        for (auto& du : radial.DU) du = derivative(du);
        continue;
      }
      auto& toroidal = result.toroidal(source, receiver);
      for (auto& w : toroidal.W) w = derivative(w);
      for (auto& dw : toroidal.DW) dw = derivative(dw);
      auto& spheroidal = result.spheroidal(source, receiver);
      for (auto j : {0, 1, 2, 3}) {
        const auto U = spheroidal.U[j];
        const auto V = spheroidal.V[j];
        spheroidal.U[j] = derivative(U);
        spheroidal.V[j] = derivative(V);
        spheroidal.DU[j] = derivative(spheroidal.DU[j]);
        spheroidal.DV[j] = derivative(spheroidal.DV[j]);
        if (!corrections) continue;
        if (full) {
          // The change in gravity at the moving instrument, phi' + U g',
          // and the tilt of its axis, (phi + g U) / r.
          spheroidal.U[j] += spheroidal.DPhi[j] + U * at.dg;
          spheroidal.V[j] += (spheroidal.Phi[j] + at.g * U) / at.r;
        } else {
          // Without the potential, the free-air term and the tilt.
          spheroidal.U[j] -= 2 * at.g * U / at.r;
          spheroidal.V[j] += at.g * U / at.r;
        }
      }
    }
  }
  return result;
}

Coefficients Combine(const Responses& responses, std::size_t sourceDepth,
                     std::size_t receiverDepth, int l,
                     const SourceWeights& weights) {
  auto result = Coefficients{};
  if (l == 0) {
    const auto& radial = responses.radial(sourceDepth, receiverDepth);
    const auto& w = weights.radial;
    result.u[0] = w[0] * radial.U[0] + w[1] * radial.U[1];
    result.du[0] = w[0] * radial.DU[0] + w[1] * radial.DU[1];
    return result;
  }
  const auto& toroidal = responses.toroidal(sourceDepth, receiverDepth);
  const auto& spheroidal = responses.spheroidal(sourceDepth, receiverDepth);
  for (auto m : OrderArray<Complex>::Orders(l)) {
    for (auto j : {0, 1}) {
      const auto weight = weights.toroidal[m][j];
      result.w[m] += weight * toroidal.W[j];
      result.dw[m] += weight * toroidal.DW[j];
    }
    for (auto j : {0, 1, 2, 3}) {
      const auto weight = weights.spheroidal[m][j];
      result.u[m] += weight * spheroidal.U[j];
      result.v[m] += weight * spheroidal.V[j];
      result.du[m] += weight * spheroidal.DU[j];
      result.dv[m] += weight * spheroidal.DV[j];
      result.phi[m] += weight * spheroidal.Phi[j];
      result.dphi[m] += weight * spheroidal.DPhi[j];
    }
  }
  return result;
}

Coefficients Weight(const Responses& responses, std::size_t sourceDepth,
                    std::size_t receiverDepth, const Context& c,
                    const ComplexMaterial& atSource, const MomentTensor& moment,
                    const Force& force) {
  return Combine(responses, sourceDepth, receiverDepth, c.Degree(),
                 PoleWeights(c, atSource, moment, force));
}

}  // namespace YSpec
