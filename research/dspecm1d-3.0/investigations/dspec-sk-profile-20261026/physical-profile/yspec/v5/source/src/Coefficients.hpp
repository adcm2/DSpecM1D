#pragma once

// The displacement coefficients of the programs: for each source, receiver
// depth, degree and frequency, the coefficients of the displacement in the
// frame with the source at the north pole, order by order. Seismograms sum
// them at their receivers; snapshots at the points of their grids.

#include <YSpec/Calculation.hpp>
#include <YSpec/Model.hpp>
#include <YSpec/RadialSolver.hpp>
#include <YSpec/Source.hpp>
#include <optional>
#include <span>
#include <vector>

#include "Equations.hpp"
#include "OrderArray.hpp"
#include "SourceWeights.hpp"
#include "TimeSeries.hpp"

namespace YSpec {

// The frequency grid of a calculation, non-dimensional.
FrequencyGrid MakeGrid(const Calculation& calculation, const Scales& scales);

// The taper of a calculation, if any, non-dimensional.
std::optional<Taper> MakeTaper(const Calculation& calculation,
                               const Scales& scales);

// i z, without a complex multiplication.
inline Complex TimesI(Complex z) { return {-z.imag(), z.real()}; }

// The SI unit of a quantity, in the model's scales.
double Unit(const Scales& scales, Quantity quantity);

// The number of time derivatives of displacement a quantity is.
int Derivatives(Quantity quantity);

// The members of a set that share a depth, and the radius of that depth.
struct DepthGroup {
  double radius;
  std::vector<std::size_t> members;
};

// What the corrections need at a receiver depth.
struct ReceiverGravity {
  double r;   // The radius.
  double g;   // Gravity.
  double dg;  // Its radial derivative.
};

// The radial solves for a set of sources and of receiver depths (m): the
// depths grouped, one solver for every one, and what the corrections need at
// each receiver depth.
struct DepthSolves {
  DepthSolves(const Model& model, const Calculation& calculation,
              std::span<const Source> sources, std::span<const double> depths);

  std::vector<DepthGroup> sourceDepths;
  std::vector<DepthGroup> receiverDepths;
  RadialSolver solver;
  std::vector<ReceiverGravity> gravity;
};

// The responses at one degree and frequency for every source depth and
// receiver depth, converted to the quantity recorded and corrected, if asked,
// for the potential and tilt.
struct Responses {
  ResponseTable<ToroidalResponse> toroidal{0, 0};
  ResponseTable<SpheroidalResponse> spheroidal{0, 0};
  ResponseTable<RadialResponse> radial{0, 0};
};

Responses Record(const DepthSolves& solves, int l, Complex omega,
                 const Calculation& calculation);

// The displacement coefficients for one source, order by order, and their
// radial derivatives (that of v is zero in a fluid), with those of the
// perturbation of the potential (zero without self-gravitation, and at l =
// 0, where it is not formed).
struct Coefficients {
  OrderArray<Complex> u;   // radial
  OrderArray<Complex> v;   // consoidal
  OrderArray<Complex> w;   // toroidal
  OrderArray<Complex> du;  // d/dr
  OrderArray<Complex> dv;
  OrderArray<Complex> dw;
  OrderArray<Complex> phi;
  OrderArray<Complex> dphi;
};

// The responses to a source's jumps, combined with its weights.
Coefficients Combine(const Responses& responses, std::size_t sourceDepth,
                     std::size_t receiverDepth, int l,
                     const SourceWeights& weights);

// With the moment tensor and the force non-dimensional.
Coefficients Weight(const Responses& responses, std::size_t sourceDepth,
                    std::size_t receiverDepth, const Context& c,
                    const ComplexMaterial& atSource, const MomentTensor& moment,
                    const Force& force);

// Calls each(source, receiverDepth, coefficients) for every source and every
// receiver depth, at one degree and frequency.
template <class Each>
void ForEachSource(const DepthSolves& solves, const Model& model,
                   const Calculation& calculation,
                   std::span<const Source> sources, int l, Complex omega,
                   Each&& each) {
  const auto responses = Record(solves, l, omega, calculation);
  const auto c = Context(model, l, omega, calculation.attenuation);
  const auto& scales = model.GetScales();
  for (auto sourceDepth :
       std::views::iota(std::size_t{0}, solves.sourceDepths.size())) {
    const auto point = solves.solver.Source(sourceDepth);
    const auto atSource = c(point.layer, point.r);
    for (auto source : solves.sourceDepths[sourceDepth].members) {
      const auto moment = sources[source].moment / scales.Moment();
      const auto force = sources[source].force / scales.Force();
      for (auto receiverDepth :
           std::views::iota(std::size_t{0}, solves.receiverDepths.size())) {
        each(source, receiverDepth,
             Weight(responses, sourceDepth, receiverDepth, c, atSource, moment,
                    force));
      }
    }
  }
}

}  // namespace YSpec
