#include <GaussQuad/All>
#include <YSpec/Kernels.hpp>
#include <YSpec/RadialSolver.hpp>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <map>
#include <mutex>
#include <numbers>
#include <ranges>
#include <stdexcept>
#include <utility>

#include "AngularFunctions.hpp"
#include "Coefficients.hpp"
#include "Geometry.hpp"
#include "KernelParts.hpp"
#include "OrderArray.hpp"
#include "Report.hpp"
#include "SourceWeights.hpp"
#include "Sphere.hpp"

namespace YSpec {

SampledKernels::SampledKernels(std::vector<std::string> names,
                               KernelSampling sampling)
    : names_{std::move(names)}, sampling_{std::move(sampling)} {
  values_.assign(names_.size() * Nodes().size() * Places().size(), 0.0);
}

namespace {

using namespace std::complex_literals;

// The spin components of the strain, as the strain seismograms have them:
// rr and the horizontal trace, of spin weight 0; r+ and r-, the components
// e_rtheta -+ i e_rphi of the shear on horizontal planes, of weights +1 and
// -1; and ++ and --, the components e_thetatheta - e_phiphi -+ 2i
// e_thetaphi of the traceless horizontal part, of weights +2 and -2.
constexpr auto Components = std::size_t{6};
constexpr auto Spin = std::array{0, 0, 1, -1, 2, -2};

// The coefficients of the spin components of the strain at one radius and
// degree, order by order, from those of the displacement. In a fluid, where
// the radial derivative of the tangential displacement is not given, r+ and
// r- are left zero.
std::array<OrderArray<Complex>, Components> StrainOf(const Coefficients& c,
                                                     int l, double r,
                                                     bool fluid) {
  const auto zeta2 = l * (l + 1.0);
  const auto zeta = std::sqrt(zeta2);
  const auto twist = l >= 2 ? zeta * std::sqrt((l - 1.0) * (l + 2)) / r : 0.0;
  auto strain = std::array<OrderArray<Complex>, Components>{};
  for (auto m : OrderArray<Complex>::Orders(l)) {
    strain[0][m] = c.du[m];
    strain[1][m] = (2.0 * c.u[m] - zeta2 * c.v[m]) / r;
    if (!fluid) {
      const auto shear = c.dv[m] - (c.v[m] - c.u[m]) / r;
      const auto turn = c.dw[m] - c.w[m] / r;
      strain[2][m] = -zeta / 2 * (shear + 1i * turn);
      strain[3][m] = zeta / 2 * (shear - 1i * turn);
    }
    strain[4][m] = twist * (c.v[m] + 1i * c.w[m]);
    strain[5][m] = twist * (c.v[m] - 1i * c.w[m]);
  }
  return strain;
}

// The coefficients of the strain of a field at every node and degree, for
// the orders |m| <= orders, laid out with the degree innermost so that the
// sums over degree read contiguous memory.
class StrainTable {
 public:
  StrainTable(std::size_t nodes, int orders, int lMax)
      : orders_{orders},
        degrees_{static_cast<std::size_t>(lMax) + 1},
        values_(nodes * Components * (2 * orders + 1) * degrees_) {}

  Complex& operator()(std::size_t node, std::size_t component, int m, int l) {
    return values_[Index(node, component, m) + static_cast<std::size_t>(l)];
  }
  const Complex* Row(std::size_t node, std::size_t component, int m) const {
    return values_.data() + Index(node, component, m);
  }
  int Orders() const { return orders_; }
  std::size_t Degrees() const { return degrees_; }

 private:
  int orders_;
  std::size_t degrees_;
  std::vector<Complex> values_;

  std::size_t Index(std::size_t node, std::size_t component, int m) const {
    return ((node * Components + component) * (2 * orders_ + 1) +
            static_cast<std::size_t>(m + orders_)) *
           degrees_;
  }
};

// The generalised Legendre functions X^N_lm at a set of colatitudes, for
// |N| <= 2 and |m| <= orders, with the degree innermost.
class LegendreTable {
 public:
  LegendreTable(int lMax, int orders, const std::vector<double>& theta)
      : orders_{orders},
        degrees_{static_cast<std::size_t>(lMax) + 1},
        values_(theta.size() * 5 * (2 * orders + 1) * degrees_) {
    const auto spin = SpinFunctions(lMax, theta);
    for (auto i : std::views::iota(std::size_t{0}, theta.size())) {
      for (auto l : std::views::iota(0, lMax + 1)) {
        for (auto m :
             std::views::iota(-std::min(l, orders), std::min(l, orders) + 1)) {
          const auto& x = spin(static_cast<Int>(i), l, m);
          for (auto n : std::views::iota(-2, 3)) {
            values_[Index(i, n, m) + static_cast<std::size_t>(l)] = x[n + 2];
          }
        }
      }
    }
  }

  const double* Row(std::size_t point, int n, int m) const {
    return values_.data() + Index(point, n, m);
  }

 private:
  int orders_;
  std::size_t degrees_;
  std::vector<double> values_;

  std::size_t Index(std::size_t point, int n, int m) const {
    return ((point * 5 + static_cast<std::size_t>(n + 2)) * (2 * orders_ + 1) +
            static_cast<std::size_t>(m + orders_)) *
           degrees_;
  }
};

// The spin components of a field's strain at a point and node: the sum over
// degree and order, with the point at colatitude row and longitude phi in the
// field's frame.
std::array<Complex, Components> Evaluate(const StrainTable& table,
                                         std::size_t node,
                                         const LegendreTable& legendre,
                                         std::size_t point, double phi) {
  const auto orders = table.Orders();
  const auto degrees = table.Degrees();
  auto result = std::array<Complex, Components>{};
  for (auto m : std::views::iota(-orders, orders + 1)) {
    const auto phase = std::polar(1.0, m * phi);
    for (auto c : std::views::iota(std::size_t{0}, Components)) {
      const auto* f = table.Row(node, c, m);
      const auto* x = legendre.Row(point, Spin[c], m);
      auto re = 0.0;
      auto im = 0.0;
      for (auto l = std::size_t{0}; l < degrees; ++l) {
        re += f[l].real() * x[l];
        im += f[l].imag() * x[l];
      }
      result[c] += Complex(re, im) * phase;
    }
  }
  return result;
}

// Where a point is in the frames of the two fields, and the angle from the
// forward field's local theta direction to the adjoint field's, by which the
// adjoint field's spin components turn into the forward field's frame.
struct PointGeometry {
  double thetaS, phiS;  // in the frame with the source at the pole
  double thetaR, phiR;  // in the frame with the receiver at the pole
  double turn;
};

// The elastic kernels of the model's own parameters accumulated at a point:
// A, C, F, L, N, qA, qF, qMu.
constexpr auto Elastic = std::size_t{8};

// Radii through every layer of the model, equally spaced from its bottom to
// its top, both included, no further apart than a spacing (m).
std::vector<KernelNode> NodesThrough(const Model& model, double spacing) {
  const auto length = model.GetScales().Length();
  auto nodes = std::vector<KernelNode>{};
  for (auto layer : std::views::iota(Int{0}, model.NumberOfLayers())) {
    const auto lower = model.GetLayer(layer).Lower();
    const auto upper = model.GetLayer(layer).Upper();
    const auto pieces = std::max(
        1,
        static_cast<int>(std::ceil((upper - lower) * length / spacing - 1e-9)));
    for (auto k : std::views::iota(0, pieces + 1)) {
      const auto r = k == pieces ? upper : lower + (upper - lower) * k / pieces;
      nodes.push_back({layer, r * length});
    }
  }
  return nodes;
}

using Vector = Sphere::Vector;

Vector Unit(const Location& location) {
  return Sphere::Unit(Sphere::Colatitude(location.latitude),
                      location.longitude * Sphere::Degree);
}

double Dot(const Vector& a, const Vector& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vector Cross(const Vector& a, const Vector& b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

Vector Normalised(const Vector& a) {
  const auto size = std::sqrt(Dot(a, a));
  return {a[0] / size, a[1] / size, a[2] / size};
}

// The frame of the path from the source to the receiver: the source, the
// direction of travel there, and the normal to the great circle, to the left
// of travel; and the length of the path (degrees).
struct PathFrame {
  Vector source, travel, normal;
  double distance;

  PathFrame(const Location& from, const Location& to)
      : source{Unit(from)}, travel{}, normal{}, distance{} {
    const auto receiver = Unit(to);
    const auto along = Dot(source, receiver);
    travel = Normalised({receiver[0] - along * source[0],
                         receiver[1] - along * source[1],
                         receiver[2] - along * source[2]});
    normal = Cross(source, travel);
    distance = std::atan2(Dot(receiver, travel), along) / Sphere::Degree;
    if (!(distance > 0 && distance < 180)) {
      throw std::invalid_argument(
          "a path needs a source and a receiver neither at one place nor "
          "antipodal");
    }
  }

  // The place at distances along the path and across it (degrees).
  KernelPlace At(double along, double across) const {
    const auto a = along * Sphere::Degree;
    const auto c = across * Sphere::Degree;
    auto x = Vector{};
    for (auto i : {0, 1, 2}) {
      x[i] = std::cos(c) * (source[i] * std::cos(a) + travel[i] * std::sin(a)) +
             std::sin(c) * normal[i];
    }
    return {.colatitude = std::acos(std::clamp(x[2], -1.0, 1.0)),
            .longitude = std::atan2(x[1], x[0]),
            .along = along,
            .across = across};
  }
};

// Values from a to b, both included, no further apart than a spacing.
std::vector<double> Steps(double a, double b, double spacing) {
  const auto pieces =
      std::max(1, static_cast<int>(std::ceil((b - a) / spacing - 1e-9)));
  auto values = std::vector<double>{};
  for (auto k : std::views::iota(0, pieces + 1)) {
    values.push_back(a + (b - a) * k / pieces);
  }
  return values;
}

void CheckSpacing(double spacing, double margin, double halfWidth) {
  if (!(spacing > 0 && margin >= 0 && halfWidth > 0)) {
    throw std::invalid_argument(
        "a section or map needs a positive spacing and width, and a margin "
        "of at least zero");
  }
}

}  // namespace

KernelSampling VolumeSampling(const Model& model, const VolumeGrid& grid) {
  if (grid.degree < 1 || !(grid.maxElement > 0)) {
    throw std::invalid_argument(
        "the grid of volume kernels needs a degree of at least one and "
        "intervals of positive size");
  }
  auto sampling = KernelSampling{.nodes = NodesThrough(model, grid.maxElement),
                                 .places = {}};
  // Gauss-Legendre colatitudes, increasing, equally spaced longitudes from
  // zero, and the two poles.
  const auto rule =
      GaussQuad::GaussLegendreQuadrature1D<double>(grid.degree + 1);
  const auto columns = 2 * grid.degree + 1;
  for (auto x : rule.Points() | std::views::reverse) {
    for (auto j : std::views::iota(0, columns)) {
      sampling.places.push_back(
          {.colatitude = std::acos(x),
           .longitude = 2 * std::numbers::pi * j / columns});
    }
  }
  sampling.places.push_back({.colatitude = 0, .longitude = 0});
  sampling.places.push_back({.colatitude = std::numbers::pi, .longitude = 0});
  return sampling;
}

KernelSampling PathSection(const Model& model, const Location& source,
                           const Location& receiver,
                           const SectionOptions& options) {
  CheckSpacing(options.spacing, options.margin, options.halfWidth);
  const auto path = PathFrame(source, receiver);
  if (path.distance + 2 * options.margin > 360) {
    throw std::invalid_argument("a section cannot go round more than once");
  }
  auto sampling = KernelSampling{
      .nodes = NodesThrough(model, options.radialSpacing), .places = {}};
  for (auto along : Steps(-options.margin, path.distance + options.margin,
                          options.spacing)) {
    sampling.places.push_back(path.At(along, 0));
  }
  return sampling;
}

KernelSampling CrossSection(const Model& model, const Location& source,
                            const Location& receiver,
                            const SectionOptions& options) {
  CheckSpacing(options.spacing, options.margin, options.halfWidth);
  const auto path = PathFrame(source, receiver);
  auto sampling = KernelSampling{
      .nodes = NodesThrough(model, options.radialSpacing), .places = {}};
  for (auto across :
       Steps(-options.halfWidth, options.halfWidth, options.spacing)) {
    sampling.places.push_back(path.At(path.distance / 2, across));
  }
  return sampling;
}

KernelSampling PathMaps(const Model& model, const Location& source,
                        const Location& receiver, const MapOptions& options) {
  CheckSpacing(options.spacing, options.margin, options.halfWidth);
  const auto path = PathFrame(source, receiver);
  const auto& scales = model.GetScales();
  const auto surface = model.SurfaceRadius() * scales.Length();
  auto sampling = KernelSampling{};
  for (auto depth : options.depths) {
    const auto r = surface - depth;
    if (!(depth >= 0 && r > 0)) {
      throw std::invalid_argument("a map's depth is outside the model");
    }
    sampling.nodes.push_back(
        {model.LayerIndex(r / scales.Length(), Side::Above), r});
  }
  if (sampling.nodes.empty()) {
    throw std::invalid_argument("maps need at least one depth");
  }
  for (auto along : Steps(-options.margin, path.distance + options.margin,
                          options.spacing)) {
    for (auto across :
         Steps(-options.halfWidth, options.halfWidth, options.spacing)) {
      sampling.places.push_back(path.At(along, across));
    }
  }
  return sampling;
}

namespace {

// A block of places of one sampling, whose functions are tabulated together.
struct Block {
  std::size_t sampling;
  std::size_t first;
  std::vector<PointGeometry> geometry;
};

}  // namespace

std::vector<std::vector<SampledKernels>> ComputeSampledKernels(
    const Model& model, const Calculation& calculation, const Source& source,
    std::span<const Measurement> measurements, ParameterSet set,
    Scaling scaling, std::span<const KernelSampling> samplings,
    const Progress& progress) {
  const auto& scales = model.GetScales();
  const auto& degrees = calculation.degrees;
  const auto radii =
      CheckKernelInputs(model, calculation, source, measurements);
  const auto length = scales.Length();
  const auto grid = MakeGrid(calculation, scales);
  const auto count = measurements.size();

  // The nodes of every sampling, each solved for once. One at the centre is
  // evaluated just off it, below any start.
  auto points = std::vector<Point>{};
  auto nodeIndex = std::vector<std::vector<std::size_t>>{};
  {
    auto seen = std::map<std::pair<Int, double>, std::size_t>{};
    for (const auto& sampling : samplings) {
      auto& indices = nodeIndex.emplace_back();
      for (const auto& node : sampling.nodes) {
        const auto r = std::max(node.radius / length, 1e-6);
        if (node.layer < 0 || node.layer >= model.NumberOfLayers()) {
          throw std::invalid_argument("a node is in no layer of the model");
        }
        const auto [at, added] =
            seen.try_emplace({node.layer, r}, points.size());
        if (added) points.push_back({node.layer, r});
        indices.push_back(at->second);
      }
    }
  }
  const auto nodeCount = points.size();
  const auto solver = KernelSolver(model, calculation, radii, points);
  const auto adjointIndex = AdjointIndex(radii);

  // The places, in blocks, and where each is seen from the source and from
  // the receiver.
  using Sphere::Degree;
  const auto& receiver = measurements.front().receiver;
  const auto thetaS = Sphere::Colatitude(source.location.latitude);
  const auto phiS = source.location.longitude * Degree;
  const auto thetaR = Sphere::Colatitude(receiver.location.latitude);
  const auto phiR = receiver.location.longitude * Degree;
  constexpr auto blockSize = std::size_t{64};
  auto blocks = std::vector<Block>{};
  for (auto s : std::views::iota(std::size_t{0}, samplings.size())) {
    const auto& places = samplings[s].places;
    for (auto first = std::size_t{0}; first < places.size();
         first += blockSize) {
      auto& block = blocks.emplace_back(Block{s, first, {}});
      const auto last = std::min(first + blockSize, places.size());
      for (auto p : std::views::iota(first, last)) {
        const auto theta = places[p].colatitude;
        const auto phi = places[p].longitude;
        const auto g = GeometryAt(theta, phi, thetaS, phiS);
        const auto r = GeometryAt(theta, phi, thetaR, phiR);
        block.geometry.push_back({.thetaS = g.distance,
                                  .phiS = g.phi,
                                  .thetaR = r.distance,
                                  .phiR = r.phi,
                                  .turn = r.rotation - g.rotation});
      }
    }
  }

  // The adjoint fields of the forces of a basis of the measurements'
  // directions, each at the pole of its own frame, where theta is south and
  // phi east, impulsive; and each measurement's spectrum.
  const auto basis = MakeBasis(measurements);
  auto forces = std::vector<Force>{};
  for (const auto& d : basis.directions) {
    forces.push_back({.r = d.up, .t = -d.north, .p = d.east});
  }
  const auto bases = forces.size();
  auto adjoint = std::vector<std::vector<Complex>>{};
  auto frequencyWeight = std::vector<double>{};
  for (const auto& measurement : measurements) {
    auto spectrum = MakeAdjointSpectrum(model, calculation, measurement);
    adjoint.push_back(std::move(spectrum.force));
    frequencyWeight = std::move(spectrum.weight);
  }

  // The model at the nodes, the source and the receiver.
  auto materials = std::vector<Material>{};
  auto fluid = std::vector<bool>{};
  for (const auto& point : points) {
    materials.push_back(model(point.layer, point.r));
    fluid.push_back(model.GetLayer(point.layer).IsFluid());
  }
  const auto atSource = solver.Source(0);
  const auto atReceiver = solver.Source(adjointIndex);
  const auto moment = source.moment / scales.Moment();
  const auto sourceForce = source.force / scales.Force();

  // The kernels of each measurement on each sampling, accumulated over
  // frequency, node by node, each block behind its own lock.
  auto total = std::vector<std::vector<std::vector<double>>>(count);
  for (auto& kernels : total) {
    for (const auto& sampling : samplings) {
      kernels.emplace_back(sampling.nodes.size() * Elastic *
                           sampling.places.size());
    }
  }
  auto locks = std::vector<std::mutex>(blocks.size());
  auto done = std::atomic<std::size_t>{0};
  auto error = std::exception_ptr{};
  // Once a frequency fails the rest are skipped, and the error is thrown
  // when the loop ends.
  auto failed = std::atomic<bool>{false};
  const auto lMax = degrees.max;

#pragma omp parallel
  {
    auto forward = StrainTable(nodeCount, 2, lMax);
    auto backward = std::vector<StrainTable>{};
    for (auto b = std::size_t{0}; b < bases; ++b) {
      backward.emplace_back(nodeCount, 1, lMax);
    }
#pragma omp for schedule(dynamic)
    for (long index = 0; index < static_cast<long>(grid.Count()); ++index) {
      if (failed) continue;
      try {
        const auto k = static_cast<std::size_t>(index);
        const auto wanted = std::ranges::any_of(
            adjoint, [k](const auto& a) { return a[k] != Complex{}; });
        if (wanted) {
          const auto omega = grid.Omega(k);
          // The strains of the forward field and of each adjoint field,
          // degree by degree.
          for (auto l = degrees.min; l <= lMax; ++l) {
            const auto c = Context(model, l, omega, calculation.attenuation);
            auto responses = Responses{};
            if (l == 0) {
              responses.radial = solver.Radial(omega);
            } else {
              responses.toroidal = solver.Toroidal(l, omega);
              responses.spheroidal = solver.Spheroidal(l, omega);
            }
            const auto fromSource = PoleWeights(
                c, c(atSource.layer, atSource.r), moment, sourceForce);
            auto fromReceiver = std::vector<SourceWeights>{};
            for (const auto& force : forces) {
              // Impulsive: the step's 1 / (i omega) undone.
              auto weights = PoleWeights(c, c(atReceiver.layer, atReceiver.r),
                                         MomentTensor{}, force);
              Scale(weights, TimesI(omega));
              fromReceiver.push_back(weights);
            }
            for (auto n : std::views::iota(std::size_t{0}, nodeCount)) {
              const auto r = points[n].r;
              const auto f = StrainOf(Combine(responses, 0, n, l, fromSource),
                                      l, r, fluid[n]);
              for (auto i : std::views::iota(std::size_t{0}, Components)) {
                for (auto m : OrderArray<Complex>::Orders(l)) {
                  forward(n, i, m, l) = f[i][m];
                }
              }
              for (auto b : std::views::iota(std::size_t{0}, bases)) {
                const auto a = StrainOf(
                    Combine(responses, adjointIndex, n, l, fromReceiver[b]), l,
                    r, fluid[n]);
                for (auto i : std::views::iota(std::size_t{0}, Components)) {
                  for (auto m : OrderArray<Complex>::Orders(std::min(l, 1))) {
                    backward[b](n, i, m, l) = a[i][m];
                  }
                }
              }
            }
          }

          // The factors of each measurement's adjoint field: its spectrum
          // along its direction.
          auto by = std::vector<std::vector<Complex>>(count);
          for (auto j : std::views::iota(std::size_t{0}, count)) {
            for (auto b : std::views::iota(std::size_t{0}, bases)) {
              by[j].push_back(adjoint[j][k] * basis.coordinates[j][b]);
            }
          }
          // The sum over frequency is of the real part weighted, each
          // modulus dispersed by its quality factor.
          const auto weight = frequencyWeight[k];
          const auto dispersion =
              Context(model, 0, omega, calculation.attenuation).Dispersion();
          auto factors = std::vector<std::array<Complex, Elastic>>{};
          for (const auto& m : materials) {
            factors.push_back({-weight * (1.0 + m.qA * dispersion),
                               -weight * (1.0 + m.qA * dispersion),
                               -weight * (1.0 + m.qF * dispersion),
                               -weight * (1.0 + m.qMu * dispersion),
                               -weight * (1.0 + m.qMu * dispersion),
                               -weight * dispersion, -weight * dispersion,
                               -weight * dispersion});
          }

          for (auto step : std::views::iota(std::size_t{0}, blocks.size())) {
            // Threads start at different blocks, to meet less at the locks.
            const auto& block = blocks[(step + k * 7) % blocks.size()];
            const auto& here = block.geometry;
            const auto places = here.size();
            const auto& nodes = nodeIndex[block.sampling];
            auto thetas = [&](auto member) {
              auto values = std::vector<double>{};
              for (const auto& g : here) values.push_back(g.*member);
              return values;
            };
            const auto sourceTable =
                LegendreTable(lMax, 2, thetas(&PointGeometry::thetaS));
            const auto receiverTable =
                LegendreTable(lMax, 1, thetas(&PointGeometry::thetaR));
            // [measurement][node][parameter][place]
            auto local = std::vector<std::vector<double>>(
                count, std::vector<double>(nodes.size() * Elastic * places));
            for (auto n : std::views::iota(std::size_t{0}, nodes.size())) {
              const auto node = nodes[n];
              const auto solid = !fluid[node];
              const auto& m = materials[node];
              const auto& f = factors[node];
              for (auto p : std::views::iota(std::size_t{0}, places)) {
                const auto& g = here[p];
                const auto turn = std::polar(1.0, g.turn);
                const auto turn2 = turn * turn;
                const auto s = Evaluate(forward, node, sourceTable, p, g.phiS);
                // The products with each adjoint field, the coefficients of
                // A, C, F, L and N in the elastic energy.
                auto e = std::vector<std::array<Complex, 5>>{};
                for (auto b : std::views::iota(std::size_t{0}, bases)) {
                  auto a =
                      Evaluate(backward[b], node, receiverTable, p, g.phiR);
                  a[2] *= std::conj(turn);
                  a[3] *= turn;
                  a[4] *= std::conj(turn2);
                  a[5] *= turn2;
                  const auto eA = s[1] * a[1];
                  e.push_back(
                      {eA, s[0] * a[0], s[0] * a[1] + s[1] * a[0],
                       solid ? 2.0 * (s[2] * a[3] + s[3] * a[2]) : Complex{},
                       solid ? (s[4] * a[5] + s[5] * a[4]) / 2.0 - eA
                             : Complex{}});
                }
                for (auto j : std::views::iota(std::size_t{0}, count)) {
                  if (adjoint[j][k] == Complex{}) continue;
                  auto x = std::array<Complex, 5>{};
                  for (auto b : std::views::iota(std::size_t{0}, bases)) {
                    for (auto q :
                         std::views::iota(std::size_t{0}, std::size_t{5})) {
                      x[q] += by[j][b] * e[b][q];
                    }
                  }
                  const auto [eA, eC, eF, eL, eN] = x;
                  const auto products = std::array<Complex, Elastic>{
                      f[0] * eA,       f[1] * eC,
                      f[2] * eF,       f[3] * eL,
                      f[4] * eN,       f[5] * (m.A * eA + m.C * eC),
                      f[6] * m.F * eF, f[7] * (m.L * eL + m.N * eN)};
                  for (auto q : std::views::iota(std::size_t{0}, Elastic)) {
                    local[j][(n * Elastic + q) * places + p] =
                        products[q].real();
                  }
                }
              }
            }
            const auto lock = std::scoped_lock(
                locks[static_cast<std::size_t>(&block - blocks.data())]);
            const auto all = samplings[block.sampling].places.size();
            for (auto j : std::views::iota(std::size_t{0}, count)) {
              auto& out = total[j][block.sampling];
              for (auto n : std::views::iota(std::size_t{0}, nodes.size())) {
                for (auto q : std::views::iota(std::size_t{0}, Elastic)) {
                  auto* to = out.data() + (n * Elastic + q) * all + block.first;
                  const auto* from =
                      local[j].data() + (n * Elastic + q) * places;
                  for (auto p : std::views::iota(std::size_t{0}, places)) {
                    to[p] += from[p];
                  }
                }
              }
            }
          }
        }
        const auto finished = ++done;
        Report(progress, finished, grid.Count());
      } catch (...) {
        failed = true;
#pragma omp critical(YSpecError)
        if (!error) error = std::current_exception();
      }
    }
  }
  if (error) std::rethrow_exception(error);

  // In SI units, and in the parameters asked for, density left out.
  auto names = Names(set);
  names.erase(names.begin());
  const auto volumeUnit = length * length * length;
  const auto modulus = scales.Modulus();
  constexpr auto parameters = static_cast<std::size_t>(RadialParameter::Count);
  auto result = std::vector<std::vector<SampledKernels>>(count);
  for (auto j : std::views::iota(std::size_t{0}, count)) {
    for (auto s : std::views::iota(std::size_t{0}, samplings.size())) {
      auto& kernels = result[j].emplace_back(names, samplings[s]);
      const auto places = samplings[s].places.size();
      for (auto n :
           std::views::iota(std::size_t{0}, samplings[s].nodes.size())) {
        const auto& point = points[nodeIndex[s][n]];
        const auto values = ValuesAt(model, point.layer, point.r);
        for (auto p : std::views::iota(std::size_t{0}, places)) {
          auto k = std::array<double, parameters>{};
          for (auto q : std::views::iota(std::size_t{0}, Elastic)) {
            // After density, A to N are moduli and the rest quality factors.
            const auto unit = q < 5 ? modulus : 1.0;
            k[q + 1] = total[j][s][(n * Elastic + q) * places + p] /
                       (unit * volumeUnit);
          }
          auto [transformed, parameter] = Transform(k, values, set);
          for (auto i : std::views::iota(std::size_t{0}, names.size())) {
            auto value = transformed[i + 1];
            if (scaling == Scaling::Relative) value *= parameter[i + 1];
            kernels.Values(i, n)[p] = value;
          }
        }
      }
    }
  }
  return result;
}

SampledKernels ComputeVolumeKernels(
    const Model& model, const Calculation& calculation, const Source& source,
    const Measurement& measurement, ParameterSet set, Scaling scaling,
    const VolumeGrid& grid, const Progress& progress) {
  const auto sampling = VolumeSampling(model, grid);
  return std::move(ComputeSampledKernels(
                       model, calculation, source, std::span{&measurement, 1},
                       set, scaling, std::span{&sampling, 1}, progress)
                       .front()
                       .front());
}

}  // namespace YSpec
