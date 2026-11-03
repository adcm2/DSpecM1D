#include <gtest/gtest.h>

#include <GSHTrans/Wigner.hpp>
#include <GaussQuad/All>
#include <YSpec/Kernels.hpp>
#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
#include <string>
#include <vector>

#include "Coefficients.hpp"
#include "Geometry.hpp"
#include "KernelParts.hpp"
#include "SourceWeights.hpp"
#include "Sphere.hpp"
#include "TimeSeries.hpp"

namespace {

using namespace YSpec;
using P = RadialParameter;

const auto kPrem = std::string(YSPEC_EXAMPLES_DIR) + "/prem.200";

Calculation Small() {
  auto c = Calculation{};
  c.degrees = {0, 16};
  c.band = {0.2e-3, 3e-3};
  c.sampling = {1800, 60};
  c.taper = Taper{0.3e-3, 0.6e-3, 2.5e-3, 3e-3};
  c.quantity = Quantity::Velocity;
  c.solver.tolerance = 1e-9;
  return c;
}

// The integral over each sphere of the kernels at a node is the radial
// kernel there. The grid resolves twice the largest degree, the band of the
// products of the two fields, so its quadrature is exact; the transverse
// component makes the turn between the fields' frames matter.
TEST(VolumeKernels, IntegrateOverSpheresToTheRadialKernels) {
  const auto model = Model::ReadDeck(kPrem);
  const auto calculation = Small();
  const auto source = Source{{10, 20},
                             30e3,
                             MomentTensor{.rr = 1e19,
                                          .rt = 2e19,
                                          .rp = -1e19,
                                          .tt = -3e18,
                                          .tp = 5e18,
                                          .pp = -7e18},
                             Force{}};
  const auto receiver = Receiver{{-15, 50}, 3e3};
  for (auto direction :
       {TransverseDirection(source.location, receiver.location),
        Along(Vertical)}) {
    const auto measurement =
        Traveltime(model, calculation, source, receiver, direction,
                   {.start = 600, .end = 1500, .taper = 100});
    const auto degree = 2 * calculation.degrees.max;
    const auto volume = ComputeVolumeKernels(
        model, calculation, source, measurement, ParameterSet::Love,
        Scaling::Absolute, {.degree = degree, .maxElement = 500e3});

    // The radial kernels at the same nodes.
    const auto& scales = model.GetScales();
    const auto length = scales.Length();
    auto points = std::vector<Point>{};
    for (const auto& node : volume.Nodes()) {
      points.push_back({node.layer, std::max(node.radius / length, 1e-6)});
    }
    const auto radial = AccumulateRadial(model, calculation, source,
                                         std::span{&measurement, 1}, points)
                            .front();

    // The quadrature over the sphere: Gauss-Legendre in cos(theta), and the
    // trapezium rule in phi.
    const auto rule = GaussQuad::GaussLegendreQuadrature1D<double>(degree + 1);
    // The grid's places, colatitude by colatitude, then the poles.
    const auto columns = static_cast<std::size_t>(2 * degree + 1);
    const auto rows = static_cast<std::size_t>(degree + 1);
    ASSERT_EQ(volume.Places().size(), rows * columns + 2);
    ASSERT_EQ(rows, static_cast<std::size_t>(degree + 1));
    // The names of the Love set, density first.
    for (auto j = std::size_t{0}; j < volume.Names().size(); ++j) {
      const auto parameter = static_cast<P>(j + 1);
      const auto unit = j < 5 ? scales.Modulus() : 1.0;
      auto largest = 0.0;
      auto worst = 0.0;
      for (auto n = std::size_t{0}; n < points.size(); ++n) {
        const auto values = volume.Values(j, n);
        auto integral = 0.0;
        for (auto i = std::size_t{0}; i < rows; ++i) {
          // The colatitudes increase, the rule's points in cos(theta) too.
          const auto w = rule.Weights()[rows - 1 - i];
          for (auto c = std::size_t{0}; c < columns; ++c) {
            integral += w * 2 * std::numbers::pi /
                        static_cast<double>(columns) * values[i * columns + c];
          }
        }
        const auto expected =
            radial.points[n][static_cast<std::size_t>(parameter)] /
            (unit * length * length * length);
        largest = std::max(largest, std::abs(expected));
        worst = std::max(worst, std::abs(integral - expected));
      }
      ASSERT_GT(largest, 0) << volume.Names()[j];
      EXPECT_LT(worst, 1e-8 * largest) << volume.Names()[j];
    }
  }
}

// For an explosion and a vertical component, with the source and receiver
// on the equator, the kernels are symmetric about it.
TEST(VolumeKernels, AreSymmetricAboutThePlaneOfSourceAndReceiver) {
  const auto model = Model::ReadDeck(kPrem);
  const auto calculation = Small();
  const auto source = Source{
      {0, 0}, 30e3, MomentTensor{.rr = 1e19, .tt = 1e19, .pp = 1e19}, Force{}};
  const auto receiver = Receiver{{0, 60}, 3e3};
  const auto measurement =
      Traveltime(model, calculation, source, receiver, Along(Vertical),
                 {.start = 600, .end = 1500, .taper = 100});
  const auto volume = ComputeVolumeKernels(
      model, calculation, source, measurement, ParameterSet::Velocities,
      Scaling::Relative, {.degree = 20, .maxElement = 1000e3});
  const auto columns = std::size_t{41};
  const auto rows = std::size_t{21};
  // Against the largest of each kernel anywhere: at the sea surface, where
  // the pressure vanishes, the kernels of vp and qkappa are rounding errors.
  for (auto j = std::size_t{0}; j < volume.Names().size(); ++j) {
    auto largest = 0.0;
    for (auto n = std::size_t{0}; n < volume.Nodes().size(); ++n) {
      for (auto v : volume.Values(j, n))
        largest = std::max(largest, std::abs(v));
    }
    ASSERT_GT(largest, 0) << volume.Names()[j];
    for (auto n = std::size_t{0}; n < volume.Nodes().size(); ++n) {
      const auto values = volume.Values(j, n);
      for (auto i = std::size_t{0}; i < rows; ++i) {
        for (auto c = std::size_t{0}; c < columns; ++c) {
          EXPECT_NEAR(values[i * columns + c],
                      values[(rows - 1 - i) * columns + c], 1e-9 * largest)
              << volume.Names()[j] << ", node " << n;
        }
      }
      // And the poles alike.
      EXPECT_NEAR(values[rows * columns], values[rows * columns + 1],
                  1e-9 * largest);
    }
  }
}

// At a few points, against the kernels computed in the source's frame alone:
// the adjoint field there has every order, its force being away from the
// pole, so the two fields meet without a turn between frames.
TEST(VolumeKernels, AgreeWithTheAdjointFieldInTheSourcesFrame) {
  using namespace std::complex_literals;
  const auto model = Model::ReadDeck(kPrem);
  auto calculation = Small();
  calculation.degrees = {0, 10};
  const auto source = Source{{10, 20},
                             30e3,
                             MomentTensor{.rr = 1e19,
                                          .rt = 2e19,
                                          .rp = -1e19,
                                          .tt = -3e18,
                                          .tp = 5e18,
                                          .pp = -7e18},
                             Force{}};
  const auto receiver = Receiver{{-15, 50}, 3e3};
  const auto direction =
      TransverseDirection(source.location, receiver.location);
  const auto measurement =
      Traveltime(model, calculation, source, receiver, direction,
                 {.start = 600, .end = 1500, .taper = 100});
  const auto volume = ComputeVolumeKernels(
      model, calculation, source, measurement, ParameterSet::Love,
      Scaling::Absolute, {.degree = 12, .maxElement = 1500e3});

  // Some nodes in the mantle and some points of the grid.
  const auto& scales = model.GetScales();
  const auto length = scales.Length();
  auto nodes = std::vector<std::size_t>{};
  for (auto n = std::size_t{0}; n < volume.Nodes().size(); ++n) {
    const auto r = volume.Nodes()[n].radius;
    if (r > 3600e3 && r < 6300e3) nodes.push_back(n);
  }
  ASSERT_GE(nodes.size(), 2u);
  const auto columns = std::size_t{25};
  const auto picked = std::vector<std::size_t>{3, 40, 77, 120, 151, 200};
  auto thetas = std::vector<double>{};
  auto phis = std::vector<double>{};
  using Sphere::Degree;
  const auto thetaS = Sphere::Colatitude(source.location.latitude);
  const auto phiS = source.location.longitude * Degree;
  for (auto p : picked) {
    const auto& place = volume.Places()[p];
    ASSERT_EQ(p % columns,
              static_cast<std::size_t>(std::round(place.longitude * columns /
                                                  (2 * std::numbers::pi))));
    const auto g = GeometryAt(place.colatitude, place.longitude, thetaS, phiS);
    thetas.push_back(g.distance);
    phis.push_back(g.phi);
  }

  // The solves, the force in the source's frame, and the functions there.
  const auto lMax = calculation.degrees.max;
  const auto rs = model.SurfaceRadius() - source.depth / length;
  const auto rr = model.SurfaceRadius() - receiver.depth / length;
  auto points = std::vector<Point>{};
  for (auto n : nodes) {
    points.push_back(
        {volume.Nodes()[n].layer, volume.Nodes()[n].radius / length});
  }
  const auto solver = RadialSolver(model, {rs, rr}, points, calculation.gravity,
                                   calculation.attenuation, calculation.solver);
  const auto at = Geometry(receiver.location, source.location);
  const auto f = InSourceFrame(direction, at.rotation);
  const auto plus = Complex(f.t, -f.p);
  const auto minus = Complex(f.t, f.p);
  using Table = GSHTrans::Wigner<double, GSHTrans::All, GSHTrans::All,
                                 GSHTrans::Multiple>;
  const auto atReceiver = Table(lMax, lMax, 1, std::vector{at.distance});
  const auto atPoints = Table(lMax, lMax, 2, thetas);
  const auto [adjoint, weight] =
      MakeAdjointSpectrum(model, calculation, measurement);
  const auto grid = MakeGrid(calculation, scales);
  const auto moment = source.moment / scales.Moment();

  // The kernels A, C, F, L, N, qA, qF, qMu at each node and point.
  auto kernels =
      std::vector<std::array<double, 8>>(nodes.size() * picked.size());
  for (auto k = std::size_t{0}; k < grid.Count(); ++k) {
    if (adjoint[k] == Complex{}) continue;
    const auto omega = grid.Omega(k);
    // The spin components of both strains, rr, h, r+, r-, ++, --.
    using Strain = std::array<Complex, 6>;
    auto forward = std::vector<Strain>(kernels.size());
    auto backward = std::vector<Strain>(kernels.size());
    for (auto l = 0; l <= lMax; ++l) {
      const auto c = Context(model, l, omega, calculation.attenuation);
      const auto toroidal = l > 0 ? solver.Toroidal(l, omega)
                                  : ResponseTable<ToroidalResponse>(0, 0);
      const auto spheroidal = l > 0 ? solver.Spheroidal(l, omega)
                                    : ResponseTable<SpheroidalResponse>(0, 0);
      const auto radial =
          l == 0 ? solver.Radial(omega) : ResponseTable<RadialResponse>(0, 0);
      const auto weights =
          PoleWeights(c, c(solver.Source(0).layer, rs), moment, Force{});
      const auto zeta2 = l * (l + 1.0);
      const auto zeta = std::sqrt(zeta2);
      for (auto i = std::size_t{0}; i < nodes.size(); ++i) {
        const auto r = points[i].r;
        const auto twist =
            l >= 2 ? zeta * std::sqrt((l - 1.0) * (l + 2)) / r : 0.0;
        // The strain of displacement coefficients, as a sum over the points.
        auto add = [&](std::vector<Strain>& out, int m, Complex u, Complex v,
                       Complex w, Complex du, Complex dv, Complex dw) {
          const auto shear = dv - (v - u) / r;
          const auto turn = dw - w / r;
          const auto strain = Strain{du,
                                     (2.0 * u - zeta2 * v) / r,
                                     -zeta / 2 * (shear + 1i * turn),
                                     zeta / 2 * (shear - 1i * turn),
                                     twist * (v + 1i * w),
                                     twist * (v - 1i * w)};
          constexpr auto spin = std::array{0, 0, 1, -1, 2, -2};
          for (auto q = std::size_t{0}; q < picked.size(); ++q) {
            const auto e = std::polar(1.0, m * phis[q]);
            for (auto s = std::size_t{0}; s < 6; ++s) {
              if (std::abs(spin[s]) > l) continue;
              out[i * picked.size() + q][s] +=
                  strain[s] * atPoints[spin[s], q][l][m] * e;
            }
          }
        };
        // The forward field, of orders |m| <= 2.
        for (auto m : OrderArray<Complex>::Orders(l)) {
          auto u = Complex{}, v = Complex{}, w = Complex{}, du = Complex{},
               dv = Complex{}, dw = Complex{};
          if (l == 0) {
            u = weights.radial[0] * radial(0, i).U[0] +
                weights.radial[1] * radial(0, i).U[1];
            du = weights.radial[0] * radial(0, i).DU[0] +
                 weights.radial[1] * radial(0, i).DU[1];
          } else {
            for (auto j : {0, 1, 2, 3}) {
              const auto x = weights.spheroidal[m][j];
              u += x * spheroidal(0, i).U[j];
              v += x * spheroidal(0, i).V[j];
              du += x * spheroidal(0, i).DU[j];
              dv += x * spheroidal(0, i).DV[j];
            }
            for (auto j : {0, 1}) {
              w += weights.toroidal[m][j] * toroidal(0, i).W[j];
              dw += weights.toroidal[m][j] * toroidal(0, i).DW[j];
            }
          }
          add(forward, m, u, v, w, du, dv, dw);
        }
        // The adjoint field, of every order: the force's weights in the
        // jumps of P, S and T, impulsive, with the spectrum A.
        for (auto m = -l; m <= l; ++m) {
          const auto e = std::polar(1.0, -m * at.phi) / (rr * rr) * adjoint[k];
          // The views of spin 1 and -1 start at degree 1.
          const auto xm = l > 0 ? atReceiver[-1, 0][l][m] : 0.0;
          const auto x0 = atReceiver[0, 0][l][m];
          const auto xp = l > 0 ? atReceiver[1, 0][l][m] : 0.0;
          const auto p = -rr * f.r * x0 * e;
          const auto s = -0.5 * rr * (minus * xm - plus * xp) * e;
          const auto t = -0.5i * rr * (plus * xp + minus * xm) * e;
          auto u = Complex{}, v = Complex{}, w = Complex{}, du = Complex{},
               dv = Complex{}, dw = Complex{};
          if (l == 0) {
            u = p * radial(1, i).U[1];
            du = p * radial(1, i).DU[1];
          } else {
            u = p * spheroidal(1, i).U[2] + s * spheroidal(1, i).U[3];
            v = p * spheroidal(1, i).V[2] + s * spheroidal(1, i).V[3];
            du = p * spheroidal(1, i).DU[2] + s * spheroidal(1, i).DU[3];
            dv = p * spheroidal(1, i).DV[2] + s * spheroidal(1, i).DV[3];
            w = t * toroidal(1, i).W[1];
            dw = t * toroidal(1, i).DW[1];
          }
          add(backward, m, u, v, w, du, dv, dw);
        }
      }
    }
    const auto dispersion =
        Context(model, 0, omega, calculation.attenuation).Dispersion();
    for (auto i = std::size_t{0}; i < nodes.size(); ++i) {
      const auto m = model(points[i].layer, points[i].r);
      for (auto q = std::size_t{0}; q < picked.size(); ++q) {
        const auto& a = forward[i * picked.size() + q];
        const auto& b = backward[i * picked.size() + q];
        const auto eC = a[0] * b[0];
        const auto eA = a[1] * b[1];
        const auto eF = a[0] * b[1] + a[1] * b[0];
        const auto eL = 2.0 * (a[2] * b[3] + a[3] * b[2]);
        const auto eN = (a[4] * b[5] + a[5] * b[4]) / 2.0 - eA;
        const auto w = -weight[k];
        auto& out = kernels[i * picked.size() + q];
        out[0] += (w * (1.0 + m.qA * dispersion) * eA).real();
        out[1] += (w * (1.0 + m.qA * dispersion) * eC).real();
        out[2] += (w * (1.0 + m.qF * dispersion) * eF).real();
        out[3] += (w * (1.0 + m.qMu * dispersion) * eL).real();
        out[4] += (w * (1.0 + m.qMu * dispersion) * eN).real();
        out[5] += (w * dispersion * (m.A * eA + m.C * eC)).real();
        out[6] += (w * dispersion * m.F * eF).real();
        out[7] += (w * dispersion * (m.L * eL + m.N * eN)).real();
      }
    }
  }

  // In SI units, against the volume's.
  const auto cube = length * length * length;
  for (auto j = std::size_t{0}; j < 8; ++j) {
    const auto unit = (j < 5 ? scales.Modulus() : 1.0) * cube;
    auto largest = 0.0;
    for (const auto& k : kernels) largest = std::max(largest, std::abs(k[j]));
    ASSERT_GT(largest, 0) << volume.Names()[j];
    for (auto i = std::size_t{0}; i < nodes.size(); ++i) {
      for (auto q = std::size_t{0}; q < picked.size(); ++q) {
        EXPECT_NEAR(volume.Values(j, nodes[i])[picked[q]] * unit,
                    kernels[i * picked.size() + q][j], 1e-8 * largest)
            << volume.Names()[j] << ", node " << nodes[i] << ", point "
            << picked[q];
      }
    }
  }
}

// Kernels computed together, sharing the forward field and the adjoint
// fields of a basis of the directions, are those computed one by one.
TEST(VolumeKernels, AreTheSameComputedTogether) {
  const auto model = Model::ReadDeck(kPrem);
  auto calculation = Small();
  calculation.degrees = {0, 10};
  const auto source = Source{{10, 20},
                             30e3,
                             MomentTensor{.rr = 1e19,
                                          .rt = 2e19,
                                          .rp = -1e19,
                                          .tt = -3e18,
                                          .tp = 5e18,
                                          .pp = -7e18},
                             Force{}};
  const auto receiver = Receiver{{-15, 50}, 3e3};
  const auto window = Window{.start = 600, .end = 1500, .taper = 100};
  // Three measurements in two directions.
  const auto measurements = std::vector<Measurement>{
      Traveltime(model, calculation, source, receiver, Along(Vertical), window),
      Traveltime(model, calculation, source, receiver,
                 TransverseDirection(source.location, receiver.location),
                 {.start = 700, .end = 1200, .taper = 100}),
      Traveltime(model, calculation, source, receiver, Along(Vertical),
                 {.start = 900, .end = 1500, .taper = 100})};
  const auto samplings = std::vector<KernelSampling>{
      PathSection(model, source.location, receiver.location,
                  {.margin = 5, .spacing = 5, .radialSpacing = 1000e3}),
      PathMaps(
          model, source.location, receiver.location,
          {.depths = {200e3}, .margin = 5, .halfWidth = 10, .spacing = 5})};
  const auto together = ComputeSampledKernels(
      model, calculation, source, measurements, ParameterSet::Velocities,
      Scaling::Relative, samplings);
  const auto radial =
      ComputeRadialKernels(model, calculation, source, measurements,
                           {.order = 3, .maxElement = 500e3});
  for (auto j = std::size_t{0}; j < measurements.size(); ++j) {
    const auto alone = ComputeSampledKernels(
        model, calculation, source, std::span{&measurements[j], 1},
        ParameterSet::Velocities, Scaling::Relative, samplings);
    for (auto s = std::size_t{0}; s < samplings.size(); ++s) {
      const auto& a = together[j][s];
      const auto& b = alone[0][s];
      for (auto i = std::size_t{0}; i < a.Names().size(); ++i) {
        // Against the largest of the kernel: where it vanishes, as vp's does
        // at the sea surface, what is left is rounding.
        auto largest = 0.0;
        for (auto n = std::size_t{0}; n < b.Nodes().size(); ++n) {
          for (auto v : b.Values(i, n))
            largest = std::max(largest, std::abs(v));
        }
        for (auto n = std::size_t{0}; n < a.Nodes().size(); ++n) {
          const auto x = a.Values(i, n);
          const auto y = b.Values(i, n);
          for (auto p = std::size_t{0}; p < x.size(); ++p) {
            EXPECT_NEAR(x[p], y[p], 1e-10 * largest)
                << "measurement " << j << ", " << a.Names()[i];
          }
        }
      }
    }
    const auto one =
        ComputeRadialKernels(model, calculation, source, measurements[j],
                             {.order = 3, .maxElement = 500e3});
    for (auto p = std::size_t{0}; p < static_cast<std::size_t>(P::Count); ++p) {
      const auto parameter = static_cast<P>(p);
      for (auto i = std::size_t{0}; i < one.Size(); ++i) {
        EXPECT_NEAR(radial[j][parameter][i], one[parameter][i],
                    1e-10 * std::abs(one[parameter][i]) + 1e-300);
      }
      EXPECT_NEAR(radial[j].AtSource(parameter), one.AtSource(parameter),
                  1e-10 * std::abs(one.AtSource(parameter)) + 1e-300);
    }
  }
}

// The sections and maps lie where they should: the path through the source
// and the receiver, the section across it at its middle.
TEST(VolumeKernels, LayOutSectionsAndMapsByThePath) {
  const auto model = Model::ReadDeck(kPrem);
  const auto source = Location{10, 20};
  const auto receiver = Location{-15, 50};
  const auto distance = [](const KernelPlace& a, double theta, double phi) {
    return GeometryAt(a.colatitude, a.longitude, theta, phi).distance /
           Sphere::Degree;
  };
  const auto thetaS = Sphere::Colatitude(source.latitude);
  const auto phiS = source.longitude * Sphere::Degree;
  const auto thetaR = Sphere::Colatitude(receiver.latitude);
  const auto phiR = receiver.longitude * Sphere::Degree;
  const auto delta =
      GeometryAt(thetaR, phiR, thetaS, phiS).distance / Sphere::Degree;
  const auto path =
      PathSection(model, source, receiver, {.margin = 10, .spacing = 0.5});
  for (const auto& place : path.places) {
    EXPECT_NEAR(distance(place, thetaS, phiS), std::abs(place.along), 1e-9);
    EXPECT_NEAR(distance(place, thetaR, phiR), std::abs(place.along - delta),
                1e-9);
    EXPECT_EQ(place.across, 0);
  }
  EXPECT_NEAR(path.places.front().along, -10, 1e-12);
  EXPECT_NEAR(path.places.back().along, delta + 10, 1e-12);
  // Radii through every layer, both ends of each.
  EXPECT_EQ(path.nodes.front().radius, 0);
  EXPECT_EQ(path.nodes.back().layer, model.NumberOfLayers() - 1);

  const auto across =
      CrossSection(model, source, receiver, {.halfWidth = 20, .spacing = 1});
  for (const auto& place : across.places) {
    // Equally far from source and receiver, at the distance from the
    // middle of the path that a right angle there gives.
    const auto toSource = distance(place, thetaS, phiS);
    EXPECT_NEAR(toSource, distance(place, thetaR, phiR), 1e-9);
    const auto c = std::cos(place.across * Sphere::Degree) *
                   std::cos(delta / 2 * Sphere::Degree);
    EXPECT_NEAR(toSource, std::acos(c) / Sphere::Degree, 1e-9);
  }

  const auto maps = PathMaps(model, source, receiver,
                             {.depths = {100e3, 670e3}, .spacing = 2});
  ASSERT_EQ(maps.nodes.size(), 2u);
  // At a discontinuity, the layer above.
  const auto length = model.GetScales().Length();
  EXPECT_EQ(maps.nodes[1].layer,
            model.LayerIndex((maps.nodes[1].radius) / length, Side::Above));
  EXPECT_THROW((void)PathSection(model, source, source), std::invalid_argument);
}

}  // namespace
