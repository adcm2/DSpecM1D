#pragma once

// The kinds of seismograms, Displacement and Strain: what each records, as
// sums of the displacement coefficients against the angular functions in the
// frame with the source at the pole, and the turn to the receiver's frame.
// They are kept apart from the computation so that the sums can be tested on
// fields given by hand.

#include <YSpec/Receiver.hpp>
#include <YSpec/Seismograms.hpp>
#include <YSpec/Source.hpp>
#include <YSpec/Synthetics.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <ranges>
#include <span>
#include <vector>

#include "AngularFunctions.hpp"
#include "Coefficients.hpp"
#include "Geometry.hpp"
#include "OrderArray.hpp"

namespace YSpec {

// Where each receiver is seen from each source, and the angular functions
// there (a Table: AngularFunctions or SpinFunctions). These depend only on
// where the receiver is on the surface, so receivers at one place and several
// depths, as in a snapshot, share them.
template <class Table>
class Views {
 public:
  Views(std::span<const Source> sources, std::span<const Receiver> receivers,
        int maxDegree)
      : place_(receivers.size()) {
    using namespace std::complex_literals;
    auto places = std::vector<Location>{};
    for (auto k : std::views::iota(std::size_t{0}, receivers.size())) {
      const auto& location = receivers[k].location;
      const auto found = std::ranges::find_if(places, [&](const Location& p) {
        return p.latitude == location.latitude &&
               p.longitude == location.longitude;
      });
      place_[k] =
          static_cast<std::size_t>(std::distance(places.begin(), found));
      if (found == places.end()) places.push_back(location);
    }
    places_ = places.size();
    auto distances = std::vector<double>{};
    for (const auto& source : sources) {
      for (const auto& place : places) {
        auto view =
            View{.geometry = Geometry(place, source.location), .phase = {}};
        for (auto m : OrderArray<Complex>::Orders(2)) {
          view.phase[m] = std::exp(1i * (m * view.geometry.phi));
        }
        views_.push_back(view);
        distances.push_back(view.geometry.distance);
      }
    }
    table_ = Table(maxDegree, distances);
  }

  // A receiver as seen from a source.
  struct View {
    ReceiverGeometry geometry;
    OrderArray<Complex> phase;  // exp(i m phi)
  };

  const View& Of(std::size_t source, std::size_t receiver) const {
    return views_[Index(source, receiver)];
  }
  const auto& Functions(std::size_t source, std::size_t receiver, int l,
                        int m) const {
    return table_(static_cast<Int>(Index(source, receiver)), l, m);
  }

 private:
  std::vector<std::size_t> place_;  // the place of each receiver
  std::size_t places_ = 0;
  std::vector<View> views_;
  Table table_{0, {}};

  std::size_t Index(std::size_t source, std::size_t receiver) const {
    return source * places_ + place_[receiver];
  }
};

// What seismograms of displacement record: the sums are the displacement in
// the frame with the source at the pole, (r, theta, phi), rotated at the end
// to (up, north, east).
struct Displacement {
  static constexpr std::size_t Components = 3;
  using Table = AngularFunctions;
  using Sums = std::array<Complex, 3>;
  using Output = Seismograms;
  using ComponentType = Component;

  static void Add(Sums& sums, const Coefficients& coefficients, int l, double,
                  bool, const Views<Table>& views, std::size_t source,
                  std::size_t receiver) {
    const auto& phase = views.Of(source, receiver).phase;
    for (auto m : OrderArray<Complex>::Orders(l)) {
      const auto& f = views.Functions(source, receiver, l, m);
      const auto e = phase[m];
      const auto u = coefficients.u[m];
      const auto v = coefficients.v[m];
      const auto w = coefficients.w[m];
      sums[0] += u * f.value * e;
      sums[1] += (v * f.dTheta + TimesI(w) * f.mOverSin) * e;
      sums[2] += (TimesI(v) * f.mOverSin - w * f.dTheta) * e;
    }
  }

  // With rotation the angle that turns (theta, phi) to (south, east).
  static Sums ToReceiver(const Sums& u, double rotation) {
    const auto cosine = std::cos(rotation);
    const auto sine = std::sin(rotation);
    return {u[0], -u[1] * cosine + u[2] * sine, u[1] * sine + u[2] * cosine};
  }

  static double Unit(const Scales& scales, Quantity quantity) {
    return YSpec::Unit(scales, quantity);
  }
};

// What seismograms of strain record. The sums are the strain in the frame
// with the source at the pole, in its spin components: rr, the horizontal
// trace, the components r+ and r- of the vector (rtheta, rphi), and the
// components ++ and -- of the traceless horizontal part, with a+ = atheta -
// i aphi, a- = atheta + i aphi, and t+- = ttheta,theta - tphi,phi -+ 2i
// ttheta,phi. With
//
//   u = U Y r + V grad1 Y - W r x grad1 Y,
//
// the spin components of the displacement are u+- = -+zeta (V +- iW) X^+-1,
// and those of the strain
//
//   rr = U' X^0,   trace = (2U - zeta^2 V) X^0 / r,
//   r+- = -+(zeta/2) (V' - V/r + U/r +- i (W' - W/r)) X^+-1,
//   +-+- = zeta sqrt((l-1)(l+2)) (V +- iW) X^+-2 / r,
//
// each times exp(i m phi). They are turned at the end into the components
// (rr, theta theta, phi phi, r theta, r phi, theta phi), and rotated to
// (south, east).
struct Strain {
  static constexpr std::size_t Components = 6;
  using Table = SpinFunctions;
  using Sums = std::array<Complex, 6>;
  using Output = Strains;
  using ComponentType = StrainComponent;

  // In a fluid the radial derivative of the tangential displacement is not
  // given, and r+ and r- are left zero.
  static void Add(Sums& sums, const Coefficients& coefficients, int l, double r,
                  bool fluid, const Views<Table>& views, std::size_t source,
                  std::size_t receiver) {
    const auto& phase = views.Of(source, receiver).phase;
    const auto zeta2 = l * (l + 1.0);
    const auto zeta = std::sqrt(zeta2);
    const auto twist = l >= 2 ? zeta * std::sqrt((l - 1.0) * (l + 2)) / r : 0.0;
    for (auto m : OrderArray<Complex>::Orders(l)) {
      // X^-2, ..., X^2.
      const auto& x = views.Functions(source, receiver, l, m);
      const auto e = phase[m];
      const auto u = coefficients.u[m];
      const auto v = coefficients.v[m];
      const auto w = coefficients.w[m];
      const auto shear = coefficients.dv[m] - (v - u) / r;
      const auto turn = coefficients.dw[m] - w / r;
      sums[0] += coefficients.du[m] * x[2] * e;
      sums[1] += (2.0 * u - zeta2 * v) / r * x[2] * e;
      if (!fluid) {
        sums[2] += -zeta / 2 * (shear + TimesI(turn)) * x[3] * e;
        sums[3] += zeta / 2 * (shear - TimesI(turn)) * x[1] * e;
      }
      sums[4] += twist * (v + TimesI(w)) * x[4] * e;
      sums[5] += twist * (v - TimesI(w)) * x[0] * e;
    }
  }

  static Sums ToReceiver(const Sums& spin, double rotation) {
    // The components in the frame with the source at the pole.
    const auto rr = spin[0];
    const auto trace = spin[1];
    const auto rTheta = (spin[2] + spin[3]) / 2.0;
    const auto rPhi = TimesI(spin[2] - spin[3]) / 2.0;
    const auto difference = (spin[4] + spin[5]) / 2.0;
    const auto thetaPhi = TimesI(spin[4] - spin[5]) / 4.0;
    const auto thetaTheta = (trace + difference) / 2.0;
    const auto phiPhi = (trace - difference) / 2.0;
    // Rotated by R, which turns (theta, phi) to (south, east): a vector to
    // R a and a tensor to R t R^T.
    const auto c = std::cos(rotation);
    const auto s = std::sin(rotation);
    const auto south = c * rTheta - s * rPhi;
    const auto east = s * rTheta + c * rPhi;
    const auto ss = c * c * thetaTheta - 2 * c * s * thetaPhi + s * s * phiPhi;
    const auto ee = s * s * thetaTheta + 2 * c * s * thetaPhi + c * c * phiPhi;
    const auto se = c * s * (thetaTheta - phiPhi) + (c * c - s * s) * thetaPhi;
    return {rr, ss, ee, south, east, se};
  }

  // Strain is displacement over length.
  static double Unit(const Scales& scales, Quantity quantity) {
    return YSpec::Unit(scales, quantity) / scales.Length();
  }
};

}  // namespace YSpec
