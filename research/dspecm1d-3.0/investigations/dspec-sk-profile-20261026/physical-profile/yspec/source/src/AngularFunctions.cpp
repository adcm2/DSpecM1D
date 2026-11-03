#include "AngularFunctions.hpp"

#include <GSHTrans/Wigner.hpp>
#include <algorithm>
#include <cmath>
#include <ranges>
#include "../../profile.hpp"

namespace YSpec {

AngularFunctions::AngularFunctions(int lMax, std::vector<double> theta)
    : lMax_{lMax}, count_{std::ssize(theta)} {
  values_.assign(static_cast<std::size_t>(count_) * (lMax + 1) * 5, {});
  // GSHTrans refuses a double table above a degree at which a recursion for
  // high orders would underflow; long double reaches far higher.
  if (lMax <= GSHTrans::MaxSafeDegree<double>()) {
    Fill<double>(theta);
  } else {
    Fill<long double>(theta);
  }
}

template <class Real>
void AngularFunctions::Fill(const std::vector<double>& theta) {
  using Table =
      GSHTrans::Wigner<Real, GSHTrans::All, GSHTrans::All, GSHTrans::Multiple>;
  const auto nMax = lMax_ >= 1 ? 1 : 0;
  const auto angles = std::vector<Real>(theta.begin(), theta.end());
  const auto wigner = Table(lMax_, 2, nMax, angles);
  for (auto i : std::views::iota(Int{0}, count_)) {
    for (auto l : std::views::iota(0, lMax_ + 1)) {
      const auto zeta = std::sqrt(static_cast<double>(l) * (l + 1));
      const auto mMax = std::min(l, 2);
      for (auto m : std::views::iota(-mMax, mMax + 1)) {
        auto& functions = Value(i, l, m);
        functions.value = static_cast<double>(wigner[0, i][l][m]);
        if (l == 0) continue;
        const auto minus = static_cast<double>(wigner[-1, i][l][m]);
        const auto plus = static_cast<double>(wigner[1, i][l][m]);
        functions.dTheta = zeta / 2 * (minus - plus);
        functions.mOverSin = -zeta / 2 * (plus + minus);
      }
    }
  }
}

SpinFunctions::SpinFunctions(int lMax, std::vector<double> theta)
    : lMax_{lMax}, count_{std::ssize(theta)} {
  values_.assign(static_cast<std::size_t>(count_) * (lMax + 1) * 5, {});
  if (lMax <= GSHTrans::MaxSafeDegree<double>()) {
    Fill<double>(theta);
  } else {
    Fill<long double>(theta);
  }
}

template <class Real>
void SpinFunctions::Fill(const std::vector<double>& theta) {
  using Table =
      GSHTrans::Wigner<Real, GSHTrans::All, GSHTrans::All, GSHTrans::Multiple>;
  const auto nMax = std::min(lMax_, 2);
  const auto angles = std::vector<Real>(theta.begin(), theta.end());
  const auto wigner = [&] {
    mathprofile::Scope scope("angular.wigner_function_table_build");
    return Table(lMax_, 2, nMax, angles);
  }();
  mathprofile::count("wigner_table_angles", theta.size());
  mathprofile::count("wigner_table_max_degree", lMax_);
  mathprofile::Scope copyScope("angular.spin_legendre_table_population");
  for (auto i : std::views::iota(Int{0}, count_)) {
    for (auto l : std::views::iota(0, lMax_ + 1)) {
      const auto mMax = std::min(l, 2);
      for (auto m : std::views::iota(-mMax, mMax + 1)) {
        auto& values = values_[Index(i, l, m)];
        for (auto n :
             std::views::iota(-std::min(l, nMax), std::min(l, nMax) + 1)) {
          values[n + 2] = static_cast<double>(wigner[n, i][l][m]);
        }
      }
    }
  }
}

}  // namespace YSpec
