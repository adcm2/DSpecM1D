#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <numbers>
#include <vector>

#include "CheckAdditionTheorem.hpp"
#include "CheckLegendre.hpp"
#include "CheckWignerBoundary.hpp"
#include "CheckWignerConvention.hpp"

// The Wigner d-function tables: the value convention, the boundary recursions
// against their closed forms, the upper-index-zero row against
// std::sph_legendre, the addition theorem, and the transform-major and
// reflected layouts against the block layout. The checks themselves live in
// the Check*.hpp headers.

namespace {

void CheckSingleUpperIndexAccess(std::ptrdiff_t n) {
  using namespace GSHTrans;

  constexpr std::ptrdiff_t lMax = 5;
  constexpr std::ptrdiff_t mMax = 3;

  auto singleAngle = Wigner<double, All, Single, Single>(lMax, mMax, n, 0.7);
  auto explicitSingleAngle = singleAngle[n, 0];
  for (auto l : singleAngle.Degrees()) {
    for (auto m : explicitSingleAngle[l].Orders()) {
      EXPECT_DOUBLE_EQ(singleAngle[l][m], explicitSingleAngle[l][m]);
    }
  }

  const auto angles = std::array{0.2, 0.7, 1.3};
  auto multipleAngles =
      Wigner<double, All, Single, Multiple>(lMax, mMax, n, angles);
  for (auto iTheta : multipleAngles.AngleIndices()) {
    auto implicitView = multipleAngles[iTheta];
    auto explicitView = multipleAngles[n, iTheta];
    for (auto l : explicitView.Degrees()) {
      for (auto m : explicitView[l].Orders()) {
        EXPECT_DOUBLE_EQ(implicitView[l][m], explicitView[l][m]);
      }
    }
  }
}

}  // namespace

// Pin the value convention: stored values are sqrt((2l+1)/(4 pi)) d^l_{Nm},
// with the upper index first, per Dahlen & Tromp (1998) eq. (C.115).
TEST(Wigner, CheckConventionDouble) {
  EXPECT_EQ(CheckWignerConvention<double>(), 0);
}

TEST(Wigner, CheckConventionLongDouble) {
  EXPECT_EQ(CheckWignerConvention<long double>(), 0);
}

// The seed row and the boundary orders come from recursions; the closed
// forms in WignerDetails are the definition they answer to.
TEST(Wigner, CheckBoundaryRecursionDouble) {
  EXPECT_LT(CheckWignerBoundary<double>(),
            CheckWignerBoundaryTolerance<double>());
}

TEST(Wigner, CheckBoundaryRecursionLongDouble) {
  EXPECT_LT(CheckWignerBoundary<long double>(),
            CheckWignerBoundaryTolerance<long double>());
}

// Compare values for n = 0 to the std library function.
TEST(Wigner, CheckLegendreDouble) {
  int i = CheckLegendre<double>();
  EXPECT_EQ(i, 0);
}

TEST(Wigner, CheckLegendreLongDouble) {
  int i = CheckLegendre<long double>();
  EXPECT_EQ(i, 0);
}

// Check the addition theorem is satisfied.
TEST(Wigner, CheckAdditionTheoremDouble) {
  int i = CheckAdditionTheorem<double>();
  EXPECT_EQ(i, 0);
}

TEST(Wigner, CheckAdditionTheoremLongDouble) {
  int i = CheckAdditionTheorem<long double>();
  EXPECT_EQ(i, 0);
}

TEST(Wigner, SinglePositiveUpperIndexAccess) { CheckSingleUpperIndexAccess(2); }

TEST(Wigner, SingleNegativeUpperIndexAccess) {
  CheckSingleUpperIndexAccess(-2);
}

TEST(Wigner, SingleMaximumUpperIndexAccess) {
  CheckSingleUpperIndexAccess(5);
  CheckSingleUpperIndexAccess(-5);
}

// -- The transform-major layout.
//
// The claim to establish is that [n][m][l][theta] holds the same values
// as [n][theta][(l, m)], value for value and bit for bit. Bit-identity is the
// right standard rather than a tolerance: both layouts run the same recursion
// through WignerDetails::ComputeBlock with the same seeds and the same
// evaluation order, so any difference at all would be a layout fault rather
// than an arithmetic one -- and the whole purpose of the layout is that it
// changes nothing about the values.
//
// The test computes the index into the matrix itself instead of asking the
// class for it, so that it checks the documented layout rather than agreeing
// with an accessor that could be wrong the same way twice.
namespace {

template <typename Real, typename MRange, typename NRange>
void CheckTransformMajorAgreesWithWigner(std::ptrdiff_t lMax,
                                         std::ptrdiff_t mMax,
                                         std::ptrdiff_t nMax) {
  using namespace GSHTrans;
  using Int = std::ptrdiff_t;

  // Angles chosen away from the poles and unevenly, so that no accidental
  // symmetry of the sample points can hide an index that is transposed.
  auto angles = std::vector<Real>{};
  const auto nTheta = Int{7};
  for (auto i = Int{0}; i < nTheta; i++) {
    angles.push_back(static_cast<Real>(0.17) +
                     static_cast<Real>(i) * static_cast<Real>(0.39));
  }

  const auto table =
      Wigner<Real, MRange, NRange, Multiple>(lMax, mMax, nMax, angles);
  const auto matrices =
      WignerMatrices<Real, MRange, NRange>(lMax, mMax, nMax, angles);

  ASSERT_EQ(matrices.NumberOfAngles(), nTheta);
  ASSERT_EQ(matrices.MaxDegree(), lMax);
  ASSERT_EQ(matrices.MaxOrder(), mMax);

  // Every value the matrix layout holds is present in the block layout, and
  // equal. Walked from the matrix side, since that is the layout under test.
  auto matrixValues = std::size_t{0};
  for (auto n : matrices.UpperIndices()) {
    for (auto m : matrices.Orders()) {
      const auto block = matrices[n, m];
      const auto lMin = matrices.MinDegree(n, m);
      EXPECT_EQ(block.size(), static_cast<std::size_t>(
                                  matrices.NumberOfDegrees(n, m) * nTheta));

      for (auto l : matrices.Degrees(n, m)) {
        for (auto iTheta = Int{0}; iTheta < nTheta; iTheta++) {
          const auto fromMatrix = block[(l - lMin) * nTheta + iTheta];
          const auto fromTable = table[n, iTheta][l][m];
          EXPECT_EQ(fromMatrix, fromTable)
              << "n = " << n << ", m = " << m << ", l = " << l
              << ", iTheta = " << iTheta;
          matrixValues++;
        }
      }
    }
  }

  // And nothing is missing: the two layouts hold the same number of values.
  // This is the half the walk above cannot see, and it is the statement that
  // the transposed triangle really is the same triangle.
  auto tableValues = std::size_t{0};
  for (auto n : table.UpperIndices()) {
    for ([[maybe_unused]] auto iTheta : table.AngleIndices()) {
      tableValues +=
          static_cast<std::size_t>(GSHIndices<MRange>(lMax, mMax, n).Size());
    }
  }
  EXPECT_EQ(matrixValues, tableValues);
}

}  // namespace

TEST(WignerMatrices, AgreesWithTheBlockLayout) {
  using namespace GSHTrans;
  CheckTransformMajorAgreesWithWigner<double, All, All>(8, 8, 2);
}

TEST(WignerMatrices, AgreesWithTheBlockLayoutLongDouble) {
  using namespace GSHTrans;
  CheckTransformMajorAgreesWithWigner<long double, All, All>(6, 6, 2);
}

// mMax below lMax truncates the orders, which changes both layouts' shapes in
// different places -- the block loses columns at high degree, the matrix set
// loses whole matrices. That they still agree is the check that the identity
// is not an artefact of the square case.
TEST(WignerMatrices, AgreesWhenOrdersAreTruncated) {
  using namespace GSHTrans;
  CheckTransformMajorAgreesWithWigner<double, All, All>(9, 4, 3);
}

// The reduced m >= 0 storage of a real scalar grid, which is the one case
// where MinOrder() is zero rather than -mMax.
TEST(WignerMatrices, AgreesForNonNegativeOrders) {
  using namespace GSHTrans;
  CheckTransformMajorAgreesWithWigner<double, NonNegative, All>(7, 7, 0);
}

// The matrix at (n, m) starts at degree max(|n|, |m|) and its height falls
// linearly in |m|. That is the load imbalance a parallel kernel has to
// divide work for, so it is worth pinning as a property rather than leaving
// it implied by the agreement test.
TEST(WignerMatrices, MatrixHeightFallsWithOrder) {
  using namespace GSHTrans;
  using Int = std::ptrdiff_t;

  constexpr Int lMax = 10;
  const auto angles = std::vector<double>{0.3, 0.9, 1.7};
  const auto matrices = WignerMatrices<double, All, All>(lMax, lMax, 2, angles);

  for (auto n : matrices.UpperIndices()) {
    for (auto m : matrices.Orders()) {
      EXPECT_EQ(matrices.MinDegree(n, m), std::max(std::abs(n), std::abs(m)));
      EXPECT_EQ(matrices.NumberOfDegrees(n, m),
                lMax - std::max(std::abs(n), std::abs(m)) + 1);
    }
  }

  // The tallest matrix is at m = 0 and the shortest at |m| = mMax.
  EXPECT_EQ(matrices.NumberOfDegrees(0, 0), lMax + 1);
  EXPECT_EQ(matrices.NumberOfDegrees(0, lMax), 1);
  EXPECT_EQ(matrices.NumberOfDegrees(0, -lMax), 1);
}

// -- The reflected layout.
//
// D&T (C.118) in this library's stored values reads
//
//     d^l_{nm}(pi - theta) = (-1)^{l+n} d^l_{n,-m}(theta)
//
// so a table over colatitudes symmetric about pi/2 need store only the
// non-negative orders. These tests pin the relation itself, the halving it
// licenses, and the refusal when the angles do not support it -- because a
// table that is quietly the wrong values for half its orders is the failure
// mode here, and it would show up nowhere else until a transform was wrong.
namespace {

auto SymmetricAngles(std::ptrdiff_t nTheta) {
  // Symmetric about pi/2 by construction rather than by quadrature, so the
  // test does not depend on GaussQuad.
  auto theta = std::vector<double>{};
  for (auto i = std::ptrdiff_t{0}; i < nTheta; i++) {
    theta.push_back(std::numbers::pi_v<double> *
                    (static_cast<double>(i) + 0.5) /
                    static_cast<double>(nTheta));
  }
  return theta;
}

}  // namespace

TEST(WignerMatrices, ReflectionRecoversTheNegativeOrders) {
  using namespace GSHTrans;
  using Int = std::ptrdiff_t;

  constexpr Int lMax = 9;
  constexpr Int nMax = 2;
  const auto theta = SymmetricAngles(8);
  const auto nTheta = static_cast<Int>(theta.size());

  const auto full =
      WignerMatrices<double, All, All>::Full(lMax, lMax, nMax, theta);
  const auto half =
      WignerMatrices<double, All, All>::Reflected(lMax, lMax, nMax, theta);

  EXPECT_FALSE(full.IsReflected());
  EXPECT_TRUE(half.IsReflected());
  EXPECT_EQ(full.MinOrder(), -lMax);
  EXPECT_EQ(half.MinOrder(), 0);

  for (auto n : half.UpperIndices()) {
    for (auto m : half.Orders()) {
      ASSERT_GE(m, 0);
      const auto stored = half[n, m];
      const auto reference = full[n, m];
      ASSERT_EQ(stored.size(), reference.size());
      // What is kept is kept exactly.
      for (std::size_t i = 0; i < stored.size(); i++) {
        EXPECT_EQ(stored[i], reference[i]) << "n " << n << " m " << m;
      }

      if (m == 0) continue;

      // And what is dropped is recoverable: the matrix at -m is this one with
      // its columns reversed and the sign applied to row l.
      const auto lMin = half.MinDegree(n, m);
      ASSERT_EQ(full.MinDegree(n, -m), lMin);
      const auto negative = full[n, -m];
      for (auto l : half.Degrees(n, m)) {
        const auto sign = WignerMatrices<double, All, All>::Sign(l, n);
        for (auto iTheta = Int{0}; iTheta < nTheta; iTheta++) {
          const auto mirror = nTheta - 1 - iTheta;
          EXPECT_NEAR(negative[(l - lMin) * nTheta + iTheta],
                      sign * stored[(l - lMin) * nTheta + mirror], 1e-14)
              << "n " << n << " m " << m << " l " << l << " i " << iTheta;
        }
      }
    }
  }
}

TEST(WignerMatrices, ReflectedStorageIsHalfTheOrders) {
  using namespace GSHTrans;
  using Int = std::ptrdiff_t;

  constexpr Int lMax = 12;
  const auto theta = SymmetricAngles(10);

  auto count = [&](const auto& table) {
    std::size_t total = 0;
    for (auto n : table.UpperIndices()) {
      for (auto m : table.Orders()) total += table[n, m].size();
    }
    return total;
  };

  const auto full =
      WignerMatrices<double, All, All>::Full(lMax, lMax, 2, theta);
  const auto half =
      WignerMatrices<double, All, All>::Reflected(lMax, lMax, 2, theta);

  // Not exactly half: order zero is its own reflection and is stored once
  // either way. So the saving is (total - zeroth) / 2, and stating it that
  // way is the check that nothing else was dropped or duplicated.
  std::size_t zeroth = 0;
  for (auto n : full.UpperIndices()) zeroth += full[n, 0].size();
  EXPECT_EQ(count(half), (count(full) - zeroth) / 2 + zeroth);
}

TEST(WignerMatrices, ReflectedRefusesUnsymmetricAngles) {
  using namespace GSHTrans;
  // Symmetric about pi/2 is the whole premise; without it half the table
  // would be quietly wrong.
  const auto skewed = std::vector<double>{0.3, 0.9, 1.4, 2.0};
  EXPECT_THROW((WignerMatrices<double, All, All>::Reflected(6, 6, 1, skewed)),
               std::invalid_argument);
  // The unreflected layout takes any angles at all.
  EXPECT_NO_THROW((WignerMatrices<double, All, All>::Full(6, 6, 1, skewed)));
}
