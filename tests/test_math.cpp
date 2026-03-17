// tests/test_math.cpp
#include <gtest/gtest.h>
#include <Eigen/Dense>

// #include "YourRealHeader.h" // You would include your library headers here

TEST(MathTests, BasicEigenAddition) {
  Eigen::Matrix2d a;
  a << 1.0, 2.0, 3.0, 4.0;

  Eigen::Matrix2d b;
  b << 1.0, 1.0, 1.0, 1.0;

  Eigen::Matrix2d result = a + b;

  // 1. Exact comparison (Good for integers or exact math)
  EXPECT_DOUBLE_EQ(result(0, 0), 2.0);

  // 2. Near comparison (CRITICAL for floating point/FFT math)
  // Checks if the difference is within 1e-9
  EXPECT_NEAR(result(1, 1), 5.0, 1e-9);

  // 3. Testing whole Eigen matrices at once using Eigen's built-in fuzzy
  // compare
  Eigen::Matrix2d expected;
  expected << 2.0, 3.0, 4.0, 5.0;

  EXPECT_TRUE(result.isApprox(expected, 1e-9));
}