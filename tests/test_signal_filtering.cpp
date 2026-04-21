#include <gtest/gtest.h>
#include <DSpecM1D/src/NormClass.h>
#include <DSpecM1D/src/SignalFiltering.h>
#include <SpectraSolver/FF>
#include "test_utils.h"

namespace {

SpectraSolver::FreqFull
makeFilterFreq() {
  prem_norm<double> norm;
  return SpectraSolver::FreqFull(0.1, 1.0, 0.1, 0.2, 0.8, 1.0, 1.0, 1.0, 1.0,
                                 0.05, 0.0, 1.0, 1, norm.TimeNorm());
}

}   // namespace

TEST(SignalFilteringTests, RejectsInvalidPassCount) {
  auto freq = makeFilterFreq();
  Eigen::MatrixXcd raw = Eigen::MatrixXcd::Zero(3, 8);

  DSpecM::FilterOptions options;
  options.passes = 0;

  EXPECT_THROW(DSpecM::applyFilter(raw, freq, options), std::invalid_argument);
}

TEST(SignalFilteringTests, ZeroFrequencyInputProducesZeroOutputs) {
  auto freq = makeFilterFreq();
  Eigen::MatrixXcd raw = Eigen::MatrixXcd::Zero(3, 8);

  const auto filtered = DSpecM::applyFilter(raw, freq);
  EXPECT_EQ(filtered.timeSeries.rows(), 3);
  EXPECT_EQ(filtered.frequencySeries.rows(), 3);
  EXPECT_TRUE(filtered.timeSeries.isZero(1e-12));
  EXPECT_TRUE(filtered.frequencySeries.isZero(1e-12));
}
