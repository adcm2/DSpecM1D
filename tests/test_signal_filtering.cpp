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

TEST(SignalFilteringTests, FilterOptionsExposeExpectedDefaults) {
  DSpecM::FilterOptions options;
  EXPECT_DOUBLE_EQ(options.preTaper, 0.05);
  EXPECT_DOUBLE_EQ(options.finalTaper, 0.5);
  EXPECT_EQ(options.passes, 1);
  EXPECT_FALSE(options.enforceRealSignal);
}
