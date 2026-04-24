#include <gtest/gtest.h>
#include <DSpecM1D/src/NormClass.h>
#include <DSpecM1D/src/SignalFiltering.h>
#include <DSpecM1D/FrequencyTools>
#include "test_utils.h"

namespace {

SpectraSolver::FreqFull
makeFilterFreq() {
  prem_norm<double> norm;
  return SpectraSolver::FreqFull(0.1, 1.0, 0.1, 0.2, 0.8, 1.0, 1.0, 1.0, 1.0,
                                 0.05, 0.0, 1.0, 1, norm.TimeNorm());
}

Eigen::MatrixXcd
makeComplexFrequencyInput(const SpectraSolver::FreqFull& freq) {
  Eigen::MatrixXcd out(3, freq.nt() / 2 + 1);
  for (int row = 0; row < out.rows(); ++row) {
    for (int col = 0; col < out.cols(); ++col) {
      const double omega = static_cast<double>(col) * 0.02;
      out(row, col) =
          std::complex<double>((row + 1) * std::cos(omega),
                               0.2 * (row + 1) * std::sin(omega));
    }
  }
  return out;
}

Eigen::MatrixXd
makeRealTimeInput(const SpectraSolver::FreqFull& freq) {
  Eigen::MatrixXd out(3, freq.nt());
  for (int row = 0; row < out.rows(); ++row) {
    for (int col = 0; col < out.cols(); ++col) {
      const double time = static_cast<double>(col) * freq.dt();
      out(row, col) = (row + 1) * std::sin(2.0 * M_PI * 0.0005 * time) +
                      0.1 * std::cos(2.0 * M_PI * 0.001 * time);
    }
  }
  return out;
}

}   // namespace

TEST(SignalFilteringTests, RejectsInvalidPassCount) {
  auto freq = makeFilterFreq();
  Eigen::MatrixXcd raw = Eigen::MatrixXcd::Zero(3, freq.nt() / 2 + 1);

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

TEST(SignalFilteringTests, ZeroFrequencyInputStaysZeroAcrossMultiplePasses) {
  auto freq = makeFilterFreq();
  Eigen::MatrixXcd raw = Eigen::MatrixXcd::Zero(3, freq.nt() / 2 + 1);

  DSpecM::FilterOptions options;
  options.passes = 2;
  options.enforceRealSignal = true;

  const auto result = DSpecM::applyFilter(raw, freq, options);
  ASSERT_EQ(result.timeSeries.rows(), 3);
  ASSERT_EQ(result.timeSeries.cols(), freq.nt());
  ASSERT_EQ(result.frequencySeries.rows(), 3);
  ASSERT_EQ(result.frequencySeries.cols(), freq.nt0() / 2 + 1);
  EXPECT_TRUE(result.timeSeries.isZero(0.0));
  EXPECT_TRUE(result.frequencySeries.isZero(0.0));
}

TEST(SignalFilteringTests, RejectsFrequencyInputWithUnexpectedColumnCount) {
  auto freq = makeFilterFreq();
  Eigen::MatrixXcd raw = Eigen::MatrixXcd::Zero(3, freq.nt() / 2);

  EXPECT_THROW(DSpecM::applyFilter(raw, freq), std::invalid_argument);
}

TEST(SignalFilteringTests, RejectsTimeInputWithUnexpectedColumnCount) {
  auto freq = makeFilterFreq();
  Eigen::MatrixXd raw = Eigen::MatrixXd::Zero(3, freq.nt() - 1);

  EXPECT_THROW(DSpecM::applyFilter(raw, freq), std::invalid_argument);
}

TEST(SignalFilteringTests, TimeInputPreservesDimensionsAndProducesFiniteOutput) {
  auto freq = makeFilterFreq();
  const auto raw = makeRealTimeInput(freq);

  DSpecM::FilterOptions options;
  options.passes = 2;

  const auto result = DSpecM::applyFilter(raw, freq, options);
  ASSERT_EQ(result.timeSeries.rows(), raw.rows());
  ASSERT_EQ(result.timeSeries.cols(), freq.nt());
  ASSERT_EQ(result.frequencySeries.rows(), raw.rows());
  ASSERT_EQ(result.frequencySeries.cols(), freq.nt0() / 2 + 1);
  EXPECT_TRUE(result.timeSeries.array().isFinite().all());
  EXPECT_TRUE(result.frequencySeries.real().array().isFinite().all());
  EXPECT_TRUE(result.frequencySeries.imag().array().isFinite().all());
  EXPECT_GT(result.timeSeries.cwiseAbs().maxCoeff(), 0.0);
}

TEST(SignalFilteringTests,
     FrequencyInputPreservesDimensionsAndSupportsEnforceRealSignal) {
  auto freq = makeFilterFreq();
  const auto raw = makeComplexFrequencyInput(freq);

  DSpecM::FilterOptions options;
  options.enforceRealSignal = true;
  options.passes = 2;

  const auto result = DSpecM::applyFilter(raw, freq, options);
  ASSERT_EQ(result.timeSeries.rows(), raw.rows());
  ASSERT_EQ(result.timeSeries.cols(), freq.nt());
  ASSERT_EQ(result.frequencySeries.rows(), raw.rows());
  ASSERT_EQ(result.frequencySeries.cols(), freq.nt0() / 2 + 1);
  EXPECT_TRUE(result.timeSeries.array().isFinite().all());
  EXPECT_TRUE(result.frequencySeries.real().array().isFinite().all());
  EXPECT_TRUE(result.frequencySeries.imag().array().isFinite().all());
  EXPECT_GT(result.frequencySeries.cwiseAbs().maxCoeff(), 0.0);
}
