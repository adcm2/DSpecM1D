#include <gtest/gtest.h>
#include <algorithm>
#include <fstream>
#include <filesystem>
#include <iterator>
#include <string>
#include <vector>
#include <DSpecM1D/src/NormClass.h>
#include <DSpecM1D/src/OutputWriters.h>
#include <DSpecM1D/FrequencyTools>
#include "test_utils.h"

namespace {

SpectraSolver::FreqFull
makeFreq() {
  prem_norm<double> norm;
  return SpectraSolver::FreqFull(0.1, 1.0, 0.1, 0.2, 0.8, 1.0, 1.0, 1.0, 1.0,
                                 0.05, 0.0, 1.0, 1, norm.TimeNorm());
}

std::vector<std::string>
readLines(const std::filesystem::path &path) {
  std::ifstream in(path);
  std::vector<std::string> lines;
  for (std::string line; std::getline(in, line);) {
    lines.push_back(line);
  }
  return lines;
}

std::size_t
countChar(const std::string &text, char needle) {
  return static_cast<std::size_t>(std::count(text.begin(), text.end(), needle));
}

Eigen::MatrixXcd
makeComplexMatrix(int cols, double base) {
  Eigen::MatrixXcd out(3, cols);
  for (int row = 0; row < out.rows(); ++row) {
    for (int col = 0; col < out.cols(); ++col) {
      out(row, col) = std::complex<double>(base + row, 10.0 * base + col);
    }
  }
  return out;
}

Eigen::MatrixXd
makeRealMatrix(int cols, double base) {
  Eigen::MatrixXd out(3, cols);
  for (int row = 0; row < out.rows(); ++row) {
    for (int col = 0; col < out.cols(); ++col) {
      out(row, col) = base + row + 0.1 * col;
    }
  }
  return out;
}

InputParametersNew
makeParamsNew(const std::filesystem::path &dir) {
  const auto path = DSpecMTest::writeFile(
      dir / "params.txt",
      DSpecMTest::makeParameterText(DSpecMTest::modelPath().string()));
  return InputParametersNew(path.string());
}

}   // namespace

TEST(OutputWriterTests, FrequencyWriterTruncatesToSmallestInput) {
  DSpecMTest::TempDir temp;
  const auto path = temp.path() / "freq.out";
  auto freq = makeFreq();
  std::vector<double> angularFrequency{1.0, 2.0, 3.0};

  Eigen::MatrixXcd primary(3, 3);
  Eigen::MatrixXcd reference(3, 2);
  primary.setConstant(std::complex<double>(1.0, 0.0));
  reference.setConstant(std::complex<double>(2.0, 0.0));

  DSpecM::writeFrequencyComparison(path.string(), angularFrequency, freq, 1.0,
                                   primary, reference, 10, 6);

  const auto lines = readLines(path);
  ASSERT_EQ(lines.size(), 2u);
  EXPECT_NE(lines.front().find("159.154"), std::string::npos);
}

TEST(OutputWriterTests, ThreeWayFrequencyWriterIncludesAllReferenceBlocks) {
  DSpecMTest::TempDir temp;
  const auto path = temp.path() / "freq_three_way.out";
  auto freq = makeFreq();
  std::vector<double> angularFrequency{1.0, 2.0, 3.0};

  const auto primary = makeComplexMatrix(3, 1.0);
  const auto reference = makeComplexMatrix(2, 2.0);
  const auto secondary = makeComplexMatrix(1, 3.0);

  DSpecM::writeFrequencyComparison(path.string(), angularFrequency, freq, 1.0,
                                   primary, reference, secondary, 10, 6);

  const auto lines = readLines(path);
  ASSERT_EQ(lines.size(), 1u);
  EXPECT_EQ(countChar(lines.front(), ';'), 27u);
  EXPECT_NE(lines.front().find("3.000000"), std::string::npos);
  EXPECT_NE(lines.front().find("30.000000"), std::string::npos);
}

TEST(OutputWriterTests, FrequencyWriterThrowsWhenParentDirectoryIsMissing) {
  DSpecMTest::TempDir temp;
  const auto path = temp.path() / "missing" / "freq.out";
  auto freq = makeFreq();
  std::vector<double> angularFrequency{1.0};

  const auto primary = makeComplexMatrix(1, 1.0);
  const auto reference = makeComplexMatrix(1, 2.0);

  EXPECT_THROW(
      DSpecM::writeFrequencyComparison(path.string(), angularFrequency, freq,
                                       1.0, primary, reference),
      std::runtime_error);
}

TEST(OutputWriterTests, FrequencyWriterConvenienceOverloadUsesInputParametersNew) {
  DSpecMTest::TempDir temp;
  const auto path = temp.path() / "freq_params_new.out";
  auto paramsNew = makeParamsNew(temp.path());

  const auto cols =
      static_cast<int>(paramsNew.freqFull().w().size());
  const auto primary = makeComplexMatrix(cols, 1.0);
  const auto reference = makeComplexMatrix(cols, 2.0);

  DSpecM::writeFrequencyComparison(path.string(), paramsNew, primary, reference,
                                   -paramsNew.freqFull().i2() + 1, 6);

  const auto lines = readLines(path);
  ASSERT_EQ(lines.size(), 1u);
  EXPECT_EQ(countChar(lines.front(), ';'), 18u);
}

TEST(OutputWriterTests,
     ThreeWayFrequencyConvenienceOverloadUsesInputParametersNew) {
  DSpecMTest::TempDir temp;
  const auto path = temp.path() / "freq_params_new_three_way.out";
  auto paramsNew = makeParamsNew(temp.path());

  const auto cols =
      static_cast<int>(paramsNew.freqFull().w().size());
  const auto primary = makeComplexMatrix(cols, 1.0);
  const auto reference = makeComplexMatrix(cols, 2.0);
  const auto secondary = makeComplexMatrix(cols, 3.0);

  DSpecM::writeFrequencyComparison(path.string(), paramsNew, primary, reference,
                                   secondary, -paramsNew.freqFull().i2() + 1,
                                   6);

  const auto lines = readLines(path);
  ASSERT_EQ(lines.size(), 1u);
  EXPECT_EQ(countChar(lines.front(), ';'), 27u);
}

TEST(OutputWriterTests, TimeWriterStopsAtRequestedOutputWindow) {
  DSpecMTest::TempDir temp;
  const auto path = temp.path() / "time.out";
  auto freq = makeFreq();

  Eigen::MatrixXd primary = Eigen::MatrixXd::Zero(3, 4);
  Eigen::MatrixXd reference = Eigen::MatrixXd::Ones(3, 5);

  DSpecM::writeTimeComparison(path.string(), freq, 1.0, 0.02, primary,
                              reference, 3);

  const auto lines = readLines(path);
  ASSERT_FALSE(lines.empty());
  EXPECT_EQ(lines.size(), 4u);
  EXPECT_TRUE(lines.front().rfind("0.000;", 0) == 0);
}

TEST(OutputWriterTests, TimeWriterTruncatesToShortestInputAndKeepsLayout) {
  DSpecMTest::TempDir temp;
  const auto path = temp.path() / "time_truncate.out";
  auto freq = makeFreq();

  const auto primary = makeRealMatrix(3, 1.0);
  const auto reference = makeRealMatrix(2, 2.0);

  DSpecM::writeTimeComparison(path.string(), freq, 1.0, 10.0, primary,
                              reference, 4);

  const auto lines = readLines(path);
  ASSERT_EQ(lines.size(), 2u);
  EXPECT_EQ(countChar(lines.front(), ';'), 6u);
}

TEST(OutputWriterTests, TimeWriterThrowsWhenParentDirectoryIsMissing) {
  DSpecMTest::TempDir temp;
  const auto path = temp.path() / "missing" / "time.out";
  auto freq = makeFreq();

  const auto primary = makeRealMatrix(1, 1.0);
  const auto reference = makeRealMatrix(1, 2.0);

  EXPECT_THROW(DSpecM::writeTimeComparison(path.string(), freq, 1.0, 1.0,
                                           primary, reference),
               std::runtime_error);
}

TEST(OutputWriterTests, TimeWriterConvenienceOverloadUsesInputParametersNew) {
  DSpecMTest::TempDir temp;
  const auto path = temp.path() / "time_params_new.out";
  auto paramsNew = makeParamsNew(temp.path());

  const auto primary = makeRealMatrix(2, 1.0);
  const auto reference = makeRealMatrix(3, 2.0);

  DSpecM::writeTimeComparison(path.string(), paramsNew, primary, reference, 4);

  const auto lines = readLines(path);
  ASSERT_FALSE(lines.empty());
  EXPECT_EQ(countChar(lines.front(), ';'), 6u);
}
