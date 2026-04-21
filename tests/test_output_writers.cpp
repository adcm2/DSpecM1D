#include <gtest/gtest.h>
#include <fstream>
#include <string>
#include <vector>
#include <DSpecM1D/src/NormClass.h>
#include <DSpecM1D/src/OutputWriters.h>
#include <SpectraSolver/FF>
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
  EXPECT_LE(lines.size(), 3u);
}
