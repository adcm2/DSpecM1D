#include <gtest/gtest.h>
#include <DSpecM1D/src/ReferenceSeriesIO.h>
#include "test_utils.h"

TEST(ReferenceSeriesIOTests, LoadYSpecTimeSeriesPadsAndTruncates) {
  DSpecMTest::TempDir temp;
  const auto path = DSpecMTest::writeFile(
      temp.path() / "yspec.txt",
      "0 1 2 3\n"
      "1 4 5 6\n");

  const auto matrix = DSpecM::loadYSpecTimeSeries(path.string(), 4);
  ASSERT_EQ(matrix.rows(), 3);
  ASSERT_EQ(matrix.cols(), 4);
  EXPECT_DOUBLE_EQ(matrix(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(matrix(2, 1), 6.0);
  EXPECT_DOUBLE_EQ(matrix(1, 3), 0.0);
}

TEST(ReferenceSeriesIOTests, LoadSpecnmTimeSeriesParsesSemicolonFormat) {
  DSpecMTest::TempDir temp;
  const auto path = DSpecMTest::writeFile(
      temp.path() / "specnm.txt",
      "0;1;2;3\n"
      "1;4;5;6\n");

  const auto matrix = DSpecM::loadSpecnmTimeSeries(path.string(), 2);
  ASSERT_EQ(matrix.rows(), 3);
  ASSERT_EQ(matrix.cols(), 2);
  EXPECT_DOUBLE_EQ(matrix(0, 1), 4.0);
  EXPECT_DOUBLE_EQ(matrix(2, 0), 3.0);
}

TEST(ReferenceSeriesIOTests, LoadMineosTimeSeriesAppliesScale) {
  DSpecMTest::TempDir temp;
  const auto z = DSpecMTest::writeFile(temp.path() / "z.txt", "0 1\n1 2\n");
  const auto n = DSpecMTest::writeFile(temp.path() / "n.txt", "0 3\n1 4\n");
  const auto e = DSpecMTest::writeFile(temp.path() / "e.txt", "0 5\n1 6\n");

  const auto matrix =
      DSpecM::loadMineosTimeSeries(z.string(), n.string(), e.string(), 3, 2.0);
  ASSERT_EQ(matrix.rows(), 3);
  ASSERT_EQ(matrix.cols(), 3);
  EXPECT_DOUBLE_EQ(matrix(0, 0), 2.0);
  EXPECT_DOUBLE_EQ(matrix(1, 1), 8.0);
  EXPECT_DOUBLE_EQ(matrix(2, 2), 0.0);
}
