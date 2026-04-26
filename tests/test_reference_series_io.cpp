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

TEST(ReferenceSeriesIOTests, LoadSpecnmTimeSeriesPadsShortInputs) {
  DSpecMTest::TempDir temp;
  const auto path = DSpecMTest::writeFile(temp.path() / "specnm.txt",
                                          "0;1;2;3\n");

  const auto matrix = DSpecM::loadSpecnmTimeSeries(path.string(), 3);
  ASSERT_EQ(matrix.rows(), 3);
  ASSERT_EQ(matrix.cols(), 3);
  EXPECT_DOUBLE_EQ(matrix(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(matrix(1, 0), 2.0);
  EXPECT_DOUBLE_EQ(matrix(2, 0), 3.0);
  EXPECT_DOUBLE_EQ(matrix(0, 2), 0.0);
  EXPECT_DOUBLE_EQ(matrix(2, 2), 0.0);
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

TEST(ReferenceSeriesIOTests, LoadMineosTimeSeriesTruncatesToShortestComponent) {
  DSpecMTest::TempDir temp;
  const auto z = DSpecMTest::writeFile(
      temp.path() / "z.txt",
      "0 1\n"
      "1 2\n"
      "2 3\n");
  const auto n = DSpecMTest::writeFile(
      temp.path() / "n.txt",
      "0 4\n"
      "1 5\n");
  const auto e = DSpecMTest::writeFile(
      temp.path() / "e.txt",
      "0 6\n"
      "1 7\n"
      "2 8\n"
      "3 9\n");

  const auto matrix =
      DSpecM::loadMineosTimeSeries(z.string(), n.string(), e.string(), 5, 1.0);
  ASSERT_EQ(matrix.rows(), 3);
  ASSERT_EQ(matrix.cols(), 5);
  EXPECT_DOUBLE_EQ(matrix(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(matrix(1, 1), 5.0);
  EXPECT_DOUBLE_EQ(matrix(2, 1), 7.0);
  EXPECT_DOUBLE_EQ(matrix(0, 2), 0.0);
  EXPECT_DOUBLE_EQ(matrix(2, 4), 0.0);
}

TEST(ReferenceSeriesIOTests, LoadReferenceTimeSeriesBundlesYspecAndMineos) {
  DSpecMTest::TempDir temp;
  const auto yspec = DSpecMTest::writeFile(
      temp.path() / "yspec.txt",
      "0 1 2 3\n"
      "1 4 5 6\n");
  const auto z = DSpecMTest::writeFile(temp.path() / "z.txt", "0 7\n1 8\n");
  const auto n = DSpecMTest::writeFile(temp.path() / "n.txt", "0 9\n1 10\n");
  const auto e =
      DSpecMTest::writeFile(temp.path() / "e.txt", "0 11\n1 12\n");

  const auto bundle = DSpecM::loadReferenceTimeSeries(
      yspec.string(), z.string(), n.string(), e.string(), 3, 0.5);

  ASSERT_EQ(bundle.yspecTime.rows(), 3);
  ASSERT_EQ(bundle.yspecTime.cols(), 3);
  EXPECT_DOUBLE_EQ(bundle.yspecTime(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(bundle.yspecTime(2, 1), 6.0);
  EXPECT_DOUBLE_EQ(bundle.yspecTime(1, 2), 0.0);

  ASSERT_EQ(bundle.mineosTime.rows(), 3);
  ASSERT_EQ(bundle.mineosTime.cols(), 3);
  EXPECT_DOUBLE_EQ(bundle.mineosTime(0, 0), 3.5);
  EXPECT_DOUBLE_EQ(bundle.mineosTime(1, 1), 5.0);
  EXPECT_DOUBLE_EQ(bundle.mineosTime(2, 2), 0.0);
}

TEST(ReferenceSeriesIOTests, LoadersReturnEmptyMatricesForZeroRequestedColumns) {
  DSpecMTest::TempDir temp;
  const auto yspec = DSpecMTest::writeFile(temp.path() / "yspec.txt",
                                           "0 1 2 3\n");
  const auto specnm = DSpecMTest::writeFile(temp.path() / "specnm.txt",
                                            "0;4;5;6\n");
  const auto z = DSpecMTest::writeFile(temp.path() / "z.txt", "0 1\n");
  const auto n = DSpecMTest::writeFile(temp.path() / "n.txt", "0 2\n");
  const auto e = DSpecMTest::writeFile(temp.path() / "e.txt", "0 3\n");

  const auto yspecMatrix = DSpecM::loadYSpecTimeSeries(yspec.string(), 0);
  const auto specnmMatrix = DSpecM::loadSpecnmTimeSeries(specnm.string(), 0);
  const auto mineosMatrix =
      DSpecM::loadMineosTimeSeries(z.string(), n.string(), e.string(), 0);

  EXPECT_EQ(yspecMatrix.rows(), 3);
  EXPECT_EQ(yspecMatrix.cols(), 0);
  EXPECT_EQ(specnmMatrix.rows(), 3);
  EXPECT_EQ(specnmMatrix.cols(), 0);
  EXPECT_EQ(mineosMatrix.rows(), 3);
  EXPECT_EQ(mineosMatrix.cols(), 0);
}
