#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include <stdexcept>
#include <DSpecM1D/src/InputParser.h>
#include "test_utils.h"

TEST(InputParserTests, GetNextValueLineSkipsCommentsAndQuotes) {
  DSpecMTest::TempDir temp;
  auto path = DSpecMTest::writeFile(
      temp.path() / "quoted.txt",
      "# comment\n\n   \"models/prem.200.noatten.txt\"\n");

  std::ifstream in(path);
  ASSERT_TRUE(in.is_open());
  EXPECT_EQ(get_next_value_line(in), "models/prem.200.noatten.txt");
}

TEST(InputParserTests, ReadRequiredScalarRejectsTrailingTokens) {
  DSpecMTest::TempDir temp;
  auto path = DSpecMTest::writeFile(temp.path() / "bad_scalar.txt", "10 extra\n");

  std::ifstream in(path);
  ASSERT_TRUE(in.is_open());
  EXPECT_THROW(read_required_scalar<int>(in, "lmax"), std::runtime_error);
}

TEST(InputParserTests, ReadRequiredLatLonParsesPair) {
  DSpecMTest::TempDir temp;
  auto path =
      DSpecMTest::writeFile(temp.path() / "latlon.txt", "51.5 -0.12\n");

  std::ifstream in(path);
  ASSERT_TRUE(in.is_open());
  const auto [lat, lon] = read_required_lat_lon(in, "receiver_lat_lon");
  EXPECT_DOUBLE_EQ(lat, 51.5);
  EXPECT_DOUBLE_EQ(lon, -0.12);
}

TEST(InputParserTests, InputParametersParsesKnownGoodFile) {
  DSpecMTest::TempDir temp;
  auto path = DSpecMTest::writeFile(
      temp.path() / "params.txt",
      DSpecMTest::makeParameterText(
          DSpecMTest::modelPath().string(), 2, 1e-4));

  InputParameters params(path.string());
  EXPECT_EQ(params.output_type(), 2);
  EXPECT_EQ(params.num_receivers(), 2);
  EXPECT_EQ(params.lmax(), 6);
  EXPECT_NEAR(params.relative_error(), 1e-4, 1e-12);
}

TEST(InputParserTests, InputParametersRejectsOutOfRangeLatitude) {
  DSpecMTest::TempDir temp;
  std::string bad = DSpecMTest::makeParameterText(
      DSpecMTest::modelPath().string(), 0, 1e-5);
  const std::string needle = "45.0 90.0";
  bad.replace(bad.find(needle), needle.size(), "95.0 90.0");
  auto path = DSpecMTest::writeFile(temp.path() / "bad_params.txt", bad);

  EXPECT_THROW(InputParameters(path.string()), std::runtime_error);
}
