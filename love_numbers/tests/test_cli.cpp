#include <gtest/gtest.h>

#include <array>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#include <DSpecM1D/LoveNumbers>

namespace {

struct CommandResult {
  int status;
  std::string standard_output;
  std::string standard_error;
};

struct OutputRow {
  int degree;
  std::array<double, 8> values;
};

class TemporaryDirectory {
public:
  TemporaryDirectory() {
    std::mt19937_64 random(
        std::random_device{}() ^
        std::chrono::steady_clock::now().time_since_epoch().count());
    const std::filesystem::path base =
        std::filesystem::temp_directory_path();
    for (int attempt = 0; attempt < 32; ++attempt) {
      path_ = base / ("dspecm1d-love-cli-" +
                      std::to_string(random()));
      std::error_code error;
      if (std::filesystem::create_directory(path_, error)) {
        return;
      }
    }
    throw std::runtime_error("Could not create a temporary directory.");
  }

  ~TemporaryDirectory() {
    std::error_code error;
    std::filesystem::remove_all(path_, error);
  }

  const std::filesystem::path &path() const { return path_; }

private:
  std::filesystem::path path_;
};

std::string readText(const std::filesystem::path &path) {
  std::ifstream input(path);
  return std::string(
      std::istreambuf_iterator<char>(input),
      std::istreambuf_iterator<char>());
}

std::string shellQuote(const std::string &value) {
  std::string quoted = "'";
  for (const char character : value) {
    quoted += character == '\'' ? "'\\''" : std::string(1, character);
  }
  return quoted + "'";
}

CommandResult runCli(const std::vector<std::string> &arguments,
                     const TemporaryDirectory &temporary_directory) {
  const std::filesystem::path standard_output =
      temporary_directory.path() / "stdout.txt";
  const std::filesystem::path standard_error =
      temporary_directory.path() / "stderr.txt";

  std::ostringstream command;
  command << shellQuote(DSPECM1D_LOVE_NUMBERS_CLI);
  for (const std::string &argument : arguments) {
    command << ' ' << shellQuote(argument);
  }
  command << " > " << shellQuote(standard_output.string())
          << " 2> " << shellQuote(standard_error.string());

  const int status = std::system(command.str().c_str());
  return {
      .status = status,
      .standard_output = readText(standard_output),
      .standard_error = readText(standard_error),
  };
}

std::vector<OutputRow> parseRows(const std::string &text) {
  std::vector<OutputRow> rows;
  std::istringstream input(text);
  for (std::string line; std::getline(input, line);) {
    if (line.empty() || line.front() == '#') {
      continue;
    }

    std::istringstream row_input(line);
    OutputRow row;
    if (!(row_input >> row.degree)) {
      throw std::runtime_error("Could not read output degree.");
    }
    for (double &value : row.values) {
      if (!(row_input >> value)) {
        throw std::runtime_error("Output row has fewer than nine columns.");
      }
    }
    std::string extra;
    if (row_input >> extra) {
      throw std::runtime_error("Output row has more than nine columns.");
    }
    rows.push_back(row);
  }
  return rows;
}

TEST(LoveNumbersCliTests, ExplicitOptionsWriteFileAndMatchPublicResults) {
  TemporaryDirectory temporary_directory;
  const std::filesystem::path output =
      temporary_directory.path() / "love.txt";
  const CommandResult command = runCli(
      {DSPECM1D_LOVE_NUMBERS_CLI_MODEL, "2", output.string(),
       "2", "0.4"},
      temporary_directory);

  ASSERT_EQ(command.status, 0) << command.standard_error;
  EXPECT_TRUE(command.standard_output.empty());
  const std::string text = readText(output);
  EXPECT_NE(text.find("# polynomial_order: 2"), std::string::npos);
  EXPECT_NE(
      text.find("# maximum_radial_step: 4.0000000000000002e-01"),
      std::string::npos);

  const std::vector<OutputRow> rows = parseRows(text);
  ASSERT_EQ(rows.size(), 3);
  const EarthModels::ModelInput<double> model(
      DSPECM1D_LOVE_NUMBERS_CLI_MODEL);
  const std::vector<DSpecM1D::LoveNumbers::DegreeResult> expected =
      DSpecM1D::LoveNumbers::calculate(
          model,
          {.maximum_degree = 2,
           .polynomial_order = 2,
           .maximum_radial_step = 0.4});

  for (int degree = 0; degree <= 2; ++degree) {
    EXPECT_EQ(rows[degree].degree, degree);
    const auto &result = expected[degree];
    const std::array<double, 8> values{
        result.h_u,      result.k_u,      result.h_phi, result.k_phi,
        result.h_load(), result.k_load(), result.h_t,   result.k_t,
    };
    for (std::size_t column = 0; column < values.size(); ++column) {
      EXPECT_DOUBLE_EQ(rows[degree].values[column], values[column]);
    }
  }
}

TEST(LoveNumbersCliTests, DefaultsWriteStandardOutput) {
  TemporaryDirectory temporary_directory;
  const CommandResult command = runCli(
      {DSPECM1D_LOVE_NUMBERS_CLI_MODEL, "0", "-"},
      temporary_directory);

  ASSERT_EQ(command.status, 0) << command.standard_error;
  EXPECT_NE(
      command.standard_output.find("# polynomial_order: 6"),
      std::string::npos);
  EXPECT_NE(
      command.standard_output.find(
          "# maximum_radial_step: 1.0000000000000000e-02"),
      std::string::npos);
  const std::vector<OutputRow> rows =
      parseRows(command.standard_output);
  ASSERT_EQ(rows.size(), 1);
  EXPECT_EQ(rows.front().degree, 0);
}

TEST(LoveNumbersCliTests, InvalidArgumentsAndMissingModelFail) {
  TemporaryDirectory temporary_directory;
  const CommandResult usage =
      runCli({}, temporary_directory);
  EXPECT_NE(usage.status, 0);
  EXPECT_NE(usage.standard_error.find("usage: dspecm1d-love"),
            std::string::npos);

  const CommandResult invalid = runCli(
      {DSPECM1D_LOVE_NUMBERS_CLI_MODEL, "not-a-degree", "-"},
      temporary_directory);
  EXPECT_NE(invalid.status, 0);
  EXPECT_NE(invalid.standard_error.find("LMAX must be an integer"),
            std::string::npos);

  const CommandResult missing = runCli(
      {(temporary_directory.path() / "missing-model.txt").string(),
       "0", "-", "2", "0.4"},
      temporary_directory);
  EXPECT_NE(missing.status, 0);
  EXPECT_NE(missing.standard_error.find("dspecm1d-love:"),
            std::string::npos);
}

TEST(LoveNumbersCliTests, SurfaceFluidIsRejectedClearly) {
  TemporaryDirectory temporary_directory;
  const CommandResult command = runCli(
      {DSPECM1D_LOVE_NUMBERS_CLI_OCEAN_MODEL, "0", "-", "2", "0.4"},
      temporary_directory);

  EXPECT_NE(command.status, 0);
  EXPECT_NE(
      command.standard_error.find(
          "Surface fluids and oceans are not supported"),
      std::string::npos);
}

}   // namespace
