#include <gtest/gtest.h>

#include <YSpec/Instrument.hpp>
#include <YSpec/Synthetics.hpp>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using namespace YSpec;
namespace fs = std::filesystem;

const auto kPrem = std::string(YSPEC_EXAMPLES_DIR) + "/prem.200";

fs::path Write(const std::string& name, const std::string& text) {
  const auto directory = fs::temp_directory_path() / "yspec-instrument-tests";
  fs::create_directories(directory);
  const auto path = directory / name;
  std::ofstream(path) << text;
  return path;
}

TEST(Instrument, ReadsASacPoleZeroFile) {
  // In the manner of the IRIS DMC's files, with a zero left unlisted.
  const auto instrument = ReadSacPoleZero(Write("sts.pz", R"(* *******
* NETWORK   (KNETWK): IU
* INPUT UNIT       : M/S
* OUTPUT UNIT      : COUNTS
ZEROS 2
  +1.0e+00  -2.0e+00
POLES 2
  -3.7e-02  +3.7e-02
  -3.7e-02  -3.7e-02
CONSTANT  +4.5e+09
)"));
  ASSERT_EQ(instrument.zeros.size(), 2u);
  EXPECT_EQ(instrument.zeros[0], Complex(1, -2));
  EXPECT_EQ(instrument.zeros[1], Complex(0, 0));
  ASSERT_EQ(instrument.poles.size(), 2u);
  EXPECT_EQ(instrument.poles[1], Complex(-3.7e-2, -3.7e-2));
  EXPECT_EQ(instrument.constant, 4.5e9);
  EXPECT_EQ(instrument.input, Quantity::Velocity);
}

TEST(Instrument, TakesDisplacementWhenNoUnitIsGiven) {
  const auto instrument =
      ReadSacPoleZero(Write("bare.pz", "zeros 1\npoles 0\nconstant 2\n"));
  EXPECT_EQ(instrument.input, Quantity::Displacement);
  EXPECT_EQ(instrument.zeros.size(), 1u);
  EXPECT_TRUE(instrument.poles.empty());
}

TEST(Instrument, RefusesFaultyFiles) {
  for (const auto& text :
       {"", "* nothing but a comment\n", "ZEROS 1\n1 0\n2 0\n",
        "POLES 1\n-1 0\nPOLES 1\n-2 0\n", "ZEROS 1\nhello\n", "1 0\n",
        "* INPUT UNIT : NM\nCONSTANT 1\n", "CONSTANT\n", "ZEROS -1\n"}) {
    EXPECT_THROW((void)ReadSacPoleZero(Write("faulty.pz", text)),
                 std::runtime_error)
        << text;
  }
  EXPECT_THROW((void)ReadSacPoleZero("no/such/file.pz"), std::runtime_error);
}

TEST(Instrument, RespondsAsItsPolesAndZerosSay) {
  // A first-order low pass, a / (s + a).
  const auto a = 2.0;
  const auto lowPass = Instrument{.zeros = {}, .poles = {-a}, .constant = a};
  EXPECT_NEAR(std::abs(lowPass.Response(0) - 1.0), 0, 1e-15);
  EXPECT_NEAR(std::abs(lowPass.Response(a)), 1 / std::sqrt(2.0), 1e-15);
  EXPECT_NEAR(std::arg(lowPass.Response(a)), -std::numbers::pi / 4, 1e-15);
}

// Seismograms for one source and two receivers, with an instrument at both.
Seismograms Record(Quantity quantity,
                   const std::optional<Instrument>& instrument) {
  const auto model = Model::ReadDeck(kPrem);
  auto calculation = Calculation{};
  calculation.degrees = {0, 30};
  calculation.band = {0.2e-3, 5e-3};
  calculation.sampling = {3600, 60};
  calculation.quantity = quantity;
  const auto source =
      Source{{10, 20}, 30e3, {.rr = 1e18, .rt = 2e18, .tp = -1e18}};
  auto receivers = std::vector<Receiver>{{{-20, 60}, 3e3}, {{40, 100}, 3e3}};
  for (auto& receiver : receivers) receiver.instrument = instrument;
  return ComputeSeismograms(model, calculation, std::span{&source, 1},
                            receivers);
}

void ExpectSame(const Seismograms& a, const Seismograms& b) {
  for (auto k : {0, 1}) {
    for (auto c : {Vertical, North, East}) {
      const auto trace = b.Trace(0, k, c);
      const auto peak = std::ranges::max(
          trace | std::views::transform([](double x) { return std::abs(x); }));
      ASSERT_GT(peak, 0);
      for (auto i = std::size_t{0}; i < b.Samples(); ++i) {
        EXPECT_NEAR(a(0, k, c, i), b(0, k, c, i), 1e-10 * peak)
            << "receiver " << k << ", component " << c << ", sample " << i;
      }
    }
  }
}

TEST(Instrument, ThatDifferentiatesRecordsVelocity) {
  // H(s) = s, responding to displacement.
  const auto differentiator =
      Instrument{.zeros = {0.0}, .poles = {}, .constant = 1};
  ExpectSame(Record(Quantity::Displacement, differentiator),
             Record(Quantity::Velocity, std::nullopt));
}

TEST(Instrument, RespondsToItsOwnInputWhateverIsComputed) {
  // H(s) = 1, responding to displacement, in a calculation of velocity or
  // acceleration, records displacement.
  const auto unit = Instrument{};
  const auto displacement = Record(Quantity::Displacement, std::nullopt);
  ExpectSame(Record(Quantity::Velocity, unit), displacement);
  ExpectSame(Record(Quantity::Acceleration, unit), displacement);
  // And one responding to velocity records velocity.
  const auto velocity = Instrument{
      .zeros = {}, .poles = {}, .constant = 1, .input = Quantity::Velocity};
  ExpectSame(Record(Quantity::Displacement, velocity),
             Record(Quantity::Velocity, std::nullopt));
}

}  // namespace
