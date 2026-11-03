#include <gtest/gtest.h>

#include <YSpec/Synthetics.hpp>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace {

using namespace YSpec;

const auto kPrem = std::string(YSPEC_EXAMPLES_DIR) + "/prem.200";

Calculation Small() {
  auto c = Calculation{};
  c.degrees = {0, 40};
  c.band = {0.2e-3, 5e-3};
  c.sampling = {3600, 60};
  return c;
}

const auto kMoment = MomentTensor{
    .rr = 1e18, .rt = 1e19, .rp = 5e18, .tt = -2e18, .tp = 3e18, .pp = 1e18};

double Peak(const Seismograms& s, std::size_t source, std::size_t receiver,
            Component c) {
  const auto trace = s.Trace(source, receiver, c);
  return std::ranges::max(
      trace | std::views::transform([](double x) { return std::abs(x); }));
}

TEST(Synthetics, AreContinuousAtTheEpicentre) {
  const auto model = Model::ReadDeck(kPrem);
  const auto source = Source{{20, 30}, 30e3, kMoment};
  const auto d = 0.001;
  auto receivers = std::vector<Receiver>{};
  for (auto [lat, lon] : {std::pair{20.0, 30.0},
                          {20 + d, 30.0},
                          {20 - d, 30.0},
                          {20.0, 30 + d},
                          {20.0, 30 - d}}) {
    receivers.push_back({{lat, lon}, 3e3});
  }
  const auto s =
      ComputeSeismograms(model, Small(), std::span{&source, 1}, receivers);
  for (auto c : {Vertical, North, East}) {
    const auto peak = Peak(s, 0, 0, c);
    ASSERT_GT(peak, 0);
    for (auto k : {1, 2, 3, 4}) {
      for (auto i = std::size_t{0}; i < s.Samples(); ++i) {
        EXPECT_NEAR(s(0, k, c, i), s(0, 0, c, i), 1e-2 * peak)
            << "component " << c << ", receiver " << k << ", sample " << i;
      }
    }
  }
}

TEST(Synthetics, ScaleWithTheMomentTensor) {
  const auto model = Model::ReadDeck(kPrem);
  const auto receivers =
      std::vector<Receiver>{{{-10, 60}, 3e3}, {{50, 100}, 3e3}};
  auto sources = std::vector<Source>{{{20, 30}, 30e3, kMoment}};
  const auto a = ComputeSeismograms(model, Small(), sources, receivers);
  sources[0].moment = {3 * kMoment.rr, 3 * kMoment.rt, 3 * kMoment.rp,
                       3 * kMoment.tt, 3 * kMoment.tp, 3 * kMoment.pp};
  const auto b = ComputeSeismograms(model, Small(), sources, receivers);
  for (auto k : {0, 1}) {
    for (auto c : {Vertical, North, East}) {
      for (auto i = std::size_t{0}; i < a.Samples(); ++i) {
        EXPECT_NEAR(b(0, k, c, i), 3 * a(0, k, c, i),
                    1e-12 * std::abs(b(0, k, c, i)) + 1e-30);
      }
    }
  }
}

TEST(Synthetics, GiveEachSourceWhatItGivesAlone) {
  // Several sources, at two depths and in different places, computed
  // together and one at a time. Together the integrations are shared and
  // stop at every depth, so the steps differ: they agree to the tolerance.
  const auto model = Model::ReadDeck(kPrem);
  auto calculation = Small();
  calculation.solver.tolerance = 1e-9;
  const auto receivers =
      std::vector<Receiver>{{{-10, 60}, 3e3}, {{50, 100}, 3e3}, {{5, 5}, 0}};
  const auto sources = std::vector<Source>{
      {{20, 30}, 30e3, kMoment},
      {{-40, 120}, 30e3, MomentTensor{.rr = 2e18, .tp = 1e18}},
      {{0, 0}, 200e3, MomentTensor{.rt = 1e18, .pp = -1e18}}};
  // The receiver at the surface is in the ocean.
  const auto together =
      ComputeSeismograms(model, calculation, sources, receivers);
  for (auto s : {0, 1, 2}) {
    const auto alone = ComputeSeismograms(model, calculation,
                                          std::span{&sources[s], 1}, receivers);
    for (auto k : {0, 1, 2}) {
      for (auto c : {Vertical, North, East}) {
        const auto peak = Peak(alone, 0, k, c);
        for (auto i = std::size_t{0}; i < alone.Samples(); ++i) {
          EXPECT_NEAR(together(s, k, c, i), alone(0, k, c, i), 1e-6 * peak);
        }
      }
    }
  }
}

TEST(Synthetics, RefuseASourceInAFluid) {
  const auto model = Model::ReadDeck(kPrem);
  const auto receivers = std::vector<Receiver>{{{10, 10}, 3e3}};
  for (auto depth : {1e3, 3000e3}) {
    const auto source = Source{{0, 0}, depth, kMoment};
    EXPECT_THROW((void)ComputeSeismograms(model, Small(), std::span{&source, 1},
                                          receivers),
                 std::invalid_argument)
        << depth;
  }
}

TEST(Synthetics, ServeReceiversAtSeveralDepthsTogether) {
  // As for a snapshot: receivers beneath the same places at several depths,
  // computed together and one depth at a time.
  const auto model = Model::ReadDeck(kPrem);
  auto calculation = Small();
  calculation.solver.tolerance = 1e-9;
  const auto sources = std::vector<Source>{{{20, 30}, 30e3, kMoment}};
  auto receivers = std::vector<Receiver>{};
  for (auto depth : {3e3, 100e3, 400e3, 1000e3}) {
    for (auto place : {Location{20, 50}, Location{-30, 100}}) {
      receivers.push_back({place, depth});
    }
  }
  const auto together =
      ComputeSeismograms(model, calculation, sources, receivers);
  // One place per depth suffices for the comparison: the places share the
  // depth's integration, so the second adds nothing but a full solve.
  for (auto k = std::size_t{0}; k < receivers.size(); k += 2) {
    const auto alone = ComputeSeismograms(model, calculation, sources,
                                          std::span{&receivers[k], 1});
    for (auto c : {Vertical, North, East}) {
      const auto peak = Peak(alone, 0, 0, c);
      for (auto i = std::size_t{0}; i < alone.Samples(); ++i) {
        EXPECT_NEAR(together(0, k, c, i), alone(0, 0, c, i), 1e-6 * peak)
            << "receiver " << k << " component " << c;
      }
    }
  }
}

}  // namespace
