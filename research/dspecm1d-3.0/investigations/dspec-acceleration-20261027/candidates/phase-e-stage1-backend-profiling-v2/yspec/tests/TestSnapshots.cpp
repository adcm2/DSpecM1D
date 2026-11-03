#include <gtest/gtest.h>

#include <YSpec/Location.hpp>
#include <YSpec/Snapshots.hpp>
#include <YSpec/Synthetics.hpp>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

#include "Geometry.hpp"

namespace {

using namespace YSpec;
constexpr auto pi = std::numbers::pi;

const auto kPrem = std::string(YSPEC_EXAMPLES_DIR) + "/prem.200";

Calculation Small() {
  auto c = Calculation{};
  c.degrees = {0, 40};
  c.band = {0.2e-3, 5e-3};
  c.sampling = {3600, 60};
  c.taper = Taper{0.1e-3, 0.3e-3, 4e-3, 5e-3};
  return c;
}

TEST(Snapshots, HoldTheSamplesOfSeismogramsAtTheSamePoints) {
  const auto model = Model::ReadDeck(kPrem);
  auto calculation = Small();
  calculation.quantity = Quantity::Velocity;
  const auto sources = std::vector<Source>{
      {{10, 20}, 30e3, {.rr = 1e18, .rt = 2e18, .tp = -1e18}},
      {{-30, 100}, 200e3, {.rp = 1e18, .tt = 1e18, .pp = -1e18}}};
  // On the sea surface, beneath the sea floor, in the mantle and in the
  // outer core.
  const auto grid = SnapshotGrid{.depths = {0, 3e3, 400e3, 4000e3},
                                 .distances = {0, 30, 75, 150},
                                 .azimuths = {40, 220}};
  const auto samples = std::vector<std::size_t>{10, 20, 30};
  auto times = std::vector<double>{};
  for (auto n : samples) times.push_back(60.0 * static_cast<double>(n));
  const auto snapshots =
      ComputeSnapshots(model, calculation, sources, grid, times);

  // A receiver at each point of the grid around each source. Computed
  // together, the seismograms share the snapshots' radial solves exactly.
  auto receivers = std::vector<Receiver>{};
  for (const auto& source : sources) {
    for (auto depth : grid.depths) {
      for (auto azimuth : grid.azimuths) {
        for (auto distance : grid.distances) {
          const auto place =
              GreatCircle(source.location, azimuth, distance, distance, 1);
          receivers.push_back({place[0], depth});
        }
      }
    }
  }
  const auto seismograms =
      ComputeSeismograms(model, calculation, sources, receivers);
  const auto points = receivers.size() / 2;

  for (auto s : {0, 1}) {
    const auto& source = sources[s];
    auto peak = 0.0;
    for (auto k = s * points; k < (s + 1) * points; ++k) {
      for (auto c : {Vertical, North, East}) {
        for (auto x : seismograms.Trace(s, k, c)) {
          peak = std::max(peak, std::abs(x));
        }
      }
    }
    ASSERT_GT(peak, 0);

    auto k = s * points;
    for (auto d = std::size_t{0}; d < grid.depths.size(); ++d) {
      for (auto a : {0, 1}) {
        for (auto i : {0, 1, 2, 3}) {
          // The azimuth of the direction away from the source at the point.
          const auto beta =
              grid.distances[i] == 0
                  ? grid.azimuths[a] * pi / 180
                  : pi - Geometry(receivers[k].location, source.location)
                             .rotation;
          for (auto t : {0, 1, 2}) {
            const auto n = samples[t];
            const auto away = snapshots(s, t, d, i, a, Away);
            const auto transverse = snapshots(s, t, d, i, a, Transverse);
            const auto where = ::testing::Message()
                               << "source " << s << ", depth " << d
                               << ", azimuth " << a << ", distance " << i
                               << ", time " << t;
            EXPECT_NEAR(snapshots(s, t, d, i, a, Up),
                        seismograms(s, k, Vertical, n), 1e-9 * peak)
                << where;
            EXPECT_NEAR(away * std::cos(beta) + transverse * std::sin(beta),
                        seismograms(s, k, North, n), 1e-9 * peak)
                << where;
            EXPECT_NEAR(away * std::sin(beta) - transverse * std::cos(beta),
                        seismograms(s, k, East, n), 1e-9 * peak)
                << where;
          }
          ++k;
        }
      }
    }
  }
}

TEST(Snapshots, RefuseWhatTheyCannotCompute) {
  const auto model = Model::ReadDeck(kPrem);
  const auto calculation = Small();
  const auto sources = std::vector<Source>{{{0, 0}, 30e3, {.rr = 1e18}}};
  const auto grid =
      SnapshotGrid{.depths = {100e3}, .distances = {10, 20}, .azimuths = {0}};
  const auto times = std::vector<double>{600};
  auto refuses = [&](const SnapshotGrid& g, const std::vector<double>& t) {
    EXPECT_THROW((void)ComputeSnapshots(model, calculation, sources, g, t),
                 std::invalid_argument);
  };
  refuses(grid, {4000});
  refuses(grid, {});
  auto deep = grid;
  deep.depths = {7000e3};
  refuses(deep, times);
  auto far = grid;
  far.distances.assign(1, 190.0);
  refuses(far, times);
  auto empty = grid;
  empty.azimuths.clear();
  refuses(empty, times);
}

}  // namespace
