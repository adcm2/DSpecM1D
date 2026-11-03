#pragma once

/**
 * @file Snapshots.hpp
 * @brief The wavefield at chosen times, on a grid around each source.
 *
 * @details The grid is in the frame of each source: points at a set of
 * depths, distances from the source and azimuths at it. In that frame only
 * orders @f$|m| \le 2@f$ are excited, so at each depth and distance the field
 * is a sum over five orders of functions of the distance, each times
 * @f$e^{im\phi}@f$. The cost therefore grows with the numbers of depths and
 * distances, and the azimuths come almost free. A cross-section through a
 * source is a grid with two azimuths, @f$a@f$ and @f$a + 180@f$; a map is a
 * grid with one depth and many azimuths.
 *
 * The values at the chosen times are summed directly from the spectra, with
 * the damping undone, as the seismograms are by a fast Fourier transform: a
 * snapshot at the time of a sample of a seismogram at the same point holds
 * that sample.
 *
 * The points may lie anywhere but at the centre: in the core, or in an
 * ocean, as well as in the solid region of the sources.
 */

#include <cstddef>
#include <span>
#include <utility>
#include <vector>

#include "Calculation.hpp"
#include "Model.hpp"
#include "Source.hpp"

namespace YSpec {

/** @brief The points of a snapshot, in the frame of a source. */
struct SnapshotGrid {
  std::vector<double> depths;     ///< The depths (m).
  std::vector<double> distances;  ///< The distances from the source (degrees).
  std::vector<double> azimuths;   ///< At the source, clockwise from north.
};

/**
 * @brief The components of a snapshot: up, away from the source along the
 * great circle from it, and transverse to that, @f$\hat{up} \times
 * \hat{away}@f$, which is 90 degrees anticlockwise from away seen from above.
 */
enum SnapshotComponent { Up = 0, Away = 1, Transverse = 2 };

/**
 * @brief Snapshots for a set of sources, at a set of times, on a grid, in SI
 * units of the quantity computed.
 */
class Snapshots {
 public:
  /** @brief Zero snapshots of a given size. */
  Snapshots(std::size_t sources, std::vector<double> times, SnapshotGrid grid)
      : sources_{sources},
        times_{std::move(times)},
        grid_{std::move(grid)},
        data_(sources_ * times_.size() * grid_.depths.size() *
              grid_.distances.size() * grid_.azimuths.size() * 3) {}

  /** @brief The number of sources. */
  std::size_t Sources() const { return sources_; }
  /** @brief The times after the origin (s). */
  const std::vector<double>& Times() const { return times_; }
  /** @brief The grid. */
  const SnapshotGrid& Grid() const { return grid_; }

  /** @brief A value: indices of the source, time, depth, distance, azimuth. */
  double& operator()(std::size_t source, std::size_t time, std::size_t depth,
                     std::size_t distance, std::size_t azimuth,
                     SnapshotComponent component) {
    return data_[Offset(source, time, depth, distance, azimuth, component)];
  }
  /** @brief A value: indices of the source, time, depth, distance, azimuth. */
  double operator()(std::size_t source, std::size_t time, std::size_t depth,
                    std::size_t distance, std::size_t azimuth,
                    SnapshotComponent component) const {
    return data_[Offset(source, time, depth, distance, azimuth, component)];
  }

 private:
  std::size_t sources_;
  std::vector<double> times_;
  SnapshotGrid grid_;
  std::vector<double> data_;

  std::size_t Offset(std::size_t source, std::size_t time, std::size_t depth,
                     std::size_t distance, std::size_t azimuth,
                     SnapshotComponent component) const {
    auto offset = source;
    offset = offset * times_.size() + time;
    offset = offset * grid_.depths.size() + depth;
    offset = offset * grid_.azimuths.size() + azimuth;
    offset = offset * grid_.distances.size() + distance;
    return offset * 3 + component;
  }
};

/**
 * @brief Computes snapshots of the wavefield.
 * @param model The earth model.
 * @param calculation What to compute. The sampling's duration sets the
 * damping and the spacing of the frequencies, as for seismograms; its time
 * step only needs to resolve the band.
 * @param sources The sources.
 * @param grid The points, the same in the frame of each source.
 * @param times The times after the origin (s), between zero and the duration.
 * @param progress Called as frequencies are completed, if given.
 * @throws std::invalid_argument if the grid is empty or out of range, a time
 * is out of range, or the sources do not lie in one solid region.
 */
Snapshots ComputeSnapshots(const Model& model, const Calculation& calculation,
                           std::span<const Source> sources,
                           const SnapshotGrid& grid,
                           std::span<const double> times,
                           const Progress& progress = {});

}  // namespace YSpec
