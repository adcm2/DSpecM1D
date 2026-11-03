#pragma once

/**
 * @file Input.hpp
 * @brief Reads the input files of the programs.
 *
 * @details The files are strict JSON (RFC 8259). The programs share most of
 * their members: those of the yspec program, which computes seismograms, are
 * shown here, with the units in their names; "output" and "receivers" are its
 * own, and the rest are shared.
 *
 * @code{.json}
 * {
 *   "model": "prem.200",
 *   "attenuation": true,
 *   "gravity": "full",
 *   "output": {
 *     "prefix": "yspec.out",
 *     "quantity": "velocity",
 *     "potential_and_tilt_corrections": false
 *   },
 *   "degrees": { "min": 0, "max": 600 },
 *   "frequencies_mHz": { "min": 0.2, "max": 40 },
 *   "time_series": { "length_min": 60, "time_step_s": 1.0 },
 *   "filter_mHz": [1.0, 5.0, 38, 40],
 *   "sources": [
 *     { "name": "event", "depth_km": 29.5, "latitude": 36.87,
 *       "longitude": 35.58,
 *       "moment_tensor_Nm": { "rr": 4.22e17, "rt": 3.08e17, "rp": 6.99e17,
 *                             "tt": -2.976e18, "tp": 6.54e17, "pp": 2.55e18 }
 *     }
 *   ],
 *   "receivers": [
 *     { "depth_km": 3, "locations": [[0, 4], [0, 6]],
 *       "instrument": "ANMO.pz" },
 *     { "depth_km": [3, 100],
 *       "profile": { "from": [0, 0], "to": [0, 90], "count": 10 } },
 *     { "depth_km": 3,
 *       "profile": { "from": [0, 0], "azimuth_deg": 45,
 *                    "distances_deg": [10, 170], "count": 17 } }
 *   ],
 *   "solver": { "tolerance": 1e-5, "start_decay": 11 }
 * }
 * @endcode
 *
 * "gravity" is one of "none", "cowling" and "full"; "quantity" one of
 * "displacement", "velocity" and "acceleration". The moment tensor is in
 * spherical coordinates at the source (r up, theta south, phi east). A source
 * may have instead, or as well, a force, "force_N": { "r": 0, "t": 0,
 * "p": 1e12 }, in the same coordinates, whose members default to zero.
 * Locations are [latitude, longitude] in geographic degrees. The filter is a
 * raised-cosine taper rising between its first two corners and falling
 * between its last two; it may be left out. "attenuation", "gravity", the
 * corrections, "degrees.min", the filter, the sources' names and "solver"
 * have defaults (true, "full", false, 0, none, their numbers from 1, and the
 * solver's defaults).
 *
 * "receivers" is one group or an array of them. A group has one depth or
 * several, and either a list of locations or a profile along a great circle
 * (see GreatCircle), from one point to another or from a point at an azimuth
 * between two distances, with a number of points equally spaced. The
 * receivers are numbered group by group and, within a group, depth by depth.
 * A group's optional instrument is a SAC pole-zero file (see
 * ReadSacPoleZero).
 *
 * The yspec-snapshot program has "output" and "snapshots" instead:
 *
 * @code{.json}
 * {
 *   "output": { "prefix": "snapshot", "quantity": "displacement" },
 *   "snapshots": {
 *     "times_s": [300, 600, 900],
 *     "depths_km": { "min": 25, "max": 2850, "count": 114 },
 *     "distances_deg": { "min": 0, "max": 180, "count": 361 },
 *     "azimuths_deg": [30, 210]
 *   }
 * }
 * @endcode
 *
 * The yspec-kernels program has one source, and "output", "receiver",
 * "measurements" and "kernels" instead:
 *
 * @code{.json}
 * {
 *   "output": { "prefix": "kernels", "quantity": "velocity" },
 *   "receiver": { "depth_km": 3, "distance_deg": 60, "azimuth_deg": 30 },
 *   "measurements": [
 *     { "name": "S", "component": "T",
 *       "traveltime": { "window_s": [1000, 1200], "taper_s": 50 } },
 *     { "name": "P", "component": "Z",
 *       "traveltime": { "window_s": [580, 640], "taper_s": 20 } }
 *   ],
 *   "kernels": {
 *     "parameters": "velocities",
 *     "scaling": "relative",
 *     "mesh": { "order": 5, "max_element_km": 50 },
 *     "volume": { "degree": 64, "max_element_km": 100 },
 *     "sections": { "margin_deg": 15, "half_width_deg": 15,
 *                   "spacing_deg": 0.2, "radial_spacing_km": 25 },
 *     "maps": { "depths_km": [600, 1500], "margin_deg": 10,
 *               "half_width_deg": 15, "spacing_deg": 0.5 }
 *   }
 * }
 * @endcode
 *
 * The receiver is at a "location", [latitude, longitude], or at a distance
 * and azimuth from the source, and may have an "instrument". Each
 * measurement, with a name to label its files (by default its number), is
 * the traveltime of a component, one of "Z", "N", "E", "R" (radial, away from
 * the source) and "T" (transverse, 90 degrees anticlockwise from radial seen
 * from above), measured by cross-correlation in a window with cosine tapers
 * of the given length inside it (default none). The measurements share the
 * receiver, and with it most of the work. "parameters" is one of "love",
 * "moduli", "velocities" (the default) and "anisotropic" (see ParameterSet),
 * and "scaling" "absolute" (the default) or "relative"; the mesh has the
 * defaults of RadialMesh. Three-dimensional kernels, for the same parameters
 * without density, are computed on whatever of "volume" (see VolumeGrid),
 * "sections" (see SectionOptions) and "maps" (see MapOptions) is given, each
 * with the defaults of its structure.
 *
 * Each of the times, depths, distances and azimuths is one number, an array
 * of them, or a range of "count" values equally spaced from "min" to "max",
 * both included. The grid is in the frame of each source: see
 * ComputeSnapshots. "time_series.time_step_s" may be left out, by either
 * program: it is then a quarter of the period of the highest frequency.
 *
 * Relative paths are relative to the directory of the input file. A member
 * not listed here is an error, so that a misspelt one is not silently
 * ignored.
 */

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "Calculation.hpp"
#include "Kernels.hpp"
#include "Receiver.hpp"
#include "Snapshots.hpp"
#include "Source.hpp"

namespace YSpec {

/** @brief The input of the yspec program, in SI units. */
struct SeismogramInput {
  std::filesystem::path model;      ///< The model deck.
  Calculation calculation;          ///< What to compute.
  std::vector<Source> sources;      ///< The sources.
  std::vector<std::string> names;   ///< The names of the sources.
  std::vector<Receiver> receivers;  ///< The receivers.
  std::filesystem::path prefix;     ///< The prefix of the output files.
};

/**
 * @brief Reads the input of the yspec program.
 * @throws std::runtime_error naming the file, the member and the fault.
 */
SeismogramInput ReadSeismogramInput(const std::filesystem::path& path);

/** @brief The input of the yspec-snapshot program, in SI units. */
struct SnapshotInput {
  std::filesystem::path model;     ///< The model deck.
  Calculation calculation;         ///< What to compute.
  std::vector<Source> sources;     ///< The sources.
  std::vector<std::string> names;  ///< The names of the sources.
  SnapshotGrid grid;               ///< The points, around each source.
  std::vector<double> times;       ///< The times after the origin (s).
  std::filesystem::path prefix;    ///< The prefix of the output files.
};

/**
 * @brief Reads the input of the yspec-snapshot program.
 * @throws std::runtime_error naming the file, the member and the fault.
 */
SnapshotInput ReadSnapshotInput(const std::filesystem::path& path);

/** @brief The input of the yspec-kernels program, in SI units. */
/** @brief A traveltime on a component, as yspec-kernels reads it. */
struct KernelMeasurement {
  std::string name;       ///< Its name, which labels its files.
  std::string component;  ///< "Z", "N", "E", "R" or "T".
  Direction direction;    ///< The direction of that component.
  Window window;          ///< The window of the traveltime.
};

/** @brief The input of the yspec-kernels program, in SI units. */
struct KernelInput {
  std::filesystem::path model;  ///< The model deck.
  Calculation calculation;      ///< How the seismogram is computed.
  Source source;                ///< The source.
  std::string name;             ///< The name of the source.
  Receiver receiver;            ///< The receiver.
  std::vector<KernelMeasurement> measurements;         ///< The measurements.
  ParameterSet parameters = ParameterSet::Velocities;  ///< The kernels wanted.
  Scaling scaling = Scaling::Absolute;  ///< Absolute or relative.
  RadialMesh mesh;                      ///< How finely.
  /// The grid of three-dimensional kernels through the earth, if wanted.
  std::optional<VolumeGrid> volume;
  /// The sections along and across the path, if wanted.
  std::optional<SectionOptions> sections;
  /// The maps around the path, if wanted.
  std::optional<MapOptions> maps;
  std::filesystem::path prefix;  ///< The prefix of the output files.
};

/**
 * @brief Reads the input of the yspec-kernels program.
 * @throws std::runtime_error naming the file, the member and the fault.
 */
KernelInput ReadKernelInput(const std::filesystem::path& path);

}  // namespace YSpec
