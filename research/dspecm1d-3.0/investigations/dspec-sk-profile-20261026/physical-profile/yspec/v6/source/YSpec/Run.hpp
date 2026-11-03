#pragma once

/**
 * @file Run.hpp
 * @brief Calculations described by input files, from end to end: what the
 * programs do.
 */

#include <filesystem>
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "Input.hpp"
#include "Kernels.hpp"
#include "Seismograms.hpp"
#include "Snapshots.hpp"
#include "Synthetics.hpp"

namespace YSpec {

/**
 * @brief Reads the model an input names and computes its seismograms.
 * @param input The calculation.
 * @param progress Called as frequencies are completed, if given.
 */
Seismograms ComputeSeismograms(const SeismogramInput& input,
                               const Progress& progress = {});

/**
 * @brief Writes seismograms as text, one file per source and receiver,
 * named after the input's output prefix: `prefix.1`, `prefix.2`, ... in the
 * order of the receivers when there is one source, and `prefix.name.1`, ...,
 * with the source's name, when there are several. Each line holds the time
 * after the origin (s) and the vertical, north and east components.
 * @param seismograms The seismograms.
 * @param input The calculation they are for, which gives the prefix and the
 * names of the sources.
 * @throws std::runtime_error if a file cannot be written.
 */
void WriteSeismograms(const Seismograms& seismograms,
                      const SeismogramInput& input);

/**
 * @brief Reads an input file, computes the seismograms it describes and
 * writes them.
 * @param path The input file.
 * @param progress Called as frequencies are completed, if given.
 * @return The seismograms.
 */
Seismograms RunSeismograms(const std::filesystem::path& path,
                           const Progress& progress = {});

/**
 * @brief Reads the model an input names and computes its snapshots.
 * @param input The calculation.
 * @param progress Called as frequencies are completed, if given.
 */
Snapshots ComputeSnapshots(const SnapshotInput& input,
                           const Progress& progress = {});

/**
 * @brief Writes snapshots as text, one file per source and time, named
 * `prefix.1`, `prefix.2`, ... in the order of the times when there is one
 * source, and `prefix.name.1`, ..., with the source's name, when there are
 * several. After two comment lines starting with '#', each line holds a point
 * of the grid, depth by depth, azimuth by azimuth and distance by distance:
 * its depth (km), distance (degrees) and azimuth (degrees) from the source,
 * its latitude and longitude, and the up, away and transverse components (see
 * SnapshotComponent).
 * @param snapshots The snapshots.
 * @param input The calculation they are for.
 * @throws std::runtime_error if a file cannot be written.
 */
void WriteSnapshots(const Snapshots& snapshots, const SnapshotInput& input);

/**
 * @brief Reads an input file, computes the snapshots it describes and writes
 * them.
 * @param path The input file.
 * @param progress Called as frequencies are completed, if given.
 * @return The snapshots.
 */
Snapshots RunSnapshots(const std::filesystem::path& path,
                       const Progress& progress = {});

/** @brief The kernels of one of the measurements of yspec-kernels. */
struct MeasurementKernels {
  std::string name;         ///< Its name.
  Measurement measurement;  ///< The traveltime, by its adjoint source.
  RadialKernels kernels;    ///< Its radial kernels in the model's parameters.
  NamedKernels named;       ///< And in those asked for.
  std::optional<SampledKernels> volume;  ///< Through the earth, if asked for.
  std::optional<SampledKernels> path;    ///< On the section along the path.
  std::optional<SampledKernels> across;  ///< On the section across it.
  std::optional<SampledKernels> maps;    ///< On maps around it.
};

/** @brief What the yspec-kernels program computes. */
struct KernelResult {
  Seismograms seismogram;  ///< The seismogram at the receiver.
  std::vector<MeasurementKernels> measurements;  ///< The kernels of each.
};

/** @brief A progress report for each stage of a calculation, by its name. */
using StageProgress = std::function<Progress(std::string_view stage)>;

/**
 * @brief Reads the model an input names and computes the kernels of its
 * traveltimes.
 * @param input The calculation.
 * @param progress Gives a report for each stage, as frequencies are done.
 */
KernelResult ComputeKernels(const KernelInput& input,
                            const StageProgress& progress = {});

/**
 * @brief Whether this YSpec writes model files, and so three-dimensional
 * kernels through the earth: it does when built with netCDF
 * (YSPEC_WITH_NETCDF).
 */
bool WritesModelFiles();

/**
 * @brief Writes three-dimensional kernels on VolumeSampling as a model file
 * of planetmodel (schema planetmodel.model/1), on the skeleton of the
 * model's discontinuities, with a field `K_X` for each parameter, read
 * linearly between the nodes in radius and in angle.
 * @param kernels The kernels.
 * @param input What they are of, which the file's metadata record.
 * @param measurement Which measurement they are of.
 * @param path The file.
 * @throws std::runtime_error if this YSpec does not write model files, and
 * planetmodel::Error if the file cannot be written.
 */
void WriteVolumeKernels(const SampledKernels& kernels, const KernelInput& input,
                        const KernelMeasurement& measurement,
                        const std::filesystem::path& path);

/**
 * @brief Writes the kernels of each measurement and the measurement as text,
 * with the prefix, or with the prefix and the measurement's name when there
 * are several: `.kernels`, whose header gives the measurement, the units and
 * the point parts at the source (on a line beginning `# at_source`), and
 * whose lines give the radius and depth (km), the weight (m^3) and each
 * kernel, in the order of the column names on the last line of the header;
 * `.measurement`, whose lines give the time (s), the seismogram measured, the
 * window and the adjoint source; with sections, `.sections`, whose lines
 * give the section (`path` or `across`), the layer, the depth (km), the
 * distances along and across the path and the latitude and longitude
 * (degrees), and each kernel; with maps, `.maps`, the same without the
 * section; and through the earth, `.nc` (see WriteVolumeKernels).
 * @throws std::runtime_error if a file cannot be written.
 */
void WriteKernels(const KernelResult& result, const KernelInput& input);

/**
 * @brief Reads an input file, computes the kernels it describes and writes
 * them.
 */
KernelResult RunKernels(const std::filesystem::path& path,
                        const StageProgress& progress = {});

}  // namespace YSpec
