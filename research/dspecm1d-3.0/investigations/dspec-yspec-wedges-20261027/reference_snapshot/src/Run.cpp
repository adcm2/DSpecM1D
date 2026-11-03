#include <YSpec/Model.hpp>
#include <YSpec/Run.hpp>
#include <format>
#include <fstream>
#include <limits>
#include <map>
#include <numbers>
#include <span>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

#include "Geometry.hpp"
#include "Sphere.hpp"

#ifdef PLANETMODEL_WITH_NETCDF
#include <planetmodel/modelfile.hpp>
#endif

namespace YSpec {

Seismograms ComputeSeismograms(const SeismogramInput& input,
                               const Progress& progress) {
  const auto model = Model::ReadDeck(input.model);
  return ComputeSeismograms(model, input.calculation, input.sources,
                            input.receivers, progress);
}

namespace {

// The file of the seismograms of a source at a receiver.
std::filesystem::path SeismogramPath(const SeismogramInput& input,
                                     std::size_t sources, std::size_t source,
                                     std::size_t receiver) {
  auto path = input.prefix;
  path += sources == 1
              ? std::format(".{}", receiver + 1)
              : std::format(".{}.{}", input.names.at(source), receiver + 1);
  return path;
}

// A table of the seismogram files, one line for each source and receiver:
// where the receiver is, and where it is seen from the source.
void WriteReceivers(const SeismogramInput& input, std::size_t sources) {
  auto path = input.prefix;
  path += ".receivers";
  auto file = std::ofstream(path);
  if (!file) throw std::runtime_error("cannot write " + path.string());
  file << "# file source latitude longitude depth_km distance_deg "
          "azimuth_deg back_azimuth_deg\n";
  for (auto source = std::size_t{0}; source < sources; ++source) {
    for (auto receiver = std::size_t{0}; receiver < input.receivers.size();
         ++receiver) {
      const auto& at = input.receivers[receiver];
      const auto view = Geometry(at.location, input.sources[source].location);
      const auto degree = std::numbers::pi / 180;
      file << std::format(
          "{} {} {:10.4f}{:10.4f}{:10.3f}{:10.4f}{:10.4f}{:10.4f}\n",
          SeismogramPath(input, sources, source, receiver).filename().string(),
          input.names.at(source), at.location.latitude, at.location.longitude,
          at.depth / 1e3, view.distance / degree,
          (std::numbers::pi - view.phi) / degree,
          (2 * std::numbers::pi - view.rotation) / degree);
    }
  }
  if (!file) throw std::runtime_error("cannot write " + path.string());
}

}  // namespace

void WriteSeismograms(const Seismograms& seismograms,
                      const SeismogramInput& input) {
  WriteReceivers(input, seismograms.Sources());
  for (auto source = std::size_t{0}; source < seismograms.Sources(); ++source) {
    for (auto receiver = std::size_t{0}; receiver < seismograms.Receivers();
         ++receiver) {
      const auto path =
          SeismogramPath(input, seismograms.Sources(), source, receiver);
      auto file = std::ofstream(path);
      if (!file) throw std::runtime_error("cannot write " + path.string());
      const auto vertical = seismograms.Trace(source, receiver, Vertical);
      const auto north = seismograms.Trace(source, receiver, North);
      const auto east = seismograms.Trace(source, receiver, East);
      for (auto i = std::size_t{0}; i < seismograms.Samples(); ++i) {
        file << std::format("{:15.6E}{:15.6E}{:15.6E}{:15.6E}\n",
                            static_cast<double>(i) * seismograms.TimeStep(),
                            vertical[i], north[i], east[i]);
      }
      if (!file) throw std::runtime_error("cannot write " + path.string());
    }
  }
}

Seismograms RunSeismograms(const std::filesystem::path& path,
                           const Progress& progress) {
  const auto input = ReadSeismogramInput(path);
  auto seismograms = ComputeSeismograms(input, progress);
  WriteSeismograms(seismograms, input);
  return seismograms;
}

Snapshots ComputeSnapshots(const SnapshotInput& input,
                           const Progress& progress) {
  const auto model = Model::ReadDeck(input.model);
  return ComputeSnapshots(model, input.calculation, input.sources, input.grid,
                          input.times, progress);
}

void WriteSnapshots(const Snapshots& snapshots, const SnapshotInput& input) {
  const auto& grid = snapshots.Grid();
  const auto quantity = [&] {
    switch (input.calculation.quantity) {
      case Quantity::Displacement:
        return "displacement (m)";
      case Quantity::Velocity:
        return "velocity (m/s)";
      case Quantity::Acceleration:
        return "acceleration (m/s^2)";
    }
    return "";
  }();
  for (auto source = std::size_t{0}; source < snapshots.Sources(); ++source) {
    // Where the points are, around this source.
    auto places = std::vector<Location>{};
    for (auto azimuth : grid.azimuths) {
      for (auto distance : grid.distances) {
        places.push_back(GreatCircle(input.sources[source].location, azimuth,
                                     distance, distance, 1)
                             .front());
      }
    }
    for (auto time = std::size_t{0}; time < snapshots.Times().size(); ++time) {
      auto path = input.prefix;
      path += snapshots.Sources() == 1
                  ? std::format(".{}", time + 1)
                  : std::format(".{}.{}", input.names.at(source), time + 1);
      auto file = std::ofstream(path);
      if (!file) throw std::runtime_error("cannot write " + path.string());
      file << std::format("# source {}, {} s after the origin, {}\n",
                          input.names.at(source), snapshots.Times()[time],
                          quantity);
      file << "# depth_km distance_deg azimuth_deg latitude longitude up away "
              "transverse\n";
      for (auto depth = std::size_t{0}; depth < grid.depths.size(); ++depth) {
        for (auto a = std::size_t{0}; a < grid.azimuths.size(); ++a) {
          for (auto i = std::size_t{0}; i < grid.distances.size(); ++i) {
            const auto& place = places[a * grid.distances.size() + i];
            file << std::format(
                "{:10.3f}{:10.3f}{:10.3f}{:10.4f}{:10.4f}{:15.6E}{:15.6E}"
                "{:15.6E}\n",
                grid.depths[depth] / 1e3, grid.distances[i], grid.azimuths[a],
                place.latitude, place.longitude,
                snapshots(source, time, depth, i, a, Up),
                snapshots(source, time, depth, i, a, Away),
                snapshots(source, time, depth, i, a, Transverse));
          }
        }
      }
      if (!file) throw std::runtime_error("cannot write " + path.string());
    }
  }
}

Snapshots RunSnapshots(const std::filesystem::path& path,
                       const Progress& progress) {
  const auto input = ReadSnapshotInput(path);
  auto snapshots = ComputeSnapshots(input, progress);
  WriteSnapshots(snapshots, input);
  return snapshots;
}

KernelResult ComputeKernels(const KernelInput& input,
                            const StageProgress& progress) {
  auto stage = [&](std::string_view name) {
    return progress ? progress(name) : Progress{};
  };
  const auto model = Model::ReadDeck(input.model);
  const auto& calculation = input.calculation;
  auto result =
      KernelResult{.seismogram = ComputeSeismograms(
                       model, calculation, std::span{&input.source, 1},
                       std::span{&input.receiver, 1}, stage("seismogram")),
                   .measurements = {}};
  auto directions = std::vector<Direction>{};
  auto windows = std::vector<Window>{};
  for (const auto& m : input.measurements) {
    directions.push_back(m.direction);
    windows.push_back(m.window);
  }
  auto measurements = Traveltimes(model, calculation, input.source,
                                  input.receiver, directions, windows);
  auto radial =
      ComputeRadialKernels(model, calculation, input.source, measurements,
                           input.mesh, stage("radial kernels"));

  // The three-dimensional kernels, on every sampling at once.
  auto samplings = std::vector<KernelSampling>{};
  const auto& from = input.source.location;
  const auto& to = input.receiver.location;
  if (input.volume) samplings.push_back(VolumeSampling(model, *input.volume));
  if (input.sections) {
    samplings.push_back(PathSection(model, from, to, *input.sections));
    samplings.push_back(CrossSection(model, from, to, *input.sections));
  }
  if (input.maps) samplings.push_back(PathMaps(model, from, to, *input.maps));
  auto sampled = std::vector<std::vector<SampledKernels>>(measurements.size());
  if (!samplings.empty()) {
    sampled = ComputeSampledKernels(
        model, calculation, input.source, measurements, input.parameters,
        input.scaling, samplings, stage("three-dimensional kernels"));
  }

  for (auto j : std::views::iota(std::size_t{0}, measurements.size())) {
    auto named =
        Reparameterise(radial[j], model, input.parameters, input.scaling);
    auto next = sampled[j].begin();
    auto take = [&](bool wanted) {
      return wanted ? std::optional{std::move(*next++)} : std::nullopt;
    };
    auto volume = take(input.volume.has_value());
    auto path = take(input.sections.has_value());
    auto across = take(input.sections.has_value());
    auto maps = take(input.maps.has_value());
    result.measurements.push_back({.name = input.measurements[j].name,
                                   .measurement = std::move(measurements[j]),
                                   .kernels = std::move(radial[j]),
                                   .named = std::move(named),
                                   .volume = std::move(volume),
                                   .path = std::move(path),
                                   .across = std::move(across),
                                   .maps = std::move(maps)});
  }
  return result;
}

bool WritesModelFiles() {
#ifdef PLANETMODEL_WITH_NETCDF
  return true;
#else
  return false;
#endif
}

namespace {

#ifdef PLANETMODEL_WITH_NETCDF
// The dimensions of a parameter, as exponents of mass, length and time, which
// a model file records.
std::array<int, 3> DimensionsOf(std::string_view name) {
  if (name == "rho") return {1, -3, 0};
  for (auto modulus : {"A", "C", "F", "L", "N", "kappa", "mu"}) {
    if (name == modulus) return {1, -1, -2};
  }
  if (name.starts_with('v')) return {0, 1, -1};
  return {0, 0, 0};
}
#endif

std::string_view QuantityName(Quantity quantity) {
  switch (quantity) {
    case Quantity::Displacement:
      return "displacement";
    case Quantity::Velocity:
      return "velocity";
    case Quantity::Acceleration:
      return "acceleration";
  }
  return "";
}

// What a measurement is, in words.
std::string Described(const KernelInput& input, const KernelMeasurement& m) {
  return std::format(
      "the cross-correlation traveltime (s) of {} {} at the receiver, in a "
      "window from {} s to {} s",
      m.component, QuantityName(input.calculation.quantity), m.window.start,
      m.window.end);
}

}  // namespace

void WriteVolumeKernels(const SampledKernels& kernels, const KernelInput& input,
                        const KernelMeasurement& measurement,
                        const std::filesystem::path& path) {
#ifndef PLANETMODEL_WITH_NETCDF
  (void)kernels;
  (void)input;
  (void)measurement;
  throw std::runtime_error(
      "this YSpec was built without netCDF, and does not write the model "
      "file " +
      path.string() + ": build it with YSPEC_WITH_NETCDF on");
#else
  namespace pm = planetmodel;
  if (!input.volume) {
    throw std::invalid_argument("the kernels are not on a volume grid");
  }
  const auto degree = input.volume->degree;
  const auto cells = static_cast<std::size_t>((degree + 1) * (2 * degree + 1));
  const auto nodes = kernels.Nodes();
  if (kernels.Places().size() != cells + 2) {
    throw std::invalid_argument("the kernels are not on the volume grid");
  }
  // The layers' boundaries, and the elements between the nodes, which are
  // the files' radial nodes, two to an element.
  const auto layers = nodes.back().layer + 1;
  auto boundaries = std::vector<double>{nodes.front().radius};
  auto edges = std::vector<double>{};
  for (auto n : std::views::iota(std::size_t{0}, nodes.size())) {
    if (edges.empty() || nodes[n].radius > edges.back()) {
      edges.push_back(nodes[n].radius);
    }
    if (n + 1 == nodes.size() || nodes[n + 1].layer != nodes[n].layer) {
      boundaries.push_back(nodes[n].radius);
    }
  }
  auto file = pm::ModelFile(pm::Skeleton(boundaries));
  file.SetName("yspec-kernels");
  file.SetMeta("code", std::string("yspec-kernels ") + YSPEC_VERSION);
  file.SetMeta("measurement", Described(input, measurement));
  file.SetMeta("scaling",
               std::string(input.scaling == Scaling::Relative ? "relative"
                                                              : "absolute"));
  file.SetMeta("source_latitude", input.source.location.latitude);
  file.SetMeta("source_longitude", input.source.location.longitude);
  file.SetMeta("source_depth", input.source.depth);
  file.SetMeta("receiver_latitude", input.receiver.location.latitude);
  file.SetMeta("receiver_longitude", input.receiver.location.longitude);
  file.SetMeta("receiver_depth", input.receiver.depth);
  file.AddPart(0, static_cast<int>(layers) - 1, edges, 2,
               pm::GaussLegendre(degree));

  const auto relative = input.scaling == Scaling::Relative;
  for (auto j : std::views::iota(std::size_t{0}, kernels.Names().size())) {
    const auto& name = kernels.Names()[j];
    // Per unit volume of a traveltime, per unit of the parameter if
    // absolute.
    const auto [mass, length, time] =
        relative ? std::array{0, 0, 0} : DimensionsOf(name);
    const auto spec = pm::FieldSpec{
        .character = pm::SCALAR,
        .dimensions = pm::Dimensions{.mass = -mass,
                                     .length = -3 - length,
                                     .time = 1 - time,
                                     .temperature = 0},
        .meaning = std::format(
            "the derivative of the traveltime with respect to {}{}, per unit "
            "volume",
            relative ? "the logarithm of " : "", name)};
    auto sources = std::map<int, pm::LayerSource>{};
    for (auto layer : std::views::iota(Int{0}, layers)) {
      const auto points = file.PointsOf(static_cast<int>(layer));
      auto values = std::vector<double>(points.Size());
      for (auto i : std::views::iota(std::size_t{0}, points.radii.size())) {
        // The node of this layer at the radius.
        const auto node = static_cast<std::size_t>(
            std::ranges::min_element(
                nodes,
                [&](const auto& a, const auto& b) {
                  auto distance = [&](const auto& n) {
                    return n.layer == layer
                               ? std::abs(n.radius - points.radii[i])
                               : std::numeric_limits<double>::infinity();
                  };
                  return distance(a) < distance(b);
                }) -
            nodes.begin());
        const auto at = kernels.Values(j, node);
        std::ranges::copy(at.first(cells), values.begin() + i * cells);
        values[points.radii.size() * cells + 2 * i] = at[cells];
        values[points.radii.size() * cells + 2 * i + 1] = at[cells + 1];
      }
      sources.emplace(static_cast<int>(layer),
                      pm::LayerSource::LateralValues(std::move(values),
                                                     pm::Frame::Spherical));
    }
    file.AddField(
        "K_" + name, spec, sources,
        {.radial = pm::RadialRule::Linear, .lateral = pm::LateralRule::Linear});
  }
  pm::WriteModelFile(path, file);
#endif
}

namespace {

// Writes kernels on a section or a map, a line for each node and place,
// after the section's name if it has one.
void WriteSampled(std::ofstream& file, const SampledKernels& kernels,
                  std::string_view section, double surface) {
  const auto nodes = kernels.Nodes();
  const auto places = kernels.Places();
  for (auto n : std::views::iota(std::size_t{0}, nodes.size())) {
    for (auto p : std::views::iota(std::size_t{0}, places.size())) {
      const auto& place = places[p];
      if (!section.empty()) file << section << " ";
      file << std::format(
          "{:3}{:12.4f}{:11.4f}{:11.4f}{:11.5f}{:12.5f}", nodes[n].layer,
          (surface - nodes[n].radius) / 1e3, place.along, place.across,
          Sphere::Latitude(place.colatitude), place.longitude / Sphere::Degree);
      for (auto j : std::views::iota(std::size_t{0}, kernels.Names().size())) {
        file << std::format("{:14.5E}", kernels.Values(j, n)[p]);
      }
      file << "\n";
    }
  }
}

}  // namespace

void WriteKernels(const KernelResult& result, const KernelInput& input) {
  const auto model = Model::ReadDeck(input.model);
  const auto surface = model.SurfaceRadius() * model.GetScales().Length();
  const auto relative = input.scaling == Scaling::Relative;
  const auto several = input.measurements.size() > 1;
  for (auto j : std::views::iota(std::size_t{0}, result.measurements.size())) {
    const auto& kernel = result.measurements[j];
    const auto& m = input.measurements[j];
    const auto measured = Described(input, m);
    auto open = [&](std::string_view suffix) {
      auto path = input.prefix;
      if (several) path += "." + m.name;
      path += suffix;
      auto file = std::ofstream(path);
      if (!file) throw std::runtime_error("cannot write " + path.string());
      return std::pair{std::move(file), path};
    };
    auto check = [](const std::ofstream& file,
                    const std::filesystem::path& path) {
      if (!file) throw std::runtime_error("cannot write " + path.string());
    };

    {
      auto [file, path] = open(".measurement");
      file << "# yspec-kernels: " << measured << "\n";
      file << "# time_s seismogram window adjoint_source\n";
      const auto& d = m.direction;
      const auto& u = result.seismogram;
      const auto& adjoint = kernel.measurement.adjointSource;
      for (auto n = std::size_t{0}; n < u.Samples(); ++n) {
        const auto t = static_cast<double>(n) * u.TimeStep();
        file << std::format("{:12.3f}{:15.6E}{:10.4f}{:15.6E}\n", t,
                            d.up * u(0, 0, Vertical, n) +
                                d.north * u(0, 0, North, n) +
                                d.east * u(0, 0, East, n),
                            m.window(t), adjoint[n]);
      }
      check(file, path);
    }

    {
      auto [file, path] = open(".kernels");
      const auto& named = kernel.named;
      file << "# yspec-kernels: radial kernels of " << measured << "\n";
      file << std::format(
          "# dJ = integral of sum over X of K_X d{}X r^2 dr, plus sum over X "
          "of P_X d{}X at the source; K in s m^-3{}, P in s{}\n",
          relative ? "ln " : "", relative ? "ln " : "",
          relative ? "" : " per unit of X", relative ? "" : " per unit of X");
      file << std::format(
          "# the point parts P_X at the source, at {:.3f} km:\n",
          kernel.kernels.SourceRadius() / 1e3);
      file << "# at_source";
      for (auto value : named.atSource) file << std::format(" {:.6E}", value);
      file << "\n# radius_km depth_km weight_m3";
      for (const auto& name : named.names) file << " " << name;
      file << "\n";
      const auto& radial = kernel.kernels;
      for (auto i = std::size_t{0}; i < radial.Size(); ++i) {
        const auto r = radial.Radii()[i];
        file << std::format("{:15.6f}{:15.6f}{:15.6E}", r / 1e3,
                            (surface - r) / 1e3, radial.Weights()[i]);
        for (const auto& values : named.values) {
          file << std::format("{:15.6E}", values[i]);
        }
        file << "\n";
      }
      check(file, path);
    }

    auto header = [&](std::ofstream& file, const SampledKernels& kernels,
                      std::string_view what, std::string_view first) {
      file << "# yspec-kernels: " << what << " of the kernels of " << measured
           << "\n";
      file << std::format("# per unit volume of {}X, in s m^-3{}\n",
                          relative ? "ln " : "",
                          relative ? "" : " per unit of X");
      file << "# " << first
           << "layer depth_km along_deg across_deg latitude longitude";
      for (const auto& name : kernels.Names()) file << " " << name;
      file << "\n";
    };
    if (kernel.path && kernel.across) {
      auto [file, path] = open(".sections");
      header(file, *kernel.path, "sections along and across the path",
             "section ");
      WriteSampled(file, *kernel.path, "path", surface);
      WriteSampled(file, *kernel.across, "across", surface);
      check(file, path);
    }
    if (kernel.maps) {
      auto [file, path] = open(".maps");
      header(file, *kernel.maps, "maps", "");
      WriteSampled(file, *kernel.maps, "", surface);
      check(file, path);
    }
    if (kernel.volume) {
      auto path = input.prefix;
      if (several) path += "." + m.name;
      path += ".nc";
      WriteVolumeKernels(*kernel.volume, input, m, path);
    }
  }
}

KernelResult RunKernels(const std::filesystem::path& path,
                        const StageProgress& progress) {
  const auto input = ReadKernelInput(path);
  if (input.volume && !WritesModelFiles()) {
    throw std::runtime_error(
        "three-dimensional kernels through the earth are written as model "
        "files, and this YSpec was built without netCDF: build it with "
        "YSPEC_WITH_NETCDF on, or ask for sections and maps");
  }
  auto result = ComputeKernels(input, progress);
  WriteKernels(result, input);
  return result;
}

}  // namespace YSpec
