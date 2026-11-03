// The Python module yspec._core: the library as it is, in SI units, with
// numpy arrays for its results. The package around it, python/yspec, gives
// the interface that is documented; this one stays close to the C++, so that
// each binding is a line or two.
//
// The calculations release the interpreter's lock, so that the library's
// threads run while Python waits. Their progress comes back through a Python
// callable, called with the lock taken again, which also checks for signals:
// Ctrl-C stops a calculation at the next frequency to finish.

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/filesystem.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <YSpec/YSpec.hpp>
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <optional>
#include <ranges>
#include <span>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace nb = nanobind;
using namespace nb::literals;

namespace {

template <class T, std::size_t N>
using Array = nb::ndarray<nb::numpy, T, nb::ndim<N>>;

// A numpy array holding a copy of the values.
template <class T = double, std::ranges::sized_range Range>
Array<T, 1> Copy(const Range& values) {
  const auto size = static_cast<std::size_t>(std::ranges::size(values));
  auto* data = new T[size];
  std::ranges::copy(values, data);
  auto owner =
      nb::capsule(data, [](void* p) noexcept { delete[] static_cast<T*>(p); });
  return Array<T, 1>(data, {size}, owner);
}

// The library's progress from a Python callable, or None. The callable is
// not owned: the argument it came in keeps it alive while the calculation
// runs.
YSpec::Progress Reporting(nb::handle callback) {
  return [callback](std::size_t done, std::size_t total) {
    auto gil = nb::gil_scoped_acquire{};
    if (PyErr_CheckSignals() != 0) throw nb::python_error();
    if (!callback.is_none()) callback(done, total);
  };
}

// The same for the stages of the kernels, the callable taking the stage too.
YSpec::StageProgress StageReporting(nb::handle callback) {
  return [callback](std::string_view stage) -> YSpec::Progress {
    return [callback, stage = std::string(stage)](std::size_t done,
                                                  std::size_t total) {
      auto gil = nb::gil_scoped_acquire{};
      if (PyErr_CheckSignals() != 0) throw nb::python_error();
      if (!callback.is_none()) callback(stage, done, total);
    };
  };
}

// The traces as an array (source, receiver, component, sample) that views
// their storage.
template <class Traces>
Array<double, 4> TraceView(Traces& traces, std::size_t components) {
  auto* data = traces.Trace(0, 0, {}).data();
  return Array<double, 4>(data, {traces.Sources(), traces.Receivers(),
                                 components, traces.Samples()});
}

// The model's material at a radius (m) in a layer, in SI units.
nb::dict MaterialAt(const YSpec::Model& model, YSpec::Int layer,
                    double radius) {
  const auto& scales = model.GetScales();
  const auto m = model(layer, radius / scales.Length());
  const auto modulus = scales.Modulus();
  auto values = nb::dict{};
  values["rho"] = m.rho * scales.Density();
  values["A"] = m.A * modulus;
  values["C"] = m.C * modulus;
  values["F"] = m.F * modulus;
  values["L"] = m.L * modulus;
  values["N"] = m.N * modulus;
  values["kappa"] = m.kappa * modulus;
  values["mu"] = m.mu * modulus;
  values["q_kappa"] = m.qKappa;
  values["q_mu"] = m.qMu;
  values["g"] = m.g * scales.Acceleration();
  return values;
}

// The nodes and places of a sampling, as arrays, for a class T whose nodes
// and places the functions given return.
template <class T, class Nodes, class Places>
void BindSamplingArrays(nb::class_<T>& cls, Nodes nodes, Places places) {
  auto node = [nodes](auto member) {
    return [nodes, member](const T& s) {
      using Value = std::remove_cvref_t<decltype(nodes(s).front().*member)>;
      return Copy<Value>(
          nodes(s) |
          std::views::transform([member](const auto& n) { return n.*member; }));
    };
  };
  auto place = [places](auto member) {
    return [places, member](const T& s) {
      return Copy(places(s) | std::views::transform([member](const auto& p) {
                    return p.*member;
                  }));
    };
  };
  constexpr auto copy = nb::rv_policy::move;
  cls.def_prop_ro("layers", node(&YSpec::KernelNode::layer), copy)
      .def_prop_ro("radii", node(&YSpec::KernelNode::radius), copy)
      .def_prop_ro("colatitudes", place(&YSpec::KernelPlace::colatitude), copy)
      .def_prop_ro("longitudes", place(&YSpec::KernelPlace::longitude), copy)
      .def_prop_ro("along", place(&YSpec::KernelPlace::along), copy)
      .def_prop_ro("across", place(&YSpec::KernelPlace::across), copy);
}

}  // namespace

NB_MODULE(_core, m) {
  using namespace YSpec;
  m.doc() = "YSpec's library, in SI units; see the yspec package.";
  m.attr("__version__") = YSPEC_VERSION;

  m.def("threads", [] {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
  });
  m.def(
      "set_threads",
      [](int count) {
#ifdef _OPENMP
        omp_set_num_threads(count);
#else
        static_cast<void>(count);
#endif
      },
      "count"_a);
  m.def("writes_model_files", &WritesModelFiles);

  // Places, sources and receivers.
  nb::class_<Location>(m, "Location")
      .def(
          "__init__",
          [](Location* self, double latitude, double longitude) {
            new (self) Location{latitude, longitude};
          },
          "latitude"_a = 0.0, "longitude"_a = 0.0)
      .def_rw("latitude", &Location::latitude)
      .def_rw("longitude", &Location::longitude);
  m.def("great_circle",
        nb::overload_cast<const Location&, const Location&, std::size_t>(
            &GreatCircle),
        "start"_a, "end"_a, "count"_a);
  m.def("great_circle",
        nb::overload_cast<const Location&, double, double, double, std::size_t>(
            &GreatCircle),
        "start"_a, "azimuth"_a, "min_distance"_a, "max_distance"_a, "count"_a);

  nb::class_<MomentTensor>(m, "MomentTensor")
      .def(
          "__init__",
          [](MomentTensor* self, double rr, double rt, double rp, double tt,
             double tp,
             double pp) { new (self) MomentTensor{rr, rt, rp, tt, tp, pp}; },
          "rr"_a = 0.0, "rt"_a = 0.0, "rp"_a = 0.0, "tt"_a = 0.0, "tp"_a = 0.0,
          "pp"_a = 0.0)
      .def_rw("rr", &MomentTensor::rr)
      .def_rw("rt", &MomentTensor::rt)
      .def_rw("rp", &MomentTensor::rp)
      .def_rw("tt", &MomentTensor::tt)
      .def_rw("tp", &MomentTensor::tp)
      .def_rw("pp", &MomentTensor::pp);
  nb::class_<Force>(m, "Force")
      .def(
          "__init__",
          [](Force* self, double r, double t, double p) {
            new (self) Force{r, t, p};
          },
          "r"_a = 0.0, "t"_a = 0.0, "p"_a = 0.0)
      .def_rw("r", &Force::r)
      .def_rw("t", &Force::t)
      .def_rw("p", &Force::p);
  nb::class_<Source>(m, "Source")
      .def(
          "__init__",
          [](Source* self, const Location& location, double depth,
             const MomentTensor& moment, const Force& force) {
            new (self) Source{location, depth, moment, force};
          },
          "location"_a, "depth"_a, "moment"_a = MomentTensor{},
          "force"_a = Force{})
      .def_rw("location", &Source::location)
      .def_rw("depth", &Source::depth)
      .def_rw("moment", &Source::moment)
      .def_rw("force", &Source::force);

  nb::enum_<Quantity>(m, "Quantity")
      .value("displacement", Quantity::Displacement)
      .value("velocity", Quantity::Velocity)
      .value("acceleration", Quantity::Acceleration);
  nb::class_<Instrument>(m, "Instrument")
      .def(
          "__init__",
          [](Instrument* self, std::vector<Complex> zeros,
             std::vector<Complex> poles, double constant, Quantity input) {
            new (self)
                Instrument{std::move(zeros), std::move(poles), constant, input};
          },
          "zeros"_a, "poles"_a, "constant"_a = 1.0,
          "input"_a = Quantity::Displacement)
      .def_rw("zeros", &Instrument::zeros)
      .def_rw("poles", &Instrument::poles)
      .def_rw("constant", &Instrument::constant)
      .def_rw("input", &Instrument::input)
      .def("response", &Instrument::Response, "omega"_a);
  m.def("read_sac_pole_zero", &ReadSacPoleZero, "path"_a);
  nb::class_<Receiver>(m, "Receiver")
      .def(
          "__init__",
          [](Receiver* self, const Location& location, double depth,
             std::optional<Instrument> instrument) {
            new (self) Receiver{location, depth, std::move(instrument)};
          },
          "location"_a, "depth"_a = 0.0, "instrument"_a = nb::none())
      .def_rw("location", &Receiver::location)
      .def_rw("depth", &Receiver::depth)
      .def_rw("instrument", &Receiver::instrument);

  // What is computed.
  nb::enum_<Gravity>(m, "Gravity")
      .value("none", Gravity::None)
      .value("cowling", Gravity::Cowling)
      .value("full", Gravity::Full);
  nb::class_<Degrees>(m, "Degrees")
      .def(
          "__init__",
          [](Degrees* self, int min, int max) { new (self) Degrees{min, max}; },
          "min"_a = 0, "max"_a = 0)
      .def_rw("min", &Degrees::min)
      .def_rw("max", &Degrees::max);
  nb::class_<FrequencyBand>(m, "FrequencyBand")
      .def(
          "__init__",
          [](FrequencyBand* self, double min, double max) {
            new (self) FrequencyBand{min, max};
          },
          "min"_a = 0.0, "max"_a = 0.0)
      .def_rw("min", &FrequencyBand::min)
      .def_rw("max", &FrequencyBand::max);
  nb::class_<Sampling>(m, "Sampling")
      .def(
          "__init__",
          [](Sampling* self, double duration, double timeStep) {
            new (self) Sampling{duration, timeStep};
          },
          "duration"_a = 0.0, "time_step"_a = 0.0)
      .def_rw("duration", &Sampling::duration)
      .def_rw("time_step", &Sampling::timeStep);
  nb::class_<Taper>(m, "Taper")
      .def(
          "__init__",
          [](Taper* self, double f1, double f2, double f3, double f4) {
            new (self) Taper{f1, f2, f3, f4};
          },
          "f1"_a, "f2"_a, "f3"_a, "f4"_a)
      .def_rw("f1", &Taper::f1)
      .def_rw("f2", &Taper::f2)
      .def_rw("f3", &Taper::f3)
      .def_rw("f4", &Taper::f4);
  nb::class_<SolverOptions>(m, "SolverOptions")
      .def(
          "__init__",
          [](SolverOptions* self, double tolerance, double startDecay,
             double stepFactor) {
            new (self) SolverOptions{tolerance, startDecay, stepFactor};
          },
          "tolerance"_a = 1e-5, "start_decay"_a = 11.0, "step_factor"_a = 0.25)
      .def_rw("tolerance", &SolverOptions::tolerance)
      .def_rw("start_decay", &SolverOptions::startDecay)
      .def_rw("step_factor", &SolverOptions::stepFactor);
  nb::class_<Calculation>(m, "Calculation")
      .def(nb::init<>())
      .def_rw("gravity", &Calculation::gravity)
      .def_rw("attenuation", &Calculation::attenuation)
      .def_rw("quantity", &Calculation::quantity)
      .def_rw("corrections", &Calculation::corrections)
      .def_rw("degrees", &Calculation::degrees)
      .def_rw("band", &Calculation::band)
      .def_rw("sampling", &Calculation::sampling)
      .def_rw("taper", &Calculation::taper)
      .def_rw("solver", &Calculation::solver);

  // The model.
  nb::enum_<Side>(m, "Side")
      .value("below", Side::Below)
      .value("above", Side::Above);
  nb::class_<Model>(m, "Model")
      .def_static(
          "read_deck",
          [](const std::filesystem::path& path) {
            return Model::ReadDeck(path);
          },
          "path"_a)
      .def_static(
          "read_deck_text",
          [](const std::string& text) { return Model::ReadDeckText(text); },
          "text"_a)
      .def_prop_ro("number_of_layers", &Model::NumberOfLayers)
      .def_prop_ro("surface_radius",
                   [](const Model& model) {
                     return model.SurfaceRadius() * model.GetScales().Length();
                   })
      .def(
          "layer_bounds",
          [](const Model& model, Int layer) {
            const auto& l = model.GetLayer(layer);
            const auto length = model.GetScales().Length();
            return std::pair{l.Lower() * length, l.Upper() * length};
          },
          "layer"_a)
      .def(
          "is_fluid",
          [](const Model& model, Int layer) {
            return model.GetLayer(layer).IsFluid();
          },
          "layer"_a)
      .def(
          "layer_index",
          [](const Model& model, double radius, Side side) {
            return model.LayerIndex(radius / model.GetScales().Length(), side);
          },
          "radius"_a, "side"_a = Side::Below)
      .def("material", &MaterialAt, "layer"_a, "radius"_a);

  // Seismograms and strains.
  nb::enum_<Component>(m, "Component")
      .value("vertical", Vertical)
      .value("north", North)
      .value("east", East);
  nb::class_<Seismograms>(m, "Seismograms")
      .def_prop_ro("sources", &Seismograms::Sources)
      .def_prop_ro("receivers", &Seismograms::Receivers)
      .def_prop_ro("samples", &Seismograms::Samples)
      .def_prop_ro("time_step", &Seismograms::TimeStep)
      .def_prop_ro(
          "data", [](Seismograms& s) { return TraceView(s, 3); },
          nb::rv_policy::reference_internal);
  nb::class_<Strains>(m, "Strains")
      .def_prop_ro("sources", &Strains::Sources)
      .def_prop_ro("receivers", &Strains::Receivers)
      .def_prop_ro("samples", &Strains::Samples)
      .def_prop_ro("time_step", &Strains::TimeStep)
      .def_prop_ro(
          "data", [](Strains& s) { return TraceView(s, 6); },
          nb::rv_policy::reference_internal);
  m.def(
      "compute_seismograms",
      [](const Model& model, const Calculation& calculation,
         const std::vector<Source>& sources,
         const std::vector<Receiver>& receivers, nb::handle progress) {
        const auto reporter = Reporting(progress);
        auto release = nb::gil_scoped_release{};
        return ComputeSeismograms(model, calculation, sources, receivers,
                                  reporter);
      },
      "model"_a, "calculation"_a, "sources"_a, "receivers"_a,
      "progress"_a = nb::none());
  m.def(
      "compute_strains",
      [](const Model& model, const Calculation& calculation,
         const std::vector<Source>& sources,
         const std::vector<Receiver>& receivers, nb::handle progress) {
        const auto reporter = Reporting(progress);
        auto release = nb::gil_scoped_release{};
        return ComputeStrains(model, calculation, sources, receivers, reporter);
      },
      "model"_a, "calculation"_a, "sources"_a, "receivers"_a,
      "progress"_a = nb::none());

  // Snapshots.
  nb::class_<SnapshotGrid>(m, "SnapshotGrid")
      .def(
          "__init__",
          [](SnapshotGrid* self, std::vector<double> depths,
             std::vector<double> distances, std::vector<double> azimuths) {
            new (self) SnapshotGrid{std::move(depths), std::move(distances),
                                    std::move(azimuths)};
          },
          "depths"_a, "distances"_a, "azimuths"_a)
      .def_rw("depths", &SnapshotGrid::depths)
      .def_rw("distances", &SnapshotGrid::distances)
      .def_rw("azimuths", &SnapshotGrid::azimuths);
  nb::class_<Snapshots>(m, "Snapshots")
      .def_prop_ro("sources", &Snapshots::Sources)
      .def_prop_ro("times", &Snapshots::Times)
      .def_prop_ro("grid", &Snapshots::Grid)
      .def_prop_ro(
          "data",
          // (source, time, depth, azimuth, distance, component)
          [](Snapshots& s) {
            const auto& grid = s.Grid();
            const auto shape = std::array{
                s.Sources(),          s.Times().size(),      grid.depths.size(),
                grid.azimuths.size(), grid.distances.size(), std::size_t{3}};
            const auto empty =
                std::ranges::find(shape, std::size_t{0}) != shape.end();
            auto* data = empty ? nullptr : &s(0, 0, 0, 0, 0, Up);
            return Array<double, 6>(data, 6, shape.data());
          },
          nb::rv_policy::reference_internal);
  m.def(
      "compute_snapshots",
      [](const Model& model, const Calculation& calculation,
         const std::vector<Source>& sources, const SnapshotGrid& grid,
         const std::vector<double>& times, nb::handle progress) {
        const auto reporter = Reporting(progress);
        auto release = nb::gil_scoped_release{};
        return ComputeSnapshots(model, calculation, sources, grid, times,
                                reporter);
      },
      "model"_a, "calculation"_a, "sources"_a, "grid"_a, "times"_a,
      "progress"_a = nb::none());

  // Measurements.
  nb::class_<Window>(m, "Window")
      .def(
          "__init__",
          [](Window* self, double start, double end, double taper) {
            new (self) Window{start, end, taper};
          },
          "start"_a, "end"_a, "taper"_a = 0.0)
      .def_rw("start", &Window::start)
      .def_rw("end", &Window::end)
      .def_rw("taper", &Window::taper)
      .def("__call__", &Window::operator(), "t"_a);
  nb::class_<Direction>(m, "Direction")
      .def(
          "__init__",
          [](Direction* self, double up, double north, double east) {
            new (self) Direction{up, north, east};
          },
          "up"_a = 0.0, "north"_a = 0.0, "east"_a = 0.0)
      .def_rw("up", &Direction::up)
      .def_rw("north", &Direction::north)
      .def_rw("east", &Direction::east);
  m.def("along", &Along, "component"_a);
  m.def("radial_direction", &RadialDirection, "source"_a, "receiver"_a);
  m.def("transverse_direction", &TransverseDirection, "source"_a, "receiver"_a);
  nb::class_<Measurement>(m, "Measurement")
      .def(
          "__init__",
          [](Measurement* self, const Receiver& receiver,
             const Direction& direction, std::vector<double> adjointSource) {
            new (self)
                Measurement{receiver, direction, std::move(adjointSource)};
          },
          "receiver"_a, "direction"_a, "adjoint_source"_a)
      .def_rw("receiver", &Measurement::receiver)
      .def_rw("direction", &Measurement::direction)
      .def_prop_rw(
          "adjoint_source",
          [](const Measurement& s) { return Copy(s.adjointSource); },
          [](Measurement& s, std::vector<double> values) {
            s.adjointSource = std::move(values);
          },
          nb::rv_policy::move);
  m.def(
      "traveltime",
      [](const Model& model, const Calculation& calculation,
         const Source& source, const Receiver& receiver,
         const Direction& direction, const Window& window) {
        auto release = nb::gil_scoped_release{};
        return Traveltime(model, calculation, source, receiver, direction,
                          window);
      },
      "model"_a, "calculation"_a, "source"_a, "receiver"_a, "direction"_a,
      "window"_a);
  m.def(
      "traveltimes",
      [](const Model& model, const Calculation& calculation,
         const Source& source, const Receiver& receiver,
         const std::vector<Direction>& directions,
         const std::vector<Window>& windows) {
        auto release = nb::gil_scoped_release{};
        return Traveltimes(model, calculation, source, receiver, directions,
                           windows);
      },
      "model"_a, "calculation"_a, "source"_a, "receiver"_a, "directions"_a,
      "windows"_a);

  // Radial kernels.
  nb::enum_<RadialParameter>(m, "RadialParameter")
      .value("rho", RadialParameter::Rho)
      .value("A", RadialParameter::A)
      .value("C", RadialParameter::C)
      .value("F", RadialParameter::F)
      .value("L", RadialParameter::L)
      .value("N", RadialParameter::N)
      .value("qA", RadialParameter::QA)
      .value("qF", RadialParameter::QF)
      .value("qMu", RadialParameter::QMu);
  nb::class_<RadialMesh>(m, "RadialMesh")
      .def(
          "__init__",
          [](RadialMesh* self, int order, double maxElement) {
            new (self) RadialMesh{order, maxElement};
          },
          "order"_a = 5, "max_element"_a = 50e3)
      .def_rw("order", &RadialMesh::order)
      .def_rw("max_element", &RadialMesh::maxElement);
  nb::class_<RadialKernels>(m, "RadialKernels")
      .def_prop_ro(
          "radii", [](const RadialKernels& k) { return Copy(k.Radii()); },
          nb::rv_policy::move)
      .def_prop_ro(
          "weights", [](const RadialKernels& k) { return Copy(k.Weights()); },
          nb::rv_policy::move)
      .def_prop_ro("source_radius", &RadialKernels::SourceRadius)
      .def(
          "values",
          [](const RadialKernels& k, RadialParameter p) { return Copy(k[p]); },
          "parameter"_a)
      .def(
          "at_source",
          [](const RadialKernels& k, RadialParameter p) {
            return k.AtSource(p);
          },
          "parameter"_a);
  m.def(
      "compute_radial_kernels",
      [](const Model& model, const Calculation& calculation,
         const Source& source, const std::vector<Measurement>& measurements,
         const RadialMesh& mesh, nb::handle progress) {
        const auto reporter = Reporting(progress);
        auto release = nb::gil_scoped_release{};
        return ComputeRadialKernels(model, calculation, source, measurements,
                                    mesh, reporter);
      },
      "model"_a, "calculation"_a, "source"_a, "measurements"_a,
      "mesh"_a = RadialMesh{}, "progress"_a = nb::none());
  nb::enum_<ParameterSet>(m, "ParameterSet")
      .value("love", ParameterSet::Love)
      .value("moduli", ParameterSet::Moduli)
      .value("velocities", ParameterSet::Velocities)
      .value("anisotropic", ParameterSet::Anisotropic);
  nb::enum_<Scaling>(m, "Scaling")
      .value("absolute", Scaling::Absolute)
      .value("relative", Scaling::Relative);
  nb::class_<NamedKernels>(m, "NamedKernels")
      .def_ro("names", &NamedKernels::names)
      .def_prop_ro("values",
                   [](const NamedKernels& k) {
                     auto values = nb::list{};
                     for (const auto& v : k.values) values.append(Copy(v));
                     return values;
                   })
      .def_ro("at_source", &NamedKernels::atSource);
  m.def("reparameterise", &Reparameterise, "kernels"_a, "model"_a,
        "parameters"_a, "scaling"_a);

  // Kernels through the earth.
  auto sampling = nb::class_<KernelSampling>(m, "KernelSampling");
  sampling
      .def_prop_ro("node_count",
                   [](const KernelSampling& s) { return s.nodes.size(); })
      .def_prop_ro("place_count",
                   [](const KernelSampling& s) { return s.places.size(); });
  BindSamplingArrays(
      sampling, [](const KernelSampling& s) -> const auto& { return s.nodes; },
      [](const KernelSampling& s) -> const auto& { return s.places; });
  nb::class_<VolumeGrid>(m, "VolumeGrid")
      .def(
          "__init__",
          [](VolumeGrid* self, int degree, double maxElement) {
            new (self) VolumeGrid{degree, maxElement};
          },
          "degree"_a = 64, "max_element"_a = 100e3)
      .def_rw("degree", &VolumeGrid::degree)
      .def_rw("max_element", &VolumeGrid::maxElement);
  nb::class_<SectionOptions>(m, "SectionOptions")
      .def(
          "__init__",
          [](SectionOptions* self, double margin, double halfWidth,
             double spacing, double radialSpacing) {
            new (self)
                SectionOptions{margin, halfWidth, spacing, radialSpacing};
          },
          "margin"_a = 15.0, "half_width"_a = 15.0, "spacing"_a = 0.2,
          "radial_spacing"_a = 25e3)
      .def_rw("margin", &SectionOptions::margin)
      .def_rw("half_width", &SectionOptions::halfWidth)
      .def_rw("spacing", &SectionOptions::spacing)
      .def_rw("radial_spacing", &SectionOptions::radialSpacing);
  nb::class_<MapOptions>(m, "MapOptions")
      .def(
          "__init__",
          [](MapOptions* self, std::vector<double> depths, double margin,
             double halfWidth, double spacing) {
            new (self)
                MapOptions{std::move(depths), margin, halfWidth, spacing};
          },
          "depths"_a, "margin"_a = 10.0, "half_width"_a = 15.0,
          "spacing"_a = 0.5)
      .def_rw("depths", &MapOptions::depths)
      .def_rw("margin", &MapOptions::margin)
      .def_rw("half_width", &MapOptions::halfWidth)
      .def_rw("spacing", &MapOptions::spacing);
  m.def("volume_sampling", &VolumeSampling, "model"_a, "grid"_a = VolumeGrid{});
  m.def("path_section", &PathSection, "model"_a, "source"_a, "receiver"_a,
        "options"_a = SectionOptions{});
  m.def("cross_section", &CrossSection, "model"_a, "source"_a, "receiver"_a,
        "options"_a = SectionOptions{});
  m.def("path_maps", &PathMaps, "model"_a, "source"_a, "receiver"_a,
        "options"_a);
  auto sampled = nb::class_<SampledKernels>(m, "SampledKernels");
  sampled.def_prop_ro("names", &SampledKernels::Names)
      .def_prop_ro(
          "values",
          // (parameter, node, place)
          [](SampledKernels& k) {
            const auto shape = std::array{k.Names().size(), k.Nodes().size(),
                                          k.Places().size()};
            const auto empty =
                std::ranges::find(shape, std::size_t{0}) != shape.end();
            auto* data = empty ? nullptr : k.Values(0, 0).data();
            return Array<double, 3>(data, 3, shape.data());
          },
          nb::rv_policy::reference_internal);
  BindSamplingArrays(
      sampled, [](const SampledKernels& k) { return k.Nodes(); },
      [](const SampledKernels& k) { return k.Places(); });
  m.def(
      "compute_sampled_kernels",
      [](const Model& model, const Calculation& calculation,
         const Source& source, const std::vector<Measurement>& measurements,
         ParameterSet parameters, Scaling scaling,
         const std::vector<KernelSampling>& samplings, nb::handle progress) {
        const auto reporter = Reporting(progress);
        auto release = nb::gil_scoped_release{};
        return ComputeSampledKernels(model, calculation, source, measurements,
                                     parameters, scaling, samplings, reporter);
      },
      "model"_a, "calculation"_a, "source"_a, "measurements"_a, "parameters"_a,
      "scaling"_a, "samplings"_a, "progress"_a = nb::none());

  // The input files and the runs the programs make of them.
  nb::class_<SeismogramInput>(m, "SeismogramInput")
      .def_ro("model", &SeismogramInput::model)
      .def_ro("calculation", &SeismogramInput::calculation)
      .def_ro("sources", &SeismogramInput::sources)
      .def_ro("names", &SeismogramInput::names)
      .def_ro("receivers", &SeismogramInput::receivers)
      .def_ro("prefix", &SeismogramInput::prefix);
  m.def("read_seismogram_input", &ReadSeismogramInput, "path"_a);
  m.def(
      "run_seismograms",
      [](const SeismogramInput& input, nb::handle progress) {
        const auto reporter = Reporting(progress);
        auto release = nb::gil_scoped_release{};
        auto seismograms = ComputeSeismograms(input, reporter);
        WriteSeismograms(seismograms, input);
        return seismograms;
      },
      "input"_a, "progress"_a = nb::none());

  nb::class_<SnapshotInput>(m, "SnapshotInput")
      .def_ro("model", &SnapshotInput::model)
      .def_ro("calculation", &SnapshotInput::calculation)
      .def_ro("sources", &SnapshotInput::sources)
      .def_ro("names", &SnapshotInput::names)
      .def_ro("grid", &SnapshotInput::grid)
      .def_ro("times", &SnapshotInput::times)
      .def_ro("prefix", &SnapshotInput::prefix);
  m.def("read_snapshot_input", &ReadSnapshotInput, "path"_a);
  m.def(
      "run_snapshots",
      [](const SnapshotInput& input, nb::handle progress) {
        const auto reporter = Reporting(progress);
        auto release = nb::gil_scoped_release{};
        auto snapshots = ComputeSnapshots(input, reporter);
        WriteSnapshots(snapshots, input);
        return snapshots;
      },
      "input"_a, "progress"_a = nb::none());

  nb::class_<KernelMeasurement>(m, "KernelMeasurement")
      .def_ro("name", &KernelMeasurement::name)
      .def_ro("component", &KernelMeasurement::component)
      .def_ro("direction", &KernelMeasurement::direction)
      .def_ro("window", &KernelMeasurement::window);
  nb::class_<KernelInput>(m, "KernelInput")
      .def_ro("model", &KernelInput::model)
      .def_ro("calculation", &KernelInput::calculation)
      .def_ro("source", &KernelInput::source)
      .def_ro("name", &KernelInput::name)
      .def_ro("receiver", &KernelInput::receiver)
      .def_ro("measurements", &KernelInput::measurements)
      .def_ro("parameters", &KernelInput::parameters)
      .def_ro("scaling", &KernelInput::scaling)
      .def_ro("mesh", &KernelInput::mesh)
      .def_ro("volume", &KernelInput::volume)
      .def_ro("sections", &KernelInput::sections)
      .def_ro("maps", &KernelInput::maps)
      .def_ro("prefix", &KernelInput::prefix);
  m.def("read_kernel_input", &ReadKernelInput, "path"_a);
  nb::class_<MeasurementKernels>(m, "MeasurementKernels")
      .def_ro("name", &MeasurementKernels::name)
      .def_ro("measurement", &MeasurementKernels::measurement)
      .def_ro("kernels", &MeasurementKernels::kernels)
      .def_ro("named", &MeasurementKernels::named)
      .def_ro("volume", &MeasurementKernels::volume)
      .def_ro("path", &MeasurementKernels::path)
      .def_ro("across", &MeasurementKernels::across)
      .def_ro("maps", &MeasurementKernels::maps);
  nb::class_<KernelResult>(m, "KernelResult")
      .def_ro("seismogram", &KernelResult::seismogram)
      .def_ro("measurements", &KernelResult::measurements);
  m.def(
      "run_kernels",
      [](const KernelInput& input, nb::handle progress) {
        const auto reporter = StageReporting(progress);
        auto release = nb::gil_scoped_release{};
        auto result = ComputeKernels(input, reporter);
        WriteKernels(result, input);
        return result;
      },
      "input"_a, "progress"_a = nb::none());
}
