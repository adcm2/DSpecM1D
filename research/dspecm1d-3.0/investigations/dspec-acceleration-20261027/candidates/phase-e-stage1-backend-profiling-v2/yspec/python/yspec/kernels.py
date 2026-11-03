"""Traveltimes and their sensitivity kernels, by the adjoint method.

A traveltime is measured by cross-correlation of a synthetic seismogram with
itself in a window, and is given to the kernels by its adjoint source. Its
radial kernels give its sensitivity to spherically symmetric perturbations of
the model; its three-dimensional kernels, to laterally varying perturbations
of the elastic parameters and quality factors at fixed density, sampled
through the earth, on sections about the path or on maps.
"""

from dataclasses import dataclass, replace

import numpy as np

from . import _core
from .calculation import _choice
from .places import Receiver, Source, _location
from .progress import reporter

_PARAMETERS = {"love": _core.ParameterSet.love,
               "moduli": _core.ParameterSet.moduli,
               "velocities": _core.ParameterSet.velocities,
               "anisotropic": _core.ParameterSet.anisotropic}
_SCALING = {"absolute": _core.Scaling.absolute,
            "relative": _core.Scaling.relative}
_GEOCENTRIC_RATIO = 0.993277


def _direction(component, source, receiver):
    if component == "Z":
        return _core.along(_core.Component.vertical)
    if component == "N":
        return _core.along(_core.Component.north)
    if component == "E":
        return _core.along(_core.Component.east)
    if component == "R":
        return _core.radial_direction(_location(source), _location(receiver))
    if component == "T":
        return _core.transverse_direction(_location(source),
                                          _location(receiver))
    raise ValueError("component should be Z, N, E, R or T")


@dataclass
class Traveltime:
    """A traveltime to measure: of the seismogram of a component (Z, N, E, or
    R and T, radial and transverse to the source), by cross-correlation in a
    window (s after the origin time) with cosine tapers of `taper_s` inside
    it. A window of a period or two about one phase keeps its kernel
    simple."""

    component: str
    window_s: tuple[float, float]
    taper_s: float = 0.0
    name: str = ""


class Measurement:
    """A traveltime measured on the synthetic: what the kernels need of it,
    its adjoint source, sampled at `time_s`."""

    def __init__(self, core, traveltime, source, receiver, time_s):
        self._core = core
        self.traveltime = traveltime
        self.name = traveltime.name
        self.source = source
        self.receiver = receiver
        self.time_s = time_s

    @property
    def adjoint_source(self):
        return self._core.adjoint_source

    def __repr__(self):
        start, end = self.traveltime.window_s
        return (f"Measurement({self.name or 'traveltime'} on "
                f"{self.traveltime.component}, {start:g} to {end:g} s)")


def traveltimes(model, calculation, source, receiver, measured):
    """Measures traveltimes at a receiver, from one synthetic: `measured` a
    Traveltime, or several, giving a Measurement, or a list of them."""
    single = isinstance(measured, Traveltime)
    measured = [measured] if single else list(measured)
    if not isinstance(source, Source) or not isinstance(receiver, Receiver):
        raise TypeError("traveltimes need one Source and one Receiver")
    measured = [t if t.name else replace(t, name=f"{i + 1}")
                for i, t in enumerate(measured)]
    directions = [_direction(t.component, source, receiver) for t in measured]
    windows = [_core.Window(*t.window_s, t.taper_s) for t in measured]
    calc = calculation._core()
    cores = _core.traveltimes(model._model, calc, source._core(),
                              receiver._core(), directions, windows)
    measurements = [
        Measurement(c, t, source, receiver,
                    np.arange(len(c.adjoint_source)) * calc.sampling.time_step)
        for c, t in zip(cores, measured)]
    return measurements[0] if single else measurements


def traveltime(model, calculation, source, receiver, component, window_s,
               taper_s=0.0):
    """Measures one traveltime: see Traveltime."""
    return traveltimes(model, calculation, source, receiver,
                       Traveltime(component, tuple(window_s), taper_s))


def _measurements(measurements):
    single = isinstance(measurements, Measurement)
    measurements = [measurements] if single else list(measurements)
    if not measurements:
        raise ValueError("there should be at least one measurement")
    source = measurements[0].source
    return single, measurements, source


class RadialKernels:
    """The radial kernels of a measurement, at the nodes of a mesh in
    radius: `kernels[name]` per unit volume (s m^-3, or per unit volume of
    the logarithm when relative), `weight` the weight of each node in the
    integral over the earth (m^3), so that the change in the traveltime is
    the sum over nodes of weight times kernel times perturbation, plus the
    point parts `at_source[name]` times the perturbations at the source."""

    def __init__(self, core, named, model, parameters, scaling, measurement):
        self.measurement = measurement
        self.parameters = parameters
        self.scaling = scaling
        self.radius_km = core.radii / 1e3
        self.depth_km = model.surface_radius_km - self.radius_km
        self.weight = core.weights
        self.source_radius_km = core.source_radius / 1e3
        self.names = list(named.names)
        self.kernels = dict(zip(self.names, named.values))
        self.at_source = dict(zip(self.names, named.at_source))

    def __getitem__(self, name):
        return self.kernels[name]

    def __repr__(self):
        return (f"RadialKernels({', '.join(self.names)} at "
                f"{len(self.radius_km)} nodes, {self.scaling})")


def radial_kernels(model, calculation, measurements, parameters="velocities",
                   scaling="absolute", mesh_order=5, max_element_km=50.0,
                   progress=None):
    """The radial kernels of measurements made at one receiver, which share
    the forward field and the adjoint fields: RadialKernels for one
    Measurement, a list for several.

    Parameters
    ----------
    parameters
        "love" (density, A, C, F, L, N and their inverse quality factors),
        "moduli" (density, the bulk and shear moduli and their inverse quality
        factors), "velocities" (density at fixed velocities, the isotropic
        speeds and the inverse quality factors) or "anisotropic" (density at
        fixed velocities, the five speeds and eta of transverse isotropy, and
        the inverse quality factors).
    scaling
        "absolute", for perturbations of the parameters, or "relative", for
        perturbations of their logarithms.
    mesh_order, max_element_km
        Gauss-Legendre nodes in elements between the discontinuities, the
        source and the receiver, no larger than asked.
    """
    single, measurements, source = _measurements(measurements)
    chosen = _choice("parameters", parameters, _PARAMETERS)
    scaled = _choice("scaling", scaling, _SCALING)
    cores = _core.compute_radial_kernels(
        model._model, calculation._core(), source._core(),
        [m._core for m in measurements],
        _core.RadialMesh(int(mesh_order), max_element_km * 1e3),
        reporter(progress, "yspec-kernels: radial kernels"))
    kernels = [RadialKernels(c, _core.reparameterise(c, model._model, chosen,
                                                     scaled),
                             model, parameters, scaling, m)
               for c, m in zip(cores, measurements)]
    return kernels[0] if single else kernels


class Sampling:
    """Where three-dimensional kernels are evaluated: nodes in radius at
    places on the sphere."""

    def __init__(self, kind, core):
        self.kind = kind
        self._core = core

    def __repr__(self):
        return (f"Sampling({self.kind}: {self._core.node_count} nodes x "
                f"{self._core.place_count} places)")


def volume(model, degree=64, max_element_km=100.0):
    """The whole earth: at nodes in each layer no further apart than
    `max_element_km`, both sides of each discontinuity included, on a
    Gauss-Legendre grid resolving degree `degree` and at the poles."""
    return Sampling("volume", _core.volume_sampling(
        model._model, _core.VolumeGrid(int(degree), max_element_km * 1e3)))


def _section_options(margin_deg, half_width_deg, spacing_deg,
                     radial_spacing_km):
    return _core.SectionOptions(margin_deg, half_width_deg, spacing_deg,
                                radial_spacing_km * 1e3)


def path_section(model, source, receiver, margin_deg=15.0,
                 half_width_deg=15.0, spacing_deg=0.2,
                 radial_spacing_km=25.0):
    """The section along the great circle through the source and the
    receiver, from `margin_deg` before the source to as far beyond the
    receiver, at places `spacing_deg` apart and nodes in each layer no
    further apart than `radial_spacing_km`."""
    return Sampling("path", _core.path_section(
        model._model, _location(source), _location(receiver),
        _section_options(margin_deg, half_width_deg, spacing_deg,
                         radial_spacing_km)))


def cross_section(model, source, receiver, margin_deg=15.0,
                  half_width_deg=15.0, spacing_deg=0.2,
                  radial_spacing_km=25.0):
    """The section across the path at its middle, `half_width_deg` to each
    side, spaced as path_section's."""
    return Sampling("across", _core.cross_section(
        model._model, _location(source), _location(receiver),
        _section_options(margin_deg, half_width_deg, spacing_deg,
                         radial_spacing_km)))


def path_maps(model, source, receiver, depths_km, margin_deg=10.0,
              half_width_deg=15.0, spacing_deg=0.5):
    """Maps at depths (km) about the path, from `margin_deg` before the
    source to as far beyond the receiver and `half_width_deg` to each side,
    at places `spacing_deg` apart. A map at a discontinuity is on its upper
    side."""
    options = _core.MapOptions([d * 1e3 for d in depths_km], margin_deg,
                               half_width_deg, spacing_deg)
    return Sampling("maps", _core.path_maps(
        model._model, _location(source), _location(receiver), options))


class SampledKernels:
    """Three-dimensional kernels on a sampling: `values[parameter, node,
    place]`, per unit volume (s m^-3, or of the logarithm when relative), with
    each node's layer, radius and depth and each place's latitude and
    longitude (degrees), and its distances along and across the path
    (degrees) on sections and maps."""

    def __init__(self, core, model, measurement, sampling):
        self._core = core
        self.measurement = measurement
        self.sampling = sampling
        self.names = list(core.names)
        self.values = core.values
        self.layers = core.layers
        self.radius_km = core.radii / 1e3
        self.depth_km = model.surface_radius_km - self.radius_km
        colatitudes = core.colatitudes
        self.latitude = np.degrees(np.arctan(
            np.tan(np.pi / 2 - colatitudes) / _GEOCENTRIC_RATIO))
        self.longitude = np.degrees(core.longitudes)
        self.along_deg = core.along
        self.across_deg = core.across

    def __getitem__(self, name):
        """The kernel of a parameter, `[node, place]`."""
        return self.values[self.names.index(name)]

    def __repr__(self):
        return (f"SampledKernels({', '.join(self.names)} on "
                f"{self.sampling!r})")


def sampled_kernels(model, calculation, measurements, samplings,
                    parameters="velocities", scaling="absolute",
                    progress=None):
    """The three-dimensional kernels of measurements made at one receiver,
    at fixed density, on samplings (volume, path_section, cross_section,
    path_maps). For one Measurement and one Sampling, SampledKernels; for
    several of either, lists of them, by measurement and then by sampling.
    The parameters and scaling are as radial_kernels takes them, density
    left out."""
    single, measurements, source = _measurements(measurements)
    one_sampling = isinstance(samplings, Sampling)
    samplings = [samplings] if one_sampling else list(samplings)
    cores = _core.compute_sampled_kernels(
        model._model, calculation._core(), source._core(),
        [m._core for m in measurements],
        _choice("parameters", parameters, _PARAMETERS),
        _choice("scaling", scaling, _SCALING), [s._core for s in samplings],
        reporter(progress, "yspec-kernels: three-dimensional kernels"))
    result = [[SampledKernels(c, model, m, s) for c, s in zip(row, samplings)]
              for row, m in zip(cores, measurements)]
    if one_sampling:
        result = [row[0] for row in result]
    return result[0] if single else result

