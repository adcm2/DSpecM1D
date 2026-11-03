"""Sources, receivers and their instruments."""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from . import _core
from .calculation import _QUANTITY, _choice

_MOMENT = ("rr", "rt", "rp", "tt", "tp", "pp")
_FORCE = ("r", "t", "p")


def _check_place(latitude, depth_km):
    if not -90 <= latitude <= 90:
        raise ValueError("latitude should lie between -90 and 90 degrees")
    if not depth_km >= 0:
        raise ValueError("depth_km should not be negative")


@dataclass
class Source:
    """A point source: a moment tensor, a force, or both.

    Parameters
    ----------
    latitude, longitude
        The epicentre (degrees; geographic latitude).
    depth_km
        The depth below the surface (km).
    moment_tensor
        The components rr, rt, rp, tt, tp and pp (N m), r up, t south and p
        east, as a mapping.
    force
        The components r, t and p (N), as a mapping; those left out are zero.
    name
        Its name, which labels its output.
    """

    latitude: float
    longitude: float
    depth_km: float
    moment_tensor: Mapping[str, float] | None = None
    force: Mapping[str, float] | None = None
    name: str = ""

    def _core(self):
        _check_place(self.latitude, self.depth_km)
        if self.moment_tensor is None and self.force is None:
            raise ValueError("a source needs a moment tensor, a force or both")
        moment = _core.MomentTensor()
        if self.moment_tensor is not None:
            unknown = set(self.moment_tensor) - set(_MOMENT)
            missing = set(_MOMENT) - set(self.moment_tensor)
            if unknown or missing:
                raise ValueError("moment_tensor should have exactly "
                                 f"{', '.join(_MOMENT)}")
            moment = _core.MomentTensor(**{k: float(v) for k, v in
                                           self.moment_tensor.items()})
        force = _core.Force()
        if self.force is not None:
            if set(self.force) - set(_FORCE):
                raise ValueError(f"force should have only {', '.join(_FORCE)}")
            force = _core.Force(**{k: float(v) for k, v in
                                   self.force.items()})
        return _core.Source(_core.Location(self.latitude, self.longitude),
                            self.depth_km * 1e3, moment, force)


@dataclass
class Instrument:
    """A seismometer, by its poles and zeros (rad/s) and constant: its
    response to the `input` motion at angular frequency w is constant times
    the product of (iw - zero) over the product of (iw - pole)."""

    zeros: Sequence[complex]
    poles: Sequence[complex]
    constant: float = 1.0
    input: str = "displacement"

    @classmethod
    def from_sac(cls, path):
        """The instrument of a SAC pole-zero file."""
        core = _core.read_sac_pole_zero(path)
        names = {v: k for k, v in _QUANTITY.items()}
        return cls(list(core.zeros), list(core.poles), core.constant,
                   names[core.input])

    def _core(self):
        return _core.Instrument(
            [complex(z) for z in self.zeros], [complex(p) for p in self.poles],
            float(self.constant), _choice("input", self.input, _QUANTITY))


@dataclass
class Receiver:
    """A receiver.

    Parameters
    ----------
    latitude, longitude
        Its place (degrees; geographic latitude).
    depth_km
        Its depth below the surface (km). It may be anywhere but the centre:
        in the ocean, the crust or the core.
    instrument
        A seismometer whose output is recorded, or None for the ground motion.
    name
        Its name, which labels its output.
    """

    latitude: float
    longitude: float
    depth_km: float = 0.0
    instrument: Instrument | None = None
    name: str = ""

    def _core(self):
        _check_place(self.latitude, self.depth_km)
        return _core.Receiver(
            _core.Location(self.latitude, self.longitude), self.depth_km * 1e3,
            None if self.instrument is None else self.instrument._core())


def _location(place):
    """A place as the library has it, from a Source, a Receiver or a pair of
    a latitude and a longitude."""
    if isinstance(place, (Source, Receiver)):
        return _core.Location(place.latitude, place.longitude)
    latitude, longitude = place
    return _core.Location(latitude, longitude)


def great_circle(start, azimuth_deg, distances_deg, depth_km=0.0,
                 instrument=None):
    """Receivers along a great circle: at the distances given (degrees) from
    `start`, a Source or a (latitude, longitude), in the direction of
    `azimuth_deg` there (clockwise from north)."""
    origin = _location(start)
    receivers = []
    for i, distance in enumerate(distances_deg):
        place = _core.great_circle(origin, float(azimuth_deg), float(distance),
                                   float(distance), 1)[0]
        receivers.append(Receiver(place.latitude, place.longitude, depth_km,
                                  instrument, name=f"{i + 1}"))
    return receivers


def _sources(sources):
    """One source or several, as a list."""
    return [sources] if isinstance(sources, Source) else list(sources)


def _receivers(receivers):
    return [receivers] if isinstance(receivers, Receiver) else list(receivers)

