"""Seismograms and strains at receivers."""

import numpy as np

from . import _core
from .places import _location, _receivers, _sources
from .progress import reporter


class Seismograms:
    """Seismograms: `data[source, receiver, component, sample]`, the
    components up (Z), north (N) and east (E), in SI units of the quantity
    computed, or in those of a receiver's instrument; `time_s` the time of
    each sample after the origin time."""

    components = ("Z", "N", "E")

    def __init__(self, core, sources, receivers, quantity):
        self._core = core
        self.sources = sources
        self.receivers = receivers
        self.quantity = quantity
        self.data = core.data
        self.time_step_s = core.time_step
        self.time_s = np.arange(core.samples) * core.time_step

    def __repr__(self):
        return (f"Seismograms({len(self.sources)} sources x "
                f"{len(self.receivers)} receivers x 3 components x "
                f"{len(self.time_s)} samples, {self.quantity})")

    def rotated(self):
        """The data with the horizontal components turned to radial (R, away
        from the source) and transverse (T): `[source, receiver, component,
        sample]` with the components Z, R and T."""
        out = np.empty_like(self.data)
        out[:, :, 0] = self.data[:, :, 0]
        for i, source in enumerate(self.sources):
            for j, receiver in enumerate(self.receivers):
                a, b = _location(source), _location(receiver)
                n, e = self.data[i, j, 1], self.data[i, j, 2]
                for k, d in ((1, _core.radial_direction(a, b)),
                             (2, _core.transverse_direction(a, b))):
                    out[i, j, k] = d.north * n + d.east * e
        return out

    def trace(self, component="Z", receiver=0, source=0):
        """One trace: of component Z, N, E, R or T, at a receiver from a
        source, by index."""
        if component in self.components:
            return self.data[source, receiver, self.components.index(component)]
        if component in ("R", "T"):
            a = _location(self.sources[source])
            b = _location(self.receivers[receiver])
            d = (_core.radial_direction(a, b) if component == "R" else
                 _core.transverse_direction(a, b))
            return (d.north * self.data[source, receiver, 1] +
                    d.east * self.data[source, receiver, 2])
        raise ValueError("component should be Z, N, E, R or T")

    def to_obspy(self, origin_time=None, rotate=False):
        """The seismograms as an ObsPy Stream, which needs ObsPy: a trace for
        each source, receiver and component, its station the receiver's name
        (or its number), its location the source's number and its channel
        LX followed by the component. `origin_time` is the source's (by
        default 1970-01-01). With `rotate`, the components are Z, R and T."""
        from obspy import Stream, Trace, UTCDateTime

        start = UTCDateTime(0) if origin_time is None else UTCDateTime(
            origin_time)
        data = self.rotated() if rotate else self.data
        names = ("Z", "R", "T") if rotate else self.components
        stream = Stream()
        for i in range(len(self.sources)):
            for j, receiver in enumerate(self.receivers):
                for k, component in enumerate(names):
                    trace = Trace(np.array(data[i, j, k]))
                    trace.stats.delta = self.time_step_s
                    trace.stats.starttime = start
                    trace.stats.network = "YS"
                    trace.stats.station = receiver.name or f"{j + 1}"
                    trace.stats.location = f"{i + 1:02d}"
                    trace.stats.channel = f"LX{component}"
                    stream.append(trace)
        return stream


class Strains:
    """Strains: `data[source, receiver, component, sample]`, the components
    rr, tt, pp, rt, rp and tp of the strain tensor, r up, t south and p east,
    and `time_s` the time of each sample after the origin time."""

    components = ("rr", "tt", "pp", "rt", "rp", "tp")

    def __init__(self, core, sources, receivers):
        self._core = core
        self.sources = sources
        self.receivers = receivers
        self.data = core.data
        self.time_step_s = core.time_step
        self.time_s = np.arange(core.samples) * core.time_step


def seismograms(model, calculation, sources, receivers, progress=None):
    """Synthetic seismograms at receivers from sources.

    Parameters
    ----------
    model
        A Model.
    calculation
        A Calculation.
    sources
        A Source, or a sequence of them.
    receivers
        A Receiver, or a sequence of them.
    progress
        None, True to print progress, or a callable taking (done, total).
    """
    sources, receivers = _sources(sources), _receivers(receivers)
    core = _core.compute_seismograms(
        model._model, calculation._core(), [s._core() for s in sources],
        [r._core() for r in receivers], reporter(progress, "yspec"))
    return Seismograms(core, sources, receivers, calculation.quantity)


def strains(model, calculation, sources, receivers, progress=None):
    """The strain tensors at receivers from sources, as `seismograms` takes
    them."""
    sources, receivers = _sources(sources), _receivers(receivers)
    core = _core.compute_strains(
        model._model, calculation._core(), [s._core() for s in sources],
        [r._core() for r in receivers], reporter(progress, "yspec"))
    return Strains(core, sources, receivers)
