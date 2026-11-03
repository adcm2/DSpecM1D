"""Snapshots of the wavefield on a grid around each source."""

import numpy as np

from . import _core
from .places import _sources
from .progress import reporter


class Snapshots:
    """The wavefield at chosen times: `data[source, time, depth, azimuth,
    distance, component]`, the components up, away from the source along the
    great circle, and transverse to it (90 degrees anticlockwise from away,
    seen from above), in SI units of the quantity computed."""

    components = ("up", "away", "transverse")

    def __init__(self, core, sources, times_s, depths_km, distances_deg,
                 azimuths_deg, quantity):
        self._core = core
        self.sources = sources
        self.times_s = times_s
        self.depths_km = depths_km
        self.distances_deg = distances_deg
        self.azimuths_deg = azimuths_deg
        self.quantity = quantity
        self.data = core.data

    def __repr__(self):
        shape = " x ".join(str(n) for n in self.data.shape)
        return f"Snapshots({shape}, {self.quantity})"


def snapshots(model, calculation, sources, times_s, depths_km, distances_deg,
              azimuths_deg, progress=None):
    """The wavefield at times (s) on a grid in the frame of each source: at
    depths (km), distances from it (degrees) and azimuths at it (degrees,
    clockwise from north). Two azimuths 180 degrees apart make a
    cross-section through the source; one depth and many azimuths, a map
    around it. The points may lie anywhere but at the centre.

    The time series is computed only as far as `calculation.length_min`, so
    the times should lie within it."""
    sources = _sources(sources)
    arrays = [np.atleast_1d(np.asarray(v, float))
              for v in (times_s, depths_km, distances_deg, azimuths_deg)]
    times, depths, distances, azimuths = arrays
    grid = _core.SnapshotGrid(list(depths * 1e3), list(distances),
                              list(azimuths))
    core = _core.compute_snapshots(
        model._model, calculation._core(), [s._core() for s in sources], grid,
        list(times), reporter(progress, "yspec-snapshot"))
    return Snapshots(core, sources, times, depths, distances, azimuths,
                     calculation.quantity)
