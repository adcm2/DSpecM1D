"""The output of a run, found from its input file."""

import json
import math
import pathlib
import re
from dataclasses import dataclass

import numpy as np
from planetmodel.units import Dimensions

#: The dimensions of each quantity a run may record.
DIMENSIONS = {
    "displacement": Dimensions(length=1),
    "velocity": Dimensions(length=1, time=-1),
    "acceleration": Dimensions(length=1, time=-2),
}

#: The SI unit of each, for labels.
UNITS = {"displacement": "m", "velocity": "m/s", "acceleration": "m/s$^2$"}


def values(spec):
    """The values a number, an array or a range {min, max, count} gives, as
    the programs read them."""
    if isinstance(spec, (int, float)):
        return [float(spec)]
    if isinstance(spec, list):
        return [float(v) for v in spec]
    count = spec["count"]
    if count == 1:
        return [float(spec["min"])]
    return list(np.linspace(spec["min"], spec["max"], count))


def source_names(sources):
    """The names of the sources, numbered from 1 where not given."""
    return [s.get("name", str(i + 1)) for i, s in enumerate(sources)]


def _choose(names, source, path):
    name = names[0] if source is None else source
    if name not in names:
        raise ValueError(f"there is no source {name} in {path}")
    return name


@dataclass(frozen=True)
class Receiver:
    """A receiver as the table beside the seismograms gives it, as seen
    from one source."""

    number: int
    file: pathlib.Path
    latitude: float
    longitude: float
    depth: float          # km
    distance: float       # degrees
    azimuth: float        # degrees, at the source
    back_azimuth: float   # degrees, of the source at the receiver
    units: str

    def trace(self, component):
        """The times (s) and one component: Z (up), N, E, or R and T,
        radial (away from the source) and transverse (90 degrees
        anticlockwise from radial, seen from above, as in the snapshots).
        """
        t, z, n, e = np.loadtxt(self.file).T
        back = math.radians(self.back_azimuth)
        components = {
            "Z": lambda: z,
            "N": lambda: n,
            "E": lambda: e,
            "R": lambda: -n * math.cos(back) - e * math.sin(back),
            "T": lambda: -n * math.sin(back) + e * math.cos(back),
        }
        if component not in components:
            raise ValueError(f"unknown component {component}: use Z, N, E, R "
                             "or T")
        return t, components[component]()


class SeismogramRun:
    """A run of yspec: its input, and the receivers of one source."""

    def __init__(self, path, source=None):
        self.path = pathlib.Path(path)
        self.input = json.loads(self.path.read_text())
        output = self.input["output"]
        prefix = self.path.parent / output["prefix"]
        self.quantity = output.get("quantity", "displacement")
        self.source = _choose(source_names(self.input["sources"]), source,
                              path)
        table = prefix.with_name(prefix.name + ".receivers")
        rows = [line.split() for line in table.read_text().splitlines()
                if not line.startswith("#")]
        rows = [row for row in rows if row[1] == self.source]
        self.receivers = []
        for number, (row, instrument_unit) in enumerate(
                zip(rows, self._instruments()), start=1):
            file, _, latitude, longitude, depth, distance, azimuth, back = row
            self.receivers.append(Receiver(
                number=number, file=prefix.parent / file,
                latitude=float(latitude), longitude=float(longitude),
                depth=float(depth),
                distance=float(distance), azimuth=float(azimuth),
                back_azimuth=float(back),
                units=(instrument_unit or "instrument output (unit unspecified)"
                       if instrument_unit is not False
                       else UNITS[self.quantity])))

    def receiver(self, number):
        """A receiver, numbered from 1 as the files are."""
        if not 1 <= number <= len(self.receivers):
            raise ValueError(f"there is no receiver {number} in {self.path}")
        return self.receivers[number - 1]

    def _instruments(self):
        """Whether each receiver has an instrument, group by group, depth by
        depth, place by place, as the program numbers them."""
        groups = self.input["receivers"]
        if isinstance(groups, dict):
            groups = [groups]
        for group in groups:
            depths = group["depth_km"]
            depths = len(depths) if isinstance(depths, list) else 1
            places = (len(group["locations"]) if "locations" in group
                      else group["profile"]["count"])
            output_unit = False
            if "instrument" in group:
                instrument = self.path.parent / group["instrument"]
                try:
                    text = instrument.read_text(errors="replace")
                except OSError:
                    text = ""
                match = re.search(r"^\s*\*\s*OUTPUT\s+UNIT\s*:\s*(\S.*?)\s*$",
                                  text, flags=re.IGNORECASE | re.MULTILINE)
                output_unit = match.group(1).strip() if match else None
            for _ in range(depths * places):
                yield output_unit


class SnapshotRun:
    """A run of yspec-snapshot: its grid and the fields at each time, for
    one source, as arrays (depth, azimuth, distance)."""

    COLUMNS = {"up": 5, "away": 6, "transverse": 7}

    def __init__(self, path, source=None):
        self.path = pathlib.Path(path)
        self.input = json.loads(self.path.read_text())
        prefix = self.path.parent / self.input["output"]["prefix"]
        self.quantity = self.input["output"].get("quantity", "displacement")
        names = source_names(self.input["sources"])
        self.source = _choose(names, source, path)
        self.hypocentre = self.input["sources"][names.index(self.source)]
        self.model = self.path.parent / self.input["model"]
        self.times = values(self.input["snapshots"]["times_s"])
        self._data = []
        for k in range(1, len(self.times) + 1):
            suffix = f".{k}" if len(names) == 1 else f".{self.source}.{k}"
            self._data.append(np.loadtxt(prefix.with_name(prefix.name +
                                                          suffix)))
        first = self._data[0]
        self.depths = np.unique(first[:, 0])
        self.distances = np.unique(first[:, 1])
        self.azimuths = np.array(list(dict.fromkeys(first[:, 2])))
        shape = (len(self.depths), len(self.azimuths), len(self.distances))
        coordinates = first[:, 3:5].reshape(shape + (2,))[0]
        self.latitude = coordinates[..., 0]
        self.longitude = coordinates[..., 1]

    def field(self, time_index, component):
        """One component at one time, (depth, azimuth, distance)."""
        shape = (len(self.depths), len(self.azimuths), len(self.distances))
        return self._data[time_index][:, self.COLUMNS[component]].reshape(shape)


def _along(start, azimuth, distance):
    """The geographic latitude and longitude (degrees) at a distance and
    azimuth (degrees) from a point, along a great circle of the sphere of
    geocentric latitudes, as YSpec places receivers."""
    ratio = 0.993277
    latitude = math.atan(ratio * math.tan(math.radians(start[0])))
    longitude = math.radians(start[1])
    a = math.radians(azimuth)
    d = math.radians(distance)
    end = math.asin(math.sin(latitude) * math.cos(d) +
                    math.cos(latitude) * math.sin(d) * math.cos(a))
    turn = math.atan2(math.sin(a) * math.sin(d) * math.cos(latitude),
                      math.cos(d) - math.sin(latitude) * math.sin(end))
    return (math.degrees(math.atan(math.tan(end) / ratio)),
            math.degrees(longitude + turn))


def _distance(a, b):
    """The angle (degrees) between two geographic places, on the sphere of
    their geocentric latitudes."""
    def unit(place):
        latitude = math.atan(0.993277 * math.tan(math.radians(place[0])))
        longitude = math.radians(place[1])
        return np.array([math.cos(latitude) * math.cos(longitude),
                         math.cos(latitude) * math.sin(longitude),
                         math.sin(latitude)])
    x, y = unit(a), unit(b)
    return math.degrees(math.atan2(np.linalg.norm(np.cross(x, y)),
                                   np.dot(x, y)))


@dataclass(frozen=True)
class SectionLayer:
    """One layer of a section or a stack of maps: depths (km) by places."""

    index: int
    depth: np.ndarray         # (n,), km
    values: dict              # name -> (n, m)


@dataclass(frozen=True)
class Sampled:
    """Kernels sampled on a section or on maps: the places' distances along
    and across the path (degrees), their latitudes and longitudes, and the
    kernels layer by layer."""

    along: np.ndarray
    across: np.ndarray
    latitude: np.ndarray
    longitude: np.ndarray
    layers: list


def _read_sampled(file, section=None):
    """The kernels in a .sections or .maps file, of one section if named."""
    header = [line for line in file.read_text().splitlines()
              if line.startswith("#")]
    columns = header[-1].lstrip("# ").split()
    rows = [line.split() for line in file.read_text().splitlines()
            if not line.startswith("#")]
    if section is not None:
        rows = [row[1:] for row in rows if row[0] == section]
        columns = columns[1:]
    table = np.array(rows, dtype=float)
    names = columns[6:]
    layer, depth = table[:, 0].astype(int), table[:, 1]
    # The places of the first node, which every node shares.
    first = (layer == layer[0]) & (depth == depth[0])
    places = int(first.sum())
    along, across, latitude, longitude = table[:places, 2:6].T
    layers = []
    for index in dict.fromkeys(layer):
        rows_here = table[layer == index]
        nodes = len(rows_here) // places
        depth = rows_here[::places, 1]
        values = {name: rows_here[:, 6 + j].reshape(nodes, places)
                  for j, name in enumerate(names)}
        layers.append(SectionLayer(int(index), depth, values))
    return Sampled(along, across, latitude, longitude, layers)


class KernelRun:
    """A run of yspec-kernels: one of its measurements, and its kernels."""

    def __init__(self, path, measurement=None):
        self.path = pathlib.Path(path)
        self.input = json.loads(self.path.read_text())
        output = self.input["output"]
        prefix = self.path.parent / output["prefix"]
        self.quantity = output.get("quantity", "displacement")
        self.model = self.path.parent / self.input["model"]
        self.hypocentre = self.input["sources"][0]
        measurements = self.input["measurements"]
        names = [m.get("name", str(i + 1))
                 for i, m in enumerate(measurements)]
        self.names = names
        name = names[0] if measurement is None else measurement
        if name not in names:
            raise ValueError(f"there is no measurement {name} in {path}: "
                             f"there are {', '.join(names)}")
        self.name = name
        chosen = measurements[names.index(name)]
        self.component = chosen["component"]
        self.window = chosen["traveltime"]["window_s"]
        if len(names) > 1:
            prefix = prefix.with_name(f"{prefix.name}.{name}")
        self.prefix = prefix
        kernels = self.input.get("kernels", {})
        self.relative = kernels.get("scaling", "absolute") == "relative"
        self.time, self.seismogram, self.taper, self.adjoint = np.loadtxt(
            self._file(".measurement")).T
        file = self._file(".kernels")
        header = [line for line in file.read_text().splitlines()
                  if line.startswith("#")]
        self.description = header[0].removeprefix("# yspec-kernels: ")
        self.at_source = [float(v) for v in next(
            line for line in header
            if line.startswith("# at_source")).split()[2:]]
        self.parameters = header[-1].split()[4:]
        table = np.loadtxt(file)
        self.radius, self.depth, self.weight = table[:, :3].T
        self.kernels = {p: table[:, 3 + j]
                        for j, p in enumerate(self.parameters)}
        self.at_source = dict(zip(self.parameters, self.at_source))
        self.surface = self.radius[-1] + self.depth[-1]
        self.source_radius = self.surface - self.hypocentre["depth_km"]
        self.volume = self._file(".nc")
        # Where the receiver is: at a location, or at a distance and azimuth
        # from the source.
        receiver = self.input["receiver"]
        self.receiver_depth = receiver["depth_km"]
        source = (self.hypocentre["latitude"], self.hypocentre["longitude"])
        if "location" in receiver:
            self.receiver_location = tuple(receiver["location"])
        else:
            self.receiver_location = _along(source, receiver["azimuth_deg"],
                                            receiver["distance_deg"])
        self.distance = _distance(source, self.receiver_location)

    def _file(self, suffix):
        return self.prefix.with_name(self.prefix.name + suffix)

    def sections(self):
        """The kernels on the sections along and across the path."""
        file = self._file(".sections")
        return (_read_sampled(file, "path"), _read_sampled(file, "across"))

    def maps(self):
        """The kernels on the maps."""
        return _read_sampled(self._file(".maps"))
