"""The programs yspec, yspec-snapshot and yspec-kernels, as the package
installs them: each reads an input file, computes what it describes and
writes it, as the C++ programs do, reporting on the standard error."""

import sys
import time

from . import _core
from .progress import Reporter, StageReporter


def _plural(n, noun):
    return f"{n} {noun}{'' if n == 1 else 's'}"


def _summary(name, path, calculation, *parts):
    degrees, band = calculation.degrees, calculation.band
    print(f"{name}: model {path}, degrees {degrees.min} to {degrees.max}, "
          f"{band.min * 1e3:g} to {band.max * 1e3:g} mHz, "
          + ", ".join([*parts, _plural(_core.threads(), "thread")]),
          file=sys.stderr, flush=True)


def _main(name, usage, run):
    args = sys.argv[1:]
    if args in (["--help"], ["-h"]):
        print(usage, end="")
        return 0
    if args == ["--version"]:
        print(f"{name} {_core.__version__}")
        return 0
    if len(args) != 1:
        print(usage, end="", file=sys.stderr)
        return 2
    try:
        start = time.perf_counter()
        written = run(args[0])
        print(f"{name}: wrote {written} in {time.perf_counter() - start:.2f} s",
              file=sys.stderr)
    except KeyboardInterrupt:
        print(f"{name}: interrupted", file=sys.stderr)
        return 130
    except Exception as error:  # noqa: BLE001 - reported as the C++ does
        print(f"{name}: {error}", file=sys.stderr)
        return 1
    return 0


def _seismograms(path):
    data = _core.read_seismogram_input(path)
    _summary("yspec", data.model, data.calculation,
             _plural(len(data.sources), "source"),
             _plural(len(data.receivers), "receiver"))
    result = _core.run_seismograms(data, Reporter("yspec"))
    return (f"{_plural(result.sources * result.receivers, 'file')} with "
            f"prefix {data.prefix}")


def _snapshots(path):
    data = _core.read_snapshot_input(path)
    grid = data.grid
    _summary("yspec-snapshot", data.model, data.calculation,
             _plural(len(data.sources), "source"),
             _plural(len(data.times), "time"),
             f"{_plural(len(grid.depths), 'depth')} x "
             f"{_plural(len(grid.distances), 'distance')} x "
             f"{_plural(len(grid.azimuths), 'azimuth')}")
    result = _core.run_snapshots(data, Reporter("yspec-snapshot"))
    return (f"{_plural(result.sources * len(result.times), 'file')} with "
            f"prefix {data.prefix}")


def _kernels(path):
    data = _core.read_kernel_input(path)
    names = ", ".join(f"{m.name} on {m.component}" for m in data.measurements)
    _summary("yspec-kernels", data.model, data.calculation,
             f"{_plural(len(data.measurements), 'traveltime')}: {names}")
    result = _core.run_kernels(data, StageReporter("yspec-kernels"))
    return (f"the kernels of {_plural(len(result.measurements), 'traveltime')}"
            f" with prefix {data.prefix}")


def seismograms():
    sys.exit(_main("yspec", """usage: yspec input.json

Computes synthetic seismograms for the calculation described in input.json
and writes one file per receiver, <prefix>.1, <prefix>.2, ...
(<prefix>.<source>.1, ... for several sources), each line holding the time
(s) and the vertical, north and east components in SI units, or in those of
the receiver's instrument. The number of threads is set by OMP_NUM_THREADS.
""", _seismograms))


def snapshots():
    sys.exit(_main("yspec-snapshot", """usage: yspec-snapshot input.json

Computes the wavefield at the times asked for in input.json, on a grid of
depths, distances and azimuths around each source, and writes one file per
time, <prefix>.1, <prefix>.2, ... (<prefix>.<source>.1, ... for several
sources), each line holding a point of the grid, its latitude and longitude,
and the up, away and transverse components in SI units. The number of threads
is set by OMP_NUM_THREADS.
""", _snapshots))


def kernels():
    sys.exit(_main("yspec-kernels", """usage: yspec-kernels input.json

Computes the sensitivity kernels of cross-correlation traveltimes on the
seismogram of one source at one receiver, as input.json describes, by the
adjoint method: radial kernels, in <prefix>.kernels, and, if asked,
three-dimensional kernels on sections along and across the path
(<prefix>.sections) and on maps around it (<prefix>.maps), with the
seismogram, the window and the adjoint source in <prefix>.measurement. With
several traveltimes each name has its own after the prefix. The number of
threads is set by OMP_NUM_THREADS.
""", _kernels))
