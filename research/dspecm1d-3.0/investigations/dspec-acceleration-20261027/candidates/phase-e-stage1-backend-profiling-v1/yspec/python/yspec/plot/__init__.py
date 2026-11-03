"""Pictures of what the YSpec programs write.

Each picture starts from the input file of a run, from which it finds the
output: for yspec, the seismogram files and the table of receivers beside
them (`SeismogramRun`); for yspec-snapshot, the file of each time
(`SnapshotRun`); for yspec-kernels, the kernels and the measurement
(`KernelRun`). The functions draw on figures they make and return, and
the command `yspec-plot` saves or shows them.
"""

from .kernels import radial_kernels
from .runs import KernelRun, SeismogramRun, SnapshotRun
from .seismograms import compare, record_section, traces
from .snapshots import snapshots

__all__ = ["SeismogramRun", "SnapshotRun", "KernelRun", "record_section",
           "traces", "compare", "snapshots", "radial_kernels"]
