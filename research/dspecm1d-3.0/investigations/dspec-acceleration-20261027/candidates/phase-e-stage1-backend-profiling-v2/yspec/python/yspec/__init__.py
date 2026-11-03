"""YSpec: synthetic seismograms, snapshots of the wavefield and adjoint
sensitivity kernels for spherically symmetric earth models, by direct radial
integration of the equations of motion in the frequency domain.

    import numpy as np
    import yspec

    model = yspec.Model.prem()
    source = yspec.Source(latitude=36.87, longitude=35.58, depth_km=29.5,
                          moment_tensor=dict(rr=4.22e17, rt=3.08e17,
                                             rp=6.99e17, tt=-2.976e18,
                                             tp=6.54e17, pp=2.55e18))
    receivers = yspec.great_circle(source, azimuth_deg=60,
                                   distances_deg=np.arange(5, 180, 5))
    calculation = yspec.Calculation(degrees=300, frequencies_mhz=(0.2, 15),
                                    length_min=60, filter_mhz=(1, 2, 12, 15),
                                    quantity="velocity")
    seismograms = yspec.seismograms(model, calculation, source, receivers,
                                    progress=True)
    seismograms.data  # [source, receiver, component, sample], Z, N and E

Lengths are in km, angles in degrees, frequencies in mHz and times in s (the
length of the seismograms in minutes), as in the input files of the programs,
which the package installs too, with yspec-plot to draw what they write.
"""

from . import _core
from .calculation import Calculation
from .kernels import (Measurement, RadialKernels, SampledKernels, Sampling,
                      Traveltime, cross_section, path_maps, path_section,
                      radial_kernels, sampled_kernels, traveltime, traveltimes,
                      volume)
from .model import Model
from .places import Instrument, Receiver, Source, great_circle
from .seismograms import Seismograms, Strains, seismograms, strains
from .snapshots import Snapshots, snapshots

__version__ = _core.__version__


def threads():
    """The number of threads the calculations use."""
    return _core.threads()


def set_threads(count):
    """Sets the number of threads the calculations use (by default, as
    OMP_NUM_THREADS says, or every core)."""
    _core.set_threads(int(count))


__all__ = [
    "Calculation", "Model", "Source", "Receiver", "Instrument",
    "great_circle", "seismograms", "strains", "Seismograms", "Strains",
    "snapshots", "Snapshots", "Traveltime", "Measurement", "traveltime",
    "traveltimes", "radial_kernels", "RadialKernels", "Sampling", "volume",
    "path_section", "cross_section", "path_maps", "sampled_kernels",
    "SampledKernels", "threads", "set_threads",
]
