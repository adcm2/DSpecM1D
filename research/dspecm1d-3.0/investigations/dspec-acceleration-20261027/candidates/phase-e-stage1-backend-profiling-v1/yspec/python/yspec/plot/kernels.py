"""Pictures of the radial kernels of a measurement."""

import numpy as np
from planetmodel.deck import MINEOS, read_deck
from planetmodel.plotting import BOUNDARY_STYLE, profile
from planetmodel.slices import ProfileSlice, Shown, SliceLayer

from .runs import UNITS
from .seismograms import COLOURS
from .snapshots import GUIDE, STAR, _phase_boundaries


def _parameter_unit(parameter):
    # Keep this list aligned with Kernels.cpp::Names(ParameterSet): these are
    # the parameter names the solver can actually write to a kernel file.
    units = {
        "rho": "kg m$^{-3}$",
        "vp": "m s$^{-1}$", "vs": "m s$^{-1}$",
        "vpv": "m s$^{-1}$", "vph": "m s$^{-1}$",
        "vsv": "m s$^{-1}$", "vsh": "m s$^{-1}$",
        "A": "Pa", "C": "Pa", "F": "Pa", "L": "Pa", "N": "Pa",
        "kappa": "Pa", "mu": "Pa",
        "eta": "dimensionless",
        "qA": "dimensionless", "qF": "dimensionless",
        "qkappa": "dimensionless", "qmu": "dimensionless",
    }
    try:
        return units[parameter]
    except KeyError as error:
        raise ValueError(f"unsupported kernel parameter: {parameter}") from error


def _kernel_unit(parameter, relative):
    """The units of a three-dimensional traveltime kernel."""
    unit = "s m$^{-3}$"
    if relative:
        return unit
    return f"{unit} per ({_parameter_unit(parameter)})"


def _profile(radius, values, boundaries, name, lowest):
    """The saved kernel nodes, split into separate profile layers at boundaries.

    Do not extrapolate near-boundary nodes to an interface: the output does
    not provide one-sided interface values, and extrapolating can imply them.
    """
    layers = []
    which = np.searchsorted(boundaries, radius)
    for index in np.unique(which):
        inside = which == index
        r = radius[inside]
        v = values[inside]
        points = np.stack([np.zeros_like(r), np.zeros_like(r), r], axis=-1)
        layers.append(SliceLayer(index=int(index), radius=r, values=v,
                                 reference=points, physical=points))
    return ProfileSlice(theta=None, phi=None, layers=tuple(layers),
                        shown=Shown(name))


def radial_kernels(run, names=None, *, deepest=None):
    """The measured seismogram and its window above, and the kernel of each
    parameter (by default all of them) against radius below, each with its
    point part at the source marked, down to `deepest` km if given."""
    import matplotlib.pyplot as plt

    names = names or run.parameters
    deck = read_deck(run.model, MINEOS)
    boundaries = deck.boundaries / 1e3
    phase, surface = _phase_boundaries(run.model)
    fig = plt.figure(figsize=(2.6 * len(names) + 1.5, 9))
    grid = fig.add_gridspec(2, len(names), height_ratios=[1, 3.2])
    top = fig.add_subplot(grid[0, :])
    t = run.time / 60
    peak = np.abs(run.seismogram).max() or 1.0
    top.plot(t, run.seismogram, lw=0.7,
             color=COLOURS.get(run.component, "k"))
    top.fill_between(t, -peak * run.taper, peak * run.taper, color="0.85",
                     zorder=0, label="window")
    top.set_xlim(t[0], t[-1])
    top.set_xlabel("time after the origin (min)")
    top.set_ylabel(f"{run.component}, {run.quantity} [{UNITS[run.quantity]}]")
    top.legend(fontsize=8, loc="upper right")
    top.set_title(run.description, fontsize=9)
    lowest = 0.0 if deepest is None else surface - deepest
    axes = []
    for j, name in enumerate(names):
        ax = fig.add_subplot(grid[1, j], sharey=axes[0] if axes else None)
        axes.append(ax)
        values = run.kernels[name]
        keep = run.radius >= lowest
        profile(ax, _profile(run.radius[keep], values[keep], boundaries,
                             name, lowest), joiners=True, boundaries=False,
                color="k", lw=0.9)
        for b in boundaries:
            if b >= lowest:
                ax.axhline(b, **BOUNDARY_STYLE)
        for b in phase:
            if b >= lowest:
                ax.axhline(b, **{k: v for k, v in GUIDE.items()
                                 if k != "zorder"})
        ax.axvline(0, color="0.6", lw=0.6)
        point = run.at_source[name]
        ax.plot(0, run.source_radius, **STAR)
        if point != 0:
            point_unit = ("s" if run.relative else
                          f"s per ({_parameter_unit(name)})")
            ax.annotate(f"P = {point:.2g} {point_unit}",
                        (0, run.source_radius), xytext=(4, -12),
                        textcoords="offset points", fontsize=7)
        largest = np.abs(values[keep]).max(initial=0.0)
        limit = 1.1 * largest if largest else 1.0
        ax.set_xlim(-limit, limit)
        ax.ticklabel_format(axis="x", style="sci", scilimits=(-2, 3))
        ax.tick_params(labelsize=8)
        ax.xaxis.get_offset_text().set_fontsize(8)
        ax.set_title(("ln " if run.relative else "") + name)
        ax.set_xlabel(_kernel_unit(name, run.relative), fontsize=8)
        if j > 0:
            ax.tick_params(labelleft=False)
    axes[0].set_ylabel("radius (km)")
    axes[0].set_ylim(lowest, surface)
    change = "δ ln X" if run.relative else "δX"
    fig.text(0.5, 0.005,
             f"★ source; P is the concentrated source term.\n"
             f"Kernel action: ∫K {change} r² dr + P {change}(source).\n"
             "Each parameter panel uses its own horizontal SI scale.\n"
             "Small plotted terms do not establish significance.",
             ha="center", va="bottom", fontsize=7)
    fig.tight_layout(rect=(0, 0.095, 1, 1))
    return fig
