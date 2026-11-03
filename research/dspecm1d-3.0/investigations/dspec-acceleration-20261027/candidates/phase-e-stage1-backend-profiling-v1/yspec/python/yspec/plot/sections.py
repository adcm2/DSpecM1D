"""Pictures of three-dimensional kernels on sections and maps about the
path from the source to the receiver, with the seismogram measured: the
banana of finite-frequency theory along the path, and its doughnut across
it, drawn as planetmodel draws sections.

Near the source and the receiver the kernels are far larger than elsewhere;
the colours are scaled to the kernel away from them, at more than `clear`
degrees, and saturate there.
"""

import math

import numpy as np
from planetmodel.plotting import section
from planetmodel.slices import SectionSlice, Shown, SliceLayer

from .runs import KernelRun, UNITS
from .kernels import _kernel_unit
from ._labels import clipping_stats
from .seismograms import COLOURS
from .snapshots import GUIDE, STAR, _phase_boundaries
from .volume import TRIANGLE


def _parameter(run, parameter):
    """The kernel drawn: as asked, or vs for transverse/shear measurements."""
    if parameter is None:
        parameter = "vs" if run.component == "T" or "S" in run.name else "vp"
    if parameter not in run.parameters:
        raise ValueError(f"there is no kernel of {parameter}: there are "
                         f"{', '.join(run.parameters)}")
    return parameter


def _slice(sampled, coordinate, parameter, name, surface):
    """A section as planetmodel's slice, drawn in a plane of its own: the
    distance along it from its start (radians), and each layer's radii (km)
    by places."""
    angle = np.radians(coordinate - coordinate[0])
    path = np.stack([np.cos(angle), np.sin(angle), np.zeros_like(angle)],
                    axis=-1)
    layers = []
    for layer in sampled.layers:
        radius = np.broadcast_to((surface - layer.depth)[:, None],
                                 layer.values[parameter].shape)
        points = radius[..., None] * path[None, :, :]
        layers.append(SliceLayer(index=layer.index, radius=np.array(radius),
                                 values=layer.values[parameter],
                                 reference=points, physical=points))
    plane = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    return SectionSlice(distance=angle, path=path, plane=plane,
                        layers=tuple(layers), shown=Shown(name))


def _limit(run, sampled, parameter, keep, clip, shallowest=150.0):
    """The colour limit: the clip percentile of the kernel's size at the
    places kept, below the shallowest depth (km), above which the crust's
    large sensitivity would set it."""
    values = np.concatenate([
        layer.values[parameter][layer.depth > shallowest][:, keep].ravel()
        for layer in sampled.layers])
    values = np.abs(values[np.isfinite(values)])
    return float(np.percentile(values, clip)) if values.size else 1.0


def _deepest(run):
    """How deep a section is drawn: through the core for a phase with a leg
    in it, otherwise to just below the core-mantle boundary."""
    return 0.0 if "K" in run.name else 3000.0


def _extent_label(run):
    """The radial extent shown in a section (km)."""
    deepest = _deepest(run)
    return ("full radius" if deepest == 0 else
            f"r ≥ {deepest:g} km")


def _draw(ax, run, sampled, coordinate, parameter, limit, boundaries,
          marks, stats):
    """A section in its plane, with the boundaries between solid and fluid
    and the marks (angle along it, radius, style) of the source and
    receiver, cropped to the depths the phase samples."""
    label = (f"K, {'ln ' if run.relative else ''}{parameter} "
             f"[{_kernel_unit(parameter, run.relative)}]")
    cut = _slice(sampled, coordinate, parameter, label, run.surface)
    section(ax, cut, physical=True, diverging=True, vmin=-limit,
            vmax=limit, interfaces=False, colorbar=True)
    maximum, fraction = stats
    ax.text(0.02, 0.02, clipping_stats(limit, maximum, fraction),
            transform=ax.transAxes, fontsize=6.5, va="bottom",
            bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"})
    half = 0.5 * cut.length
    deepest = _deepest(run)
    arc = np.linspace(-half, half, 721)
    for r in [run.surface] + [b for b in boundaries if b > deepest]:
        ax.plot(r * np.sin(arc), r * np.cos(arc), **GUIDE)
    for side in (-half, half):
        ax.plot([deepest * math.sin(side), run.surface * math.sin(side)],
                [deepest * math.cos(side), run.surface * math.cos(side)],
                **GUIDE)
    for angle, radius, style in marks:
        turn = math.radians(angle - coordinate[0]) - half
        ax.plot(radius * math.sin(turn), radius * math.cos(turn), **style)
    width = run.surface * math.sin(min(half, math.pi / 2))
    bottom = deepest * math.cos(half) if half < math.pi / 2 else -run.surface
    if deepest == 0:
        bottom = min(0.0, run.surface * math.cos(half))
    ax.set_xlim(-1.02 * width, 1.02 * width)
    ax.set_ylim(bottom - 0.02 * run.surface, 1.02 * run.surface)
    ax.set_aspect("equal")
    ax.set_axis_off()


def _trace(ax, run):
    """The seismogram measured, its window shaded."""
    t = run.time / 60
    peak = np.abs(run.seismogram).max() or 1.0
    inside = run.taper > 0
    ax.plot(t, run.seismogram / peak, lw=0.5, color="0.6")
    ax.plot(t[inside], run.seismogram[inside] / peak, lw=0.9,
            color=COLOURS.get(run.component, "k"))
    ax.fill_between(t, -run.taper, run.taper, color="0.9", zorder=0)
    ax.set_xlim(t[0], t[-1])
    ax.set_ylim(-1.05, 1.05)
    ax.set_yticks([])
    ax.set_xlabel("time (min)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.text(0.01, 0.92, f"{run.name} on {run.component}",
            transform=ax.transAxes, va="top", fontsize=9)
    ax.text(0.99, 0.92,
            f"normalised amplitude / peak; peak {peak:.3g} "
            f"{UNITS[run.quantity]}",
            transform=ax.transAxes, va="top", ha="right", fontsize=7)


def _marks(run):
    surface = run.surface
    return ([(0.0, run.source_radius, STAR),
             (run.distance, surface - run.receiver_depth, TRIANGLE)],
            [])


def _stats(sampled, parameter, limit):
    values = np.concatenate([layer.values[parameter].ravel()
                             for layer in sampled.layers])
    finite = np.abs(values[np.isfinite(values)])
    if finite.size == 0:
        return 0.0, 0.0
    return float(finite.max()), float(100.0 * np.mean(finite > limit))


def banana(run, parameter=None, *, clear=5.0, clip=99.0):
    """One measurement: its seismogram and window, and its kernel on the
    sections along and across the path."""
    import matplotlib.pyplot as plt

    parameter = _parameter(run, parameter)
    path, across = run.sections()
    boundaries, _ = _phase_boundaries(run.model)
    keep = (np.abs(path.along) > clear) & (np.abs(path.along - run.distance) >
                                           clear)
    limit = _limit(run, path, parameter, keep, clip)
    path_stats = _stats(path, parameter, limit)
    across_limit = _limit(run, across, parameter,
                          np.ones(len(across.across), bool), clip)
    across_stats = _stats(across, parameter, across_limit)
    fig = plt.figure(figsize=(14, 9))
    grid = fig.add_gridspec(2, 2, height_ratios=[1, 4],
                            width_ratios=[2.2, 1])
    _trace(fig.add_subplot(grid[0, :]), run)
    on_path, _ = _marks(run)
    ax = fig.add_subplot(grid[1, 0])
    _draw(ax, run, path, path.along, parameter, limit, boundaries, on_path,
          path_stats)
    ax.set_title(f"along the path\n{_extent_label(run)}", fontsize=9)
    ax = fig.add_subplot(grid[1, 1])
    _draw(ax, run, across, across.across, parameter, across_limit, boundaries,
          [], across_stats)
    ax.set_title(f"across it, at its middle\n{_extent_label(run)}",
                 fontsize=9)
    measured = run.description.removeprefix("radial kernels of ")
    fig.suptitle(f"The kernel of {'ln ' if run.relative else ''}{parameter} "
                 f"for {measured}; colours saturate near the source (star) "
                 f"and the receiver (triangle)", fontsize=10)
    fig.text(0.5, 0.008,
             f"Each panel has its own symmetric {clip:g}th-percentile limit; "
             f"limits exclude depths above 150 km, and the along-path limit "
             f"also excludes points within {clear:g}° of source and receiver. "
             "Section colours are interpolated for display between saved samples.",
             ha="center", va="bottom", fontsize=7)
    fig.tight_layout()
    return fig


def gallery(path, parameter=None, *, clear=5.0, clip=99.0):
    """Every measurement of a run, a row for each: its seismogram and
    window, and its kernel along the path and across it."""
    import matplotlib.pyplot as plt

    names = KernelRun(path).names
    fig = plt.figure(figsize=(16, 3.0 * len(names)))
    grid = fig.add_gridspec(len(names), 3, width_ratios=[1.1, 2.0, 0.8])
    for row, name in enumerate(names):
        run = KernelRun(path, name)
        chosen = _parameter(run, parameter)
        along, across = run.sections()
        boundaries, _ = _phase_boundaries(run.model)
        keep = (np.abs(along.along) > clear) & (
            np.abs(along.along - run.distance) > clear)
        along_limit = _limit(run, along, chosen, keep, clip)
        along_stats = _stats(along, chosen, along_limit)
        across_limit = _limit(run, across, chosen,
                              np.ones(len(across.across), bool), clip)
        across_stats = _stats(across, chosen, across_limit)
        _trace(fig.add_subplot(grid[row, 0]), run)
        on_path, _ = _marks(run)
        ax = fig.add_subplot(grid[row, 1])
        _draw(ax, run, along, along.along, chosen, along_limit, boundaries,
              on_path, along_stats)
        ax.set_title(f"{run.name}: {'ln ' if run.relative else ''}{chosen}, "
                     f"along the path\n{_extent_label(run)}", fontsize=9)
        ax = fig.add_subplot(grid[row, 2])
        _draw(ax, run, across, across.across, chosen, across_limit,
              boundaries, [], across_stats)
        ax.set_title(f"across it, at its middle\n{_extent_label(run)}",
                     fontsize=9)
    fig.suptitle("Traveltime kernels by cross-correlation in the windows "
                 "shaded; colours saturate near the source (star) and "
                 "receiver (triangle)", fontsize=11)
    fig.text(0.5, 0.003,
             f"Each panel has its own symmetric {clip:g}th-percentile limit; "
             f"limits exclude depths above 150 km, and along-path limits "
             f"also exclude points within {clear:g}° of source and receiver. "
             "Section colours are interpolated for display between saved samples.",
             ha="center", va="bottom", fontsize=7)
    fig.tight_layout()
    return fig
