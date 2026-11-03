"""Pictures of snapshots of the wavefield: sections through the earth and
maps around the source.

A grid of depths along a great circle (one azimuth, or two opposite) is a
section, drawn as planetmodel draws one, through a `SectionSlice` made from
the snapshot; a grid at one depth and many azimuths is a map, drawn in
azimuthal equidistant projection about the source with cartopy, so that its
distances and azimuths are those of the grid. Anything else is drawn as
lines against distance.
"""

import math
import textwrap
from dataclasses import replace

import numpy as np
from planetmodel.deck import MINEOS, read_deck
from planetmodel.plotting import colour_scale, section
from planetmodel.slices import SectionSlice, Shown, SliceLayer

from .runs import DIMENSIONS
from ._labels import clipping_stats

#: The guides drawn at the boundaries between solid and fluid, over the
#: picture.
GUIDE = {"color": "0.25", "lw": 0.6, "alpha": 0.7, "zorder": 3}

#: The mark of the source.
STAR = {"marker": "*", "color": "gold", "markeredgecolor": "k",
        "markersize": 11, "zorder": 4}

#: The radius of the sphere on which the maps are drawn (m).
RADIUS = 6371e3


def _unit(distance, azimuth):
    """Unit vectors at angular distances and azimuths from a source at the
    pole, as YSpec places them: longitude pi minus the azimuth."""
    phi = np.pi - azimuth
    return np.stack([np.sin(distance) * np.cos(phi),
                     np.sin(distance) * np.sin(phi),
                     np.cos(distance) * np.ones_like(phi)], axis=-1)


def _section_slice(run, field, surface, shown):
    """The snapshot on a great circle as a slice: along the second azimuth
    from its far end to the source, then along the first, radii upward."""
    delta = np.radians(run.distances)
    if len(run.azimuths) == 2:
        separation = (run.azimuths[1] - run.azimuths[0]) % 360.0
        if not math.isclose(separation, 180.0, abs_tol=1e-6):
            raise ValueError("a snapshot section needs two opposite azimuths")
    pieces = [(delta, np.full_like(delta, math.radians(run.azimuths[0])),
               field[:, 0, :])]
    if len(run.azimuths) == 2:
        pieces.insert(0, (delta[::-1],
                          np.full_like(delta, math.radians(run.azimuths[1])),
                          field[:, 1, ::-1]))
    path = np.concatenate([_unit(d, a) for d, a, _ in pieces])
    along = (np.concatenate([np.pi - pieces[0][0], np.pi + delta])
             if len(pieces) == 2 else delta)
    along = along - along[0]
    values = np.concatenate([v for _, _, v in pieces], axis=1)[::-1]
    radius = (surface - run.depths)[::-1]
    if len(pieces) == 2:
        # The two azimuths both contain the source sample when distance zero
        # is requested. Keep the raw samples in SnapshotRun, but render a
        # deliberate gap at the origin so no colour is interpolated between
        # the two component bases there.
        if math.isclose(delta[0], 0.0, abs_tol=1e-12):
            source_columns = (len(delta) - 1, len(delta))
            path = np.delete(path, source_columns, axis=0)
            along = np.delete(along, source_columns)
            values = np.delete(values, source_columns, axis=1)
        source_along = float(delta[-1])
        insert = int(np.searchsorted(along, source_along))
        path = np.insert(path, insert, np.array([0.0, 0.0, 1.0]), axis=0)
        along = np.insert(along, insert, source_along)
        values = np.insert(values, insert, np.nan, axis=1)
    grid = np.broadcast_to(radius[:, None], values.shape)
    points = grid[..., None] * path[None, :, :]
    travel = path[1] - path[0]
    travel = travel - np.dot(travel, path[0]) * path[0]
    plane = np.stack([path[0], travel / np.linalg.norm(travel)])
    # Split existing radial samples at every duplicated radius in the model
    # deck. No value is moved to an interface or interpolated to its other
    # side.
    deck = read_deck(run.model, MINEOS)
    duplicate = np.flatnonzero(np.diff(deck.radius) == 0)
    boundaries = np.unique(deck.radius[duplicate]) / 1e3
    shell = np.searchsorted(boundaries, radius, side="right")
    if boundaries.size:
        # DepthSolves gives each exact interface radius to the solver's
        # receiver constructor. It uses the solid side at solid-fluid
        # boundaries and Side::Below (the lower-radius layer) otherwise.
        # Keep each saved sample in that same shell; do not move or
        # interpolate its value to either side.
        shear = deck["vsv"] if "vsv" in deck else deck["vs"]
        for index in duplicate:
            boundary = deck.radius[index] / 1e3
            matches = np.isclose(radius, boundary, rtol=0, atol=1e-8)
            inner_fluid = shear[index] == 0
            outer_solid = shear[index + 1] != 0
            if inner_fluid and outer_solid:
                shell[matches] = np.searchsorted(boundaries, boundary,
                                                 side="right")
            else:
                shell[matches] = np.searchsorted(boundaries, boundary,
                                                 side="left")
    layers = []
    for index in np.unique(shell):
        rows = shell == index
        if not np.any(rows):
            continue
        layers.append(SliceLayer(index=int(index), radius=np.array(grid[rows]),
                                 values=values[rows],
                                 reference=points[rows],
                                 physical=points[rows]))
    return SectionSlice(distance=along, path=path, plane=plane,
                        layers=tuple(layers), shown=shown)


def _phase_boundaries(model):
    """The radii (km) of the surface and of every boundary between solid and
    fluid in a deck."""
    deck = read_deck(model, MINEOS)
    shear = deck["vsv"] if "vsv" in deck else deck["vs"]
    r = deck.radius
    radii = [r[-1]]
    for i in range(1, len(r)):
        if r[i] == r[i - 1] and (shear[i] == 0) != (shear[i - 1] == 0):
            radii.append(r[i])
    return np.array(radii) / 1e3, r[-1] / 1e3


def _draw_section(ax, run, field, surface, boundaries, shown, limits,
                  cmap):
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    from planetmodel.plotting import colour_extend

    slice_ = _section_slice(run, field, surface, shown)
    drawable = tuple(layer for layer in slice_.layers
                     if layer.radius.shape[0] >= 2)
    isolated = tuple(layer for layer in slice_.layers
                     if layer.radius.shape[0] == 1)
    meshes = []
    if drawable:
        meshes = section(ax, replace(slice_, layers=drawable), diverging=True,
                         vmin=-limits, vmax=limits, interfaces=False)
        norm, cmap = meshes[0].norm, meshes[0].cmap
    else:
        norm = Normalize(-limits, limits)
    for layer in isolated:
        angle = np.broadcast_to(slice_.distance, layer.radius.shape)
        turn = angle - 0.5 * slice_.length
        ax.scatter(turn.ravel(), layer.radius.ravel(),
                   c=layer.values.ravel(), cmap=cmap, norm=norm,
                   s=2, linewidths=0)
    if not meshes:
        meshes = [ScalarMappable(norm=norm, cmap=cmap)]
    all_radii = np.concatenate([layer.radius.ravel()
                                for layer in slice_.layers])
    if all_radii.size:
        ax.set_ylim(float(all_radii.min()), float(all_radii.max()))
    half = 0.5 * slice_.length
    arc = np.linspace(-half, half, 721)
    for r in boundaries:
        ax.plot(arc, np.full_like(arc, r), **GUIDE)
    at = 0.0 if len(run.azimuths) == 2 else -half
    ax.plot(at, surface - run.hypocentre["depth_km"], **STAR)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    return meshes[0]


def _draw_map(ax, run, field, projection, limits, cmap, coastlines):
    import cartopy.crs as ccrs

    image = field[0]
    latitude, longitude = run.latitude, run.longitude
    if not (np.allclose(latitude[0], latitude[-1], atol=1e-6, rtol=0) and
            np.allclose(np.sin(np.radians(longitude[0])),
                        np.sin(np.radians(longitude[-1])), atol=1e-6, rtol=0) and
            np.allclose(np.cos(np.radians(longitude[0])),
                        np.cos(np.radians(longitude[-1])), atol=1e-6, rtol=0)):
        latitude = np.concatenate([latitude, latitude[:1]], axis=0)
        longitude = np.concatenate([longitude, longitude[:1]], axis=0)
        image = np.concatenate([image, image[:1]], axis=0)
    mesh = ax.pcolormesh(longitude, latitude, image,
                         transform=ccrs.PlateCarree(), shading="nearest",
                         cmap=cmap, vmin=-limits, vmax=limits)
    if coastlines:
        ax.coastlines(lw=0.4, color="0.3")
    ax.plot(0, 0, transform=projection, **STAR)
    ax.set_global()
    return mesh


def snapshots(run, component="up", *, common=False, clip=99.5,
              coastlines=True):
    """The snapshots of a run at each of its times, a panel for each, each
    scaled to the `clip` percentile of its own values (or all alike, with
    `common`)."""
    import matplotlib.pyplot as plt

    shown = Shown(f"{component} {run.quantity}",
                  dimensions=DIMENSIONS[run.quantity])
    fields = [run.field(k, component) for k in range(len(run.times))]
    cmap, _, common_limit = colour_scale(fields, shown, diverging=True,
                                         clip=clip)
    columns = min(len(run.times), 3)
    rows = math.ceil(len(run.times) / columns)
    is_section = len(run.depths) > 1 and len(run.azimuths) <= 2
    is_map = len(run.depths) == 1 and len(run.azimuths) > 2
    fig = plt.figure(figsize=(4.2 * columns, 4.4 * rows + 0.8))
    if is_section:
        boundaries, surface = _phase_boundaries(run.model)
        kw = {"projection": "polar"}
    elif is_map:
        import cartopy.crs as ccrs
        globe = ccrs.Globe(ellipse=None, semimajor_axis=RADIUS,
                           semiminor_axis=RADIUS)
        projection = ccrs.AzimuthalEquidistant(
            central_longitude=run.hypocentre["longitude"],
            central_latitude=run.hypocentre["latitude"], globe=globe)
        kw = {"projection": projection}
    else:
        kw = {}
    axes = []
    mappables = []
    panel_stats = []
    for k, (time, field) in enumerate(zip(run.times, fields)):
        ax = fig.add_subplot(rows, columns, k + 1, **kw)
        axes.append(ax)
        limits = (common_limit if common else
                  colour_scale([field], shown, diverging=True, clip=clip)[2])
        finite = np.abs(field[np.isfinite(field)])
        panel_stats.append((float(finite.max()) if finite.size else 0.0,
                            float(100 * np.mean(finite > limits)) if finite.size else 0.0,
                            limits))
        if is_section:
            mappable = _draw_section(ax, run, field, surface, boundaries,
                                     shown, limits, cmap)
            ax.set_title(f"{time:g} s")
        elif is_map:
            mesh = _draw_map(ax, run, field, projection, limits, cmap,
                             coastlines)
            mappable = mesh
            ax.set_title(f"{time:g} s, {run.depths[0]:g} km deep")
        else:
            for d, depth in enumerate(run.depths):
                for a, azimuth in enumerate(run.azimuths):
                    ax.plot(run.distances, field[d, a],
                            label=f"{depth:g} km, {azimuth:g}°")
            ax.set_xlabel("distance (degrees)")
            ax.set_ylabel(shown.label)
            ax.legend(fontsize=7)
            ax.set_title(f"{time:g} s")
            mappable = None
        mappables.append(mappable)
        if mappable is not None:
            maximum, fraction, panel_limit = panel_stats[-1]
            ax.text(0.02, 0.02,
                    clipping_stats(panel_limit, maximum, fraction),
                    transform=ax.transAxes, fontsize=6.5, va="bottom",
                    bbox={"facecolor": "white", "alpha": 0.8,
                           "edgecolor": "none"})
    scaling = (("common colour scale across times" if common else
                "independent colour scale per time")
               if is_section or is_map else
               "line amplitudes shown in native units")
    display_note = (
        "At exact interfaces, samples use the solver's solid side, or the "
        "deeper side otherwise. Gouraud shading interpolates colours between "
        "saved samples; it adds no numerical resolution."
        if is_section else
        "Saved map samples are assigned to nearest-colour cells; amplitudes "
        "are not interpolated."
        if is_map else "")
    if is_section or is_map:
        fig.suptitle(f"{shown.label}, source {run.source}",
                     y=0.995 if is_section else 0.98)
        caption = textwrap.fill(f"{scaling}. {display_note}", width=58)
        fig.text(0.5, 0.955 if is_section else 0.965, caption,
                 ha="center", va="top", fontsize=8)
        # Separate the figure title, wrapped caption and panel titles on
        # sections with the interface-side note.
        fig.subplots_adjust(top=0.84 if is_section else 0.88)
    else:
        fig.suptitle(f"{shown.label}, source {run.source} ({scaling})")
    if is_section or is_map:
        from planetmodel.plotting import colour_extend

        if common:
            mappable = plt.cm.ScalarMappable(
                norm=plt.Normalize(-common_limit, common_limit), cmap=cmap)
            fig.colorbar(mappable, ax=axes, shrink=0.6, label=shown.label,
                         extend=colour_extend(fields, -common_limit,
                                              common_limit))
        else:
            for index, (ax, (_, _, limit), mappable) in enumerate(
                    zip(axes, panel_stats, mappables)):
                ax.figure.colorbar(mappable, ax=ax, shrink=0.72,
                                   label=shown.label,
                                   extend=colour_extend(
                                       [fields[index]], -limit, limit))
    return fig
