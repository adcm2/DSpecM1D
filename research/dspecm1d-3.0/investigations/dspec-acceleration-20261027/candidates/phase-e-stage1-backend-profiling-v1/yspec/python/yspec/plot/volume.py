"""Pictures of three-dimensional kernels, from the model file yspec-kernels
writes: a section through the earth along the great circle through the
source and the receiver, and maps at some depths, drawn by planetmodel."""

import math
from dataclasses import replace

import numpy as np
from planetmodel.modelfile import read_model
from planetmodel.plotting import section, sphere_map
from planetmodel.sampling import AngularGrid
from planetmodel.slices import section as section_slice
from planetmodel.slices import Shown, sphere

from .kernels import _kernel_unit
from ._labels import clipping_stats
from .snapshots import STAR

#: The mark of the receiver.
TRIANGLE = {"marker": "v", "color": "white", "markeredgecolor": "k",
            "markersize": 9, "zorder": 4}


def _unit(latitude, longitude):
    """The unit vector of a geographic latitude and longitude (degrees), at
    its geocentric latitude, as YSpec places it."""
    theta = math.pi / 2 - math.atan(0.993277 * math.tan(math.radians(latitude)))
    phi = math.radians(longitude)
    return np.array([math.sin(theta) * math.cos(phi),
                     math.sin(theta) * math.sin(phi), math.cos(theta)])


def _spherical(x):
    return math.acos(max(-1.0, min(1.0, x[2]))), math.atan2(x[1], x[0])


def _geographic_great_circle(source, receiver, samples=361):
    """Geographic coordinates of the great circle used by YSpec's field
    grid, sampled from the source to the receiver."""
    s = _unit(*source)
    r = _unit(*receiver)
    normal = np.cross(s, r)
    sine = np.linalg.norm(normal)
    cosine = float(np.dot(s, r))
    distance = math.atan2(sine, cosine)
    if sine < 1e-12 and cosine < 0:
        raise ValueError("source and receiver do not define a unique "
                         "great-circle path")
    fraction = np.linspace(0.0, 1.0, samples)
    if sine < 1e-12:
        points = np.broadcast_to(s, (samples, 3)).copy()
    else:
        normal /= sine
        points = (np.cos(fraction * distance)[:, None] * s
                  + np.sin(fraction * distance)[:, None]
                  * np.cross(normal, s))
    geocentric = np.arcsin(np.clip(points[:, 2], -1.0, 1.0))
    latitude = np.degrees(np.arctan2(np.sin(geocentric),
                                     0.993277 * np.cos(geocentric)))
    longitude = np.degrees(np.unwrap(np.arctan2(points[:, 1], points[:, 0])))
    return longitude, latitude


def _geographic_slice(slice_):
    """Adapt a model slice's geocentric grid to the geographic map CRS."""
    geocentric = 90.0 - np.degrees(slice_.grid.colatitudes)
    geographic = np.degrees(np.arctan(
        np.tan(np.radians(geocentric)) / 0.993277))
    grid = AngularGrid(np.radians(90.0 - geographic),
                       slice_.grid.longitudes, kind="custom")
    return replace(slice_, grid=grid)


def _annotate(ax, values, mesh):
    values = np.asarray(values)
    finite = np.abs(values[np.isfinite(values)])
    maximum = float(finite.max()) if finite.size else 0.0
    limit = max(abs(float(mesh.norm.vmin)), abs(float(mesh.norm.vmax)))
    fraction = float(100 * np.mean(finite > limit)) if finite.size else 0.0
    ax.text(0.02, 0.02, clipping_stats(limit, maximum, fraction),
            transform=ax.transAxes, fontsize=6.5, va="bottom",
            bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none"})


def _arc(source, receiver, margin):
    """The great circle from `margin` radians behind the source to as far
    beyond the receiver: its start (theta, phi), its azimuth there (radians
    from north toward east), its length, and the distance of the source and
    receiver along it."""
    s, r = source, receiver
    distance = math.atan2(np.linalg.norm(np.cross(s, r)), np.dot(s, r))
    normal = np.cross(s, r)
    normal /= np.linalg.norm(normal)
    across = np.cross(normal, s)
    start = s * math.cos(margin) - across * math.sin(margin)
    travel = s * math.sin(margin) + across * math.cos(margin)
    theta, phi = _spherical(start)
    north = -np.array([math.cos(theta) * math.cos(phi),
                       math.cos(theta) * math.sin(phi), -math.sin(theta)])
    east = np.array([-math.sin(phi), math.cos(phi), 0.0])
    azimuth = math.atan2(np.dot(travel, east), np.dot(travel, north))
    length = min(distance + 2 * margin, 2 * math.pi)
    return (theta, phi), azimuth, length, margin, margin + distance


def volume_kernels(run, parameter=None, depths=(100, 500, 1500, 2500),
                   margin=20.0, clip=98.0):
    """The kernel of a parameter (by default vs, or the first) on a section
    along the great circle through the source and receiver, extended by
    `margin` degrees beyond each, and on maps at the depths given (km)."""
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt

    model = read_model(run.volume)
    names = [n.removeprefix("K_") for n in model.field_names()]
    if parameter is None:
        parameter = "vs" if "vs" in names else names[0]
    if parameter not in names:
        raise ValueError(f"there is no kernel of {parameter}: there are "
                         f"{', '.join(names)}")
    field = f"K_{parameter}"
    source = run.hypocentre
    receiver = run.receiver_location
    s = _unit(source["latitude"], source["longitude"])
    r = _unit(*receiver)
    start, azimuth, length, at_source, at_receiver = _arc(
        s, r, math.radians(margin))
    surface = float(model.skeleton.boundaries[-1])

    # The section on the left, the maps two to a row beside it.
    rows = max(math.ceil(len(depths) / 2), 1)
    fig = plt.figure(figsize=(16, 4.2 * rows + 1))
    grid = fig.add_gridspec(rows, 3, width_ratios=[1.5, 1, 1])
    ax = fig.add_subplot(grid[:, 0], projection="polar")
    cut = section_slice(model, field, start, azimuth=azimuth, length=length,
                        m=721, n=60)
    shown = Shown(f"K, {'ln ' if run.relative else ''}{parameter} "
                  f"[{_kernel_unit(parameter, run.relative)}]")
    cut = replace(cut, shown=shown)
    section_meshes = section(ax, cut, scale=1e-3, diverging=True, clip=clip,
                             colorbar=True)
    _annotate(ax, np.concatenate([layer.values.ravel()
                                  for layer in cut.layers]),
              section_meshes[0])
    half = 0.5 * length
    ax.plot(at_source - half, (surface - source["depth_km"] * 1e3) / 1e3,
            **STAR)
    ax.plot(at_receiver - half, (surface - run.receiver_depth * 1e3) / 1e3,
            **TRIANGLE)
    ax.set_xticks([])
    ax.set_yticks([])
    measured = run.description.removeprefix("radial kernels of ")
    fig.suptitle(f"The kernel of {'ln ' if run.relative else ''}{parameter} "
                 f"for {measured}, on the great circle through the source "
                 f"(star) and the receiver (triangle), and at depths;\n"
                 f"colours clipped at the {clip:g}th percentile", fontsize=10)

    # The maps are of the hemisphere centred on the middle of the path.
    middle = (s + r) / np.linalg.norm(s + r)
    theta, phi = _spherical(middle)
    geocentric_middle = 90.0 - math.degrees(theta)
    geographic_middle = math.degrees(math.atan(
        math.tan(math.radians(geocentric_middle)) / 0.993277))
    projection = ccrs.Orthographic(central_longitude=math.degrees(phi),
                                   central_latitude=geographic_middle)
    for j, depth in enumerate(depths):
        ax = fig.add_subplot(grid[j // 2, 1 + j % 2], projection=projection)
        try:
            cut = sphere(model, field, depth=depth * 1e3)
        except ValueError:
            # On a discontinuity: the side above it.
            cut = sphere(model, field, depth=depth * 1e3, side="above")
        cut = _geographic_slice(replace(cut, shown=shown))
        mesh = sphere_map(ax, cut, diverging=True, clip=clip,
                          coastlines=True, colorbar=True)
        _annotate(ax, cut.values, mesh)
        path_lon, path_lat = _geographic_great_circle(
            (source["latitude"], source["longitude"]), receiver)
        ax.plot(path_lon, path_lat, color="k", lw=0.6,
                transform=ccrs.PlateCarree())
        ax.plot(source["longitude"], source["latitude"],
                transform=ccrs.PlateCarree(), **STAR)
        ax.plot(receiver[1], receiver[0], transform=ccrs.PlateCarree(),
                **TRIANGLE)
        ax.set_title(f"{depth:g} km deep", fontsize=9)
    fig.tight_layout()
    return fig
