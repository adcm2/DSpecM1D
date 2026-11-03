"""Small regressions for plot interpretation and coordinate contracts."""

import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

pytest.importorskip("matplotlib")
pytest.importorskip("planetmodel")

import matplotlib
matplotlib.use("Agg")
import numpy as np
from planetmodel.deck import MINEOS, read_deck
from planetmodel.slices import Shown

from yspec.plot.runs import DIMENSIONS, SeismogramRun, SnapshotRun
from yspec.plot.sections import _parameter
from yspec.plot.seismograms import compare, traces
from yspec.plot.snapshots import _draw_map, _section_slice, snapshots
from yspec.plot.volume import _geographic_slice


class TraceReceiver:
    def __init__(self, number=1, *, times=None, units="m", latitude=0.0):
        self.number = number
        self.file = Path("unused")
        self.units = units
        self.latitude = latitude
        self.longitude = 20.0
        self.depth = 3.0
        self.distance = 30.0
        self.azimuth = 40.0
        self.back_azimuth = 220.0
        self.times = np.array([0.0, 1.0, 2.0] if times is None else times)

    def trace(self, component):
        values = {"Z": [1.0, -0.5, 0.25], "R": [0.01, 0.02, 0.0],
                  "T": [0.0, 0.0, 0.0]}[component]
        return self.times, np.asarray(values)


class TraceRun:
    def __init__(self, *, time=None, quantity="displacement", units="m",
                 latitude=0.0):
        self.receivers = [TraceReceiver(times=time, units=units,
                                        latitude=latitude)]
        self.quantity = quantity
        self.source = "source"
        self.path = Path("run.json")

    def receiver(self, number):
        return self.receivers[number - 1]


@pytest.mark.parametrize("requested, expected", [(None, "vs"), ("vp", "vp")])
def test_transverse_shear_packet_defaults_to_vs_but_honours_explicit_parameter(
        requested, expected):
    run = SimpleNamespace(component="T", name="shear-packet",
                          parameters=("vp", "vs"))
    assert _parameter(run, requested) == expected


def test_traces_share_component_scale_by_default_and_disclose_independent():
    run = TraceRun()
    fig = traces(run, components="ZR")
    zlim = fig.axes[0].get_ylim()
    rlim = fig.axes[1].get_ylim()
    assert zlim == rlim
    assert "share a scale" in fig._suptitle.get_text()

    independent = traces(run, components="ZR", scale="independent")
    assert independent.axes[0].get_ylim() != independent.axes[1].get_ylim()
    assert "not comparable" in independent._suptitle.get_text()


@pytest.mark.parametrize(
    "second, message",
    [(TraceRun(time=[0.0, 2.0, 4.0]), "time grids"),
     (TraceRun(quantity="velocity", units="m/s"), "quantities"),
     (TraceRun(latitude=0.1), "geometry"),
     (TraceRun(units="COUNTS"), "units")])
def test_compare_refuses_incompatible_physical_traces(second, message):
    with pytest.raises(ValueError, match=message):
        compare([TraceRun(), second])


def test_instrument_plot_unit_comes_from_sac_metadata(tmp_path):
    (tmp_path / "station.pz").write_text(
        "* INPUT UNIT : M/S\n* OUTPUT UNIT : COUNTS\nZEROS 0\nPOLES 0\n")
    (tmp_path / "out.receivers").write_text(
        "# file source latitude longitude depth_km distance_deg azimuth_deg back_azimuth_deg\n"
        "out.1 1 12 30 3 20 40 220\n")
    config = {
        "output": {"prefix": "out", "quantity": "velocity"},
        "sources": [{"name": "1"}],
        "receivers": [{"depth_km": 3, "locations": [[12, 30]],
                       "instrument": "station.pz"}],
    }
    path = tmp_path / "input.json"
    path.write_text(json.dumps(config))
    run = SeismogramRun(path)
    assert run.receiver(1).units == "COUNTS"


def test_first_seismogram_plot_uses_the_selected_receivers_units(tmp_path):
    (tmp_path / "broadband.pz").write_text(
        "  *  INPUT UNIT : M/S\n   *   OUTPUT   UNIT : COUNTS  \n"
        "ZEROS 0\nPOLES 0\n")
    (tmp_path / "ground.1").write_text(
        "0 1 0 0\n1 -0.5 0 0\n")
    (tmp_path / "instrument.1").write_text(
        "0 1 0 0\n1 -0.5 0 0\n")
    (tmp_path / "out.receivers").write_text(
        "# file source latitude longitude depth_km distance_deg azimuth_deg back_azimuth_deg\n"
        "ground.1 1 12 30 3 20 40 220\n"
        "instrument.1 1 13 31 3 25 45 225\n")
    config = {
        "output": {"prefix": "out", "quantity": "displacement"},
        "sources": [{"name": "1"}],
        "receivers": [
            {"depth_km": 3, "locations": [[12, 30]]},
            {"depth_km": 3, "locations": [[13, 31]],
             "instrument": "broadband.pz"},
        ],
    }
    input_path = tmp_path / "input.json"
    input_path.write_text(json.dumps(config))
    script = Path(__file__).parents[2] / "docs" / "site" / "examples" / "plot_first_seismogram.py"
    env = {**os.environ, "PYTHONPATH": str(Path(__file__).parents[1]),
           "MPLCONFIGDIR": str(tmp_path / "mplconfig")}

    for trace, expected in (("ground.1", "Up component (m)"),
                            ("instrument.1", "Up component (COUNTS)")):
        output = tmp_path / f"{trace}.svg"
        completed = subprocess.run(
            [sys.executable, str(script), str(tmp_path / trace),
             str(output), str(input_path)], check=True, capture_output=True,
            text=True, env=env)
        assert expected in output.read_text()
        assert expected.split("(")[1].rstrip(")") in completed.stdout


class FakeMapAxes:
    def __init__(self):
        self.args = None
        self.kwargs = None

    def pcolormesh(self, *args, **kwargs):
        self.args, self.kwargs = args, kwargs
        return object()

    def coastlines(self, **kwargs):
        pass

    def plot(self, *args, **kwargs):
        pass

    def set_global(self):
        pass


def test_snapshot_map_draws_saved_geographic_coordinates():
    pytest.importorskip("cartopy")
    lat = np.array([[12.0, 13.0, 14.0], [20.0, 21.0, 22.0]])
    lon = np.array([[30.0, 40.0, 50.0], [120.0, 130.0, 140.0]])
    run = SimpleNamespace(latitude=lat, longitude=lon,
                          hypocentre={"longitude": 35.0})
    ax = FakeMapAxes()
    _draw_map(ax, run, np.ones((1, 2, 3)), projection=object(), limits=1,
              cmap="RdBu", coastlines=False)
    np.testing.assert_array_equal(ax.args[0], np.concatenate([lon, lon[:1]]))
    np.testing.assert_array_equal(ax.args[1], np.concatenate([lat, lat[:1]]))
    assert ax.kwargs["transform"].__class__.__name__ == "PlateCarree"
    assert ax.kwargs["shading"] == "nearest"


def test_snapshot_run_keeps_saved_latitude_and_longitude(tmp_path):
    coordinates = [(11.0, 21.0), (12.0, 22.0),
                   (31.0, 41.0), (32.0, 42.0)]
    rows = []
    for azimuth, pair in zip((0.0, 0.0, 180.0, 180.0), coordinates):
        distance = 0.0 if len(rows) % 2 == 0 else 30.0
        rows.append([3.0, distance, azimuth, *pair, 0.0, 0.0, 0.0])
    (tmp_path / "out.1").write_text(
        "# columns\n" + "\n".join(" ".join(map(str, row))
                                   for row in rows) + "\n")
    config = {
        "model": str(Path(__file__).parents[2] / "examples" / "prem.200"),
        "output": {"prefix": "out", "quantity": "displacement"},
        "sources": [{"name": "1", "latitude": 0.0, "longitude": 0.0,
                     "depth_km": 1.0}],
        "snapshots": {"times_s": [1], "depths_km": [3],
                      "distances_deg": [0, 30], "azimuths_deg": [0, 180]},
    }
    path = tmp_path / "input.json"
    path.write_text(json.dumps(config))
    run = SnapshotRun(path)
    np.testing.assert_array_equal(run.latitude,
                                  np.array([[11.0, 12.0], [31.0, 32.0]]))
    np.testing.assert_array_equal(run.longitude,
                                  np.array([[21.0, 22.0], [41.0, 42.0]]))


def test_snapshot_sections_reject_skewed_azimuths_and_gap_duplicate_origin():
    model = Path(__file__).parents[2] / "examples" / "prem.200"
    run = SimpleNamespace(model=model, distances=np.array([0.0, 30.0, 60.0]),
                          azimuths=np.array([0.0, 180.0]),
                          depths=np.array([100.0, 500.0, 1000.0]))
    field = np.ones((3, 2, 3))
    with pytest.raises(ValueError, match="opposite azimuths"):
        _section_slice(SimpleNamespace(**{**vars(run),
                                         "azimuths": np.array([0.0, 150.0])}),
                       field, 6371.0,
                       Shown("up", dimensions=DIMENSIONS["displacement"]))
    cut = _section_slice(run, field, 6371.0,
                         Shown("up", dimensions=DIMENSIONS["displacement"]))
    assert np.all(np.diff(cut.distance) > 0)
    assert np.unique(cut.path, axis=0).shape[0] == cut.path.shape[0]
    assert len(cut.layers) > 1
    assert any(np.isnan(layer.values).any() for layer in cut.layers)


def test_snapshot_render_keeps_solver_side_samples_at_prem_interfaces():
    model = Path(__file__).parents[2] / "examples" / "prem.200"

    class Run:
        times = np.array([150.0, 300.0, 600.0, 900.0, 1200.0, 1500.0])
        depths = np.array([2890.0, 2891.0, 2892.0,
                           5148.5, 5149.5, 5150.5])
        azimuths = np.array([0.0])
        distances = np.array([10.0, 30.0, 60.0])
        quantity = "displacement"
        source = "1"
        hypocentre = {"depth_km": 29.5}

        def field(self, index, component):
            return np.arange(18.0).reshape(6, 1, 3) + 1

    Run.model = model
    run = Run()
    field = run.field(0, "up")
    cut = _section_slice(
        run, field, 6371.0,
        Shown("up", dimensions=DIMENSIONS["displacement"]))
    # At the PREM CMB, the outer side is solid (shell 2); at the ICB,
    # the inner side is solid (shell 0), matching RadialSolver::SolidSide.
    assert {layer.index for layer in cut.layers} == {0, 1, 2}
    radii = {layer.index: layer.radius for layer in cut.layers}
    assert np.max(radii[0]) <= 1221.5
    assert np.min(radii[1]) > 1221.5 and np.max(radii[1]) < 3480.0
    assert np.min(radii[2]) >= 3480.0
    np.testing.assert_array_equal(
        np.sort(np.concatenate([layer.values.ravel()
                                for layer in cut.layers])),
        np.sort(field.reshape(-1)))
    fig = snapshots(run, component="up")
    fig.canvas.draw()
    assert len(fig.axes) == 12
    caption = fig.texts[1]
    assert "solid side" in caption.get_text()
    caption_box = caption.get_window_extent()
    assert caption_box.y1 < fig.texts[0].get_window_extent().y0
    assert all(caption_box.y0 > ax.title.get_window_extent().y1
               for ax in fig.axes if ax.get_title())


def test_snapshot_renders_when_all_depths_are_prem_interfaces():
    model = Path(__file__).parents[2] / "examples" / "prem.200"

    class Run:
        times = np.array([1.0])
        depths = np.array([2891.0, 5149.5])
        azimuths = np.array([0.0])
        distances = np.array([10.0, 30.0, 60.0])
        quantity = "displacement"
        source = "1"
        hypocentre = {"depth_km": 29.5}

        def field(self, index, component):
            return np.arange(6.0).reshape(2, 1, 3) + 1

    Run.model = model
    run = Run()
    field = run.field(0, "up")
    cut = _section_slice(
        run, field, 6371.0,
        Shown("up", dimensions=DIMENSIONS["displacement"]))
    assert {layer.index for layer in cut.layers} == {0, 2}
    np.testing.assert_array_equal(
        np.sort(np.concatenate([layer.values.ravel()
                                for layer in cut.layers])),
        np.sort(field.reshape(-1)))
    fig = snapshots(run, component="up")
    fig.canvas.draw()
    assert len(fig.axes) == 2
    caption_box = fig.texts[1].get_window_extent()
    assert caption_box.y1 < fig.texts[0].get_window_extent().y0
    assert all(caption_box.y0 > ax.title.get_window_extent().y1
               for ax in fig.axes if ax.get_title())


def test_snapshot_sections_show_units_limits_and_saturation_for_zero_fields():
    model = Path(__file__).parents[2] / "examples" / "prem.200"

    class Run:
        times = np.array([1.0])
        depths = np.array([100.0, 500.0, 1000.0])
        azimuths = np.array([0.0, 180.0])
        distances = np.array([0.0, 30.0, 60.0])
        quantity = "displacement"
        source = "1"
        hypocentre = {"depth_km": 29.5, "latitude": 36.87,
                      "longitude": 35.58}

        def field(self, index, component):
            return np.zeros((3, 2, 3))

    Run.model = model
    fig = snapshots(Run(), component="up")
    assert len(fig.axes) == 2  # one narrow panel and its quantitative bar
    assert "up displacement [m]" in fig._suptitle.get_text()
    caption = fig.texts[1]
    assert "interpolates colours" in caption.get_text()
    fig.canvas.draw()
    assert caption.get_window_extent().x0 >= 0
    assert caption.get_window_extent().x1 <= fig.bbox.width
    stats = fig.axes[0].texts[0].get_text()
    assert "slice max 0" in stats
    assert "saturated 0.00%" in stats
    assert len(stats.splitlines()) == 3
    assert all(len(line) <= 20 for line in stats.splitlines())


def test_snapshot_map_uses_noninterpolating_cells_without_wrap_warning():
    cartopy = pytest.importorskip("cartopy.crs")
    import warnings
    import matplotlib.pyplot as plt

    lat = np.array([[-60.0, -55.0], [10.0, 15.0], [70.0, 75.0]])
    lon = np.array([[0.0, 60.0], [120.0, 180.0], [240.0, 300.0]])
    run = SimpleNamespace(latitude=lat, longitude=lon,
                          hypocentre={"longitude": 0.0})
    fig = plt.figure(figsize=(4.2, 5.2))
    ax = fig.add_subplot(1, 1, 1, projection=cartopy.PlateCarree())
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        mesh = _draw_map(ax, run, np.arange(6).reshape(1, 3, 2),
                         cartopy.PlateCarree(), 5, "RdBu", coastlines=False)
        fig.canvas.draw()
    # Matplotlib converts nearest-neighbour sample assignment to flat cell
    # shading internally after deriving cell edges.
    assert mesh._shading == "flat"
    assert not any("wrapped coordinates" in str(w.message).lower()
                   for w in caught)
    plt.close(fig)


def test_parameter_unit_covers_all_solver_kernel_families():
    from yspec.plot.kernels import _kernel_unit, _parameter_unit

    families = {
        "love": {"rho", "A", "C", "F", "L", "N", "qA", "qF", "qmu"},
        "moduli": {"rho", "kappa", "mu", "qkappa", "qmu"},
        "velocities": {"rho", "vp", "vs", "qkappa", "qmu"},
        "anisotropic": {"rho", "vpv", "vph", "vsv", "vsh", "eta",
                        "qkappa", "qmu"},
    }
    expected = {
        "rho": "kg m$^{-3}$",
        "vp": "m s$^{-1}$", "vs": "m s$^{-1}$",
        "vpv": "m s$^{-1}$", "vph": "m s$^{-1}$",
        "vsv": "m s$^{-1}$", "vsh": "m s$^{-1}$",
        "A": "Pa", "C": "Pa", "F": "Pa", "L": "Pa", "N": "Pa",
        "kappa": "Pa", "mu": "Pa", "eta": "dimensionless",
        "qA": "dimensionless", "qF": "dimensionless",
        "qkappa": "dimensionless", "qmu": "dimensionless",
    }
    for names in families.values():
        for name in names:
            assert _parameter_unit(name) == expected[name]
            assert _kernel_unit(name, relative=False).endswith(
                f"per ({expected[name]})")
    with pytest.raises(ValueError, match="unsupported kernel parameter"):
        _parameter_unit("not-a-solver-parameter")


def test_radial_profiles_show_saved_nodes_only_and_split_at_boundaries():
    from yspec.plot.kernels import _profile

    radius = np.array([1.0, 2.0, 4.0, 5.0])
    values = np.array([10.0, 20.0, 40.0, 50.0])
    profile = _profile(radius, values, np.array([3.0]), "vp", 0.0)
    assert len(profile.layers) == 2
    np.testing.assert_array_equal(profile.layers[0].radius, radius[:2])
    np.testing.assert_array_equal(profile.layers[0].values, values[:2])
    np.testing.assert_array_equal(profile.layers[1].radius, radius[2:])
    np.testing.assert_array_equal(profile.layers[1].values, values[2:])
    assert all(3.0 not in layer.radius for layer in profile.layers)


def test_volume_adapter_converts_geocentric_grid_to_geographic_latitude():
    from planetmodel.sampling import AngularGrid
    from planetmodel.slices import SphereSlice

    colatitude = np.radians([30.0, 90.0, 150.0])
    grid = AngularGrid(colatitude, np.radians([0.0, 90.0]), kind="custom")
    cut = SphereSlice(1.0, grid, np.ones((3, 2)), np.zeros((3, 2, 3)),
                      np.zeros((3, 2, 3)), Shown("K"))
    geographic = _geographic_slice(cut)
    actual = 90.0 - np.degrees(geographic.grid.colatitudes)
    expected = np.degrees(np.arctan(
        np.tan(np.radians([60.0, 0.0, -60.0])) / 0.993277))
    np.testing.assert_allclose(actual, expected)
    np.testing.assert_array_equal(geographic.values, cut.values)


def test_volume_overlay_follows_the_field_great_circle_plane():
    from yspec.plot.volume import _geographic_great_circle, _unit

    source = (36.87, 35.58)
    receiver = (12.0, 80.0)
    longitude, latitude = _geographic_great_circle(source, receiver, 51)
    assert len(longitude) == len(latitude) == 51
    assert (latitude[0], longitude[0]) == pytest.approx(source)
    assert (latitude[-1], longitude[-1]) == pytest.approx(receiver)
    s, r = _unit(*source), _unit(*receiver)
    normal = np.cross(s, r)
    points = np.stack([_unit(lat, lon) for lat, lon in
                       zip(latitude, longitude)])
    np.testing.assert_allclose(points @ normal, 0.0, atol=2e-14)
    seam_lon, _ = _geographic_great_circle((10.0, 170.0), (20.0, -170.0), 51)
    assert np.max(np.abs(np.diff(seam_lon))) < 2.0



def test_clipping_statistics_stay_short_for_large_spatial_values():
    from yspec.plot._labels import clipping_stats

    label = clipping_stats(1.2345e100, 9.8765e120, 100.0)
    lines = label.splitlines()
    assert len(lines) == 3
    assert all(len(line) <= 20 for line in lines)
    assert lines[0].startswith("limit")
    assert lines[1].startswith("slice max")
    assert lines[2] == "saturated 100.00%"


def test_radial_kernel_footer_fits_narrow_single_panel_and_discloses_scales():
    import matplotlib.pyplot as plt
    from yspec.plot.kernels import radial_kernels

    model = Path(__file__).parents[2] / "examples" / "prem.200"
    radius = np.linspace(10.0, 6371.0, 80)
    run = SimpleNamespace(
        parameters=["vs"], model=model, radius=radius,
        kernels={"vs": np.sin(radius)}, relative=True,
        time=np.linspace(0.0, 2400.0, 1201),
        seismogram=np.sin(np.linspace(0.0, 4.0 * np.pi, 1201)),
        taper=np.zeros(1201), component="T", quantity="velocity",
        description="radial kernels of shear-packet", source_radius=6371.0-29.5,
        at_source={"vs": 0.0}, surface=6371.0)
    fig = radial_kernels(run)
    fig.canvas.draw()
    footer = fig.texts[0]
    assert "own horizontal SI scale" in footer.get_text()
    assert "do not establish significance" in footer.get_text()
    bounds = footer.get_window_extent()
    assert bounds.x0 >= 0 and bounds.x1 <= fig.bbox.width
    assert bounds.y0 >= 0 and bounds.y1 <= fig.bbox.height
    plt.close(fig)


def test_kernel_section_titles_disclose_the_displayed_radial_extent():
    from yspec.plot.sections import _extent_label

    assert _extent_label(SimpleNamespace(name="shear-packet")) == "r ≥ 3000 km"
    assert _extent_label(SimpleNamespace(name="PKIKP")) == "full radius"
