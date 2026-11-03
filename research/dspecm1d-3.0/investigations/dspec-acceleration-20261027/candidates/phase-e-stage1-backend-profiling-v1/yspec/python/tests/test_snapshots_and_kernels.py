import numpy as np
import pytest
from conftest import shared_members, write

import yspec
from yspec import _core


def test_snapshots_are_those_of_an_input_file(model, deck, source,
                                              calculation, tmp_path):
    times, depths = [300.0, 600.0], [10.0, 1000.0]
    distances, azimuths = [0.0, 30.0, 60.0], [20.0, 200.0]
    ours = yspec.snapshots(model, calculation, source, times, depths,
                           distances, azimuths)
    members = shared_members(deck, str(tmp_path / "snap"))
    members["snapshots"] = {"times_s": times, "depths_km": depths,
                            "distances_deg": distances,
                            "azimuths_deg": azimuths}
    theirs = _core.run_snapshots(
        _core.read_snapshot_input(write(tmp_path / "in.json", members)))
    assert ours.data.shape == (1, 2, 2, 2, 3, 3)
    np.testing.assert_allclose(ours.data, theirs.data, rtol=0,
                               atol=1e-12 * np.abs(theirs.data).max())


@pytest.fixture
def receiver():
    return yspec.Receiver(36.87, 105.0, depth_km=3)


def test_kernels_are_those_of_an_input_file(model, deck, source, calculation,
                                            receiver, tmp_path):
    measured = [yspec.Traveltime("T", (600, 1000), 100, name="S"),
                yspec.Traveltime("Z", (300, 700), 100, name="P")]
    measurements = yspec.traveltimes(model, calculation, source, receiver,
                                     measured)
    ours = yspec.radial_kernels(model, calculation, measurements,
                                scaling="relative", max_element_km=200)

    members = shared_members(deck, str(tmp_path / "k"))
    members["receiver"] = {"depth_km": 3, "location": [36.87, 105.0]}
    members["measurements"] = [
        {"name": "S", "component": "T",
         "traveltime": {"window_s": [600, 1000], "taper_s": 100}},
        {"name": "P", "component": "Z",
         "traveltime": {"window_s": [300, 700], "taper_s": 100}}]
    members["kernels"] = {"scaling": "relative",
                          "mesh": {"max_element_km": 200}}
    theirs = _core.run_kernels(
        _core.read_kernel_input(write(tmp_path / "in.json", members)))

    for kernels, their in zip(ours, theirs.measurements):
        assert kernels.names == list(their.named.names)
        np.testing.assert_allclose(kernels.measurement.adjoint_source,
                                   their.measurement.adjoint_source)
        for name, values in zip(their.named.names, their.named.values):
            np.testing.assert_allclose(kernels[name], values, rtol=0,
                                       atol=1e-12 * np.abs(values).max())
    assert (tmp_path / "k.S.kernels").exists()


def test_sampled_kernels_integrate_to_the_radial_kernels(model, source,
                                                         calculation,
                                                         receiver):
    measurement = yspec.traveltime(model, calculation, source, receiver, "T",
                                   (600, 1000), 100)
    radial = yspec.radial_kernels(model, calculation, measurement,
                                  max_element_km=500)
    sampling = yspec.volume(model, degree=8, max_element_km=500)
    sampled = yspec.sampled_kernels(model, calculation, measurement, sampling)
    assert sampled.names == ["vp", "vs", "qkappa", "qmu"]
    assert sampled.values.shape == (4, sampling._core.node_count,
                                    sampling._core.place_count)
    assert np.all(np.isfinite(sampled.values))
    assert radial["vs"].shape == radial.radius_km.shape

    sections = yspec.sampled_kernels(
        model, calculation, [measurement],
        [yspec.path_section(model, source, receiver, spacing_deg=5,
                            radial_spacing_km=500),
         yspec.path_maps(model, source, receiver, [100.0], spacing_deg=5)])
    path, maps = sections[0]
    assert np.all(np.diff(path.along_deg) > 0)
    assert np.all(np.abs(maps.depth_km - 100) < 1e-6)
    assert np.all(np.abs(maps.latitude) <= 90)
