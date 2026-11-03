import numpy as np
import pytest
from conftest import shared_members, write

import yspec
from yspec import _core


def test_seismograms_are_those_of_an_input_file(model, deck, source,
                                                calculation, tmp_path):
    receivers = [yspec.Receiver(10.0, 80.0, depth_km=3),
                 *yspec.great_circle(source, azimuth_deg=60,
                                     distances_deg=[40])]
    ours = yspec.seismograms(model, calculation, source, receivers)

    members = shared_members(deck, str(tmp_path / "out"))
    members["receivers"] = [
        {"depth_km": 3, "locations": [[10.0, 80.0]]},
        {"depth_km": 0, "profile": {"from": [36.87, 35.58],
                                    "azimuth_deg": 60,
                                    "distances_deg": [40, 40], "count": 1}},
    ]
    theirs = _core.run_seismograms(
        _core.read_seismogram_input(write(tmp_path / "in.json", members)))

    assert ours.data.shape == (1, 2, 3, len(ours.time_s))
    np.testing.assert_allclose(ours.data, theirs.data, rtol=0,
                               atol=1e-12 * np.abs(theirs.data).max())
    assert (tmp_path / "out.1").exists()


def test_rotation_keeps_the_horizontal_motion(model, source, calculation):
    receivers = yspec.great_circle(source, azimuth_deg=30,
                                   distances_deg=[20, 70])
    seismograms = yspec.seismograms(model, calculation, source, receivers)
    turned = seismograms.rotated()
    horizontal = np.hypot(seismograms.data[:, :, 1], seismograms.data[:, :, 2])
    np.testing.assert_allclose(np.hypot(turned[:, :, 1], turned[:, :, 2]),
                               horizontal, atol=1e-12 * horizontal.max())
    np.testing.assert_array_equal(turned[:, :, 0], seismograms.data[:, :, 0])
    np.testing.assert_allclose(seismograms.trace("T", receiver=1),
                               turned[0, 1, 2])


def test_strains_have_six_components(model, source, calculation):
    strains = yspec.strains(model, calculation, source,
                            yspec.Receiver(0.0, 60.0))
    assert strains.data.shape == (1, 1, 6, len(strains.time_s))
    assert np.all(np.isfinite(strains.data))


def test_progress_is_reported_and_can_stop_a_calculation(model, source,
                                                         calculation):
    reported = []
    yspec.seismograms(model, calculation, source, yspec.Receiver(0.0, 60.0),
                      progress=lambda done, total: reported.append(
                          (done, total)))
    total = reported[0][1]
    assert sorted(done for done, _ in reported) == list(range(1, total + 1))

    def stop(done, total):
        raise KeyError("stopped")

    with pytest.raises(KeyError, match="stopped"):
        yspec.seismograms(model, calculation, source,
                          yspec.Receiver(0.0, 60.0), progress=stop)


def test_faulty_calculations_are_refused(model, source):
    receiver = yspec.Receiver(0.0, 60.0)
    with pytest.raises(ValueError, match="frequencies"):
        yspec.Calculation(degrees=10, frequencies_mhz=(5, 1),
                          length_min=10)._core()
    with pytest.raises(ValueError, match="gravity"):
        yspec.Calculation(degrees=10, frequencies_mhz=(1, 5), length_min=10,
                          gravity="some")._core()
    with pytest.raises(ValueError, match="latitude"):
        yspec.Receiver(100.0, 0.0)._core()
    with pytest.raises(ValueError):
        # Below the centre of the earth.
        deep = yspec.Source(0.0, 0.0, 7000.0, moment_tensor=dict(
            rr=1, rt=0, rp=0, tt=0, tp=0, pp=0))
        yspec.seismograms(model, yspec.Calculation(10, (1, 5), 10), deep,
                          receiver)


def test_instruments_come_from_poles_and_zeros(model, source, calculation):
    instrument = yspec.Instrument(zeros=[0, 0], poles=[-0.01 + 0.01j,
                                                       -0.01 - 0.01j],
                                  constant=2.0)
    plain = yspec.seismograms(model, calculation, source,
                              yspec.Receiver(0.0, 60.0))
    recorded = yspec.seismograms(
        model, calculation, source,
        yspec.Receiver(0.0, 60.0, instrument=instrument))
    assert not np.allclose(plain.data, recorded.data)
    assert np.all(np.isfinite(recorded.data))
