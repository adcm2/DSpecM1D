"""What the tests share: a small calculation, described both through the
Python interface and as an input file, which must give the same results."""

import json
from importlib import resources

import pytest

import yspec

MOMENT = dict(rr=4.22e17, rt=3.08e17, rp=6.99e17, tt=-2.976e18, tp=6.54e17,
              pp=2.55e18)


@pytest.fixture(scope="session")
def model():
    return yspec.Model.prem()


@pytest.fixture(scope="session")
def deck():
    with resources.as_file(resources.files("yspec") / "models" /
                           "prem.200") as path:
        return str(path)


@pytest.fixture
def source():
    return yspec.Source(latitude=36.87, longitude=35.58, depth_km=29.5,
                        moment_tensor=MOMENT)


@pytest.fixture
def calculation():
    return yspec.Calculation(degrees=60, frequencies_mhz=(0.5, 5),
                             length_min=20, time_step_s=10,
                             filter_mhz=(0.5, 1, 4, 5), quantity="velocity")


def shared_members(deck, prefix):
    """The members of an input file for the calculation above."""
    return {
        "model": deck,
        "gravity": "full",
        "degrees": {"max": 60},
        "frequencies_mHz": {"min": 0.5, "max": 5},
        "time_series": {"length_min": 20, "time_step_s": 10},
        "filter_mHz": [0.5, 1, 4, 5],
        "sources": [{"depth_km": 29.5, "latitude": 36.87, "longitude": 35.58,
                     "moment_tensor_Nm": MOMENT}],
        "output": {"prefix": prefix, "quantity": "velocity"},
    }


def write(path, content):
    path.write_text(json.dumps(content))
    return str(path)
