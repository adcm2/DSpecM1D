"""A source-built Python calculation matching the documented Bolivia case."""

import numpy as np
import yspec

yspec.set_threads(2)
model = yspec.Model.prem()
calculation = yspec.Calculation(
    degrees=300,
    frequencies_mhz=(1.5, 21.0),
    length_min=300,
    time_step_s=1.0,
    filter_mhz=(1.9, 2.0, 20.0, 20.1),
    gravity="full",
    attenuation=False,
    quantity="displacement",
    tolerance=1e-5,
)
source = yspec.Source(
    latitude=-13.909919806613164,
    longitude=-67.25,
    depth_km=647.1,
    moment_tensor={
        "rr": -7.59e20,
        "rt": -2.5e21,
        "rp": 4.2e19,
        "tt": 7.75e20,
        "tp": -2.48e20,
        "pp": -1.6e19,
    },
    name="bolivia",
)
receiver = yspec.Receiver(
    latitude=80.06588632614357, longitude=0.0, depth_km=3.0)
result = yspec.seismograms(model, calculation, source, receiver, progress=True)

assert result.data.shape == (1, 1, 3, 18001)
assert result.components == ("Z", "N", "E")
assert result.time_s.shape == (18001,)
assert result.time_s[0] == 0 and result.time_s[-1] == 18000
assert np.isfinite(result.data).all()
assert np.max(np.abs(result.trace("Z"))) > 0
print("data shape [source, receiver, component, sample]:", result.data.shape)
print("components:", result.components)
print("time_s range:", result.time_s[0], result.time_s[-1])
print("vertical peak displacement (m):", np.max(np.abs(result.trace("Z"))))
