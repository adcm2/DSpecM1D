#!/usr/bin/env python3
"""Writes inputs with cross-correlation windows centred on PREM TauP ray times.

The phase names are labels for the targeted arrivals and CLI selection, not a
claim that a finite-frequency window isolates a single phase. The generated
windows can contain depth phases or other branches; inspect the reported
co-arrivals and validate each filtered waveform before interpreting its
kernel as a phase-specific sensitivity. The high-frequency gallery band and
the longer-period diffracted band are specialized packet demonstrations, not
the Bolivia comparison setup used by the first-seismogram guide.

    python make_input.py    # in an environment with yspec-seismo[obspy]

The inputs beside it were written so, and the example runs without it; this
is how to change the phases, the distances or the bands.
"""

import json
import pathlib

from obspy.taup import TauPyModel

AZIMUTH = 70.0  # degrees, clockwise from north at the source
DEPTH = 29.5    # km
RECEIVER_DEPTH = 3.0  # km; use the depth specified by the YSpec input

# Each run: its receiver's distance (degrees), the phases measured and on
# which components, the dominant period of its band (s), and what differs
# between the runs of the input.
RUNS = {
    "gallery": {
        "distance": 60.0,
        "phases": [("P", "Z"), ("PP", "Z"), ("S", "T"), ("ScS", "T"),
                   ("SS", "T")],
        "period": 25.0,
        "input": {
            # Enough degrees for the fundamental surface waves up to the top
            # of the band, l ~ 2 pi f a / c ~ 410 at 40 mHz: cut off below
            # that, they would be cut off in frequency where their branches
            # pass the last degree, and their sinc tails run through the whole
            # record, windows and all. Nothing below the band is computed.
            "degrees": {"max": 450},
            "frequencies_mHz": {"min": 7, "max": 40},
            # A time step at which the samples are just within a power of
            # two, for the fewest frequencies the damping allows.
            "time_series": {"length_min": 30, "time_step_s": 1.2},
            "filter_mHz": [8, 12, 35, 40],
            "sections": {"margin_deg": 15, "half_width_deg": 15,
                         "spacing_deg": 0.25, "radial_spacing_km": 30},
        },
    },
    # Diffracted waves lose their high frequencies along the core-mantle
    # boundary, and a longer period suits them.
    "diffracted": {
        "distance": 110.0,
        "phases": [("Pdiff", "Z"), ("Sdiff", "T")],
        "period": 40.0,
        "input": {
            "degrees": {"max": 300},
            "frequencies_mHz": {"min": 4, "max": 25},
            # The surface waves arrive after the record, at 2800 to 3500 s:
            # a time step at which the transform's period holds them, so
            # that they do not wrap round onto Pdiff.
            "time_series": {"length_min": 30, "time_step_s": 1.8},
            "filter_mHz": [5, 8, 20, 25],
            "sections": {"margin_deg": 15, "half_width_deg": 15,
                         "spacing_deg": 0.35, "radial_spacing_km": 40},
        },
    },
}

# The arrivals that may share a window, which are printed.
OTHERS = ["P", "pP", "sP", "PcP", "pPcP", "PP", "pPP", "sPP", "PPP", "PS",
          "SP", "S", "sS", "ScS", "sScS", "SS", "sSS", "SSS", "SKS", "sSKS",
          "SKKS", "ScP", "PcS", "Pdiff", "pPdiff", "sPdiff", "Sdiff",
          "sSdiff", "PKiKP", "PKIKP"]


def measure(model, distance, phase, component, period):
    """The measurement of a phase, and what else arrives in its window."""
    arrivals = model.get_travel_times(
        source_depth_in_km=DEPTH, receiver_depth_in_km=RECEIVER_DEPTH,
        distance_in_degree=distance, phase_list=[phase])
    time = min(a.time for a in arrivals)
    start = round(time - period)
    end = round(time + 1.5 * period)
    others = sorted({(round(a.time), a.name) for a in model.get_travel_times(
        source_depth_in_km=DEPTH, receiver_depth_in_km=RECEIVER_DEPTH,
        distance_in_degree=distance,
        phase_list=OTHERS) if start <= a.time <= end and a.name != phase})
    print(f"{phase:5s} on {component} at {distance:g} degrees: arrives at "
          f"{time:7.1f} s, window {start} to {end} s; with "
          + (", ".join(f"{name} {t} s" for t, name in others)
             or "nothing else"))
    return {
        "name": phase,
        "component": component,
        "traveltime": {"window_s": [start, end], "taper_s": period / 2},
    }


def main():
    model = TauPyModel("prem")
    for name, run in RUNS.items():
        settings = dict(run["input"])
        sections = settings.pop("sections")
        content = {
            "model": "../prem.200",
            "gravity": "full",
            "output": {"prefix": name, "quantity": "velocity"},
            **settings,
            "sources": [
                {"depth_km": DEPTH, "latitude": 36.87, "longitude": 35.58,
                 "moment_tensor_Nm": {"rr": 4.22e17, "rt": 3.08e17,
                                      "rp": 6.99e17, "tt": -2.976e18,
                                      "tp": 6.54e17, "pp": 2.55e18}}
            ],
            "receiver": {"depth_km": RECEIVER_DEPTH,
                         "distance_deg": run["distance"],
                         "azimuth_deg": AZIMUTH},
            "measurements": [measure(model, run["distance"], phase,
                                     component, run["period"])
                             for phase, component in run["phases"]],
            "kernels": {
                "parameters": "velocities",
                "scaling": "relative",
                "mesh": {"order": 4, "max_element_km": 100},
                "sections": sections,
            },
        }
        path = pathlib.Path(__file__).with_name(f"{name}.json")
        path.write_text(json.dumps(content, indent=2) + "\n")
        print(f"wrote {path.name}")


if __name__ == "__main__":
    main()
