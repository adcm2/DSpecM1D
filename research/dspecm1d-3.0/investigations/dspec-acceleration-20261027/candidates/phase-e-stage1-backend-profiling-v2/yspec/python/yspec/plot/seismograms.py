"""Pictures of seismograms: record sections, traces and comparisons."""

import numpy as np

#: A colour for each component.
COLOURS = {"Z": "#1f4e79", "N": "#7f2704", "E": "#2d6a4f",
           "R": "#7f2704", "T": "#2d6a4f"}


def record_section(run, component="Z", *, common=False, gain=1.0):
    """One component at every receiver against epicentral distance, each
    trace scaled to its peak (or all alike, with `common`), amplified by
    `gain` and clipped where it would reach the next."""
    import matplotlib.pyplot as plt

    receivers = sorted(run.receivers, key=lambda r: r.distance)
    traces = [r.trace(component) for r in receivers]
    distances = np.array([r.distance for r in receivers])
    spacing = np.median(np.diff(distances)) if len(receivers) > 1 else 1.0
    peak = max(np.abs(u).max() for _, u in traces) or 1.0
    fig, ax = plt.subplots(figsize=(8, 9))
    for distance, (t, u) in zip(distances, traces):
        scale = peak if common else (np.abs(u).max() or 1.0)
        u = np.clip(gain * u / scale, -1, 1)
        ax.plot(t / 60, distance + 0.9 * spacing * u, lw=0.6, color="k")
    t = traces[0][0]
    ax.set_xlim(t[0] / 60, t[-1] / 60)
    ax.set_ylim(distances[0] - spacing, distances[-1] + spacing)
    ax.set_xlabel("time after the origin (min)")
    ax.set_ylabel("epicentral distance (degrees)")
    scaling = "common scale" if common else "each trace scaled to its peak"
    if gain != 1:
        scaling += f", times {gain:g} and clipped"
    ax.set_title(f"{component}, {run.quantity}, source {run.source} "
                 f"({scaling})")
    return fig


def traces(run, numbers=None, components="ZRT", *, scale="common"):
    """The components at some receivers (by default the first eight), a row
    for each. By default the components in each receiver row share a physical
    amplitude scale; ``scale='independent'`` gives every component its own
    scale and is identified in the figure title."""
    import matplotlib.pyplot as plt

    if scale not in ("common", "independent"):
        raise ValueError("scale must be 'common' or 'independent'")
    numbers = numbers or list(range(1, min(len(run.receivers), 8) + 1))
    receivers = [run.receiver(k) for k in numbers]
    components = list(components)
    fig, axes = plt.subplots(len(receivers), len(components), sharex=True,
                             squeeze=False,
                             figsize=(3.4 * len(components) + 1.2,
                                      1.5 * len(receivers) + 1.2))
    for i, receiver in enumerate(receivers):
        row_traces = [receiver.trace(component) for component in components]
        row_peak = max((np.max(np.abs(u)) for _, u in row_traces), default=0.0)
        for j, component in enumerate(components):
            t, u = row_traces[j]
            ax = axes[i][j]
            ax.plot(t / 60, u, lw=0.7, color=COLOURS[component])
            ax.set_xlim(t[0] / 60, t[-1] / 60)
            if scale == "common":
                limit = row_peak or 1.0
                ax.set_ylim(-1.05 * limit, 1.05 * limit)
            ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 3))
            ax.tick_params(labelsize=8)
            ax.yaxis.get_offset_text().set_fontsize(8)
            if i == 0:
                ax.set_title(component)
            if j == 0:
                ax.set_ylabel(f"{receiver.number}: {receiver.depth:g} km\n"
                              f"{receiver.distance:.1f}°, {receiver.units}",
                              fontsize=8)
    for ax in axes[-1]:
        ax.set_xlabel("time after the origin (min)")
    scale_note = ("components share a scale within each receiver row"
                  if scale == "common" else
                  "each component scaled independently; amplitudes are not comparable")
    fig.suptitle(f"{run.quantity}, source {run.source} ({scale_note})")
    fig.tight_layout()
    return fig


def compare(runs, labels=None, *, receiver=1, component="Z"):
    """One receiver in several runs, overlaid, and their differences from
    the last."""
    import matplotlib.pyplot as plt

    if len(runs) < 2:
        raise ValueError("compare needs at least two runs")
    labels = labels or [run.path.stem for run in runs]
    if len(labels) != len(runs):
        raise ValueError("there must be a label for each run")
    receivers = [run.receiver(receiver) for run in runs]
    traces = [r.trace(component) for r in receivers]
    reference_run, reference_receiver = runs[-1], receivers[-1]
    reference_time = traces[-1][0]
    for run, rec, (time, _) in zip(runs[:-1], receivers[:-1], traces[:-1]):
        if run.quantity != reference_run.quantity:
            raise ValueError("runs have different physical quantities")
        if rec.units != reference_receiver.units:
            raise ValueError("receivers have different physical units")
        if not np.array_equal(time, reference_time):
            raise ValueError("runs have different physical time grids")
        geometry = (rec.latitude, rec.longitude, rec.depth, rec.distance,
                    rec.azimuth, rec.back_azimuth)
        reference_geometry = (reference_receiver.latitude,
                              reference_receiver.longitude,
                              reference_receiver.depth,
                              reference_receiver.distance,
                              reference_receiver.azimuth,
                              reference_receiver.back_azimuth)
        if not np.allclose(geometry, reference_geometry, rtol=0, atol=1e-4):
            raise ValueError("runs have different receiver geometry")
    fig, (top, bottom) = plt.subplots(2, 1, sharex=True, figsize=(9, 6))
    for (t, u), label in zip(traces, labels):
        top.plot(t / 60, u, lw=0.8, label=label)
    t, reference = traces[-1]
    for (t, u), label in zip(traces[:-1], labels[:-1]):
        bottom.plot(t / 60, u - reference, lw=0.8,
                    label=f"{label} $-$ {labels[-1]}")
    for ax in (top, bottom):
        ax.legend(fontsize=8, loc="upper right")
        ax.set_ylabel(receivers[-1].units)
    bottom.set_xlim(t[0] / 60, t[-1] / 60)
    bottom.set_xlabel("time after the origin (min)")
    last = receivers[-1]
    top.set_title(f"{component}, {runs[-1].quantity} at "
                  f"{last.distance:.1f}°, {last.depth:g} km deep")
    fig.tight_layout()
    return fig
