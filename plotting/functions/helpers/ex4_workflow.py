import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patheffects

from .ex1_workflow import print_workflow_context
from .plot_common import OKABE_ITO, loadtxt_or_exit


def output_to_screen_code_diff_ex4(path_waveforms, path_travel_times):
    print_workflow_context(
        data_file=path_waveforms,
        primary_label="Travel-time file",
        primary_range=path_travel_times,
    )


def plot_ex4(path_waveforms, path_travel_times, show=True):
    d1 = loadtxt_or_exit(path_waveforms)
    try:
        data = np.genfromtxt(path_travel_times, delimiter=";", names=True, dtype=None, encoding="utf-8")
    except FileNotFoundError:
        print(f"Error: Data file {path_travel_times} not found")
        raise SystemExit(1)

    fig, ax = plt.subplot_mosaic([["left", "right"]], figsize=(14, 10), constrained_layout=True)
    ax_p, ax_s = ax["left"], ax["right"]

    medium_size = 20
    bigger_size = 24

    colors = {
        "waveform": "#111111",
        "phase": OKABE_ITO["vermillion"],
    }

    for i in range(0, 91):
        maxvalp = np.max(np.abs(d1[:, 3 * i + 1])) / 2
        if maxvalp > 0:
            ax_p.plot(d1[:, 3 * i + 1] / maxvalp + 2 * i, d1[:, 0], color=colors["waveform"], linewidth=1.1)

    for i in range(0, 91):
        maxvals = np.max(np.abs(d1[:, 3 * i + 2])) / 2
        if maxvals > 0:
            ax_s.plot(d1[:, 3 * i + 2] / maxvals + 2 * i, d1[:, 0], color=colors["waveform"], linewidth=1.1)

    dist_min, dist_max = 0, 180
    time_min, time_max = 0, 3600

    phases_s = [
        "s", "S", "ScS", "sS", "ScsScS", "sScsScS", "ScsScsScS", "sScsScsScS",
        "sScS", "Sdiff", "SS", "sSv410sScS", "sSv660sScS", "SKS", "SSS", "SSSS", "SSSSS",
    ]
    phases_p = ["p", "P", "PcP", "Pdiff", "pP", "PKP", "PP", "pPP", "PKIKP", "pPPP", "PPP", "PcPPcP"]

    text_style = [patheffects.withStroke(linewidth=3, foreground="white")]

    def plot_phases(axes, phase_list, data_struct):
        for phase_name in phase_list:
            mask = data_struct["Phase"] == phase_name
            dist = data_struct["Distance_deg"][mask]
            time = data_struct["TravelTime_s"][mask]

            if len(dist) > 0:
                axes.plot(dist, time, color=colors["phase"], label=phase_name, linewidth=2.0)
                in_view = (
                    (dist >= dist_min)
                    & (dist <= dist_max)
                    & (time >= time_min)
                    & (time <= time_max)
                )

                if np.any(in_view):
                    idx_candidates = np.where(in_view)[0]
                    idx = idx_candidates[-1]
                    txt = axes.text(
                        dist[idx],
                        time[idx],
                        f" {phase_name}",
                        color=colors["phase"],
                        fontsize=14,
                        fontweight="bold",
                        va="center",
                        ha="left",
                    )
                    txt.set_path_effects(text_style)

    plot_phases(ax_p, phases_p, data)
    plot_phases(ax_s, phases_s, data)

    xticks_val = [0, 45, 90, 135, 180]
    yticks_val = [0, 600, 1200, 1800, 2400, 3000, 3600]
    spine_width = 2.0

    for target_ax, title in zip([ax_p, ax_s], ["Vertical (Z)", "North (T)"]):
        target_ax.set_xticks(xticks_val)
        target_ax.set_yticks(yticks_val)
        target_ax.set_xlim(dist_min, dist_max)
        target_ax.set_ylim(time_min, time_max)
        target_ax.set_xlabel("Distance (deg)", fontsize=bigger_size)
        target_ax.tick_params(axis="both", which="major", labelsize=medium_size, width=spine_width, length=6)
        for spine in target_ax.spines.values():
            spine.set_linewidth(spine_width)
        target_ax.set_title(title, fontsize=bigger_size, fontweight="bold")

    ax_p.set_ylabel("Time (s)", fontsize=bigger_size)

    if show:
        plt.show()
    return fig, (ax_p, ax_s)
