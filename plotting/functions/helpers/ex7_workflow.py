import matplotlib.pyplot as plt
import numpy as np

from .ex2_workflow import _prepare_ex2_state, output_to_screen_code_diff_ex2
from .plot_common import OKABE_ITO, apply_whitegrid_style, set_symmetric_ylim_and_ticks, style_component_axis


def _resolve_window(comparison, lowidx=400, upidx=1301, tup=None):
    nrows = comparison.data.shape[0]
    resolved_lowidx = max(0, int(lowidx))
    resolved_upidx = nrows if upidx is None else min(int(upidx), nrows)
    if resolved_upidx <= resolved_lowidx:
        raise ValueError("Invalid index window for ex7 plotting")

    resolved_tup = float(comparison.data[resolved_upidx - 1, 0]) if tup is None else float(tup)
    return resolved_lowidx, resolved_upidx, resolved_tup


def output_to_screen_code_diff_ex7(comparison, lowidx=0, upidx=None, tup=None):
    lowidx, upidx, tup = _resolve_window(comparison, lowidx=lowidx, upidx=upidx, tup=tup)
    output_to_screen_code_diff_ex2(
        comparison=comparison,
        lowidx=lowidx,
        upidx=upidx,
        tup=tup,
    )


def plot_ex7_comparison(comparison, lowidx=400, upidx=1301, tup=None, show=True):
    lowidx, upidx, tup = _resolve_window(comparison, lowidx=lowidx, upidx=upidx, tup=tup)
    s = _prepare_ex2_state(comparison, lowidx=lowidx, upidx=upidx, tup=tup)

    lwidth = 1.9
    m_size = 15
    bigger_size = 20
    colors = dict(s["colors"])
    colors["mineos"] = OKABE_ITO["orange"]
    colors["dspecm"] = OKABE_ITO["black"]

    line_styles = {
        "yspec": "-",
        "mineos": "-.",
        "specnm": "--",
        "dspecm": "-",
    }
    line_widths = {
        "yspec": lwidth,
        "mineos": lwidth,
        "specnm": lwidth,
        "dspecm": .5*lwidth,
    }
    line_alpha = {
        "yspec": 0.95,
        "mineos": 0.95,
        "specnm": 0.95,
        "dspecm": 0.95,
    }
    line_zorder = {
        "dspecm": 4,
        "yspec": 1,
        "mineos": 2,
        "specnm": 3,
    }

    fig, axes = plt.subplots(3, 1, figsize=(14, 12), sharex=True)
    apply_whitegrid_style(xtick_size=16, ytick_size=16)

    ax_data = axes[0]
    ax_data.plot(
        s["time_vector"],
        s["yspec_z"] / 1.01,
        color=colors["yspec"],
        linestyle=line_styles["yspec"],
        linewidth=line_widths["yspec"],
        alpha=line_alpha["yspec"],
        zorder=line_zorder["yspec"],
        label="YSpec",
    )
    ax_data.plot(
        s["time_vector"],
        s["mineos_z"] / 1.01,
        color=colors["mineos"],
        linestyle=line_styles["mineos"],
        linewidth=line_widths["mineos"],
        alpha=line_alpha["mineos"],
        zorder=line_zorder["mineos"],
        label="MINEOS",
    )
    ax_data.plot(
        s["time_vector"],
        s["specnm_z"] / 1.01,
        color=colors["specnm"],
        linestyle=line_styles["specnm"],
        linewidth=line_widths["specnm"],
        alpha=line_alpha["specnm"],
        zorder=line_zorder["specnm"],
        label="specnm",
    )
    ax_data.plot(
        s["time_vector"],
        s["dspecm_z"] / 1.01,
        color=colors["dspecm"],
        linestyle=line_styles["dspecm"],
        linewidth=line_widths["dspecm"],
        alpha=line_alpha["dspecm"],
        zorder=line_zorder["dspecm"],
        label="DSpecM1D",
    )

    ax_data = axes[1]
    ax_data.plot(s["time_vector"], s["yspec_n"] / 1.01, color=colors["yspec"], linestyle=line_styles["yspec"], linewidth=line_widths["yspec"], alpha=line_alpha["yspec"], zorder=line_zorder["yspec"])
    ax_data.plot(s["time_vector"], s["mineos_n"] / 1.01, color=colors["mineos"], linestyle=line_styles["mineos"], linewidth=line_widths["mineos"], alpha=line_alpha["mineos"], zorder=line_zorder["mineos"])
    ax_data.plot(s["time_vector"], s["specnm_n"] / 1.01, color=colors["specnm"], linestyle=line_styles["specnm"], linewidth=line_widths["specnm"], alpha=line_alpha["specnm"], zorder=line_zorder["specnm"])
    ax_data.plot(s["time_vector"], s["dspecm_n"] / 1.01, color=colors["dspecm"], linestyle=line_styles["dspecm"], linewidth=line_widths["dspecm"], alpha=line_alpha["dspecm"], zorder=line_zorder["dspecm"])

    ax_data = axes[2]
    ax_data.plot(s["time_vector"], s["yspec_e"] / 1.01, color=colors["yspec"], linestyle=line_styles["yspec"], linewidth=line_widths["yspec"], alpha=line_alpha["yspec"], zorder=line_zorder["yspec"])
    ax_data.plot(s["time_vector"], s["mineos_e"] / 1.01, color=colors["mineos"], linestyle=line_styles["mineos"], linewidth=line_widths["mineos"], alpha=line_alpha["mineos"], zorder=line_zorder["mineos"])
    ax_data.plot(s["time_vector"], s["specnm_e"] / 1.01, color=colors["specnm"], linestyle=line_styles["specnm"], linewidth=line_widths["specnm"], alpha=line_alpha["specnm"], zorder=line_zorder["specnm"])
    ax_data.plot(s["time_vector"], s["dspecm_e"] / 1.01, color=colors["dspecm"], linestyle=line_styles["dspecm"], linewidth=line_widths["dspecm"], alpha=line_alpha["dspecm"], zorder=line_zorder["dspecm"])

    offset = 10

    ax_data = axes[0]
    ax_data.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 1.3),
        ncol=4,
        fancybox=True,
        shadow=True,
        fontsize=bigger_size,
        frameon=True,
        edgecolor="black",
    )
    style_component_axis(ax_data, offset=offset, label_size=m_size, hide_bottom_spine=True, hide_x_ticks=True)
    z_mval = max(np.max(np.abs(s["yspec_z"])), np.max(np.abs(s["dspecm_z"])), np.max(np.abs(s["mineos_z"])))
    z_axlim = 1.0 * z_mval
    set_symmetric_ylim_and_ticks(ax_data, z_axlim)

    ax_data = axes[1]
    style_component_axis(ax_data, offset=offset, label_size=m_size, hide_bottom_spine=True, hide_x_ticks=True)
    n_mval = max(np.max(np.abs(s["yspec_n"])), np.max(np.abs(s["dspecm_n"])), np.max(np.abs(s["mineos_n"])))
    n_axlim = 1.0 * n_mval
    set_symmetric_ylim_and_ticks(ax_data, n_axlim)

    ax_data = axes[2]
    style_component_axis(ax_data, offset=offset, label_size=m_size, hide_bottom_spine=False, hide_x_ticks=False)
    e_mval = max(np.max(np.abs(s["yspec_e"])), np.max(np.abs(s["dspecm_e"])), np.max(np.abs(s["mineos_e"])))
    e_axlim = 1.0 * e_mval
    set_symmetric_ylim_and_ticks(ax_data, e_axlim)
    xticks = [
        s["tlow"],
        s["tlow"] + (s["tup"] - s["tlow"]) * 0.25,
        s["tlow"] + (s["tup"] - s["tlow"]) * 0.5,
        s["tlow"] + (s["tup"] - s["tlow"]) * 0.75,
        s["tup"],
    ]
    ax_data.set_xticks(xticks)
    ax_data.get_xticklabels()[-1].set_horizontalalignment("right")
    ax_data.get_xticklabels()[-1].set_visible(False)

    axes[2].set_xlabel("Time (s)", fontsize=bigger_size)
    for ax in axes:
        ax.set_xlim(s["tlow"], s["tup"])

    def _annotate_avg_diffs(ax, yspec_diff, mineos_diff, specnm_diff):
        ax.text(
            0.02,
            0.92,
            f"YSpec-DSpec: {yspec_diff:.2f} %",
            transform=ax.transAxes,
            fontsize=m_size,
            color=colors["yspec"],
            va="top",
        )
        ax.text(
            0.02,
            0.84,
            f"MINEOS-DSpec: {mineos_diff:.2f} %",
            transform=ax.transAxes,
            fontsize=m_size,
            color=colors["mineos"],
            va="top",
        )
        ax.text(
            0.02,
            0.76,
            f"specnm-DSpec: {specnm_diff:.2f} %",
            transform=ax.transAxes,
            fontsize=m_size,
            color=colors["specnm"],
            va="top",
        )

    mineos_av_diff_z = s["z_stats"]["mineos_vs_dspecm"].average_percent
    mineos_av_diff_n = s["n_stats"]["mineos_vs_dspecm"].average_percent
    mineos_av_diff_e = s["e_stats"]["mineos_vs_dspecm"].average_percent

    _annotate_avg_diffs(axes[0], s["yspec_av_diff_z"], mineos_av_diff_z, s["specnm_av_diff_z"])
    axes[0].text(s["tup"], -z_mval, "Z displacement", fontsize=bigger_size, color=colors["text"], ha="right")

    _annotate_avg_diffs(axes[1], s["yspec_av_diff_n"], mineos_av_diff_n, s["specnm_av_diff_n"])
    axes[1].text(s["tup"], -n_mval, "N displacement", fontsize=bigger_size, color=colors["text"], ha="right")

    _annotate_avg_diffs(axes[2], s["yspec_av_diff_e"], mineos_av_diff_e, s["specnm_av_diff_e"])
    axes[2].text(s["tup"], -e_mval, "E displacement", fontsize=bigger_size, color=colors["text"], ha="right")

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    if show:
        plt.show()
    return fig, axes
