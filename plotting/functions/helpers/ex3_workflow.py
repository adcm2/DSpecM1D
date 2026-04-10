import matplotlib.pyplot as plt
import numpy as np

from .ex1_workflow import print_workflow_context
from .plot_common import (
    apply_whitegrid_style,
    get_series_colors,
    set_symmetric_ylim_and_ticks,
    style_component_axis,
)


def _prepare_ex3_state(comparison, lowidx=400, upidx=1301):
    dmf = comparison.data
    colors = get_series_colors()

    time_vector = dmf[lowidx:upidx, 0].copy()
    dspecm_z = dmf[lowidx:upidx, 1].copy()
    dspecm_n = dmf[lowidx:upidx, 2].copy()
    dspecm_e = dmf[lowidx:upidx, 3].copy()
    yspec_z = dmf[lowidx:upidx, 4].copy()
    yspec_n = dmf[lowidx:upidx, 5].copy()
    yspec_e = dmf[lowidx:upidx, 6].copy()
    mineos_z = dmf[lowidx:upidx, 7].copy()
    mineos_n = dmf[lowidx:upidx, 8].copy()
    mineos_e = dmf[lowidx:upidx, 9].copy()
    specnm_z = dmf[lowidx:upidx, 10].copy()
    specnm_n = dmf[lowidx:upidx, 11].copy()
    specnm_e = dmf[lowidx:upidx, 12].copy()

    tlow = time_vector[0]
    tup = time_vector[-1]

    norm_z = comparison.normalization(yspec_z, dspecm_z, mineos_z, eps=1e-17)
    norm_n = comparison.normalization(yspec_n, dspecm_n, mineos_n, eps=1e-17)
    norm_e = comparison.normalization(yspec_e, dspecm_e, mineos_e, eps=1e-17)

    z_stats = comparison.summarize_pairs(
        {
            "yspec_vs_dspecm": (yspec_z, dspecm_z),
            "mineos_vs_dspecm": (mineos_z, dspecm_z),
            "specnm_vs_dspecm": (specnm_z, dspecm_z),
            "specnm_vs_yspec": (specnm_z, yspec_z),
        },
        norm_z,
    )
    n_stats = comparison.summarize_pairs(
        {
            "yspec_vs_dspecm": (yspec_n, dspecm_n),
            "mineos_vs_dspecm": (mineos_n, dspecm_n),
            "specnm_vs_dspecm": (specnm_n, dspecm_n),
            "specnm_vs_yspec": (specnm_n, yspec_n),
        },
        norm_n,
    )
    e_stats = comparison.summarize_pairs(
        {
            "yspec_vs_dspecm": (yspec_e, dspecm_e),
            "mineos_vs_dspecm": (mineos_e, dspecm_e),
            "specnm_vs_dspecm": (specnm_e, dspecm_e),
            "specnm_vs_yspec": (specnm_e, yspec_e),
        },
        norm_e,
    )

    yspec_av_diff_z = z_stats["yspec_vs_dspecm"].average_percent
    specnm_av_diff_z = z_stats["specnm_vs_dspecm"].average_percent
    yspec_av_diff_n = n_stats["yspec_vs_dspecm"].average_percent
    specnm_av_diff_n = n_stats["specnm_vs_dspecm"].average_percent
    yspec_av_diff_e = e_stats["yspec_vs_dspecm"].average_percent
    specnm_av_diff_e = e_stats["specnm_vs_dspecm"].average_percent

    l2_z = comparison.l2_relative_many(dspecm_z, {"YSpec": yspec_z, "MINEOS": mineos_z, "specnm": specnm_z})
    l2_n = comparison.l2_relative_many(dspecm_n, {"YSpec": yspec_n, "MINEOS": mineos_n, "specnm": specnm_n})
    l2_e = comparison.l2_relative_many(dspecm_e, {"YSpec": yspec_e, "MINEOS": mineos_e, "specnm": specnm_e})

    yspec_z /= norm_z
    yspec_n /= norm_n
    yspec_e /= norm_e
    dspecm_z /= norm_z
    dspecm_n /= norm_n
    dspecm_e /= norm_e
    mineos_z /= norm_z
    mineos_n /= norm_n
    mineos_e /= norm_e
    specnm_z /= norm_z
    specnm_n /= norm_n
    specnm_e /= norm_e

    return {
        "colors": colors,
        "time_vector": time_vector,
        "tlow": tlow,
        "tup": tup,
        "yspec_z": yspec_z,
        "yspec_n": yspec_n,
        "yspec_e": yspec_e,
        "dspecm_z": dspecm_z,
        "dspecm_n": dspecm_n,
        "dspecm_e": dspecm_e,
        "mineos_z": mineos_z,
        "mineos_n": mineos_n,
        "mineos_e": mineos_e,
        "specnm_z": specnm_z,
        "specnm_n": specnm_n,
        "specnm_e": specnm_e,
        "z_stats": z_stats,
        "n_stats": n_stats,
        "e_stats": e_stats,
        "l2_z": l2_z,
        "l2_n": l2_n,
        "l2_e": l2_e,
        "yspec_av_diff_z": yspec_av_diff_z,
        "specnm_av_diff_z": specnm_av_diff_z,
        "yspec_av_diff_n": yspec_av_diff_n,
        "specnm_av_diff_n": specnm_av_diff_n,
        "yspec_av_diff_e": yspec_av_diff_e,
        "specnm_av_diff_e": specnm_av_diff_e,
    }


def output_to_screen_code_diff_ex3(comparison, lowidx=400, upidx=1301):
    s = _prepare_ex3_state(comparison, lowidx=lowidx, upidx=upidx)

    print_workflow_context(
        data_file=comparison.path,
        primary_label="Index window",
        primary_range=f"[{lowidx}, {upidx})",
        extras={"Time limits": f"tlow={s['tlow']}, tup={s['tup']}"},
    )

    print(f"Time limits: tlow = {s['tlow']}, tup = {s['tup']}")
    print(
        f"Average relative difference for Z component: YSpec = {s['z_stats']['yspec_vs_dspecm'].average_percent:.2f} %, "
        f"MINEOS = {s['z_stats']['mineos_vs_dspecm'].average_percent:.2f} %, "
        f"specnm = {s['z_stats']['specnm_vs_dspecm'].average_percent:.2f} %"
    )
    print(
        f"Average relative difference for N component: YSpec = {s['n_stats']['yspec_vs_dspecm'].average_percent:.2f} %, "
        f"MINEOS = {s['n_stats']['mineos_vs_dspecm'].average_percent:.2f} %, "
        f"specnm = {s['n_stats']['specnm_vs_dspecm'].average_percent:.2f} %"
    )
    print(
        f"Average relative difference for E component: YSpec = {s['e_stats']['yspec_vs_dspecm'].average_percent:.2f} %, "
        f"MINEOS = {s['e_stats']['mineos_vs_dspecm'].average_percent:.2f} %, "
        f"specnm = {s['e_stats']['specnm_vs_dspecm'].average_percent:.2f} %"
    )
    print("\n")

    print(
        f"Max relative difference for Z component: YSpec = {s['z_stats']['yspec_vs_dspecm'].max_percent:.2f} %, "
        f"MINEOS = {s['z_stats']['mineos_vs_dspecm'].max_percent:.2f} %, "
        f"specnm = {s['z_stats']['specnm_vs_dspecm'].max_percent:.2f} %"
    )
    print(
        f"Max relative difference for N component: YSpec = {s['n_stats']['yspec_vs_dspecm'].max_percent:.2f} %, "
        f"MINEOS = {s['n_stats']['mineos_vs_dspecm'].max_percent:.2f} %, "
        f"specnm = {s['n_stats']['specnm_vs_dspecm'].max_percent:.2f} %"
    )
    print(
        f"Max relative difference for E component: YSpec = {s['e_stats']['yspec_vs_dspecm'].max_percent:.2f} %, "
        f"MINEOS = {s['e_stats']['mineos_vs_dspecm'].max_percent:.2f} %, "
        f"specnm = {s['e_stats']['specnm_vs_dspecm'].max_percent:.2f} %"
    )
    print("\n")

    print(
        f"Average relative difference between specnm and YSpec for Z component: "
        f"{s['z_stats']['specnm_vs_yspec'].average_percent:.2f} %"
    )
    print(
        f"Average relative difference between specnm and YSpec for N component: "
        f"{s['n_stats']['specnm_vs_yspec'].average_percent:.2f} %"
    )
    print(
        f"Average relative difference between specnm and YSpec for E component: "
        f"{s['e_stats']['specnm_vs_yspec'].average_percent:.2f} %"
    )

    print(
        f"L2 norm for Z component: YSpec = {s['l2_z']['YSpec']:.4e}, "
        f"MINEOS = {s['l2_z']['MINEOS']:.4e}, specnm = {s['l2_z']['specnm']:.4e}"
    )
    print(
        f"L2 norm for N component: YSpec = {s['l2_n']['YSpec']:.4e}, "
        f"MINEOS = {s['l2_n']['MINEOS']:.4e}, specnm = {s['l2_n']['specnm']:.4e}"
    )
    print(
        f"L2 norm for E component: YSpec = {s['l2_e']['YSpec']:.4e}, "
        f"MINEOS = {s['l2_e']['MINEOS']:.4e}, specnm = {s['l2_e']['specnm']:.4e}"
    )


def plot_ex3_comparison(comparison, lowidx=400, upidx=1301, show=True):
    s = _prepare_ex3_state(comparison, lowidx=lowidx, upidx=upidx)

    lwidth = 2.0
    m_size = 15
    bigger_size = 20
    colors = s["colors"]

    fig, axes = plt.subplots(3, 1, figsize=(14, 12), sharex=True)
    apply_whitegrid_style(xtick_size=16, ytick_size=16)

    ax_data = axes[0]
    ax_data.plot(s["time_vector"], s["yspec_z"], color=colors["yspec"], linestyle="-", linewidth=lwidth, label="YSpec")
    ax_data.plot(s["time_vector"], s["specnm_z"], color=colors["specnm"], linestyle="--", linewidth=lwidth, label="specnm")
    ax_data.plot(s["time_vector"], s["dspecm_z"], color=colors["dspecm"], linestyle="-.", linewidth=lwidth, label="DSpecM1D")

    ax_data = axes[1]
    ax_data.plot(s["time_vector"], s["yspec_n"], color=colors["yspec"], linestyle="-", linewidth=lwidth)
    ax_data.plot(s["time_vector"], s["specnm_n"], color=colors["specnm"], linestyle="--", linewidth=lwidth)
    ax_data.plot(s["time_vector"], s["dspecm_n"], color=colors["dspecm"], linestyle="-.", linewidth=lwidth)

    ax_data = axes[2]
    ax_data.plot(s["time_vector"], s["yspec_e"], color=colors["yspec"], linestyle="-", linewidth=lwidth)
    ax_data.plot(s["time_vector"], s["specnm_e"], color=colors["specnm"], linestyle="--", linewidth=lwidth)
    ax_data.plot(s["time_vector"], s["dspecm_e"], color=colors["dspecm"], linestyle="-.", linewidth=lwidth)

    offset = 10

    ax_data = axes[0]
    ax_data.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 1.3),
        ncol=3,
        fancybox=True,
        shadow=True,
        fontsize=bigger_size,
        frameon=True,
        edgecolor="black",
    )
    style_component_axis(ax_data, offset=offset, label_size=m_size, hide_bottom_spine=True, hide_x_ticks=True)
    z_mval = 1.0
    z_axlim = 1.1 * z_mval
    set_symmetric_ylim_and_ticks(ax_data, z_axlim)

    ax_data = axes[1]
    style_component_axis(ax_data, offset=offset, label_size=m_size, hide_bottom_spine=True, hide_x_ticks=True)
    n_mval = np.max(np.abs(s["yspec_n"]))
    n_mval = 1.0
    n_axlim = 1.1 * n_mval
    set_symmetric_ylim_and_ticks(ax_data, n_axlim)

    ax_data = axes[2]
    style_component_axis(ax_data, offset=offset, label_size=m_size, hide_bottom_spine=False, hide_x_ticks=False)
    e_mval = np.max(np.abs(s["yspec_e"]))
    e_mval = 1.0
    e_axlim = 1.1 * e_mval
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

    txval = s["tlow"]
    axes[0].text(txval, 0.8 * z_mval, f"{s['yspec_av_diff_z']:.2f} %", fontsize=m_size, color=colors["yspec"])
    axes[0].text(txval, 0.6 * z_mval, f"{s['specnm_av_diff_z']:.2f} %", fontsize=m_size, color=colors["specnm"])
    axes[0].text(s["tup"], -z_mval, "Z acceleration", fontsize=bigger_size, color=colors["text"], ha="right")

    axes[1].text(txval, 0.8 * n_mval, f"{s['yspec_av_diff_n']:.2f} %", fontsize=m_size, color=colors["yspec"])
    axes[1].text(txval, 0.6 * n_mval, f"{s['specnm_av_diff_n']:.2f} %", fontsize=m_size, color=colors["specnm"])
    axes[1].text(s["tup"], -n_mval, "N acceleration", fontsize=bigger_size, color=colors["text"], ha="right")

    axes[2].text(txval, 0.8 * e_mval, f"{s['yspec_av_diff_e']:.2f} %", fontsize=m_size, color=colors["yspec"])
    axes[2].text(txval, 0.6 * e_mval, f"{s['specnm_av_diff_e']:.2f} %", fontsize=m_size, color=colors["specnm"])
    axes[2].text(s["tup"], -e_mval, "E acceleration", fontsize=bigger_size, color=colors["text"], ha="right")

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    if show:
        plt.show()
    return fig, axes
