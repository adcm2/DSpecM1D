import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from .ex1_workflow import print_workflow_context
from .plot_common import OKABE_ITO, apply_whitegrid_style, loadtxt_or_exit


def output_to_screen_code_diff_ex6(path_mf, path_mf2, path_mf3, path_freqs, path_n2, max_freq_panels=6):
    print_workflow_context(
        data_file=path_mf,
        primary_label="Additional files",
        primary_range=f"{path_mf2}, {path_mf3}, {path_freqs}, {path_n2}",
        extras={"Max frequency panels": max_freq_panels},
    )


def plot_ex6(
    path_mf,
    path_mf2,
    path_mf3,
    path_freqs,
    path_n2,
    max_freq_panels=6,
    show=True,
):
    dmf = loadtxt_or_exit(path_mf)
    dmf2 = loadtxt_or_exit(path_mf2)
    dmf3 = loadtxt_or_exit(path_mf3)
    dfreqs = loadtxt_or_exit(path_freqs)
    dn2 = loadtxt_or_exit(path_n2)
    

    nfreqs = dfreqs.shape[0]
    if nfreqs > max_freq_panels:
        nfreqs = max_freq_panels

    fig, axes = plt.subplots(3, nfreqs + 1, figsize=(14, 15), sharey=True)
    # fig, axes = plt.subplots(3, nfreqs + 1, figsize=(14, 15))
    fig.subplots_adjust(wspace=0, hspace=0.1)
    apply_whitegrid_style()

    lw = 1.2
    m_size = 12
    bigger_size = 18

    colors = {
        "u": OKABE_ITO["blue"],
        "v": OKABE_ITO["vermillion"],
        "n2": OKABE_ITO["bluish_green"],
    }

    rad_val = dmf[:, 0]
    rad_val2 = dmf2[:, 0]
    rad_val3 = dmf3[:, 0]

    low_rad = dn2[0, 0]
    up_rad = dn2[-1, 0]
    roman_numerals = ["i)", "ii)", "iii)", "iv)", "v)", "vi)"]

    def robust_panel_xlim(*series, percentile=99.5, pad=1.1, fallback=1.0):
        finite = []
        for values in series:
            arr = np.asarray(values, dtype=float)
            arr = np.abs(arr[np.isfinite(arr)])
            if arr.size:
                finite.append(arr)

        if not finite:
            return fallback

        combined = np.concatenate(finite)
        robust = np.percentile(combined, percentile)
        max_abs = np.max(combined)

        if robust <= 0:
            robust = max_abs
        if robust <= 0:
            return fallback

        # Use a robust percentile so single-point spikes do not visually
        # collapse the rest of the panel onto x=0.
        limit = pad * robust
        return max(limit, fallback)

    def plot_row(row_idx, data_mf, rad_v):
        for i in range(nfreqs):
            ax_data = axes[row_idx, i]

            u_series = data_mf[:, 6 * i + 1]
            v_series = data_mf[:, 6 * i + 4]
            # don't need robust one if there isn't an ocean layer, but it doesn't hurt to have it as a fallback in case of any unexpected spikes in the data. The robust percentile will help ensure that the x-axis limits are set in a way that allows the main features of the plot to be visible without being overly influenced by outliers.
            # maxvalt = 1.1 * robust_panel_xlim(u_series, v_series)
            maxvalt = max(np.abs(u_series).max(), np.abs(v_series).max()) * 1.1
            ax_data.plot(u_series, rad_v, color=colors["u"], linewidth=lw, label="U")
            ax_data.plot(v_series, rad_v, color=colors["v"], linewidth=lw, label="V")

            label_text = f"{roman_numerals[i]}"
            ax_data.text(0.05, 0.95, label_text, transform=ax_data.transAxes, fontsize=m_size, verticalalignment="top", fontweight="bold")

            ax_data.axhline(low_rad, color="k", linestyle="--", linewidth=1)
            ax_data.axhline(up_rad, color="k", linestyle="--", linewidth=1)
            ax_data.axvline(0, color="k", linestyle="--", linewidth=1)

            ax_data.set_ylim([0, rad_v[-1]])
            ax_data.set_xlim([-maxvalt, maxvalt])

            if i == 0:
                y_offset = 0.015 * rad_v[-1]
                x_label = -0.92 * maxvalt
                ax_data.text(
                    x_label,
                    low_rad + y_offset,
                    "ICB",
                    fontsize=m_size,
                    fontweight="bold",
                    va="bottom",
                    ha="left",
                )
                ax_data.text(
                    x_label,
                    up_rad + y_offset,
                    "CMB",
                    fontsize=m_size,
                    fontweight="bold",
                    va="bottom",
                    ha="left",
                )

            ax_data.set_xticks([0])
            ax_data.set_xticklabels([0])

            ax_data.spines["left"].set_linewidth(1.5)
            ax_data.spines["right"].set_linewidth(1.5)
            ax_data.spines["bottom"].set_linewidth(1.5)
            ax_data.spines["top"].set_linewidth(1.5)
            ax_data.tick_params(axis="x", which="major", length=5, width=1.5, labelcolor="black", labelsize=m_size)

            if i > 0:
                ax_data.set_yticklabels([])
                ax_data.set_yticks([])
            if row_idx < 2:
                ax_data.set_xticklabels([])
                ax_data.set_xticks([])

    plot_row(0, dmf, rad_val)
    plot_row(1, dmf2, rad_val2)
    plot_row(2, dmf3, rad_val3)

    for r in range(3):
        axes[r, 0].set_ylabel("Radius (N.D.)", fontsize=m_size, fontweight="bold")

    ax_legend = axes[0, 0]
    legend_properties = {"weight": "bold", "size": bigger_size}
    legend_handles = [
        Line2D([0], [0], color=colors["u"], linewidth=lw, label="U"),
        Line2D([0], [0], color=colors["v"], linewidth=lw, label="V"),
        Line2D([0], [0], color=colors["n2"], linewidth=lw, label=r"BV freq. ($N^2$)"),
    ]
    ax_legend.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(2.5, 1.25),
        ncol=4,
        fancybox=True,
        shadow=True,
        frameon=True,
        edgecolor="black",
        prop=legend_properties,
    )
    

    for r in range(3):
        ax_n2 = axes[r, nfreqs]
        ax_n2.plot(dn2[:, 1], dn2[:, 0], color=colors["n2"], linewidth=lw, label=r"$N^2$")
        ax_n2.axvline(0, color="k", linestyle="--", linewidth=1)
        ax_n2.axhline(low_rad, color="k", linestyle="--", linewidth=1)
        ax_n2.axhline(up_rad, color="k", linestyle="--", linewidth=1)
        ax_n2.set_xticks([0])
        ax_n2.set_xticklabels([0])
        ax_n2.set_yticklabels([])
        ax_n2.set_yticks([])

        ax_n2.spines["left"].set_linewidth(1.5)
        ax_n2.spines["right"].set_linewidth(1.5)
        ax_n2.spines["bottom"].set_linewidth(1.5)
        ax_n2.spines["top"].set_linewidth(1.5)

        ax_n2.text(0.05, 0.95, "v)", transform=ax_n2.transAxes, fontsize=m_size, verticalalignment="top", fontweight="bold")
        ax_n2.tick_params(axis="x", which="major", length=5, width=1.2, labelcolor="black", labelsize=m_size)
        if r < 2:
            ax_n2.set_xticklabels([])
            ax_n2.set_xticks([])

    plt.tight_layout(pad=0.1, w_pad=0.0, h_pad=0.0)
    fig.subplots_adjust(wspace=0, hspace=0.05)

    
    if show:
        plt.show()
    return fig, axes
