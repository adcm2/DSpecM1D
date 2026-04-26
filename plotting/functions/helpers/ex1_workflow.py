import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch

from .mode_annotations import annotate_modes, compute_mode_heights, get_modes_in_range
from .plot_common import apply_whitegrid_style, get_series_colors


def print_workflow_context(data_file, primary_label, primary_range=None, extras=None):
    print(f"Loaded file: {data_file}")
    if primary_range is not None:
        print(f"{primary_label}: {primary_range}")
    if extras:
        for key, value in extras.items():
            print(f"{key}: {value}")


def output_to_screen_code_diff(comparison, fmin, fmax, zoom_regions=None, specnm_preferred_col=21):
    report = comparison.ex1_diff_report(fmin=fmin, fmax=fmax, specnm_preferred_col=specnm_preferred_col)
    extras = {"specnm column used": report.specnm_column}
    if zoom_regions:
        extras["Zoom regions"] = zoom_regions
    print_workflow_context(
        data_file=comparison.path,
        primary_label="Frequency range",
        primary_range=f"[{report.fmin:.3f}, {report.fmax:.3f}] mHz",
        extras=extras,
    )
    print(
        "YSpec vs DSpecM1D: "
        f"avg={report.yspec_vs_dspecm.average_percent:.4f} %, "
        f"max={report.yspec_vs_dspecm.max_percent:.4f} %"
    )
    print(
        "specnm vs DSpecM1D: "
        f"avg={report.specnm_vs_dspecm.average_percent:.4f} %, "
        f"max={report.specnm_vs_dspecm.max_percent:.4f} %"
    )


def _add_zoom_inset(
    fig,
    ax_main,
    frequency,
    yspec_plot,
    specnm_plot,
    dspecm_plot,
    zoom_region,
    inset_rect,
    colors,
    lwidth,
    lwidth2,
    connector_targets,
    mode_label_dy,
    alternate_labels=False,
):
    zoom_fmin, zoom_fmax = zoom_region
    x1, x2 = zoom_fmin, zoom_fmax

    idx_zmin = np.searchsorted(frequency, zoom_fmin)
    idx_zmax = np.searchsorted(frequency, zoom_fmax, side="right")

    zoom_x = frequency[idx_zmin:idx_zmax]
    zoom_y_yspec = yspec_plot[idx_zmin:idx_zmax]
    zoom_y_specnm = specnm_plot[idx_zmin:idx_zmax]
    zoom_y_dspecm = dspecm_plot[idx_zmin:idx_zmax]

    inset_ax = fig.add_axes(inset_rect)
    inset_ax.plot(zoom_x, zoom_y_yspec, color=colors["yspec"], linestyle="-", linewidth=lwidth)
    inset_ax.plot(zoom_x, zoom_y_specnm, color=colors["specnm"], linestyle="--", linewidth=lwidth)
    inset_ax.plot(zoom_x, zoom_y_dspecm, color=colors["dspecm"], linestyle="-.", linewidth=lwidth)

    zoom_ymax = max(np.max(zoom_y_yspec), np.max(zoom_y_specnm), np.max(zoom_y_dspecm)) * 1.07
    y1, y2 = 0, zoom_ymax

    inset_ax.set_xlim(zoom_fmin, zoom_fmax)
    inset_ax.set_ylim(0, zoom_ymax)
    inset_ax.tick_params(axis="x", which="major", labelsize=12, labelcolor="black", colors="black")
    inset_ax.set_yticks([])
    inset_ax.grid(False)

    for spine in inset_ax.spines.values():
        spine.set_edgecolor("black")
        spine.set_linewidth(lwidth2)

    (bx1, by1), (bx2, by2) = connector_targets
    con_a = ConnectionPatch(
        xyA=(x1, y2),
        coordsA=ax_main.transData,
        xyB=(bx1, by1),
        coordsB=inset_ax.transAxes,
        color="black",
        linestyle="-",
        linewidth=lwidth2,
    )
    fig.add_artist(con_a)

    con_b = ConnectionPatch(
        xyA=(x2, y2),
        coordsA=ax_main.transData,
        xyB=(bx2, by2),
        coordsB=inset_ax.transAxes,
        color="black",
        linestyle="-",
        linewidth=lwidth2,
    )
    fig.add_artist(con_b)

    rect = plt.Rectangle((x1, y1), x2 - x1, y2 - y1, facecolor="none", edgecolor="black", linewidth=lwidth2)
    ax_main.add_patch(rect)

    names, xvalues = get_modes_in_range(zoom_fmin, zoom_fmax)
    yvalues = compute_mode_heights(xvalues, zoom_x, [zoom_y_yspec, zoom_y_dspecm, zoom_y_specnm])
    annotate_modes(
        inset_ax,
        names,
        xvalues,
        yvalues,
        y2,
        dy=mode_label_dy,
        line_width=lwidth2,
        alternate=alternate_labels,
    )


def plot_ex1_comparison(
    comparison,
    fmin,
    fmax,
    zoom_regions=None,
    specnm_preferred_col=21,
    show=True,
):
    colors = get_series_colors()
    dataset = comparison.ex1_dataset(fmin=fmin, fmax=fmax, specnm_preferred_col=specnm_preferred_col)

    data_offset = 1e-3
    yspec_plot = dataset.yspec_z / dataset.norm_z + data_offset
    specnm_plot = dataset.specnm_z / dataset.norm_z + data_offset
    dspecm_plot = dataset.dspecm_z / dataset.norm_z + data_offset

    apply_whitegrid_style()
    fig, ax_data = plt.subplots(1, 1, figsize=(14, 12), sharex=True)

    lwidth = 1.3
    lwidth2 = 2.0
    m_size = 15
    bigger_size = 20
    offset = 10

    ax_data.plot(dataset.frequency, yspec_plot, color=colors["yspec"], linestyle="-", linewidth=lwidth, label="YSpec")
    ax_data.plot(dataset.frequency, specnm_plot, color=colors["specnm"], linestyle="--", linewidth=lwidth, label="specnm")
    ax_data.plot(dataset.frequency, dspecm_plot, color=colors["dspecm"], linestyle="-.", linewidth=lwidth, label="DSpecM1D")

    if zoom_regions:
        inset_layouts = [
            {
                "rect": [0.07, 0.45, 0.15, 0.45],
                "connectors": ((0, 0), (1, 0)),
                "dy": 0.002,
                "alternate": False,
            },
            {
                "rect": [0.25, 0.45, 0.15, 0.45],
                "connectors": ((1, 0), (1, 1)),
                "dy": 0.01,
                "alternate": True,
            },
        ]

        for i, region in enumerate(zoom_regions[:2]):
            cfg = inset_layouts[i]
            _add_zoom_inset(
                fig=fig,
                ax_main=ax_data,
                frequency=dataset.frequency,
                yspec_plot=yspec_plot,
                specnm_plot=specnm_plot,
                dspecm_plot=dspecm_plot,
                zoom_region=region,
                inset_rect=cfg["rect"],
                colors=colors,
                lwidth=lwidth,
                lwidth2=lwidth2,
                connector_targets=cfg["connectors"],
                mode_label_dy=cfg["dy"],
                alternate_labels=cfg["alternate"],
            )

    ax_data.set_xlabel("Frequency (mHz)", fontsize=bigger_size)
    ax_data.set_ylabel("Z displacement (N.D.)", fontsize=bigger_size)
    ax_data.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, 1.1),
        ncol=3,
        fancybox=True,
        shadow=True,
        fontsize=bigger_size,
        frameon=True,
        edgecolor="black",
    )

    ax_data.spines["top"].set_visible(False)
    ax_data.spines["right"].set_visible(False)
    ax_data.spines["left"].set_position(("outward", offset))
    ax_data.spines["left"].set_linewidth(1.5)
    ax_data.spines["left"].set_color("black")
    ax_data.spines["bottom"].set_position(("outward", offset))
    ax_data.spines["bottom"].set_linewidth(1.5)
    ax_data.spines["bottom"].set_color("black")

    ax_data.tick_params(
        axis="both",
        which="major",
        labelsize=m_size,
        length=10,
        width=1.5,
        colors="black",
        labelcolor="black",
    )
    ax_data.grid(False)
    ax_data.set_xlim([fmin, fmax])
    ax_data.set_ylim(0, 1.0)
    ax_data.set_yticks([0, 1.0])

    plt.tight_layout()
    if show:
        plt.show()
    return fig, ax_data
