import matplotlib.pyplot as plt
import numpy as np

OKABE_ITO = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "bluish_green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "reddish_purple": "#CC79A7",
}


def get_series_colors():
    return {
        "yspec": OKABE_ITO["blue"],
        "specnm": OKABE_ITO["bluish_green"],
        "dspecm": OKABE_ITO["vermillion"],
        "text": "#111111",
    }


def apply_whitegrid_style(xtick_size=None, ytick_size=None):
    plt.style.use("seaborn-v0_8-whitegrid")
    if xtick_size is not None:
        plt.rc("xtick", labelsize=xtick_size)
    if ytick_size is not None:
        plt.rc("ytick", labelsize=ytick_size)


def loadtxt_or_exit(path, delimiter=";"):
    try:
        return np.loadtxt(path, delimiter=delimiter)
    except FileNotFoundError:
        print(f"Error: Data file not found at '{path}'")
        raise SystemExit(1)


def style_component_axis(ax, offset, label_size, hide_bottom_spine=True, hide_x_ticks=True):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_position(("outward", offset))
    ax.spines["left"].set_linewidth(1.5)

    if hide_bottom_spine:
        ax.spines["bottom"].set_visible(False)
    else:
        ax.spines["bottom"].set_position(("outward", offset))
        ax.spines["bottom"].set_linewidth(1.5)

    if hide_x_ticks:
        ax.tick_params(axis="x", which="both", bottom=False, top=False)

    ax.tick_params(axis="both", which="major", length=10, width=1.5)
    ax.tick_params(axis="both", which="major", labelcolor="black", labelsize=label_size)


def set_symmetric_ylim_and_ticks(ax, amplitude):
    ax.set_ylim(-amplitude, amplitude)
    ax.set_yticks([-amplitude, 0, amplitude])
