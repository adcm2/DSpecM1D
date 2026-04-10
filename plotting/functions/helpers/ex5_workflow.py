import matplotlib.pyplot as plt

from .ex1_workflow import print_workflow_context
from .plot_common import OKABE_ITO, apply_whitegrid_style, loadtxt_or_exit


def output_to_screen_code_diff_ex5(path_nq5, path_nq6, path_nq4):
    print_workflow_context(
        data_file=path_nq5,
        primary_label="Additional files",
        primary_range=f"{path_nq6}, {path_nq4}",
    )


def plot_ex5(path_nq5, path_nq6, path_nq4, show=True):
    data_3 = loadtxt_or_exit(path_nq5)
    data_6 = loadtxt_or_exit(path_nq6)
    data_4 = loadtxt_or_exit(path_nq4)

    step_sizes_nq5 = data_3[:, 0]
    t_error_1 = (data_3[:, 1] + data_3[:, 2] + data_3[:, 3]) / 3.0

    step_sizes_nq6 = data_6[:, 0]
    t_error_2 = (data_6[:, 1] + data_6[:, 2] + data_6[:, 3]) / 3.0

    step_sizes_nq4 = data_4[:, 0]
    t_error_3 = (data_4[:, 1] + data_4[:, 2] + data_4[:, 3]) / 3.0

    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    apply_whitegrid_style()

    lwidth = 2.0
    bigger_size = 24

    colors = {
        "nq4": OKABE_ITO["bluish_green"],
        "nq5": OKABE_ITO["blue"],
        "nq6": OKABE_ITO["vermillion"],
    }

    ratio_4 = 0.63 / (50 * step_sizes_nq4)
    ratio_5 = 0.63 / (50 * step_sizes_nq5)
    ratio_6 = 0.63 / (50 * step_sizes_nq6)

    ax.plot(ratio_4, t_error_3, color=colors["nq4"], linestyle="-", linewidth=lwidth, label="NQ = 4")
    ax.plot(ratio_5, t_error_1, color=colors["nq5"], linestyle="-", linewidth=lwidth, label="NQ = 5")
    ax.plot(ratio_6, t_error_2, color=colors["nq6"], linestyle="-", linewidth=lwidth, label="NQ = 6")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("# elements per wavelength", fontsize=bigger_size)
    ax.set_ylabel("Average relative error (%)", fontsize=bigger_size)
    ax.tick_params(axis="both", which="major", labelsize=bigger_size - 2)
    ax.legend(fontsize=bigger_size - 2)

    plt.tight_layout()
    if show:
        plt.show()
    return fig, ax
