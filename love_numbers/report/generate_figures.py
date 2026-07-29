#!/usr/bin/env python3

import csv
from pathlib import Path

import matplotlib


matplotlib.use("Agg")
import matplotlib.pyplot as plt


REPORT_DIRECTORY = Path(__file__).resolve().parent
DATA_DIRECTORY = REPORT_DIRECTORY / "data"
FIGURE_DIRECTORY = REPORT_DIRECTORY / "figures"


def read_csv(name):
    with (DATA_DIRECTORY / name).open(encoding="utf-8") as stream:
        rows = (line for line in stream if not line.startswith("#"))
        return list(csv.DictReader(rows))


def configure():
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.5,
        "lines.linewidth": 1.4,
        "lines.markersize": 4.5,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,
    })


def save(figure, name):
    FIGURE_DIRECTORY.mkdir(parents=True, exist_ok=True)
    figure.savefig(FIGURE_DIRECTORY / name)
    plt.close(figure)


def background_gravity():
    rows = read_csv("background_gravity_sfs.csv")
    elements = [int(row["elements"]) for row in rows]
    curves = (
        ("dspecm1d_absolute_error", "DSpecM1D", "o"),
        ("gia3d_absolute_error", "gia3D", "^"),
    )
    figure, axis = plt.subplots(figsize=(5.8, 3.6))
    for column, label, marker in curves:
        axis.loglog(
            elements,
            [float(row[column]) for row in rows],
            marker=marker,
            label=label,
        )
    axis.set_xlabel("Total radial elements, $N$")
    axis.set_ylabel(r"Absolute surface-gravity error (m s$^{-2}$)")
    axis.set_title("Controlled solid--fluid--solid sphere")
    axis.legend(frameon=False)
    axis.set_xticks(elements)
    axis.set_xticklabels(elements)
    save(figure, "background_gravity_current.pdf")


def controlled_fluid():
    rows = read_csv("controlled_sfs_convergence.csv")
    figure, axis = plt.subplots(figsize=(5.8, 3.6))
    for degree, marker in ((1, "o"), (2, "s"), (10, "^")):
        selected = [
            row for row in rows if int(row["degree"]) == degree
        ]
        axis.loglog(
            [int(row["elements"]) for row in selected],
            [float(row["maximum_relative_difference"]) for row in selected],
            marker=marker,
            label=fr"$l={degree}$",
        )
    axis.set_xlabel("Total radial elements, $N$")
    axis.set_ylabel("Maximum component-wise relative difference")
    axis.set_title("DSpecM1D minus gia3D: controlled internal fluid")
    axis.legend(frameon=False)
    axis.set_xticks([8, 16, 32, 64, 128, 256])
    axis.set_xticklabels([8, 16, 32, 64, 128, 256])
    save(figure, "controlled_fluid_convergence.pdf")


def isotropic_prem():
    rows = read_csv("isotropic_prem_knot_convergence.csv")
    curves = (
        (
            "dspec_sampled_minus_gia3d_sampled",
            "DSpecM1D sampled $-$ gia3D sampled",
            "o",
        ),
        (
            "gia3d_sampled_minus_gia3d_analytic",
            "gia3D sampled $-$ gia3D analytic",
            "s",
        ),
        (
            "dspec_sampled_minus_gia3d_analytic",
            "DSpecM1D sampled $-$ gia3D analytic",
            "^",
        ),
    )
    figure, axis = plt.subplots(figsize=(5.8, 3.8))
    for direction, label, marker in curves:
        selected = [row for row in rows if row["direction"] == direction]
        selected.sort(
            key=lambda row: float(row["knot_spacing_km"]), reverse=True
        )
        axis.loglog(
            [float(row["knot_spacing_km"]) for row in selected],
            [float(row["maximum_relative_difference"]) for row in selected],
            marker=marker,
            label=label,
        )
    axis.invert_xaxis()
    axis.set_xlabel("Maximum model-knot spacing (km)")
    axis.set_ylabel("Maximum component-wise relative difference")
    axis.set_title("Isotropized no-ocean PREM, radial step 0.0025")
    axis.legend(frameon=False)
    axis.set_xticks([100.0, 50.0, 25.0, 12.5])
    axis.set_xticklabels(["100", "50", "25", "12.5"])
    save(figure, "isotropic_prem_convergence.pdf")


def elln_sensitivity():
    rows = read_csv("elln_gravitational_constant.csv")
    labels = ("Production $G$", "ELLN $G$")
    h_errors = [abs(float(row["h_relative_difference"])) for row in rows]
    k_errors = [abs(float(row["k_relative_difference"])) for row in rows]
    positions = (0.0, 1.0)
    width = 0.34

    figure, axis = plt.subplots(figsize=(5.8, 3.6))
    axis.bar(
        [position - width / 2 for position in positions],
        h_errors,
        width,
        label=r"$h_{10}$",
    )
    axis.bar(
        [position + width / 2 for position in positions],
        k_errors,
        width,
        label=r"$k_{10}$",
    )
    axis.set_yscale("log")
    axis.set_xticks(positions)
    axis.set_xticklabels(labels)
    axis.set_ylabel("Absolute relative difference from ELLN Table S7")
    axis.set_title("Degree-10 TI-PREM gravitational-constant sensitivity")
    axis.legend(frameon=False)
    save(figure, "elln_gravitational_constant.pdf")


def main():
    configure()
    background_gravity()
    controlled_fluid()
    isotropic_prem()
    elln_sensitivity()
    print(f"Wrote four vector figures to {FIGURE_DIRECTORY}")


if __name__ == "__main__":
    main()
