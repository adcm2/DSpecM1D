#!/usr/bin/env python3
"""Plot the four-column output of generate_seismogram for one receiver."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import numpy as np


COMPONENTS = ("Radial", "North", "East")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="semicolon-delimited seismogram output")
    parser.add_argument("--save", type=Path, help="save the figure at this path")
    parser.add_argument("--no-show", action="store_true", help="do not open an interactive window")
    args = parser.parse_args()

    if args.no_show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    try:
        data = np.loadtxt(args.input, delimiter=";")
    except (OSError, ValueError) as error:
        parser.error(f"cannot read {args.input}: {error}")
    if data.ndim != 2 or data.shape[1] != 4 or data.shape[0] < 2:
        parser.error("expected at least two rows and exactly four columns (time plus one receiver's three components)")
    if not np.isfinite(data).all():
        parser.error("input contains non-finite values")
    if data[0, 0] < 0 or not np.all(np.diff(data[:, 0]) > 0):
        parser.error("time column must be nonnegative and strictly increasing")
    if not np.any(np.abs(data[:, 1:]) > 0):
        parser.error("all three seismogram components are zero")

    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(9, 7), constrained_layout=True)
    for index, (axis, label) in enumerate(zip(axes, COMPONENTS, strict=True), start=1):
        axis.plot(data[:, 0], data[:, index], color="#176b70", linewidth=1.1)
        axis.set_ylabel(f"{label} (m)")
        axis.grid(True, alpha=0.25)
    axes[-1].set_xlabel("Time (s)")
    fig.suptitle("First seismogram · displacement")
    if args.save:
        args.save.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.save, dpi=160)
    if not args.no_show:
        plt.show()
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
