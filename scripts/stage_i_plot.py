#!/usr/bin/env python3
"""Create the Stage I median wall-time plot and derived speedup table."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import statistics


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("timings", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()

    records = list(csv.DictReader(args.timings.open(), delimiter="\t"))
    grouped: dict[tuple[float, str], list[float]] = {}
    for record in records:
        key = (float(record["fmax_mHz"]), record["backend"])
        grouped.setdefault(key, []).append(float(record["seconds"]))
    fmax_values = sorted({key[0] for key in grouped})
    summary = []
    for fmax in fmax_values:
        eigen = grouped[(fmax, "eigen")]
        lapack = grouped[(fmax, "lapack")]
        eigen_median = statistics.median(eigen)
        lapack_median = statistics.median(lapack)
        summary.append({
            "fmax_mHz": fmax,
            "eigen_median_s": eigen_median,
            "lapack_median_s": lapack_median,
            "eigen_min_s": min(eigen),
            "eigen_max_s": max(eigen),
            "lapack_min_s": min(lapack),
            "lapack_max_s": max(lapack),
            "speedup_eigen_over_lapack": eigen_median / lapack_median,
        })
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")

    import matplotlib.pyplot as plt
    plt.figure(figsize=(7.0, 4.5))
    plt.plot(fmax_values, [row["eigen_median_s"] for row in summary], "o-",
             label="Eigen SparseLU")
    plt.plot(fmax_values, [row["lapack_median_s"] for row in summary], "o-",
             label="LAPACK band LU")
    plt.xlabel(r"$f_{\max}$ [mHz]")
    plt.ylabel("median complete spectra wall time [s]")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.output, dpi=160)
    print(json.dumps({"figure": str(args.output), "summary": str(args.summary)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
