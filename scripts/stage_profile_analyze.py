#!/usr/bin/env python3
"""Validate and summarize the completed adaptive multi-SEM profile campaign."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


FMAX = (20.0, 35.0, 55.0)
BACKENDS = ("eigen", "lapack")
MODES = ("radial", "toroidal", "spheroidal")
CATEGORIES = (
    "base_operator_preparation",
    "start_truncation_extraction",
    "frequency_matrix_construction",
    "sparse_compression",
    "lapack_band_packing",
    "factorization",
    "solve",
    "source_receiver_projection",
)
LABELS = {
    "base_operator_preparation": "base operators",
    "start_truncation_extraction": "start/truncation",
    "frequency_matrix_construction": "dynamic A",
    "sparse_compression": "sparse compression",
    "lapack_band_packing": "band packing",
    "factorization": "factorization",
    "solve": "solve",
    "source_receiver_projection": "source/receiver/projection",
    "unclassified": "unclassified remainder",
}
COLORS = {
    "base_operator_preparation": "#4C78A8",
    "start_truncation_extraction": "#72B7B2",
    "frequency_matrix_construction": "#F58518",
    "sparse_compression": "#B8B8B8",
    "lapack_band_packing": "#E45756",
    "factorization": "#54A24B",
    "solve": "#B279A2",
    "source_receiver_projection": "#FFBF79",
    "unclassified": "#79706E",
}


def median(values: list[float]) -> float:
    return float(statistics.median(values))


def load_and_validate(raw_path: Path, metadata_path: Path) -> tuple[list[dict], dict]:
    metadata = json.loads(metadata_path.read_text())
    if metadata.get("status") != "complete":
        raise RuntimeError("profiling metadata is not complete")
    completed = {(item["fmax_mHz"], item["backend"])
                 for item in metadata.get("completed_backends", [])}
    expected_completed = {(fmax, backend) for fmax in FMAX for backend in BACKENDS}
    if completed != expected_completed:
        raise RuntimeError("metadata must contain six completed backend blocks")
    if sorted(metadata.get("completed_points", [])) != list(FMAX):
        raise RuntimeError("metadata frequency grid is incomplete")
    required_metadata = ("git", "machine", "eigen_cmake", "lapack_cmake",
                         "profiling_source_sha256", "threads")
    if any(key not in metadata for key in required_metadata):
        raise RuntimeError("profiling metadata lacks reproducibility fields")
    repo = Path(__file__).resolve().parents[1]
    source_hashes = metadata["profiling_source_sha256"]
    if len(source_hashes) != 6:
        raise RuntimeError("expected six profiling source hashes")
    for relative_path, expected_hash in source_hashes.items():
        actual_hash = hashlib.sha256((repo / relative_path).read_bytes()).hexdigest()
        if actual_hash != expected_hash:
            raise RuntimeError(f"profiling source changed since run: {relative_path}")

    rows: list[dict] = []
    with raw_path.open() as source:
        header = source.readline().rstrip("\n").split("\t")
        expected = ["fmax_mHz", "backend", "repetition", "rows", "cols",
                    "norm", "profile_json"]
        if header != expected:
            raise RuntimeError(f"unexpected raw header: {header}")
        for line_number, line in enumerate(source, 2):
            fields = line.rstrip("\n").split("\t", 6)
            if len(fields) != 7:
                raise RuntimeError(f"malformed line {line_number}")
            row = {
                "fmax": float(fields[0]), "backend": fields[1],
                "repetition": int(fields[2]), "rows": int(fields[3]),
                "cols": int(fields[4]), "norm": float(fields[5]),
                "profile": json.loads(fields[6]),
            }
            if row["fmax"] not in FMAX or row["backend"] not in BACKENDS:
                raise RuntimeError(f"unexpected campaign key on line {line_number}")
            if row["rows"] != 6 or row["cols"] != 16385:
                raise RuntimeError(f"unexpected output shape on line {line_number}")
            if not math.isfinite(row["norm"]):
                raise RuntimeError(f"non-finite norm on line {line_number}")
            profile = row["profile"]
            if not math.isfinite(profile["wall_seconds"]) or profile["wall_seconds"] <= 0:
                raise RuntimeError(f"invalid wall time on line {line_number}")
            for mode, accounting in profile["accounting"]["workers"].items():
                values = [accounting[key] for key in
                          ("total", "categorized", "unclassified")]
                if not all(math.isfinite(value) and value >= 0 for value in values):
                    raise RuntimeError(f"invalid {mode} accounting on line {line_number}")
                if abs(values[0] - values[1] - values[2]) > 1e-8 * max(1.0, values[0]):
                    raise RuntimeError(f"worker accounting mismatch on line {line_number}")
            for timings in profile["timings_seconds"].values():
                if not all(math.isfinite(value) and value >= 0
                           for value in timings.values()):
                    raise RuntimeError(f"invalid category timing on line {line_number}")
            rows.append(row)

    if len(rows) != 18:
        raise RuntimeError(f"expected 18 records, found {len(rows)}")
    for fmax in FMAX:
        for backend in BACKENDS:
            group = [row for row in rows
                     if row["fmax"] == fmax and row["backend"] == backend]
            if sorted(row["repetition"] for row in group) != [1, 2, 3]:
                raise RuntimeError(f"invalid repetitions for {fmax:g}/{backend}")
            first = group[0]["profile"]["counts"]
            if any(row["profile"]["counts"] != first for row in group[1:]):
                raise RuntimeError(f"workload counts vary for {fmax:g}/{backend}")
            systems = first["frequency_systems"]
            if first["solves"] != systems:
                raise RuntimeError(f"solve/system mismatch for {fmax:g}/{backend}")
            if first["rhs"] < systems or first["sems"] <= 0 or first["degrees"] != 1500:
                raise RuntimeError(f"invalid workload counts for {fmax:g}/{backend}")
            if backend == "eigen":
                if first["eigen_compute"] + first["eigen_factorize"] != systems:
                    raise RuntimeError(f"Eigen cadence mismatch for {fmax:g}")
                if first["lapack_factorize"] or first["band_packs"]:
                    raise RuntimeError(f"LAPACK counts in Eigen record for {fmax:g}")
            else:
                if first["lapack_factorize"] != systems or first["band_packs"] != systems:
                    raise RuntimeError(f"LAPACK pack/factor mismatch for {fmax:g}")
                if first["eigen_compute"] or first["eigen_factorize"]:
                    raise RuntimeError(f"Eigen counts in LAPACK record for {fmax:g}")
            norms = [row["norm"] for row in group]
            if max(norms) - min(norms) > 1e-12 * max(norms):
                raise RuntimeError(f"repetition norms vary for {fmax:g}/{backend}")
        eigen_counts = next(row for row in rows if row["fmax"] == fmax and
                            row["backend"] == "eigen")["profile"]["counts"]
        lapack_counts = next(row for row in rows if row["fmax"] == fmax and
                             row["backend"] == "lapack")["profile"]["counts"]
        common = ("sems", "degrees", "frequency_systems", "solves", "rhs",
                  "dimension_min", "dimension_max", "nnz_min", "nnz_max",
                  "kl_min", "kl_max", "ku_min", "ku_max")
        if any(eigen_counts[key] != lapack_counts[key] for key in common):
            raise RuntimeError(f"backend workload mismatch for {fmax:g}")
        eigen_norm = median([row["norm"] for row in rows
                             if row["fmax"] == fmax and row["backend"] == "eigen"])
        lapack_norm = median([row["norm"] for row in rows
                              if row["fmax"] == fmax and row["backend"] == "lapack"])
        if abs(eigen_norm - lapack_norm) / eigen_norm >= 1e-10:
            raise RuntimeError(f"backend norm mismatch for {fmax:g}")
    return rows, metadata


def summarize(rows: list[dict], metadata: dict, raw_path: Path,
              metadata_path: Path) -> dict:
    summary: dict = {
        "campaign": metadata["campaign"],
        "accounting_basis": (
            "wall_seconds is complete spectra wall time; category and worker "
            "values inside OpenMP are summed thread-work seconds and are not wall-additive"),
        "input_sha256": {
            raw_path.name: hashlib.sha256(raw_path.read_bytes()).hexdigest(),
            metadata_path.name: hashlib.sha256(metadata_path.read_bytes()).hexdigest(),
        },
        "groups": [], "numerical_sanity": [],
    }
    for fmax in FMAX:
        backend_norms = {}
        backend_walls = {}
        for backend in BACKENDS:
            group = [row for row in rows
                     if row["fmax"] == fmax and row["backend"] == backend]
            walls = [row["profile"]["wall_seconds"] for row in group]
            worker_totals = [row["profile"]["accounting"]["workers"]["all"]["total"]
                             for row in group]
            worker_median = median(worker_totals)
            categories = {
                category: median([row["profile"]["timings_seconds"]["all"][category]
                                  for row in group])
                for category in CATEGORIES
            }
            reconciled_unclassified = max(0.0, worker_median - sum(categories.values()))
            categories["unclassified"] = reconciled_unclassified
            mode_totals = {
                mode: median([row["profile"]["accounting"]["workers"][mode]["total"]
                              for row in group]) for mode in MODES
            }
            summary["groups"].append({
                "fmax_mHz": fmax, "backend": backend,
                "wall_seconds": {"min": min(walls), "median": median(walls),
                                 "max": max(walls), "samples": walls},
                "worker_total_median_seconds": worker_median,
                "worker_accounting_median_seconds": {
                    key: median([
                        row["profile"]["accounting"]["workers"]["all"][key]
                        for row in group])
                    for key in ("total", "categorized", "unclassified")
                },
                "worker_categories_median_seconds": categories,
                "worker_category_share_percent": {
                    key: 100.0 * value / worker_median
                    for key, value in categories.items()},
                "worker_modes_median_seconds": mode_totals,
                "worker_mode_share_percent": {
                    key: 100.0 * value / worker_median
                    for key, value in mode_totals.items()},
                "counts": group[0]["profile"]["counts"],
                "norms": [row["norm"] for row in group],
                "max_worker_accounting_residual": max(
                    abs(accounting["total"] - accounting["categorized"] -
                        accounting["unclassified"])
                    for row in group
                    for accounting in row["profile"]["accounting"]["workers"].values()),
            })
            backend_norms[backend] = median([row["norm"] for row in group])
            backend_walls[backend] = median(walls)
        summary["numerical_sanity"].append({
            "fmax_mHz": fmax,
            "eigen_norm": backend_norms["eigen"],
            "lapack_norm": backend_norms["lapack"],
            "relative_norm_difference": abs(backend_norms["eigen"] -
                                                backend_norms["lapack"]) /
                                        backend_norms["eigen"],
            "wall_speedup_eigen_over_lapack":
                backend_walls["eigen"] / backend_walls["lapack"],
        })
    return summary


def write_csv(path: Path, summary: dict) -> None:
    fields = ["fmax_mHz", "backend", "wall_min_s", "wall_median_s",
              "wall_max_s", "worker_total_s"]
    fields += [f"{category}_s" for category in (*CATEGORIES, "unclassified")]
    fields += [f"{mode}_worker_s" for mode in MODES]
    with path.open("w", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=fields)
        writer.writeheader()
        for group in summary["groups"]:
            row = {
                "fmax_mHz": group["fmax_mHz"], "backend": group["backend"],
                "wall_min_s": group["wall_seconds"]["min"],
                "wall_median_s": group["wall_seconds"]["median"],
                "wall_max_s": group["wall_seconds"]["max"],
                "worker_total_s": group["worker_total_median_seconds"],
            }
            row.update({f"{key}_s": value for key, value in
                        group["worker_categories_median_seconds"].items()})
            row.update({f"{key}_worker_s": value for key, value in
                        group["worker_modes_median_seconds"].items()})
            writer.writerow(row)


def plot_wall(path: Path, summary: dict) -> None:
    fig, axis = plt.subplots(figsize=(7.2, 4.5))
    for backend, color, marker, label in (
            ("eigen", "#4C78A8", "o", "Eigen SparseLU"),
            ("lapack", "#E45756", "s", "LAPACK band LU")):
        groups = [group for group in summary["groups"] if group["backend"] == backend]
        medians = [group["wall_seconds"]["median"] for group in groups]
        lower = [m - group["wall_seconds"]["min"] for m, group in zip(medians, groups)]
        upper = [group["wall_seconds"]["max"] - m for m, group in zip(medians, groups)]
        axis.errorbar(FMAX, medians, yerr=[lower, upper], label=label,
                      color=color, marker=marker, linewidth=2, capsize=4)
    axis.set_xlabel("$f_{max}$ [mHz]")
    axis.set_ylabel("Median complete spectra wall time [s]")
    axis.set_xticks(FMAX)
    axis.grid(axis="y", alpha=0.25)
    axis.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_breakdown(path: Path, summary: dict) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.4), sharey=True)
    order = (*CATEGORIES, "unclassified")
    for axis, backend, title in zip(
            axes, BACKENDS, ("Eigen SparseLU", "LAPACK band LU")):
        groups = [group for group in summary["groups"] if group["backend"] == backend]
        bottoms = [0.0] * len(groups)
        for category in order:
            values = [group["worker_categories_median_seconds"][category]
                      for group in groups]
            axis.bar(range(len(groups)), values, bottom=bottoms,
                     color=COLORS[category], label=LABELS[category])
            bottoms = [bottom + value for bottom, value in zip(bottoms, values)]
        axis.set_title(title)
        axis.set_xlabel("$f_{max}$ [mHz]")
        axis.set_xticks(range(len(groups)), [f"{value:g}" for value in FMAX])
        axis.grid(axis="y", alpha=0.2)
    axes[0].set_ylabel("Median aggregate worker time [thread-s]")
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.01))
    fig.suptitle("Adaptive multi-SEM worker-time breakdown\n"
                 "Thread-work seconds; categories are not wall-additive", y=0.99)
    fig.tight_layout(rect=(0, 0.14, 1, 0.94))
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=Path("results/profiling"))
    parser.add_argument("--tag", default="stage_profile_full")
    args = parser.parse_args()
    raw = args.input_dir / f"{args.tag}_profiles.jsonl"
    metadata = args.input_dir / f"{args.tag}_metadata.json"
    rows, campaign_metadata = load_and_validate(raw, metadata)
    summary = summarize(rows, campaign_metadata, raw, metadata)
    summary_path = args.input_dir / "stage_profile_summary.json"
    csv_path = args.input_dir / "stage_profile_summary.csv"
    wall_path = args.input_dir / "stage_profile_wall_time.png"
    breakdown_path = args.input_dir / "stage_profile_worker_breakdown.png"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    write_csv(csv_path, summary)
    plot_wall(wall_path, summary)
    plot_breakdown(breakdown_path, summary)
    print(json.dumps({"records": len(rows), "summary": str(summary_path),
                      "table": str(csv_path), "figures":
                      [str(wall_path), str(breakdown_path)]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
