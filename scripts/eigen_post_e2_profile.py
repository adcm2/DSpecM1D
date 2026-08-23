#!/usr/bin/env python3
"""Run the single post-E2 Eigen SparseLU profiling reassessment point."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replace_after_marker(lines: list[str], marker: str, value: str) -> None:
    for index, line in enumerate(lines):
        if marker in line:
            for next_index in range(index + 1, len(lines)):
                if lines[next_index].strip() and not lines[next_index].lstrip().startswith("#"):
                    ending = "\n" if lines[next_index].endswith("\n") else ""
                    lines[next_index] = f"  {value}{ending}"
                    return
    raise RuntimeError(f"parameter marker not found: {marker}")


def git(repo: Path, *args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=repo, text=True).strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--fmax", type=float, default=30.0)
    parser.add_argument("--threads", type=int, default=32)
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[1]
    if args.fmax != 30.0 or args.threads != 32:
        raise RuntimeError("this reassessment is fixed at fmax=30, threads=32")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    parameter = (args.build / "data/params/eigen_post_e2_profile_f30.txt").resolve()
    parameter.parent.mkdir(parents=True, exist_ok=True)
    lines = (repo / "data/params/ex3.txt").read_text().splitlines(keepends=True)
    replace_after_marker(lines, "# fmax", "30")
    replace_after_marker(lines, "# f21 filter", "27")
    replace_after_marker(lines, "# f22 filter", "30")
    parameter.write_text("".join(lines))

    environment = os.environ.copy()
    environment.update({
        "OMP_NUM_THREADS": "32", "OMP_DYNAMIC": "FALSE",
        "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
        "BLIS_NUM_THREADS": "1",
    })
    executable = (args.build / "bin/stage_profile_solver").resolve()
    command = [str(executable), str(parameter), "eigen", "1", "1"]
    completed = subprocess.run(command, cwd=args.build, env=environment,
                               check=True, capture_output=True, text=True)
    records = [line for line in completed.stdout.splitlines()
               if line.startswith("PROFILE\t")]
    if len(records) != 1:
        raise RuntimeError(f"expected one measured profile record, got {len(records)}")
    fields = records[0].split("\t", 6)
    if len(fields) != 7 or fields[1] != "eigen":
        raise RuntimeError("malformed profile record")
    profile = json.loads(fields[6])
    if not profile["wall_seconds"] > 0 or not all(
            value >= 0 for value in profile["timings_seconds"]["all"].values()):
        raise RuntimeError("invalid profile timings")

    raw = output / "eigen_post_e2_profile_raw.tsv"
    raw.write_text("fmax_mHz\tbackend\trepetition\trows\tcols\tnorm\tprofile_json\n"
                   + "30\teigen\t1\t" + "\t".join(fields[3:6]) + "\t"
                   + fields[6] + "\n")
    counts = profile["counts"]
    worker = profile["accounting"]["workers"]["all"]
    timings = profile["timings_seconds"]["all"]
    categories = {key: value for key, value in timings.items()
                  if key not in ("total_spectra", "unclassified")}
    summary = {
        "campaign": "Post-E2 Eigen SparseLU profiling reassessment",
        "status": "complete", "fmax_mHz": 30.0, "threads": 32,
        "one_untimed_warmup": True, "measured_repetitions": 1,
        "raw_sha256": sha256(raw),
        "profile": {"wall_seconds": profile["wall_seconds"],
                    "worker_total_seconds": worker["total"],
                    "worker_categories_seconds": timings,
                    "worker_category_share_percent": {
                        key: 100.0 * value / worker["total"]
                        for key, value in categories.items()},
                    "factorization_to_solve_ratio": (
                        timings["factorization"] / timings["solve"]),
                    "worker_unclassified_seconds": worker["unclassified"],
                    "counts": counts,
                    "output_shape": [int(fields[3]), int(fields[4])],
                    "norm": float(fields[5])},
    }
    summary_path = output / "eigen_post_e2_profile_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    source_paths = ["DSpecM1D/src/Profiling.h", "DSpecM1D/src/FullSpecMultiSem.h",
                    "benchmarks/stage_profile_solver.cpp", __file__.replace(str(repo) + "/", "")]
    metadata = {
        "campaign": summary["campaign"], "status": "complete",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "git": {"revision": git(repo, "rev-parse", "HEAD"),
                "status_before_run": git(repo, "status", "--short")},
        "machine": {"hostname": platform.node()},
        "configuration": "ex3/PREM-style type-4 all-mode attenuation=1, Eigen SparseLU NaturalOrdering",
        "environment": {key: environment[key] for key in (
            "OMP_NUM_THREADS", "OMP_DYNAMIC", "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS", "BLIS_NUM_THREADS")},
        "build": str(args.build.resolve()), "executable": str(executable),
        "executable_sha256": sha256(executable), "command": " ".join(command),
        "parameter_template": str((repo / "data/params/ex3.txt").resolve()),
        "parameter_template_sha256": sha256(repo / "data/params/ex3.txt"),
        "parameter_path_provenance": (
            "generated from the retained repository ex3 template into the build "
            "data/params directory so InputParametersNew resolves data/models"),
        "parameter_file": str(parameter), "parameter_sha256": sha256(parameter),
        "raw_file": str(raw), "summary_file": str(summary_path),
        "source_sha256": {path: sha256(repo / path) for path in source_paths},
        "profile_counts": counts, "profile_wall_seconds": profile["wall_seconds"],
    }
    metadata_path = output / "eigen_post_e2_profile_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    summary["metadata_sha256"] = sha256(metadata_path)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
