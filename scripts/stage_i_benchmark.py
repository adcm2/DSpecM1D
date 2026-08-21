#!/usr/bin/env python3
"""Run the Stage I Eigen-vs-LAPACK complete-spectra benchmark.

The executable performs one untimed warm-up followed by all timed complete
spectra() calls for one backend/fmax.  This driver alternates which backend is
launched first at successive fmax points, retaining one raw TSV row per
measurement plus JSON metadata.  PR1's LAPACK path is intentionally documented
as: Eigen sparse matrix -> LAPACK band packing -> band LU factorization/solve.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
from pathlib import Path
import platform
import subprocess
import sys


FINAL_FMAX = [5.0, 10.0, 20.0, 35.0, 55.0]


def command_text(command: list[str]) -> str:
    return " ".join(command)


def run_text(command: list[str], *, cwd: Path | None = None) -> str:
    return subprocess.check_output(command, cwd=cwd, text=True).strip()


def make_parameter_file(template: Path, destination: Path, fmax: float,
                        lmax: int | None) -> None:
    text = template.read_text()
    lines = text.splitlines(keepends=True)

    def replace_value_after(marker: str, value: str) -> None:
        for index, line in enumerate(lines):
            if marker in line:
                for next_index in range(index + 1, len(lines)):
                    if lines[next_index].strip() and not lines[next_index].lstrip().startswith("#"):
                        ending = "\n" if lines[next_index].endswith("\n") else ""
                        lines[next_index] = f"  {value}{ending}"
                        return
        raise RuntimeError(f"parameter marker not found: {marker}")

    replace_value_after("# fmax", f"{fmax:.12g}")
    replace_value_after("# f21 filter", f"{max(0.1, 0.9 * fmax):.12g}")
    replace_value_after("# f22 filter", f"{fmax:.12g}")
    if lmax is not None:
        replace_value_after("# lmax", str(lmax))
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("".join(lines))


def git_metadata(repo: Path) -> dict[str, str]:
    return {
        "git_revision": run_text(["git", "rev-parse", "HEAD"], cwd=repo),
        "git_status": run_text(["git", "status", "--short"], cwd=repo),
    }


def compiler_metadata() -> dict[str, str]:
    compiler = "/opt/gcc-15.2.0/bin/g++"
    return {
        "compiler": compiler,
        "compiler_version": run_text([compiler, "--version"]).splitlines()[0],
        "cmake_version": run_text(["cmake", "--version"]).splitlines()[0],
    }


def build_metadata(build: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in (build / "CMakeCache.txt").read_text().splitlines():
        if line.startswith(("CMAKE_CXX_FLAGS:STRING=",
                            "CMAKE_CXX_FLAGS_RELEASE:STRING=",
                            "CMAKE_BUILD_TYPE:STRING=",
                            "CMAKE_CXX_STANDARD:STRING=")):
            key, value = line.split("=", 1)
            values[key.split(":", 1)[0]] = value
    return values


def cpu_metadata() -> dict[str, str]:
    details: dict[str, str] = {}
    try:
        for line in run_text(["lscpu"]).splitlines():
            if ":" in line:
                key, value = line.split(":", 1)
                if key.strip() in {"Model name", "CPU(s)", "Thread(s) per core",
                                   "Core(s) per socket", "Socket(s)", "Flags"}:
                    details[key.strip()] = value.strip()
    except (OSError, subprocess.CalledProcessError):
        details["platform"] = platform.platform()
    return details


def invoke(executable: Path, parameter_file: Path, backend: str,
           first_rep: int, repetitions: int,
           env: dict[str, str]) -> list[tuple[int, float, int, int, float]]:
    command = [str(executable), str(parameter_file), backend, str(first_rep),
               str(repetitions)]
    output = subprocess.check_output(command, text=True, env=env)
    rows = [line for line in output.splitlines() if line.startswith("RESULT\t")]
    if len(rows) != repetitions:
        raise RuntimeError(f"unexpected benchmark output from {command}: {output!r}")
    parsed = []
    for row in rows:
        _, actual_backend, actual_rep, seconds, nrows, ncols, norm = row.split("\t")
        if actual_backend != backend:
            raise RuntimeError(f"result metadata mismatch: {row}")
        parsed.append((int(actual_rep), float(seconds), int(nrows), int(ncols),
                       float(norm)))
    return parsed


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--eigen-build", type=Path, required=True)
    parser.add_argument("--lapack-build", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--fmax", type=float, nargs="+", default=FINAL_FMAX)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--lmax", type=int, default=None)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--tag", default="stage_i")
    args = parser.parse_args()

    if args.repetitions < 1 or args.threads < 1:
        parser.error("--repetitions and --threads must be positive")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw_path = args.output_dir / f"{args.tag}_timings.tsv"
    metadata_path = args.output_dir / f"{args.tag}_metadata.json"
    template = args.repo / "data/params/ex3.txt"
    env = os.environ.copy()
    env.update({
        "OMP_NUM_THREADS": str(args.threads),
        "OMP_DYNAMIC": "FALSE",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "BLIS_NUM_THREADS": "1",
    })

    metadata = {
        "campaign": "DSpecM1D Stage I PR1 solver benchmark",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "configuration": "ex3/PREM-style type-4 all-mode attenuation=1",
        "measurement": "complete preferred adaptive multi-SEM SparseFSpec::spectra(paramsNew)",
        "lapack_pr1_pipeline": "Eigen sparse frequency matrix -> LAPACK band packing -> zgbtrf/zgbtrs",
        "grid_fmax_mHz": args.fmax,
        "repetitions": args.repetitions,
        "one_untimed_warmup_per_backend_fmax": True,
        "backend_order": "alternating first backend per fmax",
        "threads": {key: env[key] for key in (
            "OMP_NUM_THREADS", "OMP_DYNAMIC", "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS", "BLIS_NUM_THREADS")},
        "eigen_build": str(args.eigen_build.resolve()),
        "lapack_build": str(args.lapack_build.resolve()),
        "template": str(template.resolve()),
        "machine": {"hostname": platform.node(), **cpu_metadata()},
        "eigen_cmake_flags": build_metadata(args.eigen_build),
        "lapack_cmake_flags": build_metadata(args.lapack_build),
        **git_metadata(args.repo),
        **compiler_metadata(),
        "command": command_text(sys.argv),
    }

    with raw_path.open("w") as raw:
        raw.write("fmax_mHz\tbackend\trepetition\tseconds\trows\tcols\tnorm\torder\n")
        order = 0
        for fmax_index, fmax in enumerate(args.fmax):
            backend_order = ["eigen", "lapack"] if fmax_index % 2 == 0 else ["lapack", "eigen"]
            parameter_files: dict[str, Path] = {}
            for backend, build in (("eigen", args.eigen_build), ("lapack", args.lapack_build)):
                parameter = build / "data/params" / f"stage_i_{args.tag}_{fmax:g}.txt"
                make_parameter_file(template, parameter, fmax, args.lmax)
                parameter_files[backend] = parameter
            # Each backend/fmax process performs exactly one warm-up, then all
            # measured repetitions. Alternate which backend launches first at
            # each fmax to keep launch/order effects symmetric.
            for backend in backend_order:
                build = args.eigen_build if backend == "eigen" else args.lapack_build
                executable = build / "bin/stage_i_solver_benchmark"
                measured = invoke(executable, parameter_files[backend], backend,
                                  1, args.repetitions, env)
                for repetition, seconds, rows, cols, norm in measured:
                    order += 1
                    raw.write(f"{fmax:.12g}\t{backend}\t{repetition}\t{seconds:.17g}\t"
                              f"{rows}\t{cols}\t{norm:.17g}\t{order}\n")
                    raw.flush()
            print(f"completed fmax={fmax:g} mHz", file=sys.stderr, flush=True)

    metadata["raw_timing_file"] = str(raw_path.resolve())
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"raw": str(raw_path), "metadata": str(metadata_path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
