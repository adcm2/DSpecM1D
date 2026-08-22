#!/usr/bin/env python3
"""Run aggregate adaptive multi-SEM production profiling.

Each executable call performs one warm-up and the requested measured complete
spectra() calls.  Profile records are appended and flushed after every call so
completed backend/frequency points survive an interrupted campaign.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
from typing import Iterator


def run_text(command: list[str], cwd: Path | None = None) -> str:
    return subprocess.check_output(command, cwd=cwd, text=True).strip()


def git_metadata(repo: Path) -> dict[str, str]:
    return {
        "revision": run_text(["git", "rev-parse", "HEAD"], cwd=repo),
        "status_short": run_text(["git", "status", "--short"], cwd=repo),
    }


def cpu_metadata() -> dict[str, str]:
    wanted = {"Model name", "CPU(s)", "Socket(s)", "Core(s) per socket",
              "Thread(s) per core"}
    result: dict[str, str] = {"hostname": platform.node()}
    for line in run_text(["lscpu"]).splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        if key.strip() in wanted:
            result[key.strip()] = value.strip()
    return result


def build_metadata(build: Path) -> dict[str, str]:
    wanted = {"CMAKE_BUILD_TYPE", "CMAKE_CXX_FLAGS",
              "CMAKE_CXX_FLAGS_RELEASE", "CMAKE_CXX_STANDARD",
              "DSPECM1D_ENABLE_LAPACK_BAND_SOLVER"}
    result: dict[str, str] = {}
    for line in (build / "CMakeCache.txt").read_text().splitlines():
        if line.startswith(("#", "//")) or "=" not in line or ":" not in line:
            continue
        key_with_type, value = line.split("=", 1)
        key = key_with_type.split(":", 1)[0]
        if key in wanted:
            result[key] = value
    result.setdefault("CMAKE_CXX_STANDARD", "23")
    result["effective_release_flags"] = " ".join(filter(None, (
        result.get("CMAKE_CXX_FLAGS", ""),
        result.get("CMAKE_CXX_FLAGS_RELEASE", ""),
        "-fopenmp", "-std=gnu++23")))
    return result


def profiling_source_hashes(repo: Path) -> dict[str, str]:
    paths = [
        "DSpecM1D/src/Profiling.h",
        "DSpecM1D/src/FullSpecMultiSem.h",
        "DSpecM1D/src/LapackBandSolver.h",
        "DSpecM1D/src/LapackBandStorage.h",
        "benchmarks/stage_profile_solver.cpp",
        "scripts/stage_profile.py",
    ]
    return {
        path: hashlib.sha256((repo / path).read_bytes()).hexdigest()
        for path in paths
    }


def make_parameter_file(template: Path, destination: Path, fmax: float,
                        lmax: int | None) -> None:
    lines = template.read_text().splitlines(keepends=True)

    def replace(marker: str, value: str) -> None:
        for index, line in enumerate(lines):
            if marker in line:
                for next_index in range(index + 1, len(lines)):
                    if lines[next_index].strip() and not lines[next_index].lstrip().startswith("#"):
                        ending = "\n" if lines[next_index].endswith("\n") else ""
                        lines[next_index] = f"  {value}{ending}"
                        return
        raise RuntimeError(f"parameter marker not found: {marker}")

    replace("# fmax", f"{fmax:.12g}")
    replace("# f21 filter", f"{max(0.1, 0.9 * fmax):.12g}")
    replace("# f22 filter", f"{fmax:.12g}")
    if lmax is not None:
        replace("# lmax", str(lmax))
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("".join(lines))


def metadata(repo: Path, args: argparse.Namespace, env: dict[str, str]) -> dict:
    compiler = "/opt/gcc-15.2.0/bin/g++"
    return {
        "campaign": "DSpecM1D profiling stage",
        "status": "running",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "configuration": "ex3/PREM-style type-4 all-mode attenuation=1",
        "measurement": "complete preferred adaptive multi-SEM spectra(paramsNew)",
        "lapack_pr1_pipeline": "Eigen sparse matrix -> LAPACK band packing -> zgbtrf/zgbtrs",
        "grid_fmax_mHz": args.fmax,
        "repetitions": args.repetitions,
        "one_untimed_warmup_per_backend_fmax": True,
        "backend_order": "alternating first backend per fmax",
        "lmax": args.lmax if args.lmax is not None else 750,
        "threads": {key: env[key] for key in (
            "OMP_NUM_THREADS", "OMP_DYNAMIC", "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS", "BLIS_NUM_THREADS")},
        "eigen_build": str(args.eigen_build.resolve()),
        "lapack_build": str(args.lapack_build.resolve()),
        "git": git_metadata(repo),
        "machine": cpu_metadata(),
        "compiler": compiler,
        "compiler_version": run_text([compiler, "--version"]).splitlines()[0],
        "eigen_cmake": build_metadata(args.eigen_build),
        "lapack_cmake": build_metadata(args.lapack_build),
        "profiling_source_sha256": profiling_source_hashes(repo),
        "command": " ".join(sys.argv),
        "completed_points": [],
        "completed_backends": [],
    }


def write_metadata(path: Path, data: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def invoke(executable: Path, parameter: Path, backend: str, repetitions: int,
           env: dict[str, str]) -> Iterator[dict]:
    process = subprocess.Popen(
        [str(executable), str(parameter), backend, "1", str(repetitions)],
        stdout=subprocess.PIPE, text=True, bufsize=1, env=env)
    assert process.stdout is not None
    records = 0
    for line in process.stdout:
        line = line.rstrip("\n")
        if not line.startswith("PROFILE\t"):
            continue
        fields = line.split("\t", 6)
        if len(fields) != 7 or fields[1] != backend:
            raise RuntimeError(f"invalid profile record: {line!r}")
        records += 1
        yield {
            "backend": backend,
            "repetition": int(fields[2]),
            "rows": int(fields[3]),
            "cols": int(fields[4]),
            "norm": float(fields[5]),
            "profile": json.loads(fields[6]),
        }
    returncode = process.wait()
    if returncode != 0:
        raise subprocess.CalledProcessError(returncode, process.args)
    if records != repetitions:
        raise RuntimeError(f"expected {repetitions} profile records, got {records}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--eigen-build", type=Path, required=True)
    parser.add_argument("--lapack-build", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--fmax", type=float, nargs="+", default=[20.0, 35.0, 55.0])
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--lmax", type=int, default=None)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--tag", default="stage_profile_campaign")
    args = parser.parse_args()
    if args.repetitions < 1 or args.threads < 1:
        parser.error("--repetitions and --threads must be positive")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw_path = args.output_dir / f"{args.tag}_profiles.jsonl"
    metadata_path = args.output_dir / f"{args.tag}_metadata.json"
    template = args.repo / "data/params/ex3.txt"
    env = os.environ.copy()
    env.update({"OMP_NUM_THREADS": str(args.threads), "OMP_DYNAMIC": "FALSE",
                "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                "BLIS_NUM_THREADS": "1"})
    info = metadata(args.repo, args, env)
    info["raw_profile_file"] = str(raw_path.resolve())
    write_metadata(metadata_path, info)

    with raw_path.open("a") as raw:
        if raw.tell() == 0:
            raw.write("fmax_mHz\tbackend\trepetition\trows\tcols\tnorm\tprofile_json\n")
            raw.flush()
        for fmax_index, fmax in enumerate(args.fmax):
            backend_order = (["eigen", "lapack"] if fmax_index % 2 == 0
                             else ["lapack", "eigen"])
            parameters = {}
            for backend, build in (("eigen", args.eigen_build), ("lapack", args.lapack_build)):
                parameter = build / "data/params" / f"{args.tag}_{fmax:g}.txt"
                make_parameter_file(template, parameter, fmax, args.lmax)
                parameters[backend] = parameter
            for backend in backend_order:
                build = args.eigen_build if backend == "eigen" else args.lapack_build
                executable = build / "bin/stage_profile_solver"
                for record in invoke(executable, parameters[backend], backend,
                                     args.repetitions, env):
                    raw.write(f"{fmax:.12g}\t{backend}\t{record['repetition']}\t"
                              f"{record['rows']}\t{record['cols']}\t{record['norm']:.17g}\t"
                              f"{json.dumps(record['profile'], separators=(',', ':'))}\n")
                    raw.flush()
                info["completed_backends"].append(
                    {"fmax_mHz": fmax, "backend": backend})
                write_metadata(metadata_path, info)
            info["completed_points"].append(fmax)
            write_metadata(metadata_path, info)
            print(f"completed fmax={fmax:g} mHz", file=sys.stderr, flush=True)
    info["status"] = "complete"
    write_metadata(metadata_path, info)
    print(json.dumps({"raw": str(raw_path), "metadata": str(metadata_path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
