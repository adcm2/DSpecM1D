#!/usr/bin/env python3
"""Validate and summarize the focused Level-1 SparseLU experiment."""

import argparse
import hashlib
import json
import math
import pathlib
import statistics
import subprocess


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def finite(value, name):
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"non-finite {name}")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", type=pathlib.Path, required=True)
    parser.add_argument("--output-dir", type=pathlib.Path, required=True)
    parser.add_argument("--model", type=pathlib.Path, required=True)
    parser.add_argument("--source", type=pathlib.Path, required=True)
    args = parser.parse_args()
    lines = [line.rstrip("\n") for line in args.raw.read_text().splitlines()]
    if len(lines) < 2 or not lines[0].startswith("system\tfamily\tsize"):
        raise ValueError("invalid benchmark header")
    names = lines[0].split("\t")
    records = []
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) != len(names):
            raise ValueError("malformed benchmark row")
        row = dict(zip(names, fields))
        for key in ("analyze_s", "factorize_s", "solve_s", "discrepancy", "residual"):
            row[key] = finite(row[key], key)
            if row[key] < 0:
                raise ValueError(f"negative {key}")
        for key in ("nnz_A", "kl", "ku", "lapack_envelope", "nnz_L", "nnz_U",
                    "info_analyze", "info_factorize", "info_solve"):
            row[key] = int(row[key])
        if (row["info_analyze"] != 0 or row["info_factorize"] != 0 or
                row["info_solve"] != 0):
            raise ValueError("SparseLU factorization/solve failed")
        if row["discrepancy"] >= 1e-9 or row["residual"] >= 1e-10:
            raise ValueError("numerical validation threshold exceeded")
        row["factor_nnz"] = row["nnz_L"] + row["nnz_U"]
        row["fill_ratio"] = row["factor_nnz"] / row["nnz_A"]
        row["factor_to_lapack_envelope"] = row["factor_nnz"] / row["lapack_envelope"]
        records.append(row)
    if len(records) != 44:
        raise ValueError(f"expected 44 rows, got {len(records)}")
    configs = {(r["ordering"], r["symmetric"]) for r in records}
    if configs != {(o, s) for o in ("COLAMD", "NaturalOrdering") for s in ("0", "1")}:
        raise ValueError("configuration coverage is incomplete")
    systems = sorted({r["system"] for r in records})
    if len(systems) != 11:
        raise ValueError("system coverage is incomplete")
    by_config = {}
    for config in sorted(configs):
        subset = [r for r in records if (r["ordering"], r["symmetric"]) == config]
        by_config[f"{config[0]}_symmetric_{config[1]}"] = {
            "analyze_s_total": sum(r["analyze_s"] for r in subset),
            "factorize_s_total": sum(r["factorize_s"] for r in subset),
            "solve_s_total": sum(r["solve_s"] for r in subset),
            "factor_nnz_mean": statistics.mean(r["factor_nnz"] for r in subset),
            "fill_ratio_mean": statistics.mean(r["fill_ratio"] for r in subset),
            "factor_to_lapack_envelope_mean": statistics.mean(
                r["factor_to_lapack_envelope"] for r in subset),
            "kl_values": sorted({r["kl"] for r in subset}),
            "ku_values": sorted({r["ku"] for r in subset}),
            "lapack_envelope_mean": statistics.mean(
                r["lapack_envelope"] for r in subset),
            "max_discrepancy": max(r["discrepancy"] for r in subset),
            "max_residual": max(r["residual"] for r in subset),
        }
    pilot = {}
    for fmax, threads in ((5, 1), (30, 32)):
        key = f"f{fmax}_t{threads}"
        values = {}
        for ordering in ("default", "natural"):
            path = args.output_dir / f"pilot_{key}_{ordering}.txt"
            if not path.exists():
                raise ValueError(f"missing amended pilot artifact: {path}")
            fields = path.read_text().splitlines()[-1].split("\t")
            if len(fields) != 5 or fields[0] != "RESULT":
                raise ValueError(f"malformed pilot artifact: {path}")
            values[ordering] = {
                "seconds": finite(fields[1], "pilot wall time"),
                "rows": int(fields[2]),
                "cols": int(fields[3]),
                "norm": finite(fields[4], "pilot norm"),
                "path": str(path),
            }
        if (values["default"]["rows"], values["default"]["cols"]) != (
            values["natural"]["rows"], values["natural"]["cols"]):
            raise ValueError("pilot output shapes differ")
        values["wall_speedup_default_over_natural"] = (
            values["default"]["seconds"] / values["natural"]["seconds"])
        values["norm_difference"] = abs(
            values["default"]["norm"] - values["natural"]["norm"])
        pilot[key] = values
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary = {
        "experiment": "eigen_sparselu_e1",
        "rows": len(records),
        "systems": systems,
        "configurations": by_config,
        "records": records,
        "conclusion": (
            "NaturalOrdering reduces symbolic analysis time, lowers mean fill, and "
            "redistributes fill across spheroidal systems while having a small numerical "
            "difference relative to COLAMD. NaturalOrdering was checked in the "
            "amended two-point production pilot; this remains an experiment."
        ),
        "production_pilot": pilot,
    }
    summary_path = args.output_dir / "eigen_sparselu_e1_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    metadata = {
        "experiment": "eigen_sparselu_e1",
        "git_head": commit,
        "compiler": "/opt/gcc-15.2.0/bin/g++",
        "focused_kernel_threads": {
            "OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE",
            "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
            "BLIS_NUM_THREADS": "1",
        },
        "pilot_threads": {
            "f5_t1": {"OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE",
                      "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                      "BLIS_NUM_THREADS": "1"},
            "f30_t32": {"OMP_NUM_THREADS": "32", "OMP_DYNAMIC": "FALSE",
                        "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
                        "BLIS_NUM_THREADS": "1"},
        },
        "warmup": 1,
        "measured_repetitions": 1,
        "omega": 0.70,
        "systems": systems,
        "model": {"path": str(args.model.resolve()), "sha256": sha256(args.model)},
        "executable_source": {"path": str(args.source), "sha256": sha256(args.source)},
        "raw": {"path": str(args.raw), "sha256": sha256(args.raw)},
        "summary": {"path": str(summary_path), "sha256": sha256(summary_path)},
        "production_pilot": {
            key: {
                ordering: {"path": value[ordering]["path"],
                           "sha256": sha256(pathlib.Path(value[ordering]["path"]))}
                for ordering in ("default", "natural")
            }
            for key, value in pilot.items()
        },
        "provenance": {
            "focused_kernel": {
                "path": "/tmp/dspecm1d-backend-on-gcc15-final/bin/eigen_sparselu_e1",
                "sha256": sha256(pathlib.Path(
                    "/tmp/dspecm1d-backend-on-gcc15-final/bin/eigen_sparselu_e1")),
                "build_command": "cmake --build /tmp/dspecm1d-backend-on-gcc15-final --target eigen_sparselu_e1 -j4",
                "run_command": "OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 /tmp/dspecm1d-backend-on-gcc15-final/bin/eigen_sparselu_e1 /tmp/dspecm1d-backend-on-gcc15-final/data/models/prem.200.no.txt",
            },
            "pilot_source": {
                "path": str(pathlib.Path("benchmarks/eigen_sparselu_e1_pilot.cpp").resolve()),
                "sha256": sha256(pathlib.Path("benchmarks/eigen_sparselu_e1_pilot.cpp")),
                "note": "same source/configuration compiled twice; one binary uses committed COLAMD and one uses the one-line temporary MultiSem alias override",
            },
            "pilot_executables": {
                "default": {
                    "path": "/tmp/dspecm1d-backend-on-gcc15-final/bin/eigen_sparselu_e1_pilot",
                    "sha256": sha256(pathlib.Path("/tmp/dspecm1d-backend-on-gcc15-final/bin/eigen_sparselu_e1_pilot")),
                },
                "natural": {
                    "path": "/tmp/eigen_sparselu_e1_natural_pilot",
                    "sha256": sha256(pathlib.Path("/tmp/eigen_sparselu_e1_natural_pilot")),
                },
            },
            "pilot_build_method": {
                "default": "cmake --build /tmp/dspecm1d-backend-on-gcc15-final --target eigen_sparselu_e1_pilot -j4",
                "natural": "copy repository DSpecM1D to /tmp/e1_override; replace only Eigen::COLAMDOrdering<int> with Eigen::NaturalOrdering<int> in FullSpecMultiSem.h; compile benchmarks/eigen_sparselu_e1_pilot.cpp with GCC 15.2, the same Release flags/includes/libs, output /tmp/eigen_sparselu_e1_natural_pilot",
                "natural_command": "/opt/gcc-15.2.0/bin/g++ -DDSPECM1D_ENABLE_LAPACK_BAND_SOLVER -DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP -I/tmp/e1_override -I/tmp/dspecm1d-backend-on-gcc15-final/benchmarks -I/space/adcm2/CleanBuilds/DSpecM1D -I/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/eigen3-src -I/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/fftwpp-src -I/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/numericconcepts-src/include -I/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/gaussquad-src -I/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/gshtrans-src -I/home/adcm2/Documents/c++/DSpecM3D/build-docs/_deps/interpolation-src -O3 -DNDEBUG -fopenmp -std=gnu++23 benchmarks/eigen_sparselu_e1_pilot.cpp -o /tmp/eigen_sparselu_e1_natural_pilot /usr/lib/x86_64-linux-gnu/libfftw3.so /usr/lib/x86_64-linux-gnu/libfftw3f.so /usr/lib/x86_64-linux-gnu/libfftw3l.so /opt/gcc-15.2.0/lib64/libgomp.so /lib/x86_64-linux-gnu/libpthread.a /usr/lib/x86_64-linux-gnu/liblapacke.so -lm -ldl /usr/lib/x86_64-linux-gnu/libopenblas.so",
            },
            "pilot_parameters": {
                "f5_t1": {"path": "/tmp/dspecm1d-backend-on-gcc15-final/data/params/e1_pilot_f5.txt", "sha256": sha256(pathlib.Path("/tmp/dspecm1d-backend-on-gcc15-final/data/params/e1_pilot_f5.txt"))},
                "f30_t32": {"path": "/tmp/dspecm1d-backend-on-gcc15-final/data/params/e1_pilot_f30.txt", "sha256": sha256(pathlib.Path("/tmp/dspecm1d-backend-on-gcc15-final/data/params/e1_pilot_f30.txt"))},
            },
            "pilot_commands": {
                "f5_t1_default": "cd /tmp/dspecm1d-backend-on-gcc15-final && OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 ./bin/eigen_sparselu_e1_pilot data/params/e1_pilot_f5.txt",
                "f5_t1_natural": "cd /tmp/dspecm1d-backend-on-gcc15-final && OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 /tmp/eigen_sparselu_e1_natural_pilot data/params/e1_pilot_f5.txt",
                "f30_t32_default": "cd /tmp/dspecm1d-backend-on-gcc15-final && OMP_NUM_THREADS=32 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 ./bin/eigen_sparselu_e1_pilot data/params/e1_pilot_f30.txt",
                "f30_t32_natural": "cd /tmp/dspecm1d-backend-on-gcc15-final && OMP_NUM_THREADS=32 OMP_DYNAMIC=FALSE OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 /tmp/eigen_sparselu_e1_natural_pilot data/params/e1_pilot_f30.txt",
            },
            "natural_override": {
                "path": str((args.output_dir / "natural_ordering_override.diff").resolve()),
                "sha256": sha256(args.output_dir / "natural_ordering_override.diff"),
                "scope": "exactly one alias line in FullSpecMultiSem.h; no committed production file was changed",
            },
            "analyzer": {"path": str(pathlib.Path(__file__).resolve()), "sha256": sha256(pathlib.Path(__file__))},
        },
    }
    (args.output_dir / "eigen_sparselu_e1_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
