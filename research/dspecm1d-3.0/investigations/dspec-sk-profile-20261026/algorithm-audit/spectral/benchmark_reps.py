#!/usr/bin/env python3
"""Run isolated frozen spectral-consumer cases and retain JSONL outputs."""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import time


def main() -> int:
    if len(sys.argv) != 4:
        raise SystemExit("usage: benchmark_reps.py FROZEN_BINARY PACKET_ROOT OUTPUT_ROOT")
    binary, packet_root, output_root = map(Path, sys.argv[1:])
    output_root.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    summary = []
    for degree in (3, 6, 12):
        for frequencies in (3, 48):
            case = f"l{degree}-f{frequencies}"
            packet = packet_root / case / "solver_packet.bin"
            for repetition in (1, 2, 3):
                output_path = output_root / f"{case}-rep{repetition}.jsonl"
                command = [str(binary), str(packet), str(degree)]
                started = time.perf_counter()
                completed = subprocess.run(command, env=env, text=True,
                                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                wall = time.perf_counter() - started
                if completed.returncode:
                    raise RuntimeError(f"{case} rep{repetition} failed: {completed.stderr}")
                output_path.write_text(completed.stdout)
                records = [json.loads(line) for line in completed.stdout.splitlines() if line.strip()]
                aggregate = next(record for record in records if record.get("aggregate"))
                meta = {"case": case, "repetition": repetition, "wall_seconds": wall,
                        "command": command, "aggregate": aggregate,
                        "stderr": completed.stderr}
                (output_root / f"{case}-rep{repetition}.meta.json").write_text(
                    json.dumps(meta, indent=2) + "\n")
                summary.append(meta)
                print(json.dumps({"case": case, "repetition": repetition,
                                  "wall_seconds": wall,
                                  "plan_seconds": aggregate["global_plan_prepare_seconds"],
                                  "direct_synthesis_seconds": aggregate["direct_synthesis_seconds_sum"],
                                  "covariance_seconds": aggregate["covariance_formation_seconds_sum"] +
                                                        aggregate["covariance_evaluation_seconds_sum"],
                                  "recoupling_fft_seconds": aggregate["recoupling_seconds_sum"] +
                                                            aggregate["kernel_assembly_seconds_sum"] +
                                                            aggregate["fft_synthesis_seconds_sum"],
                                  "field_l2": aggregate["field_reconstruction_l2_max"],
                                  "kernel_l2": aggregate["frequency_covariance_kernel_l2_max"],
                                  "gaunt_l2": aggregate["gaunt_fft_l2_max"]}))
    (output_root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
