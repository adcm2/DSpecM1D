#!/usr/bin/env python3
"""Validate a frozen 1024-row trace and write its exact 721-row prefix.

This is a data-format adapter only. It does not alter the input, solver,
observable, window, or any profile/build source.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import math
import os
from pathlib import Path
import sys

EXPECTED_HEADER = ["time_s", "vertical_m", "window_weight"]
FULL_ROWS = 1024
PREFIX_ROWS = 721
SAMPLE_STEP_S = 5.0
WINDOW_START_S = 1000.0
WINDOW_END_S = 1200.0
WINDOW_TOLERANCE = 1.0e-12


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def hann_weight(time_s: float) -> float:
    if time_s < WINDOW_START_S or time_s > WINDOW_END_S:
        return 0.0
    phase = 2.0 * math.pi * (time_s - WINDOW_START_S) / (WINDOW_END_S - WINDOW_START_S)
    return 0.5 * (1.0 - math.cos(phase))


def prefix_trace(input_value: str, output_value: str, receipt_value: str) -> dict:
    input_path = Path(input_value).expanduser().resolve(strict=True)
    output_path = Path(output_value).expanduser().resolve(strict=False)
    receipt_path = Path(receipt_value).expanduser().resolve(strict=False)
    if not input_path.is_file():
        raise ValueError(f"input is not a regular file: {input_path}")
    if output_path == receipt_path:
        raise ValueError("output CSV and receipt paths must differ")
    if input_path in (output_path, receipt_path):
        raise ValueError("input path cannot also be an output path")
    for path in (output_path, receipt_path):
        if not path.parent.is_dir():
            raise ValueError(f"output parent directory does not exist: {path.parent}")
        if path.exists() or path.is_symlink():
            raise FileExistsError(f"refusing existing output: {path}")

    raw = input_path.read_bytes()
    input_before = sha256_bytes(raw)
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError(f"input is not UTF-8 CSV: {exc}") from exc
    physical_lines = text.splitlines(keepends=True)
    if len(physical_lines) != FULL_ROWS + 1:
        raise ValueError(f"expected header plus {FULL_ROWS} physical rows, found {len(physical_lines)} lines")
    try:
        rows = list(csv.reader(physical_lines, strict=True))
    except csv.Error as exc:
        raise ValueError(f"invalid CSV: {exc}") from exc
    if len(rows) != FULL_ROWS + 1:
        raise ValueError(f"expected header plus {FULL_ROWS} CSV rows, found {len(rows) - 1}")
    if rows[0] != EXPECTED_HEADER:
        raise ValueError(f"unexpected CSV header: {rows[0]!r}")

    values: list[tuple[float, float, float]] = []
    for index, row in enumerate(rows[1:]):
        if len(row) != len(EXPECTED_HEADER):
            raise ValueError(f"row {index} has {len(row)} columns; expected 3")
        try:
            time_s, vertical_m, weight = (float(field) for field in row)
        except ValueError as exc:
            raise ValueError(f"row {index} contains a nonnumeric field") from exc
        if not all(math.isfinite(value) for value in (time_s, vertical_m, weight)):
            raise ValueError(f"row {index} contains a nonfinite field")
        expected_time = index * SAMPLE_STEP_S
        if time_s != expected_time:
            raise ValueError(f"row {index} time {time_s!r} != expected {expected_time!r} s")
        expected_weight = hann_weight(time_s)
        if abs(weight - expected_weight) > WINDOW_TOLERANCE:
            raise ValueError(
                f"row {index} Hann weight {weight:.17g} != expected {expected_weight:.17g} "
                f"within {WINDOW_TOLERANCE:g}"
            )
        if index >= PREFIX_ROWS and weight != 0.0:
            raise ValueError(f"removed row {index} has nonzero Hann weight {weight:.17g}")
        values.append((time_s, vertical_m, weight))

    full_energy = math.fsum(weight * vertical_m * vertical_m for _, vertical_m, weight in values)
    prefix_energy = math.fsum(weight * vertical_m * vertical_m
                              for _, vertical_m, weight in values[:PREFIX_ROWS])
    if not (math.isfinite(full_energy) and math.isfinite(prefix_energy)):
        raise ValueError("Hann-weighted energy is nonfinite")
    if full_energy <= 0.0 or prefix_energy <= 0.0:
        raise ValueError("Hann-weighted energy must be positive")
    if full_energy != prefix_energy:
        raise ValueError(f"full/prefix Hann-weighted energy differs: {full_energy!r} != {prefix_energy!r}")

    output_bytes = "".join(physical_lines[:PREFIX_ROWS + 1]).encode("utf-8")
    output_sha = sha256_bytes(output_bytes)
    created_output = False
    created_receipt = False
    try:
        with output_path.open("xb") as output_file:
            created_output = True
            output_file.write(output_bytes)
            output_file.flush()
            os.fsync(output_file.fileno())
        input_after = sha256_bytes(input_path.read_bytes())
        if input_after != input_before:
            raise RuntimeError("input SHA-256 changed while adapter was running")
        receipt = {
            "schema": 1,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "input": str(input_path),
            "input_sha256_before": input_before,
            "input_sha256_after": input_after,
            "output_prefix_csv": str(output_path),
            "output_sha256": output_sha,
            "rows": {"input_data": FULL_ROWS, "output_data": PREFIX_ROWS,
                     "removed_data": FULL_ROWS - PREFIX_ROWS},
            "grid": {"start_s": values[0][0], "end_s": values[-1][0],
                     "step_s": SAMPLE_STEP_S, "prefix_end_s": values[PREFIX_ROWS - 1][0]},
            "window": {"kind": "Hann", "start_s": WINDOW_START_S, "end_s": WINDOW_END_S,
                       "weight_abs_tolerance": WINDOW_TOLERANCE},
            "hann_weighted_energy_m2": {
                "full": full_energy, "prefix": prefix_energy,
                "identical": full_energy == prefix_energy,
            },
            "prefix_bytes_are_original_header_and_rows": True,
        }
        receipt_bytes = (json.dumps(receipt, indent=2, sort_keys=True) + "\n").encode("utf-8")
        with receipt_path.open("xb") as receipt_file:
            created_receipt = True
            receipt_file.write(receipt_bytes)
            receipt_file.flush()
            os.fsync(receipt_file.fileno())
        return receipt
    except Exception:
        if created_receipt:
            receipt_path.unlink(missing_ok=True)
        if created_output:
            output_path.unlink(missing_ok=True)
        raise


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_full_csv", help="full 1024-row DSpec trace CSV")
    parser.add_argument("output_prefix_csv", help="fresh path for the unchanged 721-row CSV prefix")
    parser.add_argument("receipt_json", help="fresh path for the adapter validation receipt")
    args = parser.parse_args(argv)
    try:
        receipt = prefix_trace(args.input_full_csv, args.output_prefix_csv, args.receipt_json)
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"prefix_trace: {exc}", file=sys.stderr)
        return 2
    print(json.dumps({"output_prefix_csv": receipt["output_prefix_csv"],
                      "receipt_json": str(Path(args.receipt_json).expanduser().resolve()),
                      "output_sha256": receipt["output_sha256"],
                      "energy_m2": receipt["hann_weighted_energy_m2"]["full"]},
                     sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
