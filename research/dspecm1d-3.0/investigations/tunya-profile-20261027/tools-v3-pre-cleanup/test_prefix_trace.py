#!/usr/bin/env python3
"""Synthetic-only tests for the data-only prefix adapter."""
from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest

MODULE_PATH = Path(__file__).with_name("prefix_trace.py")
SPEC = importlib.util.spec_from_file_location("prefix_trace", MODULE_PATH)
prefix_trace = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(prefix_trace)


def make_csv() -> bytes:
    lines = ["time_s,vertical_m,window_weight\n"]
    for index in range(prefix_trace.FULL_ROWS):
        time_s = index * prefix_trace.SAMPLE_STEP_S
        weight = prefix_trace.hann_weight(time_s)
        vertical = 2.0 + math.sin(index / 31.0)
        lines.append(f"{time_s:.1f},{vertical:.17g},{weight:.17g}\n")
    return "".join(lines).encode("ascii")


class PrefixTraceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root / "full.csv"
        self.output = self.root / "prefix.csv"
        self.receipt = self.root / "receipt.json"

    def write_source(self, data: bytes | None = None) -> bytes:
        raw = make_csv() if data is None else data
        self.source.write_bytes(raw)
        return raw

    def run_adapter(self) -> dict:
        return prefix_trace.prefix_trace(str(self.source), str(self.output), str(self.receipt))

    def test_prefix_bytes_and_receipt(self) -> None:
        raw = self.write_source()
        receipt = self.run_adapter()
        expected = b"".join(raw.splitlines(keepends=True)[:prefix_trace.PREFIX_ROWS + 1])
        self.assertEqual(self.output.read_bytes(), expected)
        self.assertEqual(self.source.read_bytes(), raw)
        self.assertEqual(receipt["input_sha256_before"], receipt["input_sha256_after"])
        self.assertEqual(receipt["rows"], {"input_data": 1024, "output_data": 721,
                                            "removed_data": 303})
        energies = receipt["hann_weighted_energy_m2"]
        self.assertGreater(energies["full"], 0.0)
        self.assertEqual(energies["full"], energies["prefix"])
        self.assertTrue(json.loads(self.receipt.read_text())["prefix_bytes_are_original_header_and_rows"])

    def test_nonzero_removed_weight_rejected(self) -> None:
        lines = make_csv().decode("ascii").splitlines()
        fields = lines[1 + prefix_trace.PREFIX_ROWS].split(",")
        fields[2] = "1e-15"
        lines[1 + prefix_trace.PREFIX_ROWS] = ",".join(fields)
        self.write_source(("\n".join(lines) + "\n").encode("ascii"))
        with self.assertRaisesRegex(ValueError, "removed row"):
            self.run_adapter()
        self.assertFalse(self.output.exists())
        self.assertFalse(self.receipt.exists())

    def test_wrong_header_rejected(self) -> None:
        raw = make_csv().replace(b"time_s,vertical_m,window_weight", b"time,vertical_m,window_weight", 1)
        self.write_source(raw)
        with self.assertRaisesRegex(ValueError, "unexpected CSV header"):
            self.run_adapter()

    def test_irregular_grid_rejected(self) -> None:
        lines = make_csv().decode("ascii").splitlines()
        fields = lines[17].split(",")
        fields[0] = "81.0"
        lines[17] = ",".join(fields)
        self.write_source(("\n".join(lines) + "\n").encode("ascii"))
        with self.assertRaisesRegex(ValueError, "expected 80.0 s"):
            self.run_adapter()

    def test_nonfinite_value_rejected(self) -> None:
        lines = make_csv().decode("ascii").splitlines()
        fields = lines[9].split(",")
        fields[1] = "nan"
        lines[9] = ",".join(fields)
        self.write_source(("\n".join(lines) + "\n").encode("ascii"))
        with self.assertRaisesRegex(ValueError, "nonfinite"):
            self.run_adapter()

    def test_existing_output_rejected_before_writes(self) -> None:
        self.write_source()
        self.output.write_bytes(b"keep me\n")
        with self.assertRaises(FileExistsError):
            self.run_adapter()
        self.assertEqual(self.output.read_bytes(), b"keep me\n")
        self.assertFalse(self.receipt.exists())


if __name__ == "__main__":
    unittest.main()
