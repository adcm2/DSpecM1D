from __future__ import annotations

import os
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

REPO = Path(__file__).resolve().parents[1]
YSPECPP = REPO / "benchmarks" / "yspecpp"
sys.path.insert(0, str(YSPECPP))
import run_paper as paper


def snapshot(root: Path) -> dict[str, tuple[str, bytes | None]]:
    return {
        path.relative_to(root).as_posix(): (
            "directory" if path.is_dir() else "file",
            None if path.is_dir() else path.read_bytes(),
        )
        for path in sorted(root.rglob("*"))
    }


class PaperTransformTests(unittest.TestCase):
    def test_processed_fig1_spectrum_undoes_damping_once(self):
        nfft, dt, epsilon = 128, 0.25, 0.017
        df = 1 / (nfft * dt)
        time = np.arange(nfft) * dt
        frequency = np.arange(nfft // 2 + 1) * df
        physical = (
            0.7 * np.sin(2 * np.pi * frequency[1] * time + 0.2)
            + 0.3 * np.cos(2 * np.pi * frequency[2] * time - 0.4)
            + 0.2 * np.sin(2 * np.pi * frequency[3] * time + 0.9)
        )
        corners = [10.0, 20.0, 80.0, 110.0]
        damped = np.fft.rfft(physical * np.exp(-epsilon * time)) / (nfft * df)
        spectra = {(index, 0): value for index, value in enumerate(damped)}
        metadata = {
            "samples": nfft,
            "sample_interval_s": dt,
            "frequency_spacing_hz": df,
            "damping_epsilon_rad_s": epsilon,
            "output_samples": nfft,
        }

        actual_frequency, actual = paper.processed_spectrum(
            spectra, metadata, corners
        )
        duration = (nfft - 1) * dt
        pre_windowed = np.fft.rfft(
            physical * paper.hann_ends(time, duration, 0.05), n=nfft
        )
        filtered_time = np.fft.irfft(
            pre_windowed * paper.half_cosine(frequency, corners), n=nfft
        ) * (nfft * df)
        expected = np.abs(np.fft.rfft(
            filtered_time * paper.hann_ends(time, duration, 0.5), n=nfft
        ))

        np.testing.assert_array_equal(actual_frequency, frequency)
        np.testing.assert_allclose(
            actual / actual.max(),
            expected / expected.max(),
            rtol=2e-12,
            atol=2e-12,
        )
        twice_damped = filtered_time * np.exp(epsilon * time)
        wrong = np.abs(np.fft.rfft(
            twice_damped * paper.hann_ends(time, duration, 0.5), n=nfft
        ))
        self.assertGreater(
            float(np.mean(abs(actual / actual.max() - wrong / wrong.max()))),
            1e-3,
        )

    def test_paper_case_inventory_is_exactly_figures_one_to_four(self):
        self.assertEqual(set(paper.CASES), {"fig1", "fig2", "fig3", "fig4"})
        self.assertEqual(
            paper.CASES["fig4"]["dspec"]["maxstep_earth_radii"], 0.004
        )
        self.assertEqual(
            paper.CASES["fig2"]["criteria"]["mean_gate_percent_exclusive"], 1.0
        )

    def test_legacy_provenance_identity_for_all_four_cases(self):
        for case_id, case in paper.CASES.items():
            provenance = {
                "source_parameter_sha256": case["source_parameter_sha256"],
            }
            if case_id == "fig2":
                provenance.pop("source_parameter_sha256")
                provenance["moment_tensor_input_sha256"] = case[
                    "source_parameter_sha256"
                ]
            else:
                provenance["case"] = case_id
                provenance["paper_figure"] = case["paper_figure"]
            paper.validate_saved_identity(case_id, provenance)

    def test_wrong_saved_case_fails_before_creating_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            saved = root / "fig3" / "full"
            results = saved / "results"
            results.mkdir(parents=True)
            (results / "dspec.spectra.csv").write_text("i,j,real,imag\n")
            (results / "yspecpp.spectra.csv").write_text("i,j,real,imag\n")
            fig3 = paper.CASES["fig3"]
            (saved / "provenance.json").write_text(json.dumps({
                "case": "fig3",
                "paper_figure": 3,
                "source_parameter_sha256": fig3["source_parameter_sha256"],
            }))
            output = root / "requested-output"
            sentinel = root / "outside-sentinel"
            sentinel.write_text("unchanged")
            before = snapshot(root)

            result = subprocess.run(
                [
                    sys.executable,
                    "-B",
                    str(Path(paper.__file__)),
                    "--case",
                    "fig4",
                    "--saved-run",
                    str(saved),
                    "--output",
                    str(output),
                ],
                cwd=root,
                env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertIn("source hash does not match fig4", result.stderr)
            self.assertFalse(output.exists())
            self.assertEqual(sentinel.read_text(), "unchanged")
            current = snapshot(root)
            del current["outside-sentinel"]
            expected = before.copy()
            del expected["outside-sentinel"]
            self.assertEqual(current, expected)

    def test_populated_output_stops_before_any_subprocess_or_write(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            hooks = root / "hooks"
            hooks.mkdir()
            (hooks / "sitecustomize.py").write_text(
                "import os\n"
                "from pathlib import Path\n"
                "import subprocess\n"
                "def blocked(*args, **kwargs):\n"
                "    Path(os.environ['PROCESS_SENTINEL']).write_text('called')\n"
                "    raise RuntimeError('subprocess sentinel called')\n"
                "subprocess.run = blocked\n"
                "subprocess.check_output = blocked\n"
                "subprocess.Popen = blocked\n",
                encoding="utf-8",
            )
            output = root / "populated-output"
            (output / "nested").mkdir(parents=True)
            (output / ".hidden-provenance").write_text("keep hidden")
            (output / "nested" / "input.json").write_text('{"frozen": true}\n')
            before = snapshot(output)
            process_sentinel = root / "process-called"
            environment = os.environ.copy()
            environment["PYTHONDONTWRITEBYTECODE"] = "1"
            environment["PYTHONPATH"] = str(hooks)
            environment["PROCESS_SENTINEL"] = str(process_sentinel)

            result = subprocess.run(
                [
                    sys.executable,
                    "-B",
                    str(YSPECPP / "run_paper.py"),
                    "--case",
                    "fig2",
                    "--saved-run",
                    str(root / "missing-saved-run"),
                    "--output",
                    str(output),
                ],
                cwd=root,
                env=environment,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertIn("output directory must be empty", result.stderr)
            self.assertFalse(process_sentinel.exists())
            self.assertEqual(snapshot(output), before)


if __name__ == "__main__":
    unittest.main()
