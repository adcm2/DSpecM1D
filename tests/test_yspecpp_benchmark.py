from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

REPO = Path(__file__).resolve().parents[1]
YSPECPP = REPO / "benchmarks" / "yspecpp"
sys.path.insert(0, str(YSPECPP))

from run_paper_suite import legacy_processed_amplitude_spectrum


def _hann_ends(time: np.ndarray, duration: float, fraction: float) -> np.ndarray:
    edge = fraction * duration
    window = np.ones(time.size)
    if edge <= 0:
        return window
    left = time < edge
    right = time > duration - edge
    window[left] = 0.5 * (1.0 - np.cos(np.pi * time[left] / edge))
    window[right] = 0.5 * (
        1.0 - np.cos(np.pi * (duration - time[right]) / edge)
    )
    window[time >= duration] = 0.0
    return window


def _half_cosine(frequency: np.ndarray,
                 corners_mhz: list[float]) -> np.ndarray:
    f1, f2, f3, f4 = np.asarray(corners_mhz, dtype=float) / 1000.0
    weights = np.zeros(frequency.shape)
    rising = (frequency >= f1) & (frequency < f2)
    flat = (frequency >= f2) & (frequency < f3)
    falling = (frequency >= f3) & (frequency < f4)
    weights[rising] = 0.5 * (
        1.0 - np.cos(np.pi * (frequency[rising] - f1) / (f2 - f1))
    )
    weights[flat] = 1.0
    weights[falling] = 0.5 * (
        1.0 - np.cos(np.pi * (f4 - frequency[falling]) / (f4 - f3))
    )
    return weights


def _snapshot(root: Path) -> dict[str, tuple[str, bytes | None]]:
    return {
        path.relative_to(root).as_posix(): (
            "directory" if path.is_dir() else "file",
            None if path.is_dir() else path.read_bytes(),
        )
        for path in sorted(root.rglob("*"))
    }


class LegacyProcessedSpectrumTests(unittest.TestCase):
    def test_nonzero_damping_is_undone_once_before_filtering(self):
        nfft = 128
        dt = 0.25
        df = 1.0 / (nfft * dt)
        epsilon = 0.017
        output_samples = nfft
        time = np.arange(nfft, dtype=float) * dt
        frequency = np.arange(nfft // 2 + 1, dtype=float) * df
        physical = (0.7 * np.sin(2 * np.pi * frequency[1] * time + 0.2)
                    + 0.3 * np.cos(2 * np.pi * frequency[2] * time - 0.4)
                    + 0.2 * np.sin(2 * np.pi * frequency[3] * time + 0.9))
        corners = [10.0, 20.0, 80.0, 110.0]
        final_window_fraction = 0.5
        damped_bins = np.fft.rfft(physical * np.exp(-epsilon * time)) / (nfft * df)
        spectra = {(index, 0): value for index, value in enumerate(damped_bins)}
        meta = {
            "samples": nfft,
            "sample_interval_s": dt,
            "frequency_spacing_hz": df,
            "damping_epsilon_rad_s": epsilon,
            "output_samples": output_samples,
        }

        actual_frequency, actual = legacy_processed_amplitude_spectrum(
            spectra, meta, corners, final_window_fraction
        )

        duration = max((output_samples - 1) * dt, dt)
        pre_window = _hann_ends(time, duration, 0.05)
        pre_spectrum = np.fft.rfft(physical * pre_window, n=nfft)
        filtered_time = np.fft.irfft(
            pre_spectrum * _half_cosine(frequency, corners), n=nfft
        ) * nfft * df
        expected = np.abs(np.fft.rfft(
            filtered_time * _hann_ends(time, duration, final_window_fraction),
            n=nfft,
        ))

        actual_shape = actual / np.max(actual)
        expected_shape = expected / np.max(expected)
        self.assertTrue(np.array_equal(actual_frequency, frequency))
        np.testing.assert_allclose(actual_shape, expected_shape, rtol=2e-12,
                                   atol=2e-12)

        twice_damped = filtered_time * np.exp(epsilon * time)
        wrong = np.abs(np.fft.rfft(
            twice_damped * _hann_ends(time, duration, final_window_fraction),
            n=nfft,
        ))
        wrong_shape = wrong / np.max(wrong)
        self.assertGreater(float(np.mean(np.abs(actual_shape - wrong_shape))),
                           1e-3)


class OutputDirectoryGuardTests(unittest.TestCase):
    def test_populated_run_directories_stop_before_any_subprocess_or_write(self):
        scripts = (
            (YSPECPP / "run_paper2.py", ["--dspec-source", "missing-source"]),
            (YSPECPP / "run_benchmark.py", [
                "--dspec-source", "missing-source", "--yspec-source", "missing-yspec",
                "--dependency-root", "missing-dependencies",
            ]),
        )
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            hook_dir = root / "hooks"
            hook_dir.mkdir()
            (hook_dir / "sitecustomize.py").write_text(
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
            env = os.environ.copy()
            env["PYTHONDONTWRITEBYTECODE"] = "1"
            env["PYTHONPATH"] = str(hook_dir)
            for script, extra_args in scripts:
                with self.subTest(script=script.name):
                    output = root / f"{script.stem}-run"
                    (output / "nested").mkdir(parents=True)
                    (output / ".provenance").write_text("kept", encoding="utf-8")
                    (output / "nested" / "input.json").write_text(
                        '{"frozen": true}\n', encoding="utf-8"
                    )
                    before = _snapshot(output)
                    process_sentinel = root / f"{script.stem}-process-called"
                    env["PROCESS_SENTINEL"] = str(process_sentinel)
                    result = subprocess.run(
                        [sys.executable, "-B", str(script), *extra_args,
                         "--output", str(output)],
                        cwd=root,
                        env=env,
                        capture_output=True,
                        text=True,
                        check=False,
                    )
                    self.assertNotEqual(result.returncode, 0, result.stdout)
                    self.assertIn("output directory must be empty", result.stderr)
                    self.assertFalse(process_sentinel.exists())
                    self.assertEqual(_snapshot(output), before)


if __name__ == "__main__":
    unittest.main()
