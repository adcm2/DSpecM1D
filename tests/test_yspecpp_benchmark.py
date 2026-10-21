import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[1]
RUNNER = REPO / "benchmarks" / "yspecpp" / "run_benchmark.py"
sys.path.insert(0, str(RUNNER.parent))
import run_benchmark as benchmark


class OutputDirectoryGuardTests(unittest.TestCase):
    def test_report_reflects_resolved_metrics_and_escapes_provenance(self):
        summary = {
            "provenance": {
                "refinement": "fine", "jobs": 1,
                "runtime_environment": {"OMP_NUM_THREADS": "2"},
                "yspecpp_commit": "<commit>", "yspecpp_filetree_sha256": "tree",
                "yspecpp_source_archive_sha256": "archive", "dspec_source_head": "head",
                "dspec_tracked_source_diff_sha256": "diff", "case_spec_sha256": "cases",
                "dependency_shas": {"dep": "sha"}, "planetmodel_sha": "pin",
            },
            "defaults": {
                "degrees": [0, 1, 2, 3, 4], "frequencies_mHz": [0.5, 1.5],
                "duration_min": 30, "time_step_s": 60,
                "source": {"depth_km": 100, "latitude": 10, "longitude": 20,
                           "moment_tensor_Nm": {"rr": 1e15}},
                "solver": {"dspec": {"nq": 7, "maxstep": 0.00625,
                                      "relative_error": 0.001, "nskip": 2},
                           "yspec": {"tolerance": 1e-9, "start_decay": 11}},
            },
            "cases": [{"id": "solid-full-off", "model": "homogeneous.deck",
                       "receiver_depth_km": 10, "receiver": [45, 90],
                       "ocean": False, "gravity": "full", "attenuation": False}],
            "run_metadata": {"solid-full-off": {"dspec": {
                "samples": 64, "output_samples": 31,
                "frequency_spacing_hz": 0.00026,
                "damping_epsilon_rad_s": 0.0027,
            }}},
            "spectra": [{"case": "solid-full-off", "component": 0,
                         "relative_l2": 0.000002, "difference_over_peak": 0.000003,
                         "yspec_peak": 1.0, "max_abs_difference": 0.000003}],
            "degree_spectra": [
                {"case": "solid-full-off", "degree": 0, "component": 1,
                 "relative_l2": None, "difference_over_peak": None,
                 "yspec_peak": 0.0, "max_abs_difference": 0.0},
                {"case": "solid-full-off", "degree": 0, "component": 2,
                 "relative_l2": None, "difference_over_peak": None,
                 "yspec_peak": 0.0, "max_abs_difference": 0.0},
            ],
            "time_series": [{"case": "solid-full-off", "component": "Z",
                             "relative_l2": 0.000004,
                             "max_difference_over_reference_peak": 0.000005,
                             "reference_peak": 1.0, "max_abs_difference": 0.000005}],
            "degree_time_series": [
                {"case": "solid-full-off", "degree": 0, "component": "N",
                 "relative_l2": None, "max_difference_over_reference_peak": None,
                 "reference_peak": 0.0, "max_abs_difference": 0.0},
                {"case": "solid-full-off", "degree": 0, "component": "E",
                 "relative_l2": None, "max_difference_over_reference_peak": None,
                 "reference_peak": 0.0, "max_abs_difference": 7.4e-310},
            ],
        }
        report = benchmark.render_report(summary, "python3 runner.py --refinement fine")
        self.assertIn("5e−6", report)
        self.assertIn("2e-06", report)
        self.assertIn("5e-06", report)
        self.assertIn("64/31", report)
        self.assertIn("&lt;commit&gt;", report)
        self.assertNotIn("<commit>", report)
        self.assertIn("attenuation off", report)
        self.assertIn("saved spectrum", report)
        self.assertIn("solid-full-off-traces.png", report)
        self.assertIn("Spectrum max absolute (m/s/Hz)", report)
        self.assertIn("Time max absolute (m/s)", report)
        self.assertIn("undefined (zero reference)", report)
        self.assertIn("7.4e-310", report)
        self.assertIn("PASS (limit 5e-06)", report)
        self.assertIn("FAIL (limit 5e-06)", report)
        self.assertIn("subnormal nonzero difference fails", report)

    def test_populated_run_directory_stops_before_subprocesses_or_writes(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            hook = root / "hooks"
            hook.mkdir()
            (hook / "sitecustomize.py").write_text(
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
            output = root / "saved-run"
            (output / "nested").mkdir(parents=True)
            (output / ".provenance").write_text("keep", encoding="utf-8")
            (output / "nested" / "input.json").write_text(
                '{"frozen": true}\n', encoding="utf-8"
            )
            before = {
                path.relative_to(output).as_posix():
                ("directory", None) if path.is_dir() else
                ("file", path.read_bytes())
                for path in sorted(output.rglob("*"))
            }
            sentinel = root / "subprocess-called"
            env = os.environ.copy()
            env.update({
                "PYTHONDONTWRITEBYTECODE": "1",
                "PYTHONPATH": str(hook),
                "PROCESS_SENTINEL": str(sentinel),
            })
            result = subprocess.run(
                [sys.executable, "-B", str(RUNNER),
                 "--dspec-source", "missing-dspec",
                 "--yspec-source", "missing-yspec",
                 "--dependency-root", "missing-dependencies",
                 "--output", str(output)],
                cwd=root, env=env, capture_output=True, text=True, check=False,
            )
            after = {
                path.relative_to(output).as_posix():
                ("directory", None) if path.is_dir() else
                ("file", path.read_bytes())
                for path in sorted(output.rglob("*"))
            }
            self.assertNotEqual(result.returncode, 0, result.stdout)
            self.assertIn("output directory must be empty", result.stderr)
            self.assertFalse(sentinel.exists())
            self.assertEqual(after, before)


if __name__ == "__main__":
    unittest.main()
