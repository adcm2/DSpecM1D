import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import hashlib

MODULE_PATH = Path(__file__).with_name("check_stages12_stage1.py")
SPEC = importlib.util.spec_from_file_location("check_stages12_stage1", MODULE_PATH)
CHECK = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECK)


def valid_profile():
    rows = []
    for phase, exclusive in [
        ("setup", 2.0),
        ("initial_frequency_workspace_and_capture", 1.0),
        ("measurement_and_initial_output", 1.0),
        ("kernel_production_excluding_isolated_diagnostics", 2.0),
        ("kernel_output_and_result_io", 1.0),
        ("diagnostic_tensor_validation", 1.0),
    ]:
        rows.append({"phase": phase, "scope_name": phase, "execution_domain": "master_wall",
                     "inclusive_s": exclusive, "exclusive_s": exclusive, "calls": 1,
                     "logical_component": "test", "diagnostic": phase == "diagnostic_tensor_validation",
                     "parent_overlap_description": "exclusive time is partitioned"})
    rows.append({"phase": "profile_pipeline_total", "scope_name": "profile_pipeline_wall",
                 "execution_domain": "master_wall", "inclusive_s": 10.0, "exclusive_s": 2.0,
                 "calls": 1, "logical_component": "pipeline", "diagnostic": False,
                 "parent_overlap_description": "exclusive pipeline remainder"})
    rows.append({"phase": "kernel_production_excluding_isolated_diagnostics", "scope_name": "parallel_worker.sparse_solve",
                 "execution_domain": "worker_elapsed_sum", "inclusive_s": 12.0, "exclusive_s": 11.0,
                 "calls": 2, "logical_component": "modal_solver", "diagnostic": False,
                 "parent_overlap_description": "worker elapsed sum overlaps degree-region wall"})
    return {"schema": 1, "snapshot_report_overhead_s": 0.01,
            "wall_partition": {"profile_pipeline_wall_s": 10.0,
                                "master_exclusive_partition_s": 10.0,
                                "residual_s": 0.0, "tolerance_s": 0.01,
                                "worker_elapsed_excluded": True},
            "rows": rows, "counters": []}


def make_receipted_job(root, job_id, workers, across=b"same-across"):
    job_dir = root / job_id
    data = job_dir / "data"
    data.mkdir(parents=True)
    result = {"backend": "dspec", "degree": 3, "nq": 8, "refinement": 0.05,
              "workers": workers, "attenuation": True, "epsilon_perturbation": 0.0,
              "frequency_count": 3, "first_bin": 7, "last_bin": 9,
              "fft_samples": 2048, "fft_period_s": 4096.0,
              "native_side_count": 2, "along_angles": 3, "across_angles": 2,
              "gravity_max_abs_nd_error": 0.0, "time_shift_adjoint_identity": 1.0,
              "measurement_rate_norm": 1.0, "integrated_profile_kernel": 1.0,
              "cache_bytes": 100, "forward_s": 1.0, "measurement_s": 0.1,
              "kernel_s": 2.0, "tensor_validation_s": 0.0, "elapsed_s": 4.0,
              "kernel_parameter": "ln_vs", "kernel_units": "s m^-3",
              "layout": "native-side-major", "source_shell_contacts": False}
    files = {"parameters.txt": b"parameters\n", "result.json": json.dumps(result).encode(),
             "trace.csv": b"trace\n", "spectrum.csv": b"spectrum\n",
             "components.csv": b"components\n", "nodes.csv": b"nodes\n",
             "along.bin": b"along", "along-angles.csv": b"along angles\n",
             "across.bin": across, "across-angles.csv": b"across angles\n"}
    for name, contents in files.items():
        (data / name).write_bytes(contents)
    (job_dir / "job.log").write_text("complete\n")
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    outputs = ["data/" + name for name in files]
    state = {"jobs": [{"id": job_id, "status": "completed", "returncode": 0,
                       "inputs_unchanged": True, "cleanup": {"verified_gone": True},
                       "failure_reasons": [], "outputs": outputs,
                       "outputs_sha256": {rel: sha(job_dir / rel) for rel in outputs},
                       "log": "job.log", "log_sha256": sha(job_dir / "job.log")} ]}
    (root / "state.json").write_text(json.dumps(state))


class Stage1CheckerTests(unittest.TestCase):
    def test_valid_partition_ignores_inclusive_and_worker_elapsed(self):
        result = CHECK.validate_profile(valid_profile())
        self.assertEqual(result["wall_residual_s"], 0)

    def test_rejects_inclusive_sum_double_counting(self):
        profile = valid_profile()
        profile["wall_partition"]["master_exclusive_partition_s"] = 30.0
        with self.assertRaisesRegex(ValueError, "inclusive or worker"):
            CHECK.validate_profile(profile)

    def test_rejects_worker_elapsed_mixed_into_wall(self):
        profile = valid_profile()
        profile["wall_partition"]["master_exclusive_partition_s"] = 21.0
        with self.assertRaisesRegex(ValueError, "inclusive or worker"):
            CHECK.validate_profile(profile)

    def test_rejects_worker_name_tagged_as_master_wall(self):
        profile = valid_profile()
        worker = next(row for row in profile["rows"] if row["scope_name"].startswith("parallel_worker."))
        worker["execution_domain"] = "master_wall"
        with self.assertRaisesRegex(ValueError, "worker scope name and execution domain"):
            CHECK.validate_profile(profile)

    def test_rejects_diagnostic_contamination(self):
        profile = valid_profile()
        row = next(row for row in profile["rows"] if row["phase"] == "kernel_production_excluding_isolated_diagnostics")
        row["diagnostic"] = True
        with self.assertRaisesRegex(ValueError, "diagnostic phase tagging"):
            CHECK.validate_profile(profile)

    def test_rejects_missing_diagnostic_isolation(self):
        profile = valid_profile()
        profile["rows"] = [row for row in profile["rows"] if row["phase"] != "diagnostic_tensor_validation"]
        with self.assertRaisesRegex(ValueError, "missing required profile phases"):
            CHECK.validate_profile(profile)

    def test_rejects_modified_untouched_baseline_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            candidate = Path(directory) / "candidate"
            candidate.mkdir()
            (candidate / "doc_driver.cpp").write_text("allowed edit")
            (candidate / "fast_engine.hpp").write_text("modified dependency")
            import hashlib
            manifest = {"files": {
                "doc_driver.cpp": hashlib.sha256(b"original driver").hexdigest(),
                "fast_engine.hpp": hashlib.sha256(b"original engine").hexdigest(),
            }}
            manifest_path = Path(directory) / "baseline.json"
            manifest_path.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, "outside Stage 1 scope"):
                CHECK.validate_baseline(candidate, manifest_path)

    def test_allows_only_declared_candidate_edits(self):
        with tempfile.TemporaryDirectory() as directory:
            candidate = Path(directory) / "candidate"
            candidate.mkdir()
            original = b"old"
            (candidate / "doc_driver.cpp").write_bytes(b"new")
            (candidate / "fast_engine.hpp").write_bytes(original)
            import hashlib
            manifest_path = Path(directory) / "baseline.json"
            manifest_path.write_text(json.dumps({"files": {
                "doc_driver.cpp": hashlib.sha256(original).hexdigest(),
                "fast_engine.hpp": hashlib.sha256(original).hexdigest(),
            }}))
            result = CHECK.validate_baseline(candidate, manifest_path)
            self.assertEqual(result["changed_files"], ["doc_driver.cpp"])

    def test_parity_checks_prefixed_kernel_files_and_expected_worker_counts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            a, b = root / "a", root / "b"
            a.mkdir(); b.mkdir()
            make_receipted_job(a, "one", 1)
            make_receipted_job(b, "eight", 8)
            result = CHECK.validate_parity(a, "one", b, "eight", 1, 8)
            self.assertIn("across.bin", result["byte_exact_scientific_files"])

    def test_parity_rejects_across_kernel_difference(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            a, b = root / "a", root / "b"
            a.mkdir(); b.mkdir()
            make_receipted_job(a, "one", 1)
            make_receipted_job(b, "eight", 8, across=b"changed")
            with self.assertRaisesRegex(ValueError, "scientific output differs: across.bin"):
                CHECK.validate_parity(a, "one", b, "eight", 1, 8)

    def test_parity_rejects_unplanned_worker_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            a, b = root / "a", root / "b"
            a.mkdir(); b.mkdir()
            make_receipted_job(a, "one", 1)
            make_receipted_job(b, "eight", 8)
            with self.assertRaisesRegex(ValueError, "worker count does not match"):
                CHECK.validate_parity(a, "one", b, "eight", 8, 1)


if __name__ == "__main__":
    unittest.main()
