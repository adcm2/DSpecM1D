"""Lightweight guards for the 3D kernel timing report tools."""
import csv
import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace
from unittest import TestCase
from unittest.mock import patch

RUNNER = Path(__file__).parents[1] / "benchmarks/yspecpp/kernels/3d/benchmark_kernel_runtime.py"
SPEC = importlib.util.spec_from_file_location("kernel_runtime", RUNNER)
RUNTIME = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNTIME)


class KernelTimingToolTests(TestCase):
    def write_kernel_csv(self, path, changed=None, nonfinite=False):
        fields = ["objective", "modulus", "radius_m", "theta_rad", "phi_rad",
                  "kernel_per_Pa_m3"]
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for objective in ("J_Z", "J_N"):
                for modulus in RUNTIME.MODULI:
                    for index, radius in enumerate((5.0, 5.0)):
                        value = 1.0
                        if changed == (objective, modulus, index):
                            value = 1.1
                        if nonfinite and objective == "J_Z" and modulus == "A" and index == 0:
                            value = "nan"
                        writer.writerow({"objective": objective, "modulus": modulus,
                                         "radius_m": radius, "theta_rad": 0.1,
                                         "phi_rad": 0.2 + index * 0.1,
                                         "kernel_per_Pa_m3": value})

    def write_samples_csv(self, path):
        fields = ["objective", "modulus", "radius_m", "theta_rad", "phi_rad",
                  "element", "fluid", "kernel_objective_per_Pa_m3", "deltaX_Pa",
                  "dr_m", "dOmega"]
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for modulus in RUNTIME.MODULI:
                for index, radius in enumerate((5.0, 5.0)):
                    writer.writerow({"objective": "J_Z", "modulus": modulus,
                                     "radius_m": radius, "theta_rad": 0.1,
                                     "phi_rad": 0.2 + index * 0.1, "element": index,
                                     "fluid": 0, "kernel_objective_per_Pa_m3": 1.0,
                                     "deltaX_Pa": 2.0 * (1 + RUNTIME.MODULI.index(modulus)),
                                     "dr_m": 3.0, "dOmega": 4.0})

    def test_nonempty_output_is_rejected_before_external_work(self):
        from tempfile import TemporaryDirectory
        with TemporaryDirectory() as temp:
            output = Path(temp) / "existing"
            output.mkdir()
            marker = output / "keep.txt"
            marker.write_text("preserve", encoding="utf-8")
            with patch.object(sys, "argv", [str(RUNNER), "--output", str(output)]), \
                    patch.object(RUNTIME.subprocess, "run", side_effect=AssertionError("executed")) as run:
                with self.assertRaises(SystemExit):
                    RUNTIME.main()
            run.assert_not_called()
            self.assertEqual(marker.read_text(encoding="utf-8"), "preserve")
            invalid_output = Path(temp) / "invalid"
            with patch.object(sys, "argv", [str(RUNNER), "--output", str(invalid_output),
                                            "--accuracy-tolerance", "nan"]):
                with self.assertRaisesRegex(SystemExit, "accuracy tolerance"):
                    RUNTIME.main()
            self.assertFalse(invalid_output.exists())

    def test_timed_gate_and_main_first_pair_path(self):
        from tempfile import TemporaryDirectory
        with TemporaryDirectory() as temp:
            root = Path(temp)
            dspec, yspec, samples = root / "d.csv", root / "y.csv", root / "samples.csv"
            self.write_kernel_csv(dspec)
            self.write_kernel_csv(yspec)
            self.write_samples_csv(samples)
            result = RUNTIME.accuracy_gate(
                dspec, yspec, samples, 1e-3, 2,
                {"objective_scale_S": {"J_Z": 1.0, "J_N": 1.0}})
            by_modulus = {row["modulus"]: row["Kscale_per_Pa_m3"]
                          for row in result["profiles"] if row["objective"] == "J_Z"}
            self.assertEqual(by_modulus["A"], by_modulus["C"])
            self.assertEqual(by_modulus["L"], by_modulus["N"])
            self.assertNotEqual(by_modulus["A"], by_modulus["F"])
            self.assertNotEqual(by_modulus["F"], by_modulus["L"])
            self.write_kernel_csv(yspec, changed=("J_Z", "A", 0))
            scale = {"J_Z": 1.0, "J_N": 1.0}
            with self.assertRaisesRegex(ValueError, "actual timed-output accuracy gate failed"):
                RUNTIME.accuracy_gate(dspec, yspec, samples, 1e-3, 2,
                                      {"objective_scale_S": scale})

            model = root / "model.deck"
            model.write_text("homogeneous fixture\n", encoding="utf-8")
            metadata = {"solver": "DSpecM1D", "objective": "J_Z=Re(sum(conj(q)*P u))",
                        "lmax": 3, "frequency_bins_zero_based": [1, 2, 3], "nt": 128,
                        "angular_nodes": 192, "sample_count": 2,
                        "dt_s": 20.0, "df_hz": 1.0 / 2560.0, "epsilon_s_inv": 1.0 / 360.0,
                        "physical_G_SI": 6.672299999999999e-11, "gravity": "full",
                        "attenuation": False,
                        "source": {"latitude_deg": 10.0, "longitude_deg": 20.0,
                                   "depth_km": 50.0, "radius_m": 6321000.0},
                        "receiver": {"latitude_deg": 45.0, "longitude_deg": 90.0, "depth_km": 3.0},
                        "moment_tensor_Nm": [1e19, 2e18, -3e18, -5e18, 4e18, -5e18],
                        "q_nonzero": {"Z,1": [0.7, 0.2], "Z,2": [-0.25, 0.45]},
                        "time_window": {"type": "Hann", "sample_start_inclusive": 0,
                                        "sample_end_inclusive": 90, "sample_count": 128},
                        "provenance": {"model_sha256": RUNTIME.sha256(model)},
                        "render_provenance": {"numeric_inputs_sha256": {
                            "samples.csv": RUNTIME.sha256(samples)}}}
            (root / "run.json").write_text(json.dumps(metadata), encoding="utf-8")
            paths = [root / name for name in ("dspec", "yspec", "dspec-src", "yspec-src",
                                               "dspec-build", "yspec-build")]
            for path in paths:
                path.mkdir()
            executables = [root / "dspec-exe", root / "yspec-exe"]
            for path in executables:
                path.write_text("stub", encoding="utf-8")
            output = root / "timing-output"
            config = {"matched": True}

            def fake_run_one(executable, solver, args, run_dir, env, started,
                             model_hash, samples_hash):
                run_dir.mkdir()
                self.write_kernel_csv(run_dir / "kernels.csv")
                api = [0.2, 0.3] if solver == "DSpecM1D" else [0.4]
                return {"solver": solver, "api_s": api, "setup_s": 0.1,
                        "warmup_batch_s": 0.05, "first_workload_api_s": 0.05,
                        "first_model_to_arrays_s": 0.15,
                        "calculation_s": 0.6, "end_to_end_s": 0.8,
                        "peak_rss_bytes": 1024, "peak_rss_label": "process high-water RSS",
                        "sample_count": 2, "config": config}

            argv = [str(RUNNER), "--dspec", str(executables[0]), "--yspec", str(executables[1]),
                    "--dspec-source", str(paths[0]), "--yspec-source", str(paths[1]),
                    "--dspec-build", str(paths[4]), "--yspec-build", str(paths[5]),
                    "--model", str(model), "--samples", str(samples), "--output", str(output),
                    "--measurements", "2", "--repeats", "1"]
            build_record = {"compiler": "/usr/bin/c++", "compiler_version": "test compiler",
                            "selected_dependency_pins": {"PLANETMODEL_TAG": "pin"}}
            with patch.object(sys, "argv", argv), \
                    patch.object(RUNTIME, "build_info", return_value=build_record), \
                    patch.object(RUNTIME, "run_one", side_effect=fake_run_one), \
                    patch.object(RUNTIME, "assert_config_matches_tutorial"), \
                    patch.object(RUNTIME, "assert_configs_match"), \
                    patch.object(RUNTIME, "accuracy_gate", return_value={
                        "passed": True, "profiles": [], "integrals": []}) as gate:
                RUNTIME.main()
            self.assertEqual(gate.call_args.args[4], 2)
            self.assertEqual(gate.call_args.args[5]["solver"], "DSpecM1D")
            self.assertEqual(gate.call_args.args[6], metadata)
            self.assertTrue((output / "report.html").is_file())
            saved_rows = RUNTIME.read_csv(output / "timings.csv")
            dspec_row = next(row for row in saved_rows if row["solver"] == "DSpecM1D")
            self.assertEqual(json.loads(dspec_row["api_s_by_objective"]), [0.2, 0.3])
            self.assertEqual(float(dspec_row["first_workload_api_s"]), 0.05)

    def test_nonfinite_timed_array_is_rejected(self):
        from tempfile import TemporaryDirectory
        with TemporaryDirectory() as temp:
            root = Path(temp)
            dspec, yspec, samples = root / "d.csv", root / "y.csv", root / "samples.csv"
            self.write_kernel_csv(dspec)
            self.write_kernel_csv(yspec, nonfinite=True)
            self.write_samples_csv(samples)
            with self.assertRaisesRegex(ValueError, "nonfinite"):
                RUNTIME.accuracy_gate(dspec, yspec, samples, 1e-3, 2,
                                      {"objective_scale_S": {"J_Z": 1.0, "J_N": 1.0}})
            text = samples.read_text(encoding="utf-8").replace(",3.0,4.0", ",nan,4.0", 1)
            samples.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "nonfinite"):
                RUNTIME.validate_samples_csv(samples, 2)

    def test_frozen_tutorial_configuration_mismatch_is_rejected(self):
        from tempfile import TemporaryDirectory
        with TemporaryDirectory() as temp:
            path = Path(temp) / "run.json"
            metadata = {"solver": "DSpecM1D", "objective": "J_Z=Re(sum(conj(q)*P u))",
                        "lmax": 3, "frequency_bins_zero_based": [1, 2, 3],
                        "angular_nodes": 192, "sample_count": 384,
                        "nt": 128, "dt_s": 20.0, "df_hz": 1.0 / 2560.0,
                        "epsilon_s_inv": 1.0 / 360.0,
                        "physical_G_SI": 6.672299999999999e-11,
                        "gravity": "full", "attenuation": False,
                        "source": {"latitude_deg": 10.0, "longitude_deg": 20.0, "depth_km": 50.0},
                        "receiver": {"latitude_deg": 45.0, "longitude_deg": 90.0, "depth_km": 3.0},
                        "moment_tensor_Nm": [1e19, 2e18, -3e18, -5e18, 4e18, -5e18],
                        "q_nonzero": {"Z,1": [0.7, 0.2], "Z,2": [-0.25, 0.45]},
                        "time_window": {"type": "Hann", "sample_start_inclusive": 0,
                                        "sample_end_inclusive": 90, "sample_count": 128},
                        "provenance": {"model_sha256": "model"},
                        "render_provenance": {"numeric_inputs_sha256": {"samples.csv": "samples"}}}
            path.write_text(json.dumps(metadata), encoding="utf-8")
            RUNTIME.validate_tutorial_metadata(path, 3, "model", "samples")
            config = {"frequency_grid": {"nt": 128, "dt_s": 20.0,
                                         "df_hz": 1.0 / 2560.0, "epsilon_s_inv": 1.0 / 360.0,
                                         "solver_bins_zero_based": [1, 2, 3],
                                         "degrees_inclusive": [0, 3]},
                      "physics": {"gravity": "full", "attenuation": False,
                                  "physical_G": metadata["physical_G_SI"]},
                      "source": {"geocentric_latitude_deg": 10.0, "longitude_deg": 20.0,
                                 "depth_m": 50000.0,
                                 "moment_tensor_Nm": metadata["moment_tensor_Nm"]},
                      "receiver": {"geocentric_latitude_deg": 45.0, "longitude_deg": 90.0,
                                   "depth_m": 3000.0},
                      "window": {"type": "Hann", "first_sample": 0,
                                 "last_sample": 90, "nt": 128},
                      "objective": {"names": [{"name": "J_Z", "component": "Z",
                                                  "q_nonzero": {"bin_1": [0.7, 0.2],
                                                                "bin_2": [-0.25, 0.45]}}]},
                      "angular_node_count": 192, "radial_node_count": 2}
            RUNTIME.assert_config_matches_tutorial(config, metadata, 3, 1, 384)
            config["physics"]["physical_G"] = float("inf")
            with self.assertRaisesRegex(ValueError, "nonfinite"):
                RUNTIME.assert_config_matches_tutorial(config, metadata, 3, 1, 384)
            config["physics"]["physical_G"] = metadata["physical_G_SI"]
            config["source"]["depth_m"] = 40000.0
            with self.assertRaisesRegex(ValueError, "source differs"):
                RUNTIME.assert_config_matches_tutorial(config, metadata, 3, 1, 384)
            metadata["dt_s"] = 40.0
            path.write_text(json.dumps(metadata), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "metadata mismatch for dt_s"):
                RUNTIME.validate_tutorial_metadata(path, 3, "model", "samples")

    def test_replay_preserves_saved_timing_values(self):
        from tempfile import TemporaryDirectory
        with TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "saved"
            source.mkdir()
            saved = source / "timings.csv"
            saved.write_text(
                "replicate,solver,api_s,api_s_by_objective,first_workload_api_s,"
                "first_model_to_arrays_s,calculation_s,end_to_end_s,setup_s,warmup_batch_s,"
                "peak_rss_bytes,peak_rss_label,measurements,point_count\n"
                '1,DSpecM1D,2.5,"[1.1, 1.4]",2.4,2.6,2.7,3.0,0.2,2.4,1024,'
                '"process high-water RSS",2,192\n'
                '1,YSpecpp,1.5,"[1.5]",1.5,1.6,1.7,2.0,0.1,1.5,2048,'
                '"process high-water RSS",2,192\n',
                encoding="utf-8")
            original = saved.read_bytes()
            output = root / "replayed"
            RUNTIME.replay(SimpleNamespace(replay=source, output=output))
            self.assertEqual(saved.read_bytes(), original)
            self.assertTrue((output / "kernel_runtime.png").is_file())
            replayed = RUNTIME.read_csv(output / "timings.csv")
            self.assertEqual([row["replicate"] for row in replayed], ["1", "1"])
            self.assertEqual(json.loads(replayed[0]["api_s_by_objective"]), [1.1, 1.4])
            self.assertEqual(float(replayed[0]["first_workload_api_s"]), 2.4)
            summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["DSpecM1D"]["api_s"], [2.5])
            self.assertTrue((output / "report.html").is_file())
