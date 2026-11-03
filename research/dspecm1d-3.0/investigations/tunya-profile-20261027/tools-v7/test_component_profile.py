#!/usr/bin/env python3
"""Small isolated checks for snapshot provenance and profile accounting."""
import importlib.util
import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("collector", HERE / "collect_component_profile.py")
collector = importlib.util.module_from_spec(spec)
spec.loader.exec_module(collector)
spec2 = importlib.util.spec_from_file_location("finalizer", HERE / "finalize_profile.py")
finalizer = importlib.util.module_from_spec(spec2)
spec2.loader.exec_module(finalizer)


def test_profile_closure():
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "job.log"
        p.write_text('MATHPROFILE {"solver":"synthetic","objective":0,"stages":{"physical_case_total":{"inclusive_s":3.0,"exclusive_s":1.0,"calls":1},"leaf":{"inclusive_s":2.0,"exclusive_s":2.0,"calls":1}},"counters":{"n":2}}\n')
        rec, closure = collector.profile_from_log(p)
        assert rec["solver"] == "synthetic" and closure["closes_within_existing_tolerance"]
        assert closure["exclusive_sum_s"] == 3.0
        p.write_text('MATHPROFILE {"solver":"synthetic","objective":0,"stages":{"physical_case_total":{"inclusive_s":1.0,"exclusive_s":1.0,"calls":-1}}}\n')
        try:
            collector.profile_from_log(p)
        except ValueError:
            pass
        else:
            raise AssertionError("negative profile counter accepted")


def test_finalizer_snapshot_and_hash_rejection():
    with tempfile.TemporaryDirectory() as td:
        base = Path(td); task = base / "task"; task.mkdir(); out = base / "out"; out.mkdir()
        jobdir = task / "baseline-l128" / "data"; jobdir.mkdir(parents=True)
        (task / "baseline-l128" / "job.log").write_text("synthetic log\n")
        (jobdir / "values.csv").write_text("x,y\n1,2\n3,4\n")
        outputs = {"data/values.csv": finalizer.sha(jobdir / "values.csv")}
        logsha = finalizer.sha(task / "baseline-l128" / "job.log")
        state = {"task_id": "profile-06", "thread_cap": 50, "status": "completed", "jobs": [{
            "id": "baseline-l128", "status": "completed", "returncode": 0, "inputs_unchanged": True,
            "cleanup": {"verified_gone": True}, "peak_observed_threads": 1,
            "log": "job.log", "log_sha256": logsha, "outputs": list(outputs), "outputs_sha256": outputs}]}
        (task / "state.json").write_text(json.dumps(state))
        snap = finalizer.snapshot(task, out, state)
        loaded, jobs, listed = collector.validate_snapshot(snap)
        assert loaded["status"] == "completed" and "baseline-l128/data/values.csv" in listed
        summary = collector.collect(snap)
        assert summary["status"] == "INCOMPLETE"
        assert "baseline-l128" not in summary["pending_jobs"]
        (snap / "baseline-l128/data/values.csv").write_text("x,y\n1,999\n")
        try:
            collector.validate_snapshot(snap)
        except ValueError as e:
            assert "hash mismatch" in str(e)
        else:
            raise AssertionError("tampered snapshot was accepted")


def test_exact_plan_state_gate():
    assert len(collector.EXPECTED_IDS) == 26
    assert "dsv4-native-l1536" not in collector.EXPECTED_IDS
    base_job = {"status": "completed", "returncode": 0, "inputs_unchanged": True,
                "cleanup": {"verified_gone": True}, "peak_observed_threads": 2}
    state = {"jobs": [dict(base_job, id=jid) for jid in collector.EXPECTED_IDS]}
    finalizer.require_completed_state(state)
    state["jobs"].pop()
    try:
        finalizer.require_completed_state(state)
    except RuntimeError:
        pass
    else:
        raise AssertionError("incomplete 26-job state accepted")


def test_material_aware_cross_code_mapping():
    import numpy as np
    nodes = [{"element": 7, "node": 4, "radius_m": 4.5e6, "fluid": 0},
             {"element": 12, "node": 3, "radius_m": 2.5e6, "fluid": 1}]
    d = np.zeros((1, 2, 6), dtype=float)
    y = np.zeros((1, 2, 5), dtype=float)
    # In fluid, DS kappa matches YS A+C+F although all DS A..N channels are zero.
    y[0, 1, :3] = (2.0, 3.0, 4.0)
    d[0, 1, 5] = 9.0
    # In solid, paired A channels intentionally differ and must remain visible.
    d[0, 0, 0], y[0, 0, 0] = 2.0, 1.0
    result = collector.material_channel_comparisons(d, y, nodes, 1)
    assert result["all_fluid_radial_sides"]["channels"]["kappa_vs_YSpec_A_plus_C_plus_F"]["max_abs_difference"] == 0
    assert result["regular_outer_core_2200_3200_km"]["channels"]["kappa_vs_YSpec_A_plus_C_plus_F"]["max_abs_difference"] == 0
    assert result["all_solid_radial_sides"]["channels"]["A"]["max_abs_difference"] == 1
    assert result["regular_solid_4200_6200_km"]["channels"]["A"]["max_abs_difference"] == 1


def test_finalizer_cli_incomplete_and_failed_queue():
    with tempfile.TemporaryDirectory() as td:
        base = Path(td); task = base / "task"; task.mkdir(); out = base / "fresh-output"
        (task / "state.json").write_text(json.dumps({"task_id": "profile-06", "thread_cap": 50,
            "status": "completed", "jobs": []}))
        p = subprocess.run([sys.executable, str(HERE / "finalize_profile.py"), str(task), str(out),
                            "--max-wait-seconds", "0"], capture_output=True, text=True, timeout=15)
        assert p.returncode == 1 and "exact 26-job" in p.stderr
        assert (out / "finalization-failure.json").is_file()
        fail_task = base / "failed-task"; fail_task.mkdir(); fail_out = base / "failed-output"
        partial_dir = fail_task / "baseline-l128" / "data"; partial_dir.mkdir(parents=True)
        (fail_task / "baseline-l128" / "job.log").write_text("preserved partial log\n")
        (partial_dir / "values.csv").write_text("x,y\n1,2\n")
        outputs = {"data/values.csv": finalizer.sha(partial_dir / "values.csv")}
        completed_job = {"id": "baseline-l128", "status": "completed", "returncode": 0,
            "inputs_unchanged": True, "cleanup": {"verified_gone": True}, "peak_observed_threads": 2,
            "log": "job.log", "log_sha256": finalizer.sha(fail_task / "baseline-l128" / "job.log"),
            "outputs": list(outputs), "outputs_sha256": outputs}
        (fail_task / "state.json").write_text(json.dumps({"task_id": "profile-06", "thread_cap": 50,
            "status": "failed", "jobs": [completed_job, {"id": "dsv6-profiled-l128", "status": "failed"}]}))
        p = subprocess.run([sys.executable, str(HERE / "finalize_profile.py"), str(fail_task), str(fail_out),
                            "--max-wait-seconds", "0"], capture_output=True, text=True, timeout=15)
        assert p.returncode == 2
        assert json.loads((fail_out / "failure-manifest.json").read_text())["status"] == "campaign_failed"
        assert (fail_out / "snapshot/baseline-l128/data/values.csv").is_file()


if __name__ == "__main__":
    test_profile_closure()
    test_finalizer_snapshot_and_hash_rejection()
    test_exact_plan_state_gate()
    test_material_aware_cross_code_mapping()
    test_finalizer_cli_incomplete_and_failed_queue()
    print("component profile synthetic checks passed")
