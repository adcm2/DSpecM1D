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
        p.write_text('MATHPROFILE {"solver":"synthetic","objective":0,"stages":{"root_total":{"inclusive_s":3.0,"exclusive_s":1.0,"calls":1},"leaf":{"inclusive_s":2.0,"exclusive_s":2.0,"calls":1}},"counters":{"n":2}}\n')
        rec, closure = collector.profile_from_log(p)
        assert rec["solver"] == "synthetic" and closure["closes_within_1e-6_s"]
        assert closure["exclusive_sum_s"] == 3.0


def test_finalizer_snapshot_and_hash_rejection():
    with tempfile.TemporaryDirectory() as td:
        base = Path(td); task = base / "task"; task.mkdir(); out = base / "out"; out.mkdir()
        jobdir = task / "synthetic" / "data"; jobdir.mkdir(parents=True)
        (task / "synthetic" / "job.log").write_text("synthetic log\n")
        (jobdir / "values.csv").write_text("x,y\n1,2\n3,4\n")
        outputs = {"data/values.csv": finalizer.sha(jobdir / "values.csv")}
        logsha = finalizer.sha(task / "synthetic" / "job.log")
        state = {"task_id": "profile-06", "thread_cap": 50, "status": "completed", "jobs": [{
            "id": "synthetic", "status": "completed", "returncode": 0, "inputs_unchanged": True,
            "cleanup": {"verified_gone": True}, "peak_observed_threads": 1,
            "log": "job.log", "log_sha256": logsha, "outputs": list(outputs), "outputs_sha256": outputs}]}
        (task / "state.json").write_text(json.dumps(state))
        snap = finalizer.snapshot(task, out, state)
        loaded, jobs, listed = collector.validate_snapshot(snap)
        assert loaded["status"] == "completed" and "synthetic/data/values.csv" in listed
        summary = collector.collect(snap)
        assert "synthetic" not in summary["pending_jobs"]
        (snap / "synthetic/data/values.csv").write_text("x,y\n1,999\n")
        try:
            collector.validate_snapshot(snap)
        except ValueError as e:
            assert "hash mismatch" in str(e)
        else:
            raise AssertionError("tampered snapshot was accepted")


def test_finalizer_cli_completed_queue():
    with tempfile.TemporaryDirectory() as td:
        base = Path(td); task = base / "task"; task.mkdir(); out = base / "fresh-output"
        (task / "state.json").write_text(json.dumps({"task_id": "profile-06", "thread_cap": 50,
            "status": "completed", "jobs": []}))
        p = subprocess.run([sys.executable, str(HERE / "finalize_profile.py"), str(task), str(out),
                            "--max-wait-seconds", "0"], capture_output=True, text=True, timeout=15)
        assert p.returncode == 0, p.stderr
        receipt = json.loads((out / "finalization-receipt.json").read_text())
        assert receipt["status"] == "completed"
        assert (out / "component-profile-summary.json").is_file()


if __name__ == "__main__":
    test_profile_closure()
    test_finalizer_snapshot_and_hash_rejection()
    test_finalizer_cli_completed_queue()
    print("component profile synthetic checks passed")
