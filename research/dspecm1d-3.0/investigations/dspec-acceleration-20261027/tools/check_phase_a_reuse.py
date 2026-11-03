#!/usr/bin/env python3
"""Requalify the completed A4 numerical job from the preserved v3 task receipt."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import check_pair
import check_phase_a


EXPECTED_JOB_ORDER = [
    "phase-a-source-scope-gate",
    "phase-a-build-receipt-gate",
    "phase-a-tiny-contract-test",
    "A1-native-accepted-adjoint-radial-packets",
    "A1-packet-and-frozen-weight-gate",
    "A2-native-same-packet-reconstruction-replay",
    "A2-field-kernel-and-timing-gate",
    "A3-native-L3-central-fd-kernel",
    "A3-frozen-fixture-fd-and-case-gate",
    "A4-profiled-true-forward-only-35bin-control",
    "A4-forward-route-counter-and-parity-gate",
]
TASK_STATE_SHA256 = "e44abdaddc02aa9675407fa85a8a22b9cc8b47499fe94fe24d65089eca9f759f"
A4_JOB_ID = "A4-profiled-true-forward-only-35bin-control"
A4_GATE_ID = "A4-forward-route-counter-and-parity-gate"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def qualify(task_root: Path, expected_state_sha256: str) -> dict:
    state_path = task_root / "state.json"
    state_bytes = state_path.read_bytes()
    state_sha = hashlib.sha256(state_bytes).hexdigest()
    if expected_state_sha256 != TASK_STATE_SHA256 or state_sha != expected_state_sha256:
        raise ValueError("A4 reuse state hash differs from the reviewed failed v3 receipt")
    state = json.loads(state_bytes)
    if state.get("status") != "failed":
        raise ValueError("A4 reuse source task is not the preserved failed v3 task")
    jobs = state.get("jobs", [])
    if [job.get("id") for job in jobs] != EXPECTED_JOB_ORDER:
        raise ValueError("A4 reuse source task job inventory/order changed or A5 already started")
    if any(job.get("status") != "completed" or job.get("returncode") != 0 or
           job.get("inputs_unchanged") is not True or
           job.get("cleanup", {}).get("verified_gone") is not True or
           job.get("failure_reasons") for job in jobs[:-1]):
        raise ValueError("A4 reuse source contains an unqualified prior task")
    failed_gate = jobs[-1]
    if failed_gate.get("status") != "failed" or failed_gate.get("returncode") != 1 or \
       failed_gate.get("failure_reasons") != ["missing_outputs:checks.json"] or \
       failed_gate.get("inputs_unchanged") is not True or \
       failed_gate.get("cleanup", {}).get("verified_gone") is not True:
        raise ValueError("A4 reuse source does not contain the reviewed checker-only failure")
    gate_log = task_root / A4_GATE_ID / failed_gate["log"]
    if sha(gate_log) != failed_gate.get("log_sha256"):
        raise ValueError("A4 failed checker log hash mismatch")
    if "A4 profiler call/RHS counters failed" not in gate_log.read_text():
        raise ValueError("A4 failed checker log does not show the reviewed omitted-zero-counter issue")

    a4_job, artifact_hashes, verified_state_sha = check_pair.verified_job(task_root, A4_JOB_ID)
    if verified_state_sha != expected_state_sha256:
        raise ValueError("A4 completed-job receipt is not bound to the expected task state")
    data_dir = task_root / A4_JOB_ID / "data"
    result = check_phase_a.data_gate("a4", data_dir)
    return {
        "status": "pass",
        "reuse_source_task_state_sha256": verified_state_sha,
        "reuse_source_job_id": A4_JOB_ID,
        "reuse_source_job_returncode": a4_job["returncode"],
        "reuse_source_artifact_sha256": artifact_hashes,
        "rechecked_a4_gate": result,
        "known_v3_checker_failure_job_id": A4_GATE_ID,
        "known_v3_checker_failure_log_sha256": failed_gate["log_sha256"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-root", type=Path, required=True)
    parser.add_argument("--state-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = qualify(args.task_root, args.state_sha256)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
