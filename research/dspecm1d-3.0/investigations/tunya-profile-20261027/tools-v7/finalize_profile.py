#!/usr/bin/env python3
"""Boundedly snapshot completed profile-06 evidence and run its data collector."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from collect_component_profile import EXPECTED_IDS


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_state(path):
    try:
        state = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as e:
        raise RuntimeError(f"campaign state unavailable or malformed: {e}")
    if state.get("task_id") != "profile-06" or state.get("thread_cap") != 50:
        raise RuntimeError("unexpected campaign identity or thread cap")
    return state


def require_completed_state(state):
    jobs_list = state.get("jobs", [])
    jobs = {j.get("id"): j for j in jobs_list}
    if len(jobs) != len(jobs_list) or set(jobs) != set(EXPECTED_IDS):
        raise RuntimeError("completed state does not contain the exact 26-job reviewed plan")
    for jid in EXPECTED_IDS:
        job = jobs[jid]
        if job.get("status") != "completed" or job.get("returncode") != 0:
            raise RuntimeError(f"completed campaign contains a noncompleted/failed job: {jid}")
        if job.get("inputs_unchanged") is not True or job.get("cleanup", {}).get("verified_gone") is not True:
            raise RuntimeError(f"completed job lacks input/cleanup verification: {jid}")
        peak = job.get("peak_observed_threads")
        if peak is None or peak > 50:
            raise RuntimeError(f"missing or exceeded observed thread cap for {jid}: {peak}")


def snapshot(task: Path, out: Path, state):
    snap = out / "snapshot"
    snap.mkdir()
    state_hash_before = sha(task / "state.json")
    shutil.copyfile(task / "state.json", snap / "state.json")
    if sha(snap / "state.json") != state_hash_before or sha(task / "state.json") != state_hash_before:
        raise RuntimeError("campaign state changed during snapshot")
    manifest = {"state.json": state_hash_before}
    for job in state.get("jobs", []):
        if job.get("status") != "completed":
            continue
        jid = job["id"]
        if not jid or Path(jid).name != jid or jid in (".", ".."):
            raise RuntimeError(f"unsafe job id in state: {jid!r}")
        if job.get("returncode") != 0 or job.get("inputs_unchanged") is not True:
            raise RuntimeError(f"completed job lacks successful/invariant receipt: {jid}")
        if job.get("cleanup", {}).get("verified_gone") is not True:
            raise RuntimeError(f"cleanup not verified for {jid}")
        peak = job.get("peak_observed_threads")
        if peak is None or peak > 50:
            raise RuntimeError(f"missing or exceeded observed thread cap for {jid}: {peak}")
        entries = {job.get("log", "job.log"): job.get("log_sha256")}
        entries.update(job.get("outputs_sha256", {}))
        if set(job.get("outputs", [])) != set(job.get("outputs_sha256", {})):
            raise RuntimeError(f"incomplete output hash manifest for {jid}")
        for rel, digest in entries.items():
            if not digest:
                raise RuntimeError(f"missing recorded hash for {jid}/{rel}")
            relpath = Path(rel)
            if relpath.is_absolute() or ".." in relpath.parts:
                raise RuntimeError(f"unsafe evidence path in state: {jid}/{rel}")
            src = task / jid / relpath
            if not src.resolve().is_relative_to(task.resolve()) or not src.is_file() or sha(src) != digest:
                raise RuntimeError(f"source evidence missing or changed: {jid}/{rel}")
            dst = snap / jid / relpath
            if dst.exists():
                raise RuntimeError(f"refusing to overwrite snapshot file: {jid}/{rel}")
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, dst)
            if sha(dst) != digest:
                raise RuntimeError(f"copied evidence hash mismatch: {jid}/{rel}")
            manifest[f"{jid}/{rel}"] = digest
    receipt = {"utc": datetime.now(timezone.utc).isoformat(), "sha256": manifest}
    with (snap / "snapshot-receipt.json").open("x") as f:
        json.dump(receipt, f, indent=2, sort_keys=True)
        f.write("\n")
    return snap


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("campaign_task", type=Path)
    ap.add_argument("new_output_dir", type=Path)
    ap.add_argument("--max-wait-seconds", type=int, default=43200)
    args = ap.parse_args()
    task, out = args.campaign_task.resolve(), args.new_output_dir.resolve()
    if args.max_wait_seconds < 0 or args.max_wait_seconds > 43200:
        ap.error("--max-wait-seconds must be between 0 and 43200")
    if out.exists():
        raise SystemExit("output directory must be fresh")
    out.mkdir(parents=True)
    started = time.time()
    last = None
    try:
        while True:
            state = load_state(task / "state.json")
            status = state.get("status")
            if status != last:
                print(f"campaign status: {status}", flush=True)
                last = status
            if status in ("completed", "failed"):
                break
            if time.time() - started >= args.max_wait_seconds:
                raise TimeoutError("campaign did not reach terminal state before wait bound")
            time.sleep(min(30, max(1, args.max_wait_seconds - int(time.time() - started))))
        if status == "failed":
            print("campaign failed; preserving completed-job evidence", flush=True)
            snap_receipt = None
            snapshot_error = None
            try:
                snap = snapshot(task, out, state)
                snap_receipt = sha(snap / "snapshot-receipt.json")
            except Exception as e:
                snapshot_error = str(e)
            failure = {"status": "campaign_failed", "task_id": state.get("task_id"),
                       "state_sha256": sha(task / "state.json"),
                       "snapshot_receipt_sha256": snap_receipt, "snapshot_error": snapshot_error,
                       "thread_cap": state.get("thread_cap"),
                       "jobs_over_50_threads": [j.get("id") for j in state.get("jobs", [])
                           if j.get("peak_observed_threads") is not None and j["peak_observed_threads"] > 50],
                       "utc": datetime.now(timezone.utc).isoformat()}
            with (out / "failure-manifest.json").open("x") as f:
                json.dump(failure, f, indent=2, sort_keys=True); f.write("\n")
            return 2
        print("campaign completed; snapshotting verified evidence", flush=True)
        require_completed_state(state)
        snap = snapshot(task, out, state)
        collector = Path(__file__).with_name("collect_component_profile.py")
        output_json = out / "component-profile-summary.json"
        env = os.environ.copy()
        for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
                     "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "OMP_THREAD_LIMIT"):
            env[name] = "1"
        print("running data-only component collector", flush=True)
        completed = subprocess.run([sys.executable, str(collector), str(snap), str(output_json)],
                                   env=env, check=False)
        if completed.returncode:
            raise RuntimeError(f"collector failed with exit code {completed.returncode}")
        summary = json.loads(output_json.read_text())
        if summary.get("status") != "PASS" or summary.get("pending_jobs"):
            raise RuntimeError(f"collector did not produce a complete passing summary: {summary.get('status')}")
        receipt = {"status": "completed", "task_id": state.get("task_id"),
                   "snapshot_receipt_sha256": sha(snap / "snapshot-receipt.json"),
                   "summary_sha256": sha(output_json), "elapsed_seconds": time.time() - started,
                   "utc": datetime.now(timezone.utc).isoformat()}
        with (out / "finalization-receipt.json").open("x") as f:
            json.dump(receipt, f, indent=2, sort_keys=True); f.write("\n")
        return 0
    except Exception as e:
        error = {"status": "finalization_failed", "error": str(e),
                 "utc": datetime.now(timezone.utc).isoformat()}
        with (out / "finalization-failure.json").open("x") as f:
            json.dump(error, f, indent=2, sort_keys=True); f.write("\n")
        print(f"finalization failed: {e}", file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
