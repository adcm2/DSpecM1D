#!/usr/bin/env python3
"""Wait for verified transfer/build prerequisites, then exec the reviewed plan."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import time

WAIT_SECONDS = 4 * 60 * 60
POLL_SECONDS = 10
EXPECTED_VARIANTS = {
    "dsv4-native", "dsv4-profiled", "dsv6-native", "dsv6-profiled",
    "ysv6-native", "ysv6-profiled", "receiver-gate",
}
REVIEWED_TASK_ID = "profile-01"
TERMINAL_BUILD_FAILURES = {"failed", "unclean", "stopped_before_next_job", "interrupted"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"cannot read prerequisite JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise RuntimeError(f"prerequisite JSON is not an object: {path}")
    return value


def require_digest(value: str, name: str) -> str:
    digest = value.lower()
    if not re.fullmatch(r"[0-9a-f]{64}", digest):
        raise ValueError(f"{name} must be a 64-character SHA-256 hex digest")
    return digest


def check_statuses(root: Path) -> tuple[bool, str, str]:
    state_path = root / "continuation" / "build-03" / "state.json"
    transfer_path = root / "transfer" / "transfer-verification.json"
    build_status = "waiting"
    transfer_status = "waiting"
    if state_path.is_file():
        state = read_json(state_path)
        build_status = state.get("status")
        if not isinstance(build_status, str):
            raise RuntimeError(f"build state has no string status: {state_path}")
        if build_status in TERMINAL_BUILD_FAILURES or (build_status not in {"running", "completed"}):
            raise RuntimeError(f"build-03 reached terminal non-success status {build_status!r}")
    if transfer_path.is_file():
        try:
            transfer = json.loads(transfer_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            # The deployed extractor opens this marker exclusively, then writes
            # JSON directly; an observed partial write is still only pending.
            transfer_status = "writing"
        except FileNotFoundError:
            transfer_status = "waiting"
        else:
            if not isinstance(transfer, dict):
                raise RuntimeError(f"transfer verification marker is not an object: {transfer_path}")
            transfer_status = transfer.get("status", "writing")
            if transfer_status == "writing":
                pass
            elif transfer_status != "verified":
                raise RuntimeError(f"transfer verification has non-success status {transfer_status!r}")
    return build_status == "completed" and transfer_status == "verified", str(build_status), str(transfer_status)


def verify_build_receipt(root: Path, input_root: Path, plan_value: dict) -> dict:
    receipt_path = root / "continuation" / "build-03" / "build" / "build-result.json"
    report = read_json(receipt_path)
    if report.get("status") != "completed":
        raise RuntimeError(f"build receipt is not completed: {receipt_path}")
    build_state = read_json(root / "continuation" / "build-03" / "state.json")
    if build_state.get("status") != "completed":
        raise RuntimeError("build-03 supervisor state is not completed")
    build_jobs = [job for job in build_state.get("jobs", []) if isinstance(job, dict) and job.get("id") == "build"]
    if len(build_jobs) != 1:
        raise RuntimeError("build-03 supervisor state must contain exactly one build job")
    build_job = build_jobs[0]
    receipt_hash = sha256_file(receipt_path)
    if build_job.get("status") != "completed" or build_job.get("inputs_unchanged") is not True:
        raise RuntimeError("build-03 build job is not completed with unchanged inputs")
    if build_job.get("outputs_sha256", {}).get("build-result.json") != receipt_hash:
        raise RuntimeError("build-result.json hash does not match the completed supervisor state")

    build_config = read_json(input_root / "configs" / "build.json")
    configured_root_raw = build_config.get("build_root")
    if not isinstance(configured_root_raw, str) or not Path(configured_root_raw).is_absolute():
        raise RuntimeError("reviewed build config must name its absolute build_root")
    configured_root = Path(configured_root_raw).resolve(strict=False)
    build_root_raw = report.get("build_root")
    if not isinstance(build_root_raw, str) or not Path(build_root_raw).is_absolute():
        raise RuntimeError("build receipt must name its absolute build_root")
    build_root = Path(build_root_raw).resolve(strict=True)
    if build_root != configured_root:
        raise RuntimeError(f"build receipt root differs from reviewed build config: {build_root} != {configured_root}")
    config_variants = build_config.get("variants")
    if not isinstance(config_variants, list) or len(config_variants) != len(EXPECTED_VARIANTS):
        raise RuntimeError("reviewed build config must contain exactly seven variants")
    expected_outputs = {}
    for variant in config_variants:
        if not isinstance(variant, dict) or not isinstance(variant.get("id"), str):
            raise RuntimeError("malformed reviewed build variant")
        variant_id = variant["id"]
        output_name = variant.get("output_name", "program")
        if not isinstance(output_name, str) or Path(output_name).name != output_name:
            raise RuntimeError(f"invalid configured output name for {variant_id}")
        expected_outputs[variant_id] = (build_root / variant_id / output_name).resolve(strict=False)
    if len(expected_outputs) != len(config_variants) or set(expected_outputs) != EXPECTED_VARIANTS:
        raise RuntimeError("reviewed build config does not define exactly the seven expected variants")
    variants = report.get("variants")
    if not isinstance(variants, list) or len(variants) != len(EXPECTED_VARIANTS):
        raise RuntimeError("build receipt must contain exactly seven variants")
    ids = [row.get("id") for row in variants if isinstance(row, dict)]
    if len(ids) != len(variants) or set(ids) != EXPECTED_VARIANTS:
        raise RuntimeError(f"build receipt variant IDs do not match the seven reviewed outputs: {ids}")
    checked = {}
    for row in variants:
        variant_id = row["id"]
        expected = row.get("binary_sha256")
        if not isinstance(expected, str) or not re.fullmatch(r"[0-9a-f]{64}", expected):
            raise RuntimeError(f"build receipt has invalid binary SHA-256 for {variant_id}")
        command_path = build_root / variant_id / "command.json"
        command = read_json(command_path)
        if (row.get("returncode") != 0 or command.get("id") != variant_id or
                command.get("returncode") != 0 or command.get("tracked_inputs_unchanged") is not True or
                command.get("binary_sha256") != expected):
            raise RuntimeError(f"build command receipt does not match summary for {variant_id}")
        output_raw = command.get("output")
        if not isinstance(output_raw, str) or not Path(output_raw).is_absolute():
            raise RuntimeError(f"build command receipt has no absolute binary output for {variant_id}")
        binary = Path(output_raw).resolve(strict=True)
        expected_parent = (build_root / variant_id).resolve(strict=True)
        if binary.parent != expected_parent or binary != expected_outputs[variant_id]:
            raise RuntimeError(f"build output differs from reviewed build config for {variant_id}: {binary}")
        if not binary.is_file():
            raise RuntimeError(f"built output is not a file: {binary}")
        actual = sha256_file(binary)
        if actual != expected:
            raise RuntimeError(f"built binary hash changed for {variant_id}: {actual} != {expected}")
        checked[variant_id] = {"path": str(binary), "sha256": actual}
    used_build_paths = set()
    jobs = plan_value.get("jobs")
    if not isinstance(jobs, list):
        raise RuntimeError("reviewed qualification plan has no jobs list")
    for job in jobs:
        if not isinstance(job, dict):
            raise RuntimeError("malformed qualification job")
        argv = job.get("command", [])
        if not isinstance(argv, list):
            raise RuntimeError("qualification job command must be an argument list")
        inputs = job.get("inputs", [])
        if not isinstance(inputs, list):
            raise RuntimeError("qualification job inputs must be a list")
        argv_paths = [Path(token).resolve(strict=False) for token in argv
                      if isinstance(token, str) and Path(token).is_absolute()]
        if argv and argv[0] == "/usr/bin/taskset":
            if len(argv) < 4 or argv[1:3] != ["-c", "0"] or not Path(argv[3]).is_absolute():
                raise RuntimeError(f"build-backed taskset command does not use reviewed CPU-0 format: {job.get('id')}")
            executable = Path(argv[3]).resolve(strict=False)
        elif argv and isinstance(argv[0], str) and Path(argv[0]).is_absolute():
            executable = Path(argv[0]).resolve(strict=False)
        else:
            executable = None
        for token in argv + inputs:
            if isinstance(token, str) and Path(token).is_absolute():
                path = Path(token).resolve(strict=False)
                try:
                    path.relative_to(build_root)
                except ValueError:
                    continue
                if path not in set(expected_outputs.values()):
                    raise RuntimeError(f"plan references an unverified file under build_root: {path}")
        if executable in set(expected_outputs.values()):
            used_build_paths.add(executable)
    missing_plan_binaries = set(expected_outputs.values()) - used_build_paths
    if missing_plan_binaries:
        raise RuntimeError(f"reviewed plan does not reference every verified build output: {sorted(map(str, missing_plan_binaries))}")
    plan = input_root / "configs" / "qualification-and-profile-plan.json"
    checker = input_root / "tools" / "qualification_checks.py"
    supervisor = input_root / "tools" / "supervise_campaign.py"
    if not plan.is_file() or not checker.is_file() or not supervisor.is_file():
        raise RuntimeError("reviewed plan, checker, or supervisor is missing from INPUT_ROOT")
    if Path(str(plan_value.get("campaign_root", ""))).resolve(strict=False) != root:
        raise RuntimeError("qualification plan campaign_root does not match ROOT")
    return {"build_receipt": str(receipt_path), "build_root": str(build_root),
            "verified_binaries": checked, "plan_path": str(plan),
            "checker_path": str(checker), "supervisor_path": str(supervisor)}


def verify_launch(root: Path, input_root: Path, expected_checker: str,
                  expected_plan: str, task_id: str) -> dict:
    if task_id != REVIEWED_TASK_ID:
        raise RuntimeError(f"this reviewed plan supports only task ID {REVIEWED_TASK_ID!r}")
    checker = input_root / "tools" / "qualification_checks.py"
    plan = input_root / "configs" / "qualification-and-profile-plan.json"
    supervisor = input_root / "tools" / "supervise_campaign.py"
    if sha256_file(checker) != expected_checker:
        raise RuntimeError(f"checker SHA-256 does not match reviewed value: {checker}")
    if sha256_file(plan) != expected_plan:
        raise RuntimeError(f"plan SHA-256 does not match reviewed value: {plan}")
    if not supervisor.is_file():
        raise RuntimeError(f"supervisor script is missing: {supervisor}")
    task_dir = root / "continuation" / task_id
    if task_dir.exists():
        raise RuntimeError(f"refusing to reuse existing campaign task directory: {task_dir}")
    plan_value = read_json(plan)
    if Path(str(plan_value.get("campaign_root", ""))).resolve(strict=False) != root:
        raise RuntimeError("qualification plan campaign_root does not match ROOT")
    build = verify_build_receipt(root, input_root, plan_value)
    return {"task_dir": str(task_dir), **build}


def wait_for_prerequisites(root: Path, timeout_seconds: int = WAIT_SECONDS) -> None:
    deadline = time.monotonic() + timeout_seconds
    previous = None
    while True:
        ready, build_status, transfer_status = check_statuses(root)
        phase = (build_status, transfer_status)
        if phase != previous:
            print(f"wait_and_launch: build-03={build_status}; transfer={transfer_status}", flush=True)
            previous = phase
        if ready:
            return
        if time.monotonic() >= deadline:
            raise TimeoutError("build and verified transfer prerequisites did not both complete within four hours")
        time.sleep(POLL_SECONDS)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ROOT", type=Path)
    parser.add_argument("INPUT_ROOT", type=Path)
    parser.add_argument("EXPECTED_CHECKER_SHA256")
    parser.add_argument("EXPECTED_PLAN_SHA256")
    parser.add_argument("TASKID", nargs="?", default="profile-01")
    args = parser.parse_args()
    if not args.ROOT.is_absolute() or not args.INPUT_ROOT.is_absolute():
        parser.error("ROOT and INPUT_ROOT must be absolute paths")
    root = args.ROOT.resolve(strict=True)
    input_root = args.INPUT_ROOT.resolve(strict=True)
    expected_checker = require_digest(args.EXPECTED_CHECKER_SHA256, "expected checker SHA-256")
    expected_plan = require_digest(args.EXPECTED_PLAN_SHA256, "expected plan SHA-256")
    if args.TASKID != REVIEWED_TASK_ID:
        parser.error(f"this reviewed plan supports only TASKID {REVIEWED_TASK_ID}")

    wait_for_prerequisites(root)
    # Recheck both status files after the bounded wait and before inspecting/launching.
    ready, build_status, transfer_status = check_statuses(root)
    if not ready:
        raise RuntimeError(f"prerequisites changed before launch: build={build_status}; transfer={transfer_status}")
    verified = verify_launch(root, input_root, expected_checker, expected_plan, args.TASKID)
    print("wait_and_launch: prerequisites and all seven binary hashes verified; exec supervisor", flush=True)
    command = ["/usr/bin/python3", verified["supervisor_path"], verified["plan_path"],
               "--task-id", args.TASKID]
    os.execv(command[0], command)
    return 127


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"wait_and_launch: {type(exc).__name__}: {exc}", file=sys.stderr, flush=True)
        raise SystemExit(2)
