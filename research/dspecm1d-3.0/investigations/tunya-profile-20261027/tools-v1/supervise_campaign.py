#!/usr/bin/env python3
"""Run reviewed jobs sequentially with durable, identity-aware supervision.

Linux /proc polling provides a userspace descendant/thread watchdog; it is not
kernel aggregate accounting and can miss processes shorter than the poll interval.
Start detached with nohup as documented in README.md.
"""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import uuid

POLL_SECONDS = 0.10
THREAD_VARS = ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "OPENBLAS_NUM_THREADS",
               "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
               "NUMEXPR_NUM_THREADS", "GOTO_NUM_THREADS")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def proc_table() -> dict[int, tuple[int, int, int]]:
    """pid -> (ppid, starttime ticks, thread count), validating /proc identity."""
    table: dict[int, tuple[int, int, int]] = {}
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        try:
            raw = (entry / "stat").read_text()
            tail = raw[raw.rfind(")") + 2:].split()
            ppid = int(tail[1])
            start = int(tail[19])
            threads = int(tail[17])
            table[int(entry.name)] = (ppid, start, threads)
        except (OSError, ValueError, IndexError):
            continue
    return table


def descendants(table: dict[int, tuple[int, int, int]], roots: set[tuple[int, int]]) -> set[tuple[int, int]]:
    found = set(roots)
    changed = True
    while changed:
        changed = False
        parent_pids = {pid for pid, _ in found}
        for pid, (ppid, start, _) in table.items():
            ident = (pid, start)
            if ppid in parent_pids and ident not in found:
                found.add(ident); changed = True
    return found


def alive_threads(table: dict[int, tuple[int, int, int]], ids: set[tuple[int, int]]) -> tuple[int, list[dict[str, int]]]:
    total, current = 0, []
    for pid, start in ids:
        row = table.get(pid)
        if row and row[1] == start:
            total += row[2]
            current.append({"pid": pid, "starttime_ticks": start, "threads": row[2]})
    return total, current


def signal_owned(ids: set[tuple[int, int]], sig: int) -> None:
    for pid, start in sorted(ids, reverse=True):
        try:
            raw = Path(f"/proc/{pid}/stat").read_text()
            current_start = int(raw[raw.rfind(")") + 2:].split()[19])
            if current_start == start:
                os.kill(pid, sig)
        except (OSError, ValueError, IndexError):
            pass


def atomic_json(path: Path, obj: dict) -> None:
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    with tmp.open("w") as stream:
        stream.write(json.dumps(obj, indent=2, sort_keys=True) + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(tmp, path)
    dirfd = os.open(path.parent, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(dirfd)
    finally:
        os.close(dirfd)


def command_hash(command: list[str]) -> str:
    return hashlib.sha256(json.dumps(command, separators=(",", ":")).encode()).hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("plan", type=Path)
    ap.add_argument("--task-id", help="unique task folder name; defaults to a UUID")
    ap.add_argument("--resume", action="store_true", help="skip only completed jobs whose evidence verifies")
    args = ap.parse_args()
    plan = json.loads(args.plan.read_text())
    campaign = Path(plan["campaign_root"])
    if not campaign.is_absolute():
        raise ValueError("campaign_root must be absolute")
    task_id = args.task_id or time.strftime("%Y%m%dT%H%M%SZ", time.gmtime()) + "-" + uuid.uuid4().hex[:8]
    if not task_id.replace("-", "").replace("_", "").isalnum():
        raise ValueError("task id may contain only letters, digits, hyphens and underscores")
    task = campaign / "continuation" / task_id
    cap = int(plan.get("thread_cap", 50))
    if not 1 <= cap <= 50:
        raise ValueError("thread_cap must be between 1 and 50")
    jobs = plan["jobs"]
    if not jobs or len({j["id"] for j in jobs}) != len(jobs):
        raise ValueError("jobs must be nonempty with unique ids")

    campaign.mkdir(parents=True, exist_ok=True)
    (campaign / "continuation").mkdir(exist_ok=True)
    lock = (campaign / ".supervise_campaign.lock").open("a+")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise RuntimeError("another campaign supervisor holds the single-instance lock")
    if task.exists() and not args.resume:
        raise ValueError(f"task folder already exists; choose a new task id: {task}")
    if not task.exists() and args.resume:
        raise ValueError(f"cannot resume missing task: {task}")
    task.mkdir(parents=True, exist_ok=True)
    state_path = task / "state.json"
    state = json.loads(state_path.read_text()) if args.resume and state_path.exists() else {
        "schema": 1, "task_id": task_id, "plan_sha256": sha256_file(args.plan),
        "started_unix": time.time(), "thread_cap": cap, "status": "running", "jobs": []}
    if args.resume and state.get("plan_sha256") != sha256_file(args.plan):
        raise ValueError("plan hash differs from the original task")
    state["status"] = "running"
    atomic_json(state_path, state)
    env = os.environ.copy()
    env.update({str(k): str(v) for k, v in plan.get("environment", {}).items()})
    for key in THREAD_VARS:
        env[key] = "1"

    previous = {x["id"]: x for x in state["jobs"]}
    state["jobs"] = []
    for job in jobs:
        jid = job["id"]
        if not jid.replace("-", "").replace("_", "").isalnum():
            raise ValueError(f"job id may contain only letters, digits, hyphens and underscores: {jid}")
        jdir = task / jid
        cmd = [str(x) for x in job["command"]]
        cmd = [x.replace("{JOB_DIR}", str(jdir)).replace("{TASK_DIR}", str(task)) for x in cmd]
        chash = command_hash(cmd)
        old = previous.get(jid)
        if old and old.get("status") == "completed":
            if old.get("command_sha256") != chash:
                raise ValueError(f"completed job command changed: {jid}")
            log = jdir / "job.log"
            if not log.is_file() or sha256_file(log) != old.get("log_sha256"):
                raise ValueError(f"completed job log failed hash verification: {jid}")
            for rel, digest in old.get("outputs_sha256", {}).items():
                p = jdir / rel
                if not p.is_file() or sha256_file(p) != digest:
                    raise ValueError(f"completed job output failed hash verification: {jid}/{rel}")
            state["jobs"].append(old)
            atomic_json(state_path, state)
            continue
        if old:
            raise ValueError(f"job {jid} is not verified completed; preserve it and use a new task id")
        if jdir.exists():
            raise ValueError(f"job output directory already exists: {jdir}")
        jdir.mkdir()
        for rel in job.get("outputs", []):
            relpath = Path(rel)
            if relpath.is_absolute() or ".." in relpath.parts:
                raise ValueError(f"declared output must remain under its job directory: {rel}")
            p = jdir / relpath
            if p.exists():
                raise ValueError(f"declared output already exists: {p}")
            p.parent.mkdir(parents=True, exist_ok=True)
        cwd = Path(job["cwd"])
        if not cwd.is_absolute() or not cwd.is_dir():
            raise ValueError(f"cwd must be an existing absolute directory: {cwd}")
        rec = {"id": jid, "command": cmd, "command_sha256": chash, "cwd": str(cwd),
               "status": "running", "started_unix": time.time(), "log": "job.log",
               "outputs": job.get("outputs", [])}
        state["jobs"].append(rec); atomic_json(state_path, state)
        start = time.monotonic()
        proc = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, start_new_session=True, bufsize=0)
        os.set_blocking(proc.stdout.fileno(), False)
        root_stat = (Path(f"/proc/{proc.pid}/stat").read_text())
        root_start = int(root_stat[root_stat.rfind(")") + 2:].split()[19])
        owned = {(proc.pid, root_start)}
        reasons: list[str] = []
        last_write = time.monotonic()
        with (jdir / "job.log").open("wb", buffering=0) as log:
            while True:
                try:
                    chunk = os.read(proc.stdout.fileno(), 65536) if proc.stdout else b""
                except BlockingIOError:
                    chunk = None
                if chunk is None:
                    chunk = b""
                if chunk:
                    log.write(chunk)
                table = proc_table()
                owned |= descendants(table, owned)
                nthreads, thread_rows = alive_threads(table, owned | {(os.getpid(), int(Path(f"/proc/{os.getpid()}/stat").read_text().split(")",1)[1].split()[19]))})
                rss_kb = 0
                for pid, st in owned:
                    row = table.get(pid)
                    if not row or row[1] != st:
                        continue
                    try:
                        for line in Path(f"/proc/{pid}/status").read_text().splitlines():
                            if line.startswith("VmRSS:"):
                                rss_kb += int(line.split()[1]); break
                    except (OSError, ValueError, IndexError):
                        pass
                budget = job.get("memory_gib")
                if not reasons and nthreads > cap:
                    reasons.append(f"thread_limit_exceeded:{nthreads}>{cap}")
                if budget is not None and not reasons and rss_kb > float(budget) * 1024 * 1024:
                    reasons.append(f"memory_limit_exceeded:{rss_kb}KiB>{budget}GiB")
                timeout = job.get("timeout_seconds")
                if timeout is not None and not reasons and time.monotonic() - start > float(timeout):
                    reasons.append("timeout")
                if reasons:
                    signal_owned(owned, signal.SIGTERM)
                    time.sleep(0.5)
                    table = proc_table(); owned |= descendants(table, owned)
                    signal_owned(owned, signal.SIGKILL)
                    break
                if proc.poll() is not None and not chunk:
                    # A completed leader may leave workers behind; keep tracking until they exit.
                    remaining = alive_threads(proc_table(), owned)[0]
                    if remaining == 0:
                        break
                if time.monotonic() - last_write >= 2:
                    log.flush()
                    os.fsync(log.fileno())
                    rec["watchdog"] = {"elapsed_seconds": round(time.monotonic() - start, 2),
                                       "thread_count_including_supervisor": nthreads,
                                       "thread_cap": cap, "owned_processes": thread_rows,
                                       "aggregate_rss_gib_userspace": round(rss_kb / (1024 * 1024), 4)}
                    rec["updated_unix"] = time.time(); atomic_json(state_path, state)
                    last_write = time.monotonic()
                if not chunk:
                    time.sleep(POLL_SECONDS)
            log.flush()
            os.fsync(log.fileno())
        returncode = proc.wait()
        log_hash = sha256_file(jdir / "job.log")
        outputs_hashes = {}
        missing = []
        for rel in job.get("outputs", []):
            p = jdir / rel
            if p.is_file():
                outputs_hashes[rel] = sha256_file(p)
            else:
                missing.append(rel)
        rec.update({"finished_unix": time.time(), "returncode": returncode,
                    "status": "completed" if returncode == 0 and not reasons and not missing else "failed",
                    "failure_reasons": reasons + (["missing_outputs:" + ",".join(missing)] if missing else []),
                    "log_sha256": log_hash, "outputs_sha256": outputs_hashes,
                    "final_owned_processes": sorted([list(x) for x in owned])})
        atomic_json(state_path, state)
        if rec["status"] != "completed":
            state["status"] = "failed"; state["finished_unix"] = time.time()
            atomic_json(state_path, state)
            return returncode or 1
    state["status"] = "completed"; state["finished_unix"] = time.time()
    atomic_json(state_path, state)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"supervise_campaign: {exc}", file=sys.stderr)
        raise SystemExit(2)
