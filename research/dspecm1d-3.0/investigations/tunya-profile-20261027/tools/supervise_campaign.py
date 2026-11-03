#!/usr/bin/env python3
"""Run reviewed jobs sequentially with durable, identity-aware supervision.

Linux /proc polling provides a userspace descendant/thread watchdog; it is not
kernel aggregate accounting and can miss processes shorter than the poll interval.
Start detached with nohup as documented in README.md.
"""
from __future__ import annotations

import argparse
import ctypes
import errno
import fcntl
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import sys
import time
import uuid

POLL_SECONDS = 0.10
CLEANUP_TERM_SECONDS = 0.5
CLEANUP_TOTAL_SECONDS = 10.0
THREAD_VARS = ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "OPENBLAS_NUM_THREADS",
               "MKL_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
               "NUMEXPR_NUM_THREADS", "GOTO_NUM_THREADS")
REQUESTED_STOP: int | None = None


def request_stop(signum, _frame) -> None:
    global REQUESTED_STOP
    REQUESTED_STOP = signum


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def input_manifest(job: dict, cwd: Path, command: list[str], job_dir: Path) -> dict[str, str]:
    declared = job.get("inputs")
    if not isinstance(declared, list) or not declared:
        raise ValueError(f"job {job.get('id')} must declare input files for hash-bound resume")
    paths = set()
    outputs = {(job_dir / rel).resolve(strict=False) for rel in job.get("outputs", [])}
    for value in declared:
        path = Path(value)
        if not path.is_absolute():
            path = cwd / path
        path = path.resolve(strict=True)
        if path in outputs:
            raise ValueError(f"job input is also declared as an output: {path}")
        paths.add(path)
    # Include command-line executables and file-valued arguments such as taskset's target.
    executable = shutil.which(command[0]) if command and "/" not in command[0] else command[0]
    if executable:
        candidate = Path(executable)
        if candidate.is_file():
            paths.add(candidate.resolve(strict=True))
    for token in command[1:]:
        candidate = Path(token)
        if not candidate.is_absolute():
            candidate = cwd / candidate
        candidate = candidate.resolve(strict=False)
        if candidate not in outputs and candidate.is_file():
            paths.add(candidate.resolve(strict=True))
    return {str(path): sha256_file(path) for path in sorted(paths)}


def enable_child_subreaper() -> None:
    """Fail closed unless Linux adopts orphaned descendants into this process."""
    if not sys.platform.startswith("linux"):
        raise RuntimeError("Linux PR_SET_CHILD_SUBREAPER is required")
    libc = ctypes.CDLL(None, use_errno=True)
    if libc.prctl(36, 1, 0, 0, 0) != 0:  # PR_SET_CHILD_SUBREAPER
        err = ctypes.get_errno()
        raise OSError(err, os.strerror(err), "PR_SET_CHILD_SUBREAPER")
    value = ctypes.c_int()
    if libc.prctl(37, ctypes.byref(value), 0, 0, 0) != 0 or value.value != 1:  # GET
        err = ctypes.get_errno()
        raise OSError(err or errno.EIO, "could not verify child subreaper state")


def proc_table() -> tuple[dict[int, tuple[int, int, int]], list[str]]:
    """Return pid -> (ppid, starttime, thread count) plus inspection failures."""
    table: dict[int, tuple[int, int, int]] = {}
    issues: list[str] = []
    try:
        entries = list(Path("/proc").iterdir())
    except OSError as exc:
        return {}, [f"/proc enumeration: {type(exc).__name__}:{exc.errno}"]
    for entry in entries:
        if not entry.name.isdigit():
            continue
        try:
            raw = (entry / "stat").read_text()
            tail = raw[raw.rfind(")") + 2:].split()
            ppid = int(tail[1])
            start = int(tail[19])
            threads = int(tail[17])
            table[int(entry.name)] = (ppid, start, threads)
        except OSError as exc:
            if exc.errno not in (errno.ENOENT, errno.ESRCH):
                issues.append(f"{entry.name}: {type(exc).__name__}:{exc.errno}")
        except (ValueError, IndexError) as exc:
            issues.append(f"{entry.name}: malformed stat:{type(exc).__name__}")
    return table, issues


def descendants(table: dict[int, tuple[int, int, int]], roots: set[tuple[int, int]]) -> set[tuple[int, int]]:
    found = set(roots)
    changed = True
    while changed:
        changed = False
        parent_pids = {pid for pid, start in found
                       if pid in table and table[pid][1] == start}
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


def identity_state(ident: tuple[int, int]) -> str:
    pid, start = ident
    try:
        raw = Path(f"/proc/{pid}/stat").read_text()
        current_start = int(raw[raw.rfind(")") + 2:].split()[19])
    except OSError as exc:
        if exc.errno in (errno.ENOENT, errno.ESRCH):
            return "gone"
        return f"unknown:{type(exc).__name__}:{exc.errno}"
    except (ValueError, IndexError):
        return "unknown:malformed_stat"
    return "alive" if current_start == start else "reused"


def signal_identity(ident: tuple[int, int], sig: int) -> str:
    if not hasattr(os, "pidfd_open") or not hasattr(signal, "pidfd_send_signal"):
        return "error:pidfd_signaling_unavailable"
    try:
        pidfd = os.pidfd_open(ident[0], 0)
    except ProcessLookupError:
        return "gone"
    except OSError as exc:
        return f"error:{type(exc).__name__}:{exc.errno}"
    try:
        state = identity_state(ident)
        if state != "alive":
            return state
        signal.pidfd_send_signal(pidfd, sig, None, 0)
        return "sent"
    except ProcessLookupError:
        return "gone"
    except OSError as exc:
        return f"error:{type(exc).__name__}:{exc.errno}"
    finally:
        os.close(pidfd)


def reap_adopted(ids: set[tuple[int, int]], keep_pid: int) -> list[int]:
    reaped = []
    for pid, _ in ids:
        if pid == keep_pid:
            continue
        try:
            child, _ = os.waitpid(pid, os.WNOHANG)
            if child:
                reaped.append(child)
        except ChildProcessError:
            pass
        except InterruptedError:
            pass
    return reaped


def cleanup_owned(ids: set[tuple[int, int]], supervisor: tuple[int, int],
                  leader: subprocess.Popen) -> dict:
    """Terminate, reap, and verify every recorded identity gone or report unclean."""
    signal_log = []
    reaped = set()
    inspection_errors = []
    term_deadline = time.monotonic() + CLEANUP_TERM_SECONDS
    end_deadline = time.monotonic() + CLEANUP_TOTAL_SECONDS
    escalated = False
    sent_term: set[tuple[int, int]] = set()
    sent_kill: set[tuple[int, int]] = set()
    while True:
        table, issues = proc_table()
        inspection_errors.extend(x for x in issues if x not in inspection_errors)
        ids |= descendants(table, ids | {supervisor}) - {supervisor}
        statuses = {ident: identity_state(ident) for ident in ids}
        alive = [ident for ident, status in statuses.items() if status == "alive"]
        unknown = [ident for ident, status in statuses.items() if status.startswith("unknown:")]
        now = time.monotonic()
        sig = signal.SIGKILL if escalated or now >= term_deadline else signal.SIGTERM
        if sig == signal.SIGKILL:
            escalated = True
        already_sent = sent_kill if sig == signal.SIGKILL else sent_term
        for ident in alive:
            if ident in already_sent:
                continue
            outcome = signal_identity(ident, sig)
            signal_log.append({"identity": list(ident), "signal": int(sig), "result": outcome})
            already_sent.add(ident)
        # Popen owns/reaps the direct leader; adopted descendants are reaped below.
        leader.poll()
        reaped.update(reap_adopted(ids, leader.pid))
        if not alive and not unknown:
            return {"verified_gone": True, "remaining": [], "unknown": [],
                    "signal_log": signal_log, "reaped_pids": sorted(reaped),
                    "proc_inspection_errors": inspection_errors}
        if now >= end_deadline:
            final = {ident: identity_state(ident) for ident in ids}
            remaining = [list(ident) for ident, status in final.items() if status == "alive"]
            unknown_final = [{"identity": list(ident), "state": status}
                             for ident, status in final.items() if status.startswith("unknown:")]
            return {"verified_gone": not remaining and not unknown_final,
                    "remaining": remaining, "unknown": unknown_final,
                    "signal_log": signal_log, "reaped_pids": sorted(reaped),
                    "proc_inspection_errors": inspection_errors}
        time.sleep(POLL_SECONDS)


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
    signal.signal(signal.SIGTERM, request_stop)
    signal.signal(signal.SIGINT, request_stop)
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

    state_path = task / "state.json"
    prior_state = None
    if args.resume:
        if not task.is_dir():
            raise ValueError(f"cannot resume missing task: {task}")
        if not state_path.is_file():
            raise ValueError(f"cannot resume: task state is missing: {state_path}")
        try:
            prior_state = json.loads(state_path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError(f"cannot resume: task state is unreadable or malformed: {state_path}") from exc
        if not isinstance(prior_state, dict) or prior_state.get("schema") != 1:
            raise ValueError("cannot resume: unsupported or malformed state schema")
        if prior_state.get("task_id") != task_id:
            raise ValueError("cannot resume: state task_id does not match requested task id")
        if prior_state.get("plan_sha256") != sha256_file(args.plan):
            raise ValueError("plan hash differs from the original task")
        if not isinstance(prior_state.get("jobs"), list):
            raise ValueError("cannot resume: state jobs field is malformed")
        recorded = {entry.get("id"): entry for entry in prior_state["jobs"]
                    if isinstance(entry, dict) and isinstance(entry.get("id"), str)}
        if len(recorded) != len(prior_state["jobs"]):
            raise ValueError("cannot resume: state contains malformed or duplicate job records")
        for job in jobs:
            old = recorded.get(job["id"])
            if old is None:
                continue
            if old.get("status") != "completed":
                raise ValueError(f"cannot resume incomplete job {job['id']}; preserve it and use a new task id")
            jdir = task / job["id"]
            cwd = Path(job["cwd"])
            if not cwd.is_absolute() or not cwd.is_dir():
                raise ValueError(f"cwd must be an existing absolute directory: {cwd}")
            cmd = [str(x).replace("{JOB_DIR}", str(jdir)).replace("{TASK_DIR}", str(task))
                   for x in job["command"]]
            if old.get("command_sha256") != command_hash(cmd):
                raise ValueError(f"completed job command changed: {job['id']}")
            current_inputs = input_manifest(job, cwd, cmd, jdir)
            if old.get("inputs_sha256_before") != old.get("inputs_sha256_after") or old.get("inputs_sha256_after") != current_inputs:
                raise ValueError(f"completed job inputs changed: {job['id']}")
            log = jdir / "job.log"
            if not log.is_file() or sha256_file(log) != old.get("log_sha256"):
                raise ValueError(f"completed job log failed hash verification: {job['id']}")
            for rel, digest in old.get("outputs_sha256", {}).items():
                output = jdir / rel
                if not output.is_file() or sha256_file(output) != digest:
                    raise ValueError(f"completed job output failed hash verification: {job['id']}/{rel}")

    enable_child_subreaper()
    self_stat = Path(f"/proc/{os.getpid()}/stat").read_text()
    supervisor_start = int(self_stat[self_stat.rfind(")") + 2:].split()[19])
    supervisor_identity = (os.getpid(), supervisor_start)
    campaign.mkdir(parents=True, exist_ok=True)
    (campaign / "continuation").mkdir(exist_ok=True)
    lock = (campaign / ".supervise_campaign.lock").open("a+")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise RuntimeError("another campaign supervisor holds the single-instance lock")
    if task.exists() and not args.resume:
        raise ValueError(f"task folder already exists; choose a new task id: {task}")
    task.mkdir(parents=True, exist_ok=True)
    state = prior_state if args.resume else {
        "schema": 1, "task_id": task_id, "plan_sha256": sha256_file(args.plan),
        "started_unix": time.time(), "thread_cap": cap, "status": "running", "jobs": []}
    state["status"] = "running"
    atomic_json(state_path, state)
    env = os.environ.copy()
    env.update({str(k): str(v) for k, v in plan.get("environment", {}).items()})
    for key in THREAD_VARS:
        env[key] = "1"

    previous = {x["id"]: x for x in state["jobs"]}
    if not args.resume:
        state["jobs"] = []
    for job in jobs:
        if REQUESTED_STOP is not None:
            state["status"] = "stopped_before_next_job"
            state["stop_signal"] = REQUESTED_STOP
            state["finished_unix"] = time.time()
            atomic_json(state_path, state)
            return 128 + REQUESTED_STOP
        jid = job["id"]
        if not jid.replace("-", "").replace("_", "").isalnum():
            raise ValueError(f"job id may contain only letters, digits, hyphens and underscores: {jid}")
        jdir = task / jid
        cmd = [str(x) for x in job["command"]]
        cmd = [x.replace("{JOB_DIR}", str(jdir)).replace("{TASK_DIR}", str(task)) for x in cmd]
        chash = command_hash(cmd)
        cwd = Path(job["cwd"])
        if not cwd.is_absolute() or not cwd.is_dir():
            raise ValueError(f"cwd must be an existing absolute directory: {cwd}")
        current_inputs = input_manifest(job, cwd, cmd, jdir)
        old = previous.get(jid)
        if old and old.get("status") == "completed":
            if old.get("command_sha256") != chash:
                raise ValueError(f"completed job command changed: {jid}")
            if old.get("inputs_sha256_before") != old.get("inputs_sha256_after"):
                raise ValueError(f"completed job inputs changed during its original run: {jid}")
            if old.get("inputs_sha256_after") != current_inputs:
                raise ValueError(f"completed job input files changed since the original run: {jid}")
            log = jdir / "job.log"
            if not log.is_file() or sha256_file(log) != old.get("log_sha256"):
                raise ValueError(f"completed job log failed hash verification: {jid}")
            for rel, digest in old.get("outputs_sha256", {}).items():
                p = jdir / rel
                if not p.is_file() or sha256_file(p) != digest:
                    raise ValueError(f"completed job output failed hash verification: {jid}/{rel}")
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
        rec = {"id": jid, "command": cmd, "command_sha256": chash, "cwd": str(cwd),
               "status": "running", "started_unix": time.time(), "log": "job.log",
               "outputs": job.get("outputs", []), "inputs_sha256_before": current_inputs}
        state["jobs"].append(rec); atomic_json(state_path, state)
        start = time.monotonic()
        usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
        proc = subprocess.Popen(cmd, cwd=cwd, env=env, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, start_new_session=True, bufsize=0)
        os.set_blocking(proc.stdout.fileno(), False)
        root_stat = Path(f"/proc/{proc.pid}/stat").read_text()
        root_start = int(root_stat[root_stat.rfind(")") + 2:].split()[19])
        owned = {(proc.pid, root_start)}
        reasons: list[str] = []
        cleanup: dict | None = None
        last_write = time.monotonic()
        peak_threads = 0
        peak_rss_kb = 0
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
                table, proc_issues = proc_table()
                owned |= descendants(table, owned | {supervisor_identity}) - {supervisor_identity}
                nthreads, thread_rows = alive_threads(table, owned | {supervisor_identity})
                rss_kb = 0
                for pid, st in owned | {supervisor_identity}:
                    row = table.get(pid)
                    if not row or row[1] != st:
                        continue
                    try:
                        for line in Path(f"/proc/{pid}/status").read_text().splitlines():
                            if line.startswith("VmRSS:"):
                                rss_kb += int(line.split()[1]); break
                    except (OSError, ValueError, IndexError):
                        pass
                peak_threads = max(peak_threads, nthreads)
                peak_rss_kb = max(peak_rss_kb, rss_kb)
                budget = job.get("memory_gib")
                if not reasons and nthreads > cap:
                    reasons.append(f"thread_limit_exceeded:{nthreads}>{cap}")
                if budget is not None and not reasons and rss_kb > float(budget) * 1024 * 1024:
                    reasons.append(f"memory_limit_exceeded:{rss_kb}KiB>{budget}GiB")
                timeout = job.get("timeout_seconds")
                if timeout is not None and not reasons and time.monotonic() - start > float(timeout):
                    reasons.append("timeout")
                if proc_issues and not reasons:
                    reasons.append("proc_inspection_errors:" + ";".join(proc_issues[:20]))
                if REQUESTED_STOP is not None and not reasons:
                    reasons.append(f"supervisor_signal:{REQUESTED_STOP}")
                if reasons:
                    cleanup = cleanup_owned(owned, supervisor_identity, proc)
                    break
                if proc.poll() is not None and not chunk:
                    # A completed leader may leave workers behind; keep tracking until they exit.
                    remaining_table, _remaining_issues = proc_table()
                    remaining = alive_threads(remaining_table, owned)[0]
                    if remaining == 0:
                        break
                if time.monotonic() - last_write >= 2:
                    log.flush()
                    os.fsync(log.fileno())
                    rec["watchdog"] = {"elapsed_seconds": round(time.monotonic() - start, 2),
                                       "thread_count_including_supervisor": nthreads,
                                       "thread_cap": cap, "owned_processes": thread_rows,
                                       "aggregate_rss_gib_userspace": round(rss_kb / (1024 * 1024), 4),
                                       "peak_observed_threads": peak_threads,
                                       "peak_observed_rss_gib_userspace": round(peak_rss_kb / (1024 * 1024), 4)}
                    rec["updated_unix"] = time.time(); atomic_json(state_path, state)
                    last_write = time.monotonic()
                if not chunk:
                    time.sleep(POLL_SECONDS)
            log.flush()
            os.fsync(log.fileno())
        returncode = proc.wait()
        if cleanup is None:
            cleanup = cleanup_owned(owned, supervisor_identity, proc)
        elapsed_seconds = time.monotonic() - start
        usage_after = resource.getrusage(resource.RUSAGE_CHILDREN)
        try:
            after_inputs = input_manifest(job, cwd, cmd, jdir)
            input_error = None
        except (OSError, ValueError) as exc:
            after_inputs = {}
            input_error = f"{type(exc).__name__}:{exc}"
        inputs_unchanged = after_inputs == current_inputs and input_error is None
        log_hash = sha256_file(jdir / "job.log")
        outputs_hashes = {}
        missing = []
        for rel in job.get("outputs", []):
            p = jdir / rel
            if p.is_file():
                outputs_hashes[rel] = sha256_file(p)
            else:
                missing.append(rel)
        failed = returncode != 0 or reasons or missing or not cleanup["verified_gone"] or not inputs_unchanged
        failure_reasons = list(reasons)
        if missing:
            failure_reasons.append("missing_outputs:" + ",".join(missing))
        if not cleanup["verified_gone"]:
            failure_reasons.append("cleanup_not_verified")
        if not inputs_unchanged:
            failure_reasons.append("input_hash_changed_during_job" + (":" + input_error if input_error else ""))
        rec.update({"finished_unix": time.time(), "returncode": returncode,
                    "status": "completed" if not failed else ("unclean" if not cleanup["verified_gone"] else "failed"),
                    "failure_reasons": failure_reasons,
                    "log_sha256": log_hash, "outputs_sha256": outputs_hashes,
                    "inputs_sha256_after": after_inputs,
                    "inputs_unchanged": inputs_unchanged,
                    "elapsed_seconds": elapsed_seconds,
                    "user_cpu_seconds_children": max(0.0, usage_after.ru_utime - usage_before.ru_utime),
                    "system_cpu_seconds_children": max(0.0, usage_after.ru_stime - usage_before.ru_stime),
                    "peak_observed_threads": peak_threads,
                    "peak_observed_rss_gib_userspace": round(peak_rss_kb / (1024 * 1024), 4),
                    "cleanup": cleanup,
                    "final_owned_processes": sorted([list(x) for x in owned])})
        atomic_json(state_path, state)
        if failed:
            state["status"] = rec["status"]; state["finished_unix"] = time.time()
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
