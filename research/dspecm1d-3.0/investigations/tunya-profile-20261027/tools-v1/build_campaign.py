#!/usr/bin/env python3
"""Build explicit frozen-adapter variants from a JSON manifest.

This tool never edits inputs. It refuses a populated build root and writes a
command/provenance record beside each binary. Run only after independent review.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def tree_hash(path: Path) -> str:
    if path.is_file():
        return sha256_file(path)
    if not path.is_dir():
        raise ValueError(f"input path does not exist: {path}")
    h = hashlib.sha256()
    for item in sorted(path.rglob("*")):
        rel = item.relative_to(path).as_posix().encode()
        h.update(len(rel).to_bytes(8, "big")); h.update(rel)
        if item.is_symlink():
            target = os.readlink(item).encode()
            h.update(b"L" + len(target).to_bytes(8, "big") + target)
        elif item.is_file():
            h.update(b"F" + bytes.fromhex(sha256_file(item)))
        elif item.is_dir():
            h.update(b"D")
    return h.hexdigest()


def replace_tokens(argv: list[str], values: dict[str, str]) -> list[str]:
    out = []
    for token in argv:
        for key, value in values.items():
            token = token.replace("{" + key + "}", value)
        if "{" in token or "}" in token:
            raise ValueError(f"unresolved command token: {token}")
        out.append(token)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("config", type=Path)
    args = ap.parse_args()
    cfg = json.loads(args.config.read_text())
    root = Path(cfg["build_root"])
    if not root.is_absolute():
        raise ValueError("build_root must be absolute")
    if root.exists() and any(root.iterdir()):
        raise ValueError(f"build_root must be fresh and empty: {root}")
    root.mkdir(parents=True, exist_ok=True)
    jobs = int(cfg.get("compiler_jobs", 1))
    if not 1 <= jobs <= 2:
        raise ValueError("compiler_jobs must be 1 or 2")
    variants = cfg["variants"]
    if not variants or len({v["id"] for v in variants}) != len(variants):
        raise ValueError("variants must be nonempty with unique ids")
    env = os.environ.copy()
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
                "GOTO_NUM_THREADS", "OMP_THREAD_LIMIT"):
        env[key] = "1"
    env.update({str(k): str(v) for k, v in cfg.get("environment", {}).items()})
    # Numerical thread controls are campaign invariants, even if config supplies them.
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
                "GOTO_NUM_THREADS", "OMP_THREAD_LIMIT"):
        env[key] = "1"

    summary = {"schema": 1, "started_unix": time.time(), "build_root": str(root),
               "compiler_jobs": jobs, "environment_thread_limits": "all=1", "variants": []}
    (root / "build_campaign.json").write_text(json.dumps(summary, indent=2) + "\n")
    for v in variants:
        vid = v["id"]
        if not vid.replace("-", "").replace("_", "").isalnum():
            raise ValueError(f"variant id may contain only letters, digits, hyphens and underscores: {vid}")
        vdir = root / vid
        if vdir.exists():
            raise ValueError(f"variant output already exists: {vdir}")
        vdir.mkdir()
        output_name = v.get("output_name", "program")
        if Path(output_name).name != output_name or output_name in {".", ".."}:
            raise ValueError(f"output_name must be a simple filename: {output_name}")
        output = vdir / output_name
        values = {"BUILD_DIR": str(vdir), "OUTPUT": str(output), "JOBS": str(jobs)}
        command = replace_tokens(v["command"], values)
        cwd = Path(v["cwd"])
        if not cwd.is_absolute() or not cwd.is_dir():
            raise ValueError(f"variant cwd must be an existing absolute directory: {cwd}")
        tracked = [Path(p) for p in v.get("sources", []) + v.get("dependencies", [])]
        if any(not p.is_absolute() for p in tracked):
            raise ValueError(f"all source and dependency paths must be absolute ({vid})")
        fingerprints = {str(p): tree_hash(p) for p in tracked}
        record = {"id": vid, "cwd": str(cwd), "command": command,
                  "command_display": shlex.join(command), "command_sha256": hashlib.sha256(
                      json.dumps(command, separators=(",", ":")).encode()).hexdigest(),
                  "tracked_input_sha256": fingerprints, "started_unix": time.time(),
                  "output": str(output)}
        (vdir / "command.json").write_text(json.dumps(record, indent=2) + "\n")
        with (vdir / "stdout.log").open("wb") as stdout, (vdir / "stderr.log").open("wb") as stderr:
            proc = subprocess.run(command, cwd=cwd, env=env, stdout=stdout, stderr=stderr, check=False)
        after_fingerprints = {str(p): tree_hash(p) for p in tracked}
        inputs_changed = after_fingerprints != fingerprints
        record.update({"returncode": proc.returncode, "finished_unix": time.time(),
                       "stdout_sha256": sha256_file(vdir / "stdout.log"),
                       "stderr_sha256": sha256_file(vdir / "stderr.log"),
                       "binary_sha256": sha256_file(output) if output.is_file() else None,
                       "tracked_input_sha256_after": after_fingerprints,
                       "tracked_inputs_unchanged": not inputs_changed})
        (vdir / "command.json").write_text(json.dumps(record, indent=2) + "\n")
        summary["variants"].append({"id": vid, "returncode": proc.returncode,
                                    "command_sha256": record["command_sha256"],
                                    "binary_sha256": record["binary_sha256"]})
        summary["updated_unix"] = time.time()
        (root / "build_campaign.json").write_text(json.dumps(summary, indent=2) + "\n")
        if proc.returncode != 0 or inputs_changed:
            summary["status"] = "failed"
            (root / "build_campaign.json").write_text(json.dumps(summary, indent=2) + "\n")
            return proc.returncode or 1
        if not output.is_file():
            summary["status"] = "failed_missing_output"
            (root / "build_campaign.json").write_text(json.dumps(summary, indent=2) + "\n")
            return 1
    summary["status"] = "completed"
    summary["finished_unix"] = time.time()
    (root / "build_campaign.json").write_text(json.dumps(summary, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"build_campaign: {exc}", file=sys.stderr)
        raise SystemExit(2)
