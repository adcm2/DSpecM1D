#!/usr/bin/env python3
"""Run the two bounded Stage 1 PREM proxy cases outside the source tree."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import time
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("executable", type=Path)
    parser.add_argument("run_dir", type=Path)
    args = parser.parse_args()
    source_root = Path(__file__).resolve().parents[2]
    run_dir = args.run_dir.resolve()
    run_dir.mkdir(parents=True, exist_ok=False)
    input_dir = run_dir / "inputs"
    output_dir = run_dir / "outputs"
    log_dir = run_dir / "logs"
    for directory in (input_dir, output_dir, log_dir):
        directory.mkdir()
    model = source_root / "data" / "models" / "prem.200.isotropic"
    cases = {"toroidal": source_root / "validation" / "diffracted_wave" / "toroidal.params",
             "spheroidal": source_root / "validation" / "diffracted_wave" / "spheroidal.params"}
    results: dict[str, object] = {
        "schema": "dspecm1d-diffracted-wave-stage1-v1",
        "source_root": str(source_root),
        "source_head": subprocess.check_output(["git", "-C", str(source_root), "rev-parse", "HEAD"], text=True).strip(),
        "model": str(model),
        "model_sha256": sha256(model),
        "executable": str(args.executable.resolve()),
        "executable_sha256": sha256(args.executable),
        "omp_num_threads": 1,
        "cases": {},
        "interpretation": "finite/nonzero smoke only; unconverged proxy, not a scientific validation",
    }
    for name, template in cases.items():
        output_prefix = output_dir / name
        output = output_dir / f"{name}_t.out"
        materialized = input_dir / f"{name}.params"
        text = template.read_text(encoding="utf-8").replace("@MODEL@", str(model)).replace("@OUTPUT@", str(output_prefix))
        materialized.write_text(text, encoding="utf-8")
        log_path = log_dir / f"{name}.log"
        started = time.monotonic()
        env = os.environ.copy()
        env.update({"OMP_NUM_THREADS": "1", "OMP_THREAD_LIMIT": "1",
                    "OMP_DYNAMIC": "FALSE", "OPENBLAS_NUM_THREADS": "1",
                    "MKL_NUM_THREADS": "1"})
        resource_path = log_dir / f"{name}.resource.txt"
        with log_path.open("w", encoding="utf-8") as log:
            subprocess.run(["/usr/bin/busybox", "time", "-v", "-o", str(resource_path),
                            "/usr/bin/timeout", "--kill-after=5s", "300s",
                            str(args.executable.resolve()), str(materialized)],
                           cwd=args.executable.parent.parent, env=env, stdout=log,
                           stderr=subprocess.STDOUT, check=True, timeout=315)
        elapsed = time.monotonic() - started
        resource_lines = resource_path.read_text(encoding="utf-8").splitlines()
        max_rss_kb = next(int(line.split(":", 1)[1].strip()) for line in resource_lines
                          if line.startswith("\tMaximum resident set size (kbytes):"))
        checker = source_root / "validation" / "diffracted_wave" / "check_outputs.py"
        check = subprocess.run(["python3", str(checker), str(output)], text=True,
                               capture_output=True, check=True)
        results["cases"][name] = {
            "template": str(template),
            "template_sha256": sha256(template),
            "input": str(materialized),
            "input_sha256": sha256(materialized),
            "output": str(output),
            "output_sha256": sha256(output),
            "output_bytes": output.stat().st_size,
            "log": str(log_path),
            "resource_receipt": str(resource_path),
            "resource_command": "/usr/bin/busybox time -v -o <receipt> /usr/bin/timeout --kill-after=5s 300s",
            "elapsed_s": elapsed,
            "max_rss_kb": max_rss_kb,
            "checker_stdout": check.stdout.strip(),
        }
    summary = run_dir / "summary.json"
    summary.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(results, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
