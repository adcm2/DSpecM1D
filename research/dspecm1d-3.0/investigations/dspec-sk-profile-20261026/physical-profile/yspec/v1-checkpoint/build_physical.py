#!/usr/bin/env python3
"""Build the isolated YSpec physical adapter through the shared compile gate."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shlex
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
AUDIT = ROOT / "algorithm-audit"
MATH = ROOT / "math-profile" / "yspec"
PROJECT = Path("/home/adcm2/Documents/Research/Projects/YSpecpp")
CACHE = Path("/home/adcm2/Documents/Research/Projects/DSpecM1D/build/phase-d-build-cache/build-yspecpp")
GATE = AUDIT / "run_compile.py"


def output(mode: str, name: str) -> Path:
    path = HERE / "build" / mode / name
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def replace_source_and_output(command: list[str], old: str, new: str, target: Path) -> list[str]:
    args = list(command)
    try:
        args[args.index(old)] = new
    except ValueError as exc:
        raise RuntimeError(f"compile command has no source token: {old}") from exc
    args[args.index("-o") + 1] = str(target)
    return args


def commands(mode: str) -> tuple[list[list[str]], list[str], str, Path]:
    math_meta = json.loads((MATH / "build" / "v2" / f"{mode}-commands.json").read_text())
    audit_meta = json.loads((AUDIT / "yspec_build_commands.json").read_text())
    cache_meta = json.loads((CACHE / "compile_commands.json").read_text())
    coeff = next(item for item in cache_meta if Path(item["file"]).name == "Coefficients.cpp")
    coeff_args = shlex.split(coeff["command"])
    coeff_source = HERE / "source" / "src" / "Coefficients.cpp"
    coeff_obj = output(mode, "Coefficients.cpp.o")
    coeff_args = replace_source_and_output(coeff_args, coeff["file"], str(coeff_source), coeff_obj)
    coeff_args[1:1] = ["-I", str(HERE / "source" / "src")]

    kernels = next(item for item in cache_meta if Path(item["file"]).name == "Kernels.cpp")
    kernels_args = shlex.split(kernels["command"])
    kernels_source = HERE / "source" / "src" / "Kernels.cpp"
    kernels_obj = output(mode, "Kernels.cpp.o")
    kernels_args = replace_source_and_output(kernels_args, kernels["file"], str(kernels_source), kernels_obj)
    kernels_args[1:1] = ["-I", str(HERE / "source" / "src")]

    app_args = list(audit_meta["compile"])
    app_source = HERE / "physical_yspec.cpp"
    app_obj = output(mode, "physical_yspec.o")
    old_source = app_args[app_args.index("-c") + 1]
    app_args[app_args.index(old_source)] = str(app_source)
    app_args[app_args.index("-o") + 1] = str(app_obj)
    app_args[1:1] = ["-I", str(HERE / "source" / "src")]
    if mode == "native":
        app_args.append("-DMATHPROFILE_DISABLED")
        coeff_args.append("-DMATHPROFILE_DISABLED")
        kernels_args.append("-DMATHPROFILE_DISABLED")

    objdir = MATH / "build" / "v2" / mode
    link = list(math_meta["link"])
    old_app = next(i for i, token in enumerate(link) if token.endswith("/yspec_audit.o"))
    link[old_app] = str(app_obj)
    link[old_app + 1:old_app + 1] = [str(coeff_obj), str(kernels_obj)]
    binary = output(mode, "physical-yspec")
    link[link.index("-o") + 1] = str(binary)
    cwd = math_meta["link_cwd"]
    return [coeff_args, kernels_args, app_args], link, cwd, binary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("native", "instrumented"))
    parser.add_argument("--link-only", action="store_true")
    args = parser.parse_args()
    compile_commands, link, cwd, binary = commands(args.mode)
    commands_to_run = ([] if args.link_only else compile_commands) + [link]
    manifest = {"mode": args.mode, "compile": [] if args.link_only else compile_commands,
                "link": link, "link_cwd": cwd, "binary": str(binary),
                "hook_objects_from": str(MATH / "build" / "v2" / args.mode)}
    destination = HERE / "build" / f"{args.mode}-commands.json"
    destination.write_text(json.dumps(manifest, indent=2) + "\n")
    for command in commands_to_run:
        print("RUN", shlex.join(command), flush=True)
        result = subprocess.run([sys.executable, str(GATE), "--", *command], cwd=cwd if command is link else None)
        if result.returncode:
            return result.returncode
    print(f"BUILT {binary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
