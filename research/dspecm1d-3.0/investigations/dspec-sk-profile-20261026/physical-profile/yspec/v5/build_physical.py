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
ROOT = HERE.parents[2]
AUDIT = ROOT / "algorithm-audit"
MATH = ROOT / "math-profile" / "yspec"
PROJECT = Path("/home/adcm2/Documents/Research/Projects/YSpecpp")
CACHE = Path("/home/adcm2/Documents/Research/Projects/DSpecM1D/build/phase-d-build-cache/build-yspecpp")
GATE = AUDIT / "run_compile.py"
BUILD_TAG = "v5"


def output(mode: str, name: str) -> Path:
    path = HERE / "build" / BUILD_TAG / mode / name
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
    app_args = list(audit_meta["compile"])
    app_source = HERE / "physical_yspec.cpp"
    app_obj = output(mode, "physical_yspec.o")
    old_source = app_args[app_args.index("-c") + 1]
    app_args[app_args.index(old_source)] = str(app_source)
    app_args[app_args.index("-o") + 1] = str(app_obj)
    app_args[1:1] = ["-I", str(HERE / "source"), "-I", str(HERE / "source" / "src")]
    if mode == "native":
        app_args.append("-DMATHPROFILE_DISABLED")

    library_commands = []
    library_objects = []
    library_sources = [item for item in cache_meta
                       if Path(item["file"]).parent == PROJECT / "src"]
    expected = {"AngularFunctions.cpp", "Coefficients.cpp", "Input.cpp",
                "InputSections.cpp", "Instrument.cpp", "Kernels.cpp",
                "Location.cpp", "Model.cpp", "RadialSolver.cpp", "Run.cpp",
                "Snapshots.cpp", "Synthetics.cpp", "TimeSeries.cpp",
                "SampledKernels.cpp"}
    actual = {Path(item["file"]).name for item in library_sources}
    if actual != expected:
        raise RuntimeError(f"YSpec library translation-unit set differs: {sorted(actual)}")
    for item in library_sources:
        source = HERE / "source" / "src" / Path(item["file"]).name
        target = output(mode, Path(item["file"]).name + ".o")
        args = replace_source_and_output(shlex.split(item["command"]), item["file"],
                                         str(source), target)
        args[1:1] = ["-I", str(HERE / "source"), "-I", str(HERE / "source" / "src")]
        if mode == "native":
            args.append("-DMATHPROFILE_DISABLED")
        library_commands.append(args)
        library_objects.append(str(target))

    link = list(math_meta["link"])
    old_app = next(i for i, token in enumerate(link) if token.endswith("/yspec_audit.o"))
    link[old_app] = str(app_obj)
    link = [token for token in link
            if not (token.endswith(".cpp.o") or token.endswith("libYSpec.a"))]
    app_index = link.index(str(app_obj))
    link[app_index + 1:app_index + 1] = library_objects
    binary = output(mode, "physical-yspec")
    link[link.index("-o") + 1] = str(binary)
    cwd = math_meta["link_cwd"]
    return [*library_commands, app_args], link, cwd, binary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("native", "instrumented"))
    parser.add_argument("--link-only", action="store_true")
    args = parser.parse_args()
    compile_commands, link, cwd, binary = commands(args.mode)
    commands_to_run = ([] if args.link_only else compile_commands) + [link]
    manifest = {"mode": args.mode, "compile": [] if args.link_only else compile_commands,
                "link": link, "link_cwd": cwd, "binary": str(binary),
                "isolated_library_source": str(HERE / "source" / "src"),
                "archive_linked": False,
                "build_tag": BUILD_TAG}
    manifest["build_tag"] = BUILD_TAG
    destination = HERE / "build" / f"{BUILD_TAG}-{args.mode}-commands.json"
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
