#!/usr/bin/env python3
"""Build isolated YSpec math-profile and disabled-hook control binaries."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shlex
import subprocess
import sys

HERE = Path(__file__).resolve().parent
PROJECT = Path('/home/adcm2/Documents/Research/Projects/YSpecpp')
CACHE = Path('/home/adcm2/Documents/Research/Projects/DSpecM1D/build/phase-d-build-cache/build-yspecpp')
ROOT = HERE.parents[1]
GATE = ROOT / 'algorithm-audit' / 'run_compile.py'
FILES = ['AngularFunctions.cpp', 'RadialSolver.cpp', 'SampledKernels.cpp']

def output_path(mode: str, name: str) -> Path:
    path = HERE / 'build' / 'v2' / mode / name
    path.parent.mkdir(parents=True, exist_ok=True)
    return path

def replace_source_and_output(command: list[str], old: str, new: str, out: Path) -> list[str]:
    args = list(command)
    try:
        source_index = args.index(old)
        args[source_index] = new
    except ValueError:
        raise RuntimeError(f'compile command has no source token: {old}')
    out_index = args.index('-o')
    args[out_index + 1] = str(out)
    return args

def compile_commands(mode: str) -> list[list[str]]:
    disabled = mode == 'native'
    entries = json.loads((CACHE / 'compile_commands.json').read_text())
    by_name = {Path(item['file']).name: item for item in entries}
    commands: list[list[str]] = []
    for name in FILES:
        item = by_name[name]
        args = shlex.split(item['command'])
        src = HERE / 'source' / 'src' / name
        out = output_path(mode, name + '.o')
        args = replace_source_and_output(args, item['file'], str(src), out)
        if disabled:
            args.append('-DMATHPROFILE_DISABLED')
        commands.append(args)

    app_commands = json.loads((ROOT / 'algorithm-audit' / 'yspec_build_commands.json').read_text())
    args = app_commands['compile']
    src = HERE / 'yspec_audit.cpp'
    out = output_path(mode, 'yspec_audit.o')
    args = replace_source_and_output(args, args[args.index('-c') + 1], str(src), out)
    if disabled:
        args.append('-DMATHPROFILE_DISABLED')
    commands.append(args)
    return commands

def link_command(mode: str) -> tuple[list[str], Path, str]:
    metadata = json.loads((ROOT / 'algorithm-audit' / 'yspec_build_commands.json').read_text())
    args = list(metadata['link'])
    app_obj = output_path(mode, 'yspec_audit.o')
    app_index = next(i for i, item in enumerate(args) if item.endswith('/algorithm-audit/yspec_audit.o'))
    args[app_index] = str(app_obj)
    objs = [output_path(mode, name + '.o') for name in FILES]
    args[app_index + 1:app_index + 1] = [str(path) for path in objs]
    binary = output_path(mode, 'yspec-math-profile')
    out_index = args.index('-o')
    args[out_index + 1] = str(binary)
    return args, binary, metadata['link_cwd']

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['native', 'instrumented'])
    parser.add_argument('--link-only', action='store_true')
    args = parser.parse_args()
    mode = args.mode
    commands = [] if args.link_only else compile_commands(mode)
    link, binary, link_cwd = link_command(mode)
    record = {'mode': mode, 'compile': commands, 'link': link, 'link_cwd': link_cwd, 'binary': str(binary)}
    (HERE / 'build' / 'v2' / f'{mode}-commands.json').write_text(json.dumps(record, indent=2) + '\n')
    for command in commands + [link]:
        print('RUN', shlex.join(command), flush=True)
        result = subprocess.run([sys.executable, str(GATE), '--', *command], cwd=link_cwd if command is link else None)
        if result.returncode:
            return result.returncode
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
