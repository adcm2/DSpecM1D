#!/usr/bin/env python3
"""Prepare exact sequential compiler/linker commands without executing them."""
import json
from pathlib import Path
import shlex

HERE = Path(__file__).resolve().parent
PROJECT = Path('/home/adcm2/Documents/Research/Projects/YSpecpp')
CACHE = Path('/home/adcm2/Documents/Research/Projects/DSpecM1D/build/phase-d-build-cache/build-yspecpp')
ROOT = HERE.parents[3]
FILES = ['AngularFunctions.cpp', 'RadialSolver.cpp', 'SampledKernels.cpp']
OUT = HERE / 'build'
OUT.mkdir(parents=True, exist_ok=True)


def replace_source_and_output(command, old, new, output):
    args = list(command)
    args[args.index(old)] = str(new)
    args[args.index('-o') + 1] = str(output)
    return args


entries = json.loads((CACHE / 'compile_commands.json').read_text())
by_name = {Path(item['file']).name: item for item in entries}
commands = []
for name in FILES:
    item = by_name[name]
    source = HERE / 'source' / 'src' / 'cpp' / name
    obj = OUT / (name + '.o')
    commands.append(replace_source_and_output(shlex.split(item['command']), item['file'], source, obj))

app_metadata = json.loads((ROOT / 'algorithm-audit' / 'yspec_build_commands.json').read_text())
app_source = HERE / 'algorithm-audit' / 'yspec_audit.cpp'
app_obj = OUT / 'yspec_audit.o'
app_compile = app_metadata['compile']
original_app_source = app_compile[app_compile.index('-c') + 1]
commands.append(replace_source_and_output(app_compile, original_app_source, app_source, app_obj))

link_metadata = json.loads((ROOT / 'algorithm-audit' / 'yspec_build_commands.json').read_text())
link = list(link_metadata['link'])
app_index = next(i for i, item in enumerate(link) if item.endswith('/algorithm-audit/yspec_audit.o'))
link[app_index] = str(app_obj)
objects = [OUT / (name + '.o') for name in FILES]
link[app_index + 1:app_index + 1] = [str(p) for p in objects]
binary = OUT / 'yspec-math-profile-coarse'
link[link.index('-o') + 1] = str(binary)
record = {
    'variant': 'coarse-profiler',
    'description': 'Command recipe; execution status is recorded separately in implementation_status.md and parent-run records.',
    'compile': commands,
    'link': link,
    'link_cwd': link_metadata['link_cwd'],
    'binary': str(binary),
}
(HERE / 'build-commands.json').write_text(json.dumps(record, indent=2) + '\n')
print(HERE / 'build-commands.json')
