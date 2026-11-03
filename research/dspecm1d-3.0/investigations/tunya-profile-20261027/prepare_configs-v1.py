"""Generate explicit Tunya commands; no build or numerical execution."""
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
C = '/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya'
I = C + '/continuation-inputs/v1'
TASK = C + '/continuation/profile-01'
P = C + '/bootstrap/original/dspec-sk-profile-20261026/physical-profile'
D = C + '/bootstrap/dependencies'
B = '/space/adcm2/research/builds/kernel-profile-20261027'
DS = '/space/adcm2/research/projects/DSpecM1D'
YS = '/space/adcm2/research/projects/YSpecpp'
E = '/space/adcm2/CCDSM/deps/toolchain/eigen/include/eigen3'
PY = '/space/adcm2/CCDSM/deps/toolchain/venv/bin/python3'
CXX = '/opt/gcc-13.2.0/bin/g++'
env = {'LD_LIBRARY_PATH': '/opt/gcc-13.2.0/lib64', 'LD_PRELOAD': '',
       'DSPEC_PROFILE_FTZ': '0', 'PYTHONDONTWRITEBYTECODE': '1'}
common = [CXX, '-O3', '-DNDEBUG', '-std=gnu++23', '-fopenmp', '-DGSHTRANS_HAVE_INTERPOLATION']
incs = [D+'/numericconcepts/include', D+'/fftwpp', D+'/gaussquad',
        D+'/interpolation/include', D+'/gshtrans']
variants = []
for tag, source in [('dsv4', P+'/dspec/source-v4-frozen'),
                    ('dsv6', I+'/dspec/source-v6-frozen')]:
    for mode in ['native', 'profiled']:
        includes = [DS, *incs, D+'/planetmodel-dspec/cpp/include']
        argv = common + ['-DDSFAST_DISABLE_PROFILE'] + ['-I'+p for p in includes] + ['-isystem', E]
        if mode == 'native':
            argv += ['-DMATHPROFILE_DISABLED']
        argv += [source+'/physical_driver.cpp', '-o', '{OUTPUT}', '-lfftw3', '-lfftw3f', '-lfftw3l', '-llapack', '-lblas']
        variants.append({'id': tag+'-'+mode, 'cwd': source, 'output_name': 'physical-driver',
                         'command': argv, 'sources': [source],
                         'dependencies': [DS+'/DSpecM1D', *incs, D+'/planetmodel-dspec/cpp/include', E]})
source = P+'/yspec/v6'
tus = ['AngularFunctions', 'Coefficients', 'Input', 'InputSections', 'Instrument',
       'Kernels', 'Location', 'Model', 'RadialSolver', 'Run', 'Snapshots',
       'Synthetics', 'TimeSeries', 'SampledKernels']
for mode in ['native', 'profiled']:
    includes = [source+'/source', source+'/source/src', YS, YS+'/src', *incs,
                D+'/planetmodel-yspec/cpp/include']
    argv = common + ['-DPLANETMODEL_VERSION="1.2.5"', '-DPLANETMODEL_WITH_NETCDF=1',
                     '-DYSPEC_VERSION="0.1.0"', '-fPIC', '-Wall', '-Wextra']
    argv += ['-I'+p for p in includes] + ['-isystem', E]
    if mode == 'native':
        argv += ['-DMATHPROFILE_DISABLED']
    argv += [source+'/source/src/'+tu+'.cpp' for tu in tus]
    argv += [source+'/physical_yspec.cpp', '-o', '{OUTPUT}', '-lnetcdf', '-lfftw3', '-lfftw3f', '-lfftw3l', '-ldl', '-lm']
    variants.append({'id': 'ysv6-'+mode, 'cwd': source, 'output_name': 'physical-yspec',
                     'command': argv, 'sources': [source+'/source', source+'/physical_yspec.cpp', source+'/profile.hpp'],
                     'dependencies': [YS+'/YSpec', YS+'/src', *incs, D+'/planetmodel-yspec/cpp/include', E]})
gate = I+'/dspec/validation/receiver_gate.cpp'
variants.append({'id': 'receiver-gate', 'cwd': I+'/dspec/validation', 'output_name': 'receiver-gate',
                 'command': common+['-I'+D+'/gshtrans', '-I'+D+'/numericconcepts/include', '-isystem', E,
                                    gate, '-o', '{OUTPUT}'],
                 'sources': [gate], 'dependencies': [D+'/gshtrans', D+'/numericconcepts/include', E]})
out = HERE/'configs'
out.mkdir(exist_ok=True)
def save(name, value):
    path = out/name
    if path.exists():
        raise RuntimeError('refusing to replace config: '+str(path))
    path.write_text(json.dumps(value, indent=2)+'\n')
save('build.json', {'build_root': B, 'environment': env, 'variants': variants})
save('build-plan.json', {'campaign_root': C, 'thread_cap': 50, 'environment': env,
    'jobs': [{'id': 'build', 'cwd': I,
              'command': ['/usr/bin/python3', I+'/build_receipt.py', I+'/tools/build_campaign.py', I+'/configs/build.json', '{JOB_DIR}/build-result.json'],
              'inputs': [I+'/build_receipt.py', I+'/tools/build_campaign.py', I+'/configs/build.json', CXX],
              'outputs': ['build-result.json'], 'timeout_seconds': 14400, 'memory_gib': 24}]})

jobs = []
def add(jid, command, inputs, outputs, timeout=7200, memory=24):
    inputs = [p.replace('{TASK_DIR}', TASK) for p in inputs]
    jobs.append({'id': jid, 'cwd': I, 'command': command, 'inputs': inputs,
                 'outputs': outputs, 'timeout_seconds': timeout, 'memory_gib': memory})
pin = ['/usr/bin/taskset', '-c', '0']
model = DS+'/data/models/PREM500'
gatebin = B+'/receiver-gate/receiver-gate'
add('receiver-gate', pin+[gatebin], [gatebin], [], 3600)
def ds(tag, mode, degree, probe=False):
    jid = tag+'-'+mode+'-l'+str(degree)
    binary = B+'/'+tag+'-'+mode+'/physical-driver'
    command = pin+[binary, model, '{JOB_DIR}/data', 'kernel-probe' if probe else 'kernel',
                   '--lmax', str(degree), '--nq', '6', '--radial-refinement', '0.05',
                   '--radial-traces', '0', '--dc', 'include', '--contacts', 'off',
                   '--pole-only', 'on', '--forward-only', 'on', '--max-output-gib', '8']
    if probe:
        command += ['--probe-angle-count', '64']
    outputs = ['data/result.json', 'data/trace.csv', 'data/spectrum.csv', 'data/radial_nodes.csv',
               'data/kernel_radial_nodes.csv', 'data/kernel_probe.bin' if probe else 'data/kernel.bin',
               'data/kernel_probe.json' if probe else 'data/kernel_domain.json']
    add(jid, command, [binary, model], outputs, 14400 if degree>384 else 7200)
    return '{TASK_DIR}/'+jid+'/data'
def ys(mode, degree, baseline, waveform=False, probe=False):
    binary = B+'/ysv6-'+mode+'/physical-yspec'
    jid = 'ysv6-'+mode+('-wave-l' if waveform else '-l')+str(degree)
    if waveform:
        args = ['waveform', model, '{JOB_DIR}/data', str(degree)]
        inputs = [binary, model]
        outputs = ['data/result.json', 'data/trace.csv', 'data/spectrum.csv']
    else:
        args = ['kernel-probe' if probe else 'kernel', model, baseline+'/trace.csv',
                baseline+'/radial_nodes.csv', '{JOB_DIR}/data', str(degree)]
        inputs = [binary, model, baseline+'/trace.csv', baseline+'/radial_nodes.csv']
        outputs = ['data/result.json', 'data/kernel_probe.bin' if probe else 'data/kernels.bin']
    args += ['1e-9', '0.125', '--include-dc', '--fft-padding', '1']
    if not waveform:
        args += ['--nodes-per-element', '6', '--kernel-memory-budget-mib', '8192', '--kernel-disk-budget-mib', '8192']
    if probe:
        args += ['--probe-angle-count', '64']
    add(jid, pin+[binary]+args, inputs, outputs, 14400 if degree>384 else 7200)
for tag in ['dsv4', 'dsv6']:
    for mode in ['native', 'profiled']:
        ds(tag, mode, 3)
baseline = '{TASK_DIR}/dsv6-native-l3/data'
for mode in ['native', 'profiled']:
    ys(mode, 3, baseline)
check = I+'/tools/qualification_checks.py'
check_inputs = [check]
for prior in jobs:
    check_inputs += [TASK+'/'+prior['id']+'/'+rel for rel in prior['outputs']]
add('qualification-checks', [PY, check, '{TASK_DIR}', '{JOB_DIR}/qualification.json'],
    check_inputs, ['qualification.json'], 600)
for degree in [128, 384]:
    for tag in ['dsv4', 'dsv6']:
        for mode in ['native', 'profiled']:
            ds(tag, mode, degree, True)
    base = '{TASK_DIR}/dsv6-native-l'+str(degree)+'/data'
    for mode in ['native', 'profiled']:
        ys(mode, degree, base, waveform=True)
        ys(mode, degree, base, probe=True)
for mode in ['native', 'profiled']:
    ds('dsv6', mode, 1536, True)
base = '{TASK_DIR}/dsv6-native-l1536/data'
for mode in ['native', 'profiled']:
    ys(mode, 1536, base, waveform=True)
    ys(mode, 1536, base, probe=True)
save('qualification-and-profile-plan.json', {'campaign_root': C, 'thread_cap': 50,
    'environment': env, 'jobs': jobs,
    'measurement_note': 'One numerical thread, serial jobs, CPU 0. Degree>3 kernels use 64-angle diagnostics only, never a full global kernel claim. YSpec forward waveform and kernel timings remain separate.'})
print('Wrote', len(variants), 'build variants and', len(jobs), 'qualification/profile jobs.')
