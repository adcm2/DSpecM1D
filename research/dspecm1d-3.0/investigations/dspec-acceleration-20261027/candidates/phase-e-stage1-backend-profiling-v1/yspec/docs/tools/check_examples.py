#!/usr/bin/env python3
"""Run the small downloadable manual fixtures and compile its library example."""
import argparse
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def run(command, env=None):
    print('+', ' '.join(map(str, command)), flush=True)
    subprocess.run(command, check=True, cwd=ROOT, env=env)


def finite_rows(path, columns, skip=0):
    lines = path.read_text().splitlines()[skip:]
    if not lines:
        raise RuntimeError(f'{path}: no data rows')
    for number, line in enumerate(lines, skip + 1):
        values = line.split()
        if len(values) != columns:
            raise RuntimeError(f'{path}:{number}: expected {columns} columns, got {len(values)}')
        if not all(map(lambda value: re.fullmatch(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?', value), values)):
            raise RuntimeError(f'{path}:{number}: expected finite numeric values')
        if not all(math.isfinite(float(value)) for value in values):
            raise RuntimeError(f'{path}:{number}: non-finite value')
    return lines


def numeric_rows(path, columns=None, tagged=False):
    rows = []
    for number, line in enumerate(path.read_text().splitlines(), 1):
        if not line or line.startswith('#'):
            continue
        fields = line.split()
        if tagged:
            if fields[0] not in ('path', 'across'):
                raise RuntimeError(f'{path}:{number}: expected path/across section tag')
            fields = fields[1:]
        if columns is not None and len(fields) != columns:
            raise RuntimeError(f'{path}:{number}: expected {columns} numeric columns, got {len(fields)}')
        try:
            values = [float(value) for value in fields]
        except ValueError as error:
            raise RuntimeError(f'{path}:{number}: expected numeric data') from error
        if not all(map(math.isfinite, values)):
            raise RuntimeError(f'{path}:{number}: non-finite value')
        rows.append(values)
    if not rows:
        raise RuntimeError(f'{path}: no data rows')
    return rows


def copy_kernel_fixture(directory, fixture):
    directory.mkdir()
    source_path = ROOT / f'docs/site/examples/{fixture}'
    shutil.copy2(source_path, directory / 'input.json')
    shutil.copy2(ROOT / 'examples/prem.200', directory / 'prem.200')
    source = json.loads(source_path.read_text())['sources'][0]
    return float(source['latitude']), float(source['longitude'])


def check_sampled_grid_contract(section_file, map_file, source_location,
                                probe_truncation=False):
    expected_sampled_parameters = ['vp', 'vs', 'qkappa', 'qmu']
    section_lines = section_file.read_text().splitlines()
    map_lines = map_file.read_text().splitlines()
    expected_unit_line = '# per unit volume of ln X, in s m^-3'
    if len(section_lines) < 3 or section_lines[1] != expected_unit_line:
        raise RuntimeError(f'{section_file}: expected relative sampled-kernel units')
    if len(map_lines) < 3 or map_lines[1] != expected_unit_line:
        raise RuntimeError(f'{map_file}: expected relative sampled-kernel units')

    section_header = next((line.split() for line in section_lines
                           if line.startswith('# section layer depth_km')), None)
    if not section_header:
        raise RuntimeError(f'{section_file}: missing path/across column header')
    sampled_parameters = len(section_header) - 8  # '# section' + six coordinates
    if section_header[8:] != expected_sampled_parameters:
        raise RuntimeError(f'{section_file}: expected sampled parameters '
                           f'{expected_sampled_parameters}, got {section_header[8:]}')
    sampled_parameters = len(expected_sampled_parameters)

    def validate_sections(path):
        lines = path.read_text().splitlines()
        rows = numeric_rows(path, columns=6 + sampled_parameters, tagged=True)
        tags = [line.split()[0] for line in lines
                if line and not line.startswith('#')]
        if len(rows) != 1260:
            raise RuntimeError(f'{path}: expected 1260 sampled section rows, got {len(rows)}')
        if len(tags) != 1260 or tags.count('path') != 1085 or tags.count('across') != 175:
            raise RuntimeError(f'{path}: expected 1085 path and 175 across rows')
        if not any(value != 0 for row in rows for value in row[6:]):
            raise RuntimeError(f'{path}: expected nonzero sampled kernels')

        def check_nodes(tag, coordinates, expected_node_count):
            selected = [row for row, row_tag in zip(rows, tags) if row_tag == tag]
            if len(selected) != len(coordinates) * expected_node_count:
                raise RuntimeError(f'{path}: {tag} grid has an unexpected shape')
            node_keys = []
            for start in range(0, len(selected), len(coordinates)):
                block = selected[start:start + len(coordinates)]
                key = (block[0][0], block[0][1])  # layer, depth_km
                if any((row[0], row[1]) != key for row in block):
                    raise RuntimeError(f'{path}: {tag} rows are not grouped by radial node')
                actual = [(row[2], row[3]) for row in block]
                if any(not (math.isclose(a[0], e[0], abs_tol=1e-4) and
                            math.isclose(a[1], e[1], abs_tol=1e-4))
                       for a, e in zip(actual, coordinates)):
                    raise RuntimeError(f'{path}: {tag} grid coordinates/order do not match fixture')
                node_keys.append(key)
            if len(set(node_keys)) != expected_node_count:
                raise RuntimeError(f'{path}: {tag} radial node sequence is incomplete or duplicated')
            for layer, depth in node_keys:
                if not float(layer).is_integer() or not 0 <= layer <= 12:
                    raise RuntimeError(f'{path}: {tag} contains an invalid model layer')
                if not 0 <= depth <= 6371:
                    raise RuntimeError(f'{path}: {tag} contains a depth outside PREM')
            return selected

        path_coordinates = [(float(along), 0.0) for along in range(0, 61, 2)]
        across_coordinates = [(30.0, float(across)) for across in range(-4, 5, 2)]
        path_rows = check_nodes('path', path_coordinates, 35)
        across_rows = check_nodes('across', across_coordinates, 35)
        source_lat, source_lon = source_location
        for row in path_rows:
            if row[2] == 0.0:
                if not math.isclose(row[4], source_lat, rel_tol=0, abs_tol=5.1e-6):
                    raise RuntimeError(f'{path}: path origin latitude differs from fixture source')
                if not math.isclose(row[5], source_lon, rel_tol=0, abs_tol=5.1e-6):
                    raise RuntimeError(f'{path}: path origin longitude differs from fixture source')
        for row in rows:
            if not -90 <= row[4] <= 90 or not 0 <= row[5] < 360:
                raise RuntimeError(f'{path}: latitude/longitude outside physical ranges')
        return path_rows, across_rows

    def validate_maps(path):
        header = next((line.split() for line in map_lines
                       if line.startswith('# layer depth_km')), None)
        if not header or header[7:] != expected_sampled_parameters:
            raise RuntimeError(f'{path}: expected sampled parameters '
                               f'{expected_sampled_parameters}')
        rows = numeric_rows(path, columns=6 + sampled_parameters)
        if len(rows) != 310:
            raise RuntimeError(f'{path}: expected 310 sampled map rows, got {len(rows)}')
        expected_coordinates = [(float(along), float(across))
                               for along in range(0, 61, 2)
                               for across in range(-4, 5, 2)]
        for start, depth in ((0, 100.0), (155, 500.0)):
            block = rows[start:start + 155]
            if len(block) != 155 or any(row[1] != depth for row in block):
                raise RuntimeError(f'{path}: expected 155 rows at map depth {depth:g} km')
            if any(not float(row[0]).is_integer() or not 0 <= row[0] <= 12
                   for row in block):
                raise RuntimeError(f'{path}: map contains an invalid model layer')
            actual = [(row[2], row[3]) for row in block]
            if any(not (math.isclose(a[0], e[0], abs_tol=1e-4) and
                        math.isclose(a[1], e[1], abs_tol=1e-4))
                   for a, e in zip(actual, expected_coordinates)):
                raise RuntimeError(f'{path}: map grid coordinates/order do not match fixture')
        if not any(row[6] != 0 for row in rows):
            raise RuntimeError(f'{path}: expected nonzero sampled kernels')
        for row in rows:
            if not 0 <= row[1] <= 6371 or not -90 <= row[4] <= 90 or not 0 <= row[5] < 360:
                raise RuntimeError(f'{path}: invalid map depth or latitude/longitude')

    validate_sections(section_file)
    validate_maps(map_file)
    if probe_truncation:
        probe = section_file.with_name(section_file.name + '.contract-probe')

        def expect_rejection(label, content, check, message):
            probe.write_text(content)
            try:
                try:
                    check(probe)
                except RuntimeError as error:
                    if message not in str(error):
                        raise RuntimeError(f'{label} probe failed for an unexpected reason: {error}') from error
                else:
                    raise RuntimeError(f'{label} probe was not rejected')
            finally:
                probe.unlink(missing_ok=True)

        expect_rejection(
            'two-row truncation', '\n'.join(section_lines[:-2]) + '\n',
            validate_sections, 'expected 1260 sampled section rows, got 1258')
        bad_tag_lines = section_lines.copy()
        for index in range(len(bad_tag_lines) - 1, -1, -1):
            if bad_tag_lines[index] and not bad_tag_lines[index].startswith('#'):
                fields = bad_tag_lines[index].split()
                fields[0] = 'other'
                bad_tag_lines[index] = ' '.join(fields)
                break
        expect_rejection(
            'section-tag', '\n'.join(bad_tag_lines) + '\n',
            lambda path: check_sampled_grid_contract(path, map_file, source_location),
            'expected path/across section tag')
        bad_unit_lines = section_lines.copy()
        bad_unit_lines[1] = '# per unit volume of X, in s m^-3'
        expect_rejection(
            'sampled-units', '\n'.join(bad_unit_lines) + '\n',
            lambda path: check_sampled_grid_contract(path, map_file, source_location),
            'expected relative sampled-kernel units')
        bad_coordinate_lines = section_lines.copy()
        for index, line in enumerate(bad_coordinate_lines):
            if line and not line.startswith('#'):
                fields = line.split()
                fields[3] = '1.0000'  # first path point's along coordinate
                bad_coordinate_lines[index] = ' '.join(fields)
                break
        expect_rejection(
            'sampled-coordinate', '\n'.join(bad_coordinate_lines) + '\n',
            lambda path: check_sampled_grid_contract(path, map_file, source_location),
            'grid coordinates/order do not match fixture')


def check_kernel_outputs(directory, prefix, source_location, probe_truncation=False):
    measurement = directory / f'{prefix}.measurement'
    data = numeric_rows(measurement, 4)
    if len(data) != 1201 or not any(row[1] != 0 for row in data):
        raise RuntimeError(f'{measurement}: expected 1,201 finite samples and a nonzero seismogram')
    if not any(row[2] != 0 for row in data) or not any(row[3] != 0 for row in data):
        raise RuntimeError(f'{measurement}: expected a nonzero window and adjoint source')

    radial = directory / f'{prefix}.kernels'
    headers = [line for line in radial.read_text().splitlines() if line.startswith('#')]
    at_source = next((line.split() for line in headers
                      if line.startswith('# at_source ')), None)
    if not at_source:
        raise RuntimeError(f'{radial}: missing source-point P_X terms')
    columns = next((line.split() for line in reversed(headers)
                    if line.startswith('# radius_km depth_km weight_m3')), None)
    if not columns or len(columns) < 4:
        raise RuntimeError(f'{radial}: missing radial kernel parameter columns')
    expected_radial_parameters = ['rho', 'vp', 'vs', 'qkappa', 'qmu']
    if columns[4:] != expected_radial_parameters:
        raise RuntimeError(f'{radial}: expected radial parameters '
                           f'{expected_radial_parameters}, got {columns[4:]}')
    parameter_count = len(columns) - 4
    try:
        point_terms = [float(value) for value in at_source[2:]]
    except ValueError as error:
        raise RuntimeError(f'{radial}: invalid source-point values') from error
    if len(point_terms) != parameter_count or not all(map(math.isfinite, point_terms)):
        raise RuntimeError(f'{radial}: source-point P_X count/values do not match the parameter columns')
    values = numeric_rows(radial, len(columns) - 1)
    if not any(value != 0 for row in values for value in row[3:]):
        raise RuntimeError(f'{radial}: expected nonzero radial kernel values')

    section_file = directory / f'{prefix}.sections'
    map_file = directory / f'{prefix}.maps'
    check_sampled_grid_contract(section_file, map_file, source_location,
                                probe_truncation=probe_truncation)


def check_png(path):
    if not path.is_file() or path.stat().st_size < 1000 or path.read_bytes()[:8] != b'\x89PNG\r\n\x1a\n':
        raise RuntimeError(f'{path}: missing or invalid PNG plot')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, required=True, help='configured Release/OpenMP CMake build directory')
    parser.add_argument('--python', type=Path, required=True, help='Python executable with this checkout installed using the plot extra')
    args = parser.parse_args()
    build = args.build.resolve()
    python = args.python.absolute()  # Preserve the venv executable symlink.
    for binary in ('bin/yspec', 'bin/yspec-snapshot', 'bin/yspec-kernels', 'libYSpec.a'):
        if not (build / binary).is_file():
            parser.error(f'missing {build / binary}; build yspec, yspec-snapshot and YSpec first')
    cache = (build / 'CMakeCache.txt').read_text()
    compiler = re.search(r'^CMAKE_CXX_COMPILER:(?:FILEPATH|STRING|UNINITIALIZED)=(.+)$', cache, re.M)
    openmp = re.search(r'^YSPEC_WITH_OPENMP:BOOL=(ON|OFF)$', cache, re.M)
    if not compiler or not openmp or openmp.group(1) != 'ON':
        parser.error('this focused example check needs the Release/OpenMP build used by docs CI')
    env = os.environ | {'OMP_NUM_THREADS': '2'}
    with tempfile.TemporaryDirectory(prefix='yspec-doc-examples-', dir=build) as temporary:
        work = Path(temporary)
        fixtures = {
            'snapshot': 'snapshot-small',
            'workflows': 'workflows',
            'first-seismogram': 'first-seismogram',
        }
        for name, fixture in fixtures.items():
            (work / name).mkdir()
            shutil.copy2(ROOT / f'docs/site/examples/{fixture}.json', work / name / 'input.json')
            shutil.copy2(ROOT / 'examples/prem.200', work / name / 'prem.200')
        shutil.copy2(ROOT / 'docs/site/examples/station.pz', work / 'workflows/station.pz')
        run([build / 'bin/yspec-snapshot', work / 'snapshot/input.json'], env)
        snapshot = work / 'snapshot/snapshot.1'
        rows = finite_rows(snapshot, 8, skip=2)
        if len(rows) != 4:
            raise RuntimeError(f'{snapshot}: expected four grid rows, got {len(rows)}')
        first_run = work / 'first-seismogram'
        run([build / 'bin/yspec', first_run / 'input.json'], env)
        first_trace = first_run / 'first-seismogram.out.1'
        # The traces are the files with a numbered suffix; beside them yspec
        # writes first-seismogram.out.receivers, a table of each trace's
        # receiver, checked below.
        trace_files = sorted(path for path in first_run.glob('first-seismogram.out.*')
                             if path.suffix[1:].isdigit())
        if trace_files != [first_trace] or not first_trace.is_file():
            raise RuntimeError('first-seismogram run did not create exactly first-seismogram.out.1')
        receivers = first_run / 'first-seismogram.out.receivers'
        if not receivers.is_file():
            raise RuntimeError(f'{receivers}: not written')
        rows = receivers.read_text().splitlines()
        if len(rows) != 2 or not rows[0].startswith('#') or rows[1].split()[0] != first_trace.name:
            raise RuntimeError(f'{receivers}: expected a header and one row naming {first_trace.name}')
        rows = finite_rows(first_trace, 4)
        if len(rows) != 18001:
            raise RuntimeError(f'{first_trace}: expected 18,001 time samples, got {len(rows)}')
        times = [float(row.split()[0]) for row in rows]
        if any(not math.isclose(time, float(index), rel_tol=0.0, abs_tol=1e-9)
               for index, time in enumerate(times)):
            raise RuntimeError(f'{first_trace}: expected times from 0 to 18,000 s at 1 s intervals')
        if not any(float(row.split()[1]) != 0.0 for row in rows):
            raise RuntimeError(f'{first_trace}: expected a nonzero vertical displacement trace')
        run([build / 'bin/yspec', work / 'workflows/input.json'], env)
        traces = sorted((work / 'workflows').glob('multi.event-*.*'))
        expected = [work / 'workflows' / f'multi.{source}.{receiver}'
                    for source in ('event-a', 'event-b') for receiver in range(1, 5)]
        if sorted(traces) != sorted(expected):
            raise RuntimeError('multi-source run did not create the eight documented filenames')
        for trace in expected:
            rows = finite_rows(trace, 4)
            if len(rows) != 31:
                raise RuntimeError(f'{trace}: expected 31 time samples, got {len(rows)}')

        # Check the documented Python API and package entry points from the
        # exact checkout installed by the caller.
        python_env = env | {'MPLBACKEND': 'Agg'}
        python_script = ROOT / 'docs/site/examples/python-small-calculation.py'
        run([python, python_script], python_env)
        python_seismo = work / 'python-seismogram'
        python_seismo.mkdir()
        shutil.copy2(ROOT / 'docs/site/examples/first-seismogram.json', python_seismo / 'input.json')
        shutil.copy2(ROOT / 'examples/prem.200', python_seismo / 'prem.200')
        run([python.parent / 'yspec', python_seismo / 'input.json'], python_env)
        python_trace = python_seismo / 'first-seismogram.out.1'
        rows = finite_rows(python_trace, 4)
        if len(rows) != 18001:
            raise RuntimeError(f'{python_trace}: expected 18,001 time samples, got {len(rows)}')
        times = [float(row.split()[0]) for row in rows]
        if any(not math.isclose(time, float(index), rel_tol=0.0, abs_tol=1e-9)
               for index, time in enumerate(times)):
            raise RuntimeError(f'{python_trace}: expected times from 0 to 18,000 s at 1 s intervals')
        if not any(float(row.split()[1]) != 0 for row in rows):
            raise RuntimeError(f'{python_trace}: expected a nonzero vertical trace')

        python_snapshot = work / 'python-snapshot'
        python_snapshot.mkdir()
        # The snapshot package wrapper is exercised on the separate snapshot fixture.
        shutil.copy2(ROOT / 'docs/site/examples/snapshot-small.json', python_snapshot / 'input.json')
        shutil.copy2(ROOT / 'examples/prem.200', python_snapshot / 'prem.200')
        run([python.parent / 'yspec-snapshot', python_snapshot / 'input.json'], python_env)
        if len(finite_rows(python_snapshot / 'snapshot.1', 8, skip=2)) != 4:
            raise RuntimeError('Python snapshot wrapper did not produce the four-row fixture')

        cpp_kernels = work / 'cpp-kernels'
        cpp_source = copy_kernel_fixture(cpp_kernels, 'kernels-small.json')
        run([build / 'bin/yspec-kernels', cpp_kernels / 'input.json'], env)
        check_kernel_outputs(cpp_kernels, 'kernel-small', cpp_source,
                             probe_truncation=True)

        python_kernels = work / 'python-kernels'
        python_source = copy_kernel_fixture(python_kernels, 'kernels-small.json')
        run([python.parent / 'yspec-kernels', python_kernels / 'input.json'], python_env)
        check_kernel_outputs(python_kernels, 'kernel-small', python_source)

        plotter = python.parent / 'yspec-plot'
        run([plotter, 'traces', python_seismo / 'input.json', '--receivers', '1',
             '--components', 'ZRT', '-o', python_seismo / 'traces.png'], python_env)
        run([plotter, 'kernels', python_kernels / 'input.json', '--measurement', 'shear-packet',
             '--deepest', '1500', '-o', python_kernels / 'radial.png'], python_env)
        run([plotter, 'banana', python_kernels / 'input.json', '--measurement', 'shear-packet',
             '--parameter', 'vs', '-o', python_kernels / 'section.png'], python_env)
        for image in (python_seismo / 'traces.png', python_kernels / 'radial.png',
                      python_kernels / 'section.png'):
            check_png(image)

        # The Python source build deliberately disables NetCDF. Confirm that
        # volume writing is refused there, while the C++ build writes and
        # exposes the expected model variables when its option is enabled.
        volume_off = work / 'python-volume-off'
        copy_kernel_fixture(volume_off, 'kernels-volume-small.json')
        volume_input = json.loads((volume_off / 'input.json').read_text())
        measurement_names = [item['name'] for item in volume_input['measurements']]
        if measurement_names != ['execution-smoke']:
            raise RuntimeError('the cheap volume fixture must remain explicitly execution-only')
        refused = subprocess.run([python.parent / 'yspec-kernels',
                                  volume_off / 'input.json'], cwd=ROOT,
                                 env=python_env, capture_output=True, text=True)
        if refused.returncode == 0 or 'built without netCDF' not in refused.stderr:
            raise RuntimeError('Python package volume route did not report its NetCDF-disabled boundary')

        netcdf = re.search(r'^YSPEC_WITH_NETCDF:BOOL=(ON|OFF)$', cache, re.M)
        if not netcdf:
            raise RuntimeError('cannot determine C++ build NetCDF setting')
        if netcdf.group(1) == 'ON':
            volume_on = work / 'cpp-volume-on'
            copy_kernel_fixture(volume_on, 'kernels-volume-small.json')
            run([build / 'bin/yspec-kernels', volume_on / 'input.json'], env)
            ncfile = volume_on / 'kernel-volume-small.nc'
            inspect_netcdf = ("import netCDF4,numpy as np,sys; f=netCDF4.Dataset(sys.argv[1]); "
                "g=next(iter(f.groups.values())); "
                "assert {k:len(f.dimensions[k]) for k in ('layer','interface','bound')} == {'layer':13,'interface':13,'bound':2}; "
                "assert {k:len(g.dimensions[k]) for k in ('node','element','latitude','longitude','bound','pole')} == {'node':26,'element':13,'latitude':13,'longitude':25,'bound':2,'pole':2}; "
                "names=('K_vs','K_vp','K_qmu','K_qkappa'); "
                "assert all(n in g.variables and n+'_poles' in g.variables for n in names); "
                "assert all(g.variables[n].dimensions == ('node','latitude','longitude') for n in names); "
                "assert all(g.variables[n+'_poles'].dimensions == ('node','pole') for n in names); "
                "assert all(g.variables[n].units == 'm-3 s' for n in names); "
                "assert all(np.isfinite(g.variables[n][:]).all() and np.any(g.variables[n][:] != 0) for n in names); f.close()")
            run([python, '-c', inspect_netcdf, str(ncfile)], python_env)

        executable = work / 'library-workflows'
        command = [compiler.group(1), '-std=c++23', '-O3', '-I', str(ROOT),
                   str(ROOT / 'docs/site/examples/library-workflows.cpp'),
                   str(build / 'libYSpec.a'), '-lfftw3', '-lfftw3f', '-lfftw3l',
                   '-fopenmp', '-o', str(executable)]
        run(command)
        run([executable, ROOT / 'examples/prem.200'], env)
    print('Checked C++ and Python first-seismogram/snapshot/kernel routes, bounded Python results, radial/section/map/measurement outputs, three headless plots, NetCDF build boundaries, and the strain/kernel library consumer.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
