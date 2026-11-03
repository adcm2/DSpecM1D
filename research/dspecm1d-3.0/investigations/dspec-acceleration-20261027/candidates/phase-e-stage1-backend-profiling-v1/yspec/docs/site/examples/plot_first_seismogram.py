#!/usr/bin/env python3
"""Plot the vertical trace from a yspec four-column output file.

Usage: plot_first_seismogram.py TRACE OUTPUT [INPUT_JSON]
The JSON is optional for compatibility; when supplied, it labels the quantity
and derives the plotted receiver distance from the run's receiver table.
"""
import json
import math
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if len(sys.argv) not in (3, 4):
    raise SystemExit(__doc__.splitlines()[2])
source, destination = map(Path, sys.argv[1:3])
input_path = Path(sys.argv[3]) if len(sys.argv) == 4 else None
rows = []
for line in source.read_text().splitlines():
    fields = line.split()
    if fields:
        rows.append(tuple(map(float, fields)))
if not rows or any(len(row) != 4 for row in rows):
    raise SystemExit('expected four numeric columns per output row')
if any(not math.isfinite(value) for row in rows for value in row):
    raise SystemExit('output contains a non-finite value')
times, vertical = [r[0] for r in rows], [r[1] for r in rows]
if max(map(abs, vertical)) == 0:
    raise SystemExit('vertical trace is identically zero')

quantity = None
unit = None
distance = None
if input_path is not None:
    settings = json.loads(input_path.read_text())
    quantity = settings.get('output', {}).get('quantity', 'displacement')
    unit = {'displacement': 'm', 'velocity': 'm/s',
            'acceleration': 'm/s²'}.get(quantity)
    if unit is None:
        raise SystemExit(f'unsupported output quantity: {quantity}')
    output_prefix = input_path.parent / settings['output']['prefix']
    receiver_table = output_prefix.with_name(output_prefix.name + '.receivers')
    if receiver_table.exists():
        lines = [line.split() for line in receiver_table.read_text().splitlines()
                 if line.strip()]
        header = lines[0][1:] if lines and lines[0][0] == '#' else []
        if 'file' not in header or 'source' not in header:
            raise SystemExit(f'{receiver_table}: missing receiver table header')
        file_column, source_column = header.index('file'), header.index('source')
        distance_column = (header.index('distance_deg')
                           if 'distance_deg' in header else None)
        receiver_rows = lines[1:]
        matching = [row for row in receiver_rows
                    if len(row) > max(file_column, source_column) and
                    (output_prefix.parent / row[file_column]).resolve() ==
                    source.resolve()]
        if len(matching) != 1:
            raise SystemExit(f'{source}: expected one matching receiver row in '
                             f'{receiver_table}, found {len(matching)}')
        selected = matching[0]
        source_name = selected[source_column]
        source_rows = [row for row in receiver_rows
                       if len(row) > source_column and
                       row[source_column] == source_name]
        receiver_index = source_rows.index(selected)
        if distance_column is not None:
            distance = float(selected[distance_column])

        receiver_groups = settings.get('receivers', [])
        if isinstance(receiver_groups, dict):
            receiver_groups = [receiver_groups]
        units = []
        for group in receiver_groups:
            depths = group['depth_km']
            depths = len(depths) if isinstance(depths, list) else 1
            places = (len(group['locations']) if 'locations' in group
                      else group['profile']['count'])
            instrument_unit = False
            if 'instrument' in group:
                instrument_path = input_path.parent / group['instrument']
                try:
                    instrument = instrument_path.read_text(errors='replace')
                except OSError:
                    instrument = ''
                match = re.search(
                    r'^\s*\*\s*OUTPUT\s+UNIT\s*:\s*(\S.*?)\s*$',
                    instrument, flags=re.IGNORECASE | re.MULTILINE)
                instrument_unit = match.group(1).strip() if match else None
            units.extend([instrument_unit] * (depths * places))
        if receiver_index >= len(units):
            raise SystemExit(f'{source}: receiver row does not match the input '
                             'receiver groups')
        instrument_unit = units[receiver_index]
        unit = ((instrument_unit or 'instrument output (unit unspecified)')
                if instrument_unit is not False else unit)
    elif any(group.get('instrument') for group in
             ([settings.get('receivers', [])]
              if isinstance(settings.get('receivers', []), dict)
              else settings.get('receivers', []))):
        raise SystemExit(f'{receiver_table}: required to resolve receiver units')
else:
    receiver_table = source.with_name(
        f'{source.name.rsplit(".", 1)[0]}.receivers')
    if receiver_table.exists():
        lines = [line.split() for line in receiver_table.read_text().splitlines()
                 if line.strip()]
        if lines and lines[0][0] == '#':
            header = lines[0][1:]
            if 'file' in header and 'distance_deg' in header:
                file_column, distance_column = (header.index('file'),
                                                header.index('distance_deg'))
                match = next((row for row in lines[1:]
                              if len(row) > max(file_column, distance_column) and
                              Path(row[file_column]).name == source.name), None)
                if match is not None:
                    distance = float(match[distance_column])

quantity_label = quantity or 'vertical component'
title = quantity_label.capitalize()
if distance is not None:
    title += f' · {distance:g}° epicentral distance'
ylabel = 'Up component' + (f' ({unit})' if unit else '')
matplotlib.rcParams['svg.hashsalt'] = 'yspecpp-first-seismogram'
fig, ax = plt.subplots(figsize=(8.2, 3.6), constrained_layout=True)
ax.plot(times, vertical, color='#176b70', linewidth=1.1)
ax.set(title=title, xlabel='Time after origin (s)', ylabel=ylabel)
ax.grid(True, color='#e3e9e7', linewidth=.7)
ax.margins(x=0)
destination.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(destination, format=destination.suffix.lstrip('.') or 'svg',
            metadata={'Creator': 'YSpecpp tutorial plot script', 'Date': None})
svg = destination.read_text()
destination.write_text(''.join(
    line.rstrip('\r\n').rstrip(' \t') + line[len(line.rstrip('\r\n')):]
    for line in svg.splitlines(keepends=True)
))
peak_unit = f' {unit}' if unit else ' (output unit unspecified)'
print(f'{len(rows)} samples; vertical peak absolute amplitude '
      f'{max(map(abs, vertical)):.9g}{peak_unit}; wrote {destination}')
