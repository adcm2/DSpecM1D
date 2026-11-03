"""Observe active node RAM and pause at most one owned numerical service.

Uses MemTotal-MemAvailable, then subtracts only explicitly evictable ZFS ARC
data above c_min. Missing or invalid ARC counters disable pause decisions.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'resource_coordination_20261101_v1' / 'global_ram_guard_v4'
SLICE = 'ccdsmextended20261031.slice'
EXPECTED_PAUSED_UNIT = 'ccdsm-sp12-coupling6-64workers-v2.service'
ARCSTATS = Path('/proc/spl/kstat/zfs/arcstats')
PAUSE_BYTES = 825_000_000_000
URGENT_BYTES = 875_000_000_000
BOUNDARY_BYTES = 900_000_000_000


def atomic(path: Path, value: dict) -> None:
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    os.replace(temporary, path)


def memory_metrics(meminfo_text: str, arcstats_text: str) -> dict:
    values = {}
    for line in meminfo_text.splitlines():
        key, raw = line.split(':', 1)
        values[key] = int(raw.strip().split()[0]) * 1024
    total, free, available = values['MemTotal'], values['MemFree'], values['MemAvailable']
    if not (0 <= free <= total and 0 <= available <= total):
        raise ValueError('Invalid MemTotal/MemFree/MemAvailable values')
    arc = {}
    for line in arcstats_text.splitlines():
        fields = line.split()
        if len(fields) < 2:
            continue
        key = fields[0]
        if key in {'size', 'c_min', 'mru_evictable_data', 'mfu_evictable_data'}:
            if len(fields) != 3 or fields[1] != '4' or key in arc:
                raise ValueError(f'Invalid ZFS ARC counter format/type {key}')
            try:
                arc[key] = int(fields[2], 10)
            except ValueError as exc:
                raise ValueError(f'Invalid ZFS ARC counter {key}') from exc
    required = {'size', 'c_min', 'mru_evictable_data', 'mfu_evictable_data'}
    if required - arc.keys():
        raise ValueError('Required ZFS ARC counters are missing: ' + ','.join(sorted(required - arc.keys())))
    if any(value < 0 for value in arc.values()) or arc['c_min'] > arc['size']:
        raise ValueError('Invalid ZFS ARC counter range')
    linux_used = max(0, total - available)
    arc_headroom = max(arc['size'] - arc['c_min'], 0)
    evictable_arc_data = arc['mru_evictable_data'] + arc['mfu_evictable_data']
    arc_subtraction = min(arc_headroom, evictable_arc_data)
    active = max(0, linux_used - arc_subtraction)
    return dict(host_memory_total_bytes=total,
                host_memory_free_bytes=free,
                host_memory_available_estimate_bytes=available,
                host_occupied_cache_inclusive_bytes=total - free,
                host_used_after_linux_available_discount_bytes=linux_used,
                linux_available_estimate_bytes=available,
                zfs_arc_size_bytes=arc['size'],
                zfs_arc_c_min_bytes=arc['c_min'],
                zfs_mru_evictable_data_bytes=arc['mru_evictable_data'],
                zfs_mfu_evictable_data_bytes=arc['mfu_evictable_data'],
                zfs_arc_evictable_data_reported_bytes=evictable_arc_data,
                zfs_arc_above_c_min_bytes=arc_headroom,
                zfs_arc_subtraction_capped_bytes=arc_subtraction,
                active_ram_estimate_bytes=active)


def reconcile_pause(used: int, paused_unit: str | None, pause_armed: bool,
                    frozen_units: set[str], allow_rearm: bool = True) -> dict:
    """Pure guard state transition; all freezes/adoptions are explicit events."""
    actions = []
    if paused_unit and paused_unit not in frozen_units:
        actions.append('MANUAL_THAW_OR_UNIT_EXIT_OBSERVED')
        paused_unit = None
        pause_armed = False
    if len(frozen_units) > 1:
        return dict(paused_unit=paused_unit, pause_armed=False,
                    should_pause=False, conflict='MULTIPLE_OWNED_UNITS_ALREADY_FROZEN', actions=actions)
    if not paused_unit and frozen_units:
        only = next(iter(frozen_units))
        if only == EXPECTED_PAUSED_UNIT:
            paused_unit = only
            pause_armed = False
            actions.append('ADOPT_EXISTING_FROZEN_COUPLING_UNIT')
        else:
            return dict(paused_unit=None, pause_armed=False, should_pause=False,
                        conflict='UNEXPECTED_OWNED_UNIT_ALREADY_FROZEN', actions=actions)
    if allow_rearm and not paused_unit and not frozen_units and not pause_armed and used < PAUSE_BYTES:
        pause_armed = True
        actions.append('GUARD_REARMED_BELOW_PAUSE_THRESHOLD')
    should_pause = bool(allow_rearm and used >= PAUSE_BYTES and pause_armed and not paused_unit and not frozen_units)
    return dict(paused_unit=paused_unit, pause_armed=pause_armed,
                should_pause=should_pause, conflict=None, actions=actions)


def research_cgroup() -> Path:
    relative = Path('/proc/self/cgroup').read_text().strip().split('::', 1)[1]
    here = Path('/sys/fs/cgroup') / relative.lstrip('/')
    # The monitor is outside the research slice; use the user-manager layout
    # that contains the shared slice, matching the existing v2 deployment.
    return here.parent.parent / SLICE


def owned_units() -> list[dict]:
    research = research_cgroup()
    if not research.is_dir():
        raise RuntimeError('Expected owned shared research slice is absent')
    units = []
    registrations = ((research, ('ccdsm-sp12-', 'dspecpilot3d-modal-projection-')),
                     (research.parent / 'dspecscienceangular20261103.slice',
                      ('dspec-science-angular-v1-',)))
    for root, prefixes in registrations:
        # The new slice may not exist before its first admitted service.
        if not root.is_dir():
            continue
        for path in root.iterdir():
            if not path.name.endswith('.service'):
                continue
            if not path.name.startswith(prefixes):
                continue
            if not path.is_dir():
                continue
            try:
                pids = int((path / 'pids.current').read_text())
                memory = int((path / 'memory.current').read_text())
            except (FileNotFoundError, ValueError):
                # Unit removal during a stop is normal; reevaluate next sample.
                continue
            try:
                freeze_text = (path / 'cgroup.freeze').read_text().strip()
                if freeze_text not in ('0', '1'):
                    raise ValueError('Unexpected cgroup.freeze value')
                frozen = freeze_text == '1'
                freeze_error = None
            except (FileNotFoundError, ValueError, OSError) as exc:
                frozen = None
                freeze_error = f'{type(exc).__name__}: {exc}'
            units.append(dict(name=path.name, path=str(path), pids=pids,
                              memory_bytes=memory, frozen=frozen, freeze_error=freeze_error))
    return units

def main() -> None:
    OUT.mkdir(exist_ok=True)
    peak = 0
    paused_unit = None
    pause_armed = True
    actions = []
    started = time.time()
    while True:
        units = owned_units()
        frozen_units = {unit['name'] for unit in units if unit['frozen'] and unit['pids'] > 0}
        freeze_unknown = [unit for unit in units if unit['frozen'] is None and unit['pids'] > 0]
        try:
            metrics = memory_metrics(Path('/proc/meminfo').read_text(), ARCSTATS.read_text())
            used = metrics['active_ram_estimate_bytes']
            metric_error = None
        except (OSError, ValueError, KeyError) as exc:
            # Preserve/adopt the existing pause state without making a pressure
            # decision from a partial metric.
            used = 0
            metric_error = f'{type(exc).__name__}: {exc}'
            metrics = {}
        peak = max(peak, used)
        if freeze_unknown:
            transition = dict(paused_unit=paused_unit, pause_armed=pause_armed,
                              should_pause=False, conflict='FREEZE_STATE_UNAVAILABLE', actions=[])
        elif metric_error:
            transition = reconcile_pause(used, paused_unit, pause_armed, frozen_units,
                                         allow_rearm=False)
        else:
            transition = reconcile_pause(used, paused_unit, pause_armed, frozen_units)
        paused_unit, pause_armed = transition['paused_unit'], transition['pause_armed']
        for action in transition['actions']:
            actions.append(dict(action=action, time_unix=time.time(), active_ram_estimate_bytes=used))
        if transition['conflict']:
            actions.append(dict(action=transition['conflict'], frozen_units=sorted(frozen_units),
                                time_unix=time.time(), active_ram_estimate_bytes=used))
        if transition['should_pause']:
            candidates = [unit for unit in units if unit['pids'] > 0 and not unit['frozen']]
            victim = max(candidates, key=lambda row: row['memory_bytes'], default=None)
            if victim:
                command = ['systemctl', '--user', 'freeze', victim['name']]
                result = subprocess.run(command, text=True, capture_output=True, timeout=8)
                confirmed = False
                try:
                    confirmed = (Path(victim['path']) / 'cgroup.freeze').read_text().strip() == '1'
                except OSError:
                    confirmed = False
                actions.append(dict(action='FREEZE_LARGEST_OWNED_UNIT', command=command,
                    returncode=result.returncode, confirmed_frozen=confirmed,
                    stdout=result.stdout[-2000:], stderr=result.stderr[-2000:],
                    unit_memory_bytes=victim['memory_bytes'], unit=victim['name'],
                    active_ram_estimate_bytes=used, time_unix=time.time()))
                if result.returncode == 0 and confirmed:
                    paused_unit = victim['name']
                    pause_armed = False
        actions = actions[-32:]
        status = ('METRIC_UNAVAILABLE' if metric_error else
                  'FREEZE_STATE_UNAVAILABLE' if freeze_unknown else
                  'URGENT_NODE_PRESSURE' if used >= URGENT_BYTES else
                  'SINGLE_JOB_PAUSED' if paused_unit else
                  'PAUSE_CONFLICT_REQUIRES_OPERATOR' if transition['conflict'] else 'MONITORING')
        state = dict(schema='tunya-global-memory-observer-v4', started_unix=started,
            measured_unix=time.time(), status=status,
            metric='max(0, MemTotal - MemAvailable - min(max(ARC.size - ARC.c_min,0), ARC.mru_evictable_data + ARC.mfu_evictable_data)); excludes only reported evictable ARC data and Linux MemAvailable estimate',
            active_ram_metric_valid=metric_error is None,
            active_ram_metric_error=metric_error,
            **metrics, active_ram_peak_observed_bytes=peak,
            pause_bytes=PAUSE_BYTES, urgent_bytes=URGENT_BYTES,
            user_global_boundary_bytes=BOUNDARY_BYTES,
            owned_work_ceiling_bytes=600_000_000_000,
            paused_unit=paused_unit, pause_armed=pause_armed,
            existing_frozen_owned_units=sorted(frozen_units),
            freeze_state_unknown_units=[dict(name=u['name'], error=u['freeze_error']) for u in freeze_unknown],
            largest_active_owned_unit=max((u for u in units if u['pids'] > 0 and not u['frozen']),
                key=lambda row: row['memory_bytes'], default={}).get('name'),
            largest_active_owned_memory_bytes=max((u['memory_bytes'] for u in units if u['pids'] > 0 and not u['frozen']), default=0),
            action_policy='Thresholds apply to active RAM estimate, not cache-inclusive occupancy. Pause only the largest active owned numerical service at 825GB; preserve existing pause state; require pressure to fall below threshold before rearming after a manual thaw or unit exit. Invalid/missing ARC stats disable pressure actions.',
            prohibited_actions=['automatic thaw', 'kill', 'shared-slice freeze', 'whole-slice stop'],
            actions=actions,
            limitation='User-level guard cannot impose a machine-wide kernel cap on other users.')
        atomic(OUT / 'state.json', state)
        if actions or status == 'URGENT_NODE_PRESSURE':
            atomic(OUT / 'resource_action.json', state)
        time.sleep(0.5)


if __name__ == '__main__':
    main()
