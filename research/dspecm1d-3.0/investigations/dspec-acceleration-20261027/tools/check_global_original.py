#!/usr/bin/env python3
"""Compare a direct original job with a full-grid global-spool candidate."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import check_global_cache
import check_global_pair as reviewed
import check_pair


def _provenance() -> dict[str, Any]:
    helper = Path(reviewed.__file__)
    return {
        'reviewed_helper': helper.name,
        'reviewed_helper_sha256': check_pair.sha(helper),
        'timing_key_adapter': {
            'source_key': 'serial_grid_and_context_setup_wall_s',
            'compatibility_key': 'serial_context_setup_wall_s',
            'method': 'in-memory alias; original timing fields retained unchanged',
        },
    }


def assess(reference_root: Path, reference_job: str,
           candidate_root: Path, candidate_job: str) -> dict[str, Any]:
    original_validator = reviewed.validate_worker_timings
    adaptations: list[dict[str, Any]] = []

    def adapted_validator(manifest: dict) -> dict:
        before = dict(manifest.get('timings', {}))
        adapted = check_global_cache._final_timing_adapter(manifest)
        after = adapted.get('timings', {})
        for key, value in before.items():
            if key not in after or after[key] != value:
                raise ValueError('timing adapter altered original timing field: ' + key)
        source = before.get('serial_grid_and_context_setup_wall_s')
        compatibility = after.get('serial_context_setup_wall_s')
        if source is None or compatibility != source:
            raise ValueError('final global setup timing alias is missing or differs')
        if 'serial_context_setup_wall_s' in before:
            adaptations.append({'applied': False, 'original_fields_preserved': True})
        else:
            adaptations.append({'applied': True, 'original_fields_preserved': True,
                                'source_value_s': source,
                                'compatibility_value_s': compatibility})
        return original_validator(adapted)

    reviewed.validate_worker_timings = adapted_validator
    try:
        report = reviewed.assess(reference_root, reference_job,
                                 candidate_root, candidate_job,
                                 probe64=False, global_reference=False)
    finally:
        reviewed.validate_worker_timings = original_validator
    report['timing_adaptation'] = adaptations
    report.update(_provenance())
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('reference-root', 'reference-job', 'candidate-root', 'candidate-job', 'output'):
        parser.add_argument('--' + name, required=True)
    args = parser.parse_args()
    try:
        report = assess(Path(args.reference_root), args.reference_job,
                        Path(args.candidate_root), args.candidate_job)
    except Exception as exc:
        report = {'status': 'FAIL', 'error_type': type(exc).__name__, 'error': str(exc),
                  'timing_adaptation': [], **_provenance()}
    with Path(args.output).open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    return 0 if report['status'] == 'PASS' else 1


if __name__ == '__main__':
    raise SystemExit(main())
