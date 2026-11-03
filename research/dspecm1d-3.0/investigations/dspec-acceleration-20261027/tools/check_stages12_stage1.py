#!/usr/bin/env python3
"""Validate Stage 1 profile accounting and candidate baseline scope."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path

REQUIRED_PHASES = {
    "setup",
    "initial_frequency_workspace_and_capture",
    "measurement_and_initial_output",
    "kernel_production_excluding_isolated_diagnostics",
    "kernel_output_and_result_io",
    "diagnostic_tensor_validation",
    "profile_pipeline_total",
}
FROZEN_BASELINE_HASHES = {'additional_dependencies_sha256.json': '9c8ca72c9f40fd6b2f18aab27c176f649089f575a3e91cd6d008421906a61a2f', 'case/kernels-small.json': '7b0efae842bff4b414b89eaecf209a986d6e918e0641e147cea0881d46009cc0', 'case/prem.200': '4c3c53ca04e785ff4fd853cf211fe94e353ef70d6052a526122d0f279a3519ea', 'doc_driver.cpp': '220000be65112ea2f895b51110652d2deac1bd298f736573fb31a720fd41e532', 'engine_headers/SensitivityKernels.h': '2fa57eb58a5d937aaae592a0711999f84fae505a12e1936d0958d995850e6b8c', 'engine_headers/Windowed3DFields.h': '084dbd99117a8b26ab4eaa6df06ddd42fc38bd14d73db9555474caa80c61ef49', 'engine_headers/Windowed3DKernels.h': 'a852dead24ba03c1c4139ef926b807a88551ce2c2d77f443ddf149cf150d5506', 'engine_headers/profile.hpp': '2aff92d25d28434fc07302b73397c790a270ab76b5c858e1309785f8a4a09bc1', 'fast_engine.hpp': '3bb252b3e94f3365f02a5dd1822030ad27a0ef764210f3311381ed50386ca73f', 'global/global_gsh_tile.hpp': '0a884177fc688b6a78d95c04209571b90329124f7d7437af5a063887484671e0', 'global/global_kernel_consumer.hpp': 'cfe88cf6b02800152bd3fe6f500e718392b66cc4e9533ee27c4157b2027741f6', 'global/pole_response_spool.hpp': '5e8868439771cb502cb4e3005039f6986e9d0d6ad23a58b168a7a61c2c9152dc', 'implementation_status.md': '3ccddb7738879fbaa823ddb86c389460c49f195fea531c17762c1dfca347744c', 'matched_mesh_model.hpp': '75ab1363d6914bb9e0f893da694dd92e728004018e2786112fbff30c835329cf', 'physical_driver.cpp': 'fbc6f8ffd740ca2a9d727de0317567c8d105e813d31e795341475116d4fbb7e2', 'profile.hpp': '4158ccba5c70b6db31d0aea5aad5755a2a2ffa3700820e74e27e1a92edfdf43a', 'profile_equilibrated_sparse_lu.hpp': '8dd2b4c7be5eb42d35d49d462eb0d0b04a36e5264660df4126abcdac7ee03047', 'solvers/batch_response_cache.hpp': 'c547f7e6f24bfd2cf42edb54eb163387cca42749e9dc677531e01d9246897d85', 'solvers/solver_backends.hpp': '0a55afc975d6a1fdb8a26689a872530b8ee8cd736c044e07e5c990d526b758c4', 'stage-to-math-map.json': '90dcc32f75388ed389b125a48576d336f375f9ce73899b81fc66837fb18dc307'}
ALLOWED_CHANGED = {
    "doc_driver.cpp",
    "profile.hpp",
    "profile_equilibrated_sparse_lu.hpp",
    "engine_headers/Windowed3DFields.h",
    "engine_headers/SensitivityKernels.h",
    "solvers/batch_response_cache.hpp",
    "stage-to-math-map.json",
    "implementation_status.md",
}
REQUIRED_ROW_FIELDS = {
    "phase", "scope_name", "execution_domain", "inclusive_s", "exclusive_s",
    "calls", "logical_component", "diagnostic", "parent_overlap_description",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_profile(record: dict) -> dict:
    if record.get("schema") != 1:
        raise ValueError("unsupported profile accounting schema")
    rows = record.get("rows")
    if not isinstance(rows, list) or not rows:
        raise ValueError("missing profile rows")
    for row in rows:
        missing = REQUIRED_ROW_FIELDS - row.keys()
        if missing:
            raise ValueError(f"profile row is missing fields: {sorted(missing)}")
        if row["execution_domain"] not in {"master_wall", "worker_elapsed_sum"}:
            raise ValueError("wrong execution domain")
        worker_name = row["scope_name"].startswith("parallel_worker.")
        worker_domain = row["execution_domain"] == "worker_elapsed_sum"
        if worker_name != worker_domain:
            raise ValueError("worker scope name and execution domain disagree")
        if not isinstance(row["diagnostic"], bool):
            raise ValueError("diagnostic must be boolean")
        for field in ("inclusive_s", "exclusive_s"):
            value = row[field]
            if not isinstance(value, (float, int)) or not math.isfinite(value) or value < 0:
                raise ValueError(f"invalid {field} for {row['scope_name']}")
        if row["exclusive_s"] > row["inclusive_s"] + 1e-12:
            raise ValueError(f"exclusive duration exceeds inclusive duration for {row['scope_name']}")
        if not isinstance(row["calls"], int) or row["calls"] < 0:
            raise ValueError("invalid scope call count")
        if row["diagnostic"] != (row["phase"] == "diagnostic_tensor_validation"):
            raise ValueError("diagnostic phase tagging is inconsistent")
        if row["execution_domain"] == "worker_elapsed_sum" and "overlaps" not in row["parent_overlap_description"]:
            raise ValueError("worker overlap description is missing")
    phases = {row["phase"] for row in rows}
    if not REQUIRED_PHASES.issubset(phases):
        raise ValueError(f"missing required profile phases: {sorted(REQUIRED_PHASES - phases)}")
    total_rows = [row for row in rows if row["phase"] == "profile_pipeline_total" and row["scope_name"] == "profile_pipeline_wall"]
    if len(total_rows) != 1:
        raise ValueError("missing or duplicate pipeline total row")
    if any(row["phase"] == "kernel_production_excluding_isolated_diagnostics" and row["diagnostic"] for row in rows):
        raise ValueError("diagnostic timings contaminate kernel production phase")
    wall = record.get("wall_partition", {})
    for name in ("profile_pipeline_wall_s", "master_exclusive_partition_s", "residual_s", "tolerance_s"):
        if not isinstance(wall.get(name), (float, int)) or not math.isfinite(wall[name]):
            raise ValueError(f"invalid wall partition field: {name}")
    if wall.get("worker_elapsed_excluded") is not True:
        raise ValueError("worker elapsed sums must be explicitly excluded from wall totals")
    recomputed = sum(row["exclusive_s"] for row in rows if row["execution_domain"] == "master_wall")
    if abs(recomputed - wall["master_exclusive_partition_s"]) > 1e-9:
        raise ValueError("wall partition includes inclusive or worker durations, or omits exclusive entries")
    residual = recomputed - wall["profile_pipeline_wall_s"]
    tolerance = max(0.01, 0.001 * wall["profile_pipeline_wall_s"])
    if abs(residual - wall["residual_s"]) > 1e-9 or abs(residual) > tolerance:
        raise ValueError("exclusive wall partition does not close to profile pipeline wall")
    counters = record.get("counters")
    if not isinstance(counters, list):
        raise ValueError("missing profile counters")
    for counter in counters:
        if not {"phase", "name", "value", "diagnostic"}.issubset(counter):
            raise ValueError("malformed counter row")
        if not isinstance(counter["value"], int) or counter["value"] < 0:
            raise ValueError("invalid counter value")
        if counter["diagnostic"] != (counter["phase"] == "diagnostic_tensor_validation"):
            raise ValueError("diagnostic counter tagging is inconsistent")
    return {"profile_rows": len(rows), "profile_counters": len(counters),
            "pipeline_wall_s": wall["profile_pipeline_wall_s"],
            "master_exclusive_partition_s": recomputed, "wall_residual_s": residual,
            "wall_tolerance_s": tolerance}


def validate_baseline(candidate: Path, manifest_path: Path | None = None) -> dict:
    manifest = json.loads(manifest_path.read_text()) if manifest_path else {"files": FROZEN_BASELINE_HASHES}
    expected = manifest.get("files")
    if not isinstance(expected, dict) or not expected:
        raise ValueError("baseline candidate manifest has no file hashes")
    actual_files = {path.relative_to(candidate).as_posix() for path in candidate.rglob("*") if path.is_file()}
    if actual_files != set(expected):
        raise ValueError("candidate file set differs from frozen baseline")
    changed = []
    for relative, expected_hash in expected.items():
        actual_hash = sha256(candidate / relative)
        if actual_hash == expected_hash:
            continue
        changed.append(relative)
        if relative not in ALLOWED_CHANGED:
            raise ValueError(f"baseline hash changed outside Stage 1 scope: {relative}")
    return {"baseline_files": len(expected), "changed_files": sorted(changed),
            "manifest_sha256": sha256(manifest_path) if manifest_path else "embedded-in-checker"}


def validate_parity(candidate_root: Path, candidate_job: str,
                    reference_root: Path, reference_job: str,
                    expected_candidate_workers: int,
                    expected_reference_workers: int) -> dict:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import check_pair as receipts
    candidate, candidate_hashes, candidate_state_hash = receipts.verified_job(candidate_root, candidate_job)
    reference, reference_hashes, reference_state_hash = receipts.verified_job(reference_root, reference_job)
    candidate_data = candidate_root / candidate_job / "data"
    reference_data = reference_root / reference_job / "data"
    candidate_result = json.loads((candidate_data / "result.json").read_text())
    reference_result = json.loads((reference_data / "result.json").read_text())
    if candidate_result.get("across_angles") != reference_result.get("across_angles"):
        raise ValueError("section angle metadata differs")
    exact_files = ["parameters.txt", "trace.csv", "spectrum.csv", "components.csv", "nodes.csv",
                   "along.bin", "along-angles.csv"]
    for section in ("along", "across"):
        for suffix in (".bin", "-angles.csv"):
            name = section + suffix
            if section == "along" and name in {"along.bin", "along-angles.csv"}:
                continue
            candidate_has = "data/" + name in candidate["outputs"]
            reference_has = "data/" + name in reference["outputs"]
            if candidate_has != reference_has:
                raise ValueError(f"scientific output presence differs: {name}")
            if candidate_has:
                exact_files.append(name)
    if candidate_result.get("across_angles", 0) > 0 and not any(name.startswith("across") for name in exact_files):
        raise ValueError("across-section scientific outputs are required by the result metadata")
    for name in exact_files:
        if "data/" + name not in candidate["outputs"] or "data/" + name not in reference["outputs"]:
            raise ValueError(f"scientific output missing from one parity job: {name}")
        if sha256(candidate_data / name) != sha256(reference_data / name):
            raise ValueError(f"scientific output differs: {name}")
    timing_fields = {"forward_s", "measurement_s", "kernel_s", "tensor_validation_s", "elapsed_s", "workers"}
    if candidate_result.get("workers") != expected_candidate_workers or reference_result.get("workers") != expected_reference_workers:
        raise ValueError("worker count does not match the parity plan")
    candidate_meta = {key: value for key, value in candidate_result.items() if key not in timing_fields}
    reference_meta = {key: value for key, value in reference_result.items() if key not in timing_fields}
    if candidate_meta != reference_meta:
        raise ValueError("scientific result metadata differs")
    return {"candidate_job": candidate_job, "reference_job": reference_job,
            "byte_exact_scientific_files": exact_files,
            "candidate_verified_state_sha256": candidate_state_hash,
            "reference_verified_state_sha256": reference_state_hash,
            "candidate_verified_artifact_sha256": candidate_hashes,
            "reference_verified_artifact_sha256": reference_hashes,
            "candidate_workers": expected_candidate_workers,
            "reference_workers": expected_reference_workers,
            "scientific_metadata_identical": True}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", type=Path)
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--baseline-manifest", type=Path)
    parser.add_argument("--parity-candidate-root", type=Path)
    parser.add_argument("--parity-candidate-job")
    parser.add_argument("--parity-reference-root", type=Path)
    parser.add_argument("--parity-reference-job")
    parser.add_argument("--expected-candidate-workers", type=int)
    parser.add_argument("--expected-reference-workers", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.baseline_manifest and not args.candidate:
        parser.error("--baseline-manifest requires --candidate")
    result = {"schema": 1, "status": "pass", "checks": {}}
    if args.profile:
        result["checks"]["profile_accounting"] = validate_profile(json.loads(args.profile.read_text()))
        result["checks"]["profile_sha256"] = sha256(args.profile)
    if args.candidate:
        result["checks"]["baseline_scope"] = validate_baseline(args.candidate, args.baseline_manifest)
    parity_values = [args.parity_candidate_root, args.parity_candidate_job,
                     args.parity_reference_root, args.parity_reference_job,
                     args.expected_candidate_workers,args.expected_reference_workers]
    if any(value is not None for value in parity_values):
        if any(value is None for value in parity_values):
            parser.error("all four parity arguments are required")
        result["checks"]["scientific_parity"] = validate_parity(
            args.parity_candidate_root, args.parity_candidate_job,
            args.parity_reference_root, args.parity_reference_job,
            args.expected_candidate_workers,args.expected_reference_workers)
    if not result["checks"]:
        parser.error("at least one check is required")
    output = json.dumps(result, indent=2, sort_keys=True)
    if args.output:
        with args.output.open("x") as stream:
            stream.write(output + "\n")
    else:
        print(output)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise SystemExit(f"Stage 1 check failed: {exc}")
