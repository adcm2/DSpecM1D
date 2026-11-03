#!/usr/bin/env python3
"""Stage 2 frequency-route profile, parity, invariant, and candidate-scope gates."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path

BASELINE = Path(__file__).resolve().parents[1] / "analyses/stages12-stage2-v1-baseline-candidate-sha256.json"
FROZEN_BASELINE_SHA256 = {'physical_driver.cpp': 'fbc6f8ffd740ca2a9d727de0317567c8d105e813d31e795341475116d4fbb7e2', 'profile_equilibrated_sparse_lu.hpp': 'cad08d7bad6e71a87f9205f0e533cc6ee2106dbe766497c7992baf1e50cd8324', 'stage-to-math-map.json': 'e8deeef0c7354491c2deca3e0950139fc877a02c43fe6489ce063278390fbd29', 'matched_mesh_model.hpp': '75ab1363d6914bb9e0f893da694dd92e728004018e2786112fbff30c835329cf', 'implementation_status.md': 'f8e7aca596933628cd2b32b3a38b685681241e284853fb620d07f83d39c3004b', 'fast_engine.hpp': '3bb252b3e94f3365f02a5dd1822030ad27a0ef764210f3311381ed50386ca73f', 'additional_dependencies_sha256.json': '9c8ca72c9f40fd6b2f18aab27c176f649089f575a3e91cd6d008421906a61a2f', 'doc_driver.cpp': '465e2203e39e55820ac17555737a6e07ebe7f90ba696c38356387c0dc2f10256', 'profile.hpp': 'ece6b2df5b8ed5b7d97de732ff619bb09d9482c35b99fccf756d5be3a9991724', 'engine_headers/Windowed3DKernels.h': 'a852dead24ba03c1c4139ef926b807a88551ce2c2d77f443ddf149cf150d5506', 'engine_headers/Windowed3DFields.h': '4d9f3d2b4f6fb6540b89a0bc74d7f6d1ee2e47f2cb2bcaa21b36624a6ba0f432', 'engine_headers/SensitivityKernels.h': '2fa57eb58a5d937aaae592a0711999f84fae505a12e1936d0958d995850e6b8c', 'engine_headers/profile.hpp': '2aff92d25d28434fc07302b73397c790a270ab76b5c858e1309785f8a4a09bc1', 'solvers/batch_response_cache.hpp': '576c80b7d342fbe1f9e41fe667e6a1db94d59825c85364cc0cf06e806ccb7750', 'solvers/solver_backends.hpp': '0a55afc975d6a1fdb8a26689a872530b8ee8cd736c044e07e5c990d526b758c4', 'global/global_gsh_tile.hpp': '0a884177fc688b6a78d95c04209571b90329124f7d7437af5a063887484671e0', 'global/pole_response_spool.hpp': '5e8868439771cb502cb4e3005039f6986e9d0d6ad23a58b168a7a61c2c9152dc', 'global/global_kernel_consumer.hpp': 'cfe88cf6b02800152bd3fe6f500e718392b66cc4e9533ee27c4157b2027741f6', 'case/kernels-small.json': '7b0efae842bff4b414b89eaecf209a986d6e918e0641e147cea0881d46009cc0', 'case/prem.200': '4c3c53ca04e785ff4fd853cf211fe94e353ef70d6052a526122d0f279a3519ea'}
FROZEN_BASELINE_MANIFEST_SHA256 = 'f9e29601f520c253f618f2d7ca386adcf3263e785401ccc4e1ec2bc7267fe7c4'
ALLOWED_MODIFIED = {"doc_driver.cpp", "profile.hpp", "stage-to-math-map.json", "implementation_status.md"}
ADDED = {"frequency_scheduler.hpp", "frequency_scheduler_test.cpp"}
REQUIRED_ROWS = {"setup", "initial_frequency_workspace_and_capture",
                 "measurement_and_initial_output",
                 "kernel_frequency_parallel_pass_and_ordered_reduction_including_diagnostics",
                 "kernel_output_and_result_io", "diagnostic_tensor_validation",
                 "profile_pipeline_total"}
REQUIRED_ROW_FIELDS = {"phase", "scope_name", "execution_domain", "inclusive_s",
                      "exclusive_s", "calls", "logical_component", "diagnostic",
                      "parent_overlap_description"}
SCIENTIFIC_FILES = ["parameters.txt", "trace.csv", "spectrum.csv", "components.csv",
                    "nodes.csv", "along.bin", "along-angles.csv", "across.bin",
                    "across-angles.csv"]
TIMING_METADATA = {"forward_s", "measurement_s", "kernel_s", "tensor_validation_s",
                   "elapsed_s", "kernel_parallel_pass_wall_s", "reduction_wall_s",
                   "diagnostic_worker_elapsed_s", "frequency_contribution_buffer_bytes",
                   "workers", "schedule_mode", "frequency_workers", "solve_workers",
                   "kernel_timing_semantics"}
FROZEN_MEASURE_FUNCTION_SHA256 = "7b5aa2a57fb90d4ea858159bdcda5312524bb9d687127683e673a185183cce34"
FROZEN_KERNEL_FORMULA_SHA256 = "99687637e94756b9342cb758570a5e165d3fea087bd674dd53532e753412e295"
EXPECTED_COUNTERS = {
    "sparse_factorization_calls": 14070,
    "solver_forward_rhs_columns": 42070,
    "solver_unit_basis_rhs_columns": 21000,
    "solver_weighted_adjoint_rhs_columns": 70,
    "batch_response_cache_misses": 14000,
    "batch_response_cache_hits": 14000,
    "response_cache_radial_l0_bypass_calls": 70,
    "zero_weighted_adjoint_rhs_skips": 70,
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()



def validate_build(build_result_path: Path, build_config_path: Path,
                   variant_paths: dict[str, Path]) -> dict:
    result = json.loads(build_result_path.read_text())
    config = json.loads(build_config_path.read_text())
    if result.get("status") != "completed" or result.get("build_root") != config.get("build_root"):
        raise ValueError("build result is not a completed receipt for the configured build root")
    config_ids = {v["id"] for v in config.get("variants", [])}
    result_variants = {v["id"]:v for v in result.get("variants", [])}
    if set(variant_paths) != config_ids or set(result_variants) != config_ids:
        raise ValueError("build variant inventory differs from config or supplied binaries")
    verified={}
    for variant_id,binary in variant_paths.items():
        record=result_variants[variant_id]
        command_path=binary.parent/"command.json"
        command=json.loads(command_path.read_text())
        if record.get("returncode") != 0 or command.get("returncode") != 0 or not command.get("tracked_inputs_unchanged"):
            raise ValueError(f"build variant failed or changed tracked inputs: {variant_id}")
        digest=sha256(binary)
        if digest != record.get("binary_sha256") or digest != command.get("binary_sha256"):
            raise ValueError(f"build binary checksum mismatch: {variant_id}")
        if record.get("command_sha256") != command.get("command_sha256"):
            raise ValueError(f"build command checksum differs: {variant_id}")
        verified[variant_id]={"binary":str(binary),"binary_sha256":digest,
                              "command_sha256":command["command_sha256"]}
    return {"build_root":config["build_root"],"variants":verified}


def validate_candidate(candidate: Path, manifest_path: Path | None = None) -> dict:
    if manifest_path is not None:
        manifest = json.loads(manifest_path.read_text())
        baseline = manifest["files"]
        manifest_hash = sha256(manifest_path)
    else:
        baseline = FROZEN_BASELINE_SHA256
        manifest_hash = FROZEN_BASELINE_MANIFEST_SHA256
    expected_files = set(baseline) | ADDED
    actual_files = {p.relative_to(candidate).as_posix() for p in candidate.rglob("*") if p.is_file()}
    if actual_files != expected_files:
        raise ValueError(f"candidate file inventory differs: missing={sorted(expected_files-actual_files)} extra={sorted(actual_files-expected_files)}")
    changed = []
    for rel, expected_hash in baseline.items():
        if sha256(candidate / rel) != expected_hash:
            changed.append(rel)
            if rel not in ALLOWED_MODIFIED:
                raise ValueError(f"frozen Stage 1 candidate changed outside Stage 2 scope: {rel}")
    return {"baseline_file_count": len(baseline), "added_files": sorted(ADDED),
            "modified_baseline_files": sorted(changed),
            "baseline_manifest_sha256": manifest_hash}


def validate_result(result: dict, expected_schedule_mode: str,
                    expected_solve_workers: int, expected_workers: int,
                    expected_frequency_workers: int) -> dict:
    expected = {
        "schedule_mode": expected_schedule_mode,
        "solve_workers": expected_solve_workers,
        "workers": expected_workers,
        "frequency_workers": expected_frequency_workers,
    }
    if expected_schedule_mode not in ("degree", "frequency"):
        raise ValueError("expected schedule mode must be degree or frequency")
    for key, value in expected.items():
        actual = result.get(key)
        if key != "schedule_mode":
            if not isinstance(value, int) or isinstance(value, bool):
                raise ValueError(f"invalid expected worker count: {key}")
            if not isinstance(actual, int) or isinstance(actual, bool):
                raise ValueError(f"invalid result worker count: {key}")
        if actual != value:
            raise ValueError(f"result {key} differs from the plan")
    if expected_schedule_mode == "frequency":
        if not 1 <= expected_frequency_workers <= 8:
            raise ValueError("frequency_workers is outside the [1,8] bound")
        if expected_solve_workers != 1 or expected_workers != 1:
            raise ValueError("frequency mode must use one inner solve worker")
        if "diagnostics" not in result.get("kernel_timing_semantics", "") or "reduction" not in result.get("kernel_timing_semantics", ""):
            raise ValueError("kernel timing semantics do not include diagnostic and reduction work")
        for key in ("kernel_parallel_pass_wall_s", "reduction_wall_s", "diagnostic_worker_elapsed_s", "kernel_s", "elapsed_s"):
            value = result.get(key)
            if not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                raise ValueError(f"invalid frequency timing field: {key}")
        buffer_bytes = result.get("frequency_contribution_buffer_bytes")
        if not isinstance(buffer_bytes, int) or isinstance(buffer_bytes, bool) or buffer_bytes < 0:
            raise ValueError("invalid frequency contribution buffer byte count")
        if buffer_bytes > 12634160:
            raise ValueError("frequency contribution buffer exceeds the audited bound")
    else:
        if (expected_solve_workers, expected_workers, expected_frequency_workers) != (8, 8, 1):
            raise ValueError("degree anchor expectation must be solve=8, workers=8, frequency_workers=1")
        if result.get("kernel_timing_semantics") != "legacy degree kernel_s excludes tensor validation":
            raise ValueError("degree anchor does not use the legacy kernel timing label")
        for key in ("kernel_parallel_pass_wall_s", "reduction_wall_s", "diagnostic_worker_elapsed_s", "frequency_contribution_buffer_bytes"):
            value = result.get(key)
            if value != 0:
                raise ValueError(f"degree anchor must have zero frequency-only field: {key}")
        for key in ("kernel_s", "elapsed_s"):
            value = result.get(key)
            if not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                raise ValueError(f"invalid degree timing field: {key}")
    if result.get("tensor_packed_relative_l2", math.inf) > 1e-11:
        raise ValueError("tensor/packed regression exceeds 1e-11")
    if abs(result.get("time_shift_adjoint_identity", 0.0) - 1.0) > 1e-10:
        raise ValueError("time-shift adjoint identity failed")
    if result.get("gravity_max_abs_nd_error", math.inf) > 1e-10:
        raise ValueError("gravity audit failed")
    return {"schedule_mode": expected_schedule_mode,
            "workers": expected_workers,
            "solve_workers": expected_solve_workers,
            "frequency_workers": expected_frequency_workers,
            "contribution_buffer_bytes": result["frequency_contribution_buffer_bytes"]}



def validate_source_contract(candidate: Path) -> dict:
    source = (candidate / "doc_driver.cpp").read_text()
    measure = source[source.index("Measured measure("):source.index("// Complex weak-form contractions")]
    formula = source[source.index("double kernelVs("):source.index("int main(")]
    if hashlib.sha256(measure.encode()).hexdigest() != FROZEN_MEASURE_FUNCTION_SHA256:
        raise ValueError("measurement/FFT function changed from frozen Stage 1 source")
    if hashlib.sha256(formula.encode()).hexdigest() != FROZEN_KERNEL_FORMULA_SHA256:
        raise ValueError("kernel formula changed from frozen Stage 1 source")
    required = [
        "DSPEC_FREQUENCY_WORKERS", "runIndexedWithContext(grid.Count(),frequencyWorkers",
        "prepareFrequencyWorkspace(sem,source,parameters,z,Eigen::Vector3cd::Zero(),&cache,true,false,nullptr,solveWorkers",
        "prepareFrequencyWorkspace(sem,source,parameters,z,weights,&cache,true,false,nullptr,solveWorkers",
        "std::vector<std::vector<Plan>> plans", "plans.at(static_cast<std::size_t>(omp_get_thread_num())).at(section)",
        "contribution[n*sampling.places.size()+a]=kernelVs(",
        "for(std::size_t k=0;k<grid.Count();++k)for(std::size_t section=0;section<samplings.size();++section)",
        "destination[n*samplings[section].places.size()+a]+=contribution[n*samplings[section].places.size()+a]",
        "double kernelParallelPassSeconds=0,reductionSeconds=0,diagnosticWorkerElapsedSeconds=0;std::size_t contributionBufferBytes=0;"
    ]
    missing = [needle for needle in required if needle not in source]
    if missing:
        raise ValueError(f"frequency ownership/reduction contract missing: {missing}")
    if source.count("&cache,true,false,nullptr,solveWorkers") != 4:
        raise ValueError("both passes must use the stable shared response cache and one solve thread")
    if "values[section][n*sampling.places.size()+a]+=kernelVs" not in source:
        raise ValueError("legacy degree route accumulation is missing")
    return {"measurement_function_sha256": FROZEN_MEASURE_FUNCTION_SHA256,
            "kernel_formula_sha256": FROZEN_KERNEL_FORMULA_SHA256,
            "shared_response_cache_passes": 2, "frequency_contribution_reduction": "ascending k, section, angle, native node"}


def validate_profile(record: dict, require_l200_invariants: bool = False) -> dict:
    if record.get("schema") != 1 or record.get("frequency_route") is not True:
        raise ValueError("profile is not a Stage 2 frequency profile")
    if record.get("diagnostic_worker_elapsed_is_not_subtracted_from_wall") is not True:
        raise ValueError("diagnostic worker elapsed must not be subtracted from wall time")
    diagnostic_partition=record.get("worker_diagnostic_partition", {})
    if diagnostic_partition.get("production_rows_exclude_diagnostic_overlay") is not True or diagnostic_partition.get("raw_inclusive_umbrellas_preserve_diagnostics") is not True or diagnostic_partition.get("diagnostic_overlay_phase") != "diagnostic_tensor_validation":
        raise ValueError("worker diagnostic overlay/raw umbrella semantics are not explicit")
    rows = record.get("rows")
    if not isinstance(rows, list) or not rows:
        raise ValueError("missing profile rows")
    for row in rows:
        missing = REQUIRED_ROW_FIELDS - row.keys()
        if missing:
            raise ValueError(f"profile row lacks fields: {sorted(missing)}")
        name = row["scope_name"]
        worker = name.startswith(("frequency_worker.", "parallel_worker."))
        if worker != (row["execution_domain"] == "worker_elapsed_sum"):
            raise ValueError(f"worker scope/domain mismatch: {name}")
        if name.startswith("parallel_worker."):
            raise ValueError("frequency route unexpectedly includes degree-worker scopes")
        for field in ("inclusive_s", "exclusive_s"):
            value = row[field]
            if not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                raise ValueError(f"invalid {field}: {name}")
        if row["exclusive_s"] > row["inclusive_s"] + 1e-12:
            raise ValueError(f"exclusive exceeds inclusive: {name}")
        if not isinstance(row["calls"], int) or row["calls"] < 0:
            raise ValueError(f"invalid calls: {name}")
        if row["diagnostic"] != (row["phase"] == "diagnostic_tensor_validation"):
            raise ValueError(f"diagnostic row/phase mismatch: {name}")
        if row["logical_component"] == "diagnostic_validation" and row["phase"] != "diagnostic_tensor_validation":
            raise ValueError(f"production scope was classified as pure diagnostic: {name}")
        if row["execution_domain"] == "worker_elapsed_sum" and "overlaps" not in row["parent_overlap_description"]:
            raise ValueError(f"missing worker overlap description: {name}")
    phases = {row["phase"] for row in rows}
    if not REQUIRED_ROWS.issubset(phases):
        raise ValueError(f"missing profile phase(s): {sorted(REQUIRED_ROWS-phases)}")
    if not any(row["scope_name"] == "driver_frequency_parallel_kernel_pass_including_diagnostics" and row["execution_domain"] == "master_wall" for row in rows):
        raise ValueError("missing frequency parallel-pass wall scope")
    if not any(row["scope_name"] == "driver_frequency_ordered_kernel_contribution_reduction" and row["execution_domain"] == "master_wall" for row in rows):
        raise ValueError("missing parent ordered-reduction scope")
    wall = record.get("wall_partition", {})
    for key in ("profile_pipeline_wall_s", "master_exclusive_partition_s", "residual_s", "tolerance_s"):
        if not isinstance(wall.get(key), (int, float)) or not math.isfinite(wall[key]):
            raise ValueError(f"invalid wall partition {key}")
    if wall.get("worker_elapsed_excluded") is not True:
        raise ValueError("worker elapsed sums are not explicitly excluded from wall")
    recomputed = sum(row["exclusive_s"] for row in rows if row["execution_domain"] == "master_wall")
    if abs(recomputed-wall["master_exclusive_partition_s"]) > 1e-9:
        raise ValueError("master wall partition includes worker or inclusive time")
    residual = recomputed-wall["profile_pipeline_wall_s"]
    if abs(residual-wall["residual_s"]) > 1e-9 or abs(residual) > max(.01, .001*wall["profile_pipeline_wall_s"]):
        raise ValueError("master exclusive wall partition does not close")
    counters = record.get("counters")
    if not isinstance(counters, list):
        raise ValueError("missing counters")
    sums: dict[str, int] = {}
    for row in counters:
        if not {"phase", "name", "value", "diagnostic"}.issubset(row):
            raise ValueError("malformed counter row")
        if not isinstance(row["value"], int) or row["value"] < 0:
            raise ValueError("invalid counter value")
        if row["diagnostic"] != (row["phase"] == "diagnostic_tensor_validation"):
            raise ValueError("diagnostic counter/phase mismatch")
        if row["phase"] != "diagnostic_tensor_validation":
            sums[row["name"]] = sums.get(row["name"], 0) + row["value"]
    if require_l200_invariants:
        for name, expected in EXPECTED_COUNTERS.items():
            if sums.get(name, 0) != expected:
                raise ValueError(f"Stage 1 logical counter invariant differs for {name}: {sums.get(name, 0)} != {expected}")
    return {"profile_rows": len(rows), "profile_counters": len(counters),
            "pipeline_wall_s": wall["profile_pipeline_wall_s"],
            "master_exclusive_partition_s": recomputed, "wall_residual_s": residual,
            "worker_elapsed_excluded": True,
            "logical_counter_invariants": EXPECTED_COUNTERS if require_l200_invariants else "not-required-for-low-degree-control"}


def validate_parity(candidate_data: Path, reference_data: Path,
                    expected_schedule_mode: str, expected_solve_workers: int,
                    expected_workers: int, expected_frequency_workers: int) -> dict:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import check_pair as receipts
    candidate_job_dir, reference_job_dir = candidate_data.parent, reference_data.parent
    if candidate_data.name != "data" or reference_data.name != "data":
        raise ValueError("parity paths must point at a task data directory")
    candidate, candidate_hashes, candidate_state_hash = receipts.verified_job(candidate_job_dir.parent, candidate_job_dir.name)
    reference, reference_hashes, reference_state_hash = receipts.verified_job(reference_job_dir.parent, reference_job_dir.name)
    candidate_result = json.loads((candidate_data / "result.json").read_text())
    reference_result = json.loads((reference_data / "result.json").read_text())
    validate_result(candidate_result, expected_schedule_mode, expected_solve_workers,
                     expected_workers, expected_frequency_workers)
    if candidate_result.get("across_angles") != reference_result.get("across_angles"):
        raise ValueError("section angle metadata differs")
    common = ["parameters.txt", "trace.csv", "spectrum.csv", "components.csv", "nodes.csv",
              "along.bin", "along-angles.csv"]
    exact_files = list(common)
    across_expected = candidate_result.get("across_angles", 0) > 0
    for name in ("across.bin", "across-angles.csv"):
        candidate_has = "data/" + name in candidate["outputs"]
        reference_has = "data/" + name in reference["outputs"]
        if candidate_has != reference_has:
            raise ValueError(f"scientific output presence differs: {name}")
        if across_expected and not candidate_has:
            raise ValueError(f"across-section output required by metadata: {name}")
        if not across_expected and candidate_has:
            raise ValueError(f"unexpected across-section output: {name}")
        if candidate_has:
            exact_files.append(name)
    for name in common:
        if "data/" + name not in candidate["outputs"] or "data/" + name not in reference["outputs"]:
            raise ValueError(f"required parity output missing from receipt: {name}")
    for name in exact_files:
        if sha256(candidate_data / name) != sha256(reference_data / name):
            raise ValueError(f"scientific output differs: {name}")
    candidate_meta = {k:v for k,v in candidate_result.items() if k not in TIMING_METADATA}
    reference_meta = {k:v for k,v in reference_result.items() if k not in TIMING_METADATA}
    if candidate_meta != reference_meta:
        raise ValueError("scientific result metadata differs")
    return {"byte_exact_scientific_files": exact_files,
            "candidate_schedule_mode": candidate_result["schedule_mode"],
            "candidate_workers": candidate_result["workers"],
            "candidate_solve_workers": candidate_result["solve_workers"],
            "candidate_frequency_workers": candidate_result["frequency_workers"],
            "reference_schedule_mode": reference_result.get("schedule_mode", "stage1-degree"),
            "candidate_verified_state_sha256": candidate_state_hash,
            "reference_verified_state_sha256": reference_state_hash,
            "candidate_verified_artifact_sha256": candidate_hashes,
            "reference_verified_artifact_sha256": reference_hashes}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate")
    parser.add_argument("--build-result")
    parser.add_argument("--build-config")
    parser.add_argument("--variant", action="append", default=[], help="variant_id=/absolute/path/to/binary")
    parser.add_argument("--profile")
    parser.add_argument("--result")
    parser.add_argument("--expected-schedule-mode", choices=("degree", "frequency"))
    parser.add_argument("--expected-solve-workers", type=int)
    parser.add_argument("--expected-workers", type=int)
    parser.add_argument("--expected-frequency-workers", type=int)
    parser.add_argument("--require-l200-invariants", action="store_true")
    parser.add_argument("--parity-candidate-data")
    parser.add_argument("--parity-reference-data")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    report = {}
    if args.build_result and args.build_config:
        variants={}
        for item in args.variant:
            variant_id,separator,path=item.partition("=")
            if not separator or variant_id in variants:raise ValueError("variant arguments must be unique ID=PATH pairs")
            variants[variant_id]=Path(path)
        report["build"] = validate_build(Path(args.build_result),Path(args.build_config),variants)
    if args.candidate:
        report["candidate"] = validate_candidate(Path(args.candidate))
        report["source_contract"] = validate_source_contract(Path(args.candidate))
    if args.profile:
        report["profile"] = validate_profile(json.loads(Path(args.profile).read_text()), args.require_l200_invariants)
    expected = (args.expected_schedule_mode, args.expected_solve_workers,
                args.expected_workers, args.expected_frequency_workers)
    if args.result or (args.parity_candidate_data and args.parity_reference_data):
        if any(value is None for value in expected):
            raise ValueError("result/parity checks require explicit expected schedule and worker counts")
    if args.result:
        report["result"] = validate_result(json.loads(Path(args.result).read_text()), *expected)
    if args.parity_candidate_data and args.parity_reference_data:
        report["parity"] = validate_parity(Path(args.parity_candidate_data), Path(args.parity_reference_data), *expected)
    Path(args.output).write_text(json.dumps(report, indent=2, sort_keys=True)+"\n")
    print(json.dumps(report, sort_keys=True))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
