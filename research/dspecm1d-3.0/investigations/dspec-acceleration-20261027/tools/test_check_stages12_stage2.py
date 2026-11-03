#!/usr/bin/env python3
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import hashlib

MODULE = Path(__file__).with_name("check_stages12_stage2.py")
SPEC = importlib.util.spec_from_file_location("check_stages12_stage2", MODULE)
CHECK = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CHECK)


def result(workers):
    return {"schedule_mode":"frequency", "frequency_workers":workers, "solve_workers":1,
            "workers":1, "kernel_timing_semantics":"kernel_s includes diagnostics and reduction",
            "kernel_parallel_pass_wall_s":1.1, "reduction_wall_s":0.1,
            "diagnostic_worker_elapsed_s":0.3, "kernel_s":2.0, "elapsed_s":4.0,
            "frequency_contribution_buffer_bytes":1000, "tensor_packed_relative_l2":1e-12,
            "time_shift_adjoint_identity":1.0, "gravity_max_abs_nd_error":0.0,
            "scientific":"same", "across_angles":1}


def degree_result():
    value=result(1)
    value.update(schedule_mode="degree", workers=8, solve_workers=8,
        kernel_timing_semantics="legacy degree kernel_s excludes tensor validation",
        kernel_parallel_pass_wall_s=0, reduction_wall_s=0,
        diagnostic_worker_elapsed_s=0, frequency_contribution_buffer_bytes=0)
    return value


def profile():
    rows=[]
    def add(phase,name,domain="master_wall",inc=1.0,exc=1.0,diagnostic=False):
        rows.append({"phase":phase,"scope_name":name,"execution_domain":domain,
            "inclusive_s":inc,"exclusive_s":exc,"calls":1,"logical_component":"test",
            "diagnostic":diagnostic,"parent_overlap_description":"worker elapsed sum overlaps master pass wall" if domain=="worker_elapsed_sum" else "exclusive partition"})
    add("setup","setup",exc=1)
    add("initial_frequency_workspace_and_capture","initial",exc=1)
    add("measurement_and_initial_output","measurement",exc=1)
    add("kernel_frequency_parallel_pass_and_ordered_reduction_including_diagnostics","driver_frequency_parallel_kernel_pass_including_diagnostics",inc=5.1,exc=4.9)
    add("kernel_frequency_parallel_pass_and_ordered_reduction_including_diagnostics","driver_frequency_ordered_kernel_contribution_reduction",exc=.1)
    add("kernel_frequency_parallel_pass_and_ordered_reduction_including_diagnostics","frequency_worker.kernel", "worker_elapsed_sum",5,3)
    add("kernel_frequency_parallel_pass_and_ordered_reduction_including_diagnostics","frequency_worker.diagnostic_child", "worker_elapsed_sum",2,1)
    add("kernel_output_and_result_io","output",exc=1)
    add("diagnostic_tensor_validation","frequency_worker.diagnostic_child", "worker_elapsed_sum",2,1,True)
    add("profile_pipeline_total","profile_pipeline_wall",inc=10,exc=1)
    counters=[]
    for name,value in CHECK.EXPECTED_COUNTERS.items():
        counters.append({"phase":"initial_frequency_workspace_and_capture","name":name,"value":value,"diagnostic":False})
    return {"schema":1,"frequency_route":True,"diagnostic_worker_elapsed_is_not_subtracted_from_wall":True,
        "worker_diagnostic_partition":{"production_rows_exclude_diagnostic_overlay":True,"raw_inclusive_umbrellas_preserve_diagnostics":True,"diagnostic_overlay_phase":"diagnostic_tensor_validation"},
        "wall_partition":{"profile_pipeline_wall_s":10,"master_exclusive_partition_s":10,
            "residual_s":0,"tolerance_s":.01,"worker_elapsed_excluded":True},"rows":rows,"counters":counters}


def write_job(root,job_id,workers,across_angles=1,mode="frequency"):
    job=root/job_id; data=job/"data"; data.mkdir(parents=True)
    payload=degree_result() if mode=="degree" else result(workers); payload["across_angles"]=across_angles
    files={name:("same:"+name).encode() for name in CHECK.SCIENTIFIC_FILES if across_angles or not name.startswith("across")}
    files["result.json"]=json.dumps(payload).encode()
    for name,contents in files.items():(data/name).write_bytes(contents)
    (job/"job.log").write_text("complete\n")
    digest=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
    outputs=["data/"+name for name in files]
    state={"jobs":[{"id":job_id,"status":"completed","returncode":0,"inputs_unchanged":True,
        "cleanup":{"verified_gone":True},"failure_reasons":[],"outputs":outputs,
        "outputs_sha256":{rel:digest(job/rel) for rel in outputs},"log":"job.log","log_sha256":digest(job/"job.log")}]}
    (root/"state.json").write_text(json.dumps(state))


class Stage2CheckerTests(unittest.TestCase):
    def test_frequency_result_and_parity_ignores_scheduler_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); a=root/"a"; b=root/"b"; a.mkdir(); b.mkdir()
            write_job(a,"eight",8)
            write_job(b,"one",1)
            CHECK.validate_result(result(8),"frequency",1,1,8)
            report=CHECK.validate_parity(a/"eight"/"data",b/"one"/"data","frequency",1,1,8)
            self.assertEqual(len(report["byte_exact_scientific_files"]),9)
            self.assertIn("candidate_verified_state_sha256",report)

    def test_build_receipt_binds_each_binary_and_command(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); build=root/"build"; build.mkdir()
            variants=[]; result_variants=[]; supplied={}
            for name in ("native","profile"):
                directory=build/name; directory.mkdir(); binary=directory/"driver"; binary.write_bytes(name.encode())
                digest=hashlib.sha256(binary.read_bytes()).hexdigest(); command_hash=hashlib.sha256((name+" command").encode()).hexdigest()
                (directory/"command.json").write_text(json.dumps({"returncode":0,"tracked_inputs_unchanged":True,"binary_sha256":digest,"command_sha256":command_hash}))
                variants.append({"id":name}); result_variants.append({"id":name,"returncode":0,"binary_sha256":digest,"command_sha256":command_hash}); supplied[name]=binary
            config=root/"config.json"; config.write_text(json.dumps({"build_root":str(build),"variants":variants}))
            result_path=root/"result.json"; result_path.write_text(json.dumps({"status":"completed","build_root":str(build),"variants":result_variants}))
            self.assertEqual(set(CHECK.validate_build(result_path,config,supplied)["variants"]),{"native","profile"})

    def test_rejects_nested_solve_workers(self):
        bad=result(8); bad["solve_workers"]=2
        with self.assertRaisesRegex(ValueError,"solve_workers differs from the plan"):
            CHECK.validate_result(bad,"frequency",1,1,8)

    def test_rejects_wrong_frequency_count(self):
        with self.assertRaisesRegex(ValueError,"differs from the plan"):
            CHECK.validate_result(result(1),"frequency",1,1,8)

    def test_frequency_profile_closes_master_wall_and_separates_workers(self):
        report=CHECK.validate_profile(profile(),require_l200_invariants=True)
        self.assertEqual(report["wall_residual_s"],0)
        self.assertTrue(report["worker_elapsed_excluded"])

    def test_rejects_frequency_worker_mislabeled_master(self):
        bad=profile(); next(r for r in bad["rows"] if r["scope_name"]=="frequency_worker.kernel")["execution_domain"]="master_wall"
        with self.assertRaisesRegex(ValueError,"worker scope/domain mismatch"):
            CHECK.validate_profile(bad)

    def test_rejects_degree_workers_on_frequency_route(self):
        bad=profile(); bad["rows"][5]["scope_name"]="parallel_worker.kernel"
        with self.assertRaisesRegex(ValueError,"worker scope/domain mismatch|degree-worker"):
            CHECK.validate_profile(bad)

    def test_rejects_production_scope_classified_as_pure_diagnostic(self):
        bad=profile()
        row=next(r for r in bad["rows"] if r["phase"]=="initial_frequency_workspace_and_capture")
        row["logical_component"]="diagnostic_validation"
        with self.assertRaisesRegex(ValueError,"production scope was classified as pure diagnostic"):
            CHECK.validate_profile(bad)

    def test_rejects_diagnostic_subtraction_from_wall(self):
        bad=profile(); bad["diagnostic_worker_elapsed_is_not_subtracted_from_wall"]=False
        with self.assertRaisesRegex(ValueError,"must not be subtracted"):
            CHECK.validate_profile(bad)

    def test_rejects_scientific_difference_even_when_schedule_differs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); a=root/"a"; b=root/"b"; a.mkdir(); b.mkdir()
            write_job(a,"eight",8); write_job(b,"one",1)
            (b/"one"/"data"/"across.bin").write_bytes(b"tampered")
            state_path=b/"state.json"; state=json.loads(state_path.read_text()); job=state["jobs"][0]
            job["outputs_sha256"]["data/across.bin"]=hashlib.sha256(b"tampered").hexdigest(); state_path.write_text(json.dumps(state))
            with self.assertRaisesRegex(ValueError,"scientific output differs: across.bin"):
                CHECK.validate_parity(a/"eight"/"data",b/"one"/"data","frequency",1,1,8)

    def test_receipt_verification_rejects_changed_output_hash(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); a=root/"a"; b=root/"b"; a.mkdir(); b.mkdir()
            write_job(a,"eight",8); write_job(b,"one",1)
            (a/"eight"/"data"/"trace.csv").write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError,"hash|checksum|digest|output"):
                CHECK.validate_parity(a/"eight"/"data",b/"one"/"data","frequency",1,1,8)

    def test_parity_accepts_l3_without_across_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); a=root/"a"; b=root/"b"; a.mkdir(); b.mkdir()
            write_job(a,"eight",8,across_angles=0); write_job(b,"one",1,across_angles=0)
            report=CHECK.validate_parity(a/"eight"/"data",b/"one"/"data","frequency",1,1,8)
            self.assertNotIn("across.bin",report["byte_exact_scientific_files"])

    def test_completed_degree_anchor_parity_requires_explicit_degree_contract(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); a=root/"a"; b=root/"b"; a.mkdir(); b.mkdir()
            write_job(a,"anchor",1,mode="degree"); write_job(b,"stage1",1)
            report=CHECK.validate_parity(a/"anchor"/"data",b/"stage1"/"data","degree",8,8,1)
            self.assertEqual(report["candidate_schedule_mode"],"degree")
            self.assertEqual(report["candidate_workers"],8)
            self.assertEqual(len(report["byte_exact_scientific_files"]),9)

    def test_degree_anchor_rejects_wrong_mode_solve_workers_cli_workers_or_frequency_workers(self):
        anchor=degree_result()
        expectations=[("frequency",8,8,1),("degree",1,8,1),
                      ("degree",8,1,1),("degree",8,8,8)]
        for expected in expectations:
            with self.subTest(expected=expected), self.assertRaisesRegex(ValueError,"differs from the plan"):
                CHECK.validate_result(anchor,*expected)

    def test_preserved_degree8_anchor_result_accepts_only_the_degree_contract(self):
        anchor_path=(Path(__file__).resolve().parents[1]/"receipts"/
            "stages12-stage2-v2-measure-failed/raw-task/degree8-anchor/data/result.json")
        self.assertTrue(anchor_path.is_file(),"preserved Stage 2 degree8 anchor receipt is required")
        anchor=json.loads(anchor_path.read_text())
        self.assertEqual(CHECK.validate_result(anchor,"degree",8,8,1)["schedule_mode"],"degree")
        with self.assertRaisesRegex(ValueError,"schedule_mode differs from the plan"):
            CHECK.validate_result(anchor,"frequency",1,1,1)

    def test_frequency_runs_reject_degree_expectation_and_inner_workers(self):
        with self.assertRaisesRegex(ValueError,"schedule_mode differs from the plan"):
            CHECK.validate_result(result(8),"degree",8,8,1)
        bad=result(8);bad["solve_workers"]=2
        with self.assertRaisesRegex(ValueError,"solve_workers differs from the plan"):
            CHECK.validate_result(bad,"frequency",1,1,8)

    def test_rejects_duplicate_diagnostic_worker_row_in_production_phase(self):
        bad=profile()
        duplicate=dict(next(r for r in bad["rows"] if r["scope_name"]=="frequency_worker.diagnostic_child"))
        duplicate["phase"]="kernel_frequency_parallel_pass_and_ordered_reduction_including_diagnostics"
        duplicate["diagnostic"]=True
        bad["rows"].append(duplicate)
        with self.assertRaisesRegex(ValueError,"diagnostic row/phase mismatch"):
            CHECK.validate_profile(bad)


if __name__ == "__main__":
    unittest.main()
