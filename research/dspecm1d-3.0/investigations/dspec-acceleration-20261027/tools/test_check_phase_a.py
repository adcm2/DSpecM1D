import json
import math
import struct
import tempfile
import unittest
from pathlib import Path
import csv

import check_phase_a as checker


class PhaseACheckerTests(unittest.TestCase):
    A5_SCOPES={
        "phasea_standard_radial_static_operator_and_receiver_setup",
        "phasea_standard_radial_parallel_region_wall","phasea_standard_radial_worker_region",
        "phasea_standard_radial_frequency_source_treatment","phasea_standard_radial_dynamic_operator_assembly",
        "phasea_standard_radial_receiver_lift_and_projection",
        "phasea_standard_toroidal_parallel_region_wall","phasea_standard_toroidal_worker_region",
        "phasea_standard_toroidal_degree_static_operators_and_coefficients",
        "phasea_standard_toroidal_cutoff_policy_selection","phasea_standard_toroidal_frequency_source_treatment",
        "phasea_standard_toroidal_dynamic_operator_assembly","phasea_standard_toroidal_receiver_lift_and_projection",
        "phasea_standard_spheroidal_cowling_cutoff_boundary_selection",
        "phasea_standard_spheroidal_parallel_region_wall","phasea_standard_spheroidal_worker_region",
        "phasea_standard_spheroidal_degree_angular_and_source_coefficients",
        "phasea_standard_spheroidal_formulation_static_operators_and_receiver_basis",
        "phasea_standard_spheroidal_frequency_start_index_policy_selection",
        "phasea_standard_spheroidal_frequency_source_treatment",
        "phasea_standard_spheroidal_dynamic_operator_assembly",
        "phasea_standard_spheroidal_receiver_lift_and_projection",
    }
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.root=Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def write(self,name,data):
        path=self.root/name
        path.write_text(json.dumps(data,allow_nan=True))
        return path

    def test_a4_requires_real_forward_work_and_zero_adjoint_routes(self):
        self.write("a4-summary.json",{"bins":35,"workers":8,"forward_only_route_calls":35,
            "receiver_relative_l2":0.0,"unit_hh_rhs_columns":0,"weighted_hh_rhs_columns":0,
            "receiver_qr_scopes":0})
        profile={"scope_rows":[],"counters":[{"name":"true_forward_only_route","value":35},
            {"name":"solver_forward_rhs_columns","value":10},
            {"name":"solver_unit_basis_rhs_columns","value":0},
            {"name":"solver_weighted_adjoint_rhs_columns","value":0}]}
        self.write("a4-profile.json",profile)
        (self.root/"a4-components.csv").write_text(
            "bin,up_re,up_im,north_re,north_im,east_re,east_im\n"+
            "".join(f"{b},1,0,0,1,1,1\n" for b in range(7,42)))
        self.assertEqual(checker.data_gate("a4",self.root)["status"],"pass")
        profile["counters"].append({"name":"full_weighted_modal_adjoint_rhs_columns","value":1})
        self.write("a4-profile.json",profile)
        with self.assertRaisesRegex(ValueError,"profiler call/RHS counters"):
            checker.data_gate("a4",self.root)
        profile["counters"].pop()
        profile["counters"][1]["value"]=0
        self.write("a4-profile.json",profile)
        with self.assertRaisesRegex(ValueError,"profiler call/RHS counters"):
            checker.data_gate("a4",self.root)
        profile["counters"][1]["value"]=10
        self.write("a4-profile.json",profile)
        (self.root/"a4-components.csv").write_text(
            "bin,up_re,up_im,north_re,north_im,east_re,east_im\n7,1,0,0,1,1,1\n")
        with self.assertRaisesRegex(ValueError,"exactly 35 bins"):
            checker.data_gate("a4",self.root)

    def test_a4_rejects_nan_even_when_ordering_gate_would_not(self):
        self.write("a4-summary.json",{"bins":35,"workers":8,"forward_only_route_calls":35,
            "receiver_relative_l2":math.nan,"unit_hh_rhs_columns":0,"weighted_hh_rhs_columns":0,
            "receiver_qr_scopes":0})
        self.write("a4-profile.json",{"scope_rows":[],"counters":[]})
        with self.assertRaisesRegex(ValueError,"non-finite"):
            checker.data_gate("a4",self.root)

    def test_packet_schema_checks_padded_bytes_and_every_radial_value(self):
        header=struct.pack("<8sIIQiI4d",b"DSPHAR1\0",1,1000,1,9,0,1.0,0.1,6.0e6,500.0)
        record=bytearray(1000)
        path=self.root/"radial-bin9.bin"
        path.write_bytes(header+record)
        self.assertEqual(checker.packet(path,9,1)["record_bytes"],1000)
        struct.pack_into("<d",record,232,float("nan"))
        path.write_bytes(header+record)
        with self.assertRaisesRegex(ValueError,"non-finite field"):
            checker.packet(path,9,1)

    def test_a3_requires_fd_and_frozen_input_identity(self):
        summary={"degree":3,"nq":8,"refinement":0.05,"frequency_count":35,
            "radial_sides":157,"angular_points":91,"fd_relative_error":1e-8,
            "gravity_max_abs_nd_error":0.0,"time_shift_adjoint_identity":1.0,
            "time_shift_adjoint_identity_source":"frozen L3 baseline result.json",
            "source_integrity":"same frozen +/- traces and same L3 model/mesh/weights",
            "common_vs_frozen_ds_relative":1e-9}
        self.write("a3-summary.json",summary)
        self.assertEqual(checker.data_gate("a3",self.root)["status"],"pass")
        summary["fd_relative_error"]=2.1e-4
        self.write("a3-summary.json",summary)
        with self.assertRaisesRegex(ValueError,"finite-difference"):
            checker.data_gate("a3",self.root)

    def test_a2_uses_approved_all_side_aggregate_groups(self):
        groups=[]
        for b in (9,25,40):
            for s in (0,1):
                for role in (0,1):
                    groups.append({"bin":b,"section":s,"role":role,
                        "records":(135 if s==0 else 19)*293,"relative_l2":2e-12})
        diagnostic_groups=[{"bin":9,"section":0,"role":0,"fluid":1,
            "records":1,"relative_l2":0.5}]
        self.write("a2-summary.json",{"bins":[9,25,40],"reference_records":3*2*293*154,
            "field_gate":1e-11,"max_reference_vs_common_field_l2":2e-12,
            "max_same_packet_ds_vs_common_field_l2":2e-12,
            "approved_field_groups":groups,"approved_ds_common_field_groups":groups,
            "field_groups":diagnostic_groups,"ds_common_field_groups":diagnostic_groups,
            "kernel_products":[{} for _ in range(6)],"radial_solves":0,
            "old_reference_route":"same-workspace PACKED oracle, not Stage 2 TENSOR timing"})
        (self.root/"a2-field-comparison.csv").write_text("x\n")
        (self.root/"a2-kernel-products.csv").write_text("x\n")
        self.assertEqual(checker.data_gate("a2",self.root)["status"],"pass")
        broken=list(groups);broken.pop()
        summary=json.loads((self.root/"a2-summary.json").read_text())
        summary["approved_field_groups"]=broken
        self.write("a2-summary.json",summary)
        with self.assertRaisesRegex(ValueError,"approved aggregate groups are incomplete"):
            checker.data_gate("a2",self.root)

    def test_a5_requires_default_cadence_and_exact_counter_set(self):
        summary={"bins":35,"workers":8,"worker_domain":"standard_degree_parallel",
            "nskip":10,"prepared_single_sem_policy":True,
            "untruncated_factorizations":14035,"untruncated_forward_rhs_columns":42035,
            "untruncated_adjoint_rhs_columns":0,"truncated_factorizations":14035,
            "truncated_forward_rhs_columns":42035,"truncated_adjoint_rhs_columns":0,
            "untruncated_receiver_relative_l2":1e-12,"full_pipeline_speed_claim":False,
            "same_mesh_and_forcing":True,
            "instrumented_elapsed_includes_cut_index_capture_and_allocation":True,
            "instrumented_elapsed_includes_cut_index_csv_serialization":False,
            "untruncated_cut_index_rows":14035,"truncated_cut_index_rows":14035}
        self.write("a5-summary.json",summary)
        header="family,degree,bin,formulation,start_dof,retained_rows,full_rows,nskip,policy\n"
        for mode,policy,start in (("untruncated","untruncated",0),("truncated","source_aware",1)):
            rows=[]
            for b in range(7,42):
                rows.append(f"radial,0,{b},standard,{start},{10-start},10,10,{policy}")
                for degree in range(1,201):
                    rows.append(f"toroidal,{degree},{b},standard,{start},{10-start},10,10,{policy}")
                    rows.append(f"spheroidal,{degree},{b},full,{start},{10-start},10,10,{policy}")
            (self.root/f"a5-{mode}-cut-indices.csv").write_text(header+"\n".join(rows)+"\n")
            (self.root/f"a5-{mode}-components.csv").write_text(
                "bin,up_re,up_im,north_re,north_im,east_re,east_im\n"+
                "".join(f"{b},1,0,0,1,1,1\n" for b in range(7,42)))
        for name in ("a5-untruncated-profile.json","a5-truncated-profile.json"):
            rows=[]
            for n in self.A5_SCOPES:
                worker=(n.endswith("worker_region") or any(t in n for t in
                    ("degree_static_","frequency_source_treatment","dynamic_operator_assembly",
                     "receiver_lift_and_projection","cutoff_policy_selection",
                     "start_index_policy_selection","formulation_static_operators")))
                domain="worker_elapsed_sum" if worker else "master_wall"
                rows.append({"name":n,"execution_domain":domain})
            rows.append({"name":"standard_worker.phasea_standard_radial_worker_region","execution_domain":"worker_elapsed_sum"})
            self.write(name,{"scope_rows":rows,"counters":[
                {"name":"sparse_factorization_calls","value":14035},
                {"name":"forward_rhs_columns","value":42035},
                {"name":"adjoint_rhs_columns","value":0}]})
        self.assertEqual(checker.data_gate("a5",self.root)["status"],"pass")
        raw_path=self.root/"a5-untruncated-profile.json"
        raw=json.loads(raw_path.read_text())
        raw["counters"][0]["value"]=14034
        self.write("a5-untruncated-profile.json",raw)
        with self.assertRaisesRegex(ValueError,"raw/summary counter mismatch"):
            checker.data_gate("a5",self.root)
        summary["nskip"]=1
        self.write("a5-summary.json",summary)
        with self.assertRaisesRegex(ValueError,"policy binding"):
            checker.data_gate("a5",self.root)

    def test_a5_rejects_duplicate_raw_counter_names(self):
        summary={"bins":35,"workers":8,"worker_domain":"standard_degree_parallel",
            "nskip":10,"prepared_single_sem_policy":True,
            "untruncated_factorizations":14035,"untruncated_forward_rhs_columns":42035,
            "untruncated_adjoint_rhs_columns":0,"truncated_factorizations":14035,
            "truncated_forward_rhs_columns":42035,"truncated_adjoint_rhs_columns":0,
            "untruncated_receiver_relative_l2":0.0,"full_pipeline_speed_claim":False,
            "same_mesh_and_forcing":True,
            "instrumented_elapsed_includes_cut_index_capture_and_allocation":True,
            "instrumented_elapsed_includes_cut_index_csv_serialization":False,
            "untruncated_cut_index_rows":14035,"truncated_cut_index_rows":14035}
        self.write("a5-summary.json",summary)
        for mode,policy in (("untruncated","untruncated"),("truncated","source_aware")):
            rows=[]
            for b in range(7,42):
                rows.append(f"radial,0,{b},standard,0,10,10,10,{policy}")
                for degree in range(1,201):
                    rows.append(f"toroidal,{degree},{b},standard,0,10,10,10,{policy}")
                    rows.append(f"spheroidal,{degree},{b},full,0,10,10,10,{policy}")
            (self.root/f"a5-{mode}-cut-indices.csv").write_text(
                "family,degree,bin,formulation,start_dof,retained_rows,full_rows,nskip,policy\n"+"\n".join(rows)+"\n")
            (self.root/f"a5-{mode}-components.csv").write_text(
                "bin,up_re,up_im,north_re,north_im,east_re,east_im\n"+
                "".join(f"{b},1,0,0,1,1,1\n" for b in range(7,42)))
            scope_names=sorted(self.A5_SCOPES)
            scope_rows=[]
            for n in scope_names:
                worker=(n.endswith("worker_region") or "parallel_region_wall" not in n)
                scope_rows.append({"name":n,"execution_domain":"worker_elapsed_sum" if worker else "master_wall"})
            self.write(f"a5-{mode}-profile.json",{"scope_rows":scope_rows,"counters":[
                {"name":"sparse_factorization_calls","value":14035},
                {"name":"forward_rhs_columns","value":42035},
                {"name":"adjoint_rhs_columns","value":0},
                {"name":"forward_rhs_columns","value":42035}]})
        with self.assertRaisesRegex(ValueError,"duplicate raw counter"):
            checker.data_gate("a5",self.root)


if __name__=="__main__":
    unittest.main()
