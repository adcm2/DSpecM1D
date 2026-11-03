#!/usr/bin/env python3
"""Receipt-bound Phase A source/build/output gates; never launches calculations."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

SEED = {
    "additional_dependencies_sha256.json":"9c8ca72c9f40fd6b2f18aab27c176f649089f575a3e91cd6d008421906a61a2f",
    "case/kernels-small.json":"7b0efae842bff4b414b89eaecf209a986d6e918e0641e147cea0881d46009cc0",
    "case/prem.200":"4c3c53ca04e785ff4fd853cf211fe94e353ef70d6052a526122d0f279a3519ea",
    "doc_driver.cpp":"0014948bbaad9cc445fb6afa13c4bb3b6f4ea1d3e53b2b65f90e99d68b09c981",
    "engine_headers/SensitivityKernels.h":"2fa57eb58a5d937aaae592a0711999f84fae505a12e1936d0958d995850e6b8c",
    "engine_headers/Windowed3DFields.h":"4d9f3d2b4f6fb6540b89a0bc74d7f6d1ee2e47f2cb2bcaa21b36624a6ba0f432",
    "engine_headers/Windowed3DKernels.h":"a852dead24ba03c1c4139ef926b807a88551ce2c2d77f443ddf149cf150d5506",
    "engine_headers/profile.hpp":"2aff92d25d28434fc07302b73397c790a270ab76b5c858e1309785f8a4a09bc1",
    "fast_engine.hpp":"3bb252b3e94f3365f02a5dd1822030ad27a0ef764210f3311381ed50386ca73f",
    "frequency_scheduler.hpp":"d954acaf93f2334caa655caa82559e71251619c99339a08a5220536992d9fd60",
    "frequency_scheduler_test.cpp":"0d7a422c5cd84000e82fc76b227cc5fc35dde7f0866a8214f1cff6a8ba52b714",
    "global/global_gsh_tile.hpp":"0a884177fc688b6a78d95c04209571b90329124f7d7437af5a063887484671e0",
    "global/global_kernel_consumer.hpp":"cfe88cf6b02800152bd3fe6f500e718392b66cc4e9533ee27c4157b2027741f6",
    "global/pole_response_spool.hpp":"5e8868439771cb502cb4e3005039f6986e9d0d6ad23a58b168a7a61c2c9152dc",
    "implementation_status.md":"48383ae3b40c0adbb1fc259d0af1dbe8d14f11d86ac6fef81476b6d41e2f337f",
    "matched_mesh_model.hpp":"75ab1363d6914bb9e0f893da694dd92e728004018e2786112fbff30c835329cf",
    "physical_driver.cpp":"fbc6f8ffd740ca2a9d727de0317567c8d105e813d31e795341475116d4fbb7e2",
    "profile.hpp":"8e1ddc414cd200d9858581622c3cf9a04b26ddca26aeff673b718a75e48a191d",
    "profile_equilibrated_sparse_lu.hpp":"cad08d7bad6e71a87f9205f0e533cc6ee2106dbe766497c7992baf1e50cd8324",
    "solvers/batch_response_cache.hpp":"576c80b7d342fbe1f9e41fe667e6a1db94d59825c85364cc0cf06e806ccb7750",
    "solvers/solver_backends.hpp":"0a55afc975d6a1fdb8a26689a872530b8ee8cd736c044e07e5c990d526b758c4",
    "stage-to-math-map.json":"13d447b87568f14f311be1611f101b0ef8ed5d47a2feb4fb22f3720fdbc5bb02",
}
ADDED = {
    "phase_a_contract.hpp", "phase_a_contract_test.cpp", "phase_a_driver.cpp",
    "phase_a_reconstruction.hpp", "phase_a_standard_forward.hpp",
    "provenance/pinned_formula_sources/README.md",
    "provenance/pinned_formula_sources/FullSpecMultiSem.h",
    "provenance/pinned_formula_sources/FullSpecSingleSem.h",
    "provenance/pinned_formula_sources/SampledKernels.cpp",
    "provenance/pinned_formula_sources/StartElement.h",
    "provenance/pinned_formula_sources/source-manifest.json",
}
EXPECTED_DEPENDENCY_TREES = {
    "/space/adcm2/research/projects/DSpecM1D/DSpecM1D":"06bf110dfe7691c5f311a933157523e5b4546f34134038307a7917c733068ea8",
    "/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya/bootstrap/dependencies/numericconcepts/include":"cd8f874d886a47cb89864fa739392936076b8aafaabe2dd82406ad976e5d0096",
    "/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya/bootstrap/dependencies/fftwpp":"d6ed41f835aa44edd7952e2c2e8f1fbb80f49e4cefd7b93277deb2cace8489b5",
    "/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya/bootstrap/dependencies/gaussquad":"b554b36e9832c65e2c1c43eadfaf4d5e48add121db43bcc5b5926a156636af8f",
    "/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya/bootstrap/dependencies/interpolation/include":"0ea2c712d7a0fa2ac169272b73cb0593ef04c2f314a29de600622237556b707f",
    "/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya/bootstrap/dependencies/gshtrans":"04ecaaafa43b03022ba1c7ca277fc83c3804bbad2409cc22ed54085860386417",
    "/space/adcm2/research/runs/dspec-yspec-kernels/2026-10-27-tunya/bootstrap/dependencies/planetmodel-dspec/cpp/include":"363c310d71e1e4ed19daa0118aa9c07b49aab58a125c409117b6691ae4368abf",
    "/space/adcm2/research/projects/CCDSM/deps/toolchain/eigen/include/eigen3":"6d67d4cbc351e6f66d12298110d0e63188bb3b1b66d70c1c27a37390ef8809ba",
    "/space/adcm2/research/projects/YSpecpp/src":"26f29edb53e94fda4859f3bb962b49b33cd3c8da8e1d3433a3037a6c7608e692",
    "/space/adcm2/research/projects/YSpecpp/YSpec":"eb023ebfd0953a9b60d5be3b4d7c652ff320f573755b4cbd4913ca3118615853",
    "/space/adcm2/research/builds/documentation-wedges-20261027-01/yspec-library/libyspec-doc.so":"f891e44ed9d9ca8aaa871fcbb3d9287526aee90dc33bf4e58867abc14e9a6d8f",
}
EXPECTED_CANDIDATE_TREE = "4c94847f0a46ed96fe166165cc993ee349e9ec0ad0531fd564f89bd58bd49b4e"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def finite_tree(value, where="root"):
    if isinstance(value, bool) or value is None or isinstance(value, str):
        return
    if isinstance(value, (int, float)):
        require(math.isfinite(float(value)), f"non-finite numeric value at {where}")
    elif isinstance(value, list):
        for i, child in enumerate(value): finite_tree(child, f"{where}[{i}]")
    elif isinstance(value, dict):
        for key, child in value.items(): finite_tree(child, f"{where}.{key}")


def read_json(path: Path):
    data=json.loads(path.read_text())
    finite_tree(data,str(path))
    return data


def candidate_gate(root: Path):
    actual={str(p.relative_to(root)) for p in root.rglob("*") if p.is_file()}
    allowed=set(SEED)|ADDED
    require(actual==allowed, f"candidate file inventory differs; missing={sorted(allowed-actual)}, extra={sorted(actual-allowed)}")
    mismatches=[]
    for rel,digest in SEED.items():
        if rel in {"implementation_status.md","engine_headers/Windowed3DFields.h"}: continue
        if sha(root/rel)!=digest: mismatches.append(rel)
    require(not mismatches, f"unapproved changes to frozen seed files: {mismatches}")
    hook=(root/"engine_headers/Windowed3DFields.h").read_text()
    require("visitProjectedRadialModes" in hook and "phaseAngularCoordinates" in hook,
            "read-only projected radial export hook is missing")
    return {"status":"pass","seed_file_count":len(SEED),"allowed_new_files":len(ADDED),
            "modified_seed_files":["engine_headers/Windowed3DFields.h","implementation_status.md"]}


def build_gate(result_path: Path, config_path: Path, build_root: Path):
    result=read_json(result_path);config=read_json(config_path)
    require(result.get("status")=="completed","three-variant build did not complete")
    variants={v["id"]:v for v in result.get("variants",[])}
    require(set(variants)=={"phase-a-native","phase-a-profile","phase-a-tiny"},"build result variant inventory mismatch")
    rows=[];native_command=profile_command=None;tracked_sources={}
    configured={v["id"]:v for v in config.get("variants",[])}
    for name in ("phase-a-native","phase-a-profile","phase-a-tiny"):
        record=variants[name]
        require(record.get("returncode")==0 and record.get("binary_sha256"),f"{name} failed or has no binary digest")
        binary=build_root/name/ ("phase-a-contract-test" if name.endswith("tiny") else "phase-a-driver")
        command_path=binary.parent/"command.json"
        command=read_json(command_path)
        require(command.get("returncode")==0 and command.get("tracked_inputs_unchanged") is True,
                f"{name} failed command/input stability gate")
        tracked=command.get("tracked_input_sha256",{});after=command.get("tracked_input_sha256_after",{})
        source_path=configured[name]["sources"][0]
        require(source_path in tracked and tracked[source_path]==after.get(source_path),
                f"{name} source tree changed during build or is not recorded")
        for path,digest in EXPECTED_DEPENDENCY_TREES.items():
            require(tracked.get(path)==digest and after.get(path)==digest,
                    f"{name} dependency tree differs from accepted preflight: {path}")
        tracked_sources[name]=tracked[source_path]
        require(sha(binary)==record["binary_sha256"]==command.get("binary_sha256"),f"{name} binary digest mismatch")
        if name=="phase-a-native":
            require("-DMATHPROFILE_DISABLED" in command["command"],"native build does not disable profiling")
            native_command=command["command"]
        if name=="phase-a-profile":
            require("-DMATHPROFILE_DISABLED" not in command["command"],"profile build has profiling disabled")
            profile_command=command["command"]
        rows.append({"id":name,"binary_sha256":record["binary_sha256"],
                     "command_sha256":command.get("command_sha256"),
                     "tracked_inputs_unchanged":True})
    require(native_command is not None and profile_command is not None,"native/profile command records missing")
    require(len(set(tracked_sources.values()))==1,"native/profile/tiny builds did not consume one frozen candidate tree")
    require(set(tracked_sources.values())=={EXPECTED_CANDIDATE_TREE},
            "build variants did not consume the approved staged candidate tree")
    def norm(command):
        out=[];i=0
        while i<len(command):
            if command[i]=="-DMATHPROFILE_DISABLED": i+=1;continue
            if command[i]=="-o": out.extend(["-o","{OUTPUT}"]);i+=2;continue
            out.append(command[i]);i+=1
        return out
    require(norm(native_command)==norm(profile_command),"native/profile commands differ beyond MATHPROFILE_DISABLED")
    return {"status":"pass","variants":rows,"candidate_source_tree_sha256":tracked_sources["phase-a-native"],
            "dependency_trees_match_accepted_preflight":True,
            "native_profile_difference":"MATHPROFILE_DISABLED only"}


def packet(path: Path, expected_bin: int, expected_count: int):
    b=path.read_bytes()
    require(len(b)>=64,"packet shorter than fixed header")
    magic,schema,record_bytes,count,bin_number,reserved,zr,zi,length,time=__import__('struct').unpack("<8sIIQiI4d",b[:64])
    require(magic==b"DSPHAR1\0" and schema==1 and record_bytes>0,"packet header magic/schema/record size mismatch")
    require(bin_number==expected_bin and count==expected_count,"packet bin or record count mismatch")
    expected_record_bytes=1000 if path.name.startswith("radial-") else 160
    require(record_bytes==expected_record_bytes,"packet padded record size differs from the frozen schema")
    require(math.isfinite(zr) and math.isfinite(zi) and math.isfinite(length) and math.isfinite(time),"packet header has non-finite scales/frequency")
    require(length>0 and time>0 and len(b)==64+record_bytes*count,"packet size/count binding failed")
    import struct
    double_offsets=(range(32,72,8),range(72,216,8),range(232,1000,8)) if expected_record_bytes==1000 else (range(40,160,8),)
    for record_index in range(count):
        base=64+record_index*record_bytes
        for offsets in double_offsets:
            for offset in offsets:
                value=struct.unpack_from("<d",b,base+offset)[0]
                require(math.isfinite(value),f"packet contains a non-finite field at record {record_index}")
    return {"bin":bin_number,"record_count":count,"record_bytes":record_bytes,"sha256":sha(path)}


def finite_csv(path: Path, skip_columns=()):
    with path.open(newline="") as stream:
        reader=csv.DictReader(stream)
        require(reader.fieldnames is not None,"CSV header missing")
        for line,row in enumerate(reader,2):
            for name,value in row.items():
                if name in skip_columns: continue
                require(value not in (None,""),f"empty CSV value at {path}:{line}:{name}")
                try: number=float(value)
                except ValueError: raise ValueError(f"non-numeric CSV field at {path}:{line}:{name}")
                require(math.isfinite(number),f"non-finite CSV value at {path}:{line}:{name}")


def approved_a2_groups(d, field):
    groups=d.get(field)
    require(isinstance(groups,list) and len(groups)==12,
            f"A2 approved aggregate groups are incomplete: {field}")
    observed={}
    for row in groups:
        key=(row.get("bin"),row.get("section"),row.get("role"))
        require(key not in observed,"A2 duplicate approved aggregate group")
        expected=(19 if key[1]==1 else 135)*293 if key[1] in (0,1) else -1
        require(key[0] in (9,25,40) and key[1] in (0,1) and key[2] in (0,1) and
                row.get("records")==expected,
                "A2 approved aggregate group coverage mismatch")
        value=row.get("relative_l2")
        require(isinstance(value,(int,float)) and math.isfinite(value) and value>=0,
                "A2 approved aggregate relative error is missing or non-finite")
        observed[key]=value
    require(set(observed)=={(b,s,r) for b in (9,25,40) for s in (0,1) for r in (0,1)},
            "A2 approved aggregate group keys are incomplete")
    return max(observed.values())


def component_csv(path: Path):
    expected_header=["bin","up_re","up_im","north_re","north_im","east_re","east_im"]
    with path.open(newline="") as stream:
        reader=csv.DictReader(stream)
        require(reader.fieldnames==expected_header,f"component CSV schema mismatch: {path}")
        rows=list(reader)
    require(len(rows)==35,f"component CSV must contain exactly 35 bins: {path}")
    bins=[]
    for row in rows:
        require(None not in row and all(value is not None for value in row.values()),
                f"component CSV has an extra or missing column: {path}")
        try: bins.append(int(row["bin"]))
        except ValueError: raise ValueError(f"component CSV has a non-integer bin: {path}")
        for name in expected_header[1:]:
            try: value=float(row[name])
            except ValueError: raise ValueError(f"component CSV has a non-numeric field: {path}:{name}")
            require(math.isfinite(value),f"component CSV has a non-finite field: {path}:{name}")
    require(bins==list(range(7,42)),f"component CSV bin inventory/order mismatch: {path}")
    return {"bins":len(rows),"sha256":sha(path)}


def profile_counter_map(profile, label):
    counters=profile.get("counters")
    require(isinstance(counters,list),f"{label} raw counters missing")
    result={}
    for row in counters:
        name=row.get("name")
        require(isinstance(name,str) and name not in result,
                f"{label} has a missing or duplicate raw counter name")
        value=row.get("value")
        require(isinstance(value,int) and not isinstance(value,bool) and value>=0,
                f"{label} has an invalid raw counter value for {name}")
        result[name]=value
    return result


def data_gate(mode: str, root: Path):
    if mode=="a1":
        rows=[]
        for bin_number in (9,25,40):
            rows.append(packet(root/f"radial-bin{bin_number}.bin",bin_number,201*293))
            rows.append(packet(root/f"reference-bin{bin_number}.bin",bin_number,2*293*(135+19)))
        metadata=read_json(root/"packet-metadata.json")
        require(metadata.get("frequency_solves")==3 and metadata.get("native_sides")==293 and
                metadata.get("along_points")==135 and metadata.get("across_points")==19,
                "A1 packet metadata inventory mismatch")
        finite_csv(root/"active_bins.csv")
        return {"status":"pass","packets":rows}
    if mode=="a2":
        d=read_json(root/"a2-summary.json")
        accepted_reference=approved_a2_groups(d,"approved_field_groups")
        accepted_ds=approved_a2_groups(d,"approved_ds_common_field_groups")
        require(d.get("bins")==[9,25,40] and d.get("reference_records")==3*2*293*154,
                "A2 coverage count mismatch")
        require(d.get("field_gate")==1e-11 and
                d.get("max_reference_vs_common_field_l2")==accepted_reference and accepted_reference<=1e-11 and
                d.get("max_same_packet_ds_vs_common_field_l2")==accepted_ds and accepted_ds<=1e-11,
                "A2 complex field parity gate failed")
        require(len(d.get("field_groups",[]))>0 and len(d.get("ds_common_field_groups",[]))>0 and
                len(d.get("kernel_products",[]))==6,"A2 field/kernel groups are incomplete")
        require(d.get("radial_solves")==0 and "not Stage 2 TENSOR timing" in d.get("old_reference_route",""),
                "A2 replay scope is mislabeled")
        finite_csv(root/"a2-field-comparison.csv")
        finite_csv(root/"a2-kernel-products.csv")
        return {"status":"pass","max_reference_vs_common_field_l2":accepted_reference,
                "max_same_packet_ds_vs_common_field_l2":accepted_ds,
                "kernel_groups":len(d["kernel_products"])}
    if mode=="a3":
        d=read_json(root/"a3-summary.json")
        require(d.get("degree")==3 and d.get("nq")==8 and d.get("refinement")==0.05 and
                d.get("frequency_count")==35 and d.get("radial_sides")==157 and
                d.get("angular_points")==91,"A3 mesh/quadrature coverage mismatch")
        require(d.get("fd_relative_error",float("inf"))<=2e-4 and
                d.get("gravity_max_abs_nd_error",float("inf"))<=1e-10 and
                abs(d.get("time_shift_adjoint_identity",0)-1)<=1e-10,
                "A3 finite-difference/gravity/time-identity gate failed")
        require(d.get("source_integrity")=="same frozen +/- traces and same L3 model/mesh/weights",
                "A3 frozen-input qualification label missing")
        require(d.get("time_shift_adjoint_identity_source")=="frozen L3 baseline result.json",
                "A3 time-shift identity is not identified as frozen evidence")
        return {"status":"pass","fd_relative_error":d["fd_relative_error"],
                "common_vs_frozen_ds_relative":d["common_vs_frozen_ds_relative"]}
    if mode=="a4":
        d=read_json(root/"a4-summary.json");p=read_json(root/"a4-profile.json")
        component_csv(root/"a4-components.csv")
        require(d.get("bins")==35 and d.get("workers")==8 and d.get("forward_only_route_calls")==35,
                "A4 forward-only route coverage mismatch")
        require(d.get("receiver_relative_l2",float("inf"))<=1e-11 and
                d.get("unit_hh_rhs_columns")==0 and d.get("weighted_hh_rhs_columns")==0 and
                d.get("receiver_qr_scopes")==0,"A4 parity/no-adjoint gate failed")
        counters=profile_counter_map(p,"A4")
        require(counters.get("true_forward_only_route")==35 and counters.get("solver_forward_rhs_columns",0)>0 and
                counters.get("solver_unit_basis_rhs_columns")==0 and
                counters.get("solver_weighted_adjoint_rhs_columns")==0 and
                counters.get("full_weighted_modal_adjoint_rhs_columns",0)==0,
                "A4 profiler call/RHS counters failed")
        require(not any(any(x in r["name"] for x in ("minimal_unit_basis_HH_solve","full_weighted_HH_solve","receiver_basis_and_QR_reconstruction","receiver_basis_assembly_and_qr")) for r in p.get("scope_rows",[])),
                "A4 entered unit/weighted HH or receiver QR scopes")
        return {"status":"pass","relative_l2":d["receiver_relative_l2"],"counters":counters}
    if mode=="a5":
        d=read_json(root/"a5-summary.json")
        require(d.get("bins")==35 and d.get("workers")==8 and d.get("worker_domain")=="standard_degree_parallel" and
                d.get("nskip")==10 and d.get("prepared_single_sem_policy") is True,
                "A5 standard prepared-SEM policy binding mismatch")
        for mode_name in ("untruncated","truncated"):
            require(d.get(mode_name+"_factorizations")==14035 and
                    d.get(mode_name+"_forward_rhs_columns")==42035 and
                    d.get(mode_name+"_adjoint_rhs_columns")==0,
                    f"A5 {mode_name} counters mismatch")
        require(d.get("untruncated_receiver_relative_l2",float("inf"))<=1e-11 and
                d.get("full_pipeline_speed_claim") is False and d.get("same_mesh_and_forcing") is True and
                d.get("instrumented_elapsed_includes_cut_index_capture_and_allocation") is True and
                d.get("instrumented_elapsed_includes_cut_index_csv_serialization") is False,
                "A5 matched receiver/scope gate failed")
        cut_rows={}
        profiles=[read_json(root/"a5-untruncated-profile.json"),read_json(root/"a5-truncated-profile.json")]
        for mode_name,policy in (("untruncated","untruncated"),("truncated","source_aware")):
            component_csv(root/f"a5-{mode_name}-components.csv")
            with (root/f"a5-{mode_name}-cut-indices.csv").open(newline="") as stream:
                rows=list(csv.DictReader(stream))
            expected_keys=set()
            for bin_number in range(7,42):
                expected_keys.add(("radial",0,bin_number,"standard"))
                for degree in range(1,201):
                    expected_keys.add(("toroidal",degree,bin_number,"standard"))
                    expected_keys.add(("spheroidal",degree,bin_number,"full"))
            require(len(rows)==14035 and len(rows)==d.get(mode_name+"_cut_index_rows"),
                    f"A5 {mode_name} cut-index row count mismatch")
            seen_families=set()
            seen_keys=set()
            for row in rows:
                start=int(row["start_dof"]);retained=int(row["retained_rows"]);full=int(row["full_rows"])
                require(row["policy"]==policy and int(row["nskip"])==10 and start>=0 and full>0 and
                        retained==full-start and retained>0,f"A5 {mode_name} invalid selected cut record")
                if mode_name=="untruncated":require(start==0,f"A5 {mode_name} start DOF is not zero")
                seen_families.add(row["family"])
                key=(row["family"],int(row["degree"]),int(row["bin"]),row["formulation"])
                require(key not in seen_keys,f"A5 {mode_name} duplicate cut-index key")
                seen_keys.add(key)
            require(seen_families=={"radial","toroidal","spheroidal"},f"A5 {mode_name} family cut records incomplete")
            require(seen_keys==expected_keys,f"A5 {mode_name} cut-index key inventory mismatch")
            cut_rows[mode_name]=len(rows)
        require(all(p.get("scope_rows") and p.get("counters") for p in profiles),"A5 worker profile evidence missing")
        for mode_name,profile in zip(("untruncated","truncated"),profiles):
            raw=profile_counter_map(profile,f"A5 {mode_name}")
            expected={"sparse_factorization_calls":14035,"forward_rhs_columns":42035,"adjoint_rhs_columns":0}
            for counter,value in expected.items():
                actual=raw.get(counter,0)
                require(actual==value and d.get(f"{mode_name}_"+({"sparse_factorization_calls":"factorizations","forward_rhs_columns":"forward_rhs_columns","adjoint_rhs_columns":"adjoint_rhs_columns"}[counter]))==actual,
                        f"A5 {mode_name} raw/summary counter mismatch for {counter}")
        required_scopes={
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
        for profile in profiles:
            names={row["name"].removeprefix("standard_worker.") for row in profile["scope_rows"]}
            require(required_scopes<=names,"A5 named component/worker/wall timing coverage incomplete")
            for row in profile["scope_rows"]:
                local=row["name"].removeprefix("standard_worker.")
                expected_domain=("worker_elapsed_sum" if local.endswith("worker_region") or
                    any(token in local for token in ("degree_static_","frequency_source_treatment",
                        "dynamic_operator_assembly","receiver_lift_and_projection",
                        "cutoff_policy_selection","start_index_policy_selection",
                        "formulation_static_operators")) else None)
                if "parallel_region_wall" in local or "cowling_cutoff_boundary_selection" in local:
                    expected_domain="master_wall"
                if expected_domain is not None:
                    require(row["execution_domain"]==expected_domain,
                            f"A5 scope domain mismatch for {row['name']}")
        return {"status":"pass","untruncated_relative_l2":d["untruncated_receiver_relative_l2"],
                "nskip":d["nskip"],"counter_modes":2,"cut_index_rows":cut_rows}
    raise ValueError(f"unknown data gate {mode}")


def main():
    p=argparse.ArgumentParser()
    p.add_argument("mode",choices=["candidate","build","a1","a2","a3","a4","a5"])
    p.add_argument("--candidate",type=Path);p.add_argument("--data-dir",type=Path)
    p.add_argument("--build-result",type=Path);p.add_argument("--build-config",type=Path)
    p.add_argument("--build-root",type=Path);p.add_argument("--output",type=Path,required=True)
    p.add_argument("--input-sha256",action="append",default=[],metavar="PATH=SHA256")
    a=p.parse_args()
    bound=[]
    for item in a.input_sha256:
        path_text,digest=item.rsplit("=",1);path=Path(path_text)
        require(path.is_file() and sha(path)==digest,f"frozen input SHA256 mismatch: {path}")
        bound.append({"path":str(path),"sha256":digest})
    if a.mode=="candidate": result=candidate_gate(a.candidate)
    elif a.mode=="build": result=build_gate(a.build_result,a.build_config,a.build_root)
    else: result=data_gate(a.mode,a.data_dir)
    result["input_bindings"]=bound
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")

if __name__=="__main__":
    main()
