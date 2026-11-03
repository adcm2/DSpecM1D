#!/usr/bin/env python3
"""Validate corrected Phase B three-bin mapper replay and qualified reuse."""
import argparse, csv, hashlib, json, math, pathlib, struct

EXPECTED_Q0 = {"nq5-nominal": 308, "nq5-finer": 2618, "nq6-nominal": 308, "nq6-finer": 0}
AFFECTED = ("nq5-nominal", "nq5-finer", "nq6-nominal")
BINS = (9, 25, 40)


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def read_json(path):
    return json.loads(pathlib.Path(path).read_text())


def read_csv(path):
    with pathlib.Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def sha(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def finite(value, where):
    x = float(value)
    need(math.isfinite(x), f"nonfinite {where}")
    return x


def rel_l2(err2, norm2):
    return math.sqrt(err2 / norm2) if norm2 else None


def key(row):
    return (int(row["bin"]), int(row["section"]), int(row["reference_element"]),
            int(row["reference_node"]), int(row["section_point"]))


def packet_key(row):
    return (int(row["bin"]), int(row["section"]), int(row["element"]),
            int(row["node"]), int(row["section_point"]))


def role_values(row, role, side):
    return [complex(finite(row[f"{role}_{side}{c}_re"], role),
                    finite(row[f"{role}_{side}{c}_im"], role)) for c in range(6)]


def grouped(items, compare_product=True):
    groups = {}
    for row, left, right, lp, rp, status in items:
        group = f"{row['section']}|{row['location_class']}|fluid={row['fluid']}|{status}"
        g = groups.setdefault(group, {"count": 0, "ferr": 0.0, "fnorm": 0.0,
             "aerr": 0.0, "anorm": 0.0, "kerr": 0.0, "knorm": 0.0})
        g["count"] += 1
        for x, y in zip(left[0], right[0]):
            g["ferr"] += abs(x-y)**2; g["fnorm"] += abs(y)**2
        for x, y in zip(left[1], right[1]):
            g["aerr"] += abs(x-y)**2; g["anorm"] += abs(y)**2
        if compare_product:
            g["kerr"] += (lp-rp)**2; g["knorm"] += rp*rp
    result = {}
    for group, g in groups.items():
        result[group] = {"count": g["count"],
            "forward_complex_relative_l2": rel_l2(g["ferr"], g["fnorm"]),
            "forward_complex_absolute_l2": math.sqrt(g["ferr"]), "forward_reference_norm2": g["fnorm"],
            "adjoint_complex_relative_l2": rel_l2(g["aerr"], g["anorm"]),
            "adjoint_complex_absolute_l2": math.sqrt(g["aerr"]), "adjoint_reference_norm2": g["anorm"],
            "kernel_product_relative_l2": rel_l2(g["kerr"], g["knorm"]) if compare_product else None,
            "kernel_product_absolute_error": math.sqrt(g["kerr"]), "kernel_product_reference_norm2": g["knorm"] if compare_product else None}
    return result


def actual_cuts(profile, grid, bin_id, mesh_elements, length_scale):
    gridrow = next(row for row in grid if int(row["bin"]) == bin_id)
    zr = struct.pack(">d", float(gridrow["z_real"])).hex()
    zi = struct.pack(">d", float(gridrow["z_imag"])).hex()
    found = {}
    for row in profile.get("counters", []):
        name = row["name"].removeprefix("frequency_worker.")
        tail = f"_zr{zr}_zi{zi}_actual_cut_element"
        calls_tail = f"_zr{zr}_zi{zi}_actual_cut_element_calls"
        radius_tail = f"_zr{zr}_zi{zi}_actual_cut_radius_nd_scaled_1e15"
        if not name.startswith("phase_b_solve_"):
            continue
        if name.endswith(calls_tail):
            fd = name[len("phase_b_solve_"):-len(calls_tail)]
            found.setdefault(fd, {})["calls"] = int(row["value"])
        elif name.endswith(radius_tail):
            fd = name[len("phase_b_solve_"):-len(radius_tail)]
            found.setdefault(fd, {})["radius_m"] = finite(row["value"], "actual cut radius") * length_scale / 1e15
        elif name.endswith(tail):
            fd = name[len("phase_b_solve_"):-len(tail)]
            found.setdefault(fd, {})["element"] = int(row["value"])
    need(found, f"missing actual per-family/degree cuts at bin {bin_id}")
    cuts = []
    for fd, value in sorted(found.items()):
        need(value.get("calls") == 1, f"cut {fd} bin {bin_id} does not have one solve")
        elem = value.get("element", 0)  # zero-valued counters may be omitted by the profile writer
        need(0 <= elem < len(mesh_elements), f"invalid cut element {elem} for {fd}")
        derived_radius = 0.0 if elem == 0 else finite(mesh_elements[elem]["lower_radius_m"], "mesh cut radius")
        counter_radius = value.get("radius_m")
        if counter_radius is not None:
            need(abs(counter_radius-derived_radius) <= 1e-3, f"quantized cut-radius counter disagrees with mesh lower edge for {fd}")
        if elem == 0:
            need(counter_radius in (None, 0.0), f"full-domain fallback cut radius should be zero for {fd}")
        cuts.append({"family_degree": fd, "actual_cut_element": elem,
                     "actual_cut_lower_radius_m": derived_radius,
                     "cut_radius_source": "mesh lower edge, quantized counter cross-check" if counter_radius is not None else "mesh lower edge (profile radius counter absent)"})
    return cuts


def mapper_q0_count(rows, mesh_rows):
    mesh = {int(x["element"]): x for x in mesh_rows}
    lower = []
    for row in mesh_rows:
        lower.append((int(row["layer"]), int(row["fluid"]), finite(row["lower_radius_m"], "mesh radius")))
    affected = set()
    for row in rows:
        if int(row["reference_node"]) != 0 or row["location_class"] != "bulk":
            continue
        radius = finite(row["radius_m"], "reference radius")
        elem = int(row["candidate_element"])
        if elem not in mesh:
            raise RuntimeError(f"candidate element {elem} absent from mesh inventory")
        upper = finite(mesh[elem]["upper_radius_m"], "mesh upper radius")
        same_material_lower = any(layer == int(row["layer"]) and fluid == int(row["fluid"])
                                  and abs(lo-radius) <= 1e-3 for layer, fluid, lo in lower)
        if abs(upper-radius) <= 1e-3 and same_material_lower:
            affected.add(key(row))
    return len(affected)


def state_job(state, job_id):
    jobs = {j["id"]: j for j in state["jobs"]}
    need(job_id in jobs, f"missing supervisor job {job_id}")
    return jobs[job_id]


def verify_state_outputs(state, root, job_id, names):
    job = state_job(state, job_id)
    need(job["status"] == "completed" and job["returncode"] == 0 and not job.get("failure_reasons"),
         f"job {job_id} did not complete successfully")
    need(job.get("inputs_unchanged") is True and job["cleanup"]["verified_gone"],
         f"job {job_id} input/cleanup receipt failed")
    for name in names:
        digest = job.get("outputs_sha256", {}).get(name)
        need(digest is not None, f"{job_id} missing output hash {name}")
        need(sha(root/job_id/name) == digest, f"output hash mismatch {job_id}/{name}")
    return job


def compare_file_rows(full, restricted):
    a = {key(x): x for x in full}; b = {key(x): x for x in restricted}
    need(a.keys() == b.keys(), "full/restricted common-side inventory differs")
    for k in a:
        for field in ("layer", "radius_m", "fluid", "location_class", "candidate_element"):
            need(a[k][field] == b[k][field], f"full/restricted metadata differs at {k}: {field}")
    return a, b


def build_metrics(full_rows, restricted_rows, cuts, products):
    a, b = compare_file_rows(full_rows, restricted_rows)
    full_ref = []; restricted_ref = []; delta = []; omission = {"affected": 0, "unaffected": 0}
    changed_mapping = 0; max_applicable_cut_radius = 0.0
    for k, x in a.items():
        y = b[k]
        applicable = [c["actual_cut_lower_radius_m"] for c in cuts
                      if x["fluid"] != "1" or "toroidal" not in c["family_degree"]]
        radius_cut = max(applicable, default=0.0)
        max_applicable_cut_radius = max(max_applicable_cut_radius, radius_cut)
        status = "affected_below_at_least_one_applicable_cut" if finite(x["radius_m"], "radius") < radius_cut else "outside_applicable_cut_radii"
        omission["affected" if status.startswith("affected") else "unaffected"] += 1
        fleft = (role_values(x, "forward", "candidate"), role_values(x, "adjoint", "candidate"))
        fright = (role_values(x, "forward", "ref"), role_values(x, "adjoint", "ref"))
        rleft = (role_values(y, "forward", "candidate"), role_values(y, "adjoint", "candidate"))
        rright = (role_values(y, "forward", "ref"), role_values(y, "adjoint", "ref"))
        fprod = finite(x["kernel_vs"], "full kernel product")
        rprod = finite(y["kernel_vs"], "restricted kernel product")
        packet = {"bin": x["bin"], "section": x["section"], "element": x["reference_element"],
                  "node": x["reference_node"], "section_point": x["section_point"]}
        pkey = packet_key(packet)
        need(pkey in products, f"missing frozen NQ10 product at {k}")
        refprod = products[pkey]
        full_ref.append((x, fleft, fright, fprod, refprod, status))
        restricted_ref.append((y, rleft, rright, rprod, refprod, status))
        delta.append((y, rleft, fleft, rprod, fprod, status))
        if int(x["candidate_element"]) != int(y["candidate_element"]):
            changed_mapping += 1
    return {"requested_samples_by_actual_cut_status": omission,
        "full_domain_vs_frozen_NQ10": grouped(full_ref),
        "restricted_domain_vs_frozen_NQ10": grouped(restricted_ref),
        "restricted_minus_full_same_frozen_forcing": grouped(delta),
        "full_vs_restricted_mapping_element_changes": changed_mapping,
        "maximum_applicable_cut_radius_across_samples_m": max_applicable_cut_radius}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task-root", type=pathlib.Path, required=True)
    ap.add_argument("--old-task-root", type=pathlib.Path, required=True)
    ap.add_argument("--old-state", type=pathlib.Path, required=True)
    ap.add_argument("--frozen-trace", type=pathlib.Path, required=True)
    ap.add_argument("--packet-dir", type=pathlib.Path, required=True)
    ap.add_argument("--nq10-products", type=pathlib.Path, required=True)
    ap.add_argument("--out", type=pathlib.Path, required=True)
    args = ap.parse_args()
    new_state = read_json(args.task_root/"state.json")
    old_state = read_json(args.old_state)
    need(new_state["task_id"] in str(args.task_root), "new task root does not identify its supervisor task")
    need(old_state["status"] == "failed", "expected preserved original Phase B task with checker-only terminal failure")
    products = {packet_key(r): finite(r["packed_reference_kernel_vs"], "frozen NQ10 kernel product")
                for r in read_csv(args.nq10_products)}
    trace_hash = sha(args.frozen_trace)
    packet_hashes = {f"reference-bin{b}.bin": sha(args.packet_dir/f"reference-bin{b}.bin") for b in BINS}
    expected_elements = {"nq5-nominal": 70, "nq5-finer": 132, "nq6-nominal": 53, "nq6-finer": 96}
    report = {"schema": 1, "status": "pass", "scope": "corrected three-bin common-side and truncation replay",
              "fresh_replays": {}, "qualified_saved_reuse": {}, "original_waveform_and_cost_preserved": {},
              "frozen_forcing": {"trace_sha256": trace_hash, "packet_sha256": packet_hashes},
              "accuracy_tolerances_applied": False,
              "restriction_interpretation": "principal H_R=P H P^T with eliminated cross-cut terms and zero-extension; not a physical-boundary claim",
              "omission_interpretation": "affected means radius below at least one actual applicable family/degree cut; it does not mean the summed field is zero"}
    common_replay_names = ["data/phase-b-diagnostic-replay.json", "data/profile.json", "data/frequency-grid.csv",
                           "data/mesh-elements.csv", "data/parameters.txt"]
    for mesh in AFFECTED:
        nq = int(mesh[2]); finer = "finer" in mesh
        mesh_result = {"mesh": mesh, "bins": {}, "replay_only": True, "waveform_solve_performed": False}
        fullroot = args.task_root/f"{mesh}-full-replay"/"data"
        # Exact job names are part of the controlled plan.
        full_id = f"{mesh}-full-replay"; restricted_id = f"{mesh}-restricted-replay"
        fullroot = args.task_root/full_id/"data"; restrroot = args.task_root/restricted_id/"data"
        common_names = common_replay_names + [f"data/common-fields-bin{b}.csv" for b in BINS] + [f"data/kernel-workspace-receiver-bin{b}.csv" for b in BINS]
        fulljob = verify_state_outputs(new_state, args.task_root, full_id, common_names)
        restrjob = verify_state_outputs(new_state, args.task_root, restricted_id, common_names)
        for item in (args.frozen_trace, *(args.packet_dir/f"reference-bin{b}.bin" for b in BINS)):
            path = str(item)
            expected_hash = sha(item)
            need(fulljob["inputs_sha256_before"].get(path) == restrjob["inputs_sha256_before"].get(path) == expected_hash, f"full/restricted frozen input hash mismatch: {path}")
        md_full = read_json(fullroot/"phase-b-diagnostic-replay.json")
        md_restr = read_json(restrroot/"phase-b-diagnostic-replay.json")
        for md, restricted in ((md_full, False), (md_restr, True)):
            need(md["status"] == "replay_only" and md["waveform_solve_performed"] is False, "replay metadata mislabels a waveform run")
            need(md["degree"] == 200 and md["nq"] == nq and md["element_count"] == expected_elements[mesh], f"replay mesh metadata mismatch: {mesh}")
            need(md["restricted"] is restricted and md["bins"] == list(BINS), "replay restriction/bin metadata mismatch")
            need(md["relative_error"] == 1e-5, "replay InputParameters relative_error changed")
            grid_rows = read_csv(fullroot/"frequency-grid.csv")
            need(md["mesh_selection_zmax_real"] == finite(grid_rows[-1]["z_real"], "mesh zmax") and md["mesh_selection_zmax_imag"] == finite(grid_rows[-1]["z_imag"], "mesh zmax"), "mesh selection did not use the retained full grid last complex frequency")
        need(sha(fullroot/"frequency-grid.csv") == sha(restrroot/"frequency-grid.csv"), "paired frequency grid differs")
        need(sha(fullroot/"mesh-elements.csv") == sha(restrroot/"mesh-elements.csv"), "paired mesh geometry differs")
        old_mesh_root = args.old_task_root/f"{mesh}-full"/"data"
        old_mesh_result = read_json(old_mesh_root/"result.json")
        need(old_mesh_result["element_count"] == expected_elements[mesh] and old_mesh_result["nq"] == nq and old_mesh_result["finer_mesh"] is finer, f"preserved mesh selection metadata mismatch: {mesh}")
        need(old_mesh_result["relative_error"] == md_full["relative_error"] and old_mesh_result["zmax_real"] == md_full["mesh_selection_zmax_real"] and old_mesh_result["zmax_imag"] == md_full["mesh_selection_zmax_imag"], f"replay mesh recipe differs from preserved Phase B run: {mesh}")
        need(sha(fullroot/"frequency-grid.csv") == sha(old_mesh_root/"frequency-grid.csv"), f"replay retained frequency grid differs from original {mesh}")
        need(sha(fullroot/"mesh-elements.csv") == sha(old_mesh_root/"mesh-elements.csv"), f"replay mesh geometry differs from original {mesh}")
        mesh_rows = read_csv(restrroot/"mesh-elements.csv")
        prof = read_json(restrroot/"profile.json")
        counters = {r["name"].removeprefix("frequency_worker."): int(r["value"]) for r in prof.get("counters", [])}
        need(counters.get("phase_b_restriction_principal_boundary_zero_extension", 0) > 0, "restricted replay did not execute principal zero-extension")
        need(counters.get("phase_b_restriction_omitted_dofs", 0) > 0, "restricted replay omitted no DOFs")
        hh_calls = sum(int(r.get("calls", 0)) for r in prof.get("rows", [])
                       if r.get("scope_name", "").removeprefix("frequency_worker.") in ("minimal_unit_basis_HH_solve", "full_weighted_HH_solve"))
        need(hh_calls > 0, "restricted replay has no profiled H^H solve calls")
        restriction_counters = {k: v for k, v in counters.items() if k.startswith("phase_b_restriction_")}
        length_scale_m = max(finite(row["upper_radius_m"], "surface radius") for row in mesh_rows)
        cut_rows = {str(b): actual_cuts(prof, read_csv(restrroot/"frequency-grid.csv"), b, mesh_rows, length_scale_m) for b in BINS}
        # Verify the one-sided mapper correction changed exactly the saved artificial-q0 population.
        map_changes = {}
        for b in BINS:
            old_full_rows = read_csv(args.old_task_root/f"{mesh}-full"/"data"/f"common-fields-bin{b}.csv")
            new_full_rows = read_csv(fullroot/f"common-fields-bin{b}.csv")
            old_by_key = {key(r): r for r in old_full_rows}; new_by_key = {key(r): r for r in new_full_rows}
            need(old_by_key.keys() == new_by_key.keys(), f"mapper replay location inventory changed: {mesh} bin{b}")
            changed = [k for k in old_by_key if old_by_key[k]["candidate_element"] != new_by_key[k]["candidate_element"]]
            expected = mapper_q0_count(old_full_rows, read_csv(args.old_task_root/f"{mesh}-full"/"data"/"mesh-elements.csv"))
            need(expected == EXPECTED_Q0[mesh], f"saved q0 mapper impact for {mesh} is {expected}, expected {EXPECTED_Q0[mesh]}")
            need(len(changed) == expected, f"corrected full element mapping changed {len(changed)} rows, expected {expected}, {mesh} bin{b}")
            need(all(old_by_key[k]["location_class"] == "bulk" and int(old_by_key[k]["reference_node"]) == 0 for k in changed), "mapper element change occurred outside affected artificial q0 bulk rows")
            old_trunc_rows = read_csv(args.old_task_root/f"{mesh}-restricted"/"data"/f"common-fields-bin{b}.csv")
            new_trunc_rows = read_csv(restrroot/f"common-fields-bin{b}.csv")
            old_trunc_by_key = {key(r): r for r in old_trunc_rows}; new_trunc_by_key = {key(r): r for r in new_trunc_rows}
            need(old_trunc_by_key.keys() == new_trunc_by_key.keys(), f"restricted mapper replay location inventory changed: {mesh} bin{b}")
            changed_trunc = [k for k in old_trunc_by_key if old_trunc_by_key[k]["candidate_element"] != new_trunc_by_key[k]["candidate_element"]]
            need(len(changed_trunc) == expected, f"corrected restricted element mapping changed {len(changed_trunc)} rows, expected {expected}, {mesh} bin{b}")
            map_changes[str(b)] = {"full_changed_location_count": len(changed), "restricted_changed_location_count": len(changed_trunc), "recomputed_saved_q0_impact_count": expected}
            full_rows = new_full_rows
            trunc_rows = read_csv(restrroot/f"common-fields-bin{b}.csv")
            metrics = build_metrics(full_rows, trunc_rows, cut_rows[str(b)], products)
            expected_rows = 293*(135+19)
            need(len(full_rows) == len(trunc_rows) == expected_rows, f"common side inventory mismatch {mesh} bin{b}")
            receiver_rows = read_csv(fullroot/f"kernel-workspace-receiver-bin{b}.csv")
            receiver_trunc = read_csv(restrroot/f"kernel-workspace-receiver-bin{b}.csv")
            need(len(receiver_rows) == len(receiver_trunc) == 1, f"receiver spectrum inventory mismatch {mesh} bin{b}")
            rsq = esq = 0.0
            for col in ("up_re", "up_im", "north_re", "north_im", "east_re", "east_im"):
                x = finite(receiver_trunc[0][col], "truncated receiver spectrum")
                y = finite(receiver_rows[0][col], "full receiver spectrum")
                esq += (x-y)**2; rsq += y*y
            metrics["restricted_kernel_workspace_receiver_vs_full_relative_l2"] = rel_l2(esq, rsq)
            metrics["actual_cuts_by_family_degree"] = cut_rows[str(b)]
            metrics["restricted_operator_diagnostic_counters"] = restriction_counters
            metrics["profiled_kernel_pass_HH_scope_calls"] = hh_calls
            metrics["mapper_correction_changed_rows"] = map_changes[str(b)]
            mesh_result["bins"][str(b)] = metrics
        old_state_mesh = []
        old_jobs_by_mode = {}
        for jid in (f"{mesh}-own", f"{mesh}-full", f"{mesh}-restricted"):
            old_outputs = ["data/result.json", "data/trace.csv", "data/components.csv", "data/frequency-grid.csv", "data/mesh-elements.csv"]
            if jid.endswith(("-full", "-restricted")):
                old_outputs += ["data/profile.json", "data/frequency-grid.csv", "data/mesh-elements.csv",
                    "data/common-fields-bin9.csv", "data/common-fields-bin25.csv", "data/common-fields-bin40.csv",
                    "data/kernel-workspace-receiver-bin9.csv", "data/kernel-workspace-receiver-bin25.csv", "data/kernel-workspace-receiver-bin40.csv"]
            old = verify_state_outputs(old_state, args.old_task_root, jid, old_outputs)
            old_jobs_by_mode[jid.rsplit("-",1)[-1]] = old
            old_state_mesh.append({"job_id": jid, "result_sha256": old["outputs_sha256"]["data/result.json"],
                "trace_sha256": old["outputs_sha256"]["data/trace.csv"], "components_sha256": old["outputs_sha256"]["data/components.csv"],
                "elapsed_s": read_json(args.old_task_root/jid/"data/result.json").get("elapsed_s")})
        for mode in ("full", "restricted"):
            oldjob = old_jobs_by_mode[mode]
            for path, expected_hash in ((args.frozen_trace, trace_hash), *[(args.packet_dir/f"reference-bin{b}.bin", packet_hashes[f"reference-bin{b}.bin"]) for b in BINS]):
                need(oldjob["inputs_sha256_before"].get(str(path)) == expected_hash, f"preserved original {mode} forcing differs for {mesh}: {path}")
        report["fresh_replays"][mesh] = mesh_result
        report["original_waveform_and_cost_preserved"][mesh] = old_state_mesh
    # The saved NQ6 finer control is a negative mapper-impact control and is reused only after full qualification.
    mesh = "nq6-finer"; old_full = args.old_task_root/f"{mesh}-full"/"data"; old_restricted = args.old_task_root/f"{mesh}-restricted"/"data"
    original = read_csv(old_full/"common-fields-bin9.csv"); q0 = mapper_q0_count(original, read_csv(old_full/"mesh-elements.csv"))
    need(q0 == EXPECTED_Q0[mesh] == 0, f"NQ6 finer no-impact qualification failed: {q0}")
    old_full_job = verify_state_outputs(old_state, args.old_task_root, f"{mesh}-full",
        ["data/result.json", "data/trace.csv", "data/components.csv", "data/profile.json", "data/frequency-grid.csv", "data/mesh-elements.csv",
         *(f"data/common-fields-bin{b}.csv" for b in BINS), *(f"data/kernel-workspace-receiver-bin{b}.csv" for b in BINS)])
    old_restr_job = verify_state_outputs(old_state, args.old_task_root, f"{mesh}-restricted",
        ["data/result.json", "data/trace.csv", "data/components.csv", "data/profile.json", "data/frequency-grid.csv", "data/mesh-elements.csv",
         *(f"data/common-fields-bin{b}.csv" for b in BINS), *(f"data/kernel-workspace-receiver-bin{b}.csv" for b in BINS)])
    for input_path, expected_hash in ((args.frozen_trace, trace_hash), *[(args.packet_dir/f"reference-bin{b}.bin", packet_hashes[f"reference-bin{b}.bin"]) for b in BINS]):
        need(old_full_job["inputs_sha256_before"].get(str(input_path)) == old_restr_job["inputs_sha256_before"].get(str(input_path)) == expected_hash, f"original NQ6 finer paired forcing differs: {input_path}")
    oldmeta = read_json(old_full/"result.json")
    need(oldmeta["element_count"] == expected_elements[mesh] and oldmeta["nq"] == 6 and oldmeta["finer_mesh"] is True, "NQ6 finer saved mesh metadata mismatch")
    need(sha(old_full/"frequency-grid.csv") == sha(old_restricted/"frequency-grid.csv") and sha(old_full/"mesh-elements.csv") == sha(old_restricted/"mesh-elements.csv"), "NQ6 finer full/restricted mesh mismatch")
    old_restr_profile = read_json(old_restricted/"profile.json")
    old_mesh_rows = read_csv(old_restricted/"mesh-elements.csv")
    old_grid_rows = read_csv(old_restricted/"frequency-grid.csv")
    old_length_scale = max(finite(row["upper_radius_m"], "surface radius") for row in old_mesh_rows)
    reused_bins = {}
    for b in BINS:
        cuts = actual_cuts(old_restr_profile, old_grid_rows, b, old_mesh_rows, old_length_scale)
        full_rows = read_csv(old_full/f"common-fields-bin{b}.csv")
        trunc_rows = read_csv(old_restricted/f"common-fields-bin{b}.csv")
        need(len(full_rows) == len(trunc_rows) == 293*(135+19), f"reused NQ6 finer location inventory mismatch bin{b}")
        metrics = build_metrics(full_rows, trunc_rows, cuts, products)
        rx = read_csv(old_full/f"kernel-workspace-receiver-bin{b}.csv")
        ry = read_csv(old_restricted/f"kernel-workspace-receiver-bin{b}.csv")
        esq = nsq = 0.0
        for col in ("up_re", "up_im", "north_re", "north_im", "east_re", "east_im"):
            x=finite(ry[0][col],"reused restricted receiver"); y=finite(rx[0][col],"reused full receiver")
            esq+=(x-y)**2;nsq+=y*y
        hh = sum(int(r.get("calls",0)) for r in old_restr_profile.get("rows",[])
                 if r.get("phase") == "kernel_frequency_parallel_pass_and_ordered_reduction_including_diagnostics"
                 and r.get("scope_name","").removeprefix("frequency_worker.") in ("minimal_unit_basis_HH_solve","full_weighted_HH_solve"))
        need(hh>0,"reused NQ6 finer restriction has no profiled H^H calls")
        metrics["restricted_kernel_workspace_receiver_vs_full_relative_l2"] = rel_l2(esq,nsq)
        metrics["actual_cuts_by_family_degree"] = cuts
        metrics["profiled_kernel_pass_HH_scope_calls"] = hh
        reused_bins[str(b)] = metrics
    report["qualified_saved_reuse"][mesh] = {"status": "reused_hash_verified_original_three_bin_rows_no_new_solve",
        "bins": reused_bins,
        "zero_q0_mapper_impact_count_recomputed": q0,
        "full_result_sha256": old_full_job["outputs_sha256"]["data/result.json"],
        "restricted_result_sha256": old_restr_job["outputs_sha256"]["data/result.json"],
        "full_grid_sha256": old_full_job["outputs_sha256"]["data/frequency-grid.csv"],
        "full_mesh_sha256": old_full_job["outputs_sha256"]["data/mesh-elements.csv"],
        "trace_input_sha256": old_full_job["inputs_sha256_before"][str(args.frozen_trace)],
        "packet_input_sha256": {f"bin{b}": old_full_job["inputs_sha256_before"][str(args.packet_dir/f"reference-bin{b}.bin")] for b in BINS},
        "full_restricted_common_receiver_profile_sha256": {side: {name: old_state_job["outputs_sha256"][name] for name in [*(f"data/common-fields-bin{b}.csv" for b in BINS), *(f"data/kernel-workspace-receiver-bin{b}.csv" for b in BINS), "data/profile.json"]}
            for side, old_state_job in (("full", old_full_job), ("restricted", old_restr_job))}}
    # NQ5 and NQ6 self-refinement compare corrected nominal rows against corrected/reused finer rows.
    for nq in (5, 6):
        nominal = f"nq{nq}-nominal"; finer = f"nq{nq}-finer"
        nfull = args.task_root/f"{nominal}-full-replay"/"data" if nominal in AFFECTED else None
        ffull = args.task_root/f"{finer}-full-replay"/"data" if finer in AFFECTED else args.old_task_root/f"{finer}-full"/"data"
        if nfull is None:
            nfull = args.task_root/f"{nominal}-full-replay"/"data"
        self_ref = {}
        for b in BINS:
            nr = {key(x): x for x in read_csv(nfull/f"common-fields-bin{b}.csv")}
            fr = {key(x): x for x in read_csv(ffull/f"common-fields-bin{b}.csv")}
            need(nr.keys() == fr.keys(), f"self-refinement side inventory mismatch NQ{nq} bin{b}")
            items=[]
            for k, x in nr.items():
                y=fr[k]
                for f in ("layer","radius_m","fluid","location_class"):
                    need(x[f]==y[f],f"self-refinement metadata changed NQ{nq} bin{b}: {f}")
                left=(role_values(x,"forward","candidate"),role_values(x,"adjoint","candidate")); right=(role_values(y,"forward","candidate"),role_values(y,"adjoint","candidate"))
                items.append((x,left,right,finite(x["kernel_vs"],"nominal product"),finite(y["kernel_vs"],"finer product"),"mesh_self_refinement"))
            self_ref[str(b)] = grouped(items)
        report.setdefault("nominal_vs_finer_common_field_and_kernel_self_refinement",{})[f"nq{nq}"] = {"nominal":"corrected replay","finer":"corrected replay" if finer in AFFECTED else "reused hash-qualified NQ6 finer original; q0 impact zero","bins":self_ref}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True)+"\n")
    print(json.dumps({"status": report["status"], "fresh_replay_meshes": list(report["fresh_replays"]), "reused_meshes": list(report["qualified_saved_reuse"]), "out": str(args.out)}))


if __name__ == "__main__":
    main()
