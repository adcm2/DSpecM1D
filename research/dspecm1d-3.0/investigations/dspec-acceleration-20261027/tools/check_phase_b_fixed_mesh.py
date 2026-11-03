#!/usr/bin/env python3
"""Validate fixed-mesh Phase B controls and summarize saved complex outputs."""
import csv, json, math, pathlib, struct, sys

if len(sys.argv) != 5:
    raise SystemExit("usage: check_phase_b_fixed_mesh.py TASK_ROOT FROZEN_NQ10_COMPONENTS FROZEN_NQ10_PRODUCTS OUTPUT_JSON")
root, frozen_components, frozen_products, output = map(pathlib.Path, sys.argv[1:])
meshes = [(5, "nominal"), (5, "finer"), (6, "nominal"), (6, "finer")]


def data_dir(mesh, mode):
    return root / f"nq{mesh[0]}-{mesh[1]}-{mode}" / "data"


def read_json(path):
    return json.loads(path.read_text())


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def rel_l2(err2, ref2):
    return math.sqrt(err2 / max(ref2, 1e-300))


def finite(value, where):
    x = float(value)
    if not math.isfinite(x):
        raise RuntimeError(f"nonfinite value at {where}")
    return x


def waveform_stats(path, reference):
    a, b = read_csv(path), read_csv(reference)
    if len(a) != len(b):
        raise RuntimeError(f"waveform component row count mismatch: {path}")
    keys = ("up_re", "up_im", "north_re", "north_im", "east_re", "east_im")
    err = norm = 0.0
    for x, y in zip(a, b):
        if x["bin"] != y["bin"]:
            raise RuntimeError("waveform bin ordering mismatch")
        for key in keys:
            d, v = finite(x[key], path)-finite(y[key], reference), finite(y[key], reference)
            err += d*d; norm += v*v
    return {"bins": len(a), "relative_l2_to_reference": rel_l2(err, norm)}


def location_key(row):
    return (int(row["bin"]), int(row["section"]), int(row["reference_element"]),
            int(row["reference_node"]), int(row["section_point"]))


def packet_key(row):
    return (int(row["bin"]), int(row["section"]), int(row["element"]),
            int(row["node"]), int(row["section_point"]))


def field_values(row, role, part):
    return [complex(finite(row[f"{role}_{part}{c}_re"], "common field"),
                    finite(row[f"{role}_{part}{c}_im"], "common field")) for c in range(6)]


def grouped_error(items):
    groups = {}
    for row, candidate, reference, product_candidate, product_reference, state in items:
        group_key = f"{row['location_class']}|fluid={row['fluid']}|{state}"
        g = groups.setdefault(group_key, {"count": 0, "ferr": 0.0, "fnorm": 0.0,
            "aerr": 0.0, "anorm": 0.0, "kerr": 0.0, "knorm": 0.0,
            "fluid_kernel_nonzero_count": 0})
        g["count"] += 1
        for c, r in zip(candidate[0], reference[0]): g["ferr"] += abs(c-r)**2; g["fnorm"] += abs(r)**2
        for c, r in zip(candidate[1], reference[1]): g["aerr"] += abs(c-r)**2; g["anorm"] += abs(r)**2
        if product_reference is not None:
            g["kerr"] += (product_candidate-product_reference)**2
            g["knorm"] += product_reference**2
        if row["fluid"] == "1" and product_candidate != 0.0:
            g["fluid_kernel_nonzero_count"] += 1
    return {k: {"count": g["count"],
        "forward_complex_relative_l2": rel_l2(g["ferr"], g["fnorm"]),
        "adjoint_complex_relative_l2": rel_l2(g["aerr"], g["anorm"]),
        "kernel_product_relative_l2": rel_l2(g["kerr"], g["knorm"]) if g["knorm"] else None,
        "fluid_kernel_nonzero_count": g["fluid_kernel_nonzero_count"]} for k, g in groups.items()}


def compare_candidate_files(path_a, path_b):
    a = {location_key(x): x for x in read_csv(path_a)}
    b = {location_key(x): x for x in read_csv(path_b)}
    if a.keys() != b.keys(): raise RuntimeError(f"common-location inventory mismatch: {path_a} vs {path_b}")
    groups = {}
    for key, x in a.items():
        y = b[key]
        if any(x[k] != y[k] for k in ("layer", "radius_m", "fluid", "location_class")):
            raise RuntimeError("common physical side metadata mismatch")
        group_key = f"{x['location_class']}|fluid={x['fluid']}"
        g = groups.setdefault(group_key, {"count": 0, "ferr": 0.0, "fnorm": 0.0,
            "aerr": 0.0, "anorm": 0.0, "kerr": 0.0, "knorm": 0.0})
        g["count"] += 1
        for role, ek, nk in (("forward", "ferr", "fnorm"), ("adjoint", "aerr", "anorm")):
            va, vb = field_values(x, role, "candidate"), field_values(y, role, "candidate")
            for u, v in zip(va, vb): g[ek] += abs(u-v)**2; g[nk] += abs(u)**2
        ka, kb = finite(x["kernel_vs"], path_a), finite(y["kernel_vs"], path_b)
        g["kerr"] += (ka-kb)**2; g["knorm"] += ka*ka
    return {k: {"count": g["count"], "forward_complex_relative_l2": rel_l2(g["ferr"], g["fnorm"]),
        "adjoint_complex_relative_l2": rel_l2(g["aerr"], g["anorm"]),
        "kernel_product_relative_l2": rel_l2(g["kerr"], g["knorm"])} for k, g in groups.items()}


def cuts_for_bin(profile, grid_rows, mesh_elements, bin_id):
    grid = next(x for x in grid_rows if int(x["bin"]) == bin_id)
    zre, zim = float(grid["z_real"]), float(grid["z_imag"])
    zr = struct.pack(">d", zre).hex(); zi = struct.pack(">d", zim).hex()
    counters = {}
    for row in profile.get("counters", []):
        if row.get("phase") != "kernel_frequency_parallel_pass_and_ordered_reduction_including_diagnostics": continue
        name = row["name"].removeprefix("frequency_worker.")
        prefix = "phase_b_solve_"
        suffix = f"_zr{zr}_zi{zi}_actual_cut_element"
        calls_suffix = f"_zr{zr}_zi{zi}_actual_cut_element_calls"
        if name.startswith(prefix) and name.endswith(suffix):
            family_degree = name[len(prefix):-len(suffix)]
            counters.setdefault(family_degree, {})["cut_sum"] = int(row["value"])
        elif name.startswith(prefix) and name.endswith(calls_suffix):
            family_degree = name[len(prefix):-len(calls_suffix)]
            counters.setdefault(family_degree, {})["calls"] = int(row["value"])
    result = []
    for family_degree, values in sorted(counters.items()):
        calls, cut_sum = values.get("calls", 0), values.get("cut_sum", 0)
        if calls != 1:
            raise RuntimeError(f"cut key {family_degree} bin{bin_id} has {calls} solves; cannot map exact cut")
        cut = cut_sum
        if cut < 0 or cut > max(int(x["element"]) for x in mesh_elements):
            raise RuntimeError(f"invalid cut element {cut} for {family_degree}")
        radius = 0.0 if cut == 0 else finite(mesh_elements[cut]["lower_radius_m"], "mesh cut radius")
        result.append({"family_degree": family_degree, "cut_element": cut, "cut_lower_radius_m": radius})
    if not result: raise RuntimeError(f"no exact restricted cut inventory for bin{bin_id}")
    return result


summary = {"schema": 1, "phase": "B", "mesh_rows": [], "accuracy_tolerances_applied": False,
    "interpretation": "Measured errors are descriptive; no continuous-accuracy cutoff is imposed.",
    "truncation_boundary": "Principal H_R=P H P^T with eliminated cross-cut entries and zero-extension output; not a physical-boundary claim.",
    "omission_label": "A requested common sample is marked affected when its radius is below at least one actual per-family/degree cut radius; this does not imply its summed field is zero."}
ref_products = {packet_key(x): finite(x["packed_reference_kernel_vs"], "NQ10 kernel product") for x in read_csv(frozen_products)}
base_components = data_dir((5, "nominal"), "own") / "components.csv"
for mesh in meshes:
    label = f"nq{mesh[0]}-{mesh[1]}"; own=data_dir(mesh,"own"); full=data_dir(mesh,"full"); trunc=data_dir(mesh,"restricted")
    own_result, full_result, trunc_result = (read_json(x/"result.json") for x in (own,full,trunc))
    own_trace, full_trace = read_csv(own/"trace.csv"), read_csv(full/"trace.csv")
    if len(own_trace)!=len(full_trace): raise RuntimeError(f"own/full trace sample mismatch {label}")
    initial_pair = waveform_stats(trunc/"components.csv", full/"components.csv")
    own_full_pair = waveform_stats(full/"components.csv", own/"components.csv")
    if "own waveform-derived" not in own_result["kernel_adjoint_weight_source"]: raise RuntimeError(f"own weight label missing {label}")
    if "frozen accepted" not in full_result["kernel_adjoint_weight_source"]: raise RuntimeError(f"frozen weight label missing {label}")
    row={"mesh":label,"own_weights_source":own_result["kernel_adjoint_weight_source"],
         "frozen_control_weights_source":full_result["kernel_adjoint_weight_source"],
         "truncated_weights_source":trunc_result["kernel_adjoint_weight_source"],
         "own_trace_samples":len(own_trace),"n_elements":full_result["element_count"],
         "native_sides":full_result["native_side_count"],"nominal_h":full_result["nominal_h"],
         "finer_mesh":full_result["finer_mesh"],"frequency_count":full_result["frequency_count"],
         "first_bin":full_result["first_bin"],"last_bin":full_result["last_bin"],
         "time_shift_adjoint_identity":full_result["time_shift_adjoint_identity"],
         "elapsed_s_full":full_result["elapsed_s"],"elapsed_s_restricted":trunc_result["elapsed_s"],
         "full_vs_restricted_initial_waveform_relative_l2":initial_pair["relative_l2_to_reference"],
         "frozen_control_vs_own_run_initial_waveform_relative_l2":own_full_pair["relative_l2_to_reference"],"common_bins":{}}
    profile=read_json(trunc/"profile.json");mesh_elements=read_csv(trunc/"mesh-elements.csv")
    if len(mesh_elements)!=int(trunc_result["element_count"]): raise RuntimeError(f"mesh element inventory mismatch {label}")
    element_by_id={int(x["element"]):x for x in mesh_elements}
    if sorted(element_by_id)!=list(range(len(mesh_elements))): raise RuntimeError("mesh element IDs are not contiguous")
    grid=read_csv(trunc/"frequency-grid.csv")
    per_solve={x["name"].removeprefix("frequency_worker."):int(x["value"]) for x in profile.get("counters",[])
               if x.get("phase")=="kernel_frequency_parallel_pass_and_ordered_reduction_including_diagnostics"
               and x["name"].removeprefix("frequency_worker.").startswith("phase_b_solve_")}
    aggregate={x["name"].removeprefix("frequency_worker."):int(x["value"]) for x in profile.get("counters",[])
               if x["name"].removeprefix("frequency_worker.").startswith("phase_b_restriction_")}
    hh_calls=sum(x.get("calls",0) for x in profile.get("rows",[])
        if x.get("phase")=="kernel_frequency_parallel_pass_and_ordered_reduction_including_diagnostics"
        and x.get("scope_name","").removeprefix("frequency_worker.") in ("minimal_unit_basis_HH_solve","full_weighted_HH_solve"))
    row["restricted_solver_counters_quantized"]=aggregate
    row["per_solve_cut_counters_quantized"]=per_solve
    row["actual_restricted_HH_solve_calls"]=hh_calls
    if aggregate.get("phase_b_restriction_omitted_dofs",0)<=0 or hh_calls<=0:
        raise RuntimeError(f"no demonstrated restricted H/H^H work for {label}")
    for bin_id in (9,25,40):
        fpath=full/f"common-fields-bin{bin_id}.csv"; tpath=trunc/f"common-fields-bin{bin_id}.csv"
        frows,trows=read_csv(fpath),read_csv(tpath)
        expected=293*(135+19)
        if len(frows)!=expected or len(trows)!=expected: raise RuntimeError(f"common-side count mismatch {label} bin{bin_id}")
        cuts=cuts_for_bin(profile,grid,mesh_elements,bin_id)
        radii=[x["cut_lower_radius_m"] for x in cuts]
        low,high=min(radii),max(radii)
        full_ref_items=[];restricted_ref_items=[];truncation_delta_items=[]
        omission_counts={"affected_by_at_least_one_relevant_cut":0,"unaffected_by_relevant_cuts":0}
        for x,y in zip(frows,trows):
            kx,ky=location_key(x),location_key(y)
            if kx!=ky: raise RuntimeError(f"full/restricted side identity mismatch {label} bin{bin_id}")
            if any(x[m]!=y[m] for m in ("layer","radius_m","fluid","location_class")): raise RuntimeError("paired common metadata changed")
            relevant_cuts=[c for c in cuts if x["fluid"]!="1" or "toroidal" not in c["family_degree"]]
            relevant_high=max(c["cut_lower_radius_m"] for c in relevant_cuts)
            affected=float(x["radius_m"])<relevant_high
            state="affected_by_at_least_one_relevant_family_degree_cut" if affected else "outside_all_relevant_cut_radii"
            omission_counts["affected_by_at_least_one_relevant_cut" if affected else "unaffected_by_relevant_cuts"]+=1
            reference_product=ref_products.get(packet_key({"bin":x["bin"],"section":x["section"],"element":x["reference_element"],"node":x["reference_node"],"section_point":x["section_point"]}))
            if reference_product is None: raise RuntimeError(f"missing NQ10 reference kernel product at {label} bin{bin_id} location {kx}")
            trunc_candidate=(field_values(y,"forward","candidate"),field_values(y,"adjoint","candidate"))
            trunc_reference=(field_values(y,"forward","ref"),field_values(y,"adjoint","ref"))
            full_candidate=(field_values(x,"forward","candidate"),field_values(x,"adjoint","candidate"))
            full_reference=(field_values(x,"forward","ref"),field_values(x,"adjoint","ref"))
            trunc_product=finite(y["kernel_vs"],tpath);full_product=finite(x["kernel_vs"],fpath)
            full_ref_items.append((x,full_candidate,full_reference,full_product,reference_product,state))
            restricted_ref_items.append((y,trunc_candidate,trunc_reference,trunc_product,reference_product,state))
            truncation_delta_items.append((y,trunc_candidate,full_candidate,trunc_product,full_product,state))
        row["common_bins"][str(bin_id)]={"one_sided_location_rows":len(frows),
            "exact_cut_radii_by_family_degree":cuts,"cut_radius_range_m":{"minimum":low,"maximum":high},
            "requested_samples_by_cut_status":omission_counts,
            "full_domain_vs_frozen_NQ10":grouped_error(full_ref_items),
            "restricted_domain_vs_frozen_NQ10":grouped_error(restricted_ref_items),
            "restricted_minus_full_same_forcing":grouped_error(truncation_delta_items),
            "restricted_kernel_workspace_receiver_vs_full_relative_l2":waveform_stats(trunc/f"kernel-workspace-receiver-bin{bin_id}.csv",full/f"kernel-workspace-receiver-bin{bin_id}.csv")["relative_l2_to_reference"]}
    summary["mesh_rows"].append(row)
for nq in (5,6):
    summary.setdefault("own_waveform_nominal_vs_finer",{})[f"nq{nq}"]=waveform_stats(data_dir((nq,"finer"),"own")/"components.csv",data_dir((nq,"nominal"),"own")/"components.csv")
for mesh in meshes:
    summary.setdefault("own_waveform_vs_frozen_nq10",{})[f"nq{mesh[0]}-{mesh[1]}"]=waveform_stats(data_dir(mesh,"own")/"components.csv",frozen_components)
for nq in (5,6):
    for bin_id in (9,25,40):
        summary.setdefault("common_nominal_vs_finer",{}).setdefault(f"nq{nq}",{})[str(bin_id)]=compare_candidate_files(
            data_dir((nq,"nominal"),"full")/f"common-fields-bin{bin_id}.csv",
            data_dir((nq,"finer"),"full")/f"common-fields-bin{bin_id}.csv")
summary["no_accuracy_threshold_applied"]=True
output.parent.mkdir(parents=True,exist_ok=True)
output.write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
