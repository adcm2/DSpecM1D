#!/usr/bin/env python3

import argparse
import csv
import math
from pathlib import Path
import subprocess
import sys


sys.dont_write_bytecode = True


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--elln-archive", type=Path)
    source.add_argument("--elln-source-dir", type=Path)
    return parser.parse_args()


def run(command, working_directory):
    completed = subprocess.run(
        [str(item) for item in command],
        cwd=working_directory,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        sys.stderr.write(completed.stdout)
        sys.stderr.write(completed.stderr)
        raise RuntimeError(
            f"Benchmark command failed: {' '.join(map(str, command))}"
        )
    return completed.stdout


def values(line):
    parsed = {}
    for field in line.split():
        if "=" in field:
            name, value = field.split("=", 1)
            parsed[name] = value
    return parsed


def number(text):
    value = float(text)
    if not math.isfinite(value):
        raise RuntimeError(f"Non-finite benchmark value: {text}")
    return value


def write_csv(path, comments, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        for comment in comments:
            stream.write(f"# {comment}\n")
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            formatted = {}
            for field in fields:
                value = row[field]
                formatted[field] = (
                    f"{value:.17e}" if isinstance(value, float) else value
                )
            writer.writerow(formatted)


def parse_controlled(text, model):
    gravity = {}
    components = {}
    residuals = {}
    reciprocity = {}
    for line in text.splitlines():
        fields = values(line)
        if line.startswith("gravity "):
            key = (int(fields["N"]), fields["implementation"])
            gravity[key] = abs(number(fields["signed_difference"]))
        elif line.startswith("component "):
            key = (int(fields["N"]), int(fields["l"]))
            components.setdefault(key, []).append(
                abs(number(fields["relative_difference"]))
            )
        elif line.startswith("residual "):
            key = (int(fields["N"]), int(fields["l"]))
            residuals.setdefault(key, {"dspec": [], "gia3d": []})
            residuals[key]["dspec"].append(number(fields["dspec"]))
            residuals[key]["gia3d"].append(number(fields["gia3d"]))
        elif line.startswith("reciprocity "):
            key = (int(fields["N"]), int(fields["l"]))
            reciprocity[key] = {
                "dspec": abs(number(fields["dspec"])),
                "gia3d": abs(number(fields["gia3d"])),
            }

    element_counts = sorted({key[0] for key in components})
    convergence = []
    for element_count in element_counts:
        for degree in (1, 2, 10):
            key = (element_count, degree)
            convergence.append({
                "model": model,
                "elements": element_count,
                "degree": degree,
                "maximum_relative_difference": max(components[key]),
                "maximum_dspec_residual": max(residuals[key]["dspec"]),
                "maximum_gia3d_residual": max(residuals[key]["gia3d"]),
                "absolute_dspec_reciprocity_error":
                    reciprocity[key]["dspec"],
                "absolute_gia3d_reciprocity_error":
                    reciprocity[key]["gia3d"],
            })
    return gravity, convergence


def parse_isotropic(text):
    comparisons = []
    diagnostics = []
    for line in text.splitlines():
        fields = values(line)
        if line.startswith("comparison ") and "l" in fields:
            comparisons.append({
                "direction": fields["direction"],
                "knot_spacing": fields["spacing"],
                "radial_step": fields["step"],
                "degree": int(fields["l"]),
                "component": fields["name"],
                "signed_difference": number(fields["signed_difference"]),
                "relative_difference": number(fields["relative_difference"]),
            })
        elif line.startswith("solve "):
            diagnostics.append({
                "knot_spacing": fields["spacing"],
                "radial_step": fields["step"],
                "implementation": fields["implementation"],
                "degree": int(fields["l"]),
                "dofs": int(fields["dofs"]),
                "maximum_residual": max(
                    number(fields["displacement_residual"]),
                    number(fields["gravitational_residual"]),
                    number(fields["tide_residual"]),
                ),
                "absolute_reciprocity_error":
                    abs(number(fields["reciprocity_error"])),
            })
    return comparisons, diagnostics


def parse_elln(text):
    cases = {}
    table = {}
    sensitivity = []
    for line in text.splitlines():
        fields = values(line)
        if line.startswith("gravity_case "):
            label = fields["label"]
            cases[label] = {
                "case": label,
                "gravitational_constant": number(fields["G"]),
                "h_load": number(fields["h_load"]),
                "k_load": number(fields["k_load"]),
                "surface_gravity": number(fields["surface_gravity"]),
                "h_10": number(fields["h_10"]),
                "k_10": number(fields["k_10"]),
                "elements": int(fields["elements"]),
                "radial_nodes": int(fields["radial_nodes"]),
                "dofs": int(fields["dofs"]),
                "maximum_solve_residual": number(fields["solve_residual"]),
                "absolute_reciprocity_error":
                    abs(number(fields["reciprocity_error"])),
            }
        elif line.startswith("table_s7 "):
            table[(fields["label"], fields["name"])] = {
                "signed": number(fields["signed_difference"]),
                "relative": number(fields["relative_difference"]),
            }
        elif line.startswith("sensitivity "):
            sensitivity.append({
                "component": fields["name"],
                "relative_gravitational_constant_change":
                    number(fields["relative_G_change"]),
                "response_change": number(fields["response_change"]),
                "change_per_relative_gravitational_constant":
                    number(fields["change_per_relative_G"]),
                "fraction_signed_discrepancy_explained": number(
                    fields[
                        "fraction_original_signed_discrepancy_explained"
                    ]
                ),
            })

    rows = []
    for label in ("default_G", "elln_G"):
        row = cases[label]
        row.update({
            "h_signed_difference": table[(label, "h_10")]["signed"],
            "h_relative_difference": table[(label, "h_10")]["relative"],
            "k_signed_difference": table[(label, "k_10")]["signed"],
            "k_relative_difference": table[(label, "k_10")]["relative"],
        })
        rows.append(row)
    return rows, sensitivity


def main():
    options = arguments()
    root = options.source_root.resolve()
    build = options.build_dir.resolve()
    output = options.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    validation = root / "love_numbers" / "validation"
    binary = build / "bin"

    controlled = {}
    for model, dspec, gia3d in (
        ("homogeneous", "homogeneous_sphere.dspec",
         "homogeneous_sphere.gia3d"),
        ("solid-fluid-solid", "solid_fluid_solid.dspec",
         "solid_fluid_solid.gia3d"),
        ("central-fluid", "central_fluid.dspec",
         "central_fluid.gia3d"),
    ):
        command = [
            sys.executable,
            validation / "compare_controlled_sphere.py",
            binary / "dspecm1d_controlled_love_numbers",
            binary / "dspecm1d_gia3d_controlled_reference",
            validation / "data" / dspec,
            validation / "data" / gia3d,
        ]
        if model != "homogeneous":
            command.append(model)
        controlled[model] = parse_controlled(run(command, root), model)

    gravity_rows = []
    for elements in (8, 16, 32, 64, 128, 256):
        gravity = controlled["solid-fluid-solid"][0]
        gravity_rows.append({
            "elements": elements,
            "legacy_mesh_model_absolute_error":
                gravity[(elements, "dspec_legacy")],
            "corrected_love_number_absolute_error":
                gravity[(elements, "dspec")],
            "gia3d_absolute_error": gravity[(elements, "gia3d")],
        })
    write_csv(
        output / "background_gravity_sfs.csv",
        (
            "model=controlled solid-fluid-solid sphere",
            "quantity=absolute surface-gravity error relative to the "
            "piecewise analytic density integral; units=m s^-2",
            "resolution=five GLL nodes per element; total element count N",
            "source=compare_controlled_sphere.py and its two controlled "
            "validation executables",
        ),
        (
            "elements",
            "legacy_mesh_model_absolute_error",
            "corrected_love_number_absolute_error",
            "gia3d_absolute_error",
        ),
        gravity_rows,
    )

    sfs_rows = controlled["solid-fluid-solid"][1]
    write_csv(
        output / "controlled_sfs_convergence.csv",
        (
            "model=controlled solid-fluid-solid sphere",
            "quantity=maximum absolute relative DSpecM1D-gia3D difference "
            "over h_u,k_u,h_phi,k_phi,h_load,k_load,h_t,k_t",
            "resolution=five GLL nodes per element; total element count N",
            "source=compare_controlled_sphere.py",
        ),
        tuple(sfs_rows[0]),
        sfs_rows,
    )

    summary_rows = []
    for model in ("homogeneous", "central-fluid", "solid-fluid-solid"):
        rows = controlled[model][1]
        finest = max(row["elements"] for row in rows)
        summary_rows.extend(
            row for row in rows if row["elements"] == finest
        )
    write_csv(
        output / "controlled_finest_summary.csv",
        (
            "models=homogeneous sphere, central fluid, and "
            "solid-fluid-solid controlled comparisons",
            "quantity=finest-mesh component difference, three-column "
            "residuals, and reciprocity diagnostics",
            "source=compare_controlled_sphere.py",
        ),
        tuple(summary_rows[0]),
        summary_rows,
    )

    isotropic_command = [
        sys.executable,
        validation / "compare_isotropic_prem.py",
        binary / "dspecm1d_controlled_love_numbers",
        binary / "dspecm1d_gia3d_controlled_reference",
        binary / "dspecm1d_gia3d_export_isotropic_prem",
        build / "love_numbers" / "validation" / "isotropic-prem-models",
    ]
    isotropic, isotropic_diagnostics = parse_isotropic(
        run(isotropic_command, validation)
    )
    directions = (
        "dspec_sampled_minus_gia3d_sampled",
        "gia3d_sampled_minus_gia3d_analytic",
        "dspec_sampled_minus_gia3d_analytic",
    )
    spacing_values = {"100km": 100.0, "50km": 50.0,
                      "25km": 25.0, "12.5km": 12.5}
    isotropic_convergence = []
    for spacing, spacing_km in spacing_values.items():
        for direction in directions:
            selected = [
                abs(row["relative_difference"])
                for row in isotropic
                if row["knot_spacing"] == spacing
                and row["radial_step"] == "0.0025"
                and row["direction"] == direction
            ]
            isotropic_convergence.append({
                "knot_spacing_km": spacing_km,
                "direction": direction,
                "maximum_relative_difference": max(selected),
            })
    write_csv(
        output / "isotropic_prem_knot_convergence.csv",
        (
            "model=pinned gia3D analytic no-ocean PREM and paired "
            "isotropized sampled decks",
            "quantity=maximum absolute relative difference over degrees "
            "1,2,10,50 and eight response components",
            "resolution=radial step 0.0025; five GLL nodes",
            "source=compare_isotropic_prem.py",
        ),
        tuple(isotropic_convergence[0]),
        isotropic_convergence,
    )
    isotropic_finest = [
        row for row in isotropic
        if row["knot_spacing"] == "12.5km"
        and row["radial_step"] == "0.0025"
    ]
    write_csv(
        output / "isotropic_prem_finest.csv",
        (
            "model=12.5 km paired isotropized PREM decks and pinned "
            "analytic gia3D PREM",
            "quantity=signed and relative response differences at radial "
            "step 0.0025",
            "source=compare_isotropic_prem.py",
        ),
        tuple(isotropic_finest[0]),
        isotropic_finest,
    )
    finest_diagnostics = [
        row for row in isotropic_diagnostics
        if row["knot_spacing"] == "12.5km"
        and row["radial_step"] == "0.0025"
    ]
    write_csv(
        output / "isotropic_prem_diagnostics.csv",
        (
            "model=12.5 km paired isotropized PREM decks and analytic PREM",
            "quantity=DOF count, maximum three-column solve residual, and "
            "absolute reciprocity error",
            "source=compare_isotropic_prem.py",
        ),
        tuple(finest_diagnostics[0]),
        finest_diagnostics,
    )

    elln_model_dir = (
        build / "love_numbers" / "validation" / "elln-ti-prem-models"
    )
    elln_command = [
        sys.executable,
        validation / "elln" / "compare_gravitational_constants.py",
        binary / "dspecm1d_controlled_love_numbers",
        binary / "dspecm1d_controlled_love_numbers_elln_g",
        elln_model_dir / "elln_ti_prem_6.25km.dspec",
        elln_model_dir / "elln_ti_prem_manifest.json",
        validation / "elln" / "official_ti_prem_table_s7.txt",
        "--work-dir",
        build / "love_numbers" / "validation" / "report-elln-audit",
    ]
    if options.elln_archive:
        elln_command.extend(["--archive", options.elln_archive.resolve()])
    else:
        elln_command.extend(
            ["--source-dir", options.elln_source_dir.resolve()]
        )
    elln_rows, sensitivity = parse_elln(run(elln_command, root))
    write_csv(
        output / "elln_gravitational_constant.csv",
        (
            "model=official ELLN TI PREM sampled at 6.25 km",
            "reference=ELLN Table S7 h_10=-1.451586, k_10=-0.07066685",
            "resolution=polynomial degree 6; maximum radial step 0.00125",
            "quantity=dimensional load responses, conventional responses, "
            "and diagnostics for production and ELLN G",
            "source=compare_gravitational_constants.py",
        ),
        tuple(elln_rows[0]),
        elln_rows,
    )
    write_csv(
        output / "elln_sensitivity_summary.csv",
        (
            "model=official ELLN TI PREM sampled at 6.25 km",
            "quantity=response sensitivity and fraction of signed Table S7 "
            "discrepancy explained by changing G",
            "source=compare_gravitational_constants.py",
        ),
        tuple(sensitivity[0]),
        sensitivity,
    )
    print(f"Wrote benchmark datasets to {output}")


if __name__ == "__main__":
    main()
