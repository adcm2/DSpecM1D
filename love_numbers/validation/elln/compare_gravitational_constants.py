#!/usr/bin/env python3

import argparse
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from compare_controlled_sphere import relative_difference, run_driver
from validate_elln import (
    audit_source,
    prepare_archive,
    read_numeric_rows,
)


DEFAULT_G = 6.67230e-11
ELLN_G = 6.67384e-11
DEGREE = 10
RADIAL_STEP = "0.00125"
OFFICIAL_H = -1.451586
OFFICIAL_K = -0.07066685
COEFFICIENTS = ("rho", "A", "C", "F", "L", "N")


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("default_executable")
    parser.add_argument("elln_g_executable")
    parser.add_argument("model", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("official_table", type=Path)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--archive", type=Path)
    source.add_argument("--source-dir", type=Path)
    parser.add_argument("--work-dir", type=Path, required=True)
    return parser.parse_args()


def read_deck(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    rows = tuple(tuple(map(float, line.split())) for line in lines[3:])
    if len(rows) != int(lines[2].split()[0]):
        raise RuntimeError("DSpecM1D deck row count is inconsistent")

    layers = []
    first = 0
    for index in range(1, len(rows)):
        if rows[index][0] == rows[index - 1][0]:
            layers.append(rows[first:index])
            first = index
    layers.append(rows[first:])
    if len(layers) != 83:
        raise RuntimeError(f"Expected 83 exported regions, got {len(layers)}")
    return tuple(layers)


def read_degree_ten_reference(path):
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if fields and not line.startswith("#") and int(fields[0]) == DEGREE:
            return -float(fields[1]), -float(fields[3])
    raise RuntimeError("Table S7 does not contain degree 10")


def exported_values(row):
    _, density, vpv, vsv, _, _, vph, vsh, eta = row
    a = density * vph * vph
    c = density * vpv * vpv
    l = density * vsv * vsv
    n = density * vsh * vsh
    f = eta * (a - 2.0 * l)
    return density, a, c, f, l, n


def core_values(row):
    _, _, density, bulk_modulus = row
    return density, bulk_modulus, bulk_modulus, bulk_modulus, 0.0, 0.0


def mantle_values(lower_row, upper_row, radius):
    lower = lower_row[1]
    upper = upper_row[1]
    equivalent_radius = (
        2.0 * (upper**3 - lower**3)
        / (3.0 * (upper**2 - lower**2))
    )
    density = upper_row[2] * equivalent_radius / radius
    a, c, l, n, f = upper_row[3:]
    return density, a, c, f, l, n


def report_audit_point(label, row, expected):
    actual = exported_values(row)
    print(f"audit_point label={label} radius={row[0]:.17e}")
    errors = {}
    for name, value, reference in zip(COEFFICIENTS, actual, expected):
        signed = value - reference
        relative = relative_difference(value, reference)
        errors[name] = abs(relative) if reference != 0.0 else abs(signed)
        print(
            f"audit_value label={label} name={name} "
            f"exported={value:.17e} direct={reference:.17e} "
            f"signed_difference={signed:.17e} "
            f"relative_difference={relative:.17e}"
        )
    return errors


def audit_model(model, source_directory):
    paths = audit_source(source_directory)
    core = read_numeric_rows(paths["EarthCore26.txt"], 4)
    mantle = read_numeric_rows(paths["EarthMantleTI56.txt"], 8)
    layers = read_deck(model)
    maximum_errors = {name: 0.0 for name in COEFFICIENTS}

    samples = (
        ("fluid_core", layers[13][len(layers[13]) // 2],
         core_values(core[13])),
        ("lower_mantle", layers[37][len(layers[37]) // 2],
         mantle_values(mantle[10], mantle[11],
                        layers[37][len(layers[37]) // 2][0])),
        ("upper_mantle", layers[72][len(layers[72]) // 2],
         mantle_values(mantle[45], mantle[46],
                        layers[72][len(layers[72]) // 2][0])),
    )
    for label, row, expected in samples:
        for name, error in report_audit_point(
                label, row, expected).items():
            maximum_errors[name] = max(maximum_errors[name], error)

    mantle_boundary = 40
    boundary_samples = (
        ("cmb_below", layers[26][-1], core_values(core[26])),
        ("cmb_above", layers[27][0],
         mantle_values(mantle[0], mantle[1], layers[27][0][0])),
        ("mantle_boundary_below", layers[26 + mantle_boundary][-1],
         mantle_values(
             mantle[mantle_boundary - 1], mantle[mantle_boundary],
             layers[26 + mantle_boundary][-1][0])),
        ("mantle_boundary_above", layers[27 + mantle_boundary][0],
         mantle_values(
             mantle[mantle_boundary], mantle[mantle_boundary + 1],
             layers[27 + mantle_boundary][0][0])),
    )
    for label, row, expected in boundary_samples:
        for name, error in report_audit_point(
                label, row, expected).items():
            maximum_errors[name] = max(maximum_errors[name], error)

    print(
        "audit_convention core=constant_properties_per_region "
        "mantle=upper_row_properties rho=rho_upper*r_bar/r "
        "boundaries=duplicated_with_each_side_owned_by_its_region "
        "units=r_m,rho_kg_m-3,elastic_coefficients_Pa"
    )
    print(
        "audit_maximum_errors "
        + " ".join(
            f"{name}={maximum_errors[name]:.17e}"
            for name in COEFFICIENTS
        )
    )
    return maximum_errors


def conventional(result, gravity, radius, gravitational_constant):
    factor = (
        (2 * DEGREE + 1)
        / (4.0 * math.pi * gravitational_constant * radius)
    )
    return (
        factor * gravity * result[4],
        -(1.0 + factor * result[5]),
    )


def run_case(label, executable, model, gravitational_constant):
    run = run_driver([
        executable, "--configured", str(model), "6", RADIAL_STEP,
        str(DEGREE),
    ])
    radius, gravity, elements, radial_nodes = run["metadata"]
    result = run["results"][DEGREE]
    diagnostic = run["diagnostics"][DEGREE]
    h, k = conventional(result, gravity, radius, gravitational_constant)
    residual = max(diagnostic[1:4])
    print(
        f"gravity_case label={label} G={gravitational_constant:.17e} "
        f"h_load={result[4]:.17e} k_load={result[5]:.17e} "
        f"surface_gravity={gravity:.17e} h_10={h:.17e} k_10={k:.17e} "
        f"elements={int(elements)} radial_nodes={int(radial_nodes)} "
        f"dofs={int(diagnostic[0])} solve_residual={residual:.17e} "
        f"reciprocity_error={diagnostic[4]:.17e}"
    )
    for name, value, reference in (
            ("h_10", h, OFFICIAL_H), ("k_10", k, OFFICIAL_K)):
        print(
            f"table_s7 label={label} name={name} "
            f"value={value:.17e} reference={reference:.17e} "
            f"signed_difference={value-reference:.17e} "
            f"relative_difference="
            f"{relative_difference(value, reference):.17e}"
        )
    return {"h": h, "k": k}


def main():
    arguments = parse_arguments()
    manifest = json.loads(arguments.manifest.read_text(encoding="utf-8"))
    if read_degree_ten_reference(arguments.official_table) != (
            OFFICIAL_H, OFFICIAL_K):
        raise RuntimeError("The audited Table S7 degree-10 row changed")
    if manifest["models"][-1]["file"] != arguments.model.name:
        raise RuntimeError("Stage 17 requires the finest exported model")

    default = run_case(
        "default_G", arguments.default_executable, arguments.model, DEFAULT_G
    )
    elln = run_case(
        "elln_G", arguments.elln_g_executable, arguments.model, ELLN_G
    )
    relative_g_change = (ELLN_G - DEFAULT_G) / DEFAULT_G
    for name, reference in (("h", OFFICIAL_H), ("k", OFFICIAL_K)):
        original = default[name] - reference
        residual = elln[name] - reference
        fraction = (original - residual) / original
        sensitivity = (elln[name] - default[name]) / relative_g_change
        print(
            f"sensitivity name={name} relative_G_change="
            f"{relative_g_change:.17e} response_change="
            f"{elln[name]-default[name]:.17e} "
            f"change_per_relative_G={sensitivity:.17e} "
            f"fraction_original_signed_discrepancy_explained="
            f"{fraction:.17e}"
        )

    arguments.work_dir.mkdir(parents=True, exist_ok=True)
    if arguments.archive:
        source_directory = prepare_archive(
            arguments.archive.resolve(), arguments.work_dir
        )
    else:
        source_directory = arguments.source_dir.resolve()
    audit_model(arguments.model, source_directory)
    print(
        "classification G_effect=majority_of_signed_discrepancy "
        "model_export=no_demonstrated_transcription_error_at_audited_points "
        "residual=ELLN_representation_tabulation_or_unresolved_inter_code"
    )


if __name__ == "__main__":
    main()
