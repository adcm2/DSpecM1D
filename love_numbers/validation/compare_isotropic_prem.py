#!/usr/bin/env python3

import hashlib
import math
import pathlib
import sys

from compare_controlled_sphere import (
    COMPONENTS,
    relative_difference,
    run_driver,
)


DEGREES = ("1", "2", "10", "50")
KNOT_SPACINGS = (
    ("100km", 100_000.0),
    ("50km", 50_000.0),
    ("25km", 25_000.0),
    ("12.5km", 12_500.0),
)
MESH_STEPS = (
    ("0.01", 100),
    ("0.005", 200),
    ("0.0025", 400),
)
LAYER_RADII = (
    0.0,
    1_221_500.0,
    3_480_000.0,
    3_630_000.0,
    5_600_000.0,
    5_701_000.0,
    5_771_000.0,
    5_971_000.0,
    6_151_000.0,
    6_291_000.0,
    6_346_600.0,
    6_356_000.0,
    6_368_000.0,
)


def checksum(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_exported_decks(dspec_path, gia3d_path, maximum_spacing):
    with dspec_path.open(encoding="utf-8") as stream:
        dspec_lines = stream.readlines()
    with gia3d_path.open(encoding="utf-8") as stream:
        gia3d_lines = stream.readlines()

    if len(dspec_lines) < 4:
        raise RuntimeError(f"incomplete DSpecM1D deck: {dspec_path}")
    declared_knots = int(dspec_lines[2].split()[0])
    dspec_rows = tuple(
        tuple(map(float, line.split())) for line in dspec_lines[3:]
    )
    gia3d_rows = tuple(
        tuple(map(float, line.split())) for line in gia3d_lines
    )
    if declared_knots != len(dspec_rows):
        raise RuntimeError(f"incorrect knot count in {dspec_path}")
    if any(len(row) != 9 for row in dspec_rows):
        raise RuntimeError(f"invalid DSpecM1D row in {dspec_path}")
    if any(len(row) != 4 for row in gia3d_rows):
        raise RuntimeError(f"invalid gia3D row in {gia3d_path}")
    if tuple(row[:4] for row in dspec_rows) != gia3d_rows:
        raise RuntimeError("paired isotropic PREM deck values differ")
    for row in dspec_rows:
        if (
            row[6] != row[2]
            or row[7] != row[3]
            or row[8] != 1.0
            or not math.isfinite(row[4])
            or not math.isfinite(row[5])
        ):
            raise RuntimeError("invalid isotropic DSpecM1D columns")

    layers = []
    first = 0
    for knot in range(1, len(gia3d_rows)):
        if gia3d_rows[knot][0] == gia3d_rows[knot - 1][0]:
            layers.append(gia3d_rows[first:knot])
            first = knot
    layers.append(gia3d_rows[first:])
    if len(layers) != 12:
        raise RuntimeError(f"exported PREM has {len(layers)} layers")
    if tuple(layer[0][0] for layer in layers) != LAYER_RADII[:-1]:
        raise RuntimeError("exported PREM lower layer radii differ")
    if tuple(layer[-1][0] for layer in layers) != LAYER_RADII[1:]:
        raise RuntimeError("exported PREM upper layer radii differ")
    for layer in layers:
        if len(layer) < 4:
            raise RuntimeError("exported PREM layer has fewer than four knots")
        for lower, upper in zip(layer, layer[1:]):
            if upper[0] - lower[0] > maximum_spacing:
                raise RuntimeError("exported PREM exceeds maximum knot spacing")
    if layers[-1][-1][3] == 0.0:
        raise RuntimeError("exported PREM surface is fluid")
    return tuple(map(len, layers))


def report_difference(direction, context, left, right):
    print(
        f"comparison direction={direction} {context} "
        f"signed_difference={left-right:.17e} "
        f"relative_difference={relative_difference(left, right):.17e}"
    )


def expected_dimension(layer_elements, degree):
    radial_nodes = 4 * sum(layer_elements) + 1
    solid_nodes = (
        4 * layer_elements[0]
        + 1
        + 4 * sum(layer_elements[2:])
        + 1
    )
    return radial_nodes + 2 * solid_nodes - (degree == 1)


def verify_run(label, run, layer_elements):
    radius, gravity, elements, radial_nodes = run["metadata"]
    if run["layers"] != layer_elements:
        raise RuntimeError(
            f"{label}: layer elements {run['layers']} != {layer_elements}"
        )
    if elements != sum(layer_elements) or radial_nodes != 4 * elements + 1:
        raise RuntimeError(f"{label}: inconsistent mesh metadata")
    if not math.isclose(radius, LAYER_RADII[-1], abs_tol=1.0e-8):
        raise RuntimeError(f"{label}: surface radius is {radius}")
    return radius, gravity, elements, radial_nodes


def report_convergence(records):
    spacing = KNOT_SPACINGS[-1][0]
    for coarse, fine in zip(MESH_STEPS, MESH_STEPS[1:]):
        coarse_record = records[(spacing, coarse[0])]
        fine_record = records[(spacing, fine[0])]
        for implementation in ("analytic", "gia3d", "dspec"):
            for degree in map(int, DEGREES):
                for component, coarse_value, fine_value in zip(
                    COMPONENTS,
                    coarse_record[implementation]["results"][degree],
                    fine_record[implementation]["results"][degree],
                ):
                    print(
                        f"mesh_convergence spacing={spacing} "
                        f"step={coarse[0]}->{fine[0]} "
                        f"implementation={implementation} l={degree} "
                        f"name={component} "
                        f"signed_change={fine_value-coarse_value:.17e} "
                        f"relative_change="
                        f"{relative_difference(fine_value, coarse_value):.17e}"
                    )
        for degree in map(int, DEGREES):
            for component, dspec_coarse, gia_coarse, dspec_fine, gia_fine in zip(
                COMPONENTS,
                coarse_record["dspec"]["results"][degree],
                coarse_record["gia3d"]["results"][degree],
                fine_record["dspec"]["results"][degree],
                fine_record["gia3d"]["results"][degree],
            ):
                coarse_error = dspec_coarse - gia_coarse
                fine_error = dspec_fine - gia_fine
                print(
                    f"mesh_error_convergence spacing={spacing} "
                    f"step={coarse[0]}->{fine[0]} l={degree} "
                    f"name={component} coarse={coarse_error:.17e} "
                    f"fine={fine_error:.17e} "
                    f"signed_change={fine_error-coarse_error:.17e}"
                )

    step = MESH_STEPS[-1][0]
    for coarse, fine in zip(KNOT_SPACINGS, KNOT_SPACINGS[1:]):
        coarse_record = records[(coarse[0], step)]
        fine_record = records[(fine[0], step)]
        for implementation in ("gia3d", "dspec"):
            for degree in map(int, DEGREES):
                for component, coarse_value, fine_value, analytic in zip(
                    COMPONENTS,
                    coarse_record[implementation]["results"][degree],
                    fine_record[implementation]["results"][degree],
                    fine_record["analytic"]["results"][degree],
                ):
                    coarse_error = coarse_value - analytic
                    fine_error = fine_value - analytic
                    print(
                        f"knot_convergence step={step} "
                        f"spacing={coarse[0]}->{fine[0]} "
                        f"implementation={implementation} l={degree} "
                        f"name={component} "
                        f"coarse_analytic_error={coarse_error:.17e} "
                        f"fine_analytic_error={fine_error:.17e} "
                        f"signed_change={fine_error-coarse_error:.17e}"
                    )


def main():
    if len(sys.argv) != 5:
        raise RuntimeError(
            "usage: compare_isotropic_prem.py "
            "dspec_exe gia3d_exe exporter generated_model_directory"
        )
    dspec_executable, gia3d_executable, _, model_directory = sys.argv[1:]
    model_directory = pathlib.Path(model_directory)

    models = {}
    for label, spacing in KNOT_SPACINGS:
        dspec_path = model_directory / f"isotropic_prem_{label}.dspec"
        gia3d_path = model_directory / f"isotropic_prem_{label}.gia3d"
        layer_knots = read_exported_decks(dspec_path, gia3d_path, spacing)
        models[label] = (dspec_path, gia3d_path)
        print(
            f"model spacing={label} layers={len(layer_knots)} "
            f"layer_knots={','.join(map(str, layer_knots))} "
            f"total_knots={sum(layer_knots)} "
            f"surface_radius={LAYER_RADII[-1]:.17e} "
            f"dspec_sha256={checksum(dspec_path)} "
            f"gia3d_sha256={checksum(gia3d_path)}"
        )

    analytic_runs = {}
    records = {}
    for step_label, requested_elements in MESH_STEPS:
        analytic_runs[step_label] = run_driver(
            [
                gia3d_executable,
                "analytic-prem",
                step_label,
                *DEGREES,
            ]
        )

        for spacing_label, _ in KNOT_SPACINGS:
            dspec_path, gia3d_path = models[spacing_label]
            gia3d = run_driver(
                [
                    gia3d_executable,
                    str(gia3d_path),
                    step_label,
                    *DEGREES,
                ]
            )
            dspec = run_driver(
                [
                    dspec_executable,
                    "--selected",
                    str(dspec_path),
                    str(requested_elements),
                    *DEGREES,
                ]
            )
            analytic = analytic_runs[step_label]
            if gia3d["layers"] != analytic["layers"]:
                raise RuntimeError(
                    f"spacing={spacing_label}, step={step_label}: "
                    "analytic and sampled gia3D layer element counts differ"
                )
            records[(spacing_label, step_label)] = {
                "analytic": analytic,
                "gia3d": gia3d,
                "dspec": dspec,
            }

            metadata = {}
            for implementation, run in (
                ("analytic", analytic),
                ("gia3d", gia3d),
                ("dspec", dspec),
            ):
                metadata[implementation] = verify_run(
                    f"{spacing_label}/{step_label}/{implementation}",
                    run,
                    run["layers"],
                )
            print(
                f"mesh spacing={spacing_label} step={step_label} "
                f"analytic_layer_elements="
                f"{','.join(map(str, analytic['layers']))} "
                f"gia3d_layer_elements="
                f"{','.join(map(str, gia3d['layers']))} "
                f"dspec_layer_elements="
                f"{','.join(map(str, dspec['layers']))} "
                f"analytic_elements={metadata['analytic'][2]} "
                f"gia3d_elements={metadata['gia3d'][2]} "
                f"dspec_elements={metadata['dspec'][2]} "
                f"analytic_radial_nodes={metadata['analytic'][3]} "
                f"gia3d_radial_nodes={metadata['gia3d'][3]} "
                f"dspec_radial_nodes={metadata['dspec'][3]}"
            )
            for quantity, field in (("surface_radius", 0), ("surface_gravity", 1)):
                report_difference(
                    "dspec_sampled_minus_gia3d_sampled",
                    f"spacing={spacing_label} step={step_label} "
                    f"quantity={quantity}",
                    metadata["dspec"][field],
                    metadata["gia3d"][field],
                )
                report_difference(
                    "gia3d_sampled_minus_gia3d_analytic",
                    f"spacing={spacing_label} step={step_label} "
                    f"quantity={quantity}",
                    metadata["gia3d"][field],
                    metadata["analytic"][field],
                )
                report_difference(
                    "dspec_sampled_minus_gia3d_analytic",
                    f"spacing={spacing_label} step={step_label} "
                    f"quantity={quantity}",
                    metadata["dspec"][field],
                    metadata["analytic"][field],
                )

            for degree in map(int, DEGREES):
                for implementation, run in (
                    ("analytic", analytic),
                    ("gia3d", gia3d),
                    ("dspec", dspec),
                ):
                    diagnostic = run["diagnostics"][degree]
                    dimension = expected_dimension(run["layers"], degree)
                    if diagnostic[0] != dimension:
                        raise RuntimeError(
                            f"{implementation} l={degree}: "
                            f"DOFs {diagnostic[0]} != {dimension}"
                        )
                    print(
                        f"solve spacing={spacing_label} step={step_label} "
                        f"implementation={implementation} l={degree} "
                        f"dofs={diagnostic[0]} "
                        f"displacement_residual={diagnostic[1]:.17e} "
                        f"gravitational_residual={diagnostic[2]:.17e} "
                        f"tide_residual={diagnostic[3]:.17e} "
                        f"reciprocity_error={diagnostic[4]:.17e}"
                    )

                for component, analytic_value, gia3d_value, dspec_value in zip(
                    COMPONENTS,
                    analytic["results"][degree],
                    gia3d["results"][degree],
                    dspec["results"][degree],
                ):
                    context = (
                        f"spacing={spacing_label} step={step_label} "
                        f"l={degree} name={component}"
                    )
                    report_difference(
                        "dspec_sampled_minus_gia3d_sampled",
                        context,
                        dspec_value,
                        gia3d_value,
                    )
                    report_difference(
                        "gia3d_sampled_minus_gia3d_analytic",
                        context,
                        gia3d_value,
                        analytic_value,
                    )
                    report_difference(
                        "dspec_sampled_minus_gia3d_analytic",
                        context,
                        dspec_value,
                        analytic_value,
                    )

            for implementation, run in (
                ("analytic", analytic),
                ("gia3d", gia3d),
                ("dspec", dspec),
            ):
                for component, column in (
                    ("k_u", 1),
                    ("h_phi", 2),
                    ("k_phi", 3),
                    ("k_load", 5),
                    ("k_t", 7),
                ):
                    value = run["results"][1][column]
                    print(
                        f"degree_one_zero spacing={spacing_label} "
                        f"step={step_label} implementation={implementation} "
                        f"component={component} value={value:.17e}"
                    )
                    if value != 0.0:
                        raise RuntimeError(
                            f"{implementation} degree-one {component} "
                            "is not exactly zero"
                        )

    report_convergence(records)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
