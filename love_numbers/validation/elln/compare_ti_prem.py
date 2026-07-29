#!/usr/bin/env python3

import hashlib
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from compare_controlled_sphere import relative_difference, run_driver
from validate_elln import read_reference


DSPEC_GRAVITATIONAL_CONSTANT = 6.6723e-11
ELLN_GRAVITATIONAL_CONSTANT = 6.67384e-11
KNOT_SPACINGS = ("100km", "50km", "25km", "12.5km", "6.25km")
MESH_STEPS = ("0.005", "0.0025", "0.00125")


def checksum(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_checksums(path):
    checksums = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        digest, name = line.split()
        checksums[name] = digest
    return checksums


def read_model(path, maximum_spacing):
    lines = path.read_text(encoding="utf-8").splitlines()
    declared_knots = int(lines[2].split()[0])
    rows = tuple(tuple(map(float, line.split())) for line in lines[3:])
    if len(rows) != declared_knots or any(len(row) != 9 for row in rows):
        raise RuntimeError(f"Invalid DSpecM1D deck: {path}")

    layers = []
    first = 0
    for knot in range(1, len(rows)):
        if rows[knot][0] == rows[knot - 1][0]:
            layers.append(rows[first:knot])
            first = knot
    layers.append(rows[first:])
    if len(layers) != 83:
        raise RuntimeError(f"{path} has {len(layers)} physical regions")

    for layer_index, layer in enumerate(layers):
        if len(layer) < 4:
            raise RuntimeError("Each ELLN region requires at least four knots")
        for lower, upper in zip(layer, layer[1:]):
            if upper[0] - lower[0] > maximum_spacing:
                raise RuntimeError(f"{path} exceeds its knot spacing")
        for row in layer:
            if not all(math.isfinite(value) for value in row):
                raise RuntimeError(f"{path} contains a non-finite value")
            is_fluid = row[3] == 0.0 and row[7] == 0.0
            if is_fluid != (layer_index < 27):
                raise RuntimeError("ELLN fluid-core topology was not preserved")
    return tuple(map(len, layers))


def conventional_result(result, degree, gravity, radius):
    factor = (
        (2 * degree + 1)
        / (4.0 * math.pi * DSPEC_GRAVITATIONAL_CONSTANT * radius)
    )
    return (
        factor * gravity * result[4],
        -(1.0 + factor * result[5]),
    )


def report_comparison(context, name, value, reference):
    print(
        f"comparison {context} name={name} "
        f"dspec={value:.17e} official_elln={reference:.17e} "
        f"signed_difference={value-reference:.17e} "
        f"relative_difference={relative_difference(value, reference):.17e}"
    )


def expected_dimension(run, degree):
    elements = int(run["metadata"][2])
    radial_nodes = 6 * elements + 1
    solid_elements = sum(run["layers"][27:])
    solid_nodes = 6 * solid_elements + 1
    return radial_nodes + 2 * solid_nodes - (degree == 1)


def report_convergence(records):
    spacing = KNOT_SPACINGS[-1]
    for coarse, fine in zip(MESH_STEPS, MESH_STEPS[1:]):
        for name, component in (("h", 0), ("k", 1)):
            coarse_value = records[(spacing, coarse)]["conventional"][component]
            fine_value = records[(spacing, fine)]["conventional"][component]
            print(
                f"mesh_convergence spacing={spacing} step={coarse}->{fine} "
                f"l=10 name={name} "
                f"signed_change={fine_value-coarse_value:.17e} "
                f"relative_change="
                f"{relative_difference(fine_value, coarse_value):.17e}"
            )

    step = MESH_STEPS[-1]
    for coarse, fine in zip(KNOT_SPACINGS, KNOT_SPACINGS[1:]):
        for name, component in (("h", 0), ("k", 1)):
            coarse_value = records[(coarse, step)]["conventional"][component]
            fine_value = records[(fine, step)]["conventional"][component]
            print(
                f"knot_convergence step={step} spacing={coarse}->{fine} "
                f"l=10 name={name} "
                f"signed_change={fine_value-coarse_value:.17e} "
                f"relative_change="
                f"{relative_difference(fine_value, coarse_value):.17e}"
            )


def main():
    if len(sys.argv) != 6:
        raise RuntimeError(
            "usage: compare_ti_prem.py dspec_exe model_dir manifest "
            "checksums official_table"
        )
    executable = sys.argv[1]
    model_directory = Path(sys.argv[2])
    manifest = json.loads(Path(sys.argv[3]).read_text(encoding="utf-8"))
    expected_checksums = read_checksums(Path(sys.argv[4]))
    references = read_reference(Path(sys.argv[5]))

    radius = manifest["surface_radius"]
    density_moment = manifest["enclosed_density_moment"]
    expected_dspec_gravity = (
        4.0 * math.pi * DSPEC_GRAVITATIONAL_CONSTANT
        * density_moment / radius**2
    )
    relative_g_difference = (
        (DSPEC_GRAVITATIONAL_CONSTANT - ELLN_GRAVITATIONAL_CONSTANT)
        / ELLN_GRAVITATIONAL_CONSTANT
    )
    print(
        "constants "
        f"G_ELLN={ELLN_GRAVITATIONAL_CONSTANT:.17e} "
        f"G_DSpecM1D={DSPEC_GRAVITATIONAL_CONSTANT:.17e} "
        f"relative_difference={relative_g_difference:.17e}"
    )
    print(
        f"analytic_gravity elln={manifest['elln_surface_gravity']:.17e} "
        f"dspec={expected_dspec_gravity:.17e} "
        f"signed_difference="
        f"{expected_dspec_gravity-manifest['elln_surface_gravity']:.17e}"
    )

    models = {}
    for model in manifest["models"]:
        label = model["label"]
        path = model_directory / model["file"]
        actual_checksum = checksum(path)
        if actual_checksum != expected_checksums[path.name]:
            raise RuntimeError(f"Checksum mismatch for {path}")
        layer_knots = read_model(path, model["maximum_knot_spacing"])
        if layer_knots != tuple(model["layer_knots"]):
            raise RuntimeError(f"Manifest knot counts differ for {path}")
        models[label] = path
        print(
            f"model spacing={label} regions={len(layer_knots)} "
            f"total_knots={sum(layer_knots)} "
            f"surface_radius={radius:.17e} sha256={actual_checksum}"
        )

    official_h = -references[10][0]
    official_k = -references[10][2]
    official_h_one = -references[1][0]
    official_k_one = -references[1][2]
    records = {}
    for spacing in KNOT_SPACINGS:
        for step in MESH_STEPS:
            run = run_driver([
                executable, "--configured", str(models[spacing]),
                "6", step, "1", "10",
            ])
            surface_radius, gravity, elements, radial_nodes = run["metadata"]
            if surface_radius != radius:
                raise RuntimeError("DSpecM1D surface radius changed")
            if radial_nodes != 6 * elements + 1:
                raise RuntimeError("Unexpected p=6 radial-node count")

            conventional = conventional_result(
                run["results"][10], 10, gravity, surface_radius
            )
            records[(spacing, step)] = {
                "conventional": conventional,
                "run": run,
            }
            print(
                f"mesh spacing={spacing} step={step} "
                f"elements={elements} radial_nodes={radial_nodes} "
                f"surface_radius={surface_radius:.17e} "
                f"surface_gravity={gravity:.17e} "
                f"analytic_dspec_gravity={expected_dspec_gravity:.17e} "
                f"gravity_error={gravity-expected_dspec_gravity:.17e}"
            )

            for degree in (1, 10):
                diagnostic = run["diagnostics"][degree]
                dimension = expected_dimension(run, degree)
                if diagnostic[0] != dimension:
                    raise RuntimeError(
                        f"l={degree}: DOFs {diagnostic[0]} != {dimension}"
                    )
                print(
                    f"solve spacing={spacing} step={step} l={degree} "
                    f"dofs={diagnostic[0]} "
                    f"displacement_residual={diagnostic[1]:.17e} "
                    f"gravitational_residual={diagnostic[2]:.17e} "
                    f"tide_residual={diagnostic[3]:.17e} "
                    f"reciprocity_error={diagnostic[4]:.17e}"
                )

            context = f"spacing={spacing} step={step} l=10"
            report_comparison(context, "h", conventional[0], official_h)
            report_comparison(context, "k", conventional[1], official_k)

            degree_one = run["results"][1]
            degree_one_conventional = conventional_result(
                degree_one, 1, gravity, surface_radius
            )
            if degree_one[5] != 0.0 or degree_one_conventional[1] != -1.0:
                raise RuntimeError("Unexpected DSpecM1D degree-one potential")
            print(
                f"degree_one_excluded spacing={spacing} step={step} "
                f"dspec_generic_h={degree_one_conventional[0]:.17e} "
                f"dspec_generic_k={degree_one_conventional[1]:.17e} "
                f"official_elln_h={official_h_one:.17e} "
                f"official_elln_k={official_k_one:.17e} "
                "reason=frame_conventions_not_equivalent"
            )

    report_convergence(records)
    print(
        "horizontal_response excluded: DSpecM1D public API does not "
        "expose surface horizontal displacement"
    )


if __name__ == "__main__":
    main()
