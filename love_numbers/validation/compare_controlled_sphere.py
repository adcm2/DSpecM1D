#!/usr/bin/env python3

import math
import subprocess
import sys


COMPONENTS = (
    "h_u",
    "k_u",
    "h_phi",
    "k_phi",
    "h_load",
    "k_load",
    "h_t",
    "k_t",
)
DEGREES = ("1", "2", "10")
RADIUS = 6_371_000.0
GRAVITATIONAL_CONSTANT = 6.67230e-11

MODELS = {
    "homogeneous": {
        "elements": (4, 8, 16, 32),
        "layer_fractions": (1.0,),
        "solid_layers": (0,),
        "rows": (
            (0.0, 5500.0, 10000.0, 5000.0),
            (RADIUS / 3.0, 5500.0, 10000.0, 5000.0),
            (2.0 * RADIUS / 3.0, 5500.0, 10000.0, 5000.0),
            (RADIUS, 5500.0, 10000.0, 5000.0),
        ),
        "density_integral": 5500.0 / 3.0,
        "fluids": {},
        "interfaces": (),
    },
    "solid-fluid-solid": {
        "elements": (8, 16, 32),
        "layer_fractions": (0.25, 0.25, 0.5),
        "solid_layers": (0, 2),
        "rows": (
            (0.0, 7000.0, 11000.0, 5500.0),
            (RADIUS / 8.0, 7000.0, 11000.0, 5500.0),
            (RADIUS / 4.0, 7000.0, 11000.0, 5500.0),
            (RADIUS / 4.0, 6000.0, 9000.0, 0.0),
            (3.0 * RADIUS / 8.0, 5500.0, 9000.0, 0.0),
            (RADIUS / 2.0, 5000.0, 9000.0, 0.0),
            (RADIUS / 2.0, 4500.0, 9000.0, 4500.0),
            (3.0 * RADIUS / 4.0, 4500.0, 9000.0, 4500.0),
            (RADIUS, 4500.0, 9000.0, 4500.0),
        ),
        "density_integral": (
            7000.0 * (0.25**3) / 3.0
            + (
                7000.0 * (0.5**3 - 0.25**3) / 3.0
                - 1000.0 * (0.5**4 - 0.25**4)
            )
            + 4500.0 * (1.0 - 0.5**3) / 3.0
        ),
        "fluids": {
            1: (6000.0, 5000.0, -4000.0 / RADIUS),
        },
        "interfaces": (
            (RADIUS / 4.0, 6000.0, -1),
            (RADIUS / 2.0, 5000.0, 1),
        ),
    },
    "central-fluid": {
        "elements": (8, 16, 32),
        "layer_fractions": (0.5, 0.5),
        "solid_layers": (1,),
        "rows": (
            (0.0, 6000.0, 9000.0, 0.0),
            (RADIUS / 4.0, 6000.0, 9000.0, 0.0),
            (RADIUS / 2.0, 6000.0, 9000.0, 0.0),
            (RADIUS / 2.0, 4500.0, 9000.0, 4500.0),
            (3.0 * RADIUS / 4.0, 4500.0, 9000.0, 4500.0),
            (RADIUS, 4500.0, 9000.0, 4500.0),
        ),
        "density_integral": (
            6000.0 * (0.5**3) / 3.0
            + 4500.0 * (1.0 - 0.5**3) / 3.0
        ),
        "fluids": {
            0: (6000.0, 6000.0, 0.0),
        },
        "interfaces": (
            (RADIUS / 2.0, 6000.0, 1),
        ),
    },
}


def relative_difference(value, reference):
    if reference != 0.0:
        return (value - reference) / abs(reference)
    return 0.0 if value == reference else math.copysign(math.inf, value)


def read_decks(dspec_path, gia_path):
    with open(dspec_path, encoding="utf-8") as stream:
        dspec_lines = stream.readlines()[3:]
    with open(gia_path, encoding="utf-8") as stream:
        gia_lines = stream.readlines()

    dspec_rows = []
    for line in dspec_lines:
        fields = tuple(map(float, line.split()))
        if len(fields) != 9:
            raise RuntimeError(f"invalid DSpecM1D deck row: {line}")
        if (
            fields[6] != fields[2]
            or fields[7] != fields[3]
            or fields[8] != 1.0
            or fields[4] <= 0.0
            or fields[5] <= 0.0
        ):
            raise RuntimeError(f"invalid isotropic DSpecM1D row: {line}")
        dspec_rows.append(fields[:4])

    gia_rows = []
    for line in gia_lines:
        fields = tuple(map(float, line.split()))
        if len(fields) != 4:
            raise RuntimeError(f"invalid gia3D deck row: {line}")
        gia_rows.append(fields)
    return tuple(dspec_rows), tuple(gia_rows)


def parse_flag(text):
    if text in ("0", "F"):
        return False
    if text in ("1", "T"):
        return True
    raise RuntimeError(f"invalid floating-point flag: {text}")


def run_driver(command):
    completed = subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        sys.stderr.write(completed.stdout)
        sys.stderr.write(completed.stderr)
        raise RuntimeError(
            f"driver failed with exit code {completed.returncode}: "
            f"{' '.join(command)}"
        )

    metadata = None
    layer_elements = None
    floating_point_flags = None
    legacy_surface_gravity = None
    diagnostics = {}
    fluids = {}
    for line in completed.stderr.splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "metadata":
            metadata = (
                float(fields[1]),
                float(fields[2]),
                int(fields[3]),
                int(fields[4]),
            )
        elif fields[0] == "layers":
            number_of_layers = int(fields[1])
            layer_elements = tuple(map(int, fields[2:]))
            if len(layer_elements) != number_of_layers:
                raise RuntimeError(f"invalid layer row: {line}")
        elif fields[0] == "fluid":
            fluids[int(fields[1])] = tuple(map(float, fields[2:]))
        elif fields[0] == "diagnostic":
            diagnostics[int(fields[1])] = (
                int(fields[2]),
                float(fields[3]),
                float(fields[4]),
                float(fields[5]),
            )
        elif fields[0] == "floating_point_flags":
            floating_point_flags = (
                parse_flag(fields[1]),
                parse_flag(fields[2]),
            )
        elif fields[0] == "legacy_surface_gravity":
            legacy_surface_gravity = float(fields[1])

    results = {}
    for line in completed.stdout.splitlines():
        fields = line.split()
        if len(fields) != 9:
            raise RuntimeError(f"invalid result row: {line}")
        results[int(fields[0])] = tuple(map(float, fields[1:]))

    if metadata is None or layer_elements is None:
        raise RuntimeError("driver did not report complete mesh metadata")
    if floating_point_flags is None:
        raise RuntimeError("driver did not report floating-point flags")
    values = metadata[:2]
    if legacy_surface_gravity is not None:
        values += (legacy_surface_gravity,)
    values += tuple(value for row in diagnostics.values() for value in row[1:])
    values += tuple(value for row in fluids.values() for value in row)
    values += tuple(value for row in results.values() for value in row)
    if not all(map(math.isfinite, values)):
        raise RuntimeError("driver reported a non-finite value")
    if any(floating_point_flags):
        raise RuntimeError(
            "driver raised invalid or divide-by-zero floating-point flags"
        )
    return {
        "metadata": metadata,
        "layers": layer_elements,
        "fluids": fluids,
        "diagnostics": diagnostics,
        "results": results,
        "floating_point_flags": floating_point_flags,
        "legacy_surface_gravity": legacy_surface_gravity,
    }


def print_difference(prefix, dspec_value, gia_value):
    difference = dspec_value - gia_value
    relative = relative_difference(dspec_value, gia_value)
    print(
        f"{prefix} dspec={dspec_value:.17e} gia3d={gia_value:.17e} "
        f"signed_difference={difference:.17e} "
        f"relative_difference={relative:.17e}"
    )


def expected_dimension(model, layer_elements, degree):
    radial_nodes = 4 * sum(layer_elements) + 1
    solid_nodes = sum(
        4 * layer_elements[layer] + 1
        for layer in model["solid_layers"]
    )
    return radial_nodes + 2 * solid_nodes - (degree == 1)


def main():
    if len(sys.argv) not in (5, 6):
        raise RuntimeError(
            "usage: compare_controlled_sphere.py "
            "dspec_exe gia_exe dspec_deck gia_deck [model]"
        )
    (
        dspec_executable,
        gia_executable,
        dspec_model,
        gia_model,
    ) = sys.argv[1:5]
    model_name = sys.argv[5] if len(sys.argv) == 6 else "homogeneous"
    model = MODELS[model_name]

    dspec_rows, gia_rows = read_decks(dspec_model, gia_model)
    if dspec_rows != gia_rows or dspec_rows != model["rows"]:
        raise RuntimeError(f"{model_name}: paired deck values do not match")
    print(f"deck model={model_name} rows={len(dspec_rows)} verified=true")

    for radius, density, orientation in model["interfaces"]:
        print(
            f"interface model={model_name} radius={radius:.17e} "
            f"fluid_density={density:.17e} orientation={orientation:+d}"
        )
    for layer, fluid in model["fluids"].items():
        print(
            f"fluid_definition model={model_name} layer={layer} "
            f"lower_density={fluid[0]:.17e} "
            f"upper_density={fluid[1]:.17e} "
            f"density_derivative={fluid[2]:.17e}"
        )

    analytic_gravity = (
        4.0
        * math.pi
        * GRAVITATIONAL_CONSTANT
        * RADIUS
        * model["density_integral"]
    )
    mesh_results = {"dspec": {}, "gia3d": {}}

    for elements in model["elements"]:
        dspec = run_driver(
            [dspec_executable, dspec_model, str(elements), *DEGREES]
        )
        gia = run_driver(
            [
                gia_executable,
                gia_model,
                format(1.0 / elements, ".17g"),
                *DEGREES,
            ]
        )
        mesh_results["dspec"][elements] = dspec["results"]
        mesh_results["gia3d"][elements] = gia["results"]

        dspec_radius, dspec_gravity, dspec_elements, dspec_nodes = (
            dspec["metadata"]
        )
        gia_radius, gia_gravity, gia_elements, gia_nodes = gia["metadata"]
        expected_layers = tuple(
            int(elements * fraction)
            for fraction in model["layer_fractions"]
        )
        if dspec["layers"] != expected_layers:
            raise RuntimeError(
                f"N={elements}: DSpecM1D layer counts "
                f"{dspec['layers']} != {expected_layers}"
            )
        if gia["layers"] != expected_layers:
            raise RuntimeError(
                f"N={elements}: gia3D layer counts "
                f"{gia['layers']} != {expected_layers}"
            )
        if (dspec_elements, gia_elements) != (elements, elements):
            raise RuntimeError(
                f"N={elements}: actual element counts are "
                f"{dspec_elements} and {gia_elements}"
            )
        expected_nodes = 4 * elements + 1
        if (dspec_nodes, gia_nodes) != (expected_nodes, expected_nodes):
            raise RuntimeError(
                f"N={elements}: actual radial-node counts are "
                f"{dspec_nodes} and {gia_nodes}"
            )

        print_difference(
            f"model name={model_name} N={elements} "
            "quantity=surface_radius",
            dspec_radius,
            gia_radius,
        )
        print_difference(
            f"model name={model_name} N={elements} "
            "quantity=surface_gravity",
            dspec_gravity,
            gia_gravity,
        )
        for implementation, gravity in (
            ("dspec", dspec_gravity),
            ("gia3d", gia_gravity),
        ):
            print(
                f"gravity model={model_name} N={elements} "
                f"implementation={implementation} "
                f"calculated={gravity:.17e} "
                f"analytic={analytic_gravity:.17e} "
                f"signed_difference={gravity-analytic_gravity:.17e} "
                f"relative_difference="
                f"{relative_difference(gravity, analytic_gravity):.17e}"
            )
        if dspec["legacy_surface_gravity"] is not None:
            legacy_gravity = dspec["legacy_surface_gravity"]
            print(
                f"gravity model={model_name} N={elements} "
                "implementation=dspec_legacy "
                f"calculated={legacy_gravity:.17e} "
                f"analytic={analytic_gravity:.17e} "
                f"signed_difference="
                f"{legacy_gravity-analytic_gravity:.17e} "
                f"relative_difference="
                f"{relative_difference(legacy_gravity, analytic_gravity):.17e}"
            )
        print(
            f"mesh model={model_name} N={elements} "
            f"layer_elements={','.join(map(str, expected_layers))} "
            f"dspec_elements={dspec_elements} "
            f"gia3d_elements={gia_elements} "
            f"dspec_radial_nodes={dspec_nodes} "
            f"gia3d_radial_nodes={gia_nodes}"
        )
        print(
            f"floating_point model={model_name} N={elements} "
            "dspec_invalid=0 dspec_divide_by_zero=0 "
            "gia3d_invalid=0 gia3d_divide_by_zero=0"
        )

        if set(dspec["fluids"]) != set(model["fluids"]):
            raise RuntimeError("DSpecM1D fluid-layer diagnostics are incomplete")
        if set(gia["fluids"]) != set(model["fluids"]):
            raise RuntimeError("gia3D fluid-layer diagnostics are incomplete")
        for layer, expected in model["fluids"].items():
            for implementation, actual in (
                ("dspec", dspec["fluids"][layer]),
                ("gia3d", gia["fluids"][layer]),
            ):
                for quantity, value, reference in zip(
                    (
                        "lower_density",
                        "upper_density",
                        "minimum_density_derivative",
                        "maximum_density_derivative",
                    ),
                    actual,
                    (expected[0], expected[1], expected[2], expected[2]),
                ):
                    print(
                        f"fluid_sample model={model_name} N={elements} "
                        f"implementation={implementation} layer={layer} "
                        f"quantity={quantity} value={value:.17e} "
                        f"expected={reference:.17e} "
                        f"signed_difference={value-reference:.17e} "
                        f"relative_difference="
                        f"{relative_difference(value, reference):.17e}"
                    )

        for degree_text in DEGREES:
            degree = int(degree_text)
            dspec_dimension, *dspec_residuals = dspec["diagnostics"][degree]
            gia_dimension, *gia_residuals = gia["diagnostics"][degree]
            dimension = expected_dimension(model, expected_layers, degree)
            if (dspec_dimension, gia_dimension) != (dimension, dimension):
                raise RuntimeError(
                    f"N={elements}, l={degree}: dimensions are "
                    f"{dspec_dimension} and {gia_dimension}; "
                    f"expected {dimension}"
                )
            print(
                f"matrix model={model_name} N={elements} l={degree} "
                f"dspec_dofs={dspec_dimension} gia3d_dofs={gia_dimension}"
            )
            for column, dspec_residual, gia_residual in zip(
                ("displacement", "gravitational", "tide"),
                dspec_residuals,
                gia_residuals,
            ):
                print(
                    f"residual model={model_name} N={elements} "
                    f"l={degree} column={column} "
                    f"dspec={dspec_residual:.17e} "
                    f"gia3d={gia_residual:.17e}"
                )

            for component, dspec_value, gia_value in zip(
                COMPONENTS,
                dspec["results"][degree],
                gia["results"][degree],
            ):
                print_difference(
                    f"component model={model_name} N={elements} "
                    f"l={degree} name={component}",
                    dspec_value,
                    gia_value,
                )

        for implementation, results in (
            ("dspec", dspec["results"]),
            ("gia3d", gia["results"]),
        ):
            degree_one = results[1]
            for component, column in (
                ("k_u", 1),
                ("h_phi", 2),
                ("k_phi", 3),
                ("k_load", 5),
                ("k_t", 7),
            ):
                value = degree_one[column]
                print(
                    f"degree_one_zero model={model_name} N={elements} "
                    f"implementation={implementation} "
                    f"component={component} value={value:.17e}"
                )
                if value != 0.0:
                    raise RuntimeError(
                        f"N={elements}: {implementation} degree-one "
                        f"{component} is not exactly zero"
                    )

    for implementation, results_by_mesh in mesh_results.items():
        for coarse, fine in zip(model["elements"], model["elements"][1:]):
            for degree_text in DEGREES:
                degree = int(degree_text)
                for component, coarse_value, fine_value in zip(
                    COMPONENTS,
                    results_by_mesh[coarse][degree],
                    results_by_mesh[fine][degree],
                ):
                    difference = fine_value - coarse_value
                    relative = relative_difference(fine_value, coarse_value)
                    print(
                        f"convergence model={model_name} "
                        f"implementation={implementation} "
                        f"N={coarse}->{fine} l={degree} "
                        f"name={component} "
                        f"signed_change={difference:.17e} "
                        f"relative_change={relative:.17e}"
                    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
