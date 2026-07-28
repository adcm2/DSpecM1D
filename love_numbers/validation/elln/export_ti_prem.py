#!/usr/bin/env python3

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True

from validate_elln import audit_source, prepare_archive, read_numeric_rows


ELLN_GRAVITATIONAL_CONSTANT = 6.67384e-11
SPACINGS = (
    ("100km", 100_000.0),
    ("50km", 50_000.0),
    ("25km", 25_000.0),
    ("12.5km", 12_500.0),
    ("6.25km", 6_250.0),
)


def checksum(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sampled_radii(lower, upper, maximum_spacing):
    intervals = max(3, math.ceil((upper - lower) / maximum_spacing))
    return tuple(
        lower + (upper - lower) * point / intervals
        for point in range(intervals + 1)
    )


def read_layers(paths):
    mantle_rows = read_numeric_rows(paths["EarthMantleTI56.txt"], 8)
    core_rows = read_numeric_rows(paths["EarthCore26.txt"], 4)
    layers = []

    lower = 0.0
    for _, upper, density, bulk_modulus in core_rows:
        def core_values(radius, density=density,
                        bulk_modulus=bulk_modulus):
            del radius
            p_speed = math.sqrt(bulk_modulus / density)
            return density, p_speed, 0.0, p_speed, 0.0, 1.0

        layers.append(("fluid", lower, upper, core_values))
        lower = upper

    for lower_row, upper_row in zip(mantle_rows, mantle_rows[1:]):
        lower = lower_row[1]
        upper = upper_row[1]
        density_at_upper = upper_row[2]
        a, c, l, n, f = upper_row[3:]
        equivalent_radius = (
            2.0 * (upper**3 - lower**3)
            / (3.0 * (upper**2 - lower**2))
        )

        def mantle_values(
            radius, density_at_upper=density_at_upper,
            equivalent_radius=equivalent_radius,
            a=a, c=c, l=l, n=n, f=f,
        ):
            density = density_at_upper * equivalent_radius / radius
            eta = f / (a - 2.0 * l)
            return (
                density,
                math.sqrt(c / density),
                math.sqrt(l / density),
                math.sqrt(a / density),
                math.sqrt(n / density),
                eta,
            )

        layers.append(("solid", lower, upper, mantle_values))

    if len(layers) != 83:
        raise RuntimeError(f"Expected 83 physical regions, got {len(layers)}")
    return layers


def write_deck(path, layers, maximum_spacing):
    rows = []
    layer_knots = []
    for _, lower, upper, values in layers:
        radii = sampled_radii(lower, upper, maximum_spacing)
        layer_knots.append(len(radii))
        for radius in radii:
            density, vpv, vsv, vph, vsh, eta = values(radius)
            row = (
                radius, density, vpv, vsv, 1000.0, 600.0,
                vph, vsh, eta,
            )
            if not all(math.isfinite(value) for value in row):
                raise RuntimeError(f"Non-finite deck value at radius {radius}")
            rows.append(row)

    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(
            "official ELLN EarthMantleTI56 and EarthCore26 "
            "sampled for DSpecM1D validation\n"
        )
        stream.write("1 -1.0 1 2\n")
        stream.write(f"{len(rows)} 0 0\n")
        for row in rows:
            stream.write(" ".join(f"{value:.17e}" for value in row))
            stream.write("\n")
    return layer_knots


def density_moment(layers):
    total = 0.0
    for kind, lower, upper, values in layers:
        if kind == "fluid":
            density = values(lower)[0]
            total += density * (upper**3 - lower**3) / 3.0
        else:
            total += (
                values(upper)[0]
                * upper
                * (upper**2 - lower**2)
                / 2.0
            )
    return total


def read_expected_checksums(path):
    expected = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        digest, name = line.split()
        expected[name] = digest
    return expected


def parse_arguments():
    parser = argparse.ArgumentParser()
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--archive", type=Path)
    source.add_argument("--source-dir", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-checksums", type=Path)
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    arguments.output_dir.mkdir(parents=True, exist_ok=True)
    if arguments.archive:
        source_directory = prepare_archive(
            arguments.archive.resolve(),
            arguments.output_dir / "source-audit",
        )
    else:
        source_directory = arguments.source_dir.resolve()
    paths = audit_source(source_directory)
    layers = read_layers(paths)
    surface_radius = layers[-1][2]
    moment = density_moment(layers)

    models = []
    for label, spacing in SPACINGS:
        path = arguments.output_dir / f"elln_ti_prem_{label}.dspec"
        layer_knots = write_deck(path, layers, spacing)
        models.append({
            "file": path.name,
            "label": label,
            "layer_knots": layer_knots,
            "maximum_knot_spacing": spacing,
            "sha256": checksum(path),
            "total_knots": sum(layer_knots),
        })

    if arguments.expected_checksums:
        expected = read_expected_checksums(arguments.expected_checksums)
        actual = {model["file"]: model["sha256"] for model in models}
        if actual != expected:
            raise RuntimeError(
                f"Generated model checksums differ: {actual} != {expected}"
            )

    manifest = {
        "elln_gravitational_constant": ELLN_GRAVITATIONAL_CONSTANT,
        "elln_surface_gravity": (
            4.0 * math.pi * ELLN_GRAVITATIONAL_CONSTANT
            * moment / surface_radius**2
        ),
        "enclosed_density_moment": moment,
        "fluid_regions": 27,
        "models": models,
        "physical_regions": len(layers),
        "solid_regions": 56,
        "surface_radius": surface_radius,
    }
    manifest_path = arguments.output_dir / "elln_ti_prem_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    for model in models:
        print(
            f"ELLN TI PREM model {model['label']}: "
            f"knots={model['total_knots']} sha256={model['sha256']}"
        )
    print(
        f"ELLN TI PREM: regions={len(layers)} "
        f"surface_radius={surface_radius:.17e} "
        f"surface_gravity={manifest['elln_surface_gravity']:.17e}"
    )


if __name__ == "__main__":
    main()
