#!/usr/bin/env python3
"""Generate the hash-bound Study 1 physical probe grid from case geometry."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sphere_unit(theta: float, phi: float) -> tuple[float, float, float]:
    return (math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta))


def dot(a: tuple[float, ...], b: tuple[float, ...]) -> float:
    return sum(x * y for x, y in zip(a, b))


def scale(a: tuple[float, ...], x: float) -> tuple[float, ...]:
    return tuple(v * x for v in a)


def add(a: tuple[float, ...], b: tuple[float, ...]) -> tuple[float, ...]:
    return tuple(x + y for x, y in zip(a, b))


def cross(a: tuple[float, ...], b: tuple[float, ...]) -> tuple[float, float, float]:
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def normalise(a: tuple[float, ...]) -> tuple[float, float, float]:
    n = math.sqrt(dot(a, a))
    if not n:
        raise ValueError("cannot normalize a zero vector")
    return tuple(x / n for x in a)


def path_geometry(source_theta: float, source_phi: float,
                  receiver_theta: float, receiver_phi: float):
    src = sphere_unit(source_theta, source_phi)
    receiver = sphere_unit(receiver_theta, receiver_phi)
    along = dot(src, receiver)
    travel = normalise(tuple(receiver[i] - along * src[i] for i in range(3)))
    normal = cross(src, travel)
    measured_distance = math.degrees(math.atan2(dot(receiver, travel), along))

    def at(along_deg: float, across_deg: float = 0.0) -> tuple[float, float]:
        a, c = math.radians(along_deg), math.radians(across_deg)
        point = tuple(math.cos(c) * (src[i] * math.cos(a) + travel[i] * math.sin(a)) + math.sin(c) * normal[i]
                      for i in range(3))
        return math.acos(max(-1.0, min(1.0, point[2]))), math.atan2(point[1], point[0])

    return receiver, measured_distance, at


def read_path_endpoints(path: Path, expected_source_theta: float, expected_source_phi: float,
                        expected_distance_deg: float, expected_azimuth_deg: float):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows or list(rows[0]) != ["colatitude_rad", "longitude_rad", "along_deg", "across_deg", "angular_weight_sr"]:
        raise ValueError("frozen path-angle receipt has an unexpected schema")
    source_theta = float(rows[0]["colatitude_rad"])
    source_phi = float(rows[0]["longitude_rad"])
    receiver_theta = float(rows[-1]["colatitude_rad"])
    receiver_phi = float(rows[-1]["longitude_rad"])
    if abs(source_theta - expected_source_theta) > 2e-14 or abs(math.remainder(source_phi - expected_source_phi, 2 * math.pi)) > 2e-14:
        raise ValueError("frozen path-angle source does not match the case through YSpec Sphere::Colatitude")
    if abs(float(rows[0]["along_deg"])) > 1e-12 or abs(float(rows[-1]["along_deg"]) - expected_distance_deg) > 1e-10:
        raise ValueError("frozen path-angle endpoints do not match the configured source-receiver distance")
    source = sphere_unit(source_theta, source_phi)
    receiver = sphere_unit(receiver_theta, receiver_phi)
    along = dot(source, receiver)
    travel = normalise(tuple(receiver[i] - along * source[i] for i in range(3)))
    north = (-math.cos(source_theta) * math.cos(source_phi),
             -math.cos(source_theta) * math.sin(source_phi), math.sin(source_theta))
    east = (-math.sin(source_phi), math.cos(source_phi), 0.0)
    azimuth = math.degrees(math.atan2(dot(travel, east), dot(travel, north))) % 360.0
    if abs(azimuth - expected_azimuth_deg) > 1e-10:
        raise ValueError("frozen path-angle endpoint does not match configured geocentric initial azimuth")
    return source_theta, source_phi, receiver_theta, receiver_phi


def load_mesh(path: Path):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows or list(rows[0]) != ["element", "layer", "lower_radius_m", "upper_radius_m", "fluid"]:
        raise ValueError("baseline mesh receipt has an unexpected schema")
    return [{"element": int(row["element"]), "layer": int(row["layer"]),
             "lo": float(row["lower_radius_m"]), "hi": float(row["upper_radius_m"]),
             "fluid": row["fluid"] == "1"} for row in rows]


def resolve(radius_m: float, side: str, trace: str, mesh: list[dict]):
    tol = max(1e-7, 64 * math.ulp(radius_m))
    candidates = [e for e in mesh if e["fluid"] == (side == "fluid")]
    if trace == "inner":
        found = [e for e in candidates if abs(e["hi"] - radius_m) <= tol]
    elif trace == "outer":
        found = [e for e in candidates if abs(e["lo"] - radius_m) <= tol]
    else:
        found = [e for e in candidates if e["lo"] + tol < radius_m < e["hi"] - tol]
        if not found:
            found = [e for e in candidates if abs(e["hi"] - radius_m) <= tol]
    if len(found) != 1:
        raise ValueError(f"expected exactly one baseline element for r={radius_m} side={side} trace={trace}; got {len(found)}")
    return found[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--path-reference", type=Path, required=True)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    case = json.loads(args.case.read_text(encoding="utf-8"))
    source = case["sources"][0]
    source_theta = math.pi / 2 - math.atan(0.993277 * math.tan(math.radians(source["latitude"])))
    source_phi = math.radians(source["longitude"])
    source_theta, source_phi, receiver_theta, receiver_phi = read_path_endpoints(
        args.path_reference, source_theta, source_phi,
        case["receiver"]["distance_deg"], case["receiver"]["azimuth_deg"])
    _, distance_deg, at = path_geometry(source_theta, source_phi, receiver_theta, receiver_phi)
    mesh = load_mesh(args.mesh)
    surface = max(element["hi"] for element in mesh)
    source_radius = surface - source["depth_km"] * 1000.0
    receiver_radius = surface - case["receiver"]["depth_km"] * 1000.0
    boundary_offset = 0.44776119402985076
    alongs = [0.0, boundary_offset, 1.0, 5.0, 15.0, 30.0, 45.0, 55.0,
              distance_deg - boundary_offset, distance_deg]
    cross_offsets = [-4.0, -2.0, 2.0, 4.0]
    angles: list[tuple[str, float, float, str]] = []
    for i, along in enumerate(alongs):
        theta, phi = at(along)
        angles.append((f"path-{i:02d}", theta, phi, f"along={along:.15g};across=0"))
    for i, across in enumerate(cross_offsets):
        theta, phi = at(distance_deg / 2, across)
        angles.append((f"cross-{i:02d}", theta, phi,
                       f"along={distance_deg / 2:.15g};across={across:.15g}"))

    rows: list[dict[str, object]] = []
    for angle_id, theta, phi, construction in angles:
        for depth_km in (3.0, 10.0, 24.0, 29.0, 29.5, 30.0, 50.0, 100.0, 300.0, 500.0):
            radius = surface - depth_km * 1000.0
            if depth_km == 3.0:
                requests = [(11, "solid", "inner"), (12, "fluid", "outer")]
            elif depth_km == 29.5:
                requests = [(9, "solid", "inner"), (9, "solid", "outer")]
            else:
                candidates = [e for e in mesh if e["lo"] < radius < e["hi"]]
                if not candidates:
                    candidates = [e for e in mesh if abs(e["hi"] - radius) <= max(1e-7, 64 * math.ulp(radius))]
                if len(candidates) != 1:
                    raise ValueError(f"cannot resolve model layer at depth {depth_km} km")
                element = candidates[0]
                requests = [(element["layer"], "fluid" if element["fluid"] else "solid", "interior")]
            for layer, side, trace in requests:
                element = resolve(radius, side, trace, mesh)
                point_id = f"{angle_id}-d{depth_km:g}-{side}-{trace}"
                rows.append({"probe_id": point_id, "radius_m": f"{radius:.17g}",
                             "theta_rad": f"{theta:.17g}", "phi_rad": f"{phi:.17g}",
                             "layer": layer, "material_side": side, "radial_trace": trace,
                             "baseline_element_for_validation_only": element["element"],
                             "baseline_element_bounds_m_for_validation_only": [element["lo"], element["hi"]],
                             "angle_construction": construction})

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["probe_id", "radius_m", "theta_rad", "phi_rad",
                                                     "layer", "material_side", "radial_trace"],
                                lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in writer.fieldnames})
    receipt = {
        "schema": 1,
        "purpose": "shared fixed physical probe grid for Study 1; baseline element ids are validation metadata only",
        "probe_csv": args.csv.name,
        "probe_csv_sha256": sha256(args.csv),
        "probe_count": len(rows),
        "unique_angle_count": len(angles),
        "source": {"latitude_deg": source["latitude"], "longitude_deg": source["longitude"],
                   "geocentric_theta_rad": source_theta, "geocentric_phi_rad": source_phi,
                   "depth_km": source["depth_km"], "radius_m": source_radius},
        "receiver": {"geocentric_theta_rad": receiver_theta, "geocentric_phi_rad": receiver_phi,
                     "depth_km": case["receiver"]["depth_km"], "radius_m": receiver_radius,
                     "solid_side_functional_preserved": True},
        "angular_grid": {"path_along_deg": alongs, "cross_at_midpoint_deg": cross_offsets,
                         "path_distance_deg_recomputed": distance_deg,
                         "geometry_formula": "pinned Sphere::Colatitude plus SampledKernels.cpp Unit/PathFrame::At",
                         "receiver_endpoint_source": "hash-bound frozen Phase D along-angles.csv"},
        "depth_grid_km": [3, 10, 24, 29, 29.5, 30, 50, 100, 300, 500],
        "special_traces": {"receiver_radius": "layer 11 solid inner and layer 12 fluid outer",
                            "source_radius": "layer 9 solid inner and outer SEM traces"},
        "case_sha256": sha256(args.case),
        "mesh_receipt_sha256": sha256(args.mesh),
        "path_reference_sha256": sha256(args.path_reference),
        "generator_sha256": sha256(Path(__file__).resolve()),
        "mesh_receipt_provenance_elements": [
            {"probe_id": row["probe_id"], "element": row["baseline_element_for_validation_only"],
             "bounds_m": row["baseline_element_bounds_m_for_validation_only"]} for row in rows],
    }
    args.receipt.parent.mkdir(parents=True, exist_ok=True)
    args.receipt.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
