#!/usr/bin/env python3
"""Evaluate and replay matched DSpecM1D/YSpecpp source–receiver Earth slices."""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Circle
import numpy as np

MODULI = ("A", "C", "F", "L", "N")
RADIUS_M = 6_371_000.0
GRID_ANGLES = 181
GRID_RADII = 79
R_MIN_M = 100_000.0
SAMPLE_HEADER = ("objective", "modulus", "radius_m", "theta_rad", "phi_rad", "element",
                 "fluid", "kernel_objective_per_Pa_m3", "deltaX_Pa", "dr_m", "dOmega")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f"cannot read JSON {path}: {exc}") from exc


def validate_source_run(meta: dict, meta_path: Path) -> tuple[Path, int]:
    expected = {
        "objective": "J_Z=Re(sum(conj(q)*P u))",
        "nt": 128,
        "dt_s": 20,
        "frequency_bins_zero_based": [1, 2, 3],
        "gravity": "full",
        "attenuation": False,
    }
    for name, value in expected.items():
        if meta.get(name) != value:
            raise SystemExit(f"source run {name}={meta.get(name)!r}; expected {value!r}")
    lmax = int(meta.get("lmax", -1))
    if lmax not in (3, 12):
        raise SystemExit("source run must use the supported lmax 3 or 12 setup")
    if not math.isclose(float(meta.get("radius_m", 0)), RADIUS_M, rel_tol=0, abs_tol=1e-6):
        raise SystemExit("source run does not use the supported 6,371 km Earth radius")
    if meta.get("time_window", {}).get("type") != "Hann" or \
            meta.get("time_window", {}).get("sample_start_inclusive") != 0 or \
            meta.get("time_window", {}).get("sample_end_inclusive") != 90:
        raise SystemExit("source run has a different time window")
    if meta.get("q_nonzero") != {"Z,1": [0.7, 0.2], "Z,2": [-0.25, 0.45]}:
        raise SystemExit("source run has different J_Z objective weights")
    if not math.isclose(float(meta.get("epsilon_s_inv", 0)), 1.0 / 360.0,
                        rel_tol=0, abs_tol=1e-14):
        raise SystemExit("source run has a different damping rate")
    source = meta.get("source", {})
    receiver = meta.get("receiver", {})
    source_values = (source.get("latitude_deg"), source.get("longitude_deg"),
                     source.get("depth_km"), source.get("radius_m"))
    receiver_values = (receiver.get("latitude_deg"), receiver.get("longitude_deg"),
                       receiver.get("depth_km"))
    if not np.allclose(np.asarray(source_values, dtype=float),
                       (10.0, 20.0, 50.0, 6_321_000.0), rtol=0, atol=1e-8) or \
            not np.allclose(np.asarray(receiver_values, dtype=float),
                            (45.0, 90.0, 3.0), rtol=0, atol=1e-8):
        raise SystemExit("source run locations/depths differ from the timing adapters")
    if not np.allclose(np.asarray(meta.get("moment_tensor_Nm", []), dtype=float),
                       (1e19, 2e18, -3e18, -5e18, 4e18, -5e18), rtol=1e-12, atol=0):
        raise SystemExit("source run moment tensor differs from the timing adapters")
    if not math.isclose(float(meta.get("physical_G_SI", 0)), 6.6723e-11,
                        rel_tol=1e-12, abs_tol=0):
        raise SystemExit("source run physical G differs from the timing adapters")
    provenance = meta.get("provenance", {})
    model = Path(provenance.get("model", ""))
    expected_hash = provenance.get("model_sha256")
    if not model.is_file() or not expected_hash or sha256(model) != expected_hash:
        raise SystemExit("source-run model is missing or its recorded SHA-256 does not match")
    return model.resolve(), lmax


def unit_vector(latitude_deg: float, longitude_deg: float) -> np.ndarray:
    lat, lon = np.deg2rad([latitude_deg, longitude_deg])
    return np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])


def make_plane(meta: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Return radii, angles, geocentric coordinates, and source-to-receiver arc angle."""
    src, rec = meta["source"], meta["receiver"]
    source = unit_vector(float(src["latitude_deg"]), float(src["longitude_deg"]))
    receiver = unit_vector(float(rec["latitude_deg"]), float(rec["longitude_deg"]))
    cosine = float(np.clip(np.dot(source, receiver), -1.0, 1.0))
    delta = math.acos(cosine)
    tangent = receiver - cosine * source
    tangent_norm = np.linalg.norm(tangent)
    if tangent_norm < 1e-12:
        raise SystemExit("source and receiver do not define a unique great-circle plane")
    tangent /= tangent_norm
    angles = np.unique(np.r_[np.linspace(0.0, 2.0 * np.pi, GRID_ANGLES, endpoint=False),
                             np.pi, delta])
    # Preserve the exact source and receiver angular coordinates after sorting.
    angles[np.argmin(np.abs(angles))] = 0.0
    angles = np.sort(angles)
    source_radius = float(src["radius_m"])
    receiver_radius = RADIUS_M - float(rec["depth_km"]) * 1000.0
    radii = np.unique(np.r_[np.linspace(R_MIN_M, RADIUS_M, GRID_RADII),
                            source_radius, receiver_radius, RADIUS_M])
    xyz = np.cos(angles)[:, None] * source + np.sin(angles)[:, None] * tangent
    latitude = np.arcsin(np.clip(xyz[:, 2], -1.0, 1.0))
    longitude = np.arctan2(xyz[:, 1], xyz[:, 0]) % (2.0 * np.pi)
    theta = np.pi / 2.0 - latitude
    return radii, angles, np.column_stack((theta, longitude)), delta


def write_points(path: Path, radii: np.ndarray, angles: np.ndarray, coordinates: np.ndarray) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(SAMPLE_HEADER)
        for radius in radii:
            for theta, phi in coordinates:
                # Both timing adapters require 11 tokens but read geometry only. These
                # NA placeholders are not quadrature or perturbation weights.
                writer.writerow(("J_Z", "A", f"{radius:.17g}", f"{theta:.17g}",
                                 f"{phi:.17g}", -1, "NA", "NA", "NA", "NA", "NA"))


def run_adapter(executable: Path, model: Path, points: Path, out: Path,
                lmax: int, yspec: bool, angle_count: int) -> dict:
    out.mkdir()
    command = [str(executable), str(model), str(points), str(out), str(lmax), "1", "1"]
    if yspec:
        command.append(str(angle_count))
    result = subprocess.run(command, text=True, capture_output=True, timeout=150)
    if result.returncode:
        raise SystemExit(f"{executable.name} failed ({result.returncode}):\n{result.stderr[-4000:]}")
    try:
        metrics = json.loads(result.stdout.strip().splitlines()[-1])
    except (IndexError, json.JSONDecodeError) as exc:
        raise SystemExit(f"{executable.name} did not return its run metadata: {exc}") from exc
    (out / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    return metrics


def validate_adapter_metrics(metrics: dict, solver: str, lmax: int, sample_count: int) -> None:
    if metrics.get("solver", "").lower() != solver.lower() or int(metrics.get("lmax", -1)) != lmax:
        raise SystemExit(f"{solver} adapter resolved a different solver or lmax")
    count = int(metrics.get("sample_count", metrics.get("point_count", -1)))
    if count != sample_count:
        raise SystemExit(f"{solver} adapter returned {count} points; expected {sample_count}")
    config = metrics.get("config", {})
    grid = config.get("frequency_grid", {})
    if grid.get("nt") != 128 or grid.get("solver_bins_zero_based") != [1, 2, 3] or \
            not math.isclose(float(grid.get("dt_s", 0)), 20.0, rel_tol=0, abs_tol=1e-12) or \
            not math.isclose(float(grid.get("epsilon_s_inv", 0)), 1 / 360, rel_tol=0, abs_tol=1e-14):
        raise SystemExit(f"{solver} adapter resolved a different frequency grid")
    physics = config.get("physics", {})
    if physics.get("gravity") != "full" or physics.get("attenuation") is not False or \
            not math.isclose(float(physics.get("physical_G", 0)), 6.6723e-11, rel_tol=1e-12):
        raise SystemExit(f"{solver} adapter resolved different physics")
    src, rec = config.get("source", {}), config.get("receiver", {})
    source_values = (src.get("geocentric_latitude_deg"), src.get("longitude_deg"), src.get("depth_m"))
    receiver_values = (rec.get("geocentric_latitude_deg"), rec.get("longitude_deg"), rec.get("depth_m"))
    if not np.allclose(np.asarray(source_values, dtype=float), (10.0, 20.0, 50000.0), rtol=0, atol=1e-7) or \
            not np.allclose(np.asarray(receiver_values, dtype=float), (45.0, 90.0, 3000.0), rtol=0, atol=1e-7):
        raise SystemExit(f"{solver} adapter resolved different source or receiver coordinates")
    if not np.allclose(np.asarray(src.get("moment_tensor_Nm", []), dtype=float),
                       (1e19, 2e18, -3e18, -5e18, 4e18, -5e18), rtol=1e-12, atol=0):
        raise SystemExit(f"{solver} adapter resolved a different moment tensor")
    window = config.get("window", {})
    if window.get("type") != "Hann" or window.get("first_sample") != 0 or \
            window.get("last_sample") != 90:
        raise SystemExit(f"{solver} adapter resolved a different time window")
    objectives = config.get("objective", {}).get("names", [])
    if len(objectives) != 1 or objectives[0].get("name") != "J_Z" or \
            objectives[0].get("component") != "Z" or \
            objectives[0].get("q_nonzero") != {"bin_1": [0.7, 0.2], "bin_2": [-0.25, 0.45]}:
        raise SystemExit(f"{solver} adapter resolved a different J_Z objective")


def read_kernel(path: Path, expected_radii: np.ndarray, expected_angles: np.ndarray,
                coordinates: np.ndarray, angle_count: int) -> dict[str, np.ndarray]:
    fields = {modulus: np.full((len(expected_radii), angle_count), np.nan) for modulus in MODULI}
    with path.open(newline="", encoding="utf-8") as stream:
        rows = csv.DictReader(stream)
        for row in rows:
            if row["objective"] != "J_Z":
                continue
            modulus = row["modulus"]
            if modulus not in fields:
                continue
            radius = float(row["radius_m"])
            theta, phi = float(row["theta_rad"]), float(row["phi_rad"])
            ir = int(np.argmin(np.abs(expected_radii - radius)))
            ia = int(np.argmin(np.abs(coordinates[:, 0] - theta) +
                                np.abs(np.angle(np.exp(1j * (coordinates[:, 1] - phi))))))
            if not math.isclose(radius, expected_radii[ir], rel_tol=0, abs_tol=1e-6) or \
                    not np.allclose((theta, phi), coordinates[ia], rtol=0, atol=1e-12):
                raise SystemExit(f"kernel coordinates in {path} do not match points.csv")
            value = float(row["kernel_objective_per_Pa_m3"] if "kernel_objective_per_Pa_m3" in row
                          else row["kernel_per_Pa_m3"])
            if not math.isfinite(value) or math.isfinite(fields[modulus][ir, ia]):
                raise SystemExit(f"nonfinite or duplicate {modulus} result in {path}")
            fields[modulus][ir, ia] = value
    for modulus, values in fields.items():
        if not np.isfinite(values).all():
            raise SystemExit(f"{path} has incomplete {modulus} Earth-slice values")
    return fields


def geometry_edges(radii: np.ndarray, angles: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    redges = np.empty(len(radii) + 1)
    redges[0], redges[-1] = R_MIN_M, RADIUS_M
    redges[1:-1] = (radii[:-1] + radii[1:]) / 2.0
    # The first sample ring begins at 100 km; leave the unsampled centre white.
    aedges = np.empty(len(angles) + 1)
    aedges[0] = (angles[-1] - 2 * np.pi + angles[0]) / 2.0
    aedges[1:-1] = (angles[:-1] + angles[1:]) / 2.0
    aedges[-1] = aedges[0] + 2 * np.pi
    return redges, aedges


def point_xy(radius: float, angle: float) -> tuple[float, float]:
    return radius * math.sin(angle), radius * math.cos(angle)


def comparison_diagnostics(radii: np.ndarray, angles: np.ndarray, dspec: np.ndarray,
                            yspec: np.ndarray, source_radius: float,
                            receiver_radius: float, receiver_angle: float) -> dict:
    difference = dspec - yspec
    peak = float(max(np.max(np.abs(dspec)), np.max(np.abs(yspec))))
    physical_distance = np.sqrt(np.maximum(
        0.0, radii[:, None] ** 2 + source_radius ** 2 -
        2.0 * radii[:, None] * source_radius * np.cos(angles)[None, :]))
    receiver_distance = np.sqrt(np.maximum(
        0.0, radii[:, None] ** 2 + receiver_radius ** 2 -
        2.0 * radii[:, None] * receiver_radius *
        np.cos(angles[None, :] - receiver_angle)))
    masks = {
        "full_plane": np.ones(difference.shape, dtype=bool),
        "near_receiver_distance_lt_500km": receiver_distance < 500_000.0,
        "near_center_radius_lt_1000km": np.broadcast_to(
            (radii < 1_000_000.0)[:, None], difference.shape),
        "near_source_distance_lt_500km": physical_distance < 500_000.0,
        "near_source_excluding_exact_source_radius":
            (physical_distance < 500_000.0) & (radii[:, None] != source_radius),
        "exact_source_radius_all_angles": np.broadcast_to(
            (radii == source_radius)[:, None], difference.shape),
        "exact_source_point": np.broadcast_to(
            ((radii == source_radius)[:, None] & (np.abs(angles[None, :]) < 1e-14)),
            difference.shape),
        "validated_shell_4200_to_6200km": np.broadcast_to(
            ((radii >= 4_200_000.0) & (radii <= 6_200_000.0))[:, None], difference.shape),
        "below_validated_shell_lt_4200km": np.broadcast_to(
            (radii < 4_200_000.0)[:, None], difference.shape),
        "above_validated_shell_gt_6200km": np.broadcast_to(
            (radii > 6_200_000.0)[:, None], difference.shape),
        "above_validated_shell_excluding_exact_source_radius": np.broadcast_to(
            ((radii > 6_200_000.0) & (radii != source_radius))[:, None], difference.shape),
    }
    summary = {"global_solver_peak_abs": peak,
               "pointwise_relative_floor": 0.01 * peak,
               "regions": {}, "per_radius": []}
    for name, mask in masks.items():
        d, y, diff = dspec[mask], yspec[mask], difference[mask]
        if not diff.size:
            continue
        local_peak = float(max(np.max(np.abs(d)), np.max(np.abs(y))))
        relative = np.abs(diff) / np.maximum(np.maximum(np.abs(d), np.abs(y)), 0.01 * peak)
        flat = int(np.argmax(np.abs(diff)))
        indices = np.argwhere(mask)
        ir, ia = map(int, indices[flat])
        summary["regions"][name] = {
            "sample_count": int(diff.size),
            "dspec_peak_abs": float(np.max(np.abs(d))),
            "yspec_peak_abs": float(np.max(np.abs(y))),
            "max_abs_difference": float(np.max(np.abs(diff))),
            "max_difference_over_region_peak": float(np.max(np.abs(diff)) / local_peak)
                if local_peak else 0.0,
            "relative_l2": float(np.linalg.norm(diff) /
                                  max(np.linalg.norm(d), np.linalg.norm(y), 1e-300)),
            "p99_abs_difference": float(np.quantile(np.abs(diff), 0.99)),
            "max_pointwise_relative_above_floor": float(np.max(relative)),
            "largest_difference_location": {
                "radius_m": float(radii[ir]), "angle_deg": float(np.degrees(angles[ia])),
                "source_distance_m": float(physical_distance[ir, ia]),
            },
        }
    for ir, radius in enumerate(radii):
        row_diff = difference[ir]
        angle_index = int(np.argmax(np.abs(row_diff)))
        row_d, row_y = dspec[ir], yspec[ir]
        summary["per_radius"].append({
            "radius_m": float(radius),
            "relative_l2": float(np.linalg.norm(row_diff) /
                                  max(np.linalg.norm(row_d), np.linalg.norm(row_y), 1e-300)),
            "max_abs_difference": float(np.max(np.abs(row_diff))),
            "max_difference_angle_deg": float(np.degrees(angles[angle_index])),
        })
    return summary


def save_disk(out: Path, modulus: str, radii: np.ndarray, angles: np.ndarray,
              dspec: np.ndarray, yspec: np.ndarray, meta: dict) -> dict:
    solver_limit = float(max(np.nanmax(np.abs(dspec)), np.nanmax(np.abs(yspec)))) or 1.0
    difference = dspec - yspec
    diff_limit = float(np.nanmax(np.abs(difference))) or 1.0
    norm = Normalize(-solver_limit, solver_limit)
    diff_norm = Normalize(-diff_limit, diff_limit)
    redges, aedges = geometry_edges(radii, angles)
    ae, re = np.meshgrid(aedges, redges)
    x, z = re * np.sin(ae) / 1000.0, re * np.cos(ae) / 1000.0
    dplot, yplot = dspec, yspec
    delta_plot = dplot - yplot
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 5.6), constrained_layout=True)
    images = []
    for ax, values, title, cmap, scale in (
        (axes[0], dplot, "DSpecM1D", "RdBu_r", norm),
        (axes[1], yplot, "YSpecpp", "RdBu_r", norm),
        (axes[2], delta_plot, "DSpecM1D − YSpecpp", "RdBu_r", diff_norm),
    ):
        image = ax.pcolormesh(x, z, values, shading="flat", cmap=cmap, norm=scale,
                              rasterized=True)
        images.append(image)
        ax.add_patch(Circle((0, 0), RADIUS_M / 1000, fill=False, color="0.15", lw=0.8))
        ax.add_patch(Circle((0, 0), R_MIN_M / 1000, color="white", ec="none", zorder=3))
        ax.set_aspect("equal")
        ax.set_xlim(-RADIUS_M / 1000 * 1.06, RADIUS_M / 1000 * 1.06)
        ax.set_ylim(-RADIUS_M / 1000 * 1.06, RADIUS_M / 1000 * 1.06)
        ax.set_xlabel("Distance in source–receiver plane (km)")
        ax.set_title(title)
        ax.grid(color="0.4", alpha=0.15, lw=0.4)
    axes[0].set_ylabel("Distance in source–receiver plane (km)")
    src, rec = meta["source"], meta["receiver"]
    delta = float(meta["source_receiver_arc_rad"])
    source_xy = point_xy(float(src["radius_m"]), 0.0)
    receiver_xy = point_xy(RADIUS_M - float(rec["depth_km"]) * 1000.0, delta)
    for ax in axes:
        for label, xy, marker, color in (("S", source_xy, "*", "#f28e2b"),
                                         ("R", receiver_xy, "^", "#16a085")):
            ax.scatter(xy[0] / 1000, xy[1] / 1000, s=74, marker=marker, color=color,
                       edgecolor="black", linewidth=0.5, zorder=5)
            offset = (7, -16) if label == "S" else (-14, 9)
            ax.annotate(label, (xy[0] / 1000, xy[1] / 1000), xytext=offset,
                        textcoords="offset points", fontsize=8, fontweight="bold",
                        color="black")
    colorbar = fig.colorbar(images[0], ax=axes[:2].tolist(), shrink=0.84, pad=0.02)
    colorbar.set_label("J_Z kernel (objective Pa⁻¹ m⁻³)")
    diffbar = fig.colorbar(images[2], ax=axes[2], shrink=0.84, pad=0.02)
    diffbar.set_label("Raw DSpecM1D − YSpecpp (objective Pa⁻¹ m⁻³)")
    angle_deg = math.degrees(delta)
    fig.suptitle(f"{modulus} · homogeneous sphere · lmax {meta['lmax']} · spectral J_Z · Δ={angle_deg:.1f}°",
                 fontsize=13)
    src_radius = float(src["radius_m"]) / 1000.0
    rec_radius = (RADIUS_M - float(rec["depth_km"]) * 1000.0) / 1000.0
    fig.text(0.5, 0.015,
             f"★ Source r={src_radius:.0f} km   ▲ Receiver r={rec_radius:.0f} km   ·   "
             "blank centre <100 km   ·   matched-volume checks: 4,200–6,200 km; outside unvalidated",
             ha="center", fontsize=8)
    for ext in ("png", "svg", "pdf"):
        fig.savefig(out / f"earth_slice_{modulus}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return {"solver_limit": solver_limit, "difference_limit": diff_limit}


def save_oblique(out: Path, modulus: str, radii: np.ndarray, angles: np.ndarray,
                 dspec: np.ndarray, yspec: np.ndarray, source: np.ndarray,
                 tangent: np.ndarray, meta: dict) -> dict:
    xyz = np.cos(angles)[:, None] * source + np.sin(angles)[:, None] * tangent
    points = radii[:, None, None] * xyz[None, :, :]
    x, y, z = [points[:, :, i] / 1000.0 for i in range(3)]
    dlimit = float(max(np.nanmax(np.abs(dspec)), np.nanmax(np.abs(yspec)))) or 1.0
    difference = dspec - yspec
    diff_limit = float(np.nanmax(np.abs(difference))) or 1.0
    cmap = plt.get_cmap("RdBu_r")
    panels = ((dspec, "DSpecM1D", Normalize(-dlimit, dlimit),
               "J_Z kernel (objective Pa⁻¹ m⁻³)"),
              (yspec, "YSpecpp", Normalize(-dlimit, dlimit),
               "J_Z kernel (objective Pa⁻¹ m⁻³)"),
              (difference, "DSpecM1D − YSpecpp", Normalize(-diff_limit, diff_limit),
               "Raw J_Z difference (objective Pa⁻¹ m⁻³)"))
    earth = RADIUS_M / 1000.0
    plane_normal = np.cross(source, tangent)
    plane_normal /= np.linalg.norm(plane_normal)
    view_direction = math.cos(math.radians(32.0)) * plane_normal + \
        math.sin(math.radians(32.0)) * source
    view_elevation = math.degrees(math.asin(float(np.clip(view_direction[2], -1.0, 1.0))))
    view_azimuth = math.degrees(math.atan2(view_direction[1], view_direction[0]))
    edge = earth * xyz
    sphere_lat = np.linspace(-np.pi / 2, np.pi / 2, 33)
    sphere_lon = np.linspace(0, 2 * np.pi, 65)
    slon, slat = np.meshgrid(sphere_lon, sphere_lat)
    sx = earth * np.cos(slat) * np.cos(slon)
    sy = earth * np.cos(slat) * np.sin(slon)
    sz = earth * np.sin(slat)
    fig = plt.figure(figsize=(15.0, 6.0), constrained_layout=True)
    src, rec = meta["source"], meta["receiver"]
    src_xyz = unit_vector(float(src["latitude_deg"]), float(src["longitude_deg"])) * float(src["radius_m"])
    rec_xyz = unit_vector(float(rec["latitude_deg"]), float(rec["longitude_deg"])) * \
        (RADIUS_M - float(rec["depth_km"]) * 1000.0)
    for index, (values, title, norm, colorbar_label) in enumerate(panels, start=1):
        ax = fig.add_subplot(1, 3, index, projection="3d")
        ax.plot_surface(sx, sy, sz, color="0.86", alpha=0.08, linewidth=0,
                        antialiased=False, shade=False, rstride=2, cstride=2)
        # Geocentric latitude and longitude graticule gives the cut a geographic frame.
        lon_line = np.linspace(0, 2 * np.pi, 181)
        for lat_deg in (-60, -30, 0, 30, 60):
            lat = math.radians(lat_deg)
            gx = earth * math.cos(lat) * np.cos(lon_line)
            gy = earth * math.cos(lat) * np.sin(lon_line)
            gz = np.full_like(lon_line, earth * math.sin(lat))
            ax.plot(gx, gy, gz, color="0.45", alpha=0.2, lw=0.45)
        lat_line = np.linspace(-np.pi / 2, np.pi / 2, 91)
        for lon_deg in range(0, 360, 30):
            lon = math.radians(lon_deg)
            gx = earth * np.cos(lat_line) * math.cos(lon)
            gy = earth * np.cos(lat_line) * math.sin(lon)
            gz = earth * np.sin(lat_line)
            ax.plot(gx, gy, gz, color="0.45", alpha=0.16, lw=0.4)
        closed_values = np.c_[values, values[:, :1]]
        colors = cmap(norm(closed_values))
        xc, yc, zc = [np.c_[axis, axis[:, :1]] for axis in (x, y, z)]
        surface = ax.plot_surface(xc, yc, zc, facecolors=colors, rstride=1, cstride=1,
                                  linewidth=0, antialiased=False, shade=False)
        surface.set_rasterized(True)
        for label, position, marker, color in (("S", src_xyz, "*", "#f28e2b"),
                                               ("R", rec_xyz, "^", "#16a085")):
            pos = position / 1000.0
            ax.scatter(*pos, marker=marker, s=75, color=color, edgecolor="black", depthshade=False)
            ax.text(*(pos * 1.018), label, fontsize=8, fontweight="bold")
        ax.plot(edge[:, 0], edge[:, 1], edge[:, 2], color="0.15", lw=0.6)
        ax.set_box_aspect((1, 1, 1))
        ax.set_xlim(-earth, earth); ax.set_ylim(-earth, earth); ax.set_zlim(-earth, earth)
        ax.view_init(elev=view_elevation, azim=view_azimuth)
        ax.set_axis_off()
        ax.set_title(title, pad=8)
        mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        fig.colorbar(mappable, ax=ax, shrink=0.62, pad=0.03, label=colorbar_label)
    src_radius = float(meta["source"]["radius_m"]) / 1000.0
    rec_radius = (RADIUS_M - float(meta["receiver"]["depth_km"]) * 1000.0) / 1000.0
    fig.suptitle(f"{modulus} · source–receiver plane · lmax {meta['lmax']} · spectral J_Z",
                 fontsize=13)
    fig.text(0.5, 0.015,
             f"Geocentric latitude/longitude graticule, 30° spacing   ·   "
             f"★ Source r={src_radius:.0f} km   ▲ Receiver r={rec_radius:.0f} km",
             ha="center", fontsize=8)
    for ext in ("png", "svg", "pdf"):
        fig.savefig(out / f"earth_slice_oblique_{modulus}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return {"solver_limit": dlimit, "difference_limit": diff_limit}


def create_output(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)
    if any(path.iterdir()):
        raise SystemExit(f"output must be empty: {path}")


def write_report(out: Path, run: dict) -> None:
    meta = read_json(out / "source_run.json")
    source, receiver = run["source"], run["receiver"]
    delta_deg = math.degrees(float(run["source_receiver_arc_rad"]))
    source_r_km = float(source["radius_m"]) / 1000.0
    receiver_r_km = (RADIUS_M - float(receiver["depth_km"]) * 1000.0) / 1000.0
    def fmt(value: float) -> str:
        return f"{value:.3e}"
    rows = []
    for modulus in MODULI:
        regions = run["comparison_diagnostics"][modulus]["regions"]
        full = regions["full_plane"]
        shell = regions["validated_shell_4200_to_6200km"]
        near = regions["near_receiver_distance_lt_500km"]
        outer = regions["above_validated_shell_excluding_exact_source_radius"]
        rows.append("<tr><th>" + modulus + "</th><td>" + fmt(full["relative_l2"]) +
                    "</td><td>" + fmt(shell["relative_l2"]) + "</td><td>" +
                    fmt(near["relative_l2"]) + "</td><td>" + fmt(near["max_abs_difference"]) +
                    "</td><td>" + fmt(outer["relative_l2"]) + "</td></tr>")
    gallery = []
    for kind, label in (("earth_slice", "Circular section"),
                        ("earth_slice_oblique", "Oblique plane")):
        for modulus in MODULI:
            stem = f"{kind}_{modulus}"
            gallery.append(f'<figure><a href="{stem}.png"><img src="{stem}.png" alt="{label}, {modulus}"></a>'
                           f'<figcaption>{label} · {modulus} · <a href="{stem}.pdf">PDF</a> · '
                           f'<a href="{stem}.svg">SVG</a></figcaption></figure>')
    source_desc = (f"{float(source['latitude_deg']):g}°, {float(source['longitude_deg']):g}°; "
                   f"depth {float(source['depth_km']):g} km; r={source_r_km:.0f} km")
    receiver_desc = (f"{float(receiver['latitude_deg']):g}°, {float(receiver['longitude_deg']):g}°; "
                     f"depth {float(receiver['depth_km']):g} km; r={receiver_r_km:.0f} km")
    html_lines = [
        "<!doctype html>", '<html lang="en"><meta charset="utf-8">',
        "<title>Matched source–receiver Earth slices</title>",
        "<style>body{font:16px system-ui,sans-serif;max-width:1100px;margin:2em auto;padding:0 1em;color:#222}"
        "table{border-collapse:collapse}td,th{padding:.45em .7em;border:1px solid #bbb;text-align:right}"
        "figure{margin:0}figcaption{font-size:.9em} .gallery{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:1em}"
        "img{width:100%;height:auto;border:1px solid #ddd} code{overflow-wrap:anywhere}</style>",
        "<h1>Matched source–receiver Earth slices</h1>",
        f"<p>{html.escape(run['objective'])}; {html.escape(meta['kernel_units'])}. "
        f"Homogeneous spherical fixture, lmax {int(run['lmax'])}; spectral point kernels.</p>",
        "<h2>Experiment</h2><dl>",
        f"<dt>Source</dt><dd>{html.escape(source_desc)}</dd>",
        f"<dt>Receiver</dt><dd>{html.escape(receiver_desc)}</dd>",
        f"<dt>Great-circle separation</dt><dd>{delta_deg:.3f}°</dd>",
        f"<dt>Frequency bins (zero based)</dt><dd>{html.escape(str(meta['frequency_bins_zero_based']))}; "
        f"dt={float(meta['dt_s']):g} s, nt={int(meta['nt'])}</dd>",
        f"<dt>Window</dt><dd>{html.escape(meta['time_window']['type'])}, samples "
        f"{int(meta['time_window']['sample_start_inclusive'])}–{int(meta['time_window']['sample_end_inclusive'])}</dd>",
        f"<dt>Model</dt><dd>{html.escape(Path(run['model']).name)}; SHA-256 "
        f"<code>{html.escape(run['model_sha256'])}</code></dd>",
        f"<dt>Sample grid</dt><dd>{int(run['radial_count'])} radii × {int(run['angular_count'])} angles; "
        f"{int(run['sample_count']):,} shared points; r=100–6,371 km</dd>", "</dl>",
        "<h2>DSpecM1D versus YSpecpp</h2>",
        "<p>Relative L2 is ||D−Y||₂ / max(||D||₂, ||Y||₂). Max differences are raw and unmasked. "
        "The previously checked volume interval is 4,200–6,200 km; other queried values remain unvalidated. "
        "Source-shell contact is separate from these volume samples.</p>",
        "<table><thead><tr><th>Modulus</th><th>Full plane relative L2</th>"
        "<th>4,200–6,200 km relative L2</th><th>Within 500 km of receiver relative L2</th>"
        "<th>Near-receiver max |D−Y|</th><th>r&gt;6,200 km excluding exact source radius relative L2</th>"
        "</tr></thead><tbody>", *rows, "</tbody></table>",
        '<h2>Figures</h2><div class=\"gallery\">', *gallery, '</div>',
        '<h2>Data and provenance</h2><p><a href="run.json">Run diagnostics</a> · '
        '<a href="source_run.json">Source run metadata</a> · <a href="points.csv">Sample coordinates</a> · '
        '<a href="dspec/kernels.csv">DSpecM1D kernels</a> · '
        '<a href="yspec/kernels.csv">YSpecpp kernels</a></p>',
        '<h2>Figure context</h2><p>Source-centered Earth sections and source/receiver symbols follow '
        '<a href="https://academic.oup.com/gji/article/174/1/265/577094">Liu &amp; Tromp (2008), Fig. 2</a>; '
        'the volume-kernel framing is consistent with <a href="https://doi.org/10.1093/gji/ggae136">'
        'Adourian et al. (2024), Fig. 3</a>. The present panels are direct low-frequency J<sub>Z</sub> '
        'point-kernel evaluations, not ray paths or banana–doughnut sensitivity.</p>',
        "</html>",
    ]
    (out / "report.html").write_text("\n".join(html_lines) + "\n", encoding="utf-8")


def render(out: Path, run: dict) -> None:
    radii = np.asarray(run["radii_m"], dtype=float)
    angles = np.asarray(run["angles_rad"], dtype=float)
    coordinates = np.asarray(run["coordinates_theta_phi_rad"], dtype=float)
    dfields = read_kernel(out / "dspec" / "kernels.csv", radii, angles, coordinates, len(angles))
    yfields = read_kernel(out / "yspec" / "kernels.csv", radii, angles, coordinates, len(angles))
    source = unit_vector(float(run["source"]["latitude_deg"]), float(run["source"]["longitude_deg"]))
    receiver = unit_vector(float(run["receiver"]["latitude_deg"]), float(run["receiver"]["longitude_deg"]))
    delta = float(run["source_receiver_arc_rad"])
    tangent = (receiver - math.cos(delta) * source) / math.sin(delta)
    run["comparison_scope"] = (
        "The established matched-volume comparison applies only to 4,200–6,200 km. "
        "All queried values are retained; differences outside that interval are shown raw and unvalidated. "
        "Plots contain volume point kernels only; source-shell contact is separate.")
    run["comparison_diagnostics"] = {}
    for modulus in MODULI:
        receiver_radius = RADIUS_M - float(run["receiver"]["depth_km"]) * 1000.0
        run["comparison_diagnostics"][modulus] = comparison_diagnostics(
            radii, angles, dfields[modulus], yfields[modulus],
            float(run["source"]["radius_m"]), receiver_radius,
            float(run["source_receiver_arc_rad"]))
        run["panel_scales"][modulus] = save_disk(out, modulus, radii, angles,
                                                  dfields[modulus], yfields[modulus], run)
        run.setdefault("oblique_panel_scales", {})[modulus] = save_oblique(
            out, modulus, radii, angles, dfields[modulus], yfields[modulus],
            source, tangent, run)
    (out / "run.json").write_text(json.dumps(run, indent=2) + "\n", encoding="utf-8")
    write_report(out, run)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--source-run", type=Path, help="DSpec tutorial run metadata for a fresh matched slice")
    mode.add_argument("--replay", type=Path, help="saved Earth-slice run directory; render without solvers")
    parser.add_argument("--dspec-exe", type=Path, help="DSpec timing adapter executable")
    parser.add_argument("--yspec-exe", type=Path, help="YSpec timing adapter executable")
    parser.add_argument("--output", type=Path, required=True, help="new or empty output directory")
    args = parser.parse_args()
    out = args.output.resolve()
    create_output(out)
    if args.replay:
        replay = args.replay.resolve()
        run = read_json(replay / "run.json")
        for name in ("points.csv", "source_run.json", "dspec", "yspec"):
            source = replay / name
            destination = out / name
            if source.is_dir():
                shutil.copytree(source, destination)
            else:
                shutil.copy2(source, destination)
        run["replay_of"] = str(replay)
        run["render_mode"] = "saved numeric inputs; no solver calls"
        render(out, run)
        return 0
    if args.dspec_exe is None or args.yspec_exe is None:
        parser.error("fresh slice generation requires --dspec-exe and --yspec-exe")
    source_path = args.source_run.resolve()
    source_meta = read_json(source_path)
    model, lmax = validate_source_run(source_meta, source_path)
    radii, angles, coordinates, delta = make_plane(source_meta)
    points = out / "points.csv"
    write_points(points, radii, angles, coordinates)
    dspec_exe, yspec_exe = args.dspec_exe.resolve(), args.yspec_exe.resolve()
    if not dspec_exe.is_file() or not yspec_exe.is_file():
        raise SystemExit("both timing adapter executables must exist")
    dmetrics = run_adapter(dspec_exe, model, points, out / "dspec", lmax, False, len(angles))
    ymetrics = run_adapter(yspec_exe, model, points, out / "yspec", lmax, True, len(angles))
    validate_adapter_metrics(dmetrics, "DSpecM1D", lmax, len(radii) * len(angles))
    validate_adapter_metrics(ymetrics, "YSpecpp", lmax, len(radii) * len(angles))
    run = {
        "schema": "matched-earth-slice-v1", "objective": source_meta["objective"],
        "lmax": lmax, "earth_radius_m": RADIUS_M, "minimum_sample_radius_m": R_MIN_M,
        "adapter_timeout_seconds": 150,
        "source": source_meta["source"], "receiver": source_meta["receiver"],
        "source_receiver_arc_rad": delta,
        "geometry": "full geocentric source–receiver great-circle plane; no ray path is inferred",
        "radii_m": radii.tolist(), "angles_rad": angles.tolist(),
        "coordinates_theta_phi_rad": coordinates.tolist(),
        "sample_count": int(len(radii) * len(angles)),
        "angular_count": int(len(angles)), "radial_count": int(len(radii)),
        "source_run": str(source_path), "source_run_sha256": sha256(source_path),
        "model": str(model), "model_sha256": sha256(model),
        "points_sha256": sha256(points),
        "kernels_sha256": {"DSpecM1D": sha256(out / "dspec" / "kernels.csv"),
                            "YSpecpp": sha256(out / "yspec" / "kernels.csv")},
        "adapter_metrics": {"DSpecM1D": dmetrics, "YSpecpp": ymetrics},
        "panel_scales": {}, "convergence_claim": False,
        "fixture_note": "Homogeneous spherical model; no internal interface contrast is represented.",
    }
    shutil.copy2(source_path, out / "source_run.json")
    (out / "run.json").write_text(json.dumps(run, indent=2) + "\n", encoding="utf-8")
    render(out, run)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError, IndexError, subprocess.TimeoutExpired) as exc:
        print(f"earth_slices.py: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc
