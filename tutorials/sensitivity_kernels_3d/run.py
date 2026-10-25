#!/usr/bin/env python3
"""Run an installed 3D consumer and build reproducible, replayable kernel plots."""
import argparse
import csv
import html
import hashlib
import json
import math
from pathlib import Path
import shlex
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize

plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
                     "xtick.labelsize": 9, "ytick.labelsize": 9})
MODULI = ("A", "C", "F", "L", "N")


def read(path):
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise SystemExit(f"empty sample file: {path}")
    return np.array(rows, dtype=object)


def read_yspec(path):
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    required = {"objective", "modulus", "radius_m", "theta_rad", "phi_rad", "kernel_per_Pa_m3"}
    if not rows or not required.issubset(rows[0]):
        raise SystemExit(f"YSpec profile file lacks columns: {sorted(required - set(rows[0] if rows else []))}")
    return [r for r in rows if r["objective"] == "J_Z"]


def validate_yspec_match(path, metadata_path, dspec_meta, dspec_rows, dspec_samples_path):
    """Reject comparisons unless the saved YSpec run is the same experiment."""
    if metadata_path is None:
        raise SystemExit("YSpec comparison requires --yspec-meta from the same timed run")
    try:
        ymeta = json.loads(metadata_path.read_text(encoding="utf-8"))
        cfg = ymeta["config"]
        frequency = cfg["frequency_grid"]
        physics = cfg["physics"]
        source = cfg["source"]
        receiver = cfg["receiver"]
        objective = cfg["objective"]
        window = cfg["window"]
    except (OSError, KeyError, json.JSONDecodeError, TypeError) as error:
        raise SystemExit(f"invalid YSpec matched-run metadata {metadata_path}: {error}")
    if ymeta.get("schema") != "matched-kernel-timing-run-v1" or ymeta.get("solver", "").lower() != "yspecpp":
        raise SystemExit("YSpec metadata must use matched-kernel-timing-run-v1 and identify YSpecpp")
    hashes = ymeta.get("provenance", {})
    if not hashes.get("model_sha256") or not hashes.get("kernels_sha256") or not hashes.get("samples_sha256"):
        raise SystemExit("YSpec metadata must include model, kernel, and sample SHA-256 hashes")
    dspec_samples_sha256 = hashlib.sha256(dspec_samples_path.read_bytes()).hexdigest()
    if hashes["samples_sha256"] != dspec_samples_sha256:
        raise SystemExit("YSpec sampled-point hash does not match the DSpec samples.csv being plotted")
    if hashes.get("kernels_sha256") != hashlib.sha256(path.read_bytes()).hexdigest():
        raise SystemExit("YSpec profiles hash does not match the supplied profile file")
    provenance = dspec_meta.get("provenance", {})
    if provenance.get("model_sha256") and hashes["model_sha256"] != provenance["model_sha256"]:
        raise SystemExit("DSpec and YSpec model hashes differ")
    expected = {
        "nt": dspec_meta["nt"],
        "dt_s": dspec_meta["dt_s"],
        "df_hz": dspec_meta["df_hz"],
        "epsilon_s_inv": dspec_meta["epsilon_s_inv"],
        "solver_bins_zero_based": [1, 2, 3],
        "degrees_inclusive": [0, dspec_meta["lmax"]],
    }
    for name, value in expected.items():
        actual = frequency.get(name)
        if isinstance(value, float):
            matches = actual is not None and math.isclose(float(actual), value, rel_tol=1e-10, abs_tol=1e-12)
        else:
            matches = actual == value
        if not matches:
            raise SystemExit(f"YSpec frequency_grid.{name}={actual!r} does not match DSpec {value!r}")
    if physics.get("gravity") != "full" or physics.get("attenuation") is not False:
        raise SystemExit("YSpec comparison requires full gravity and attenuation disabled")
    if not math.isclose(float(physics.get("physical_G")), float(dspec_meta["physical_G_SI"]), rel_tol=1e-10):
        raise SystemExit("YSpec physical G does not match DSpec metadata")
    ds = dspec_meta["source"]
    dr = dspec_meta["receiver"]
    expected_source = (float(ds["latitude_deg"]), float(ds["longitude_deg"]), float(ds["depth_km"]) * 1000)
    actual_source = (float(source["geocentric_latitude_deg"]), float(source["longitude_deg"]), float(source["depth_m"]))
    expected_receiver = (float(dr["latitude_deg"]), float(dr["longitude_deg"]), float(dr["depth_km"]) * 1000)
    actual_receiver = (float(receiver["geocentric_latitude_deg"]), float(receiver["longitude_deg"]), float(receiver["depth_m"]))
    if not np.allclose(expected_source, actual_source, rtol=0, atol=1e-6):
        raise SystemExit("YSpec source location/depth does not match DSpec metadata")
    if not np.allclose(expected_receiver, actual_receiver, rtol=0, atol=1e-6):
        raise SystemExit("YSpec receiver location/depth does not match DSpec metadata")
    if not np.allclose(np.asarray(source["moment_tensor_Nm"], dtype=float),
                       np.asarray(dspec_meta["moment_tensor_Nm"], dtype=float), rtol=1e-12, atol=0):
        raise SystemExit("YSpec moment tensor does not match DSpec metadata")
    objectives = objective.get("names", [])
    matched_objective = next((item for item in objectives if item.get("name") == "J_Z"), None)
    expected_q = {"bin_1": [0.7, 0.2], "bin_2": [-0.25, 0.45]}
    if matched_objective is None or matched_objective.get("q_nonzero") != expected_q:
        raise SystemExit("YSpec objective/q weights do not match the DSpec J_Z objective")
    expected_window = {"type": "Hann", "first_sample": 0, "last_sample": 90, "nt": 128}
    if any(window.get(k) != v for k, v in expected_window.items()):
        raise SystemExit("YSpec time window does not match the DSpec Hann samples 0–90 of 128")
    dspec_a_rows = [row for row in dspec_rows if row["modulus"] == "A"]
    if int(cfg.get("angular_node_count", -1)) * int(cfg.get("radial_node_count", -1)) != len(dspec_a_rows):
        raise SystemExit("YSpec sample count does not match the DSpec point count")


def require_finite_samples(rows):
    required = ("radius_m", "theta_rad", "phi_rad", "kernel_objective_per_Pa_m3",
                "deltaX_Pa", "dr_m", "dOmega")
    for index, row in enumerate(rows):
        try:
            values = [float(row[field]) for field in required]
        except (KeyError, ValueError) as error:
            raise SystemExit(f"invalid sample row {index}: {error}")
        if not np.isfinite(values).all():
            raise SystemExit(f"non-finite value in sample row {index}; refusing to render or write report")


def col(rows, name):
    return np.array([r[name] for r in rows])


def key(row):
    return (row["modulus"], *(round(float(row[k]), 7) for k in ("radius_m", "theta_rad", "phi_rad")))


def pair_yspec_rows(dspec_rows, yspec_rows):
    """Pair duplicate interface coordinates by occurrence, preserving one-sided rows."""
    groups = {}
    for row in yspec_rows:
        groups.setdefault(key(row), []).append(row)
    offsets = {}
    paired = []
    for row in dspec_rows:
        k = key(row)
        index = offsets.get(k, 0)
        candidates = groups.get(k, [])
        if index >= len(candidates):
            raise SystemExit(f"YSpec rows do not match DSpec coordinate multiplicity at {k}")
        paired.append(candidates[index])
        offsets[k] = index + 1
    if any(offsets.get(k, 0) != len(rows) for k, rows in groups.items()):
        raise SystemExit("YSpec output has extra rows beyond the DSpec coordinate multiplicities")
    return paired


def grid_edges(centres, lower=None, upper=None):
    centres = np.asarray(centres, dtype=float)
    middle = (centres[:-1] + centres[1:]) / 2
    return np.concatenate(([lower if lower is not None else centres[0] - (middle[0] - centres[0])],
                           middle,
                           [upper if upper is not None else centres[-1] + (centres[-1] - middle[-1])]))


def angular_grid(rows, value_field, radius):
    row_radii = np.fromiter((float(r["radius_m"]) for r in rows), dtype=float,
                            count=len(rows))
    selected = [rows[i] for i in np.flatnonzero(np.isclose(row_radii, radius))]
    latitude = np.sort(np.unique(90.0 - np.rad2deg([float(r["theta_rad"]) for r in selected])))
    longitude = np.sort(np.unique(np.rad2deg([float(r["phi_rad"]) for r in selected])))
    lat_index = {round(value, 10): i for i, value in enumerate(latitude)}
    lon_index = {round(value, 10): i for i, value in enumerate(longitude)}
    values = np.full((len(latitude), len(longitude)), np.nan)
    for row in selected:
        lat = 90.0 - math.degrees(float(row["theta_rad"]))
        lon = math.degrees(float(row["phi_rad"]))
        values[lat_index[round(lat, 10)], lon_index[round(lon, 10)]] = float(row[value_field])
    if not np.isfinite(values).all():
        raise SystemExit("depth-map angular samples do not form a complete tensor grid")
    lon_edges = grid_edges(longitude, longitude[0] - (longitude[1] - longitude[0]) / 2,
                           longitude[-1] + (longitude[-1] - longitude[-2]) / 2)
    lat_edges = grid_edges(latitude, -90.0, 90.0)
    return longitude, latitude, lon_edges, lat_edges, values


def plot_depth_maps(data, out, compare=None):
    radii = np.sort(np.unique(col(data, "radius_m").astype(float)))
    target = float(radii[len(radii) // 2])
    for modulus in MODULI:
        d_rows = [r for r in data if r["modulus"] == modulus]
        lon, lat, lon_edges, lat_edges, dspec = angular_grid(
            d_rows, "kernel_objective_per_Pa_m3", target)
        yspec = difference = None
        limit = float(np.max(np.abs(dspec))) or 1.0
        if compare is not None:
            y_rows = pair_yspec_rows(d_rows, [r for r in compare if r["modulus"] == modulus])
            y_lon, y_lat, _, _, yspec = angular_grid(y_rows, "kernel_per_Pa_m3", target)
            if not np.allclose(lon, y_lon, rtol=0, atol=1e-8) or not np.allclose(lat, y_lat, rtol=0, atol=1e-8):
                raise SystemExit(f"YSpec {modulus} angular coordinates do not match DSpec")
            difference = dspec - yspec
            limit = max(limit, float(np.max(np.abs(yspec))) or 1.0)
        panels = 3 if compare is not None else 1
        fig, axes = plt.subplots(1, panels, figsize=(13.2 if panels == 3 else 5.8, 4.3),
                                 constrained_layout=True, squeeze=False)
        axes = axes[0]
        d_image = axes[0].pcolormesh(lon_edges, lat_edges, dspec, shading="flat",
                                     cmap="RdBu_r", vmin=-limit, vmax=limit)
        axes[0].set_title("DSpecM1D")
        images = [d_image]
        if compare is not None:
            images.append(axes[1].pcolormesh(lon_edges, lat_edges, yspec, shading="flat",
                                             cmap="RdBu_r", vmin=-limit, vmax=limit))
            diff_limit = float(np.max(np.abs(difference))) or 1.0
            diff_image = axes[2].pcolormesh(lon_edges, lat_edges, difference, shading="flat",
                                            cmap="RdBu_r", vmin=-diff_limit, vmax=diff_limit)
            images.append(diff_image)
            axes[1].set_title("YSpecpp")
            axes[2].set_title("DSpecM1D − YSpecpp")
            fig.colorbar(images[0], ax=axes[:2].tolist(), label="Kernel (objective Pa⁻¹ m⁻³); shared scale", pad=0.02)
            fig.colorbar(diff_image, ax=axes[2], label="Signed difference (objective Pa⁻¹ m⁻³)", pad=0.04)
        else:
            fig.colorbar(d_image, ax=axes[0], label="Kernel (objective Pa⁻¹ m⁻³)", pad=0.04)
        for ax in axes:
            ax.set_xlim(lon_edges[0], lon_edges[-1])
            ax.set_ylim(-90, 90)
            ax.set_xlabel("East longitude (deg)")
            ax.set_ylabel("Geocentric latitude (deg)")
            ax.set_xticks([0, 90, 180, 270, 360])
        fig.suptitle(f"{modulus} kernel at radius {target / 1000:.1f} km")
        for ext in ("png", "svg", "pdf"):
            fig.savefig(out / f"depth_map_{modulus}.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)


def bearing(lat1, lon1, lat2, lon2):
    dl = lon2 - lon1
    return math.atan2(math.sin(dl) * math.cos(lat2), math.cos(lat1)*math.sin(lat2) - math.sin(lat1)*math.cos(lat2)*math.cos(dl))


def section_grid(rows, plane, modulus):
    selected = [r for r in rows if r["plane"] == plane and r["modulus"] == modulus]
    distances = np.sort(np.unique([float(r["distance_m"]) / 1000 for r in selected]))
    radii = np.sort(np.unique([float(r["radius_m"]) / 1000 for r in selected]))
    distance_index = {round(v, 9): i for i, v in enumerate(distances)}
    radius_index = {round(v, 9): i for i, v in enumerate(radii)}
    values = np.full((len(radii), len(distances)), np.nan)
    for row in selected:
        d = float(row["distance_m"]) / 1000
        r = float(row["radius_m"]) / 1000
        values[radius_index[round(r, 9)], distance_index[round(d, 9)]] = float(row["kernel_objective_per_Pa_m3"])
    if not np.isfinite(values).all():
        raise SystemExit(f"{plane} {modulus} evaluations do not form a complete plane grid")
    return (grid_edges(distances), grid_edges(radii), values)


def sections(section_data, out):
    for modulus in MODULI:
        selected = [r for r in section_data if r["modulus"] == modulus]
        limit = float(max(abs(float(r["kernel_objective_per_Pa_m3"])) for r in selected)) or 1.0
        fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5), constrained_layout=True, sharey=True)
        image = None
        for ax, plane, title in zip(axes, ("great_circle", "cross_plane"),
                                    ("Great-circle plane", "Cross-plane")):
            x_edges, y_edges, values = section_grid(section_data, plane, modulus)
            image = ax.pcolormesh(x_edges, y_edges, values, shading="flat", cmap="RdBu_r",
                                  vmin=-limit, vmax=limit)
            ax.set_title(title)
            ax.set_xlabel("Signed distance from source (km)")
            ax.set_ylabel("Radius (km)")
        fig.colorbar(image, ax=axes.tolist(), label="Kernel (objective Pa⁻¹ m⁻³)", pad=0.03)
        fig.suptitle(f"{modulus} kernel on exact source-centred planes")
        for ext in ("png", "svg", "pdf"):
            fig.savefig(out / f"section_{modulus}.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)


def globe(data, out):
    all_radii = np.sort(np.unique(col(data, "radius_m").astype(float)))
    radii = all_radii[::4]
    if radii[-1] != all_radii[-1]:
        radii = np.append(radii, all_radii[-1])
    earth_radius = float(json.loads((out / "run.json").read_text())["radius_m"])
    wedge = (0.18, 0.88)
    for modulus in MODULI:
        rows = [r for r in data if r["modulus"] == modulus]
        scale = float(max(abs(float(r["kernel_objective_per_Pa_m3"])) for r in rows)) or 1.0
        norm = Normalize(-scale, scale)
        cmap = plt.get_cmap("RdBu_r")
        fig = plt.figure(figsize=(8.2, 7.3), constrained_layout=True)
        ax = fig.add_subplot(111, projection="3d")
        for radius in radii:
            lon, lat, _, _, values = angular_grid(rows, "kernel_objective_per_Pa_m3", float(radius))
            theta = np.deg2rad(90.0 - lat)
            phi = np.deg2rad(lon)
            t, p = np.meshgrid(theta, phi, indexing="ij")
            rr = float(radius)
            x = rr * np.sin(t) * np.cos(p)
            y = rr * np.sin(t) * np.sin(p)
            z = rr * np.cos(t)
            face_values = 0.25 * (values[:-1, :-1] + values[1:, :-1] +
                                  values[:-1, 1:] + values[1:, 1:])
            facecolors = cmap(norm(face_values))
            phi_mid = (phi[:-1] + phi[1:]) / 2
            facecolors[:, (phi_mid > wedge[0]) & (phi_mid < wedge[1]), 3] = 0.0
            ax.plot_surface(x, y, z, facecolors=facecolors, rstride=1, cstride=1,
                            linewidth=0, antialiased=False, shade=False)

        # A lightly drawn reference sphere gives the shell stack a geographic frame.
        tline = np.linspace(0, np.pi, 13)
        pline = np.linspace(0, 2 * np.pi, 25)
        tg, pg = np.meshgrid(tline, pline, indexing="ij")
        ax.plot_wireframe(earth_radius * np.sin(tg) * np.cos(pg),
                          earth_radius * np.sin(tg) * np.sin(pg),
                          earth_radius * np.cos(tg), color="0.35", alpha=0.22,
                          linewidth=0.45, rstride=2, cstride=2)
        meta = json.loads((out / "run.json").read_text())
        for label, geo, radius, color in (("S", meta["source"], float(meta["source"]["radius_m"]), "#e87522"),
                                          ("R", meta["receiver"], earth_radius - float(meta["receiver"]["depth_km"]) * 1000, "#138a72")):
            theta = math.radians(90.0 - float(geo["latitude_deg"]))
            phi = math.radians(float(geo["longitude_deg"]))
            xyz = radius * np.array([math.sin(theta)*math.cos(phi),
                                     math.sin(theta)*math.sin(phi), math.cos(theta)])
            ax.scatter(*xyz, s=54, c=color, edgecolors="black", linewidths=.5,
                       depthshade=False, label=("Source" if label == "S" else "Receiver"))
            ax.text(*(xyz * 1.025), label, color=color, fontsize=10, weight="bold")
        ax.set_box_aspect((1, 1, 1))
        ax.set_axis_off()
        ax.view_init(elev=24, azim=-55)
        ax.set_title(f"{modulus}: sampled shells {radii[0]/1000:.0f}–{radii[-1]/1000:.0f} km; open 40° wedge\nOutside sampled shells: unsampled")
        scalar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        fig.colorbar(scalar, ax=ax, label="Kernel (objective Pa⁻¹ m⁻³)", shrink=.7, pad=.02)
        ax.legend(loc="lower left", frameon=False)
        for ext in ("png", "svg", "pdf"):
            fig.savefig(out / f"globe_cutaway_{modulus}.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)


def radial_profile(rows, value_field, angular_weights=None):
    by_radius = {}
    for row in rows:
        by_radius.setdefault(float(row["radius_m"]), []).append(row)
    radii = np.asarray(sorted(by_radius), dtype=float)
    means = []
    for radius in radii:
        samples = by_radius[radius]
        values = np.asarray([float(r[value_field]) for r in samples])
        if angular_weights is None:
            weights = np.asarray([float(r["dOmega"]) for r in samples])
        else:
            weights = np.asarray([angular_weights[(round(float(r["theta_rad"]), 7),
                                                   round(float(r["phi_rad"]), 7))] for r in samples])
        means.append(float(np.sum(values * weights) / np.sum(weights)))
    return radii, np.asarray(means)


def radial_profiles(data, out, compare=None):
    dspec_angles = {}
    for row in data:
        dspec_angles.setdefault((row["modulus"], round(float(row["theta_rad"]), 7),
                                 round(float(row["phi_rad"]), 7)), float(row["dOmega"]))
    for modulus in MODULI:
        d_rows = [r for r in data if r["modulus"] == modulus]
        d_r, d_mean = radial_profile(d_rows, "kernel_objective_per_Pa_m3")
        y_r = y_mean = None
        if compare is not None:
            y_rows = pair_yspec_rows(d_rows, [r for r in compare if r["modulus"] == modulus])
            angular_weights = {}
            for row in y_rows:
                angle = (round(float(row["theta_rad"]), 7), round(float(row["phi_rad"]), 7))
                dspec_angle = (modulus, *angle)
                if dspec_angle not in dspec_angles:
                    raise SystemExit(f"YSpec {modulus} profile has an unmatched angular node")
                angular_weights[angle] = dspec_angles[dspec_angle]
            y_r, y_mean = radial_profile(y_rows, "kernel_per_Pa_m3", angular_weights)
            if not np.allclose(d_r, y_r, rtol=0, atol=1e-7):
                raise SystemExit(f"YSpec {modulus} radial coordinates do not match DSpec")
        fig, ax = plt.subplots(figsize=(5.8, 4.8), constrained_layout=True)
        ax.plot(d_mean, d_r / 1000, color="#315a83", lw=1.8, label="DSpecM1D")
        if compare is not None:
            ax.plot(y_mean, y_r / 1000, color="#bf4b3d", lw=1.2, ls="--",
                    marker="o", markersize=2.5, label="YSpecpp")
        ax.axvline(0, color="0.35", lw=.7)
        ax.grid(alpha=.2)
        ax.set_title(f"{modulus}: solid-angle weighted radial profile")
        ax.set_xlabel("Angular mean kernel (objective Pa⁻¹ m⁻³)")
        ax.set_ylabel("Radius (km)")
        if compare is not None:
            ax.legend(frameon=False)
        for ext in ("png", "svg", "pdf"):
            fig.savefig(out / f"radial_profile_{modulus}.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)


def source_split_plot(path, out):
    with path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if {r["modulus"] for r in rows} != {"C", "F", "L"}:
        raise SystemExit("source_split.csv must contain exactly C, F, and L")
    labels = ("Volume", "Contact", "Total")
    fig, ax = plt.subplots(figsize=(7.5, 4.6), constrained_layout=True)
    x = np.arange(3)
    width = 0.23
    for j, modulus in enumerate(("C", "F", "L")):
        row = next(r for r in rows if r["modulus"] == modulus)
        values = [float(row[k]) for k in ("volume_integral", "contact_integral", "total_integral")]
        ax.bar(x + (j - 1) * width, values, width, label=modulus)
    ax.axhline(0, color="0.25", lw=0.8)
    ax.set_xticks(x, labels)
    ax.set_ylabel("Contribution to J_Z (objective per unit fractional perturbation)")
    ax.set_title("Source-local volume, source-shell contact, and their sum")
    ax.legend(title="Modulus", ncols=3)
    ax.grid(axis="y", alpha=0.2)
    for ext in ("png", "svg", "pdf"):
        fig.savefig(out / f"source_split.{ext}", dpi=300)
    plt.close(fig)


def write_report(out):
    meta = json.loads((out / "run.json").read_text())
    integrals = json.loads((out / "integrals.json").read_text())
    render = meta.get("render_provenance", {})
    escape = lambda value: html.escape(str(value))
    bin_ids = meta.get("solver_bins", meta["frequency_bins_zero_based"])
    solved_hz = [int(k) * float(meta["df_hz"]) for k in bin_ids]
    frequency = (f"nt={meta['nt']}, dt={meta['dt_s']:.6g} s, df={meta['df_hz']:.8g} Hz, "
                 f"epsilon={meta['epsilon_s_inv']:.6g} s⁻¹; solved bins {bin_ids} "
                 f"({min(solved_hz):.8g}–{max(solved_hz):.8g} Hz)")
    src = meta["source"]
    rec = meta["receiver"]
    config_rows = [
        ("Objective", f"{meta['objective']} ({meta.get('objective_units', 'units recorded in metadata')})"),
        ("Solver band", f"lmax={meta['lmax']}; zero-based bins {meta['frequency_bins_zero_based']}"),
        ("Frequency grid", frequency),
        ("Time window", f"{meta['time_window']['type']}, samples {meta['time_window']['sample_start_inclusive']}–{meta['time_window']['sample_end_inclusive']} of {meta['time_window']['sample_count']}"),
        ("Source", f"{src['latitude_deg']}°, {src['longitude_deg']}°, depth {src['depth_km']} km; radius {src['radius_m']/1000:.3f} km"),
        ("Receiver", f"{rec['latitude_deg']}°, {rec['longitude_deg']}°, depth {rec['depth_km']} km"),
        ("Kernel coordinates", meta["coordinate_convention"]),
        ("Volume units", meta["kernel_units"]),
        ("Reference moduli", ", ".join(f"{m}={v:.6g} Pa" for m, v in zip(MODULI, meta["reference_moduli_Pa"]))),
        ("Model unit norms", json.dumps(meta["unit_norms"], sort_keys=True)),
        ("Physical G", f"{meta['physical_G_SI']:.8g} m³ kg⁻¹ s⁻²"),
        ("Convergence", "Bounded illustration; no convergence claim."),
    ]
    configuration = "\n".join(f"<tr><th>{escape(k)}</th><td>{escape(v)}</td></tr>" for k, v in config_rows)
    integral_rows = "\n".join(
        "<tr><td>{}</td><td>{:.7g}</td>{}</tr>".format(
            escape(row["modulus"]), row["dspec_volume_integral"],
            ("<td>{:.7g}</td><td>{:.4g}</td>".format(row["yspec_volume_integral"], row["relative_integral_error"])
             if "yspec_volume_integral" in row else "<td>—</td><td>—</td>"))
        for row in integrals["rows"])
    source_table = "<p>Source-local perturbation: {}</p>".format(escape(meta.get("source_local_profile", "not recorded")))
    split_path = out / "source_split.csv"
    if split_path.is_file():
        with split_path.open(newline="", encoding="utf-8") as f:
            split_rows = list(csv.DictReader(f))
        source_table += "<table><thead><tr><th>Modulus</th><th>Volume</th><th>Contact</th><th>Total</th></tr></thead><tbody>"
        source_table += "".join("<tr><td>{}</td><td>{:.7g}</td><td>{:.7g}</td><td>{:.7g}</td></tr>".format(
            escape(row["modulus"]), float(row["volume_integral"]), float(row["contact_integral"]),
            float(row["total_integral"])) for row in split_rows) + "</tbody></table>"
    comparison = "<p>No YSpecpp profile comparison was included.</p>"
    if any("yspec_volume_integral" in row for row in integrals["rows"]):
        comparison = ("<p>YSpecpp comparisons use matched saved volume coordinates. Integral error is "
                      "|I_D−I_Y|/max(|I_D|,|I_Y|); sample metrics are reported in integrals.json as "
                      "maximum absolute differences and volume-weighted RMS differences plus ratios to the YSpec reference norms. "
                      "The source-shell contact is excluded from YSpecpp sampled volume and shown separately above.</p>")
        yspec_inputs = render.get("yspec_comparison_inputs", {})
        if yspec_inputs:
            input_items = []
            for name, label in (("profiles", "YSpecpp profiles"), ("metadata", "YSpecpp run metadata")):
                item = yspec_inputs.get(name)
                if item:
                    path = Path(item["path"])
                    link = escape(path.resolve().as_uri())
                    input_items.append(
                        f'<li><a href="{link}">{label}</a>: {escape(str(path))} '
                        f'(SHA-256 {escape(item["sha256"])})</li>')
            if input_items:
                comparison += "<p>Comparison inputs (path and SHA-256):</p><ul>" + "".join(input_items) + "</ul>"
    figs = []
    for pattern in ("depth_map_*.png", "section_*.png", "globe_cutaway_*.png",
                    "radial_profile_*.png", "source_split.png"):
        for image_path in sorted(out.glob(pattern)):
            rel = image_path.name
            figs.append(f'<figure><a href="{escape(rel)}"><img src="{escape(rel)}" loading="lazy"></a><figcaption>{escape(rel)}</figcaption></figure>')
    prov = meta.get("provenance", {})
    source_meta = '<a href="source_run.json">original run metadata</a>' if (out / "source_run.json").exists() else '<a href="run.json">run metadata</a>'
    provenance = f"<p>Mode: {escape(render.get('mode', 'unknown'))}. {source_meta}. "
    if prov.get("model"):
        provenance += f"Model: {escape(prov['model'])} (SHA-256 {escape(prov.get('model_sha256', ''))}). "
    if prov.get("consumer"):
        provenance += f"Consumer: {escape(prov['consumer'])} (SHA-256 {escape(prov.get('consumer_sha256', ''))}). "
    provenance += ("Numeric file hashes and renderer hash are recorded in run.json. Files: "
                   "<a href=\"samples.csv\">samples.csv</a>, "
                   "<a href=\"sections.csv\">sections.csv</a>, "
                   "<a href=\"source_contacts.csv\">source_contacts.csv</a>, "
                   "<a href=\"source_split.csv\">source_split.csv</a>, "
                   "<a href=\"integrals.json\">integrals.json</a>.</p>")
    command = render.get("command", "")
    if command:
        command_html = f"<pre><code>{escape(command)}</code></pre>"
    else:
        command_html = "<p>Recreate figures with the command recorded in run.json render_provenance.</p>"
    document = f"""<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>3D kernel example report</title>
<style>body{{font:16px/1.5 system-ui,sans-serif;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#20252b}}h1,h2{{line-height:1.2}}table{{border-collapse:collapse;width:100%;margin:1rem 0 2rem}}th,td{{text-align:left;border-bottom:1px solid #ccd2d8;padding:.45rem .6rem;vertical-align:top}}th{{width:22%}}.figures{{display:grid;grid-template-columns:repeat(auto-fit,minmax(310px,1fr));gap:1rem}}figure{{margin:0;border:1px solid #ccd2d8;padding:.5rem;background:#fafbfc}}figure img{{width:100%;height:auto}}figcaption{{font-size:.88rem;color:#4a535c}}pre{{white-space:pre-wrap;background:#f1f3f5;padding:.8rem}}</style>
<body><h1>Windowed 3D material-kernel example</h1>
<p>{'Replayed from saved numeric inputs.' if render.get('mode') == 'replay' else 'Fresh bounded consumer run.'} The lmax={meta['lmax']} result is illustrative and has no convergence claim.</p>
<h2>Configuration</h2><table><tbody>{configuration}</tbody></table>
<h2>Volume integrals</h2><p>Directional perturbation integrals use the saved full-precision kernels and r² dr dΩ; units are objective per unit fractional perturbation.</p>
<table><thead><tr><th>Modulus</th><th>DSpecM1D</th><th>YSpecpp</th><th>Relative integral error</th></tr></thead><tbody>{integral_rows}</tbody></table>
<h2>Source-local volume and contact</h2>{source_table}
{comparison}
<h2>Figures</h2><div class="figures">{''.join(figs)}</div>
<h2>Reproduction and provenance</h2>{command_html}{provenance}
<p>Coordinates are geocentric colatitude and east longitude in radians in the data files. Source-shell contact uses dΩ; do not multiply by r².</p>
</body></html>"""
    (out / "report.html").write_text(document, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--consumer", type=Path)
    ap.add_argument("--model", type=Path)
    ap.add_argument("--replay", type=Path, help="existing consumer output directory; plot saved data without rerunning")
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--lmax", type=int, choices=(3, 12), default=3)
    ap.add_argument("--yspec-profiles", type=Path, help="optional full-precision profiles.csv from a matched YSpecpp kernel run")
    ap.add_argument("--yspec-meta", type=Path, help="required matched-kernel-timing-run-v1 meta.json for YSpec comparisons")
    a = ap.parse_args()
    if a.output.exists() and (not a.output.is_dir() or any(a.output.iterdir())):
        raise SystemExit(f"refusing nonempty output directory: {a.output}")
    if a.replay is not None:
        if a.consumer or a.model:
            raise SystemExit("--replay cannot be combined with --consumer or --model")
        a.output.mkdir(parents=True, exist_ok=True)
        import shutil
        source_run = a.replay / "run.json"
        if not source_run.is_file():
            raise SystemExit("replay directory must contain run.json")
        original_run_bytes = source_run.read_bytes()
        for name in ("samples.csv", "sections.csv", "source_contacts.csv", "source_split.csv", "run.json"):
            src = a.replay / name
            if src.is_file():
                shutil.copy2(src, a.output / name)
        if not (a.output / "samples.csv").is_file() or not (a.output / "sections.csv").is_file():
            raise SystemExit("replay directory must contain samples.csv and sections.csv")
        (a.output / "source_run.json").write_bytes(original_run_bytes)
        model_path = None
        consumer_path = None
    else:
        if not a.consumer or not a.model:
            raise SystemExit("provide --consumer and --model, or use --replay")
        a.output.mkdir(parents=True, exist_ok=True)
        try:
            subprocess.run([str(a.consumer.resolve()), str(a.model.resolve()),
                            str(a.output.resolve()), str(a.lmax)],
                           check=True, timeout=300)
        except subprocess.TimeoutExpired as error:
            raise SystemExit(
                "numerical consumer exceeded its 300-second limit; partial output is preserved"
            ) from error
        model_path = a.model.resolve()
        consumer_path = a.consumer.resolve()
    data = read(a.output / "samples.csv")
    require_finite_samples(data)
    section_data = read(a.output / "sections.csv")
    yspec = read_yspec(a.yspec_profiles) if a.yspec_profiles else None
    if a.yspec_meta and not a.yspec_profiles:
        raise SystemExit("--yspec-meta requires --yspec-profiles")
    run_meta = json.loads((a.output / "run.json").read_text(encoding="utf-8"))
    if yspec:
        validate_yspec_match(a.yspec_profiles, a.yspec_meta, run_meta, data,
                             a.output / "samples.csv")
    plot_depth_maps(data, a.output, yspec)
    sections(section_data, a.output)
    globe(data, a.output)
    radial_profiles(data, a.output, yspec)
    if (a.output / "source_split.csv").is_file():
        source_split_plot(a.output / "source_split.csv", a.output)
    integral_rows = []
    for modulus in MODULI:
        selected = [r for r in data if r["modulus"] == modulus]
        dspec_integral = sum(float(r["kernel_objective_per_Pa_m3"]) * float(r["deltaX_Pa"]) * float(r["radius_m"])**2 * float(r["dr_m"]) * float(r["dOmega"]) for r in selected)
        row = {"modulus": modulus, "dspec_volume_integral": dspec_integral}
        if yspec:
            y_rows = pair_yspec_rows(selected, [r for r in yspec if r["modulus"] == modulus])
            reference = np.array([float(r["kernel_per_Pa_m3"]) for r in y_rows])
            dspec_values = np.array([float(r["kernel_objective_per_Pa_m3"]) for r in selected])
            if not np.isfinite(reference).all():
                raise SystemExit(f"YSpec samples do not match the full DSpec volume grid for {modulus}")
            differences = dspec_values - reference
            reference_max = float(np.max(np.abs(reference)))
            volume_weights = np.array([float(r["radius_m"])**2 * float(r["dr_m"]) *
                                       float(r["dOmega"]) for r in selected])
            weight_sum = float(np.sum(volume_weights))
            if not np.isfinite(volume_weights).all() or weight_sum <= 0:
                raise SystemExit(f"invalid DSpec volume weights for {modulus}")
            reference_rms = float(np.sqrt(np.sum(reference**2 * volume_weights) / weight_sum))
            yspec_integral = sum(reference[i] * float(r["deltaX_Pa"]) * float(r["radius_m"])**2 * float(r["dr_m"]) * float(r["dOmega"]) for i, r in enumerate(selected))
            denominator = max(abs(dspec_integral), abs(yspec_integral))
            row["yspec_volume_integral"] = yspec_integral
            row["relative_integral_error"] = 0.0 if denominator == 0 else abs(dspec_integral-yspec_integral)/denominator
            row["max_abs_sample_difference"] = float(np.max(np.abs(differences)))
            row["max_abs_difference_over_reference_max"] = 0.0 if reference_max == 0 else row["max_abs_sample_difference"]/reference_max
            row["volume_weighted_rms_sample_difference"] = float(
                np.sqrt(np.sum(differences**2 * volume_weights) / weight_sum))
            row["rms_reference"] = reference_rms
            row["rms_difference_over_reference_rms"] = (
                0.0 if reference_rms == 0 else row["volume_weighted_rms_sample_difference"]/reference_rms)
        integral_rows.append(row)
    (a.output / "integrals.json").write_text(json.dumps({"units": "objective per unit fractional perturbation", "quadrature": "r^2 dr dOmega", "rows": integral_rows}, indent=2) + "\n")
    def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
    meta = json.loads((a.output / "run.json").read_text())
    if a.replay is None:
        meta["provenance"] = {"model": str(model_path), "model_sha256": sha(model_path),
                              "consumer": str(consumer_path), "consumer_sha256": sha(consumer_path),
                              "yspec_profiles": str(a.yspec_profiles.resolve()) if a.yspec_profiles else None,
                              "yspec_profiles_sha256": sha(a.yspec_profiles.resolve()) if a.yspec_profiles else None}
    meta["render_provenance"] = {
        "mode": "replay" if a.replay else "fresh_compute",
        "replay_directory": str(a.replay.resolve()) if a.replay else None,
        "renderer": str(Path(__file__).resolve()),
        "renderer_sha256": sha(Path(__file__).resolve()),
        "command": shlex.join([sys.executable, *sys.argv]),
        "numeric_inputs_sha256": {name: sha(a.output / name) for name in
                                  ("samples.csv", "sections.csv", "source_contacts.csv", "source_split.csv")
                                  if (a.output / name).is_file()},
        "yspec_comparison_inputs": ({
            "profiles": {"path": str(a.yspec_profiles.resolve()),
                         "sha256": sha(a.yspec_profiles.resolve())},
            "metadata": {"path": str(a.yspec_meta.resolve()),
                         "sha256": sha(a.yspec_meta.resolve())}}
            if a.yspec_profiles and a.yspec_meta else None)}
    (a.output / "run.json").write_text(json.dumps(meta, indent=2) + "\n")
    report_mode = "Replayed plots from saved numerical inputs; original run metadata is preserved in source_run.json. " if a.replay else "Fresh numerical run. "
    (a.output / "README.txt").write_text(report_mode + "Full precision samples.csv, sections.csv, source_contacts.csv, and source_split.csv are the numeric sources for the plots. Kernel signs are retained; solver maps share one zero-centred range and difference maps show a separately labelled signed range. Sections are evaluated on the exact planes. source_split.csv integrates the prescribed C/F/L source-local volume perturbation with r^2 dr dOmega and adds the separate contact term with dOmega. YSpecpp sampled volume kernels exclude source-shell contact; this comparison uses the separate DSpec source split. The lmax=12 preset is a bounded illustration, not a convergence study.\n")
    write_report(a.output)
    print(f"Wrote 300 dpi PNG, SVG, PDF, full-precision samples and provenance to {a.output}")


if __name__ == "__main__":
    main()
