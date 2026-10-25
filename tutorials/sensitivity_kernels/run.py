#!/usr/bin/env python3
"""Run the installed consumer, then make readable plots and an HTML report."""
import argparse
import csv
import hashlib
import html
import json
from pathlib import Path
import shlex
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def rows(path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--consumer", type=Path, required=True, help="installed consumer executable")
    ap.add_argument("--model", type=Path, required=True, help="1-D model file")
    ap.add_argument("--output", type=Path, required=True, help="new or empty results directory")
    args = ap.parse_args()
    if args.output.exists() and (not args.output.is_dir() or any(args.output.iterdir())):
        raise SystemExit(f"refusing nonempty output directory: {args.output}")
    args.output.mkdir(parents=True, exist_ok=True)
    subprocess.run([str(args.consumer.resolve()), str(args.model.resolve()), str(args.output.resolve())], check=True)
    samples, directions = rows(args.output / "samples.csv"), rows(args.output / "directions.csv")
    run_text = (args.output / "run.txt").read_text(encoding="utf-8")
    run_data = json.loads((args.output / "run.json").read_text(encoding="utf-8"))
    def sha256(path):
        h = hashlib.sha256()
        with path.open("rb") as f:
            for block in iter(lambda: f.read(1024 * 1024), b""):
                h.update(block)
        return h.hexdigest()
    run_data["provenance"] = {
        "model_file": str(args.model.resolve()), "model_sha256": sha256(args.model.resolve()),
        "consumer_file": str(args.consumer.resolve()), "consumer_sha256": sha256(args.consumer.resolve())}
    machine_directions, gate_rows = [], []
    for r in directions:
        vals = {k: float(v) for k, v in r.items() if k not in ("direction",)}
        pred = vals["predicted"]
        richardson = (4 * vals["fd_0.005"] - vals["fd_0.01"]) / 3
        denominator = max(abs(richardson), abs(pred))
        error = 0.0 if denominator == 0 else abs(richardson - pred) / denominator
        selected = bool(vals["selected_for_gate"])
        passed = error <= 1e-6
        machine_directions.append({"direction": r["direction"], **vals,
                                   "richardson_0.005_0.01": richardson,
                                   "relative_error": error, "gate_1e-6": passed,
                                   "selected_for_gate": selected})
        gate_rows.append((selected, passed))
    selected_pass = [passed for selected, passed in gate_rows if selected]
    run_data["directions"] = machine_directions
    run_data["directional_gate"] = {"name": "kernel vs Richardson directional derivative",
                                    "threshold": 1e-6, "all_selected_directions_pass": all(selected_pass),
                                    "selected_direction_count": len(selected_pass)}
    (args.output / "run.json").write_text(json.dumps(run_data, indent=2) + "\n", encoding="utf-8")

    radius = [float(r["radius_m"]) / 1000 for r in samples]
    rho = [float(r["density_kg_m3"]) for r in samples]
    solid = [r["solid"] == "1" for r in samples]
    c = [float(r["C_Pa"]) / 1e9 if is_solid else float("nan")
         for r, is_solid in zip(samples, solid)]
    kf = [float(r["fluid_bulk_pa"]) / 1e9 if not is_solid else float("nan")
          for r, is_solid in zip(samples, solid)]
    fig, ax = plt.subplots(3, 1, figsize=(8, 8), sharex=True, constrained_layout=True)
    ax[0].plot(radius, rho, color="#4466aa"); ax[0].set_ylabel(r"Density (kg m$^{-3}$)")
    ax[1].plot(radius, c, color="#b45a35"); ax[1].set_ylabel(r"Solid $C$ (GPa)")
    ax[2].plot(radius, kf, color="#26826e"); ax[2].set_ylabel(r"Fluid $\kappa$ (GPa)"); ax[2].set_xlabel("Radius (km)")
    for a in ax:
        for i in range(1, len(radius)):
            if solid[i] != solid[i-1]: a.axvline(radius[i], color="0.75", lw=.6)
        a.grid(alpha=.2)
    fig.suptitle("Reference model profiles on the NQ=4 SEM samples")
    fig.savefig(args.output / "profiles.png", dpi=160); plt.close(fig)

    kernels = ["density_kernel", "bulk_kernel", "C_kernel"]
    labels = ["Density kernel: J/((kg m$^{-3}$) m)", "Fluid bulk kernel: J/(Pa m)", "Solid C kernel: J/(Pa m)"]
    fig, ax = plt.subplots(3, 1, figsize=(8, 8), sharex=True, constrained_layout=True)
    for a, key, label in zip(ax, kernels, labels):
        vals = [float(r[key]) for r in samples]
        a.plot(radius, vals, lw=.9); a.set_ylabel(label, fontsize=9); a.grid(alpha=.2)
    ax[-1].set_xlabel("Radius (km)")
    fig.suptitle("Windowed measurement kernels on the one-sided SEM samples")
    fig.savefig(args.output / "kernels.png", dpi=160); plt.close(fig)

    names = [r["direction"] for r in directions]
    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    x = list(range(len(names)))
    ax.scatter(x, [float(r["predicted"]) for r in directions], marker="D", s=55,
               label="Integrated kernel prediction", color="#222222", zorder=3)
    for col, label, marker in [("fd_0.02", "Centered FD, a=0.02", "o"),
                               ("fd_0.01", "Centered FD, a=0.01", "s"),
                               ("fd_0.005", "Centered FD, a=0.005", "^")]:
        ax.scatter(x, [float(r[col]) for r in directions], label=label, marker=marker)
    ax.scatter(x, [r["richardson_0.005_0.01"] for r in machine_directions],
               label="Richardson estimate", marker="x", s=70, color="#aa3377")
    ax.set_xticks(x, names); ax.set_ylabel("Objective derivative per unit bump amplitude")
    ax.set_title("Localized fixed-mesh directional derivative check")
    ax.grid(axis="y", alpha=.25); ax.legend(fontsize=8)
    fig.savefig(args.output / "directional_check.png", dpi=160); plt.close(fig)

    body = "".join(f"<tr><td>{html.escape(r['direction'])}</td><td>{r['element']}</td><td>{float(r['r_lower_m'])/1000:.1f}–{float(r['r_upper_m'])/1000:.1f}</td>"+
                   "".join(f"<td>{float(r[k]):.4e}</td>" for k in ("predicted","fd_0.02","fd_0.01","fd_0.005"))+
                   f"<td>{machine_directions[i]['richardson_0.005_0.01']:.4e}</td><td>{machine_directions[i]['relative_error']:.3e} ({'PASS' if machine_directions[i]['gate_1e-6'] else 'FAIL'}; {'gate' if machine_directions[i]['selected_for_gate'] else 'screen'})</td></tr>" for i, r in enumerate(directions))
    command = " ".join(map(shlex.quote, (sys.executable, str(Path(__file__).resolve()),
        "--consumer", str(args.consumer.resolve()), "--model", str(args.model.resolve()),
        "--output", "/tmp/dspec-k2-example-rerun")))
    page = f"""<!doctype html><meta charset="utf-8"><title>K2 sensitivity-kernel example</title>
<style>body{{font:16px system-ui;max-width:980px;margin:2rem auto;color:#222}}img{{max-width:100%}}.scroll{{overflow-x:auto}}table{{border-collapse:collapse;white-space:nowrap}}td,th{{padding:.4rem .65rem;border:1px solid #bbb}}pre{{background:#f4f4f4;padding:.7rem;white-space:pre-wrap;overflow-wrap:anywhere}}code{{background:#f4f4f4;padding:.15rem;overflow-wrap:anywhere;word-break:break-word}}</style>
<h1>Windowed displacement sensitivity kernels</h1><p>This run compares material-kernel predictions with centered finite differences for localized solid-density and fluid-bulk bumps on a fixed SEM mesh. The selected-bin reconstruction is a windowed, selected-band measurement; it is not a broad-band research seismogram.</p>
<h2>Run configuration</h2><pre>{html.escape(run_text)}</pre><p>Model SHA-256: <code>{run_data['provenance']['model_sha256']}</code><br>Consumer SHA-256: <code>{run_data['provenance']['consumer_sha256']}</code></p>
<h2>Reproduce this run</h2><p>Using the recorded consumer and model, run the command below with a fresh output directory:</p><pre>{html.escape(command)}</pre>
<h2>Reference model profiles</h2><img src="profiles.png" alt="Density and solid/fluid bulk modulus profiles">
<h2>Kernel profiles</h2><img src="kernels.png" alt="Density, fluid bulk, and solid C kernel profiles">
<h2>Directional check</h2><img src="directional_check.png" alt="Kernel predictions and centered finite differences">
<div class="scroll"><table><thead><tr><th>Direction</th><th>Element</th><th>Radius interval (km)</th><th>Kernel</th><th>FD 0.02</th><th>FD 0.01</th><th>FD 0.005</th><th>Richardson</th><th>Relative error / gate</th></tr></thead><tbody>{body}</tbody></table></div>
<h2>Definitions and scope</h2><p>The scalar measurement is Re(sum(conj(q) y)). The measurement map is y=dt<sub>SI</sub> RFFT[W L<sub>m</sub> df̂ exp(ε̂ t̂) C2R(raw displacement spectrum)]; after the inverse real FFT and frequency-spacing factor, exponential damping is compensated once before SI length scaling and the Hann window. The window support is [0,(nt−1)dt<sub>SI</sub>]. The reported epsilon and z convention are in the run configuration. The passband and corner values are frequency-grid configuration metadata only: this engine applies no frequency taper or filter, and solves the listed bins only. Density kernels include the nonlocal change in background gravity from the perturbed density. Density kernels are objective/((kg m<sup>−3</sup>) m); fluid-bulk and solid-C kernels are objective/(Pa m). Predictions integrate kernel times physical perturbation and GLL dr weights. Centered differences are [J(+a)−J(−a)]/(2a), with a scaling the stated bump. The Richardson estimate is [4 FD(0.005)−FD(0.01)]/3. Relative error is |Richardson−prediction|/max(|Richardson|,|prediction|), defined as zero if both are zero. The diagnostic gate is 10<sup>−6</sup>. The 10 GPa fluid row is shown as an amplitude screen; the selected gate row is 100 GPa to increase the centered-difference signal. Every row remains visible in HTML and JSON. A failed selected gate leaves the evidence and makes the runner exit nonzero. This checks a fixed-mesh discrete derivative; it says nothing about mesh or bandwidth convergence.</p>
<p>Full-precision outputs: <a href="run.json">run.json</a>, <a href="run.txt">run.txt</a>, <a href="samples.csv">samples.csv</a>, and <a href="directions.csv">directions.csv</a>.</p>"""
    (args.output / "report.html").write_text(page, encoding="utf-8")
    print(f"Wrote report and plots in {args.output}")
    if not all(selected_pass):
        raise SystemExit("directional derivative gate failed; see report.html and run.json")


if __name__ == "__main__":
    main()
