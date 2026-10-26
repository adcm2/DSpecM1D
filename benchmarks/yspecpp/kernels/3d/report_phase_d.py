#!/usr/bin/env python3
"""Render a compact baseline-only report and diagnostic PNGs from saved outputs."""
from __future__ import annotations
import argparse
import html
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def write_report(out: Path) -> None:
    config = json.loads((out/'fixture.json').read_text())
    metrics = json.loads((out/'baseline_metrics.json').read_text())
    run = json.loads((out/'run.json').read_text())
    oracle = json.loads((out/'oracle/baseline.json').read_text())
    strain = metrics['cartesian_strain_fd_check']
    fd_path = out/'oracle/finite_differences.json'
    finite_differences = json.loads(fd_path.read_text()) if fd_path.is_file() else None
    if finite_differences is None:
        report_state = 'BASELINES_COMPLETE' if metrics['all_baseline_gates_pass'] else 'BASELINE_GATES_FAILED'
        report_stage = ('full campaign requested; this report contains baseline evidence only'
                        if run.get('requested_mode') == 'full' else 'baseline-only; finite differences not run')
    else:
        informative_checks = [row for support in finite_differences['supports'].values()
                              for directions in support['directions'].values()
                              for row in directions.values() if row.get('informative_kernel')]
        fd_pass = bool(informative_checks) and all(
            row.get('plateau', {}).get('informative') and row.get('plateau', {}).get('passed') and
            row.get('production_comparison', {}).get('informative') and
            row.get('production_comparison', {}).get('passed') and
            (not row.get('born_comparison', {}).get('informative') or
             row.get('born_comparison', {}).get('passed'))
            for row in informative_checks)
        report_state = 'PHASE_D_COMPLETE' if fd_pass else 'FINITE_DIFFERENCE_GATES_FAILED'
        report_stage = 'full finite-difference campaign; gate status derived from saved results'
    plots = out/'plots'
    plots.mkdir(exist_ok=True)

    spectrum = metrics['spectrum']
    labels = [f"{row['objective']} {row['component']}" for row in spectrum]
    values = [row['relative_l2'] if row['informative'] else row['scaled_absolute'] for row in spectrum]
    fig, ax = plt.subplots(figsize=(8.2, 4.4), constrained_layout=True)
    ax.bar(np.arange(len(labels)), values, color='#286a8e')
    ax.axhline(1e-5, color='#b34535', linestyle='--', linewidth=1.2, label='relative gate 1e-5')
    ax.set_yscale('log'); ax.set_ylabel('Relative L2 error; weak channels show difference / S')
    ax.set_xticks(np.arange(len(labels)), labels, rotation=25, ha='right')
    ax.set_title('Matched raw spectrum error, active bins 1–3')
    ax.grid(axis='y', alpha=.25); ax.legend(frameon=False)
    fig.savefig(plots/'baseline-spectrum-error.png', dpi=170); plt.close(fig)

    profiles = metrics['regular_profiles']
    source_local = metrics['source_local_production_born']
    labels = [f"{row['objective']} {row['modulus']}" for row in profiles]
    values = [row['Eprofile'] if row['profile_informative'] else row['RMSdiff_over_Kscale'] for row in profiles]
    fig, ax = plt.subplots(figsize=(9.4, 4.5), constrained_layout=True)
    ax.bar(np.arange(len(labels)), values, color='#4e8a62')
    ax.axhline(1e-3, color='#b34535', linestyle='--', linewidth=1.2, label='relative profile gate 1e-3')
    ax.set_yscale('log'); ax.set_ylabel('Relative weighted RMS; weak channels show RMS / Kscale')
    ax.set_xticks(np.arange(len(labels)), labels, rotation=35, ha='right')
    ax.set_title('Regular-volume kernel profile error')
    ax.grid(axis='y', alpha=.25); ax.legend(frameon=False)
    fig.savefig(plots/'regular-profile-error.png', dpi=170); plt.close(fig)

    values = [row['EJ'] if row['integrated_informative'] else row['abs_I_diff_over_S'] for row in profiles]
    fig, ax = plt.subplots(figsize=(9.4, 4.5), constrained_layout=True)
    ax.bar(np.arange(len(labels)), values, color='#9a7043')
    ax.axhline(1e-3, color='#b34535', linestyle='--', linewidth=1.2, label='relative integrated gate 1e-3')
    ax.set_yscale('log'); ax.set_ylabel('Relative integral error; weak channels show |ΔI| / S')
    ax.set_xticks(np.arange(len(labels)), labels, rotation=35, ha='right')
    ax.set_title('Integrated regular-volume kernel comparison')
    ax.grid(axis='y', alpha=.25); ax.legend(frameon=False)
    fig.savefig(plots/'regular-integral-error.png', dpi=170); plt.close(fig)

    born = metrics['oracle_born_comparisons']
    born_labels = [f"{row['objective']} {row['modulus']}" for row in born['rows']]
    born_values = [max(float(row['relative_error']), 1e-300) if row['informative']
                   else max(float(row['difference_over_S']), 1e-300) for row in born['rows']]
    fig, ax = plt.subplots(figsize=(9.4, 4.6), constrained_layout=True)
    ax.bar(np.arange(len(born_labels)), born_values, color='#9b4d67')
    ax.axhline(born['relative_gate'], color='#b34535', linestyle='--',
               linewidth=1.2, label=f"informative Born gate {born['relative_gate']:.0e}")
    ax.set_yscale('log'); ax.set_ylabel('Relative error; weak rows show absolute difference / S')
    ax.set_xticks(np.arange(len(born_labels)), born_labels, rotation=35, ha='right')
    ax.set_title('Independent regular-volume Born prediction vs DSpec kernel integral')
    ax.grid(axis='y', alpha=.25); ax.legend(frameon=False)
    fig.savefig(plots/'regular-born-error.png', dpi=170); plt.close(fig)

    bridge = oracle['baseline_spectrum_agreement']['relative_l2_by_quantity_component']
    bridge_labels = list(bridge)
    bridge_values = [max(float(bridge[key]), 1e-300) for key in bridge_labels]
    fig, ax = plt.subplots(figsize=(9.2, 4.5), constrained_layout=True)
    ax.bar(np.arange(len(bridge_labels)), bridge_values, color='#7653a6')
    ax.axhline(config['gates']['oracle_baseline_relative'], color='#b34535',
               linestyle='--', linewidth=1.2, label='oracle bridge gate 1e-9')
    ax.set_yscale('log'); ax.set_ylabel('Relative L2 error over the full 65-bin grid')
    ax.set_xticks(np.arange(len(bridge_labels)), bridge_labels, rotation=30, ha='right')
    ax.set_title('Independent Cartesian baseline spectrum bridge')
    ax.grid(axis='y', alpha=.25); ax.legend(frameon=False)
    fig.savefig(plots/'oracle-spectrum-bridge.png', dpi=170); plt.close(fig)

    gate_rows = []
    for row in spectrum:
        metric = f"{row['relative_l2']:.3g}" if row['informative'] else f"{row['scaled_absolute']:.3g} (scaled absolute)"
        gate_rows.append(f"<tr><td>{html.escape(row['objective'])}</td><td>{html.escape(row['component'])}</td><td>{metric}</td><td>{'PASS' if row['gate_pass'] else 'weak diagnostic'}</td></tr>")
    profile_rows = []
    for row in profiles:
        ep = f"{row['Eprofile']:.3g}" if row['profile_informative'] else f"{row['RMSdiff_over_Kscale']:.3g} (scaled absolute)"
        ei = f"{row['EJ']:.3g}" if row['integrated_informative'] else f"{row['abs_I_diff_over_S']:.3g} (scaled absolute)"
        profile_rows.append(f"<tr><td>{html.escape(row['objective'])}</td><td>{html.escape(row['modulus'])}</td><td>{ep}</td><td>{ei}</td><td>{'PASS' if row['profile_gate_pass'] and row['integrated_gate_pass'] else 'weak diagnostic'}</td></tr>")
    oracle_rows = []
    for key, value in bridge.items():
        oracle_rows.append(f"<tr><td>{html.escape(key)}</td><td>{float(value):.3g}</td><td>{'PASS' if float(value) <= config['gates']['oracle_baseline_relative'] else 'FAIL'}</td></tr>")
    objective_bridge_rows = []
    for objective, values in metrics['oracle_baseline']['objective_checks_normalized_by_S'].items():
        for quantity, check in values.items():
            objective_bridge_rows.append(f"<tr><td>{html.escape(objective)}</td><td>{html.escape(quantity)}</td><td>{float(check['absolute_difference_over_S']):.3g}</td><td>{'PASS' if check['gate_pass'] else 'FAIL'}</td></tr>")
    born_rows = []
    for objective, directions in oracle['regular_volume_born_m_s_per_h'].items():
        for modulus, value in directions.items():
            born_rows.append(f"<tr><td>{html.escape(objective)}</td><td>{html.escape(modulus)}</td><td>{float(value):.3g}</td></tr>")
    born_comparison_rows = []
    for row in born['rows']:
        error = (f"{row['relative_error']:.3g}" if row['informative']
                 else f"{row['difference_over_S']:.3g} (difference / S)")
        born_comparison_rows.append(
            f"<tr><td>{html.escape(row['objective'])}</td><td>{html.escape(row['modulus'])}</td>"
            f"<td>{row['regular_born_m_s_per_h']:.3g}</td><td>{row['dspec_kernel_integral_m_s']:.3g}</td>"
            f"<td>{error}</td><td>{'PASS' if row['gate_pass'] else ('weak diagnostic' if not row['informative'] else 'FAIL')}</td></tr>")
    source_local_rows = []
    for row in source_local['rows']:
        elements = ', '.join(str(x) for x in row['elements'])
        source_local_rows.append(
            f"<tr><td>{html.escape(row['objective'])}</td><td>{html.escape(row['modulus'])}</td>"
            f"<td>{row['I_volume_m_s_per_h']:.3g}</td><td>{row['I_contact_m_s_per_h']:.3g}</td>"
            f"<td>{row['Kpred_total_m_s_per_h']:.3g}</td><td>{row['Kscale_per_Pa_m3']:.3g}</td>"
            f"<td>{html.escape(elements)}</td></tr>")
    source_bridge = metrics['source_local_born_comparisons']
    source_bridge_rows = []
    for row in source_bridge['rows']:
        error = (f"{row['relative_error']:.3g}" if row['informative']
                 else f"{row['scaled_absolute_error']:.3g} (difference / S)")
        source_bridge_rows.append(
            f"<tr><td>{html.escape(row['objective'])}</td><td>{html.escape(row['modulus'])}</td>"
            f"<td>{html.escape(row['term'])}</td><td>{row['cartesian_born_m_s_per_h']:.3g}</td>"
            f"<td>{row['dspec_integral_m_s_per_h']:.3g}</td><td>{error}</td>"
            f"<td>{'PASS' if row['gate_pass'] else 'FAIL'}</td></tr>")
    source_bridge_labels = [f"{r['objective']} {r['modulus']} {r['term']}"
                            for r in source_bridge['rows']]
    source_bridge_values = [max(float(r['relative_error'] or r['scaled_absolute_error']), 1e-300)
                            for r in source_bridge['rows']]
    fig, ax = plt.subplots(figsize=(10.2, 4.8), constrained_layout=True)
    ax.bar(np.arange(len(source_bridge_labels)), source_bridge_values, color='#a56c4a')
    ax.axhline(source_bridge['relative_gate'], color='#b34535', linestyle='--',
               linewidth=1.2, label=f"Born bridge gate {source_bridge['relative_gate']:.0e}")
    ax.set_yscale('log'); ax.set_ylabel('Relative error; weak rows show difference / S')
    ax.set_xticks(np.arange(len(source_bridge_labels)), source_bridge_labels,
                  rotation=50, ha='right')
    ax.set_title('Independent source-local Born split vs DSpec production integrals')
    ax.grid(axis='y', alpha=.25); ax.legend(frameon=False)
    fig.savefig(plots/'source-local-born-bridge.png', dpi=170); plt.close(fig)

    fd_section = ''
    if finite_differences is not None:
        fd_rows, fd_plot_labels, fd_plot_values, fd_plot_pass = [], [], [], []
        split_rows, residual_rows, resource_rows = [], [], []
        for support, data in finite_differences['supports'].items():
            for modulus, objectives in data['directions'].items():
                for objective, row in objectives.items():
                    plateau = row.get('plateau') or {}
                    comparison = row.get('born_comparison') or {}
                    production_comparison = row.get('production_comparison') or {}
                    median = plateau.get('median_m_s_per_h')
                    spread = plateau.get('spread_over_max_abs', plateau.get('spread_over_S'))
                    born_error = (comparison.get('relative_error') if comparison.get('informative')
                                  else comparison.get('scaled_absolute_error'))
                    production_error = (production_comparison.get('relative_error')
                                        if production_comparison.get('informative')
                                        else production_comparison.get('scaled_absolute_error'))
                    production_kpred = row.get('production_Kpred_m_s_per_h')
                    status = row.get('informative_class', row.get('status', 'computed'))
                    if plateau.get('passed') is not None:
                        status += '; plateau PASS' if plateau['passed'] else '; plateau FAIL'
                    if production_comparison.get('passed') is not None:
                        if production_comparison.get('informative'):
                            status += ('; production relative PASS' if production_comparison['passed']
                                       else '; production relative FAIL')
                        else:
                            status += ('; weak absolute diagnostic within limit'
                                       if production_comparison['passed']
                                       else '; weak absolute diagnostic over limit')
                    if comparison.get('passed') is not None:
                        if comparison.get('informative'):
                            status += '; oracle Born relative PASS' if comparison['passed'] else '; oracle Born relative FAIL'
                        else:
                            status += '; oracle Born weak diagnostic'
                    production_error_text = ('—' if production_error is None else
                        (f'{float(production_error):.3g}' if production_comparison.get('informative')
                         else f'{float(production_error):.3g} / S (weak)'))
                    born_error_text = ('—' if born_error is None else
                        (f'{float(born_error):.3g}' if comparison.get('informative')
                         else f'{float(born_error):.3g} / S (weak)'))
                    fd_rows.append(
                        f"<tr><td>{html.escape(support)}</td><td>{html.escape(objective)}</td>"
                        f"<td>{html.escape(modulus)}</td><td>{'—' if production_kpred is None else f'{float(production_kpred):.3g}'}</td>"
                        f"<td>{'—' if median is None else f'{float(median):.3g}'}</td>"
                        f"<td>{'—' if spread is None else f'{float(spread):.3g}'}</td>"
                        f"<td>{production_error_text}</td><td>{born_error_text}</td>"
                        f"<td>{html.escape(status)}</td></tr>")
                    if median is not None and production_comparison.get('informative'):
                        fd_plot_labels.append(f"{support} {objective} {modulus}")
                        fd_plot_values.append(max(float(production_error), 1e-300))
                        fd_plot_pass.append(bool(production_comparison.get('passed')))
                    split = row.get('source_prediction_split_m_s_per_h')
                    if support == 'source_local' and split:
                        split_rows.append(
                            f"<tr><td>{html.escape(objective)}</td><td>{html.escape(modulus)}</td>"
                            f"<td>{float(split['oracle_volume']):.3g}</td><td>{float(split['production_volume']):.3g}</td>"
                            f"<td>{float(split['oracle_contact']):.3g}</td><td>{float(split['production_contact']):.3g}</td>"
                            f"<td>{float(split['production_total']):.3g}</td></tr>")
        for support, cases in finite_differences.get('solver_residuals', {}).items():
            residuals, iterations, adjoint = [], [], []
            for bins in cases.values():
                for solve in bins.values():
                    residuals.append(float(solve['true_relative_residual']))
                    iterations.append(int(solve['iterations']))
                    for objective in ('J_Z', 'J_N'):
                        if objective in solve:
                            adjoint.append(float(solve[objective]['lambda0_adjoint_residual_over_S']))
            if residuals:
                residual_rows.append(
                    f"<tr><td>{html.escape(support)}</td><td>{max(residuals):.3g}</td>"
                    f"<td>{max(iterations)}</td><td>{max(adjoint, default=0.0):.3g}</td>"
                    f"<td>{float(finite_differences['solver_controls']['true_relative_residual_gate']):.3g}</td></tr>")
        for support, resource_record in finite_differences.get('resources', {}).items():
            resource_rows.append(
                f"<tr><td>{html.escape(support)}</td><td>{int(resource_record['total_dofs']):,}</td>"
                f"<td>{int(resource_record['delta_upper_nnz']):,}</td>"
                f"<td>{resource_record['estimated_peak_bytes']/1e6:.3g} MB</td>"
                f"<td>{resource_record['memory_cap_bytes']/1e9:.3g} GB</td></tr>")
        if fd_plot_values:
            fig, ax = plt.subplots(figsize=(10.2, 4.8), constrained_layout=True)
            ax.bar(np.arange(len(fd_plot_labels)), fd_plot_values,
                   color=['#4e8a62' if passed else '#b34535' for passed in fd_plot_pass])
            ax.axhline(float(config['gates']['independent_born_relative']), color='#b34535',
                       linestyle='--', linewidth=1.2, label='production Kpred gate 1e-6')
            ax.set_yscale('log'); ax.set_ylabel('Production Kpred vs Richardson plateau error')
            ax.set_xticks(np.arange(len(fd_plot_labels)), fd_plot_labels,
                          rotation=55, ha='right')
            ax.set_title('Production kernel prediction vs Richardson plateau')
            ax.grid(axis='y', alpha=.25); ax.legend(frameon=False)
            fig.savefig(plots/'finite-difference-plateau-error.png', dpi=170); plt.close(fig)
        fd_section = (
            '<h2>Finite-difference and Richardson results</h2>'
            f"<p>{html.escape(finite_differences.get('status', 'Saved finite-difference packet'))}. "
            f"Steps h={', '.join(f'{float(h):g}' for h in finite_differences['steps_h'])}; "
            f"the finest consecutive passing Richardson triple is selected using only R and S. "
            f"The direct production Kpred comparison uses the frozen 1e-6 relative gate when both "
            f"production derivative and plateau exceed 1e-6 S; the independent Cartesian Born "
            f"comparison is retained as a separate check. "
            f"Solver controls and true-residual evidence are saved in the packet.</p>"
            '<div class="table-wrap"><table><thead><tr><th>Support</th><th>Objective</th><th>Modulus</th>'
            '<th>Production Kpred</th><th>Plateau median</th><th>Plateau spread</th>'
            '<th>Production error</th><th>Oracle Born error</th><th>Status</th>'
            f"</tr></thead><tbody>{''.join(fd_rows)}</tbody></table></div>"
            + ('<img src="plots/finite-difference-plateau-error.png" alt="Born prediction errors against Richardson plateaus">'
               if fd_plot_values else '')
            + ('<h3>Source-local volume and contact split</h3><div class="table-wrap"><table><thead><tr>'
               '<th>Objective</th><th>Modulus</th><th>Oracle volume</th><th>DSpec volume</th>'
               '<th>Oracle contact</th><th>DSpec contact</th><th>DSpec total</th></tr></thead><tbody>'
               f"{''.join(split_rows)}</tbody></table></div>" if split_rows else '')
            + ('<h3>Coupled-solve residuals and sparse resources</h3><div class="table-wrap"><table><thead><tr>'
               '<th>Support</th><th>Maximum true relative residual</th><th>Maximum GMRES iterations</th>'
               '<th>Maximum |λ₀ᴴr|/S diagnostic</th><th>Residual limit</th></tr></thead><tbody>'
               f"{''.join(residual_rows)}</tbody></table></div>" if residual_rows else '')
            + ('<div class="table-wrap"><table><thead><tr><th>Support</th><th>Total DOFs</th>'
               '<th>Upper ΔH nonzeros</th><th>Estimated peak</th><th>Memory cap</th></tr></thead><tbody>'
               f"{''.join(resource_rows)}</tbody></table></div>" if resource_rows else '')
            + '<p>Per-step objectives, centered derivatives, all Richardson values, and coupling norms are retained in '
              '<a href="oracle/finite_differences.json">finite_differences.json</a>. '
              + ('The append-only solve trace is available in '
                 '<a href="oracle/finite_difference_progress.jsonl">finite_difference_progress.jsonl</a>.'
                 if (out/'oracle/finite_difference_progress.jsonl').is_file() else '')
              + '</p>')
    coupling = oracle['coupling_structure']
    coupling_rows = []
    coupling_groups = (
        ('angular Cartesian Gram', coupling['angular_cartesian_gram'],
         'spheroidal_toroidal_block_norm', 'offdiagonal_lm_block_norm'),
        ('regular ΔH', coupling['assembled_regular_delta_h'],
         'spheroidal_toroidal_frobenius_norm', 'offdiagonal_lm_frobenius_norm'),
        ('source-local ΔH', coupling['assembled_source_local_delta_h'],
         'spheroidal_toroidal_frobenius_norm', 'offdiagonal_lm_frobenius_norm'))
    for label, group, st_key, offdiag_key in coupling_groups:
        for modulus in ('L', 'N'):
            values = group[modulus]
            nnz = (f"{int(values['assembled_delta_h_nnz']):,}"
                   if 'assembled_delta_h_nnz' in values else '—')
            coupling_rows.append(
                f"<tr><td>{html.escape(label)}</td><td>{modulus}</td>"
                f"<td>{float(values[st_key]):.3g}</td><td>{float(values[offdiag_key]):.3g}</td>"
                f"<td>{nnz}</td></tr>")
    resource = oracle['resource_estimate']
    provenance = run.get('provenance', {})
    origin = html.escape(str(provenance.get('origin_preflight_sha256', config.get('origin_preflight_sha256'))))
    report = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Phase D matched baselines</title><style>
body{{font:16px/1.5 system-ui,sans-serif;color:#17232c;margin:0;background:#f4f6f7}}main{{max-width:1080px;margin:auto;padding:24px}}
h1,h2{{line-height:1.2}}.status{{background:#fff2cf;border-left:5px solid #bd7b00;padding:14px 18px;margin:18px 0}}
.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(280px,1fr));gap:14px}}.card{{background:white;border:1px solid #d5dde1;border-radius:8px;padding:16px}}
table{{border-collapse:collapse;width:100%;min-width:760px;background:white;margin:0}}.table-wrap{{max-width:100%;overflow-x:auto;margin:10px 0 24px}}th,td{{padding:8px 10px;border-bottom:1px solid #dce2e5;text-align:left;vertical-align:top}}th{{background:#e8eef1}}
img{{max-width:100%;height:auto;border:1px solid #d5dde1;background:white}}code{{overflow-wrap:anywhere}}a{{color:#075a82}}
@media(max-width:640px){{main{{padding:14px}}table{{font-size:14px}}th,td{{padding:6px}}}}
</style></head><body><main>
<h1>Phase D matched baselines</h1>
<div class="status"><strong>{'Finite-difference packet available.' if finite_differences else 'Baseline-only evidence.'}</strong> {'The saved packet includes six ±h pairs and the frozen Richardson selection.' if finite_differences else 'This report contains no ±h finite differences or Richardson convergence result; it is not a completed Phase D campaign.'}</div>
<p>Fixture: homogeneous elastic sphere, radius {config['model']['radius_m']/1000:g} km, density {config['model']['rho_kg_m3']:g} kg/m³, Vp {config['model']['vp_m_s']:g} m/s, Vs {config['model']['vs_m_s']:g} m/s. Degrees {config['solver']['degrees_inclusive'][0]}–{config['solver']['degrees_inclusive'][1]}, active raw bins 1–3, nt={config['solver']['nt']}, dt={config['solver']['dt_s']:g} s, Hann window on samples 0–90. Objectives J_Z and J_N are separate.</p>
<p>Frequency is z_k=2πk/(nt·dt)−iε with ε={config['solver']['epsilon_s_inv']:.6g} s⁻¹. Raw displacement maps through L₀T₀f_z once; full-grid P/P* and real-RFFT DC/Nyquist weights then form the separate Z and N objectives. Raw solver spectra are zero outside bins 1–3. Objective and signal scale S use m·s. The angular perturbation is χ=Re(Y₂₁)+0.5Re(Y₃₂), using normalized complex scalar harmonics, geocentric colatitude, and east longitude.</p>
<p>The regular-volume profile is sin²(π(r−4200 km)/2000 km) on [4200,6200] km, with exact zero endpoints and zero at the source. The source-local profile is cos²(π(r−6321 km)/20 km) on [6311,6331] km. Volume kernels use objective/(Pa·m³) and physical dV=r²dr dΩ. YSpec sampled volume kernels test the regular profile only; source-contact terms are kept separate.</p>
<div class="grid"><section class="card"><h2>Scientific status</h2><p>{html.escape(report_state)}<br>Stage: {html.escape(report_stage)}<br>Baseline gates: {'PASS' if metrics['all_baseline_gates_pass'] else 'FAIL'}</p><p>Status is derived from saved baseline and finite-difference gates. See <a href="run.json">run.json for execution status and provenance</a>.</p><p>Frozen preflight SHA-256:<br><code>{origin}</code></p></section>
<section class="card"><h2>Source and receiver</h2><p>Source: 10° geocentric, 20° east, 50 km depth (r=6321 km). Receiver: 45° geocentric, 90° east, 3 km depth (r=6368 km).</p><p>YSpec uses geographic latitude conversion from fixture.json with round-trip checking.</p></section></div>
<h2>Matched spectrum gate</h2><p>Each component is compared over bins 1–3. The relative L2 limit is 1e-5 when both signals exceed 1e-6 S; weak signals are shown as difference/S and are not counted as relative passes.</p>
<div class="table-wrap"><table><thead><tr><th>Objective</th><th>Component</th><th>Error</th><th>Status</th></tr></thead><tbody>{''.join(gate_rows)}</tbody></table></div>
<img src="plots/baseline-spectrum-error.png" alt="Matched spectrum errors by objective and component">
<h2>Independent Cartesian baseline bridge</h2><p>The separate sparse modal solve uses {oracle['mode_count']} modes and {oracle['total_dofs']} total degrees of freedom. Before producing Born predictions, its raw receiver spectra and P-processed spectra are compared component by component with DSpec over the full grid. The predeclared limit is {config['gates']['oracle_baseline_relative']:.0e}. The J0 and S consistency checks report |oracle−DSpec|/S_DSpec against the same limit, so a small objective from cancellation does not set its own scale. An independent 12-case Cartesian displacement finite-difference check validates the spherical strain conversion: maximum absolute component error {strain['max_absolute_component_error']:.3g}, against {strain['tolerance']:.0e}. {'No modal perturbed solve was run.' if not finite_differences else 'The perturbed modal solves and their true residuals are included below.'}</p>
<div class="table-wrap"><table><thead><tr><th>Spectrum</th><th>Relative L2 error</th><th>Status</th></tr></thead><tbody>{''.join(oracle_rows)}</tbody></table></div>
<img src="plots/oracle-spectrum-bridge.png" alt="Independent Cartesian baseline spectrum bridge errors">
<div class="table-wrap"><table><thead><tr><th>Objective</th><th>Quantity</th><th>|difference| / S</th><th>Status</th></tr></thead><tbody>{''.join(objective_bridge_rows)}</tbody></table></div>
<p>The saved oracle packet also contains full-grid objective and S values, sparse workspace estimates, solver controls, and the regular-volume first-order predictions below. Prediction units are m·s per unit fractional perturbation.</p>
<div class="table-wrap"><table><thead><tr><th>Objective</th><th>Modulus</th><th>Linear prediction (m·s per h)</th></tr></thead><tbody>{''.join(born_rows)}</tbody></table></div>
<p>Source projection has scalar Hermitian error {float(oracle['source_projection_checks']['scalar_hermitian_error']):.3g}, vector Hermitian error {float(oracle['source_projection_checks']['vector_hermitian_error']):.3g}, and S–T source-vector coupling norm {float(oracle['source_projection_checks']['source_vector_ST_coupling_norm']):.3g}. Both the angular Gram and assembled sparse ΔH operators show S–T and off-diagonal l,m coupling for L and N:</p>
<div class="table-wrap"><table><thead><tr><th>Operator</th><th>Modulus</th><th>Spheroidal–toroidal norm</th><th>Off-diagonal l,m norm</th><th>ΔH nonzeros</th></tr></thead><tbody>{''.join(coupling_rows)}</tbody></table></div>
<p>Sparse resource estimate: {int(resource['total_dofs']):,} DOFs, {int(resource['delta_upper_nnz']):,} upper ΔH entries, estimated peak {resource['estimated_peak_bytes']/1e6:.3g} MB under the {resource['memory_cap_bytes']/1e9:.3g} GB cap. The global system is assembled and solved in sparse blocks; no dense global factorization is used.</p>
<h2>Regular-volume Born comparison</h2><p>Each independent linear prediction is compared with the DSpec sampled-kernel integral for the same perturbation and objective. Both use the source-vanishing profile over the full [4200,6200] km support. The relative gate is {born['relative_gate']:.0e}; a row is informative only when both derivative magnitudes exceed {born['informative_signal_floor_over_S']:.0e} S. Weak rows show |difference|/S as a diagnostic.</p>
<div class="table-wrap"><table><thead><tr><th>Objective</th><th>Modulus</th><th>Oracle dJ/dh</th><th>DSpec integral</th><th>Error</th><th>Status</th></tr></thead><tbody>{''.join(born_comparison_rows)}</tbody></table></div>
<img src="plots/regular-born-error.png" alt="Independent regular-volume Born prediction error by objective and modulus">
<h2>Regular-volume kernel gate</h2><p>For A, C, F, L, and N separately, δX/X₀=h·sin²(π(r−4200 km)/2000 km)·χ on 4200–6200 km, with exact zero endpoints. Profile weighted RMS and integrated derivative limits are each 1e-3. YSpec sampled volume kernels omit source shell/contact terms, so these comparisons apply only to the source-vanishing profile.</p>
<div class="table-wrap"><table><thead><tr><th>Objective</th><th>Modulus</th><th>Profile error</th><th>Integrated error</th><th>Status</th></tr></thead><tbody>{''.join(profile_rows)}</tbody></table></div>
<div class="grid"><img src="plots/regular-profile-error.png" alt="Regular profile kernel errors"><img src="plots/regular-integral-error.png" alt="Regular profile integrated derivative errors"></div>
<h2>Source-local production Born split</h2><p>DSpec point kernels are integrated over the closed support {source_local['support_m'][0]/1000:g}–{source_local['support_m'][1]/1000:g} km using <code>{html.escape(source_local['profile'])}</code>. The reported volume term includes the upper source-owning element’s local lift; the contact term is integrated separately over χ at the source. Their sum is the production Kpred diagnostic, not a finite-difference result or a pass/fail gate. All derivatives are in m·s per unit fractional h. Sampled elements include {', '.join(str(x) for x in next(r['elements'] for r in source_local['rows']))}; the source-owning element is {source_local['source_element']}.</p>
<div class="table-wrap"><table><thead><tr><th>Objective</th><th>Modulus</th><th>Volume dJ/dh</th><th>Contact dJ/dh</th><th>Total Kpred</th><th>Kscale</th><th>Elements</th></tr></thead><tbody>{''.join(source_local_rows)}</tbody></table></div>
<p>The independent Cartesian oracle decomposes the same source-local prediction into volume, contact, and total terms and compares all 18 objective/modulus/term values against DSpec production. The source-owning element lift is included in the volume term. Limit: {source_bridge['relative_gate']:.0e}; informative floor: {source_bridge['informative_signal_floor_over_S']:.0e} S.</p>
<div class="table-wrap"><table><thead><tr><th>Objective</th><th>Modulus</th><th>Term</th><th>Cartesian Born</th><th>DSpec integral</th><th>Error</th><th>Status</th></tr></thead><tbody>{''.join(source_bridge_rows)}</tbody></table></div>
<img src="plots/source-local-born-bridge.png" alt="Source-local independent Born bridge errors">
{fd_section}
<h2>Reproducibility</h2><p>Command: <code>python3 benchmarks/yspecpp/kernels/3d/run_phase_d.py --output &lt;new-or-empty-directory&gt; --full</code> for the complete baseline and finite-difference campaign, or <code>--baseline-only</code> for baseline evidence. The runner records source and dependency identities, exact build/run commands and exit statuses, build cache fingerprints, and executable hashes in run.json. Files: <a href="fixture.json">resolved fixture</a>, <a href="baseline_metrics.json">full precision baseline metrics</a>, <a href="oracle/baseline.json">independent oracle baseline</a>, {('<a href="oracle/finite_differences.json">finite-difference packet</a>, ' if finite_differences else '')}<a href="source_snapshot/manifest.json">source snapshot manifest</a>, <a href="run.json">provenance and commands</a>, <a href="dspec/">DSpec outputs</a>, <a href="yspecpp/">YSpec outputs</a>.</p>
</main></body></html>'''
    (out/'report.html').write_text(report)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', required=True, type=Path)
    args = ap.parse_args()
    write_report(args.output.expanduser().resolve())
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
