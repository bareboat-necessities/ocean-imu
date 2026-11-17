#!/usr/bin/env python3
"""Render compact tables and plots from the isolated startup experiment."""
from __future__ import annotations
import argparse
import csv
import gzip
import json
from pathlib import Path
import statistics
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('evidence', type=Path)
    args = parser.parse_args()
    root = args.evidence
    summary = json.loads((root/'summary.json').read_text())
    coverage = json.loads((root/'coverage.json').read_text())
    provenance = json.loads((root/'provenance.json').read_text())
    with gzip.open(root/'metrics.csv.gz', 'rt') as stream:
        rows = list(csv.DictReader(stream))
    methods = ('shipping', 'uniform', 'tapered', 'magfit', 'joint', 'joint_cal')
    text = ['# Startup learning: measured results', '',
            f"Base: `{provenance['base_commit']}`. Seeds: {provenance['seeds']}; "
            f"power-on offsets: {provenance['offsets_s']} s. Both pinned v1.4.0 datasets.", '',
            'These are held-out sensor seeds and startup phases on the existing wave records, '
            'not independent newly generated seas. Values below pool one endpoint error per trial; '
            'they are not live-filter RMS over a navigation run.', '',
            '## Equal-epoch handover comparison', '',
            'Every candidate row in this table uses exactly the shipping handover timestamp. '
            'A shipping startup missing the 35 s deadline is counted separately below.', '',
            '| Direction | Sensor profile | Method | N | Yaw RMS ° | Tilt RMS ° | Field direction RMS ° | '
            'Field magnitude RMS µT | Gyro bias RMS °/s |',
            '|---|---|---|---:|---:|---:|---:|---:|---:|']
    for case in ('following', 'head'):
        for profile in sorted({s['profile'] for s in summary}):
            for method in methods:
                cell = [s for s in summary if s['case'] == case and s['profile'] == profile and
                        s['method'] == method and s['schedule'] == ('handoff' if method == 'shipping' else 'matched')]
                if cell:
                    s = cell[0]
                    text.append(f"| {case} | {profile} | {method} | {s['n']} | {s['yaw_rms_deg']:.3f} | "
                                f"{s['tilt_rms_deg']:.3f} | {s['field_direction_rms_deg']:.3f} | "
                                f"{s['field_magnitude_rms_uT']:.3f} | {s['bg_rms_dps']:.4f} |")
    text += ['', '## Fixed 35-second endpoint', '',
             '| Direction | Sensor profile | Method | N | Yaw RMS / worst ° | Tilt RMS / worst ° | Gyro bias RMS °/s | '
             'Ready % | False ready % of ready |',
             '|---|---|---|---:|---:|---:|---:|---:|---:|']
    for s in summary:
        if s['schedule'] == 'fixed' and s['horizon'] == '35' and s['method'] in ('joint', 'joint_cal'):
            false = s['false_ready_fraction']
            false_text = f'{100*false:.1f}' if false is not None else '—'
            text.append(f"| {s['case']} | {s['profile']} | {s['method']} | {s['n']} | "
                        f"{s['yaw_rms_deg']:.3f} / {s['worst_yaw_deg']:.3f} | "
                        f"{s['tilt_rms_deg']:.3f} / {s['worst_tilt_deg']:.3f} | {s['bg_rms_dps']:.4f} | "
                        f"{100*s['ready_fraction']:.1f} | {false_text} |")
    text += ['', 'The diagnostic readiness target is yaw ≤2° and tilt ≤1°. '
             'Readiness uses sensor statistics; target failures are evaluated using truth afterward. '
             'Neither the flag nor a small residual certifies accuracy.', '',
             '## Deadline and latency', '',
             '| Direction | Sensor profile | Method | Trials | Shipping misses at 35 s | Shipping median s | '
             'Accurate by 35 s | Median first ready s | False first readiness |',
             '|---|---|---|---:|---:|---:|---:|---:|---:|']
    for case in ('following', 'head'):
        for profile in sorted({s['profile'] for s in summary}):
            cell = [c for c in coverage if c['case'] == case and c['profile'] == profile]
            times = [c['shipping_handover_s'] for c in cell if c['shipping_handover_s'] is not None]
            if not cell:
                continue
            median = f'{statistics.median(times):.2f}' if times else '—'
            n = len(cell)
            for method in ('joint', 'joint_cal'):
                ready = [c[f'{method}_first_ready_s'] for c in cell if c[f'{method}_first_ready_s'] is not None]
                first = f'{statistics.median(ready):.1f}' if ready else '—'
                text.append(f"| {case} | {profile} | {method} | {n} | {sum(c['shipping_deadline_missed'] for c in cell)} | "
                            f"{median} | {sum(c[f'{method}_accuracy_by_s'] is not None for c in cell)}/{n} | "
                            f"{first} | {sum(c[f'{method}_false_first_ready'] for c in cell)}/{len(ready)} |")
    solve_times = [float(r['solve_ms']) for r in rows if r['method'] == 'joint' and
                   r['schedule'] == 'fixed' and float(r['elapsed_s']) == 35]
    if solve_times:
        text += ['', f'Host joint solve at 35 s: median {statistics.median(solve_times):.1f} ms; '
                 f'maximum {max(solve_times):.1f} ms in this run. This does not qualify ESP32 execution time or memory.']
    text += ['', 'Accelerometer bias is retained as a calibration prior, so its error is unchanged. '
             'The gyroscope/field covariance is conditional and cannot be substituted for a full correlated '
             'live-filter handover covariance.', '',
             'The reproducible raw metrics and actual handover components are in `metrics.csv.gz`; '
             '`coverage.json` retains every deadline miss and first-readiness result. '
             '`provenance.json` records input and source hashes.', '']
    (root/'REPORT.md').write_text('\n'.join(text))

    if not any(s['profile'] == 'calibrated' for s in summary):
        return
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams['svg.fonttype'] = 'none'
    matplotlib.rcParams['svg.hashsalt'] = 'ocean-imu-startup-learning'
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'plots'))
    from vector_export import save_figure
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.5))
    colors = {'uniform': '#9aa2ad', 'tapered': '#7986b0', 'magfit': '#d99126',
              'joint': '#187c78', 'joint_cal': '#9d5aa0'}
    metrics = [('yaw_rms_deg', 'Yaw RMS (degrees)'), ('tilt_rms_deg', 'Tilt RMS (degrees)'),
               ('bg_rms_dps', 'Gyro-bias RMS (degrees/s)')]
    for row, case in enumerate(('following', 'head')):
        for col, (metric, label) in enumerate(metrics):
            ax = axes[row, col]
            for method, color in colors.items():
                cell = sorted((s for s in summary if s['case'] == case and s['profile'] == 'calibrated' and
                               s['method'] == method and s['schedule'] == 'fixed'), key=lambda s: float(s['horizon']))
                ax.plot([float(s['horizon']) for s in cell], [s[metric] for s in cell],
                        '-o', markersize=3, color=color, label=method)
            baseline = next((s for s in summary if s['case'] == case and s['profile'] == 'calibrated' and
                             s['method'] == 'shipping'), None)
            if baseline is not None:
                ax.scatter([baseline['mean_elapsed_s']], [baseline[metric]], marker='x', s=70, linewidth=2,
                           color='#333333', label='shipping at its handover', zorder=5)
            ax.set(xlabel='Available startup data (s)', ylabel=label, title=f'{case.capitalize()} seas')
            ax.grid(alpha=.2)
            ax.set_ylim(bottom=0)
    legend = {}
    for ax in axes.flat:
        handles, labels = ax.get_legend_handles_labels()
        legend.update(zip(labels, handles))
    fig.legend(legend.values(), legend.keys(), loc='lower center', bbox_to_anchor=(.5, .01),
               ncol=3, frameon=False)
    fig.suptitle('Isolated startup: noisy / biased IMU, calibrated magnetometer', fontsize=14)
    fig.tight_layout(rect=(0, .12, 1, .94))
    svg = root/'convergence.svg'
    save_figure(fig, svg, metadata={'Date': None})
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    save_figure(fig, root/'convergence.png', dpi=150)
    plt.close(fig)


if __name__ == '__main__':
    main()
