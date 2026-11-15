#!/usr/bin/env python3
"""Run both pinned directions, preserving every score and failing CI on any gate.

The --measure option is for reviewing a new dataset. It never changes limits.
Charts use full-rate metrics, not decimated samples. Dense per-filter traces and
manuscript charts are generated separately by the build.yml case matrix.
"""
from __future__ import annotations
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import sys

from sim_dataset import CASES, RELEASE, ROOT
from filter_sim_gates import REPLAYS


def metric_rows(report):
    rows = []
    for family, profiles in report['records'].items():
        for profile, records in profiles.items():
            for record in records:
                rows.append(dict(case=report['case'], family=family, profile=profile, **record))
    other = report['other_replays']
    for record in other['pii']:
        rows.append(dict(case=report['case'], family='pii', profile='default', **record))
    for record in other['ahrs']:
        rows.append(dict(case=report['case'], family='ahrs', profile='default', input=record['input'],
                         **{a+'_rms_deg': record[a] for a in ('roll', 'pitch', 'yaw')}))
    names = sorted(report['inputs']['input_sha256'])
    for name, value in zip(names, other['nlo_z_percent_hs']):
        rows.append(dict(case=report['case'], family='nlo', profile='default', input=name, disp_z_pct_hs=value))
    for profile, records in other['heave'].items():
        for record in records:
            rows.append(dict(case=report['case'], family='heave_'+record['method'], profile=profile,
                             input=Path(record['record']).name, disp_z_pct_hs=float(record['z_pct_hs'])))
    return rows


def publish(report, destination):
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams['svg.hashsalt'] = 'ocean-imu-sim-v1.4.0'
    import matplotlib.pyplot as plt
    import numpy as np
    rows = metric_rows(report)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with (destination/'metrics.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    records = sorted(report['inputs']['input_sha256'])
    for metric, title, unit in [('disp_z_pct_hs', 'Heave RMS', '% incident Hs'),
                                ('disp_3d_pct_refmax', '3-D displacement RMS', '% reference maximum'),
                                ('roll_rms_deg', 'Roll RMS', 'deg'),
                                ('pitch_rms_deg', 'Pitch RMS', 'deg'),
                                ('yaw_rms_deg', 'Yaw RMS', 'deg'),
                                ('accel_bias_3d_rms_mps2', 'Accelerometer bias RMS', 'm/s²'),
                                ('gyro_bias_3d_rms_radps', 'Gyroscope bias RMS', 'rad/s')]:
        chosen = [r for r in rows if r['profile'] == 'default' and metric in r]
        families = list(dict.fromkeys(r['family'] for r in chosen))
        values = np.array([[next(r[metric] for r in chosen if r['family']==f and r['input']==n)
                            for n in records] for f in families])
        fig, ax = plt.subplots(figsize=(11, max(4, .35*len(families)+1.5)))
        im = ax.imshow(values, aspect='auto', cmap='YlOrRd')
        ax.set_xticks(range(8), [n.split('_')[2]+' '+n.split('_')[3][1:]+' m' for n in records], rotation=30, ha='right')
        ax.set_yticks(range(len(families)), families)
        ax.set_title(f"{title} · {RELEASE} · {CASES[report['case']]['description']}", loc='left', pad=16)
        for i in range(len(families)):
            for j in range(8):
                ax.text(j, i, f'{values[i,j]:.3g}', ha='center', va='center', fontsize=8,
                        color='white' if values[i,j] > values.max()*.65 else '#202020')
        fig.colorbar(im, ax=ax, label=unit)
        fig.tight_layout()
        fig.savefig(destination/(metric+'.svg'), metadata={'Date': None})
        plt.close(fig)
    (destination/'README.md').write_text(
        f"# {RELEASE}: {CASES[report['case']]['description']}\n\n"
        f"{report['acceptance']}. {len(rows)} filter/profile/sea cells. "
        "Eight records per profile; displacement and bias scores use the trailing 900 s, "
        "AHRS the trailing 60 s. No gate is fitted during this run.\n\n"
        "`measurements.json` contains input and executable hashes and exit codes; "
        "`metrics.csv` contains the full-rate scores. SVGs compare production defaults "
        "across all eight incident seas. Legacy and truth/proxy frontend scores remain in the CSV.\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases', nargs='+', choices=CASES, default=list(CASES))
    parser.add_argument('--output-dir', type=Path, default=ROOT/'reports/results/sim_v1.4.0')
    parser.add_argument('--jobs', type=int, default=3)
    parser.add_argument('--build', action='store_true')
    parser.add_argument('--measure', action='store_true')
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    if args.build:
        for folder, binary, _ in REPLAYS.values():
            subprocess.run(['make', '-C', str(ROOT/'tests'/folder), binary], check=True)
    failed = False
    for case in args.cases:
        env = {**os.environ, 'OCEAN_IMU_SIM_CASE': case}
        data = ROOT/'.sim-data/cases'/case
        subprocess.run([sys.executable, str(ROOT/'tools/sim_dataset.py'), '--case', case,
                        '--dest', str(data)], env=env, check=True)
        out = args.output_dir.resolve()/case
        out.mkdir(parents=True, exist_ok=True)
        # Readouts are always kept, including numeric failures. The final exit
        # is a gate decision; --measure is never used by normal CI.
        command = [sys.executable, str(ROOT/'tools/filter_sim_gates.py'), '--run',
                   '--data-dir', str(data), '--log-dir', str(out/'logs'), '--jobs', str(args.jobs),
                   '--output', str(out/'measurements.json')]
        if args.measure:
            command.append('--measure')
        (out/'measurements.json').unlink(missing_ok=True)
        done = subprocess.run(command, env=env, check=False)
        failed |= done.returncode != 0
        if (out/'measurements.json').exists():
            report = json.loads((out/'measurements.json').read_text())
            publish(report, out)
    return int(failed)


if __name__ == '__main__':
    raise SystemExit(main())
