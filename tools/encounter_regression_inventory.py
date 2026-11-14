#!/usr/bin/env python3
"""Compare every gated OU/CARMA/TFG family on identical-binary speed controls."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from sim_dataset import ROOT, REFERENCE_CSV_SHA256, SHA256, digest
from stationary_vessel_control import ARCHIVE_SHA256

FAMILIES = ('ou2', 'ou3', 'tfg', 'tfg_legacy', 'carma', 'csm', 'tcarma', 'tcsm')
METRICS = ('disp_z_pct_hs', 'disp_z_rms_m', 'disp_3d_rms_m', 'roll_rms_deg',
           'pitch_rms_deg', 'yaw_rms_deg', 'accel_bias_3d_rms_mps2',
           'gyro_bias_3d_rms_radps')


def compare(moving, stationary):
    names = set(REFERENCE_CSV_SHA256)
    if (moving['inputs']['archive_sha256'] != SHA256 or
            moving['inputs']['input_sha256'] != REFERENCE_CSV_SHA256 or
            stationary['archive_sha256'] != ARCHIVE_SHA256 or
            set(stationary['input_sha256']) != names):
        raise ValueError('comparison requires the pinned exact-clock v1.3.1 inputs')
    rows, summaries = [], {}
    for family in FAMILIES:
        control = stationary['families'][family]
        if control['binary_sha256'] != moving['binary_sha256'][family]:
            raise ValueError(f'{family}: speed comparison requires the same executable')
        by_speed = {}
        for speed, values in [('stationary', control['records']),
                              ('moving', moving['records'][family]['default'])]:
            if len(values) != len(names) or {r['input'] for r in values} != names:
                raise ValueError(f'{family}/{speed}: missing or duplicate sea records')
            by_speed[speed] = {r['input']: r for r in values}
        paired = []
        for name in sorted(names):
            metrics = {}
            for metric in METRICS:
                before, after = (float(by_speed[s][name][metric])
                                 for s in ('stationary', 'moving'))
                if not all(math.isfinite(v) and v >= 0 for v in (before, after)):
                    raise ValueError(f'{family}/{name}: invalid {metric}')
                metrics[metric] = dict(stationary=before, moving=after, delta=after-before)
            row = dict(family=family, input=name, metrics=metrics)
            rows.append(row)
            paired.append(row)
        summaries[family] = {
            metric: {speed: dict(
                mean=sum(r['metrics'][metric][speed] for r in paired)/len(paired),
                worst=max(r['metrics'][metric][speed] for r in paired))
                     for speed in ('stationary', 'moving')}
            for metric in METRICS}
        summaries[family]['heave_regressions'] = sum(
            r['metrics']['disp_z_pct_hs']['delta'] > 0 for r in paired)
    return dict(
        role='remaining speed regressions after accepted fixes, not a before/after implementation claim',
        protocol='default production sampling; eight seas per family; final 900 s at 200 Hz; same executable at both speeds',
        limits='Single default sensor seed; fixed-sea records. No all-seed or transition robustness claim. '
               'TFG_LEGACY is the legacy regularization law, not the interval-rate sampling profile. '
               'Ungated CARMA research shapes are outside this inventory.',
        binary_sha256={f: moving['binary_sha256'][f] for f in FAMILIES},
        summaries=summaries, rows=rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--moving', type=Path,
                   default=ROOT/'reports/results/filter_sim_gates/measurements.json')
    p.add_argument('--stationary', type=Path,
                   default=ROOT/'reports/results/filter_sim_gates/stationary_control.json')
    p.add_argument('--output', type=Path,
                   default=ROOT/'reports/results/encounter_scaling/remaining_regressions.json')
    args = p.parse_args()
    report = compare(json.loads(args.moving.read_text()), json.loads(args.stationary.read_text()))
    report['measurement_sha256'] = dict(moving=digest(args.moving), stationary=digest(args.stationary))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    for family, data in report['summaries'].items():
        h = data['disp_z_pct_hs']
        print(f"{family}: mean {h['stationary']['mean']:.4f} -> {h['moving']['mean']:.4f}; "
              f"worst {h['stationary']['worst']:.4f} -> {h['moving']['worst']:.4f}% Hs; "
              f"{data['heave_regressions']}/8 seas regress")


if __name__ == '__main__':
    main()
