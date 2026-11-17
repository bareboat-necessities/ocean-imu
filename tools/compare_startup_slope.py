#!/usr/bin/env python3
"""Compare fixed and motion-scaled force slope on identical startup windows."""
import argparse
from collections import defaultdict
import csv
import gzip
import json
import math
from pathlib import Path

from sim_dataset import digest


def read_run(folder, model):
    receipt = json.loads((folder / 'provenance.json').read_text())
    if receipt.get('force_slope_model', 'fixed') != model:
        raise ValueError(f'expected {model} weighting in {folder}')
    path = folder / receipt['metrics_file']
    if digest(path) != receipt['metrics_sha256']:
        raise ValueError(f'metrics digest mismatch in {folder}')
    methods = {'calibrated': 'joint', 'residual': 'joint_cal'}
    trial_rows = {}
    with gzip.open(path, 'rt') as stream:
        for r in csv.DictReader(stream):
            if r['method'] not in ('shipping', methods[r['profile']]):
                continue
            if r['schedule'] == 'fixed' and float(r['elapsed_s']) != 35:
                continue
            key = tuple(r[k] for k in ('case', 'input', 'profile', 'schedule', 'seed', 'offset', 'epoch_s'))
            if key in trial_rows or r['valid'] != '1':
                raise ValueError(f'duplicate or invalid endpoint {key}')
            trial_rows[key] = tuple(float(r[k]) for k in ('yaw_error_deg', 'tilt_error_deg', 'bg_error_dps'))
    return receipt, trial_rows


def compare(fixed, scaled, output):
    first, baseline = read_run(fixed, 'fixed')
    second, candidate = read_run(scaled, 'scaled')
    for key in ('datasets', 'seeds', 'offsets_s', 'profiles', 'imu_hz', 'mag_hz', 'deadline_s'):
        if first[key] != second[key]:
            raise ValueError(f'unmatched input protocol: {key}')
    if not baseline or baseline.keys() != candidate.keys():
        raise ValueError('entered different startup windows or scored different epochs')
    groups = defaultdict(list)
    for key, prior in baseline.items():
        current = candidate[key]
        if key[3] == 'handoff' and prior != current:
            raise ValueError(f'shipping control changed at {key}')
        groups[key[:4]].append(prior + current)
    with output.open('w', newline='') as f:
        writer = csv.writer(f, lineterminator='\n')
        writer.writerow(['case', 'input', 'profile', 'schedule', 'n',
                         *(f'{model}_{metric}_rms' for model in ('fixed', 'scaled')
                           for metric in ('yaw_deg', 'tilt_deg', 'bg_dps'))])
        for key, values in sorted(groups.items()):
            writer.writerow([*key, len(values),
                             *(f'{math.sqrt(sum(x*x for x in column) / len(values)):.8f}'
                               for column in zip(*values))])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('fixed', 'scaled', 'output'):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    compare(args.fixed, args.scaled, args.output)
