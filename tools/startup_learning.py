#!/usr/bin/env python3
"""Paired, isolated startup learning on both pinned underway datasets."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import platform
import shlex
import subprocess

import sim_dataset as dataset

ROOT = Path(__file__).resolve().parents[1]
TEST = ROOT / 'tests/startup_learning'
METRICS = ('yaw_error_deg', 'tilt_error_deg', 'attitude_error_deg',
           'field_direction_error_deg', 'field_magnitude_error_uT', 'field_error_uT',
           'bg_error_dps', 'ba_error_mps2', 'rate_error_dps')
METHODS = ('uniform', 'tapered', 'magfit', 'joint', 'joint_cal')


def summarize(rows):
    groups = {}
    for row in rows:
        horizon = row['elapsed_s'] if row['schedule'] == 'fixed' else 'handoff'
        key = (row['case'], row['profile'], row['method'], row['schedule'], horizon)
        groups.setdefault(key, []).append(row)
    summary = []
    for key, cell in sorted(groups.items()):
        item = dict(zip(('case', 'profile', 'method', 'schedule', 'horizon'), key))
        item['n'] = len(cell)
        for metric in METRICS:
            values = [float(r[metric]) for r in cell]
            if not all(math.isfinite(v) for v in values):
                raise ValueError(f'non-finite metric {key}: {metric}')
            item[metric.replace('_error', '_rms')] = math.sqrt(sum(v*v for v in values)/len(values))
        item['ready_fraction'] = sum(int(r['ready']) for r in cell)/len(cell)
        item['mean_elapsed_s'] = sum(float(r['elapsed_s']) for r in cell)/len(cell)
        # Explicit accuracy target for false-readiness diagnostics, not a
        # claim that a model-residual gate certifies these truth errors.
        ready = [r for r in cell if int(r['ready'])]
        item['false_ready_fraction'] = (sum(abs(float(r['yaw_error_deg'])) > 2 or
                                              float(r['tilt_error_deg']) > 1 for r in ready)/len(ready)
                                        if ready else None)
        item['worst_yaw_deg'] = max(abs(float(r['yaw_error_deg'])) for r in cell)
        item['worst_tilt_deg'] = max(float(r['tilt_error_deg']) for r in cell)
        summary.append(item)
    return summary


def validate_completeness(rows, seeds, offsets, profiles):
    """Missing baseline handovers are deadline misses, never dropped success."""
    cells = {}
    for row in rows:
        key = (row['case'], row['input'], row['profile'], int(row['seed']), float(row['offset']))
        cells.setdefault(key, []).append(row)
    expected = {(case, name, profile, int(seed), float(offset))
                for case in dataset.CASES for name in dataset.REFERENCE_CSV_SHA256_BY_CASE[case]
                for profile in profiles.split(',') for seed in seeds.split(',') for offset in offsets.split(',')}
    if set(cells) != expected:
        raise ValueError('incomplete experiment cells')
    coverage = []
    for key, cell in sorted(cells.items()):
        fixed = [(r['method'], float(r['elapsed_s'])) for r in cell if r['schedule'] == 'fixed']
        wanted = {(m, float(t)) for m in METHODS for t in (10, 15, 20, 25, 30, 35)}
        if len(fixed) != len(wanted) or set(fixed) != wanted or any(r['valid'] != '1' for r in cell):
            raise ValueError(f'missing/duplicate/invalid fixed horizon in {key}')
        baseline = [r for r in cell if r['method'] == 'shipping']
        matched = [r for r in cell if r['schedule'] == 'matched']
        if len(baseline) > 1 or len(matched) != len(METHODS)*len(baseline):
            raise ValueError(f'incomplete paired handover in {key}')
        if baseline and (set(r['method'] for r in matched) != set(METHODS) or any(
                r['elapsed_s'] != baseline[0]['elapsed_s'] for r in matched)):
            raise ValueError(f'physical handover epoch mismatch in {key}')
        item = dict(zip(('case', 'input', 'profile', 'seed', 'offset_s'), key))
        item['shipping_handover_s'] = float(baseline[0]['elapsed_s']) if baseline else None
        item['shipping_deadline_missed'] = not baseline
        for method in METHODS:
            history = sorted((r for r in cell if r['method'] == method and r['schedule'] == 'fixed'),
                             key=lambda r: float(r['elapsed_s']))
            ready = [r for r in history if r['ready'] == '1']
            item[f'{method}_first_ready_s'] = float(ready[0]['elapsed_s']) if ready else None
            item[f'{method}_false_first_ready'] = bool(ready and (
                abs(float(ready[0]['yaw_error_deg'])) > 2 or float(ready[0]['tilt_error_deg']) > 1))
            # Truth-only diagnostic: first checkpoint meeting the target and
            # remaining within it through all later checkpoints, up to 35 s.
            accurate = [abs(float(r['yaw_error_deg'])) <= 2 and float(r['tilt_error_deg']) <= 1 for r in history]
            item[f'{method}_accuracy_by_s'] = next((float(history[i]['elapsed_s']) for i in range(len(history))
                                                    if all(accurate[i:])), None)
        coverage.append(item)
    return coverage


def check_accuracy(summary):
    """Directional, paired regression bars; no thresholds use scored truth at runtime."""
    checks = []
    for case in dataset.CASES:
        for profile, method, yaw_limit in (('calibrated', 'joint', 2.0), ('residual', 'joint_cal', 2.5)):
            def find(name, schedule, case=case, profile=profile):
                return next((s for s in summary if s['case'] == case and s['profile'] == profile and
                             s['method'] == name and s['schedule'] == schedule), None)
            baseline, candidate = find('shipping', 'handoff'), find(method, 'matched')
            if baseline is None or candidate is None:
                raise ValueError(f'missing paired accuracy cell {case}/{profile}')
            for metric, limit in (('yaw_rms_deg', min(yaw_limit, .85*baseline['yaw_rms_deg'])),
                                  ('tilt_rms_deg', baseline['tilt_rms_deg']),
                                  ('bg_rms_dps', baseline['bg_rms_dps'])):
                actual = candidate[metric]
                checks.append({'case': case, 'profile': profile, 'method': method,
                               'schedule': 'matched', 'metric': metric,
                               'actual': actual, 'limit': limit, 'passed': actual <= limit})
            # Paired comparisons omit shipping deadline misses. Gate every
            # startup again at 35 s so those omissions cannot hide regressions.
            endpoint = next((s for s in summary if s['case'] == case and s['profile'] == profile and
                             s['method'] == method and s['schedule'] == 'fixed' and
                             s['horizon'] == '35'), None)
            if endpoint is None:
                raise ValueError(f'missing 35-second accuracy cell {case}/{profile}')
            # The calibrated arm must meet the declared target in EVERY
            # trial. Residual-calibration yaw is still only an aggregate bar;
            # retain its individual target failures in the coverage report.
            yaw_metric = 'worst_yaw_deg' if profile == 'calibrated' else 'yaw_rms_deg'
            for metric, limit in ((yaw_metric, yaw_limit), ('worst_tilt_deg', 1.0), ('bg_rms_dps', .04)):
                actual = endpoint[metric]
                checks.append({'case': case, 'profile': profile, 'method': method,
                               'schedule': 'fixed', 'horizon_s': 35, 'metric': metric,
                               'actual': actual, 'limit': limit, 'passed': actual <= limit})
    return checks


def check_controls(summary):
    checks = []
    for case in dataset.CASES:
        for method in ('joint', 'joint_cal'):
            for profile in ('ideal', 'interference'):
                cell = next((s for s in summary if s['case'] == case and s['method'] == method and
                             s['profile'] == profile and s['schedule'] == 'fixed' and s['horizon'] == '35'), None)
                if cell is None:
                    raise ValueError(f'missing control cell {case}/{method}/{profile}')
                metrics = (('yaw_rms_deg', .25), ('tilt_rms_deg', .1)) if profile == 'ideal' else (
                    ('ready_fraction', 0.0),)
                for metric, limit in metrics:
                    checks.append({'case': case, 'method': method, 'profile': profile,
                                   'metric': metric, 'actual': cell[metric], 'limit': limit,
                                   'passed': cell[metric] <= limit})
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seeds', default='1')
    parser.add_argument('--offsets', default='0')
    parser.add_argument('--profiles', default='ideal,calibrated,residual')
    parser.add_argument('--jobs', type=int, default=2)
    parser.add_argument('--output', type=Path, default=ROOT/'reports/results/startup_learning')
    parser.add_argument('--scale-force-slope', action='store_true',
                        help='weight the joint force slope by transported force fluctuations')
    parser.add_argument('--check', action='store_true', help='require paired calibrated and residual accuracy bars')
    parser.add_argument('--check-controls', action='store_true', help='require ideal accuracy and field-step rejection')
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    subprocess.run(['make', '-C', str(TEST), '-j2'], check=True)
    args.output.mkdir(parents=True, exist_ok=True)
    paths = {}
    provenance = {}
    for case, spec in dataset.CASES.items():
        folder = TEST/'data'/case
        paths[case] = [folder/name for name in dataset.REFERENCE_CSV_SHA256_BY_CASE[case]]
        if not all(p.exists() for p in paths[case]):
            archive = ROOT/'.sim-data'/dataset.RELEASE/spec['archive']
            dataset.fetch(archive, case)
            dataset.materialize(archive, [folder], case)
        provenance[case] = dataset.input_provenance(paths[case], case)
    jobs = [(case, path) for case in paths for path in paths[case]]

    def run(job):
        case, path = job
        result = args.output/f'{case}-{path.stem}.csv'
        subprocess.run([str(TEST/'startup-learning'), '--input', str(path),
                        '--output', str(result), '--seeds', args.seeds,
                        '--offsets', args.offsets, '--profiles', args.profiles,
                        *(['--scale-force-slope'] if args.scale_force_slope else [])], check=True)
        with result.open() as stream:
            rows = [dict(case=case, **r) for r in csv.DictReader(stream)]
        result.unlink()
        return rows

    rows = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for chunk in pool.map(run, jobs):
            rows.extend(chunk)
    coverage = validate_completeness(rows, args.seeds, args.offsets, args.profiles)
    # Six checkpoints per startup (no full-rate traces). Deterministic gzip
    # keeps the reviewable evidence small even with all handover components.
    metrics_path = args.output/'metrics.csv.gz'
    with metrics_path.open('wb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding='utf-8', newline='') as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
    summary = summarize(rows)
    (args.output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    (args.output/'coverage.json').write_text(json.dumps(coverage, indent=2, allow_nan=False)+'\n')
    sources = [*TEST.glob('*.h'), *TEST.glob('*.cpp'), TEST/'Makefile', Path(__file__)]
    receipt = {
        'base_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'sources_sha256': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                           for p in sorted(sources)},
        'datasets': provenance, 'seeds': args.seeds, 'offsets_s': args.offsets,
        'force_slope_model': 'scaled' if args.scale_force_slope else 'fixed',
        'profiles': args.profiles, 'imu_hz': 200, 'mag_hz': 25, 'deadline_s': 35,
        'truth_access': 'Only sensor generation and scoring; candidate takes Sample only.',
        'metrics_file': metrics_path.name, 'metrics_sha256': dataset.digest(metrics_path),
        'executable_sha256': dataset.digest(TEST/'startup-learning'),
        'compiler': subprocess.check_output([*shlex.split(os.environ.get('CXX', 'c++')), '--version'],
                                            text=True).splitlines()[0],
        'python': platform.python_version(),
    }
    (args.output/'provenance.json').write_text(json.dumps(receipt, indent=2)+'\n')
    if args.check or args.check_controls:
        checks = (check_accuracy(summary) if args.check else []) + (check_controls(summary) if args.check_controls else [])
        (args.output/'checks.json').write_text(json.dumps(checks, indent=2)+'\n')
        if not all(c['passed'] for c in checks):
            raise ValueError(f'accuracy regression: {[c for c in checks if not c["passed"]]}')
    for item in summary:
        if item['method'] in ('shipping', 'joint') and (
                item['schedule'] != 'fixed' or item['horizon'] == '35'):
            print(f"{item['case']:9} {item['profile']:12} {item['method']:8} "
                  f"{item['schedule']:7} n={item['n']:3} t={item['mean_elapsed_s']:.2f} "
                  f"yaw={item['yaw_rms_deg']:.4f} tilt={item['tilt_rms_deg']:.4f} "
                  f"field={item['field_direction_rms_deg']:.4f} bg={item['bg_rms_dps']:.5f}")


if __name__ == '__main__':
    main()
