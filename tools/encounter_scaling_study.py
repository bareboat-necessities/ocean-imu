#!/usr/bin/env python3
"""Pair exact-clock zero-speed/5 kn replays of the PII phase correction.

Build tests/pii_observer first. Both arms use the same executable and sensor
seeds; PII_PHASE_MATCH=0 reproduces the previous schedule. Statistics always
use the full-rate final 900 s, never the decimated drawing files.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
import zipfile

import ou_evidence_provenance as provenance
from sim_dataset import (ROOT, REFERENCE_CSV_SHA256, digest, input_provenance,
                         verify_csv_clock)
from stationary_vessel_control import ARCHIVE_SHA256 as STATIONARY_ARCHIVE_SHA256

METRICS = ('disp_z_rms_m', 'disp_z_pct_hs', 'roll_rms_deg', 'pitch_rms_deg', 'yaw_rms_deg')


def verify_stationary(archive: Path, directory: Path) -> dict:
    if digest(archive) != STATIONARY_ARCHIVE_SHA256:
        raise ValueError('zero-speed input archive must be the pinned v1.3.1 release')
    import hashlib
    records = {}
    with zipfile.ZipFile(archive) as bundle:
        for name in REFERENCE_CSV_SHA256:
            path = directory / name
            records[name] = digest(path)
            if records[name] != hashlib.sha256(bundle.read(name)).hexdigest():
                raise ValueError(f'zero-speed extracted input differs: {name}')
            verify_csv_clock(path)
    return dict(release='v1.3.1', archive='sim-data-files-vessel-rao-28ft.zip',
                archive_sha256=STATIONARY_ARCHIVE_SHA256, input_sha256=records)


def run_one(job: dict) -> dict:
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(('W3D_', 'SF_', 'OU_', 'TFG_', 'PII_', 'NLO_', 'PROBE_'))}
    env.update(W3D_COLLECT_ALL_GATES='1', W3D_WRITE_TIMESERIES='0',
               W3D_VALIDATION_WINDOW_SEC='900', PII_PHASE_MATCH=job['phase_match'])
    if job['dataset'] == 'stationary':
        env['W3D_REFERENCE_MOTION'] = 'stationary'
    if job['seed'] != 'default':
        env['W3D_SEED'] = job['seed']
    with tempfile.TemporaryDirectory(prefix='encounter-pii-') as temporary:
        work = Path(temporary)
        (work / job['input'].name).symlink_to(job['input'])
        result = subprocess.run([str(job['binary'])], cwd=work, env=env,
                                capture_output=True, text=True, timeout=600, check=False)
    lines = [line for line in result.stdout.splitlines() if line.startswith('PII_METRICS ')]
    if result.returncode not in (0, 1) or len(lines) != 1:
        raise RuntimeError(f'{job}: exit {result.returncode}\n{result.stdout[-2000:]}{result.stderr[-2000:]}')
    fields = dict(item.split('=', 1) for item in lines[0].split()[1:])
    if (fields['samples'], fields['start_s'], fields['window_s']) != ('180000', '300', '900'):
        raise ValueError('incomplete replay window')
    metrics = {key: float(fields[key]) for key in METRICS}
    if not all(math.isfinite(value) for value in metrics.values()):
        raise ValueError('nonfinite scored metric')
    row = dict(dataset=job['dataset'], input=job['input'].name, seed=job['seed'],
               arm='phase_matched' if job['phase_match'] == '1' else 'previous',
               exit_code=result.returncode, metrics=metrics)
    print(row['dataset'], row['input'], row['seed'], row['arm'], metrics['disp_z_pct_hs'], flush=True)
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stationary-archive', type=Path, required=True)
    parser.add_argument('--stationary-dir', type=Path, required=True)
    parser.add_argument('--moving-dir', type=Path, required=True)
    parser.add_argument('--seeds', default='default,7,23')
    parser.add_argument('--jobs', type=int, default=3)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    seeds = args.seeds.split(',')
    if args.jobs < 1 or len(set(seeds)) != len(seeds) or any(
            seed != 'default' and (not seed.isdecimal() or int(seed) > 2**32-1) for seed in seeds):
        parser.error('positive jobs and distinct default/unsigned integer seeds required')
    folder = ROOT / 'tests/pii_observer'
    subprocess.run(['make', '-q', '-C', str(folder), 'pii_observer-adaptive'], check=True,
                   stdout=subprocess.DEVNULL)
    binary = folder / 'pii_observer-adaptive'
    initial_binary_hash = digest(binary)
    directories = dict(stationary=args.stationary_dir.resolve(), moving=args.moving_dir.resolve())
    inputs = dict(stationary=verify_stationary(args.stationary_archive, directories['stationary']),
                  moving=input_provenance(directories['moving']/name for name in REFERENCE_CSV_SHA256))
    for name in REFERENCE_CSV_SHA256:
        verify_csv_clock(directories['moving'] / name)
    jobs = [dict(binary=binary, dataset=dataset, input=directory/name, seed=seed, phase_match=phase)
            for dataset, directory in directories.items() for name in REFERENCE_CSV_SHA256
            for seed in seeds for phase in ('0', '1')]
    with ThreadPoolExecutor(args.jobs) as pool:
        rows = list(pool.map(run_one, jobs))
    if digest(binary) != initial_binary_hash:
        raise ValueError('executable changed during the paired experiment')
    roots = [folder/'pii_observer-adaptive.cpp', ROOT/'src/util/W3dSimCommon.cpp']
    sources = set(provenance.implementation_closure(roots)) | {Path(__file__).resolve(), folder/'Makefile'}
    pairs = []
    for previous in rows:
        if previous['arm'] != 'previous':
            continue
        after = next(row for row in rows if row['arm'] == 'phase_matched' and
                     all(row[key] == previous[key] for key in ('dataset', 'input', 'seed')))
        before, current = previous['metrics']['disp_z_pct_hs'], after['metrics']['disp_z_pct_hs']
        pairs.append(dict(dataset=previous['dataset'], input=previous['input'], seed=previous['seed'],
                          previous_heave_pct_hs=before, phase_matched_heave_pct_hs=current,
                          relative_change_percent=100*(current/before-1)))
    report = dict(protocol='200 Hz, 1200 s; final 900 s/180000 samples; paired sensor/initialization seeds',
                  role='controlled paired experiment; moving-data acceptance bars are not stationary/noise-seed bars',
                  seeds=seeds, inputs=inputs, binary_sha256=initial_binary_hash,
                  git_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  source_sha256={str(p.relative_to(ROOT)): digest(p) for p in sorted(sources)},
                  build_environment=provenance.environment_metadata(), rows=rows, pairs=pairs)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
