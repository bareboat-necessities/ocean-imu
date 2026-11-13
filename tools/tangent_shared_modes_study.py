#!/usr/bin/env python3
"""Reproduce paired tangent-bias research without changing production tuning.

Raw one-Hz traces stay outside git. Committed scores are computed online from
20-Hz estimates, not plotted/decimated traces. All seven variants in a row
consume one sensor history; seeds, noise, prior and evaluation windows are
explicit. No RMS-dependent tuning or promotion is performed by this script.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
BASELINE = '7712a0e046bb653839a4d582ec893d7671ac1e5f'
CASES = ['stationary', 'rocking', 'wave', 'low', 'high', 'transition', 'yaw',
         'inclination0', 'inclination87', 'reference_error', 'exact_reference',
         'initial_bias', 'initial_attitude', 'late_bias', 'late_attitude',
         'gyro_error', 'bias_drift', 'gyro_drift', 'vibration', 'colored', 'long', 'rate100', 'rate400',
         'mag100', 'mag10', 'jitter']


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def optional_digest(path):
    try:
        return digest(path)
    except OSError:
        return None


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=Path, default=ROOT/'tests/kalman_tangent_shared_modes/tangent-study')
    p.add_argument('--raw', type=Path, required=True)
    p.add_argument('--output', type=Path, default=ROOT/'reports/results/tangent_shared_modes')
    p.add_argument('--profile', choices=['development', 'acceptance'], default='acceptance')
    p.add_argument('--workers', type=int, default=2)
    p.add_argument('--cases', nargs='+', choices=CASES)
    p.add_argument('--seeds', nargs='+', type=int)
    p.add_argument('--resume', action='store_true')
    a = p.parse_args()
    a.raw.mkdir(parents=True, exist_ok=True); a.output.mkdir(parents=True, exist_ok=True)
    cases = a.cases or (['rocking', 'wave'] if a.profile == 'development' else CASES)
    seeds = a.seeds or ([7301, 7302] if a.profile == 'development' else list(range(4101, 4107)))
    jobs = [(case, seed, 'native') for case in cases for seed in seeds]
    if a.profile == 'acceptance' and 'wave' in cases:
        jobs += [('wave', seed, model) for seed in seeds for model in ('noise', 'population')]
    source_paths = [ROOT/'tools/tangent_shared_modes_study.py', ROOT/'src/lie/TangentBiasGroup.h',
                    ROOT/'src/lie/SE3Jacobians.h', ROOT/'src/kalman_common/BiasLearningGate.h']
    source_paths += sorted((ROOT/'tests/kalman_tangent_shared_modes').glob('*.cpp'))
    source_paths += sorted((ROOT/'tests/kalman_tangent_shared_modes').glob('*.h'))
    source_paths += sorted((ROOT/'src/kalman_carma_shared_modes').glob('*.h'))
    source_paths += [ROOT/'src/kalman_carma/Kalman3D_Wave_CARMA3.h', ROOT/'src/kalman_carma/Kalman3D_Wave_TFG_CARMA3.h']
    source_paths += [p for p in (ROOT/'src').rglob('*') if p.suffix in ('.h', '.cpp')]
    source_paths += [ROOT/'tests/kalman_carma/CarmaBenchTruth.h', ROOT/'tests/kalman_carma/carma3/Bench3.h', ROOT/'tests/common/CarmaModelExport.h', ROOT/'tests/kalman_tangent_shared_modes/Makefile']
    source_hashes = {str(p.relative_to(ROOT)):digest(p) for p in source_paths}
    binary_hash = digest(a.binary)
    implementation_commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    started = time.time()

    def run(job):
        case, seed, model = job
        prefix = a.raw/f'{case}_{seed}_{model}'
        duration = 180 if a.profile == 'development' else 7200 if case == 'long' else 900
        command = [str(a.binary.resolve()), '--scenario', case, '--seed', str(seed),
                   '--duration', str(duration), '--output', str(prefix)]
        if model in ('noise', 'population'): command += ['--matched-noise']
        if model == 'population': command += ['--population-prior']
        if case == 'wave': command += ['--filters', 'QSM,TFG,TAN,QSM-G,TFG-G,TAN-G,CARMA,QSM-LG,TFG-LG,TAN-LG']
        receipt_path = Path(str(prefix)+'_receipt.json')
        artifacts = {'scores_sha256': str(prefix)+'_scores.csv',
                     'trace_sha256': str(prefix)+'_trace.csv',
                     'model_sha256': str(prefix)+'_model.json'}
        if a.resume and receipt_path.exists():
            receipt = json.loads(receipt_path.read_text())
            if (receipt['binary_sha256'] != binary_hash or receipt['command'] != command
                    or any(receipt.get(key) != optional_digest(path) for key, path in artifacts.items())):
                raise RuntimeError(f'stale resume receipt {receipt_path}')
        else:
            begin = time.time()
            receipt = dict(command=command, binary_sha256=binary_hash)
            try:
                result = subprocess.run(command, capture_output=True, text=True)
                receipt.update(exit_code=result.returncode, stdout=result.stdout, stderr=result.stderr)
            except OSError as error:
                receipt.update(exit_code=None, stdout='', stderr=str(error),
                               launch_error=type(error).__name__)
            receipt['seconds'] = time.time()-begin
            # Preserve termination diagnostics before inspecting any output:
            # a crashed process may not have opened even its score file.
            receipt_path.write_text(json.dumps(receipt, indent=2)+'\n')
        receipt.update({key:optional_digest(path) for key, path in artifacts.items()})
        errors = []
        for key, path in artifacts.items():
            if receipt[key] is None:
                errors.append('unavailable output: '+Path(path).name)
        rows = []
        try:
            with Path(artifacts['scores_sha256']).open() as f:
                rows = list(csv.DictReader(f))
            required = {'scenario', 'seed', 'filter', 'start', 'end', 'sensor_hash', 'failed'}
            if not rows or any(not required.issubset(row) for row in rows):
                raise ValueError('empty or incomplete score schema')
            for row in rows:
                if row['scenario'] != case or int(row['seed']) != seed:
                    raise ValueError('score condition/seed differs from command')
                float(row['start']); float(row['end']); int(row['failed'])
        except (OSError, ValueError, TypeError, csv.Error) as error:
            errors.append('unreadable scores: '+str(error))
            rows = []
        if rows and len({r['sensor_hash'] for r in rows}) != 1:
            errors.append('unpaired sensors')
        for r in rows: r['noise_profile'] = model
        receipt['errors'] = errors
        receipt['outputs_complete'] = not errors
        receipt['failed'] = (receipt['exit_code'] != 0 or bool(errors)
                             or any(int(row['failed']) != 0 for row in rows))
        receipt_path.write_text(json.dumps(receipt, indent=2)+'\n')
        print(f'{case} {seed} {model}: exit={receipt["exit_code"]}, failed={receipt["failed"]}', flush=True)
        return rows, receipt

    rows, receipts = [], []
    with ThreadPoolExecutor(max_workers=max(1, a.workers)) as pool:
        futures = [pool.submit(run, j) for j in jobs]
        for future in as_completed(futures):
            rr, receipt = future.result(); rows.extend(rr); receipts.append(receipt)
    rows.sort(key=lambda r:(r['scenario'], int(r['seed']), r['noise_profile'], r['filter'], float(r['start']), float(r['end'])))
    with (a.output/'scores.csv').open('w') as f:
        fields = list(dict.fromkeys(key for row in rows for key in row))
        writer = csv.DictWriter(f, fieldnames=fields or ['scenario', 'seed', 'filter', 'failed', 'noise_profile'],
                                lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    after = {str(p.relative_to(ROOT)):optional_digest(p) for p in source_paths}
    source_changed = after != source_hashes
    failed_runs = sum(r['failed'] for r in receipts)
    manifest = dict(baseline_commit=BASELINE, implementation_commit=implementation_commit, profile=a.profile, seeds=seeds, cases=cases,
                    binary_sha256=binary_hash, source_sha256=source_hashes,
                    commands=sorted(receipts, key=lambda r:r['command']), seconds=time.time()-started,
                    failed_experiments=failed_runs,
                    incomplete_experiments=sum(not r['outputs_complete'] for r in receipts),
                    source_changed_during_run=source_changed,
                    python=sys.version, platform=platform.platform(),
                    compiler=subprocess.check_output(['g++', '--version'], text=True).splitlines()[0],
                    score_rate_hz=20, trace_rate_hz=1,
                    interpretation='TG(SE3) tangent-bias hybrid, not full TG(SE2(3)); no autonomous CARMA error theorem',
                    runtime='per-filter process CPU time; host proxy, not ESP32 board measurement',
                    thresholds='unchanged historical and current TFG replay bars; native replay and paired smoke are in make all')
    if source_changed:
        manifest['source_sha256_after_run'] = after
    (a.output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(f'Recorded {len(jobs)} paired experiments, {len(rows)} score rows; failed={failed_runs}')
    return 2 if failed_runs or source_changed else 0


if __name__ == '__main__': raise SystemExit(main())
