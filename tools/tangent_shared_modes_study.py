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
import os
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
        if a.resume and receipt_path.exists():
            receipt = json.loads(receipt_path.read_text())
            if (receipt['binary_sha256'] != binary_hash or receipt['command'] != command
                    or receipt['scores_sha256'] != digest(str(prefix)+'_scores.csv')
                    or receipt['trace_sha256'] != digest(str(prefix)+'_trace.csv')):
                raise RuntimeError(f'stale resume receipt {receipt_path}')
        else:
            begin = time.time()
            result = subprocess.run(command, capture_output=True, text=True)
            receipt = dict(command=command, binary_sha256=binary_hash, exit_code=result.returncode,
                           seconds=time.time()-begin, stdout=result.stdout, stderr=result.stderr,
                           scores_sha256=digest(str(prefix)+'_scores.csv'),
                           trace_sha256=digest(str(prefix)+'_trace.csv'),
                           model_sha256=digest(str(prefix)+'_model.json'))
            receipt_path.write_text(json.dumps(receipt, indent=2)+'\n')
        with Path(str(prefix)+'_scores.csv').open() as f: rows = list(csv.DictReader(f))
        if len({r['sensor_hash'] for r in rows}) != 1:
            raise RuntimeError(f'unpaired sensors: {prefix}')
        for r in rows: r['noise_profile'] = model
        print(f'{case} {seed} {model}: exit={receipt["exit_code"]}', flush=True)
        return rows, receipt

    rows, receipts = [], []
    with ThreadPoolExecutor(max_workers=max(1, a.workers)) as pool:
        futures = [pool.submit(run, j) for j in jobs]
        for future in as_completed(futures):
            rr, receipt = future.result(); rows.extend(rr); receipts.append(receipt)
    rows.sort(key=lambda r:(r['scenario'], int(r['seed']), r['noise_profile'], r['filter'], float(r['start']), float(r['end'])))
    with (a.output/'scores.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    if {str(p.relative_to(ROOT)):digest(p) for p in source_paths} != source_hashes:
        raise RuntimeError('study sources changed while the cohort was running')
    manifest = dict(baseline_commit=BASELINE, implementation_commit=implementation_commit, profile=a.profile, seeds=seeds, cases=cases,
                    binary_sha256=binary_hash, source_sha256=source_hashes,
                    commands=sorted(receipts, key=lambda r:r['command']), seconds=time.time()-started,
                    python=sys.version, platform=platform.platform(),
                    compiler=subprocess.check_output(['g++', '--version'], text=True).splitlines()[0],
                    score_rate_hz=20, trace_rate_hz=1,
                    interpretation='TG(SE3) tangent-bias hybrid, not full TG(SE2(3)); no autonomous CARMA error theorem',
                    runtime='per-filter process CPU time; host proxy, not ESP32 board measurement',
                    thresholds='unchanged TFG native replay bars; experimental smoke only is in make all')
    (a.output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(f'Completed {len(jobs)} paired sensor histories, {len(rows)} score rows')


if __name__ == '__main__': main()
