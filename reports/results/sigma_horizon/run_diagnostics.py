#!/usr/bin/env python3
"""Measure applied-sigma variability from full simulator CSVs, then discard them."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--configs', type=Path, required=True)
    parser.add_argument('--records', default='0,4')
    parser.add_argument('--seeds', default='default')
    parser.add_argument('--jobs', type=int, default=2)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.repo / 'tools'))
    from model_mismatch_ablation import FAMILIES, RECORDS
    from ou_validation import parse_validation_metrics
    configs = json.loads(args.configs.read_text())
    records = [RECORDS[int(i)] for i in args.records.split(',')]
    allowed = {'OU_SIGMA_STILL_DECAY_SEC', 'OU_SIGMA_VAR_K_PERIODS', 'OU_SIGMA_VAR_HORIZON_MIN_S',
               'OU_SIGMA_VAR_HORIZON_MAX_S'}
    for config in configs:
        if set(config['env']) - allowed:
            raise ValueError('Only the sigma-statistics horizon may change')
    columns = ['time', 'sigma_a_applied', 'accel_var_tuner', 'tau_applied',
               'R_p0_applied', 'disp_ref_z', 'disp_est_z']
    tasks = [(family, config, record, seed) for family in ('OU-II', 'OU-III')
             for config in configs for record in records for seed in args.seeds.split(',')]

    def run(task):
        family, config, record, seed = task
        binary = FAMILIES[family]
        source = (binary.parent / record.filename).resolve()
        env = dict(os.environ, W3D_WRITE_TIMESERIES='1', W3D_COLLECT_ALL_GATES='1',
                   W3D_VALIDATION_WINDOW_SEC='900', **config['env'])
        if seed != 'default':
            env.update(W3D_IMU_SEED=seed, W3D_INIT_SEED=seed)
        with tempfile.TemporaryDirectory() as name:
            work = Path(name)
            (work / record.filename).symlink_to(source)
            result = subprocess.run([str(binary), '--input', record.filename], cwd=work,
                                    env=env, capture_output=True, text=True)
            if result.returncode not in (0, 1) or 'QUALITY_GATE: PASS=' not in result.stdout:
                raise RuntimeError(result.stdout[-2000:] + result.stderr[-2000:])
            outputs = [p for p in work.glob('*.csv') if not p.is_symlink()]
            if len(outputs) != 1:
                raise ValueError('Expected one complete output series')
            frame = pd.read_csv(outputs[0], usecols=columns)
        row = {'family': family, 'config': config['name'], 'env': config['env'],
               'input': record.filename, 'seed': seed, 'exit_code': result.returncode,
               'metrics': parse_validation_metrics(result.stdout), 'windows': {}}
        for name, mask in [('acquisition', frame.time < 300), ('scored', frame.time >= 300)]:
            part = frame.loc[mask]
            stats = {}
            for col in ('sigma_a_applied', 'accel_var_tuner', 'tau_applied', 'R_p0_applied'):
                values = part[col].to_numpy()
                stats[col] = {'mean': float(np.mean(values)), 'std': float(np.std(values)),
                              'q05': float(np.quantile(values, .05)),
                              'q95': float(np.quantile(values, .95)),
                              'increment_rms': float(np.sqrt(np.mean(np.diff(values)**2)))}
            stats['sigma_floor_fraction'] = float(np.mean(part.sigma_a_applied <= .050001))
            stats['samples'] = len(part)
            row['windows'][name] = stats
        sampled = frame.iloc[::200].copy()
        for field in ('family', 'config', 'input', 'seed'):
            sampled[field] = row[field]
        return row, sampled

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows, series = [], []
    with ThreadPoolExecutor(args.jobs) as pool:
        for row, frame in pool.map(run, tasks):
            rows.append(row)
            series.append(frame)
            print(len(rows), '/', len(tasks), row['family'], row['config'], row['input'], flush=True)
    (args.output_dir / 'diagnostics.json').write_text(json.dumps(rows, indent=2) + '\n')
    pd.concat(series).to_csv(args.output_dir / 'series-1hz.csv.gz', index=False)
    (args.output_dir / 'manifest.json').write_text(json.dumps({
        'source_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=args.repo, text=True).strip(),
        'source_diff_sha256': hashlib.sha256(subprocess.check_output(['git', 'diff', 'HEAD'], cwd=args.repo)).hexdigest(),
        'binary_sha256': {f: hashlib.sha256(FAMILIES[f].read_bytes()).hexdigest() for f in ('OU-II', 'OU-III')},
        'producer_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'configs': configs, 'seeds': args.seeds, 'records': args.records,
        'statistics_use_all_200hz_samples': True, 'saved_series_sampling_hz': 1,
    }, indent=2) + '\n')


if __name__ == '__main__':
    main()
