#!/usr/bin/env python3
"""Test whether fixed-clock covariance refresh explains encounter RMS regressions.

Use rebuilt production simulators. This is a diagnostic: all scored metrics and
numeric gate failures are retained; it does not modify acceptance limits.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
FAMILIES = {
    'ou2': ('kalman_ou_ii', 'kalman_ou_ii-sim'),
    'ou3': ('kalman_ou_iii', 'kalman_ou_iii-sim'),
    'tfg': ('kalman_tfg', 'kalman_tfg-sim'),
}
METRICS = ('disp_z_pct_hs', 'disp_z_rms_m', 'disp_3d_rms_m', 'roll_rms_deg',
           'pitch_rms_deg', 'yaw_rms_deg', 'accel_bias_3d_rms_mps2',
           'gyro_bias_3d_rms_radps')
ARMS = {
    'baseline': {},
    'no_periodic_sync': {'W3D_AW_COV_SYNC': 'reconfigure', 'TFG_AW_COV_SYNC': '0'},
    'slower_commits': {'OU_ADAPT_EVERY_SECS': '0.5'},
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, default=ROOT)
    parser.add_argument('--stationary-archive', type=Path, required=True)
    parser.add_argument('--stationary-dir', type=Path, required=True)
    parser.add_argument('--moving-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=2)
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    root = args.source_root.resolve()
    sys.path.insert(0, str(root / 'tools'))
    import ou_evidence_provenance as provenance
    from encounter_scaling_study import verify_stationary
    from sim_dataset import REFERENCE_CSV_SHA256, input_provenance, verify_csv_clock

    directories = dict(stationary=args.stationary_dir.resolve(), moving=args.moving_dir.resolve())
    names = sorted(n for n in REFERENCE_CSV_SHA256 if 'jonswap' in n and
                   any(h in n for h in ('H0.270_', 'H1.500_', 'H8.500_')))
    if len(names) != 3:
        raise ValueError('expected small, medium and large JONSWAP screens')
    inputs = dict(stationary=verify_stationary(args.stationary_archive, directories['stationary']),
                  moving=input_provenance(directories['moving'] / n for n in REFERENCE_CSV_SHA256))
    for n in names:
        verify_csv_clock(directories['moving'] / n)
    binaries, sources, builds = {}, set(), {}
    for name, (directory, executable) in FAMILIES.items():
        folder = root / 'tests' / directory
        subprocess.run(['make', '-q', '-C', str(folder), executable], check=True,
                       stdout=subprocess.DEVNULL)
        binary = folder / executable
        binaries[name] = binary
        builds[name] = digest(binary)
        sources.update(provenance.implementation_closure(
            [folder / (executable + '.cpp'), root / 'src/util/W3dSimCommon.cpp']))
        sources.add(folder / 'Makefile')
    source_hashes = {str(p.relative_to(root)): digest(p) for p in sorted(sources)}
    base_env = {k: v for k, v in os.environ.items()
                if not k.startswith(('W3D_', 'SF_', 'OU_', 'TFG_', 'PII_', 'NLO_', 'PROBE_'))}
    base_env.update(W3D_COLLECT_ALL_GATES='1', W3D_WRITE_TIMESERIES='0',
                    W3D_VALIDATION_WINDOW_SEC='900')

    def run(job):
        family, dataset, name, arm = job
        env = {**base_env, **ARMS[arm]}
        if dataset == 'stationary':
            env['W3D_REFERENCE_MOTION'] = 'stationary'
        with tempfile.TemporaryDirectory(prefix='encounter-covariance-') as tmp:
            work = Path(tmp)
            (work / name).symlink_to(directories[dataset] / name)
            result = subprocess.run([str(binaries[family]), '--input', name], cwd=work,
                                    env=env, capture_output=True, text=True, timeout=600)
        lines = [line for line in result.stdout.splitlines()
                 if line.startswith('VALIDATION_METRICS ')]
        if result.returncode not in (0, 1) or len(lines) != 1:
            raise RuntimeError(f'{job}: {result.returncode}\n{result.stdout[-1000:]}{result.stderr[-1000:]}')
        fields = dict(part.split('=', 1) for part in lines[0].split()[1:])
        if (fields['input'], fields['start_s'], fields['window_s'], fields['samples']) != (
                name, '300', '900', '180000'):
            raise ValueError('wrong input or incomplete score window')
        metrics = {key: float(fields[key]) for key in METRICS}
        if not all(math.isfinite(value) for value in metrics.values()):
            raise ValueError('nonfinite RMS')
        print(family, dataset, name, arm, metrics['disp_z_pct_hs'], flush=True)
        return dict(family=family, dataset=dataset, input=name, arm=arm,
                    exit_code=result.returncode, metrics=metrics)

    jobs = [(f, d, n, a) for f in FAMILIES for d in directories for n in names for a in ARMS
            if not (f == 'tfg' and a == 'slower_commits')]
    with ThreadPoolExecutor(args.jobs) as pool:
        rows = list(pool.map(run, jobs))
    if source_hashes != {str(p.relative_to(root)): digest(p) for p in sorted(sources)}:
        raise ValueError('source changed during screen')
    if builds != {f: digest(binaries[f]) for f in FAMILIES}:
        raise ValueError('executable changed during screen')
    report = dict(
        role='non-promoting diagnostic; three sea sizes, two speeds, default sensor seed',
        protocol='1200 s at 200 Hz; full-rate final 900 s; no detrending or fitting',
        source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip(),
        source_sha256=source_hashes, binary_sha256=builds,
        driver_sha256=digest(Path(__file__)), build_environment=provenance.environment_metadata(),
        arms=ARMS, inputs=inputs, rows=rows,
        limits='Stationary gate failures are retained; moving gates do not define stationary acceptance. '
               'This screen can reject a mechanism, not establish all-sea or multi-seed improvement.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
