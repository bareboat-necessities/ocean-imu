#!/usr/bin/env python3
"""Isolate encounter heave errors with explicitly counterfactual physical records.

Original, level-with-original-translation, and level-vertical-only records retain
the identical z/vz/az history. These are causal diagnostics, not shipping gates
or a claim that removing hull motion is a deployable estimator improvement.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
import pandas as pd

import moving_vessel_diagnosis as physics
import ou_evidence_provenance as provenance
from encounter_scaling_study import verify_stationary
from filter_sim_gates import METRICS, REPLAYS
from sim_dataset import ROOT, REFERENCE_CSV_SHA256, digest, input_provenance, verify_csv_clock


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stationary-archive', type=Path, required=True)
    parser.add_argument('--stationary-dir', type=Path, required=True)
    parser.add_argument('--moving-dir', type=Path, required=True)
    parser.add_argument('--families', default='ou2,ou3,carma,tfg')
    parser.add_argument('--jobs', type=int, default=3)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    families = args.families.split(',')
    if args.jobs < 1 or len(set(families)) != len(families) or not set(families) <= {'ou2', 'ou3', 'carma', 'tfg', 'csm', 'tcsm', 'tcarma'}:
        parser.error('distinct supported Kalman families and positive jobs required')
    sources = {Path(__file__).resolve(), ROOT/'tools/moving_vessel_diagnosis.py',
               ROOT/'tools/encounter_scaling_study.py', ROOT/'tools/filter_sim_gates.py',
               ROOT/'tools/stationary_vessel_control.py', ROOT/'tools/sim_dataset.py',
               ROOT/'tools/ou_evidence_provenance.py'}
    binaries = {}
    for family in families:
        folder, binary, _ = REPLAYS[family]
        directory = ROOT/'tests'/folder
        subprocess.run(['make', '-q', '-C', str(directory), binary], check=True, stdout=subprocess.DEVNULL)
        binaries[family] = directory/binary
        sources.update(provenance.implementation_closure([directory/(binary+'.cpp'), ROOT/'src/util/W3dSimCommon.cpp']))
        sources.add(directory/'Makefile')
    source_hashes = {str(p.relative_to(ROOT)): digest(p) for p in sorted(sources)}
    binary_hashes = {k: digest(p) for k, p in binaries.items()}
    directories = {'stationary': args.stationary_dir.resolve(), 'moving': args.moving_dir.resolve()}
    inputs = {'stationary': verify_stationary(args.stationary_archive, directories['stationary']),
              'moving': input_provenance(directories['moving']/n for n in REFERENCE_CSV_SHA256)}
    names = [n for n in REFERENCE_CSV_SHA256 if 'jonswap' in n and 'H4.000' not in n]
    rows, derived = [], []

    def run(family, dataset, name, arm, noise, path):
        env = {k: v for k, v in os.environ.items() if not k.startswith(('W3D_', 'SF_', 'OU_', 'TFG_', 'PROBE_'))}
        env.update(W3D_WRITE_TIMESERIES='0', W3D_COLLECT_ALL_GATES='1', W3D_VALIDATION_WINDOW_SEC='900')
        if dataset == 'stationary' or arm == 'vertical':
            env['W3D_REFERENCE_MOTION'] = 'stationary'
        command = [str(binaries[family]), '--input', str(path)] + (['--no-noise'] if noise == 'none' else [])
        result = subprocess.run(command, cwd=path.parent, env=env, capture_output=True, text=True, timeout=600)
        lines = [line for line in result.stdout.splitlines() if line.startswith('VALIDATION_METRICS ')]
        if result.returncode not in (0, 1) or len(lines) != 1:
            raise RuntimeError(f'{command}: {result.returncode}\n{result.stdout[-2000:]}{result.stderr[-1000:]}')
        fields = dict(item.split('=', 1) for item in lines[0].split()[1:])
        if (Path(fields['input']).name, fields['samples'], fields['start_s'], fields['window_s']) != (name, '180000', '300', '900'):
            raise ValueError('wrong input or incomplete scoring window')
        metrics = {k: float(fields[k]) for k in METRICS}
        if not all(math.isfinite(v) for v in metrics.values()):
            raise ValueError('nonfinite metric')
        row = dict(family=family, dataset=dataset, input=name, arm=arm, noise=noise,
                   input_sha256=digest(path), metrics=metrics, exit_code=result.returncode,
                   gate_failures=[s for s in (result.stdout+'\n'+result.stderr).splitlines()
                                  if s.startswith(('ERROR:', 'QUALITY_GATE: PASS=0'))])
        print(family, dataset, name, arm, noise, metrics['disp_z_pct_hs'], flush=True)
        return row

    report = dict(role='physical/sensor counterfactual localization only; not deployment qualification',
                  protocol='default seeds; 1200 s at 200 Hz, final 900 s; unchanged estimator and z/vz/az history',
                  source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  source_sha256=source_hashes, binary_sha256=binary_hashes, inputs=inputs,
                  build_environment=provenance.environment_metadata(), derived_inputs=derived, rows=rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    partial = args.output.with_suffix('.partial.json')
    for dataset, directory in directories.items():
        for name in names:
            source = directory/name
            verify_csv_clock(source)
            original = pd.read_csv(source, float_precision='round_trip')
            with tempfile.TemporaryDirectory(prefix='encounter-axis-') as temporary:
                work = Path(temporary)
                paths = {'original': source}
                for arm in ('level', 'vertical'):
                    frame = original.copy()
                    if arm == 'vertical':
                        for prefix in ('disp', 'vel', 'acc'):
                            frame[[prefix+'_x', prefix+'_y']] = 0.0
                    q = np.zeros((len(frame), 4))
                    q[:, 0] = 1.0
                    physics.set_body_fields(frame, q, np.zeros((len(frame), 3)))
                    assert frame[['disp_z', 'vel_z', 'acc_z']].equals(original[['disp_z', 'vel_z', 'acc_z']])
                    folder = work/arm
                    folder.mkdir()
                    paths[arm] = folder/name
                    frame.to_csv(paths[arm], index=False, float_format='%.17g')
                    verify_csv_clock(paths[arm])
                    derived.append(dict(dataset=dataset, input=name, arm=arm, sha256=digest(paths[arm])))
                with ThreadPoolExecutor(args.jobs) as pool:
                    tasks = [pool.submit(run, family, dataset, name, arm, noise, path)
                             for family in families for arm, path in paths.items() for noise in ('normal', 'none')]
                    for task in as_completed(tasks):
                        rows.append(task.result())
                        partial.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    if binary_hashes != {k: digest(p) for k, p in binaries.items()} or source_hashes != {str(p.relative_to(ROOT)): digest(p) for p in sorted(sources)}:
        raise ValueError('source or executable changed during diagnosis')
    rows.sort(key=lambda r: tuple(r[k] for k in ('family', 'dataset', 'input', 'arm', 'noise')))
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    partial.unlink()


if __name__ == '__main__':
    main()
