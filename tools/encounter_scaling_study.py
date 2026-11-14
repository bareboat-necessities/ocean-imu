#!/usr/bin/env python3
"""Pair exact-clock zero-speed/5 kn replays of an encounter-scaling adjustment.

Build the selected family first. Both arms use the same executable and sensor
seeds. The previous arm disables the selected adjustment or uses an explicitly
source-bound previous executable for filters without a runtime ablation. Statistics
always use the full-rate final 900 s, never the decimated drawing files.
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
import zipfile

import ou_evidence_provenance as provenance
from sim_dataset import (ROOT, REFERENCE_CSV_SHA256, digest, input_provenance,
                         verify_csv_clock)
from stationary_vessel_control import ARCHIVE_SHA256 as STATIONARY_ARCHIVE_SHA256

METRICS = ('disp_z_rms_m', 'disp_z_pct_hs', 'roll_rms_deg', 'pitch_rms_deg', 'yaw_rms_deg')
FAMILIES = {
    'pii': ('pii_observer', 'pii_observer-adaptive', 'PII_METRICS',
            'PII_PHASE_MATCH', 'phase_matched'),
    'nlo': ('nlo', 'nlo-sim', 'NLO_METRICS', 'NLO_THETA_ACCEL_SCALE', 'energy_adjusted'),
    **{name: (folder, binary, 'VALIDATION_METRICS', None, 'adjusted')
       for name, folder, binary in (
           ('ou2', 'kalman_ou_ii', 'kalman_ou_ii-sim'),
           ('ou3', 'kalman_ou_iii', 'kalman_ou_iii-sim'),
           ('tfg', 'kalman_tfg', 'kalman_tfg-sim'),
           ('tfg_legacy', 'kalman_tfg', 'kalman_tfg-legacy-sim'),
           ('carma', 'kalman_carma', 'carma3-sim'),
           ('csm', 'kalman_carma_shared_modes', 'kalman_carma_shared_modes-sim'),
           ('tcarma', 'kalman_tfg_carma', 'tfg_carma3-sim'),
           ('tcsm', 'kalman_tfg_shared_modes', 'kalman_tfg_shared_modes-sim'))},
}
for _family, _folder, _binary in (
    ('csm', 'kalman_carma_shared_modes', 'kalman_carma_shared_modes-sim'),
    ('tcsm', 'kalman_tfg_shared_modes', 'kalman_tfg_shared_modes-sim'),
    ('tangent', 'kalman_tangent_shared_modes', 'tangent-replay'),
):
    FAMILIES[_family] = (_folder, _binary, 'VALIDATION_METRICS',
                         'W3D_CARMA_BAND_NORMALIZATION', 'band_normalized')


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
               W3D_VALIDATION_WINDOW_SEC='900')
    _, _, prefix, setting, current_arm = FAMILIES[job['family']]
    if job['previous'] and setting is not None and job['previous_binary'] is None:
        env[setting] = '0'
    if job['engine_rpm'] is not None:
        env['W3D_ENGINE_RPM'] = str(job['engine_rpm'])
    if job['dataset'] == 'stationary':
        env['W3D_REFERENCE_MOTION'] = 'stationary'
    if job['seed'] != 'default':
        env['W3D_SEED'] = job['seed']
    with tempfile.TemporaryDirectory(prefix=f'encounter-{job["family"]}-') as temporary:
        work = Path(temporary)
        (work / job['input'].name).symlink_to(job['input'])
        binary = job['previous_binary'] if job['previous'] and job['previous_binary'] else job['binary']
        command = [str(binary)] + (['--no-noise'] if job['noise'] == 'none' else [])
        if prefix == 'VALIDATION_METRICS':
            command += ['--input', job['input'].name]
        result = subprocess.run(command, cwd=work, env=env,
                                capture_output=True, text=True, timeout=600, check=False)
    lines = [line for line in result.stdout.splitlines() if line.startswith(prefix + ' ')]
    if result.returncode not in (0, 1) or len(lines) != 1:
        raise RuntimeError(f'{job}: exit {result.returncode}\n{result.stdout[-2000:]}{result.stderr[-2000:]}')
    fields = dict(item.split('=', 1) for item in lines[0].split()[1:])
    if fields['input'] != job['input'].name:
        raise ValueError('scored input differs from the requested record')
    if (fields['samples'], fields['start_s'], fields['window_s']) != ('180000', '300', '900'):
        raise ValueError('incomplete replay window')
    keys = METRICS + (('disp_3d_rms_m', 'accel_bias_3d_rms_mps2',
                       'gyro_bias_3d_rms_radps') if prefix == 'VALIDATION_METRICS' else ())
    metrics = {key: float(fields[key]) for key in keys}
    if not all(math.isfinite(value) for value in metrics.values()):
        raise ValueError('nonfinite scored metric')
    row = dict(dataset=job['dataset'], input=job['input'].name, seed=job['seed'],
               arm='previous' if job['previous'] else current_arm,
               exit_code=result.returncode, metrics=metrics,
               gate_failures=[s for s in (result.stdout + '\n' + result.stderr).splitlines()
                              if s.startswith('ERROR:') or s.startswith('QUALITY_GATE: PASS=0')])
    print(row['dataset'], row['input'], row['seed'], row['arm'], metrics['disp_z_pct_hs'], flush=True)
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--family', choices=FAMILIES, default='pii')
    parser.add_argument('--noise', choices=('normal', 'none'), default='normal')
    parser.add_argument('--previous-binary', type=Path)
    parser.add_argument('--previous-source-record', type=Path,
                        help='JSON with source_commit, source_sha256 and binary_sha256 captured before editing')
    parser.add_argument('--engine-rpm', type=float)
    parser.add_argument('--stationary-archive', type=Path, required=True)
    parser.add_argument('--stationary-dir', type=Path, required=True)
    parser.add_argument('--moving-dir', type=Path, required=True)
    parser.add_argument('--seeds', default='default,7,23')
    parser.add_argument('--jobs', type=int, default=3)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--resume', action='store_true', help='resume a matching .partial.json checkpoint')
    args = parser.parse_args()
    if args.engine_rpm is not None and (not math.isfinite(args.engine_rpm)
                                      or args.engine_rpm <= 0 or args.noise == 'none'):
        parser.error('engine RPM must be finite and positive, with normal sensor noise')
    if bool(args.previous_binary) != bool(args.previous_source_record):
        parser.error('previous executable and its source record must be supplied together')
    if FAMILIES[args.family][3] is None and args.previous_binary is None:
        parser.error('this family requires an explicitly source-bound previous executable')
    if args.family == 'tfg_legacy' and args.noise == 'none':
        parser.error('the legacy TFG harness does not support noise-free replay')
    previous_record = None
    if args.previous_binary is not None:
        import hashlib
        args.previous_binary = args.previous_binary.resolve()
        previous_record = json.loads(args.previous_source_record.read_text())
        if digest(args.previous_binary) != previous_record['binary_sha256']:
            raise ValueError('previous executable differs from its source record')
        previous_commit = subprocess.check_output(
            ['git', 'rev-parse', '--verify', previous_record['source_commit']+'^{commit}'],
            cwd=ROOT, text=True).strip()
        roots = [ROOT/'tests'/FAMILIES[args.family][0]/(FAMILIES[args.family][1]+'.cpp'),
                 ROOT/'src/util/W3dSimCommon.cpp']
        required = {str(p.relative_to(ROOT)) for p in provenance.implementation_closure(roots)}
        if not required <= previous_record['source_sha256'].keys():
            raise ValueError('previous source record omits replay implementation dependencies')
        for name, expected in previous_record['source_sha256'].items():
            content = subprocess.check_output(['git', 'show', previous_commit+':'+name], cwd=ROOT)
            if hashlib.sha256(content).hexdigest() != expected:
                raise ValueError(f'previous source is not the recorded commit: {name}')
    seeds = args.seeds.split(',')
    if args.jobs < 1 or len(set(seeds)) != len(seeds) or any(
            seed != 'default' and (not seed.isdecimal() or int(seed) > 2**32-1) for seed in seeds):
        parser.error('positive jobs and distinct default/unsigned integer seeds required')
    directory, executable, _, _, current_arm = FAMILIES[args.family]
    folder = ROOT / 'tests' / directory
    subprocess.run(['make', '-q', '-C', str(folder), executable], check=True,
                   stdout=subprocess.DEVNULL)
    binary = folder / executable
    initial_binary_hash = digest(binary)
    directories = dict(stationary=args.stationary_dir.resolve(), moving=args.moving_dir.resolve())
    inputs = dict(stationary=verify_stationary(args.stationary_archive, directories['stationary']),
                  moving=input_provenance(directories['moving']/name for name in REFERENCE_CSV_SHA256))
    for name in REFERENCE_CSV_SHA256:
        verify_csv_clock(directories['moving'] / name)
    roots = [folder/(executable+'.cpp'), ROOT/'src/util/W3dSimCommon.cpp']
    sources = set(provenance.implementation_closure(roots)) | {
        Path(__file__).resolve(), folder/'Makefile', ROOT/'tools/sim_dataset.py',
        ROOT/'tools/stationary_vessel_control.py', ROOT/'tools/ou_evidence_provenance.py'}
    source_hashes = {str(p.relative_to(ROOT)): digest(p) for p in sorted(sources)}
    jobs = [dict(binary=binary, dataset=dataset, input=directory/name, seed=seed,
                 previous=previous, previous_binary=args.previous_binary,
                 family=args.family, noise=args.noise, engine_rpm=args.engine_rpm)
            for dataset, directory in directories.items() for name in REFERENCE_CSV_SHA256
            for seed in seeds for previous in (True, False)]
    # Long paired studies survive an interrupted session. A checkpoint cannot
    # cross an executable, source, input, seed or experiment-setting change.
    binding = dict(family=args.family, noise=args.noise, engine_rpm=args.engine_rpm,
                   seeds=seeds, inputs=inputs, binary_sha256=initial_binary_hash,
                   previous_build=previous_record, source_sha256=source_hashes)
    partial = args.output.with_suffix('.partial.json')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    if args.resume and partial.exists():
        saved = json.loads(partial.read_text())
        if saved['binding'] != binding:
            raise ValueError('checkpoint does not match this experiment')
        rows = saved['rows']
    def row_key(row):
        return (row['dataset'], row['input'], row['seed'], row['arm'])
    done = {row_key(row) for row in rows}
    expected = {(j['dataset'], j['input'].name, j['seed'],
                 'previous' if j['previous'] else current_arm) for j in jobs}
    if len(done) != len(rows) or not done <= expected:
        raise ValueError('checkpoint has duplicate or unexpected rows')
    pending = [j for j in jobs if (j['dataset'], j['input'].name, j['seed'],
               'previous' if j['previous'] else current_arm) not in done]
    with ThreadPoolExecutor(args.jobs) as pool:
        futures = [pool.submit(run_one, j) for j in pending]
        for future in as_completed(futures):
            rows.append(future.result())
            tmp = partial.with_suffix('.tmp')
            tmp.write_text(json.dumps(dict(binding=binding, rows=rows), allow_nan=False)+'\n')
            tmp.replace(partial)
    rows.sort(key=row_key)
    if previous_record is not None and digest(args.previous_binary) != previous_record['binary_sha256']:
        raise ValueError('previous executable changed during the paired experiment')
    if digest(binary) != initial_binary_hash:
        raise ValueError('executable changed during the paired experiment')
    if source_hashes != {str(p.relative_to(ROOT)): digest(p) for p in sorted(sources)}:
        raise ValueError('replay source changed during the paired experiment')
    pairs = []
    for previous in rows:
        if previous['arm'] != 'previous':
            continue
        after = next(row for row in rows if row['arm'] == current_arm and
                     all(row[key] == previous[key] for key in ('dataset', 'input', 'seed')))
        before, current = previous['metrics']['disp_z_pct_hs'], after['metrics']['disp_z_pct_hs']
        pairs.append(dict(dataset=previous['dataset'], input=previous['input'], seed=previous['seed'],
                          previous_heave_pct_hs=before, **{current_arm+'_heave_pct_hs': current},
                          relative_change_percent=100*(current/before-1)))
    report = dict(protocol='200 Hz, 1200 s; final 900 s/180000 samples; paired sensor/initialization seeds',
                  family=args.family, noise=args.noise, engine_rpm=args.engine_rpm,
                  previous_build=previous_record,
                  role='controlled paired experiment; moving-data acceptance bars are not stationary/noise-seed bars',
                  seeds=seeds, inputs=inputs, binary_sha256=initial_binary_hash,
                  git_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  source_sha256=source_hashes,
                  build_environment=provenance.environment_metadata(), rows=rows, pairs=pairs)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    partial.unlink(missing_ok=True)


if __name__ == '__main__':
    main()
