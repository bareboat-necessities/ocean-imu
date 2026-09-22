#!/usr/bin/env python3
"""Execute paired RAO coefficient cases with hash-checked, per-case resume.

No quality gate is altered. Every attempted case, including process failures and
non-finite physical scores, is retained. The cache identity includes executable,
producer, pinned input, scoring protocol, random seeds and effective overrides.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import time

from model_mismatch_ablation import FAMILIES, RECORDS, source_commit
from ou_validation import parse_validation_metrics
from sim_dataset import input_provenance

ROOT = Path(__file__).resolve().parents[1]
METRICS = (
    'disp_z_pct_hs', 'disp_z_rms_m', 'disp_3d_rms_m', 'roll_rms_deg',
    'pitch_rms_deg', 'yaw_rms_deg', 'accel_3d_rms_mps2',
    'accel_bias_3d_rms_mps2', 'gyro_bias_3d_rms_radps',
)
ENV_PREFIXES = ('W3D_', 'TFG_', 'SF_', 'OU_II_', 'OU_III_', 'OU_SIGMA_')
EXTRA_KNOBS = {
    'W3D_AW_COV_SYNC', 'OU_SIGMA_STILL_DECAY_SEC',
    'OU_SIGMA_VAR_K_PERIODS', 'OU_SIGMA_VAR_HORIZON_MIN_S',
    'OU_SIGMA_VAR_HORIZON_MAX_S',
}


def digest_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical(data) -> bytes:
    return json.dumps(data, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()


def clean_nonfinite(data):
    if isinstance(data, dict):
        return {k: clean_nonfinite(v) for k, v in data.items()}
    if isinstance(data, (list, tuple)):
        return [clean_nonfinite(v) for v in data]
    if isinstance(data, float) and not math.isfinite(data):
        return None
    return data


def write_json(path: Path, data) -> None:
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_bytes(canonical(clean_nonfinite(data)) + b'\n')
    temporary.replace(path)


def validate_configs(configs: list[dict], binary_bytes: bytes) -> None:
    names = set()
    if not configs:
        raise ValueError('empty configuration list')
    for config in configs:
        name = config['name']
        if not re.fullmatch(r'[A-Za-z0-9_.+-]+', name) or name in names:
            raise ValueError(f'unsafe or duplicate configuration name: {name}')
        names.add(name)
        for key, value in config['env'].items():
            if (not key.startswith(('SF_', 'OU_II_', 'OU_III_', 'TFG_'))
                    and key not in EXTRA_KNOBS) or any(
                        word in key for word in ('GATE', 'LIMIT', 'SEED')):
                raise ValueError(f'disallowed override {key}')
            if not isinstance(value, str) or '\x00' in value:
                raise ValueError(f'override must be a NUL-free string: {key}')
            if key.encode() + b'\x00' not in binary_bytes:
                raise ValueError(f'{key} is not exposed by the executable')


def execute(configs: list[dict], family: str, record_indices: list[int],
            seeds: list[str], out: Path, jobs: int, cache: Path) -> list[dict]:
    if jobs < 1 or not seeds or len(set(seeds)) != len(seeds):
        raise ValueError('invalid worker count or duplicated/empty seed set')
    if len(set(record_indices)) != len(record_indices):
        raise ValueError('duplicate record indices')
    out.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)
    producer_hash = digest_bytes(Path(__file__).read_bytes())
    revision = source_commit()
    source_diff = subprocess.check_output(
        ['git', 'diff', 'HEAD', '--', 'src', 'tools', 'tests'], cwd=ROOT)
    (out / 'source.patch').write_bytes(source_diff)
    binary = FAMILIES[family]
    binary_bytes = binary.read_bytes()
    binary_hash = digest_bytes(binary_bytes)
    validate_configs(configs, binary_bytes)
    records = [RECORDS[i] for i in record_indices]
    paths = {r.filename: (binary.parent / r.filename).resolve() for r in records}
    provenance = input_provenance(paths.values())
    protocol = {'window_s': 900, 'write_timeseries': False,
                'collect_all_gates': True, 'timeout_s': 600,
                'seed_mapping': 'nondefault seeds set both W3D_INIT_SEED and W3D_IMU_SEED'}
    manifest = dict(source_commit=revision, source_dirty=bool(source_diff),
                    source_diff_sha256=digest_bytes(source_diff),
                    producer_sha256=producer_hash, binary_sha256=binary_hash,
                    simulation_provenance=provenance, configs=configs,
                    family=family, record_indices=record_indices, seeds=seeds,
                    protocol=protocol, status='executing')
    write_json(out / 'configs.json', configs)
    write_json(out / 'manifest.json', manifest)
    base_environment = {k: v for k, v in os.environ.items()
                        if not k.startswith(ENV_PREFIXES)}
    tasks = [(c, r, s) for c in configs for r in records for s in seeds]

    def run(task):
        config, record, seed = task
        stem = f"{config['name']}_{record.spectrum}_{record.hs_m}_{seed}"
        overrides = dict(W3D_WRITE_TIMESERIES='0', W3D_COLLECT_ALL_GATES='1',
                         W3D_VALIDATION_WINDOW_SEC='900')
        if seed != 'default':
            overrides.update(W3D_INIT_SEED=seed, W3D_IMU_SEED=seed)
        overrides.update(config['env'])
        identity = dict(producer_sha256=producer_hash, binary_sha256=binary_hash,
                        input_sha256=provenance['input_sha256'][record.filename],
                        family=family, input=record.filename, seed=seed,
                        effective_overrides=overrides, protocol=protocol)
        key = digest_bytes(canonical(identity))
        cached_path = cache / (key + '.json')
        cached_log = cache / (key + '.log')
        if cached_path.exists():
            envelope = json.loads(cached_path.read_text())
            payload = envelope['payload']
            if envelope['identity'] != identity or envelope['payload_sha256'] != digest_bytes(canonical(payload)):
                raise ValueError(f'cache identity/checksum mismatch: {key}')
            if not cached_log.exists() or digest_bytes(cached_log.read_bytes()) != payload['log_sha256']:
                raise ValueError(f'cache log mismatch: {key}')
            row = dict(payload, config=config['name'], env=config['env'], reused=True,
                       cache_key=key, cache_source_commit=payload['source_commit'])
            (out / (stem + '.log')).write_bytes(cached_log.read_bytes())
            write_json(out / (stem + '.json'), row)
            return row
        start = time.monotonic()
        error = None
        try:
            completed = subprocess.run(
                [str(binary), '--input', str(paths[record.filename])],
                cwd=binary.parent, env=dict(base_environment, **overrides),
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
                timeout=protocol['timeout_s'], check=False)
            log, code = completed.stdout, completed.returncode
        except subprocess.TimeoutExpired as exc:
            log = exc.stdout or ''
            if isinstance(log, bytes):
                log = log.decode(errors='replace')
            log += '\nREPLAY_EXECUTION_ERROR: timeout\n'
            code, error = -1, 'timeout'
        except OSError as exc:
            code, error, log = -1, str(exc), 'REPLAY_EXECUTION_ERROR: ' + str(exc)
        gates = re.findall(r'QUALITY_GATE: PASS=([01])', log)
        if code not in (0, 1) or not gates or int(gates[-1]) != (1 if code == 0 else 0):
            error = error or f'execution/quality-gate contract failed: exit={code}, gates={gates}'
        try:
            metrics = clean_nonfinite(parse_validation_metrics(log))
        except (ValueError, KeyError) as exc:
            metrics = {}
            error = error or f'metric parse error: {exc}'
        invalid = [m for m in METRICS
                   if not isinstance(metrics.get(m), (int, float)) or metrics[m] < 0]
        if invalid:
            error = error or 'missing/nonfinite physical metrics: ' + ','.join(invalid)
        if metrics.get('samples') != 180000 or metrics.get('window_s') != 900:
            error = error or 'incorrect validation scoring window'
        violations = []
        for line in log.splitlines():
            if line.startswith('ERROR:'):
                match = re.search(r'\(([0-9.eE+-]+)[^>]*>\s*([0-9.eE+-]+)', line)
                ratio = None
                if match and float(match[2]) != 0:
                    ratio = float(match[1]) / float(match[2])
                violations.append(dict(message=line, ratio=ratio))
        payload = dict(source_commit=revision, config=config['name'], env=config['env'],
                       family=family, input=record.filename, seed=seed,
                       exit_code=code, seconds=time.monotonic()-start,
                       violations=violations, metrics=metrics, error=error,
                       log_sha256=digest_bytes(log.encode()),
                       binary_sha256=binary_hash, input_sha256=identity['input_sha256'])
        payload = clean_nonfinite(payload)
        cached_log.write_text(log)
        write_json(cached_path, dict(identity=identity, payload=payload,
                                    payload_sha256=digest_bytes(canonical(payload))))
        row = dict(payload, cache_key=key, reused=False)
        (out / (stem + '.log')).write_text(log)
        write_json(out / (stem + '.json'), row)
        return row

    rows = []
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        for future in as_completed([pool.submit(run, task) for task in tasks]):
            row = future.result()
            rows.append(row)
            print(f"{len(rows)}/{len(tasks)} {row['config']} {row['input']} "
                  f"seed={row['seed']} exit={row['exit_code']} reused={row['reused']} "
                  f"error={row.get('error')}", flush=True)
    rows.sort(key=lambda r: (r['config'], r['input'], r['seed']))
    write_json(out / 'runs.json', rows)
    if binary.read_bytes() != binary_bytes:
        raise RuntimeError('simulator changed during replay')
    summary = []
    for config in configs:
        cases = [r for r in rows if r['config'] == config['name']]
        eligible = all(r.get('error') is None for r in cases)
        ratios = [v['ratio'] for r in cases for v in r['violations'] if v['ratio'] is not None]
        summary.append(dict(config=config['name'], eligible=eligible,
                            failed_records=sum(r['exit_code'] != 0 for r in cases),
                            execution_errors=sum(r.get('error') is not None for r in cases),
                            violations=sum(len(r['violations']) for r in cases),
                            worst_ratio=max(ratios, default=1),
                            **{m: sum(r['metrics'][m] for r in cases)/len(cases)
                               if eligible else None for m in METRICS}))
    write_json(out / 'summary.json', summary)
    manifest.update(status='complete', expected_cases=len(tasks), completed_cases=len(rows),
                    cache_hits=sum(r['reused'] for r in rows),
                    execution_errors=sum(r.get('error') is not None for r in rows),
                    runs_sha256=digest_bytes((out / 'runs.json').read_bytes()))
    write_json(out / 'manifest.json', manifest)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--configs', type=Path, required=True)
    parser.add_argument('--family', choices=FAMILIES, default='TFG')
    parser.add_argument('--records', default='0,1,2,3,4,5,6,7')
    parser.add_argument('--seeds', required=True)
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--cache-dir', type=Path)
    args = parser.parse_args()
    execute(json.loads(args.configs.read_text()), args.family,
            [int(i) for i in args.records.split(',')], args.seeds.split(','),
            args.output_dir.resolve(), args.jobs,
            (args.cache_dir or args.output_dir / 'cache').resolve())


if __name__ == '__main__':
    main()
