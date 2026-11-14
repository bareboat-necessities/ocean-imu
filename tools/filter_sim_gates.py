#!/usr/bin/env python3
"""Measure the pinned filter replays without changing their acceptance limits.

Build with make first. Run with --run --log-dir /tmp/filter-sim-gates, or
summarize existing logs with --log-dir. Logs and dense time series are not
versioned; the compact report records every scored record and input hash.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile

from sim_dataset import ROOT, REFERENCE_CSV_SHA256, input_provenance, digest


REPLAYS = {
    'ou2': ('kalman_ou_ii', 'kalman_ou_ii-sim', True),
    'ou3': ('kalman_ou_iii', 'kalman_ou_iii-sim', True),
    'tfg': ('kalman_tfg', 'kalman_tfg-sim', True),
    'carma': ('kalman_carma', 'carma3-sim', True),
    'csm': ('kalman_carma_shared_modes', 'kalman_carma_shared_modes-sim', True),
    'tcarma': ('kalman_tfg_carma', 'tfg_carma3-sim', True),
    'tcsm': ('kalman_tfg_shared_modes', 'kalman_tfg_shared_modes-sim', True),
    'tfg_legacy': ('kalman_tfg', 'kalman_tfg-legacy-sim', False),
    'ahrs': ('ahrs', 'ahrs-qmekf-sim', False),
    'heave': ('heave_only', 'heave_only-sim', False),
    'nlo': ('nlo', 'nlo-sim', False),
    'pii': ('pii_observer', 'pii_observer-adaptive', False),
}
METRICS = (
    'disp_z_pct_hs', 'disp_3d_pct_refmax', 'roll_rms_deg', 'pitch_rms_deg',
    'yaw_rms_deg', 'accel_bias_3d_rms_mps2', 'gyro_bias_3d_rms_radps',
    'disp_x_rms_m', 'disp_y_rms_m', 'disp_z_rms_m', 'disp_3d_rms_m',
)
BIAS_PATTERN = re.compile(
    r'Bias error RMS \(% of max TRUE bias\) \((acc|gyro)\): '
    r'X=([^%]+)% Y=([^%]+)% Z=([^%]+)% \|3D\|=([^%]+)%')
GATE_METRICS = {
    'err_limit_percent_z_jonswap': ('disp_z_pct_hs', 'jonswap'),
    'err_limit_percent_z_pmstokes': ('disp_z_pct_hs', 'pmstokes'),
    'err_limit_yaw_deg': ('yaw_rms_deg', None),
    'err_limit_roll_deg': ('roll_rms_deg', None),
    'err_limit_pitch_deg': ('pitch_rms_deg', None),
    'err_limit_percent_3d_jonswap': ('disp_3d_pct_refmax', 'jonswap'),
    'err_limit_percent_3d_pmstokes': ('disp_3d_pct_refmax', 'pmstokes'),
    'acc_z_bias_percent': ('acc_bias_z_percent', None),
    'bias_3d_percent': ('acc_bias_3d_percent', None),
    'gyro_bias_3d_percent': ('gyro_bias_3d_percent', None),
}


def measured_worst(rows: list[dict]) -> dict:
    """Keep normalization and wave-family selection identical to the C++ gates."""
    return {gate: max(row[metric] for row in rows if wave is None or wave in row['input'])
            for gate, (metric, wave) in GATE_METRICS.items()}


def parse_pii(text: str) -> list[dict]:
    """Score the vertical observer without inventing horizontal or bias states."""
    rows = []
    for line in text.splitlines():
        if not line.startswith('PII_METRICS '):
            continue
        fields = dict(item.split('=', 1) for item in line.split()[1:])
        if (fields['start_s'], fields['window_s'], fields['samples']) != ('300', '900', '180000'):
            raise ValueError('PII requires the trailing 900 s of a full replay')
        row = {'input': fields['input'], **{key: float(fields[key]) for key in
               ('disp_z_rms_m', 'disp_z_pct_hs', 'roll_rms_deg', 'pitch_rms_deg', 'yaw_rms_deg')}}
        if not all(math.isfinite(value) for key, value in row.items() if key != 'input'):
            raise ValueError('nonfinite PII metric')
        rows.append(row)
    if len(rows) != 8 or {r['input'] for r in rows} != set(REFERENCE_CSV_SHA256):
        raise ValueError('expected all eight unique PII records')
    return rows


def profiles_for(name: str) -> list[tuple[str, list[str]]]:
    profiles = [('default', [])]
    if REPLAYS[name][2]:
        profiles.append(('legacy', ['--regression']))
    if name == 'heave':
        profiles.extend((('proxy', ['--frontend', 'proxy']), ('truth', ['--frontend', 'truth'])))
    return profiles


def parse_w3d(text: str) -> list[dict]:
    """Require eight complete, finite, 900 s records, including bias scores."""
    rows = []
    for line in text.splitlines():
        if line.startswith('VALIDATION_METRICS '):
            fields = dict(item.split('=', 1) for item in line.split()[1:])
            if (fields['samples'], fields['window_s'], fields['start_s']) != ('180000', '900', '300'):
                raise ValueError('expected the trailing 900 s of a 1200 s, 200 Hz replay')
            rows.append({'input': fields['input'], **{key: float(fields[key]) for key in METRICS}})
        elif match := BIAS_PATTERN.fullmatch(line):
            if not rows:
                raise ValueError('bias summary without a scored record')
            kind, x, y, z, xyz = match.groups()
            rows[-1].update({f'{kind}_bias_{axis}_percent': float(value)
                             for axis, value in zip(('x', 'y', 'z', '3d'), (x, y, z, xyz))})
    if len(rows) != 8 or {r['input'] for r in rows} != set(REFERENCE_CSV_SHA256):
        raise ValueError('expected each of the eight pinned wave records exactly once')
    for row in rows:
        for key in (*METRICS, *(f'{kind}_bias_{axis}_percent'
                               for kind in ('acc', 'gyro') for axis in ('x', 'y', 'z', '3d'))):
            if key not in row or not math.isfinite(row[key]):
                raise ValueError(f"{row['input']}: missing/nonfinite {key}")
    return rows


def collect(log_dir: Path, jobs: int, *, measure: bool = False) -> None:
    """Run already built binaries; numeric failures never become new limits."""
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(('W3D_', 'SF_', 'OU_', 'HO_', 'NLO_', 'TFG_'))}
    env.update(W3D_COLLECT_ALL_GATES='1', W3D_WRITE_TIMESERIES='0',
               W3D_VALIDATION_WINDOW_SEC='900')

    def run(item):
        name, (folder, binary, _) = item
        here = ROOT / 'tests' / folder
        binary_sha = digest(here / binary)
        profiles = profiles_for(name)
        results = []
        # Only the eight scored records belong to this protocol. An isolated
        # working directory also keeps generated CSVs out of the source tree.
        with tempfile.TemporaryDirectory(prefix=f'filter-gates-{name}-') as tmp:
            work = Path(tmp)
            for record in REFERENCE_CSV_SHA256:
                (work / record).symlink_to(here / record)
            for profile, args in profiles:
                command = [str(here / binary), *args]
                with (log_dir / f'{name}-{profile}.log').open('w') as stream:
                    result = subprocess.run(command, cwd=work, env=env, stdout=stream,
                                            stderr=subprocess.STDOUT, check=False)
                results.append({'profile': profile, 'exit_code': result.returncode,
                                'log_sha256': digest(log_dir / f'{name}-{profile}.log')})
                print(f'{name}/{profile}: exit {result.returncode}', flush=True)
        return {'name': name, 'binary_sha256': binary_sha, 'profiles': results}

    log_dir.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        runs = list(pool.map(run, REPLAYS.items()))
    (log_dir / 'runs.json').write_text(json.dumps(runs, indent=2) + '\n')
    if not measure and any(p['exit_code'] for r in runs for p in r['profiles']):
        raise ValueError('a replay failed; inspect logs, do not refit failing gates')


def split_make_log(text: str) -> dict[tuple[str, str], str]:
    """Recover the same replay profiles from make all, without replaying twice."""
    if any(message in text for message in ('QUALITY_GATE: PASS=0',
                                           'ERROR: raw Z RMS above limit',
                                           'ERROR: nonfinite scored NLO heave error')):
        raise ValueError('make log contains a failed numeric replay')
    families = dict(zip(('OU_II', 'OU_III', 'TFG', 'CARMA3_A41',
                         'CARMA_SHARED_MODES_CSM', 'TFG_CARMA3_T41',
                         'TFG_CARMA_SHARED_MODES_CSM', 'TFG_LEGACY'),
                        ('ou2', 'ou3', 'tfg', 'carma', 'csm', 'tcarma', 'tcsm', 'tfg_legacy')))
    parts = {}
    block = []
    profile = 'default'

    def flush():
        names = {dict(item.split('=', 1) for item in line.split()[1:])['family']
                 for line in block if line.startswith('VALIDATION_METRICS ')}
        if not names:
            return
        if len(names) != 1 or not names <= families.keys():
            raise ValueError(f'unknown or mixed replay families: {names}')
        key = (families[next(iter(names))], profile)
        if key in parts:
            raise ValueError(f'duplicate replay profile in make log: {key}')
        parts[key] = '\n'.join(block) + '\n'

    ahrs = []
    heave = {'mahony': [], 'proxy': [], 'truth': []}
    nlo = []
    pii = []
    pii_active = False
    ahrs_heading = False
    for line in text.splitlines():
        if line.startswith('Running ') and line.endswith('/run_tests.sh'):
            flush()
            block = []
            pii_active = False
        if line.startswith('AdaptiveVerticalPIIMahony simulation starting'):
            flush()
            block = []
            pii_active = True
        if line.startswith(('GYRO_REPLAY_PROFILE ', 'REGRESSION_PROTOCOL:',
                            'TFG legacy comparison simulation starting')):
            flush()
            block = []
            profile = 'legacy' if 'legacy_interval' in line else 'default'
        block.append(line)
        if line.startswith('=== Last 60 s RMS summary'):
            ahrs_heading = True
        elif ahrs_heading and line.startswith('Angles RMS (deg):'):
            ahrs.append(line)
            ahrs_heading = False
        if line.startswith('HEAVE_ONLY '):
            fields = dict(item.split('=', 1) for item in line.split()[1:])
            heave[fields['frontend']].append(line)
        if line.startswith('Z RMS raw (%Hs):'):
            nlo.append(line)
        if line.startswith('PII_METRICS ') or (pii_active and line.startswith('QUALITY_GATE:')):
            pii.append(line)
    flush()
    parts['ahrs', 'default'] = '\n'.join(ahrs) + '\n'
    parts['nlo', 'default'] = '\n'.join(nlo) + '\n'
    parts['pii', 'default'] = '\n'.join(pii) + '\n'
    for frontend, lines in heave.items():
        parts['heave', 'default' if frontend == 'mahony' else frontend] = '\n'.join(lines) + '\n'
    expected = {(name, profile) for name in REPLAYS for profile, _ in profiles_for(name)}
    if parts.keys() != expected:
        raise ValueError(f'incomplete make log profiles: {expected - parts.keys()}')
    return parts


def import_make_log(source: Path, log_dir: Path) -> None:
    parts = split_make_log(source.read_text())
    log_dir.mkdir(parents=True, exist_ok=True)
    runs = []
    for name, (folder, binary, _) in REPLAYS.items():
        profiles = []
        for profile, _ in profiles_for(name):
            path = log_dir / f'{name}-{profile}.log'
            path.write_text(parts[name, profile])
            profiles.append({'profile': profile, 'exit_code': 0, 'log_sha256': digest(path)})
        runs.append({'name': name, 'binary_sha256': digest(ROOT / 'tests' / folder / binary),
                     'profiles': profiles})
    (log_dir / 'runs.json').write_text(json.dumps(runs, indent=2) + '\n')


def summarize(log_dir: Path, *, measure: bool = False) -> dict:
    runs = json.loads((log_dir / 'runs.json').read_text())
    if {r['name'] for r in runs} != set(REPLAYS):
        raise ValueError('incomplete filter inventory')
    records = {}
    for run in runs:
        name = run['name']
        if 'profiles' not in run or any(p['exit_code'] not in ((0, 1) if measure else (0,)) for p in run['profiles']):
            raise ValueError(f'{name}: unsuccessful replay')
        if [p['profile'] for p in run['profiles']] != [p for p, _ in profiles_for(name)]:
            raise ValueError(f'{name}: incomplete profile inventory')
        if name not in ('ahrs', 'heave', 'nlo', 'pii'):
            records[name] = {}
            for profile in run['profiles']:
                text = (log_dir / f"{name}-{profile['profile']}.log").read_text()
                count = len(re.findall(r'^QUALITY_GATE: PASS=[01] ', text, re.M)) if measure else text.count('QUALITY_GATE: PASS=1 REASON=ok')
                if count != 8 or 'QUALITY_GATE: SKIPPED' in text:
                    raise ValueError(f"{name}/{profile['profile']}: not eight passing numeric gates")
                records[name][profile['profile']] = parse_w3d(text)
    other = {}
    pii = (log_dir / 'pii-default.log').read_text()
    gate_count = len(re.findall(r'^QUALITY_GATE: PASS=[01] REASON=', pii, re.M)) if measure else pii.count('QUALITY_GATE: PASS=1 REASON=ok')
    if gate_count != 8 or 'QUALITY_GATE: SKIPPED' in pii:
        raise ValueError('PII did not complete eight gated replays')
    other['pii'] = parse_pii(pii)
    ahrs = (log_dir / 'ahrs-default.log').read_text()
    angles = re.findall(r'Angles RMS \(deg\): Roll=(\S+) Pitch=(\S+) Yaw=(\S+)', ahrs)
    if len(angles) != 8:
        raise ValueError('AHRS did not score all eight records')
    other['ahrs'] = [{'input': name, **dict(zip(('roll', 'pitch', 'yaw'), map(float, values)))}
                     for name, values in zip(sorted(REFERENCE_CSV_SHA256), angles)]
    if not all(math.isfinite(row[a]) for row in other['ahrs'] for a in ('roll', 'pitch', 'yaw')):
        raise ValueError('AHRS has nonfinite scored attitude')
    other['heave'] = {}
    for profile, count in (('default', 40), ('proxy', 40), ('truth', 32)):
        lines = (log_dir / f'heave-{profile}.log').read_text().splitlines()
        rows = [dict(item.split('=', 1) for item in line.split()[1:])
                for line in lines if line.startswith('HEAVE_ONLY ')]
        if len(rows) != count:
            raise ValueError(f'heave/{profile}: expected {count} scores')
        for method in {r['method'] for r in rows}:
            names = [Path(r['record']).name for r in rows if r['method'] == method]
            if len(names) != 8 or set(names) != set(REFERENCE_CSV_SHA256):
                raise ValueError(f'heave/{profile}/{method}: incomplete records')
        if not all(math.isfinite(float(r['z_pct_hs'])) for r in rows):
            raise ValueError(f'heave/{profile}: nonfinite RMS')
        other['heave'][profile] = [{k: row[k] for k in ('method', 'frontend', 'record', 'z_pct_hs')}
                                   for row in rows]
    nlo = (log_dir / 'nlo-default.log').read_text()
    other['nlo_z_percent_hs'] = [float(v) for v in re.findall(r'Z RMS raw \(%Hs\): (\S+)%', nlo)]
    if len(other['nlo_z_percent_hs']) != 8 or not all(math.isfinite(v) for v in other['nlo_z_percent_hs']):
        raise ValueError('NLO did not score all eight records')
    here = ROOT / 'tests/kalman_ou_iii'
    return {
        'source_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'compiler': subprocess.check_output(['g++', '--version'], text=True).splitlines()[0],
        'acceptance': 'measurement only; gate failures retained' if measure else 'all replay gates passed',
        'replay_exit_codes': {r['name']: {p['profile']: p['exit_code'] for p in r['profiles']} for r in runs},
        'protocol': '200 Hz; 1200 s; trailing 900 s; default noise/seeds; magnetometer on except NLO; AHRS trailing 60 s',
        'inputs': input_provenance([here / name for name in REFERENCE_CSV_SHA256]),
        'binary_sha256': {r['name']: r['binary_sha256'] for r in runs},
        'records': records, 'other_replays': other,
        'diagnostics': {
            'tests/kalman_carma/carma_process-sim.cpp': 'stochastic prior generator; independent spectral/derivative/variance gates in carma_process_validation.py',
            'tests/kalman_carma/carma_lever_arm-sim.cpp': 'off-CG ablation instrument; no default accuracy gate across study configurations; lever-arm numerical tests run in make all',
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--log-dir', type=Path, required=True)
    parser.add_argument('--run', action='store_true')
    parser.add_argument('--measure', action='store_true',
                        help='report a changed dataset without promoting failed gates; complete finite records still required')
    parser.add_argument('--make-log', type=Path,
                        help='import a make all log produced with W3D_VALIDATION_WINDOW_SEC=900')
    parser.add_argument('--jobs', type=int, default=3)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.run and args.make_log:
        parser.error('--run and --make-log are mutually exclusive')
    for folder, _, _ in REPLAYS.values():
        input_provenance([ROOT / 'tests' / folder / name for name in REFERENCE_CSV_SHA256])
    if args.run:
        collect(args.log_dir, args.jobs, measure=args.measure)
    if args.make_log:
        import_make_log(args.make_log, args.log_dir)
    report = summarize(args.log_dir, measure=args.measure)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(f'Wrote {args.output}')


if __name__ == '__main__':
    main()
