#!/usr/bin/env python3
"""Score engine vibration, displaced sensors and a continuous sea-state round trip."""
import argparse
import csv
import json
import math
import os
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from carma_lever_arm_study import AXES, arm_env
from ou_roundtrip_transition import (
    DIAGNOSTIC_HIGH_HEIGHT_M, DIAGNOSTIC_SOURCE_HEIGHT_M, ROUNDTRIP_DURATION_SEC,
    ROUNDTRIP_WAVE_NAME, make_roundtrip_wave, roundtrip_score_rows, roundtrip_segments,
)
import ou_validation as ov
from ou_evidence_provenance import environment_metadata, implementation_closure
from sim_dataset import CASE, REFERENCE_CSV_SHA256, ROOT, digest, input_provenance

OUT = ROOT / 'reports/results/tfg_shared_modes/degradation'
BINARY = ROOT / 'tests/kalman_tfg_shared_modes/installation-study'
METRICS = ('disp_x_rms_m', 'disp_y_rms_m', 'disp_z_rms_m', 'disp_3d_rms_m',
           'roll_rms_deg', 'pitch_rms_deg', 'yaw_rms_deg')


def score(path, extra):
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(('W3D_', 'CSM_', 'OU_', 'SF_', 'TFG_'))}
    env.update(OCEAN_IMU_SIM_CASE=CASE, W3D_IMU_SEED='101', W3D_INIT_SEED='1009', **extra)
    process = subprocess.run([BINARY, path], env=env, text=True, capture_output=True)
    if process.returncode:
        raise RuntimeError(process.stderr + process.stdout[-2000:])
    output = process.stdout
    metrics = ov.parse_validation_metrics(output)
    if (metrics['samples'] != 180000 or metrics['input'] != path.name or
            metrics['family'] != 'TFG_CARMA_SHARED_MODES_CSM' or
            not all(math.isfinite(metrics[m]) for m in METRICS)):
        raise ValueError(f'Incomplete replay on {path.name}')
    for variable, banner in [('W3D_ENGINE_RPM', 'ENGINE_VIBRATION rpm='),
                             ('W3D_IMU_LEVER_ARM_M', 'IMU_LEVER_ARM_RESULT '),
                             ('W3D_IMU_LEVER_ARM_FILTER', 'IMU_LEVER_ARM_FILTER placement=wrapper')]:
        if (variable in extra) != (banner in output):
            raise ValueError(f'Unexpected state of {variable} on {path.name}')
    return metrics


def conditions():
    yield 'baseline', 0, 'CG', 0, {}
    for rpm in (800, 1600, 2400, 3200):
        yield 'engine', rpm, 'CG', 0, {'W3D_ENGINE_RPM': str(rpm)}
    for axis, direction in AXES.items():
        for distance in (0.3, 1.0):
            for mode in ('unmodeled', 'wrapper', 'exact'):
                yield mode, 0, axis, distance, arm_env(mode, direction, distance, 15)


def run_cell(task):
    path, (mode, rpm, axis, distance, overrides) = task
    metrics = score(path, overrides)
    return dict(record=path.name, mode=mode, rpm=rpm, axis=axis, distance_m=distance,
                **{m: metrics[m] for m in METRICS})


def transition(data):
    low = ov.find_default_input(data, '1.500', '50.710')
    high = ov.find_default_input(data, '8.500', '202.839')
    columns, low_values = ov.read_wave_csv(low, ROUNDTRIP_DURATION_SEC)
    high_columns, high_values = ov.read_wave_csv(high, ROUNDTRIP_DURATION_SEC)
    if high_columns != columns:
        raise ValueError('Round-trip source columns differ')
    values = make_roundtrip_wave(ov, columns, low_values, high_values, 11,
                                 DIAGNOSTIC_HIGH_HEIGHT_M / DIAGNOSTIC_SOURCE_HEIGHT_M)
    segments = roundtrip_segments()
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / ROUNDTRIP_WAVE_NAME
        ov.write_wave_csv(path, columns, values)
        metrics = score(path, {'W3D_VALIDATION_SEGMENTS': ','.join(
            f'{name}:{start}:{end}' for name, start, end in segments)})
        return roundtrip_score_rows(metrics, segments), digest(path)


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, default=ROOT / '.sim-data/cases' / CASE)
    parser.add_argument('--output', type=Path, default=OUT / CASE)
    parser.add_argument('--jobs', type=int, default=2)
    args = parser.parse_args()
    files = [args.data_dir / name for name in sorted(REFERENCE_CSV_SHA256)]
    dataset = input_provenance(files)
    segments, surrogate_hash = transition(args.data_dir)
    cells = [(path, condition) for path in files for condition in conditions()]
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        rows = []
        for i, row in enumerate(pool.map(run_cell, cells), 1):
            rows.append(row)
            if i % 23 == 0:
                print(f'{CASE}: {i}/{len(cells)} scored', flush=True)
    args.output.mkdir(parents=True, exist_ok=True)
    write_csv(args.output / 'steady_seas.csv', rows)
    write_csv(args.output / 'roundtrip.csv', segments)
    sources = implementation_closure([BINARY.with_suffix('.cpp'), ROOT / 'src/util/W3dSimCommon.cpp'])
    sources.extend(ROOT / path for path in (
        'tests/kalman_tfg_shared_modes/Makefile', 'tests/common/HostReplayFlags.mk',
        'tools/tfg_shared_modes_degradation.py', 'tools/carma_lever_arm_study.py',
        'tools/ou_roundtrip_transition.py', 'tools/ou_validation.py', 'tools/underway_reference.py',
        'tools/sim_dataset.py'))
    build = environment_metadata()
    build['eigen'].pop('path', None)
    manifest = dict(dataset=dataset, build=build, filter='KalmanTfgCarmaSharedModes<float>',
                    wave_seed_roundtrip=11, imu_seed=101, initialization_seed=1009,
                    duration_s=1200, sample_rate_hz=200, score_last_s=900,
                    binary_sha256=digest(BINARY), roundtrip_input_sha256=surrogate_hash,
                    sources={str(p.relative_to(ROOT)): digest(p) for p in sources},
                    files={name: digest(args.output / name)
                           for name in ('steady_seas.csv', 'roundtrip.csv')})
    (args.output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    main()
