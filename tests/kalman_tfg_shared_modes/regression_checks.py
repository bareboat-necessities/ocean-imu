#!/usr/bin/env python3
"""Exercise the real default-seed replay, protocol admission and failing gates."""
from __future__ import annotations
import argparse
import csv
import json
import os
import re
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from sim_dataset import REFERENCE_CSV_SHA256, input_provenance


def main() -> None:
    """Verify default replay scores and distinguish each rejected protocol."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, default=HERE / 'kalman_tfg_shared_modes-sim')
    parser.add_argument('--data-dir', type=Path, default=HERE)
    args = parser.parse_args()
    binary = args.binary.resolve()
    data_dir = args.data_dir.resolve()
    records = [data_dir / name for name in REFERENCE_CSV_SHA256]
    inputs = input_provenance(records)
    baseline = json.loads((ROOT / 'reports/results/filter_sim_gates/limits.json').read_text())
    assert baseline['inputs'] == inputs, 'locked limits use another input dataset'
    source = (ROOT / 'tests/common/CurrentReplayLimits.h').read_text().split('namespace interval {', 1)[1]
    source = re.search(r'W3dFailureLimits tcsm\{(.*?)\};', source, re.S)[1]
    compiled_limits = {k: float(v) for k, v in re.findall(r'\.(\w+)\s*=\s*([0-9.]+)f', source)}
    assert compiled_limits == baseline['families']['tcsm']['interval_limits'], 'C++ gates differ from the locked baseline'
    env = {k: v for k, v in os.environ.items() if not k.startswith('W3D_')}
    env['W3D_COLLECT_ALL_GATES'] = '1'
    env['W3D_WRITE_TIMESERIES'] = os.environ.get('W3D_WRITE_TIMESERIES', '1')
    result = subprocess.run([str(binary), '--regression'], cwd=data_dir, env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    print(result.stdout, end='', flush=True)
    assert result.returncode == 0, f'default replay failed: {result.returncode}'
    assert result.stdout.count('QUALITY_GATE: PASS=1 REASON=ok') == 8
    assert 'QUALITY_GATE: SKIPPED' not in result.stdout
    assert 'REGRESSION_PROTOCOL: GATED=1 ' in result.stdout
    assert 'REPLAY_COMPLETE: records=8 gated=1' in result.stdout
    lines = [line for line in result.stdout.splitlines() if line.startswith('VALIDATION_METRICS ')]
    parsed = [dict(pair.split('=', 1) for pair in line.split()[1:]) for line in lines]
    assert {r['input'] for r in parsed} == set(REFERENCE_CSV_SHA256)
    assert all(r['samples'] == '180000' and r['window_s'] == '900' for r in parsed)

    # Production defaults must meet their current, tighter numeric bars.
    current = subprocess.run([str(binary)], cwd=data_dir,
                             env={**env, 'W3D_WRITE_TIMESERIES': '0'}, text=True,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    assert current.returncode == 0, current.stdout[-3000:]
    assert 'SAMPLING=production_defaults' in current.stdout
    assert 'CSM_GYRO_INPUT convention=defaults schmidt_hold=0' in current.stdout
    assert current.stdout.count('QUALITY_GATE: PASS=1 REASON=ok') == 8
    assert 'REPLAY_COMPLETE: records=8 gated=1' in current.stdout
    print(current.stdout, end='', flush=True)
    print('TFG production defaults: all eight current numeric gates pass.', flush=True)

    with tempfile.TemporaryDirectory(prefix='tfg-shared-gates-') as tmp:
        work = Path(tmp)
        diagnostic = {**env, 'W3D_WRITE_TIMESERIES': '0'}

        def run(args=(), changes=None):
            """Run the strict replay in an isolated directory with controlled overrides."""
            return subprocess.run([str(binary), '--regression', *args], cwd=work,
                                  env={**diagnostic, **(changes or {})}, text=True,
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        r = run()
        assert r.returncode == 2 and 'no scored wave records found' in r.stdout, r.stdout[-1000:]
        protocol_error = '--regression requires default variant, noise, magnetometer, seeds and excitation'
        for args, changes, error in [
                (['--no-noise'], {}, protocol_error), (['--nomag'], {}, protocol_error),
                ([], {'W3D_SEED': '7'}, protocol_error),
                ([], {'W3D_CSM_VARIANT': 'CSM_OFF'}, protocol_error),
                ([], {'W3D_CSM_GYRO_INPUT': 'instantaneous'}, '--regression forbids W3D_CSM_GYRO_INPUT'),
                ([], {'W3D_CSM_GYRO_INPUT': 'interval'}, '--regression forbids W3D_CSM_GYRO_INPUT'),
                ([], {'W3D_ENGINE_RPM': '1600'}, '--regression forbids W3D_ENGINE_RPM'),
                ([], {'W3D_REFERENCE_MOTION': 'stationary'}, '--regression forbids W3D_REFERENCE_MOTION')]:
            r = run(['--input', str(records[0]), *args], changes)
            assert r.returncode == 2 and error in r.stdout, (args, changes, r.stdout[-1000:])

        # A missing or short file must not be reported as a passing regression.
        bad = work / records[0].name
        r = run(['--input', str(bad)])
        assert r.returncode != 0, r.stdout[-1000:]
        with records[0].open() as src, bad.open('w') as dst:
            dst.write(src.readline())
            dst.write(src.readline())
            dst.write(src.readline())
        r = run(['--input', str(bad)])
        assert r.returncode == 2 and 'full 1200 s' in r.stdout, r.stdout[-1000:]

        # Corrupt only displacement truth, keeping all delivered IMU samples.
        # This forces a large RMS error and must trip the actual numeric gates.
        with records[0].open(newline='') as src, bad.open('w', newline='') as dst:
            reader = csv.DictReader(src)
            writer = csv.DictWriter(dst, fieldnames=reader.fieldnames)
            writer.writeheader()
            for row in reader:
                row['disp_z'] = str(float(row['disp_z']) + 10.0)
                writer.writerow(row)
        r = run(['--input', str(bad)])
        assert r.returncode == 1, r.stdout[-1500:]
        assert 'ERROR: Z RMS above limit' in r.stdout and 'QUALITY_GATE: PASS=0' in r.stdout
        print('TFG shared-mode regression: eight pinned records pass; corrupted truth fails numeric gates; invalid, empty and truncated protocols fail closed.')


if __name__ == '__main__':
    main()
