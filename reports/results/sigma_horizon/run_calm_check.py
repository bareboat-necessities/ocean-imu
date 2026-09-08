#!/usr/bin/env python3
"""Paired wave-to-calm-to-wave check using the existing closed-kinematics generator."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.repo / 'tools'))
    import ou_validation as ov
    from ou_roundtrip_transition import make_roundtrip_wave
    source = ov.find_default_input(args.repo / 'tests/kalman_ou_iii', '0.270', '14.047')
    columns, low = ov.read_wave_csv(source, 1200.)
    calm = low.copy()
    for i, column in enumerate(columns):
        if column != "time":
            calm[:, i] = 0.
    generated = make_roundtrip_wave(ov, columns, low, calm, 11, 1.)
    segments = (('wave', 300., 400.), ('fade_to_calm', 400., 520.),
                ('calm_recovery', 520., 560.), ('calm_settled', 640., 800.),
                ('wave_return', 1040., 1200.))
    rows = []
    os.environ['OU_SIGMA_VAR_HORIZON_MAX_S'] = '35'
    with tempfile.TemporaryDirectory() as name:
        path = Path(name) / source.name
        ov.write_wave_csv(path, columns, generated)
        input_hash = hashlib.sha256(path.read_bytes()).hexdigest()
        for family in ('OU_II', 'OU_III'):
            for decay in ('1', '5'):
                os.environ['OU_SIGMA_STILL_DECAY_SEC'] = decay
                metrics, _, _ = ov.run_simulator(
                    family, path, 900., imu_seed=101, initialization_seed=1009,
                    tuning_mode='adaptive', aw_cov_sync='periodic', segments=segments,
                    write_timeseries=False)
                rows.append({'family': family, 'decay_sec': float(decay), 'metrics': metrics})
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / 'calm-check.json').write_text(json.dumps({
        'source_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=args.repo, text=True).strip(),
        'producer_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'generated_input_sha256': input_hash, 'source_input_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'source_diff_sha256': hashlib.sha256(subprocess.check_output(['git', 'diff', 'HEAD'], cwd=args.repo)).hexdigest(),
        'wave_imu_initialization_seeds': [11, 101, 1009], 'segments': segments,
        'high_motion': 'exact_zero', 'controls': {'OU_SIGMA_VAR_HORIZON_MAX_S': '35'},
        'rows': rows}, indent=2) + '\n')


if __name__ == '__main__':
    main()
