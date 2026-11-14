#!/usr/bin/env python3
"""Replay exact-clock zero-speed controls with the current moving-data binaries.

The supplied directory must contain the eight reference names from the v1.3.1
zero-speed 28 ft archive (SHA256 below). These are comparison inputs only;
the default regression dataset remains the 5 kn archive in sim_dataset.py.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile

from filter_sim_gates import REPLAYS, parse_pii, parse_w3d, parse_heave, profiles_for
from sim_dataset import REFERENCE_CSV_SHA256, ROOT, digest, verify_csv_clock

ARCHIVE_SHA256 = '5d74a266cf72a908db163b3e4ffc8a2fe4fffc0949b2c347ecc5b60bf5f68732'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--data-dir', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=2)
    parser.add_argument('--families', default=','.join(REPLAYS))
    parser.add_argument('--update-existing', action='store_true')
    args = parser.parse_args()
    if args.jobs < 1 or digest(args.archive) != ARCHIVE_SHA256:
        parser.error('positive jobs and the pinned v1.3.1 zero-speed archive required')
    names = args.families.split(',')
    if len(set(names)) != len(names) or not set(names) <= REPLAYS.keys():
        parser.error('families must be distinct replay names')
    # Verify both the archive identity and every extracted member before replay.
    import zipfile
    inputs = {}
    with zipfile.ZipFile(args.archive) as archive:
        for name in REFERENCE_CSV_SHA256:
            path = args.data_dir.resolve() / name
            inputs[name] = digest(path)
            if hashlib.sha256(archive.read(name)).hexdigest() != inputs[name]:
                raise ValueError(f'extracted stationary input differs: {name}')
            verify_csv_clock(path)
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(('W3D_', 'SF_', 'OU_', 'TFG_', 'PII_', 'NLO_', 'HO_'))}
    env.update(W3D_REFERENCE_MOTION='stationary', W3D_COLLECT_ALL_GATES='1',
               W3D_WRITE_TIMESERIES='0', W3D_VALIDATION_WINDOW_SEC='900')

    def run(name):
        folder, binary, _ = REPLAYS[name]
        executable = ROOT / 'tests' / folder / binary
        binary_hash = digest(executable)
        with tempfile.TemporaryDirectory(prefix=f'stationary-control-{name}-') as temp:
            work = Path(temp)
            texts, exits = [], []
            # AHRS returns after a failed gate; retain every record independently.
            groups = [[n] for n in sorted(inputs)] if name == 'ahrs' else [sorted(inputs)]
            for group in groups:
                for n in group:
                    (work / n).symlink_to(args.data_dir.resolve() / n)
                profiles = profiles_for('heave') if name == 'heave' else [('default', [])]
                heave = {}
                for profile, extra in profiles:
                    result = subprocess.run([str(executable), *extra], cwd=work, env=env,
                                            text=True, stdout=subprocess.PIPE,
                                            stderr=subprocess.STDOUT, timeout=1200)
                    if result.returncode not in (0, 1):
                        raise RuntimeError(result.stdout[-4000:])
                    texts.append(result.stdout)
                    exits.append(result.returncode)
                    if name == 'heave':
                        heave[profile] = dict(records=parse_heave(result.stdout, profile),
                            exit_code=result.returncode,
                            log_sha256=hashlib.sha256(result.stdout.encode()).hexdigest())
                for p in work.iterdir():
                    p.unlink()
            text = '\n'.join(texts)
        if digest(executable) != binary_hash:
            raise ValueError(f'{name}: binary changed during control')
        if name == 'heave':
            print(f'heave: {sum(len(p["records"]) for p in heave.values())} stationary method/record controls', flush=True)
            return name, dict(binary_sha256=binary_hash, profiles=heave)
        if name == 'pii':
            rows = parse_pii(text)
        elif name == 'nlo':
            values = re.findall(r'Z RMS raw \(%Hs\): (\S+)%', text)
            if len(values) != 8:
                raise ValueError('incomplete NLO control')
            rows = [dict(input=n, disp_z_pct_hs=float(v)) for n, v in zip(sorted(inputs), values)]
        elif name == 'ahrs':
            values = re.findall(r'Angles RMS \(deg\): Roll=(\S+) Pitch=(\S+) Yaw=(\S+)', text)
            if len(values) != 8:
                raise ValueError('incomplete AHRS control')
            rows = [dict(input=n, **dict(zip(('roll', 'pitch', 'yaw'), map(float, v))))
                    for n, v in zip(sorted(inputs), values)]
        else:
            rows = parse_w3d(text)
        print(f'{name}: {len(rows)} stationary control records', flush=True)
        return name, dict(binary_sha256=binary_hash, exit_codes=exits,
                          log_sha256=hashlib.sha256(text.encode()).hexdigest(),
                          records=rows)

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        families = dict(pool.map(run, names))
    path = ROOT / 'reports/results/filter_sim_gates/stationary_control.json'
    if args.update_existing:
        previous = json.loads(path.read_text())
        if previous['archive_sha256'] != ARCHIVE_SHA256 or previous['input_sha256'] != inputs:
            raise ValueError('cannot merge controls from different inputs')
        for name, control in previous['families'].items():
            if name in names:
                continue
            folder, binary, _ = REPLAYS[name]
            if digest(ROOT / 'tests' / folder / binary) != control['binary_sha256']:
                raise ValueError(f'{name}: retained control executable differs; remeasure this family')
            families[name] = control
    report = dict(role='same-executable exact-clock zero-speed diagnostic controls',
                  release='v1.3.1', archive='sim-data-files-vessel-rao-28ft.zip',
                  archive_sha256=ARCHIVE_SHA256, input_sha256=inputs,
                  protocol='1200 s at 200 Hz; trailing 900 s (AHRS 60 s); production defaults; all finite scores retained, including moving-data gate misses',
                  source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                  families=families)
    path.write_text(
        json.dumps(report, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
