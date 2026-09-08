#!/usr/bin/env python3
"""Fetch and verify the single versioned vessel-motion dataset used by all replays."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile

RELEASE = 'v1.2.1'
ARCHIVE = 'sim-data-files-vessel-rao-28ft.zip'
SHA256 = '6d6eb92db78e97f1d2456c6b92387993bf23be14a9a6cf26943f0a22fc6fee60'
REPOSITORY = 'bareboat-necessities/oceanography-waves-lib'
URL = f'https://github.com/{REPOSITORY}/releases/download/{RELEASE}/{ARCHIVE}'
DESCRIPTION = f'oceanography-waves-lib {RELEASE}, estimated 28 ft fin-keel sailboat RAO at CG'
ROOT = Path(__file__).resolve().parents[1]


def digest(path: Path) -> str:
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


REFERENCE_CSV_SHA256 = {'wave_data_jonswap_H0.270_L14.047_A30.00_P60.00.csv': '6a8a27234cb8fe3b1ecd524d2c13a9ebd029d584652f4212023ba974bd218e03', 'wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.csv': 'cd8e8ea89123773e7d8dcf7ac8712d780528c6ff02436b0f40b4d1a6a1bb97e0', 'wave_data_jonswap_H4.000_L112.766_A30.00_P30.00.csv': '11e29a93ad1084b02af116311d194416567fa6bc20ac9569f5471f97b3d9f8bb', 'wave_data_jonswap_H8.500_L202.839_A-30.00_P72.00.csv': 'f147cba9d14b87d0e1af4d3e434196312d91295cd40b9a9d738e2919e60cdaa7', 'wave_data_pmstokes_H0.270_L14.047_A30.00_P60.00.csv': 'a1a4fff90e99a5be68f7767b2865db317e413c87df8e942440c6f6ca1520cf0c', 'wave_data_pmstokes_H1.500_L50.710_A-30.00_P120.00.csv': '3db45e87b465aed441637ba7d7d8b67c3944416a784ea9aef6ce57e9360e72e4', 'wave_data_pmstokes_H4.000_L112.766_A30.00_P30.00.csv': 'f7d56bcfe9f901727e3dd1e21e61144d666943f6a74e2c1b5f4cc13d0cf580f2', 'wave_data_pmstokes_H8.500_L202.839_A-30.00_P72.00.csv': '793c3f6fd08bef8c980da53f7b77fc2f4e60914fe1e652f47ed5a7996d9a9fa1'}


def input_provenance(paths) -> dict:
    records = {}
    for path in paths:
        path = Path(path)
        actual = digest(path)
        if REFERENCE_CSV_SHA256.get(path.name) != actual:
            raise ValueError(f'{path}: input is not a pinned v1.2.1 vessel RAO reference')
        records[path.name] = actual
    return {'release': RELEASE, 'archive': ARCHIVE, 'archive_sha256': SHA256,
            'motion_model': 'VESSEL_RAO_28FT', 'input_sha256': records}


def verify(archive: Path) -> None:
    actual = digest(archive)
    if actual != SHA256:
        raise ValueError(f'{archive}: dataset SHA-256 mismatch: expected {SHA256}, got {actual}')


def fetch(archive: Path) -> None:
    if archive.exists():
        verify(archive)
        return
    archive.parent.mkdir(parents=True, exist_ok=True)
    temporary = archive.with_suffix(archive.suffix + '.download')
    subprocess.run(['curl', '-fL', '--retry', '3', URL, '-o', str(temporary)], check=True)
    verify(temporary)
    temporary.replace(archive)


def materialize(archive: Path, destinations: list[Path]) -> None:
    """Share immutable inputs locally; CSV names stay compatible with every harness.

    Check all cached member hashes before reuse. A matching filename, an old
    surface ZIP, a stale marker, or a partial extraction cannot admit old data.
    """
    cache = ROOT / '.sim-data' / SHA256
    marker = cache / 'members.json'
    records = json.loads(marker.read_text()) if marker.exists() else {}
    with zipfile.ZipFile(archive) as bundle:
        members = bundle.namelist()
        if len(members) != 30 or any(Path(n).name != n or not n.endswith('.csv') for n in members):
            raise ValueError('unexpected dataset members')
        valid = set(records) == set(members) and all(
            (cache / n).is_file() and digest(cache / n) == records[n] for n in members)
        if not valid:
            cache.mkdir(parents=True, exist_ok=True)
            bundle.extractall(cache)
            records = {n: digest(cache / n) for n in members}
            marker.write_text(json.dumps(records, indent=2, sort_keys=True) + '\n')
    for destination in destinations:
        destination.mkdir(parents=True, exist_ok=True)
        for name in members:
            target = destination / name
            source = cache / name
            if target.is_symlink() and target.resolve() == source:
                continue
            if target.exists() or target.is_symlink():
                target.unlink()
            target.symlink_to(source)
    print(f'Verified {DESCRIPTION}: {len(members)} inputs in {len(destinations)} directories')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, default=ROOT / ARCHIVE)
    parser.add_argument('--check', action='store_true', help='verify existing archive without downloading')
    parser.add_argument('--dest', type=Path, nargs='+', default=[])
    args = parser.parse_args()
    if args.check:
        verify(args.archive)
    else:
        fetch(args.archive)
    if args.dest:
        materialize(args.archive, args.dest)
    print(f'{SHA256}  {args.archive}')


if __name__ == '__main__':
    main()
