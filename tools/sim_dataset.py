#!/usr/bin/env python3
"""Fetch and verify the single versioned vessel-motion dataset used by all replays."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile

RELEASE = 'v1.3.1'
ARCHIVE = 'sim-data-files-vessel-rao-28ft-5kn.zip'
SHA256 = 'c3505115275925e453e848248561803aeb6b9eadc4e03827c6657f4f1920a402'
REPOSITORY = 'bareboat-necessities/oceanography-waves-lib'
URL = f'https://github.com/{REPOSITORY}/releases/download/{RELEASE}/{ARCHIVE}'
DESCRIPTION = f'oceanography-waves-lib {RELEASE}, estimated 28 ft fin-keel sailboat RAO at CG, 5 kn, heading +/-3 deg with 120 s period'
ROOT = Path(__file__).resolve().parents[1]


def digest(path: Path) -> str:
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


REFERENCE_CSV_SHA256 = {'wave_data_jonswap_H0.270_L14.047_A30.00_P60.00.csv': '761ee9dd703c14363ceeed73b9027da6db6f941529dff33f49f6289237fbc0c8', 'wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.csv': '328bf5f24c6a579fc1d3b473a8df45c15a7a9f2146b85d03db14e92263113913', 'wave_data_jonswap_H4.000_L112.766_A30.00_P30.00.csv': '18d5ac472be4b88b6c61cfbde4fa7177b27d7b04c90262dcbf783e71ef0e45c7', 'wave_data_jonswap_H8.500_L202.839_A-30.00_P72.00.csv': 'dab59d4e72b3d75ba730f545813429a1d21bbd454cabd35fe56124d0e309cc4d', 'wave_data_pmstokes_H0.270_L14.047_A30.00_P60.00.csv': '0fb455a798fdb09a9d77bb75dc800e57d84ecdc823997f74adf8a0fcc55ce300', 'wave_data_pmstokes_H1.500_L50.710_A-30.00_P120.00.csv': '609339f115e4caa774bbe437bc429301e2d77b5603b5aaac785cb4aa62790970', 'wave_data_pmstokes_H4.000_L112.766_A30.00_P30.00.csv': '343ffbde75fe5e146018bd7790d6df4596e783a101b0b0d8346354ff2b1d56e4', 'wave_data_pmstokes_H8.500_L202.839_A-30.00_P72.00.csv': '7d471977f88b58c9a22ccc3e8c4e1e5a068367e0434e587e73307150c4bc3150'}


def input_provenance(paths) -> dict:
    records = {}
    for path in paths:
        path = Path(path)
        actual = digest(path)
        if REFERENCE_CSV_SHA256.get(path.name) != actual:
            raise ValueError(f'{path}: input is not a pinned {RELEASE} vessel RAO reference')
        records[path.name] = actual
    return {'release': RELEASE, 'archive': ARCHIVE, 'archive_sha256': SHA256,
            'motion_model': 'VESSEL_RAO_28FT_5KN', 'speed_kn': 5.0,
            'heading_amplitude_deg': 3.0, 'heading_period_s': 120.0,
            'reference_coordinates': 'wave displacement/velocity relative to declared cruise track; total acceleration',
            'input_sha256': records}


def bind_dataset_results(path: Path, data_dir: Path, binaries=()) -> None:
    """Invalidate append-only rows when input bytes or executables change.

    Filenames are identical across stationary and underway releases, so a
    filename/seed resume key alone would silently reuse the old experiment.
    Only this regenerable dataset CSV is replaced; synthetic results stay.
    """
    protocol = {'dataset': input_provenance(data_dir / n for n in REFERENCE_CSV_SHA256),
                'binaries': {str(Path(b).relative_to(ROOT)): digest(Path(b)) for b in binaries}}
    marker = path.with_suffix('.inputs.json')
    previous = json.loads(marker.read_text()) if marker.exists() else None
    if previous != protocol and path.exists():
        path.unlink()
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text(json.dumps(protocol, indent=2, sort_keys=True) + '\n')


def completed_dataset_runs(path: Path) -> set[tuple[str, str]]:
    """Resume only complete eight-record cells, never a partial invocation."""
    if not path.exists():
        return set()
    cells = {}
    with path.open() as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames
        rows = list(reader)
        for row in rows:
            cells.setdefault((row['family'], row['seed']), []).append(row['input'])
    complete = {key for key, names in cells.items()
                if len(names) == 8 and set(names) == set(REFERENCE_CSV_SHA256)}
    if complete != set(cells):
        with path.open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(row for row in rows if (row['family'], row['seed']) in complete)
    return complete


def dataset_result_description(path: Path) -> str:
    """Describe recorded rows by their receipt, never by today's fetch pin."""
    marker = path.with_suffix('.inputs.json')
    if marker.exists():
        recorded = json.loads(marker.read_text())['dataset']
        return f"oceanography-waves-lib {recorded['release']}, {recorded['motion_model']}"
    return 'archived stationary-vessel oceanography-waves-lib v1.2.3'


def verify(archive: Path) -> None:
    actual = digest(archive)
    if actual != SHA256:
        raise ValueError(f'{archive}: dataset SHA-256 mismatch: expected {SHA256}, got {actual}')


def verify_csv_clock(path: Path, *, first_index: int = 0, samples: int = 240000) -> None:
    """Require exact round trips of the integer-indexed 200 Hz clock."""
    count = 0
    with path.open() as stream:
        next(stream)  # header
        for count, line in enumerate(stream, 1):
            actual = float(line.partition(',')[0])
            expected = (first_index + count - 1) / 200
            if actual != expected:
                raise ValueError(f'{path}: invalid timestamp at sample {count}: {actual} != {expected}')
    if count != samples:
        raise ValueError(f'{path}: expected {samples} samples, got {count}')


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
    for name in members:
        if name.startswith('wave_data_') or name == 'wave_tracker_data.csv':
            fenton = '_fenton_' in name or name == 'wave_tracker_data.csv'
            samples = 5999 if name == 'wave_tracker_data.csv' else (239999 if fenton else 240000)
            verify_csv_clock(cache / name, first_index=2 if fenton else 0, samples=samples)
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
