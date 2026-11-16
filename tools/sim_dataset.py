#!/usr/bin/env python3
"""Fetch and verify the versioned vessel-motion cases used by all replays."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import zipfile

RELEASE = 'v1.4.0'
REPOSITORY = 'bareboat-necessities/oceanography-waves-lib'
ROOT = Path(__file__).resolve().parents[1]
CASES = {
    'following': {
        'archive': 'sim-data-files-vessel-rao-28ft-5kn.zip',
        'sha256': '95f7af50842aec50d6f1f023f8453d7a5cf1e71414c89ad59c9de13a3e63a940',
        'mean_heading_deg': 0.0,
        'description': 'northbound, following seas',
    },
    'head': {
        'archive': 'sim-data-files-vessel-rao-28ft-5kn-south.zip',
        'sha256': '402d5c8366678b40fac015d8f1eb9dc26da71b13ddf6bed1e0f0651a7e64a3be',
        'mean_heading_deg': 180.0,
        'description': 'southbound, head / bow-quartering seas',
    },
}


def selected_case(case=None):
    case = case or os.environ.get('OCEAN_IMU_SIM_CASE', 'following')
    if case not in CASES:
        raise ValueError(f'unknown simulation case {case!r}; expected following or head')
    return case


# Existing study drivers import these aliases. One process/CI matrix job owns
# one case; never change the environment after importing a study module.
CASE = selected_case()
ARCHIVE = CASES[CASE]['archive']
SHA256 = CASES[CASE]['sha256']
URL = f'https://github.com/{REPOSITORY}/releases/download/{RELEASE}/{ARCHIVE}'
DESCRIPTION = (f'oceanography-waves-lib {RELEASE}, estimated 28 ft fin-keel sailboat RAO at CG, '
               f"5 kn, {CASES[CASE]['description']}, heading +/-3 deg with 20 s period")


def digest(path: Path) -> str:
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


REFERENCE_CSV_SHA256_BY_CASE = json.loads((ROOT / 'tools/sim_dataset_manifest.json').read_text())
REFERENCE_CSV_SHA256 = REFERENCE_CSV_SHA256_BY_CASE[CASE]


def input_provenance(paths, case=None) -> dict:
    case = selected_case(case)
    spec = CASES[case]
    records = {}
    for path in paths:
        path = Path(path)
        actual = digest(path)
        if REFERENCE_CSV_SHA256_BY_CASE[case].get(path.name) != actual:
            raise ValueError(f'{path}: input is not a pinned {RELEASE}/{case} vessel RAO reference')
        records[path.name] = actual
    return {'release': RELEASE, 'case': case, 'archive': spec['archive'], 'archive_sha256': spec['sha256'],
            'motion_model': 'VESSEL_RAO_28FT_5KN', 'speed_kn': 5.0,
            'mean_heading_deg': spec['mean_heading_deg'],
            'heading_amplitude_deg': 3.0, 'heading_period_s': 20.0,
            'reference_coordinates': 'wave displacement/velocity relative to declared cruise track; total acceleration',
            'input_sha256': records}


def study_output(name: str) -> Path:
    """Optional case-scoped root for legacy dataset drivers with fixed outputs."""
    return Path(os.environ.get('OCEAN_IMU_STUDY_ROOT', ROOT / 'reports/results')) / name


def require_complete_metrics(rows, returncode: int) -> None:
    """Diagnostic variants can fail accuracy bars, but not data completeness."""
    import math
    names = [Path(row.get('input', '')).name for row in rows]
    if returncode not in (0, 1) or len(names) != 8 or set(names) != set(REFERENCE_CSV_SHA256):
        raise ValueError(f'incomplete dataset replay: exit={returncode}, records={names}')
    for row in rows:
        for key in ('disp_z_rms_m', 'roll_rms_deg', 'pitch_rms_deg', 'yaw_rms_deg'):
            if not math.isfinite(float(row[key])):
                raise ValueError(f'nonfinite {key} in {row["input"]}')


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
        return f"oceanography-waves-lib {recorded['release']}, {recorded['motion_model']}, {recorded.get('case', 'historical')}"
    return 'archived stationary-vessel oceanography-waves-lib v1.2.3'


def verify(archive: Path, case=None) -> None:
    actual = digest(archive)
    expected = CASES[selected_case(case)]['sha256']
    if actual != expected:
        raise ValueError(f'{archive}: dataset SHA-256 mismatch: expected {expected}, got {actual}')


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


def fetch(archive: Path, case=None) -> None:
    if archive.exists():
        verify(archive, case)
        return
    archive.parent.mkdir(parents=True, exist_ok=True)
    temporary = archive.with_suffix(archive.suffix + '.download')
    spec = CASES[selected_case(case)]
    url = f"https://github.com/{REPOSITORY}/releases/download/{RELEASE}/{spec['archive']}"
    subprocess.run(['curl', '-fL', '--retry', '3', url, '-o', str(temporary)], check=True)
    verify(temporary, case)
    temporary.replace(archive)


def materialize(archive: Path, destinations: list[Path], case=None) -> None:
    """Share immutable inputs locally; CSV names stay compatible with every harness.

    Check all cached member hashes before reuse. A matching filename, an old
    surface ZIP, a stale marker, or a partial extraction cannot admit old data.
    """
    case = selected_case(case)
    verify(archive, case)
    cache = ROOT / '.sim-data' / CASES[case]['sha256']
    marker = cache / 'members.json'
    records = json.loads(marker.read_text()) if marker.exists() else {}
    with zipfile.ZipFile(archive) as bundle:
        members = bundle.namelist()
        if len(members) != 30 or len(set(members)) != 30 or any(Path(n).name != n or not n.endswith('.csv') for n in members):
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
    print(f'Verified {RELEASE}/{case}: {len(members)} inputs in {len(destinations)} directories')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=(*CASES, 'all'), default=CASE)
    parser.add_argument('--archive', type=Path)
    parser.add_argument('--check', action='store_true', help='verify existing archive without downloading')
    parser.add_argument('--dest', type=Path, nargs='+', default=[])
    args = parser.parse_args()
    if args.case == 'all' and args.archive:
        parser.error('--archive requires one --case')
    for case in CASES if args.case == 'all' else (args.case,):
        spec = CASES[case]
        archive = args.archive or ROOT / '.sim-data' / RELEASE / spec['archive']
        if args.check:
            verify(archive, case)
        else:
            fetch(archive, case)
        # With --case all each destination gets disjoint following/head dirs.
        destinations = [p / case if args.case == 'all' else p for p in args.dest]
        if destinations:
            materialize(archive, destinations, case)
        print(f"{spec['sha256']}  {archive}")


if __name__ == '__main__':
    main()
