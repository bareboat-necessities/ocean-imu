#!/usr/bin/env python3
"""Stage and reproduce the source-bound, non-promoting low-wave experiments.

The recipes live with the measured results. This command never edits estimator
sources, installs a tuning default, or substitutes new acceptance limits.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / 'reports/results/encounter_scaling/low_wave_model_followup.json'


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work', type=Path, required=True, help='new empty staging directory')
    parser.add_argument('--moving-dir', type=Path, required=True)
    parser.add_argument('--stationary-dir', type=Path, required=True)
    parser.add_argument('--build', action='store_true', help='compile the staged diagnostic executables')
    parser.add_argument('--run', action='store_true', help='build and run all retained experiment specifications')
    args = parser.parse_args()
    report = json.loads(REPORT.read_text())
    work = args.work.resolve()
    if work.exists() and any(work.iterdir()):
        parser.error('--work must be empty, to exclude stale overlays and binaries')
    replacements = {'@REPO@': str(ROOT), '@WORK@': str(work),
                    '@MOVING@': str(args.moving_dir.resolve()),
                    '@STATIONARY@': str(args.stationary_dir.resolve())}

    def expand(text: str) -> str:
        for key, value in replacements.items():
            text = text.replace(key, value)
        return text

    for name, expected in report['reproduction']['source_sha256'].items():
        if digest((ROOT / name).read_bytes()) != expected:
            raise ValueError(f'experiment source changed: {name}')
    inputs = {(row['dataset'], row['input']): row['input_sha256']
              for experiment in report['experiments'] for row in experiment['rows']}
    for (dataset, name), expected in inputs.items():
        directory = args.moving_dir if dataset == 'moving' else args.stationary_dir
        if digest((directory / name).read_bytes()) != expected:
            raise ValueError(f'experiment input changed: {dataset}/{name}')
    staged = {}
    for name, recipe in report['reproduction']['files'].items():
        destination = work / name
        if not destination.resolve().is_relative_to(work):
            raise ValueError('recipe path escapes staging directory')
        if 'base_path' in recipe:
            lines = (ROOT / recipe['base_path']).read_text().splitlines(keepends=True)
            output, cursor = [], 0
            for edit in recipe['edits']:
                start, end = edit['start'], edit['end']
                if not cursor <= start <= end <= len(lines):
                    raise ValueError(f'invalid patch for {name}')
                output.extend(lines[cursor:start])
                output.append(edit['replacement'])
                cursor = end
            output.extend(lines[cursor:])
            text = ''.join(output)
        else:
            text = recipe['text']
        if digest(text.encode()) != recipe['template_sha256']:
            raise ValueError(f'experiment recipe hash mismatch: {name}')
        staged[destination] = expand(text)
    # Validate all source bindings and recipes before writing anything.
    for destination, text in staged.items():
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(text)
    print(f'Staged {len(staged)} source-bound diagnostic files in {work}')
    if args.build or args.run:
        for command in report['reproduction']['build_commands']:
            subprocess.run([expand(arg) for arg in command], cwd=ROOT, check=True)
    if args.run:
        for command in report['reproduction']['run_commands']:
            subprocess.run([expand(arg) for arg in command], cwd=ROOT, check=True)


if __name__ == '__main__':
    main()
