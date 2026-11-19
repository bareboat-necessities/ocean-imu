#!/usr/bin/env python3
"""Reuse a complete case study only when its input, source and products verify."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
from pathlib import Path
import subprocess
import tempfile
import zipfile

import sim_dataset

ROOT = Path(__file__).resolve().parents[1]
STUDY = 'carma_shared_modes_validation'
# Includes every runtime implementation and build recipe, plus the complete
# Python study environment. Article-only changes do not invalidate scored runs.
REPLAY_PATHS = ['src', 'tools', 'tests/common', '.github/actions/replay-toolchain', 'Makefile'] + [f'tests/{name}' for name in (
    'kalman_carma_shared_modes', 'kalman_ou_ii', 'kalman_ou_iii', 'kalman_tfg',
    'pii_observer', 'nlo')]
# These gates only consume production code. They do not build, drive, or
# score the article's replays.
REPLAY_PATHS += [':(exclude)tools/quality/fusion_smoke.cpp', ':(exclude)tools/quality_gates.sh',
                 ':(exclude)tools/shared_modes_article_evidence.py']


class StaleSource(ValueError):
    """An intact bundle produced with an older replay implementation."""


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate(results, root=ROOT):
    receipt = json.loads((results/'case-receipt.json').read_text())
    protocol = json.loads((results/'protocol.json').read_text())
    summary = json.loads((results/'publication_summary.json').read_text())
    if receipt['case'] != sim_dataset.CASE or receipt['status'] != 'completed':
        raise ValueError('study receipt is incomplete or belongs to another case')
    if receipt['inputs'] != protocol['dataset']:
        raise ValueError('receipt and study input identities disagree')
    if protocol['dataset']['input_sha256'] != sim_dataset.REFERENCE_CSV_SHA256:
        raise ValueError('study input bytes differ from the selected release/case')
    seeds = protocol['seed_triplets']
    if len(seeds) != 10 or len({json.dumps(s, sort_keys=True) for s in seeds}) != 10 or summary['replays'] != 1038:
        raise ValueError('article requires the complete ten-seed study')
    if receipt['source_commit'] != protocol['source_revision'] or receipt['mode'] != 'full':
        raise ValueError('receipt and study source/protocol disagree')
    with (results/'scores.csv').open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 1038 or len({r['key'] for r in rows}) != 1038:
        raise ValueError('incomplete or duplicate replay inventory')
    if any(int(r['score_samples']) != 180000 for r in rows):
        raise ValueError('shortened scoring window')
    if not {'protocol.json', 'scores.csv', 'engine_summary.csv'} <= set(summary['files_sha256']):
        raise ValueError('incomplete product hash inventory')
    for path, digest in summary['files_sha256'].items():
        if sha(results/path) != digest:
            raise ValueError(f'changed study product: {path}')
    for path, digest in protocol['source_hashes'].items():
        if not (root/path).is_file() or sha(root/path) != digest:
            raise StaleSource(f'stale study source: {path}')
    # The narrow protocol hashes alone are insufficient: check runtime headers,
    # shared helpers and build recipes against the actual producing source tree.
    changed = subprocess.check_output(['git', 'diff', '--name-only',
        protocol['source_revision'], '--', *REPLAY_PATHS], cwd=root, text=True).splitlines()
    if changed:
        raise StaleSource(f'replay sources changed: {changed}')


def produce(results, jobs):
    repository = os.environ.get('GITHUB_REPOSITORY', 'bareboat-necessities/ocean-imu')
    name = f'simulation-study-{sim_dataset.RELEASE}-{sim_dataset.CASE}-{STUDY}'
    artifacts = json.loads(subprocess.check_output(['gh', 'api', '--paginate', '--slurp',
        f'repos/{repository}/actions/artifacts?name={name}&per_page=100'], text=True))
    for page in artifacts:
        for artifact in page['artifacts']:
            if artifact['expired']:
                continue
            revision = artifact['workflow_run']['head_sha']
            if subprocess.run(['git', 'merge-base', '--is-ancestor', revision, 'HEAD'], cwd=ROOT).returncode:
                continue
            with tempfile.TemporaryDirectory(prefix='shared-mode-article-') as tmp:
                archive = Path(tmp)/'study.zip'
                with archive.open('wb') as stream:
                    subprocess.run(['gh', 'api', f'repos/{repository}/actions/artifacts/{artifact["id"]}/zip'],
                                   stdout=stream, check=True)
                with zipfile.ZipFile(archive) as z:
                    z.extractall(Path(tmp)/'unpacked')
                candidate = Path(tmp)/'unpacked'/STUDY
                try:
                    validate(candidate)
                    producing_source = json.loads((candidate/'protocol.json').read_text())['source_revision']
                    if producing_source != revision:
                        raise ValueError('artifact workflow and replay source disagree')
                except (ValueError, OSError, KeyError) as error:
                    print(f'Rejecting artifact {artifact["id"]}: {error}', flush=True)
                    continue
                results.parent.mkdir(parents=True, exist_ok=True)
                shutil.copytree(candidate, results)
                print(f'Reused verified artifact {artifact["id"]} from {revision}', flush=True)
                return
    # No eligible current bundle exists: generate the actual complete study.
    # Never replace its source revision or relabel old rows as this execution.
    subprocess.run(['python3', '-u', 'tools/simulation_studies.py', STUDY,
                    '--jobs', str(jobs), '--output-root', str(results.parents[2])], cwd=ROOT, check=True)


def acquire(results, jobs):
    if results.exists():
        try:
            validate(results)
            return
        except StaleSource as error:
            print(f'Refreshing complete study with stale sources: {error}', flush=True)
        except (ValueError, OSError, KeyError) as error:
            raise RuntimeError(f'refusing to overwrite an existing invalid bundle: {error}') from error

    # Stage the download or full replay on the same filesystem. Keep the
    # previous bundle available if retrieval, generation, or validation fails.
    results.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=results.parent) as tmp:
        candidate = Path(tmp)/sim_dataset.CASE/'studies'/STUDY
        produce(candidate, jobs)
        validate(candidate)
        previous = Path(tmp)/'previous'
        if results.exists():
            results.rename(previous)
        try:
            candidate.rename(results)
        except OSError:
            if previous.exists():
                previous.rename(results)
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jobs', type=int, default=2)
    parser.add_argument('--check', type=Path)
    args = parser.parse_args()
    if args.check:
        validate(args.check)
    else:
        acquire(ROOT/'reports/results/sim_v1.4.0'/sim_dataset.CASE/'studies'/STUDY, args.jobs)


if __name__ == '__main__':
    main()
