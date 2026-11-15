#!/usr/bin/env python3
"""Reproduce non-shipping wave-matched integral-observation experiments.

Recipes bind every overlay to the exact original and modified header bytes.
Preparation compiles into a separate directory; src/ and production defaults
are never edited. Gate failures remain in the report, including held-out draws.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
import tarfile

from encounter_scaling_study import verify_stationary
from filter_sim_gates import METRICS
import ou_evidence_provenance as provenance
from sim_dataset import ROOT, REFERENCE_CSV_SHA256, digest, input_provenance, verify_csv_clock

RECIPES = ROOT/'reports/results/encounter_scaling/spectral_pseudo_recipes.json'
FAMILIES = {
    'ou3': ('kalman_ou_iii', 'kalman_ou_iii-sim'),
    'tfg': ('kalman_tfg', 'kalman_tfg-sim'),
    'tfg_legacy': ('kalman_tfg', 'kalman_tfg-legacy-sim'),
}


def apply_recipe(recipe: dict, directory: Path, source_texts: dict | None = None) -> None:
    """Fail closed on stale source, overlapping edits, or incorrect output."""
    for name, patch in recipe['patches'].items():
        original = (ROOT/'src'/name).read_text() if source_texts is None else source_texts[name]
        if hashlib.sha256(original.encode()).hexdigest() != patch['base_sha256']:
            raise ValueError(f'recipe source changed: {name}')
        lines = original.splitlines(keepends=True)
        output, cursor = [], 0
        for edit in patch['edits']:
            start, end = edit['start'], edit['end']
            if not cursor <= start <= end <= len(lines):
                raise ValueError('invalid or overlapping recipe edits')
            output.extend(lines[cursor:start])
            output.append(edit['replacement'])
            cursor = end
        output.extend(lines[cursor:])
        data = ''.join(output).encode()
        if hashlib.sha256(data).hexdigest() != patch['overlay_sha256']:
            raise ValueError(f'recipe output mismatch: {name}')
        target = directory/name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
    for name, text in recipe['added'].items():
        target = directory/name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text)


def sources_for(family: str, overlay: Path, snapshot: Path | None = None) -> set[Path]:
    folder, name = FAMILIES[family]
    base = snapshot or ROOT
    previous = provenance.REPO_ROOT
    try:
        provenance.REPO_ROOT = base
        closure = set(provenance.implementation_closure([
            base/'tests'/folder/(name+'.cpp'), base/'src/util/W3dSimCommon.cpp',
        ]))
    finally:
        provenance.REPO_ROOT = previous
    return closure | set(overlay.rglob('*.h')) | {
        Path(__file__).resolve(), RECIPES, ROOT/'tools/encounter_scaling_study.py',
        ROOT/'tools/filter_sim_gates.py', ROOT/'tools/sim_dataset.py',
        ROOT/'tools/ou_evidence_provenance.py', ROOT/'tools/stationary_vessel_control.py',
    }


def prepare(family: str, variant: str, work: Path) -> None:
    recipes = json.loads(RECIPES.read_text())
    recipe = recipes['recipes'][variant+':'+family.replace('_legacy', '')]
    overlay = work/'include'
    if overlay.exists():
        raise ValueError('use a fresh build directory to exclude stale overlay headers')
    # Later shipping changes must not silently alter a historical baseline.
    # Archive the exact original implementation and test sources for both arms.
    snapshot = work/'source'
    snapshot.mkdir()
    archive = work/'source.tar'
    with archive.open('wb') as stream:
        subprocess.run(['git', 'archive', recipes['base_commit'], 'src', 'tests'],
                       cwd=ROOT, stdout=stream, check=True)
    with tarfile.open(archive) as stream:
        stream.extractall(snapshot, filter='data')
    archive.unlink()
    apply_recipe(recipe, overlay, {n: (snapshot/'src'/n).read_text() for n in recipe['patches']})
    folder, name = FAMILIES[family]
    eigen = ROOT/'third_party/eigen'
    if not (eigen/'Eigen/Dense').exists():
        eigen = Path('/usr/include/eigen3')
    binary = work/'replay'
    common = work/'W3dSimCommon.o'
    flags = ['g++', '-O3', '-std=c++20', '-Wall', '-Wextra', '-Wshadow', '-Wconversion',
             '-funroll-loops', '-fno-finite-math-only', '-march=native',
             '-I'+str(overlay), '-I'+str(snapshot/'src'), '-isystem', str(eigen)]
    commands = [flags+['-c', str(snapshot/'src/util/W3dSimCommon.cpp'), '-o', str(common)],
                flags+[str(snapshot/'tests'/folder/(name+'.cpp')), str(common), '-o', str(binary)]]
    sources = sources_for(family, overlay, snapshot)
    before = {str(p): digest(p) for p in sorted(sources)}
    # Serial builds are intentional: concurrent Eigen template builds exhaust RAM.
    for command in commands:
        subprocess.run(command, cwd=ROOT, check=True)
    if before != {str(p): digest(p) for p in sorted(sources)}:
        raise ValueError('source changed during compilation')
    record = dict(family=family, variant=variant, commands=commands, source_sha256=before,
                  source_commit=recipes['base_commit'],
                  binary_sha256=digest(binary), environment=provenance.environment_metadata())
    (work/'build.json').write_text(json.dumps(record, indent=2)+'\n')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--family', choices=FAMILIES, required=True)
    parser.add_argument('--variant', choices=('wave_matched', 'live_release', 'partial_schmidt', 'attitude_schmidt'), required=True)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--prepare', action='store_true', help='build a fresh overlay executable and stop')
    parser.add_argument('--moving-dir', type=Path)
    parser.add_argument('--stationary-dir', type=Path)
    parser.add_argument('--stationary-archive', type=Path)
    parser.add_argument('--seeds', default='default,7,23')
    parser.add_argument('--jobs', type=int, default=3)
    parser.add_argument('--engine-rpm', type=float)
    parser.add_argument('--regression', action='store_true', help='historical gyro profile, not the legacy TFG rS law')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    if args.variant == 'live_release' and args.family != 'ou3':
        parser.error('the recorded release refinement is OU-III only')
    if args.family == 'tfg_legacy' and args.regression:
        parser.error('the legacy-law harness has no alternate gyro profile')
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    if args.prepare:
        prepare(args.family, args.variant, work)
        return
    if not all((args.moving_dir, args.stationary_dir, args.stationary_archive, args.output)):
        parser.error('replay requires both input directories, zero-speed archive, and output')
    seeds = args.seeds.split(',')
    if args.jobs < 1 or len(seeds) != len(set(seeds)) or any(
            s != 'default' and (not s.isdecimal() or int(s) > 2**32-1) for s in seeds):
        parser.error('positive jobs and distinct unsigned seeds required')
    if args.engine_rpm is not None and (not math.isfinite(args.engine_rpm) or args.engine_rpm <= 0):
        parser.error('engine RPM must be finite and positive')
    binary, overlay = work/'replay', work/'include'
    build = json.loads((work/'build.json').read_text())
    if (build['family'], build['variant'], build['binary_sha256']) != (args.family, args.variant, digest(binary)):
        raise ValueError('build does not match the requested experiment')
    sources = sources_for(args.family, overlay, work/'source')
    initial = {str(p): digest(p) for p in sorted(sources)}
    if initial != build['source_sha256']:
        raise ValueError('source changed since the recorded build')
    directories = dict(moving=args.moving_dir.resolve(), stationary=args.stationary_dir.resolve())
    inputs = dict(moving=input_provenance(directories['moving']/n for n in REFERENCE_CSV_SHA256),
                  stationary=verify_stationary(args.stationary_archive, directories['stationary']))
    for directory in directories.values():
        for name in REFERENCE_CSV_SHA256:
            verify_csv_clock(directory/name)
    metadata = dict(role='experimental paired replay; not shipping/default acceptance or a stability receipt',
                    build=build, inputs=inputs, seeds=seeds, engine_rpm=args.engine_rpm,
                    regression=args.regression, protocol='1200 s at 200 Hz; score all 180000 samples in final 900 s')
    rows = []
    partial = args.output.with_suffix('.partial.json')
    if args.resume and partial.exists():
        saved = json.loads(partial.read_text())
        rows = saved.pop('rows')
        if saved != metadata:
            raise ValueError('checkpoint belongs to a different experiment')
    jobs = [(d, n, s, a) for d in directories for n in REFERENCE_CSV_SHA256
            for s in seeds for a in ('baseline', 'spectral_S')]
    done = {tuple(r[k] for k in ('dataset', 'input', 'seed', 'arm')) for r in rows}
    if len(done) != len(rows) or not done <= set(jobs):
        raise ValueError('duplicate or unexpected checkpoint rows')

    def run(job):
        dataset, name, seed, arm = job
        source = directories[dataset]/name
        env = {k: v for k, v in os.environ.items()
               if not k.startswith(('W3D_', 'SF_', 'OU_', 'TFG_', 'PROBE_'))}
        env.update(W3D_WRITE_TIMESERIES='0', W3D_COLLECT_ALL_GATES='1', W3D_VALIDATION_WINDOW_SEC='900')
        if dataset == 'stationary':
            env['W3D_REFERENCE_MOTION'] = 'stationary'
        if seed != 'default':
            env['W3D_SEED'] = seed
        if arm == 'spectral_S':
            env['PROBE_SPECTRAL_S'] = '1'
        if args.engine_rpm:
            env['W3D_ENGINE_RPM'] = str(args.engine_rpm)
        command = [str(binary), '--input', str(source)] + (['--regression'] if args.regression else [])
        with tempfile.TemporaryDirectory(prefix='pseudo-replay-', dir=work) as temporary:
            result = subprocess.run(command, cwd=temporary, env=env, capture_output=True, text=True, timeout=600)
        lines = [s for s in result.stdout.splitlines() if s.startswith('VALIDATION_METRICS ')]
        if len(lines) != 1 or result.returncode not in (0, 1):
            raise ValueError(f'invalid replay: {result.stdout[-2000:]}{result.stderr[-1000:]}')
        fields = dict(v.split('=', 1) for v in lines[0].split()[1:])
        if tuple(fields[k] for k in ('input', 'start_s', 'window_s', 'samples')) != (name, '300', '900', '180000'):
            raise ValueError('wrong input or incomplete scoring window')
        metrics = {k: float(fields[k]) for k in METRICS}
        if not all(math.isfinite(v) for v in metrics.values()):
            raise ValueError('nonfinite score')
        print(dataset, name, seed, arm, metrics['disp_z_pct_hs'], result.returncode, flush=True)
        return dict(dataset=dataset, input=name, seed=seed, arm=arm, metrics=metrics,
                    input_sha256=digest(source), exit_code=result.returncode,
                    gate_failures=[s for s in (result.stdout+'\n'+result.stderr).splitlines()
                                   if s.startswith(('ERROR:', 'QUALITY_GATE: PASS=0'))])

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(args.jobs) as pool:
        for future in as_completed([pool.submit(run, job) for job in jobs if job not in done]):
            rows.append(future.result())
            temporary = partial.with_suffix('.tmp')
            temporary.write_text(json.dumps(dict(**metadata, rows=rows), indent=2, allow_nan=False)+'\n')
            temporary.replace(partial)
    if initial != {str(p): digest(p) for p in sorted(sources)} or digest(binary) != build['binary_sha256']:
        raise ValueError('source or executable changed during replay')
    rows.sort(key=lambda r: tuple(r[k] for k in ('dataset', 'input', 'seed', 'arm')))
    args.output.write_text(json.dumps(dict(**metadata, rows=rows), indent=2, allow_nan=False)+'\n')
    partial.unlink()


if __name__ == '__main__':
    main()
