#!/usr/bin/env python3
"""Paired steering attribution controls; never an acceptance run or filter aid.

The released input is compared with two explicitly different physical histories:
wave translation with the original yaw, and wave translation at fixed yaw. Both
remove the declared track's complete p/v/a chain. No truth enters a filter and no
quality limit is refitted. Dense inputs and outputs belong in --work-dir.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import inspect
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

from filter_sim_gates import REPLAYS
from sim_dataset import CASES, REFERENCE_CSV_SHA256_BY_CASE, ROOT, digest
from underway_reference import track
from ou_validation import _quaternion_rotate, _world_to_body_quaternion_from_euler_deg

sys.path.insert(0, str(ROOT / "plots"))
from vector_export import save_figure

FAMILIES = tuple(name for name, (_, _, full) in REPLAYS.items() if full)
VARIANTS = ('released', 'wave_yaw', 'wave_fixed_yaw')
METHOD_VERSION = 1


def indices(columns, prefix):
    return [columns.index(prefix + axis) for axis in 'xyz']


def control(columns, data, case, fixed_yaw=False):
    """Kinematically closed counterfactual, retaining the same wave response.

    This is an attribution experiment, not a regenerated hydrodynamic RAO.
    Original body rates are instantaneous, so remove the analytic yaw-rate
    contribution without re-differencing quaternions at interval endpoints.
    """
    result = data.copy()
    times = data[:, columns.index('time')]
    position, velocity, acceleration, _ = track(times, case)
    for prefix, values in zip(('disp_', 'vel_', 'acc_'), (position, velocity, acceleration)):
        result[:, indices(columns, prefix)] -= values
    if not fixed_yaw:
        q = data[:, [columns.index('q_wb_zu_' + a) for a in ('w', 'x', 'y', 'z')]]
        q = q / np.linalg.norm(q, axis=1)[:, None]
        result[:, indices(columns, 'acc_b')] -= _quaternion_rotate(q, acceleration)
        return result

    roll = np.deg2rad(data[:, columns.index('roll_deg')])
    pitch = np.deg2rad(data[:, columns.index('pitch_deg')])
    rate = np.deg2rad(3.0) * (2 * np.pi / 20) * np.cos(2 * np.pi * times / 20)
    yaw_terms = rate[:, None] * np.column_stack((
        -np.sin(pitch), np.sin(roll) * np.cos(pitch), np.cos(roll) * np.cos(pitch)))
    result[:, indices(columns, 'gyro_')] -= yaw_terms
    result[:, columns.index('yaw_deg')] = CASES[case]['mean_heading_deg']
    q = _world_to_body_quaternion_from_euler_deg(
        np.rad2deg(roll), np.rad2deg(pitch), result[:, columns.index('yaw_deg')])
    result[:, [columns.index('q_wb_zu_' + a) for a in ('w', 'x', 'y', 'z')]] = q
    force = result[:, indices(columns, 'acc_')].copy()
    force[:, 2] += 9.80665
    result[:, indices(columns, 'acc_b')] = _quaternion_rotate(q, force)
    return result


def harmonic_amplitude(times, signal, frequency=0.05):
    design = np.column_stack((np.ones_like(times), np.sin(2*np.pi*frequency*times),
                              np.cos(2*np.pi*frequency*times)))
    fit = np.linalg.lstsq(design, signal, rcond=None)[0]
    return float(np.hypot(fit[1], fit[2]))


def score(path):
    with path.open() as stream:
        columns = stream.readline().strip().split(',')
    names = ['time'] + [prefix+axis for prefix in
        ('disp_ref_', 'disp_est_', 'acc_bias_', 'acc_bias_est_') for axis in 'xyz']
    trace = np.loadtxt(path, delimiter=',', skiprows=1, usecols=[columns.index(n) for n in names])
    trace = trace[trace[:, 0] >= 300]
    if len(trace) != 180000:
        raise ValueError(f'{path}: incomplete trailing 900 s: {len(trace)}')
    tail = dict(zip(names, trace.T))
    result = {}
    for prefix, truth, estimate in (('disp', 'disp_ref_', 'disp_est_'),
                                    ('acc_bias', 'acc_bias_', 'acc_bias_est_')):
        for axis in 'xyz':
            error = tail[estimate+axis] - tail[truth+axis]
            if not np.isfinite(error).all():
                raise ValueError(f'{path}: nonfinite {prefix}/{axis}')
            result[f'{prefix}_{axis}_rms'] = float(np.sqrt(np.mean(error**2)))
            result[f'{prefix}_{axis}_20s_amplitude'] = harmonic_amplitude(tail['time'], error)
            result[f'{prefix}_{axis}_10s_amplitude'] = harmonic_amplitude(tail['time'], error, 0.1)
    return result


def publish(report, work, output):
    """Compact full-rate metrics and uniformly decimated diagnostic traces."""
    import csv
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = report['runs']
    table = output.with_suffix('.csv')
    fields = ['case', 'input', 'family', 'variant', *rows[0]['metrics']]
    with table.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        for row in rows:
            writer.writerow({**{k: row[k] for k in fields[:4]}, **row['metrics']})
    for case, record in dict.fromkeys((row['case'], row['input']) for row in rows):
        selected = [row for row in rows if row['case'] == case and row['input'] == record
                    and row['family'] == 'ou3']
        stem = f'{output.stem}-{case}-{Path(record).stem}'
        if len(selected) != 3 or {r['variant'] for r in selected} != set(VARIANTS):
            continue
        fig, axes = plt.subplots(2, 3, figsize=(12, 5.6), sharex=True, sharey='row')
        for col, variant in enumerate(VARIANTS):
            row = next(r for r in selected if r['variant'] == variant)
            folder = work/'runs'/case/variant/'ou3'/Path(row['input']).stem
            path = next(folder.glob('w3d*.csv'))
            with path.open() as stream:
                names = stream.readline().strip().split(',')
            keep = ['time', 'disp_ref_y', 'disp_est_y', 'acc_bias_y', 'acc_bias_est_y']
            values = np.loadtxt(path, delimiter=',', skiprows=1, usecols=[names.index(n) for n in keep])
            values = values[values[:, 0] >= 1140][::10]  # 200 Hz -> uniform 20 Hz, plot only
            compact = output.parent/f'{stem}-{variant}-tail.csv'
            np.savetxt(compact, values, delimiter=',', header=','.join(keep), comments='', fmt='%.9g')
            for ax, (true_index, estimate_index) in zip(axes[:, col], ((1, 2), (3, 4))):
                ax.plot(values[:, 0], values[:, true_index], color='#263445', lw=1.1, label='True')
                ax.plot(values[:, 0], values[:, estimate_index], color='#1c78a8', lw=1.2, label='Estimated')
                ax.grid(alpha=.2)
                ax.spines[['top', 'right']].set_visible(False)
            axes[0, col].set_title({'released': 'Turning at 5 kn', 'wave_yaw': 'Waves + original yaw',
                                    'wave_fixed_yaw': 'Waves + fixed yaw'}[variant])
            axes[1, col].set_xlabel('Time [s]')
        axes[0, 0].set_ylabel('World Y displacement [m]')
        axes[1, 0].set_ylabel('Filter-body Y bias [m/s²]')
        axes[0, 0].legend(frameon=False)
        height = float(selected[0]['input'].split('_')[3][1:])
        fig.suptitle(f'OU-III · {case} seas · {height:g} m incident waves')
        fig.tight_layout()
        figure_path = output.parent/f'{stem}.svg'
        with matplotlib.rc_context({'svg.hashsalt': 'heading-kinematics'}):
            save_figure(fig, figure_path, metadata={'Date': None})
        figure_path.write_text('\n'.join(line.rstrip() for line in figure_path.read_text().splitlines())+'\n')
        plt.close(fig)


def cached_trace_valid(run_dir, receipt):
    traces = list(run_dir.glob('w3d*.csv'))
    return (len(traces) == 1 and traces[0].is_file() and
            receipt.get('trace_sha256') == digest(traces[0]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--families', nargs='+', choices=FAMILIES, default=list(FAMILIES))
    parser.add_argument('--cases', nargs='+', choices=CASES, default=list(CASES))
    parser.add_argument('--variants', nargs='+', choices=VARIANTS, default=list(VARIANTS))
    parser.add_argument('--records', nargs='+', default=[
        'wave_data_jonswap_H0.270_L14.047_A30.00_P60.00.csv'])
    parser.add_argument('--jobs', type=int, default=2)
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    work = args.work_dir.resolve()
    recipe = hashlib.sha256((inspect.getsource(control) + inspect.getsource(_quaternion_rotate) +
                            inspect.getsource(_world_to_body_quaternion_from_euler_deg) +
                            (ROOT/'tools/underway_reference.py').read_text()).encode()).hexdigest()
    for family in args.families:
        folder, binary, _ = REPLAYS[family]
        if not (ROOT/'tests'/folder/binary).is_file():
            parser.error(f'build first: make -C tests/{folder} {binary}')
    tasks, inputs = [], {}
    for case in args.cases:
        for name in args.records:
            original = ROOT/'.sim-data'/CASES[case]['sha256']/name
            if digest(original) != REFERENCE_CSV_SHA256_BY_CASE[case][name]:
                raise ValueError(f'{original}: not the pinned released input')
            with original.open() as stream:
                columns = stream.readline().strip().split(',')
            data = np.loadtxt(original, delimiter=',', skiprows=1)
            for variant in args.variants:
                dest = work/'inputs'/case/variant/name
                dest.parent.mkdir(parents=True, exist_ok=True)
                marker = dest.with_suffix('.recipe.json')
                identity = {'original_sha256': digest(original), 'case': case, 'variant': variant,
                            'recipe_sha256': recipe}
                cached = json.loads(marker.read_text()) if marker.exists() else {}
                valid = dest.exists() and all(cached.get(k) == v for k, v in identity.items())
                valid = valid and cached.get('result_sha256') == digest(dest)
                if not valid:
                    if variant == 'released':
                        if not dest.is_symlink() or dest.resolve() != original.resolve():
                            dest.unlink(missing_ok=True)
                            dest.symlink_to(original)
                    else:
                        values = control(columns, data, case, variant == 'wave_fixed_yaw')
                        np.savetxt(dest, values, delimiter=',', header=','.join(columns),
                                   comments='', fmt='%.17g')
                    marker.write_text(json.dumps({**identity, 'result_sha256': digest(dest)})+'\n')
                inputs[f'{case}/{variant}/{name}'] = digest(dest)
                tasks.extend((case, variant, name, family, dest) for family in args.families)

    def run(task):
        case, variant, name, family, source = task
        folder, binary, _ = REPLAYS[family]
        executable = ROOT/'tests'/folder/binary
        binary_hash = digest(executable)
        run_dir = work/'runs'/case/variant/family/Path(name).stem
        run_dir.mkdir(parents=True, exist_ok=True)
        link = run_dir/name
        if not link.exists():
            link.symlink_to(source)
        receipt = run_dir/'receipt.json'
        protocol = {'case': case, 'variant': variant, 'family': family, 'input': name,
                    'input_sha256': digest(source), 'binary_sha256': binary_hash}
        if receipt.exists():
            previous = json.loads(receipt.read_text())
            if (previous.get('method_version', 1) == METHOD_VERSION
                    and all(previous.get(k) == v for k, v in protocol.items())
                    and cached_trace_valid(run_dir, previous)):
                return previous
        env = {k: v for k, v in os.environ.items()
               if not k.startswith(('W3D_', 'SF_', 'OU_', 'TFG_', 'CSM_', 'CARMA_'))}
        env.update(OCEAN_IMU_SIM_CASE=case, W3D_WRITE_TIMESERIES='1',
                   W3D_COLLECT_ALL_GATES='1', W3D_VALIDATION_WINDOW_SEC='900',
                   W3D_REFERENCE_MOTION='underway' if variant == 'released' else 'stationary')
        with (run_dir/'run.log').open('w') as log:
            result = subprocess.run([str(executable)], cwd=run_dir, env=env, stdout=log,
                                    stderr=subprocess.STDOUT, timeout=1200)
        if result.returncode not in (0, 1):
            raise RuntimeError(f'{run_dir}: replay failed with {result.returncode}')
        if digest(executable) != binary_hash:
            raise RuntimeError(f'{executable}: executable changed during replay')
        outputs = [p for p in run_dir.glob('*.csv') if p.name != name]
        if len(outputs) != 1:
            raise ValueError(f'{run_dir}: expected one dense output, got {outputs}')
        protocol.update(method_version=METHOD_VERSION, exit_code=result.returncode, metrics=score(outputs[0]),
                        trace_sha256=digest(outputs[0]), log_sha256=digest(run_dir/'run.log'))
        receipt.write_text(json.dumps(protocol, indent=2, allow_nan=False)+'\n')
        print(f"{case}/{variant}/{family}: Y RMS {protocol['metrics']['disp_y_rms']:.6f} m", flush=True)
        return protocol

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        rows = list(pool.map(run, tasks))
    report = {'source_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
              'role': 'attribution controls, not regression acceptance or achievable navigation aid',
              'score': 'trailing 900 s at 200 Hz; fixed production defaults and paired noise seeds',
              'frames': {'displacement': 'geographic generator NWU',
                         'bias': 'sensor body after (x,y,z) -> (y,x,-z); not world axes'},
              'units': {'disp': 'm', 'acc_bias': 'm/s^2'},
              'method_version': METHOD_VERSION, 'control_recipe_sha256': recipe,
              'inputs': inputs, 'runs': rows}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    publish(report, work, args.output)


if __name__ == '__main__':
    main()
