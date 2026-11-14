#!/usr/bin/env python3
"""Attribute the stationary -> 5 kn replay change to its physical causes.

The v1.3.0 records differ from the stationary v1.2.3 records in two separable
ways: the declared steering track (5 kn, heading +/-3 deg over 120 s) and the
encounter-frequency change of the wave response.  This tool builds exact
counterfactual inputs, replays built executables on them, and reports the
estimator-independent physics of each input pair.

  inputs  paired input physics, no estimator (direction, Doppler, stillness,
          Wiener bound, CSV time-stamp rounding)
  arms    counterfactual records, plus links to the two measured records
  run     replay built *-sim executables on the arms; one CSV row per record
  table   factor tables (Markdown) from run CSVs

Arms.  Body specific force and angular rate are rebuilt from the world
kinematics of every counterfactual, so each arm is exactly self-consistent:
  S        v1.2.3 stationary record
  M        v1.3.0 5 kn record
  notrack  M with the declared track and heading removed (no steering)
  track    S with the same track and heading added (steering, no Doppler)
  stretch  S slowed by s = Tz(M)/Tz(S): the stationary spectral shape at the
           encounter period, same heave, accelerations scaled by 1/s^2
Truth and track are used only for inputs and scoring, never by an estimator.
The scorer takes reference motion from the input CSV in double precision and
evaluates the track at sample_index * dt, not at the rounded CSV time column.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
from scipy.signal import csd, lfilter, welch
from scipy.special import jv

# numpy < 2 (Ubuntu 24.04 CI) has only trapz.
trapezoid = getattr(np, 'trapezoid', None) or np.trapz

G = 9.80665
FS = 200.0
DT = 1.0 / FS
SPEED_MPS = 5 * 1852 / 3600
HEADING_AMPLITUDE_RAD = np.deg2rad(3.0)
HEADING_FREQUENCY_RADPS = 2 * np.pi / 120.0
WINDOW_SAMPLES = 180000  # trailing 900 s, as scored by the replay gates
# Simulated sensor of W3dSimCommon.h (process_wave_file_for_tracker).
ACC_WHITE_SIGMA = 1.51e-3 * G
ACC_BIAS_RW = 0.0005
# Shared StillnessAdapter (tuner/SeaStateFusionTunerCommon.h) and its input.
STILL_ENERGY_THRESHOLD_G2 = 8e-4
STILL_HOLD_S = 2.0
STILL_EMA_ALPHA = 0.05
TRACKER_INPUT_CUTOFF_HZ = 6.0
QUAT = ['q_wb_zu_w', 'q_wb_zu_x', 'q_wb_zu_y', 'q_wb_zu_z']  # world -> body
FAMILIES = {
    'ou2': 'kalman_ou_ii/kalman_ou_ii-sim',
    'ou3': 'kalman_ou_iii/kalman_ou_iii-sim',
    'tfg': 'kalman_tfg/kalman_tfg-sim',
    'carma': 'kalman_carma/carma3-sim',
    'csm': 'kalman_carma_shared_modes/kalman_carma_shared_modes-sim',
    'tcarma': 'kalman_tfg_carma/tfg_carma3-sim',
    'tcsm': 'kalman_tfg_shared_modes/kalman_tfg_shared_modes-sim',
}
ARMS = ('S', 'stretch', 'notrack', 'track', 'M')
UNDERWAY_ARMS = ('M', 'track')
BANDS_HZ = ((0.0, 0.1), (0.1, 0.2), (0.2, 0.35), (0.35, 10.0))
RECORD = re.compile(r'wave_data_(jonswap|pmstokes)_H([\d.]+)_L([\d.]+)_A([-\d.]+)_P[-\d.]+\.csv')


def track(times):
    """Declared v1.3.0 track: position, velocity, acceleration, heading, rate."""
    t = np.asarray(times, dtype=float)
    phase = HEADING_FREQUENCY_RADPS * t
    x = jv(0, HEADING_AMPLITUDE_RAD) * t
    y = np.zeros_like(t)
    for n in range(1, 13):
        c = 2 * jv(n, HEADING_AMPLITUDE_RAD) / (n * HEADING_FREQUENCY_RADPS)
        if n % 2 == 0:
            x = x + c * np.sin(n * phase)
        else:
            y = y + c * (1 - np.cos(n * phase))
    heading = HEADING_AMPLITUDE_RAD * np.sin(phase)
    rate = HEADING_AMPLITUDE_RAD * HEADING_FREQUENCY_RADPS * np.cos(phase)
    zero = np.zeros_like(t)
    position = SPEED_MPS * np.column_stack((x, y, zero))
    velocity = SPEED_MPS * np.column_stack((np.cos(heading), np.sin(heading), zero))
    acceleration = SPEED_MPS * np.column_stack((-rate * np.sin(heading), rate * np.cos(heading), zero))
    return position, velocity, acceleration, heading, rate


def sample_times(n: int) -> np.ndarray:
    """Exact sample times; the CSV time column has six significant digits."""
    return np.arange(n) * DT


def qmul(a, b):
    aw, ax, ay, az = a.T
    bw, bx, by, bz = b.T
    return np.stack([aw * bw - ax * bx - ay * by - az * bz,
                     aw * bx + ax * bw + ay * bz - az * by,
                     aw * by - ax * bz + ay * bw + az * bx,
                     aw * bz + ax * by - ay * bx + az * bw], axis=1)


def rotation(q):
    """Rotation matrices of unit quaternions [w, x, y, z]."""
    w, x, y, z = q.T
    return np.stack([
        np.stack([1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)], -1),
        np.stack([2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)], -1),
        np.stack([2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)], -1)], 1)


def euler_deg(q_wb):
    """Nautical roll, pitch, yaw (ZYX) of body -> world."""
    c = np.transpose(rotation(q_wb), (0, 2, 1))
    roll = np.degrees(np.arctan2(c[:, 2, 1], c[:, 2, 2]))
    pitch = np.degrees(-np.arcsin(np.clip(c[:, 2, 0], -1.0, 1.0)))
    yaw = np.degrees(np.arctan2(c[:, 1, 0], c[:, 0, 0]))
    return roll, pitch, yaw


def axes(prefix):
    return [f'{prefix}_{axis}' for axis in 'xyz']


def set_body_fields(frame: pd.DataFrame, q_wb: np.ndarray, gyro: np.ndarray) -> None:
    """Specific force f_b = C_wb (a + g z) and attitude columns from world truth."""
    acc = frame[axes('acc')].to_numpy()
    frame[['acc_bx', 'acc_by', 'acc_bz']] = np.einsum('nij,nj->ni', rotation(q_wb), acc + [0.0, 0.0, G])
    frame[['gyro_x', 'gyro_y', 'gyro_z']] = gyro
    frame[QUAT] = q_wb
    frame['roll_deg'], frame['pitch_deg'], frame['yaw_deg'] = euler_deg(q_wb)


def retrack(frame: pd.DataFrame, sign: int) -> pd.DataFrame:
    """Add (sign=+1) or remove (sign=-1) the declared track and heading.

    Heading is a world-z rotation: C_bw' = Rz(sign psi) C_bw, so
    q_wb' = q_wb (x) qz(-sign psi) and w_b' = w_b + sign psi_dot C_wb z.
    """
    out = frame.copy()
    position, velocity, acceleration, heading, rate = track(sample_times(len(frame)))
    for prefix, value in zip(('disp', 'vel', 'acc'), (position, velocity, acceleration)):
        out[axes(prefix)] = frame[axes(prefix)].to_numpy() + sign * value
    q = frame[QUAT].to_numpy()
    half = -sign * heading / 2
    q_new = qmul(q, np.column_stack((np.cos(half), 0 * half, 0 * half, np.sin(half))))
    body_z = rotation(q)[:, :, 2]
    gyro = frame[['gyro_x', 'gyro_y', 'gyro_z']].to_numpy() + sign * rate[:, None] * body_z
    set_body_fields(out, q_new, gyro)
    return out


def stretch(frame: pd.DataFrame, factor: float) -> pd.DataFrame:
    """Same record slowed by `factor`: x'(t) = x(t / factor), derivatives scaled."""
    t = sample_times(len(frame))
    source = t / factor
    out = frame.copy()
    for prefix, power in (('disp', 0), ('vel', 1), ('acc', 2)):
        for column in axes(prefix):
            out[column] = CubicSpline(t, frame[column].to_numpy())(source) / factor ** power
    q = np.column_stack([CubicSpline(t, frame[c].to_numpy())(source) for c in QUAT])
    q /= np.linalg.norm(q, axis=1, keepdims=True)
    gyro = np.column_stack([CubicSpline(t, frame[c].to_numpy())(source)
                            for c in ('gyro_x', 'gyro_y', 'gyro_z')]) / factor
    set_body_fields(out, q, gyro)
    return out


def spectrum(x, fs=FS, nperseg=2 ** 15):
    f, p = welch(np.asarray(x) - np.mean(x), fs=fs, nperseg=min(nperseg, len(x)))
    return f[1:], p[1:]


def zero_crossing_period(z) -> float:
    f, p = spectrum(z)
    return float(np.sqrt(trapezoid(p, f) / trapezoid(f * f * p, f)))


def stillness(acc_z, seed: int = 1, decays=(1.0, 5.0)) -> dict:
    """Replica of the shared stillness detector on noisy vertical acceleration."""
    a = np.asarray(acc_z) + np.random.default_rng(seed).normal(0.0, ACC_WHITE_SIGMA, len(acc_z))
    pole = np.exp(-2 * np.pi * TRACKER_INPUT_CUTOFF_HZ * DT)
    energy = lfilter([STILL_EMA_ALPHA], [1, STILL_EMA_ALPHA - 1],
                     (lfilter([1 - pole], [1, -pole], a) / G) ** 2)
    low = energy < STILL_ENERGY_THRESHOLD_G2
    # Continuous below-threshold time, capped at 60 s as in the detector.
    run = np.zeros(len(low))
    count = 0
    for i, flag in enumerate(low):
        count = count + 1 if flag else 0
        run[i] = min(count * DT, 60.0)
    still = run > STILL_HOLD_S
    window = slice(-WINDOW_SAMPLES, None)
    result = {'still_fraction': float(still[window].mean())}
    for decay in decays:
        attenuation = np.where(still, np.exp(-run / decay), 1.0)
        result[f'variance_factor_decay_{decay:g}s'] = float(attenuation[window].mean())
    return result


def wiener_bound(z) -> float:
    """Non-causal Wiener heave error from the simulated accelerometer (m RMS)."""
    f, s = spectrum(z, nperseg=2 ** 16)
    w = 2 * np.pi * f
    noise = (2 * ACC_WHITE_SIGMA ** 2 / FS + 2 * ACC_BIAS_RW ** 2 / w ** 2) / w ** 4
    return float(np.sqrt(trapezoid(s * noise / (s + noise), f)))


def time_stamp_artifact(frame: pd.DataFrame) -> float:
    """RMS x-track error from evaluating the track at the CSV time column."""
    rounded = track(frame['time'].to_numpy())[0][:, 0]
    exact = track(sample_times(len(frame)))[0][:, 0]
    return float(np.sqrt(np.mean((rounded - exact)[-WINDOW_SAMPLES:] ** 2)))


def body_residuals(frame: pd.DataFrame) -> dict:
    """Largest IMU-field departures from the record's own world kinematics."""
    q = frame[QUAT].to_numpy()
    q /= np.linalg.norm(q, axis=1, keepdims=True)
    force = np.einsum('nij,nj->ni', rotation(q), frame[axes('acc')].to_numpy() + [0.0, 0.0, G])
    rate = 2 * qmul(q, np.gradient(q, DT, axis=0) * np.array([1, -1, -1, -1]))[:, 1:]
    inner = slice(5, -5)
    return {'specific_force_residual_mps2': float(np.abs(frame[['acc_bx', 'acc_by', 'acc_bz']].to_numpy()
                                                          - force)[inner].max()),
            'angular_rate_residual_radps': float(np.abs(frame[['gyro_x', 'gyro_y', 'gyro_z']].to_numpy()
                                                         - rate)[inner].max())}


def read_record(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=np.float64)


def records(folder: Path, only: str = ''):
    return sorted(p for p in folder.glob('wave_data_*.csv') if RECORD.fullmatch(p.name) and only in p.name)


def input_physics(stationary: Path, moving: Path) -> dict:
    name = stationary.name
    _, height, length, angle = RECORD.fullmatch(name).groups()
    k = 2 * np.pi / float(length)
    omega = np.sqrt(G * k)
    s_frame, m_frame = read_record(stationary), read_record(moving)
    row = {'nominal_height_m': float(height),
           'deep_water_period_s': 2 * np.pi / omega,
           'encounter_period_s': 2 * np.pi / abs(omega - k * SPEED_MPS * np.cos(np.radians(float(angle))))}
    window = slice(-WINDOW_SAMPLES, None)
    for tag, frame, underway in (('stationary', s_frame, False), ('moving', m_frame, True)):
        z = frame['disp_z'].to_numpy()[window]
        vx = frame['vel_x'].to_numpy()
        if underway:
            vx = vx - track(sample_times(len(frame)))[1][:, 0]
        f, p = spectrum(z)
        row[tag] = {
            'heave_rms_m': float(z.std()),
            'vertical_acc_rms_mps2': float(frame['acc_z'].to_numpy()[window].std()),
            'heave_tz_s': zero_crossing_period(z),
            'heave_variance_below_0p1hz': float(trapezoid(p[f < 0.1], f[f < 0.1]) / trapezoid(p, f)),
            # Elevation is in phase with orbital velocity along propagation.
            'propagates_with_track': bool(np.corrcoef(z, vx[window])[0, 1] > 0),
            'wiener_bound_pct_nominal_h': 100 * wiener_bound(z) / float(height),
            **stillness(frame['acc_z'].to_numpy()),
            **body_residuals(frame),
        }
    row['csv_time_artifact_x_rms_m'] = time_stamp_artifact(m_frame)
    row['stretch_factor'] = row['moving']['heave_tz_s'] / row['stationary']['heave_tz_s']
    return row


def command_inputs(args) -> None:
    result = {}
    for path in records(args.stationary, args.only):
        result[path.name] = input_physics(path, args.moving / path.name)
        print(path.name, json.dumps(result[path.name]), flush=True)
    if args.json:
        args.json.write_text(json.dumps(result, indent=1, sort_keys=True) + '\n')


def command_arms(args) -> None:
    for path in records(args.stationary, args.only):
        s_frame, m_frame = read_record(path), read_record(args.moving / path.name)
        window = slice(-WINDOW_SAMPLES, None)
        factor = (zero_crossing_period(m_frame['disp_z'].to_numpy()[window]) /
                  zero_crossing_period(s_frame['disp_z'].to_numpy()[window]))
        built = {'notrack': retrack(m_frame, -1), 'track': retrack(s_frame, +1),
                 'stretch': stretch(s_frame, factor)}
        for arm, source in (('S', path), ('M', args.moving / path.name)):
            link = args.out / arm / path.name
            link.parent.mkdir(parents=True, exist_ok=True)
            if link.is_symlink() or link.exists():
                link.unlink()
            link.symlink_to(source.resolve())
        for arm, frame in built.items():
            (args.out / arm).mkdir(parents=True, exist_ok=True)
            frame.to_csv(args.out / arm / path.name, index=False, float_format='%.10g')
        print(f'{path.name}: stretch factor {factor:.4f}', flush=True)


SERIES = ['disp_est_x', 'disp_est_y', 'disp_est_z', 'sigma_a_applied', 'tau_applied', 'freq_tracker_hz']
HARNESS = ('roll_rms_deg', 'pitch_rms_deg', 'yaw_rms_deg', 'accel_bias_3d_rms_mps2', 'disp_z_pct_hs')


def score(record: Path, series: Path, underway: bool, nominal_height: float) -> dict:
    """Score the estimator time series against double-precision input truth."""
    truth = pd.read_csv(record, usecols=axes('disp'), dtype=np.float64).to_numpy()
    if underway:
        truth = truth - track(sample_times(len(truth)))[0]
    est = pd.read_csv(series, usecols=SERIES, dtype=np.float64)
    if len(est) != len(truth):
        raise ValueError(f'{series}: {len(est)} estimates for {len(truth)} samples')
    window = slice(-WINDOW_SAMPLES, None)
    error = est[axes('disp_est')].to_numpy()[window] - truth[window]
    row = {'z_rms_m': float(np.sqrt(np.mean(error[:, 2] ** 2))),
           'xy_rms_m': float(np.sqrt(np.mean(error[:, 0] ** 2) + np.mean(error[:, 1] ** 2)))}
    row['z_pct_h'] = 100 * row['z_rms_m'] / nominal_height
    # Heave error and energy-weighted gain by band (20 Hz is ample below 1 Hz).
    zt, ze = truth[window, 2][::10], est['disp_est_z'].to_numpy()[window][::10]
    f, pe = welch(ze - zt - np.mean(ze - zt), fs=FS / 10, nperseg=4096)
    _, pt = welch(zt - zt.mean(), fs=FS / 10, nperseg=4096)
    _, pte = csd(zt - zt.mean(), ze - ze.mean(), fs=FS / 10, nperseg=4096)
    for lo, hi in BANDS_HZ:
        band = (f >= lo) & (f < hi)
        tag = f'{lo:g}-{hi:g}hz'
        row[f'z_error_share_{tag}'] = float(trapezoid(pe[band], f[band]) / trapezoid(pe, f))
        row[f'z_truth_share_{tag}'] = float(trapezoid(pt[band], f[band]) / trapezoid(pt, f))
        weight = pt[band].sum()
        row[f'z_gain_{tag}'] = float(abs((pte[band]).sum() / weight)) if weight > 0 else float('nan')
    for column in ('sigma_a_applied', 'tau_applied', 'freq_tracker_hz'):
        row[f'{column}_median'] = float(np.median(est[column].to_numpy()[window]))
    return row


def run_one(job: dict) -> dict:
    env = {k: v for k, v in os.environ.items() if not k.startswith(('W3D_', 'SF_', 'OU_', 'HO_', 'NLO_', 'TFG_'))}
    # Stationary reference coordinates in every arm: the scorer below removes
    # the track itself, from exact sample times, for binaries with or
    # without the underway convention.
    env.update(W3D_COLLECT_ALL_GATES='1', W3D_WRITE_TIMESERIES='1', W3D_VALIDATION_WINDOW_SEC='900',
               W3D_REFERENCE_MOTION='stationary')
    env.update(job['env'])
    record = job['record']
    with tempfile.TemporaryDirectory(prefix='moving-diagnosis-') as tmp:
        work = Path(tmp)
        (work / record.name).symlink_to(record.resolve())
        result = subprocess.run([str(job['binary']), *job['extra'], '--input', record.name], cwd=work, env=env,
                                capture_output=True, text=True, check=False)
        metrics = {}
        for line in result.stdout.splitlines():
            if line.startswith('VALIDATION_METRICS '):
                metrics = dict(item.split('=', 1) for item in line.split()[1:])
        series = list(work.glob('w3d_*.csv'))
        if len(series) != 1 or not metrics:
            raise RuntimeError(f"{job['family']}/{job['arm']}/{record.name}: exit {result.returncode}, "
                               f'{len(series)} series\n{result.stdout[-2000:]}{result.stderr[-2000:]}')
        height = float(RECORD.fullmatch(record.name).group(2))
        row = {'family': job['family'], 'arm': job['arm'], 'record': record.name, 'label': job['label'],
               'exit_code': result.returncode,
               **{f'harness_{k}': float(metrics[k]) for k in HARNESS},
               **score(record, series[0], job['arm'] in UNDERWAY_ARMS, height)}
    print(row['label'], row['family'], row['arm'], record.name, f"{row['z_pct_h']:.4f}", flush=True)
    return row


def command_run(args) -> None:
    jobs = []
    for family in args.families.split(','):
        for arm in args.arms.split(','):
            for record in records(args.arms_dir / arm, args.only):
                jobs.append({'family': family, 'arm': arm, 'record': record, 'label': args.label,
                             'binary': args.tests / FAMILIES[family], 'extra': args.extra.split(),
                             'env': dict(item.split('=', 1) for item in args.env)})
    with ThreadPoolExecutor(args.jobs) as pool:
        rows = list(pool.map(run_one, jobs))
    frame = pd.DataFrame(rows)
    if args.out.exists() and args.append:
        frame = pd.concat([pd.read_csv(args.out), frame], ignore_index=True)
    frame.to_csv(args.out, index=False, float_format='%.6g')


def markdown(pivot: pd.DataFrame) -> str:
    names = list(pivot.index.names)
    lines = ['| ' + ' | '.join(names + list(pivot.columns)) + ' |',
             '|' + ' --- |' * len(names) + ' ---: |' * len(pivot.columns)]
    for key, values in pivot.iterrows():
        key = key if isinstance(key, tuple) else (key,)
        lines.append('| ' + ' | '.join([*map(str, key), *(f'{v:.4g}' for v in values)]) + ' |')
    return '\n'.join(lines)


def command_table(args) -> None:
    frame = pd.concat([pd.read_csv(p) for p in args.csv], ignore_index=True)
    frame['record'] = frame['record'].str.extract(r'wave_data_(\w+?_H[\d.]+)_')[0]
    for metric in args.metrics.split(','):
        pivot = frame.pivot_table(index=['label', 'family', 'record'], columns='arm', values=metric)
        pivot = pivot[[a for a in ARMS if a in pivot.columns]]
        print(f'\n#### {metric}\n')
        print(markdown(pivot))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='command', required=True)
    for name in ('inputs', 'arms'):
        p = sub.add_parser(name)
        p.add_argument('--stationary', type=Path, required=True, help='unpacked v1.2.3 records')
        p.add_argument('--moving', type=Path, required=True, help='unpacked v1.3.0 records')
        p.add_argument('--only', default='', help='record-name substring filter')
        if name == 'inputs':
            p.add_argument('--json', type=Path)
        else:
            p.add_argument('--out', type=Path, required=True)
    p = sub.add_parser('run')
    p.add_argument('--tests', type=Path, required=True, help='tests/ directory with built *-sim executables')
    p.add_argument('--arms-dir', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--families', default='ou2,ou3,tfg,carma,csm')
    p.add_argument('--arms', default=','.join(ARMS))
    p.add_argument('--only', default='')
    p.add_argument('--label', default='shipping')
    p.add_argument('--extra', default='', help='extra executable arguments, e.g. --no-noise')
    p.add_argument('--env', action='append', default=[], help='KEY=VALUE for diagnostic runs')
    p.add_argument('--jobs', type=int, default=2)
    p.add_argument('--append', action='store_true')
    p = sub.add_parser('table')
    p.add_argument('csv', type=Path, nargs='+')
    p.add_argument('--metrics', default='z_pct_h,xy_rms_m,harness_accel_bias_3d_rms_mps2')
    args = parser.parse_args()
    {'inputs': command_inputs, 'arms': command_arms, 'run': command_run, 'table': command_table}[args.command](args)


if __name__ == '__main__':
    main()
