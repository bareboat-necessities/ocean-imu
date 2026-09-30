"""Persistent tilt/BA compatibility orbit on the literal shipping loop.

Non-promoting finite diagnostic for the O2 obstruction of
`docs/ou3-proof-research-state.md` ("Persistent compatibility orbit: literal
reachability test"). Exact compatibility on a complete word requires the
applied nominal specific force to satisfy `f_hat=a_hat-g || b` at every
accelerometer epoch (zero-BA form), i.e. the nominal AW must carry
`P_b a_hat = P_b g`, whose norm is `g B_h/|B| >= g/5 = 1.96133 m/s^2` on the
admitted field geometry. That is the complement of Corollary A*. The audit asks
whether the unchanged OU-III execution can hold its nominal AW on that value,
or on its 16-s window mean, under admitted MARINE MOTION / IMU BIAS inputs.

The unchanged `SeaStateFusion_OU_III` runs once per history; the driver only
reads state after each sample, so no observer parity is needed. Admitted
histories reuse the exact rational envelope of the AW tracking audit along
the fixed unit direction `n=P_b g/|P_b g|` (plus an orthogonal channel along
`b`). DC probes and the inadmissible prefix of the continuation case are
linear-response probes and are recorded as not admitted. The float replay is
neither an enclosure nor a source-uniform bound.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile

import numpy as np

from .ag_readout_source_diagnostic import HEADER, REPO
from .aw_tracking_source_diagnostic import admissibility

DRIVER = REPO/'tools/stability/compat_orbit_source.cpp'
FUSION = Path('src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h')
EIGEN = '/usr/include/eigen3'
STEP = .005
G = 9.80665
B_NORM = 75.0
T0, RAMP, T_END, T_REC = 160, 20, 460, 300
WINDOW_S = 16
ROCKING = ('0.05', '0.1', '0.03', '0.07')    # roll/pitch amplitude [rad], frequency [Hz]
F_SYNC = F(200, 21)                          # literal 21-sample AW covariance sync cycle
BIAS_MAX = 0.22516660498395405               # constants.json imu_bias.B_a_mps2
PHASES = tuple(str(F(i, 10)) for i in range(10))
# (name, triangle amplitude, frequency): jerk-limited triangles phase-locked to
# the sync cycle, its harmonic and third subharmonic, and to the magnetic
# (8-sample) and default S-update cadences.
PHASE_SCAN = (('sync', '2.6', str(F_SYNC)), ('sync-harmonic', '1.3', str(2*F_SYNC)),
              ('sync-third', '7.5', str(F_SYNC/3)), ('mag-cadence', '0.97', '25'),
              ('mag-half', '1.95', '25/2'), ('S-cadence', '3.3', '22/3'))
# Tuner states: a swell along b sets tau (and T_S); the triangle shares jerk.
TUNER_SCAN = tuple((f'swell-{f2}', '2.5', f2, ph) for f2 in ('2/25', '3/25', '1/5')
                   for ph in ('3/10', '4/5'))
# Combined adversary: sync triangle, displacement-bounded slow wave along n and
# a maximal constant body accelerometer bias along n.
COMBINED = tuple((frac, ph) for frac in ('1/5', '7/25') for ph in ('3/10', '4/5'))
COMBINED_SLOW = ('3/10', '1/32')
DC_PROBES = ('1', '3')
CONTINUATION = ('2.6', 212.0, '4/5')        # DC prefix, switch time, triangle phase
FLIP_RAD = 1.0


def field(frac):
    bh = float(F(frac))*B_NORM
    return np.array([bh, 0.0, math.sqrt(B_NORM**2-bh**2)])


def directions(bv):
    b = bv/np.linalg.norm(bv)
    g = np.array([0.0, 0.0, G])
    v = g-(g@b)*b
    return v/np.linalg.norm(v), b, float(np.linalg.norm(v))


def _tri(x):
    x = x-np.floor(x)
    return np.where(x < .5, 1-4*x, 4*x-3)


def _tri_u(x):
    x = x-np.floor(x)
    return np.where(x < .5, x-2*x*x, 2*x*x-3*x+1)


def _tri_w(x):
    x = x-np.floor(x)
    w = np.where(x < .5, x*x/2-2*x**3/3, 1/24+(2*x**3/3-1.5*x*x+x)-(2/24-.375+.5))
    return w-1/48


def _envelope(t, t0=T0, tr=RAMP):
    x = np.clip((t-t0)/tr, 0, 1)
    inside = (x > 0) & (x < 1)
    e0 = x**3*(10-15*x+6*x*x)
    e1 = np.where(inside, 30*x*x*(1-x)**2/tr, 0)
    e2 = np.where(inside, 60*x*(1-x)*(1-2*x)/tr**2, 0)
    return e0, e1, e2


def admitted_scalar(t, a1, f1, phase, a2='0', f2='1'):
    """p=env d, d''=a1 tri(f1 t+phase)+a2 cos(2 pi f2 t) (as in the AW audit)."""
    a1, f1, phase, a2, f2 = (float(F(x)) for x in (a1, f1, phase, a2, f2))
    w2 = 2*math.pi*f2
    d = a1/f1**2*_tri_w(f1*t+phase)-a2/w2**2*np.cos(w2*t)
    v = a1/f1*_tri_u(f1*t+phase)+a2/w2*np.sin(w2*t)
    a = a1*_tri(f1*t+phase)+a2*np.cos(w2*t)
    e0, e1, e2 = _envelope(t)
    return e0*a+2*e1*v+e2*d


def smooth(t, t0, tr):
    x = np.clip((t-t0)/tr, 0, 1)
    return x**3*(10-15*x+6*x*x)


def compile_driver(tmp, eigen):
    exe = Path(tmp)/'compat.bin'
    subprocess.run(['g++', '-O2', '-std=c++20', '-I'+str(REPO/'src'), '-isystem', str(eigen),
                    str(DRIVER), '-o', str(exe)], check=True)
    return exe


def run(exe, tmp, tag, acc, bv, bias=(0.0, 0.0, 0.0), trec=T_REC):
    path = Path(tmp)/f'{tag}.txt'
    np.savetxt(path, acc, fmt='%.12g')
    out = subprocess.run([str(exe), str(path), *ROCKING, *(repr(float(x)) for x in bv), '8',
                          str(trec), *(repr(float(x)) for x in bias)],
                         check=True, capture_output=True, text=True).stdout.splitlines()
    path.unlink()
    rows = np.array([[float(x) for x in line.split()] for line in out[:-1]])
    return rows, int(out[-1].split()[1])


def sig(x, digits=4):
    return float(f'{x:.{digits}g}')


def ratios(rows, n, required):
    an = rows[:, 1:4]@n
    width = int(round(WINDOW_S/STEP))
    windows = np.convolve(an, np.ones(width)/width, 'valid')
    return {'long_mean_ratio': sig(an.mean()/required),
            'window_16s_max_abs_ratio': sig(np.abs(windows).max()/required),
            'pointwise_max_ratio': sig(an.max()/required)}


def diagnostic(eigen=EIGEN):
    t = np.arange(1, int(round(T_END/STEP))+1)*STEP
    record = {'qualification': 'OU3_COMPAT_ORBIT_SOURCE_FEASIBILITY_V1',
              'source_uniform_certificate': False, 'theorem_closed': False,
              'driver_sha256': hashlib.sha256(DRIVER.read_bytes()).hexdigest(),
              'header_sha256': hashlib.sha256((REPO/HEADER).read_bytes()).hexdigest(),
              'fusion_sha256': hashlib.sha256((REPO/FUSION).read_bytes()).hexdigest(),
              'rocking_roll_pitch': list(ROCKING), 'window_s': WINDOW_S,
              'record_interval_s': [T_REC, T_END]}
    with tempfile.TemporaryDirectory() as tmp:
        exe = compile_driver(tmp, eigen)
        bv = field('1/5')
        n, b, required = directions(bv)
        record['required_transverse_aw_mps2'] = sig(required, 6)
        jobs = []
        for name, a1, f1 in PHASE_SCAN:
            for ph in PHASES:
                jobs.append(('phase', name, a1, f1, ph, '0', '1', '1/5', False))
        for name, a1, f2, ph in TUNER_SCAN:
            jobs.append(('tuner', name, a1, str(F_SYNC), ph, '1', f2, '1/5', False))
        for frac, ph in COMBINED:
            jobs.append(('combined', f'{frac}-{ph}', '2.55', str(F_SYNC), ph,
                         *COMBINED_SLOW, frac, True))

        def admitted_job(i_job):
            i, (group, name, a1, f1, ph, a2, f2, frac, slow_on_n) = i_job
            bv = field(frac)
            n, b, required = directions(bv)
            if group == 'tuner':
                # Triangle along n, swell along b (orthogonal channel).
                acc = np.outer(admitted_scalar(t, a1, f1, ph), n) \
                    + np.outer(admitted_scalar(t, '0', '1', '0', a2, f2), b)
            else:
                acc = np.outer(admitted_scalar(t, a1, f1, ph, a2, f2), n)
            # Both channels are charged to one direction: conservative, since
            # |x n+z b|<=|x|+|z| bounds every derivative of the sum.
            motion = admissibility(a1, f1, a2, f2, '0', '1', ph)
            bias = BIAS_MAX*n if slow_on_n else np.zeros(3)
            rows, active = run(exe, tmp, f'a{i}', acc, bv, bias)
            row = {'group': group, 'arguments': [a1, f1, ph, a2, f2, frac],
                   'accelerometer_bias_body_mps2': float(np.linalg.norm(bias)),
                   'marine_motion': motion, 'a21_active_step': active,
                   'tau_aw_end': sig(rows[-1, 13]), 'pseudo_period_end_s': sig(rows[-1, 14])}
            row.update(ratios(rows, n, required))
            return f'{group}:{name}:{ph}', row

        with ThreadPoolExecutor(4) as pool:
            admitted = dict(pool.map(admitted_job, enumerate(jobs)))
        record['admitted'] = admitted

        probes = {}
        for d in DC_PROBES:
            acc = np.outer(float(d)*smooth(t, 200, 20), n)
            rows, _ = run(exe, tmp, f'dc{d}', acc, bv, trec=190)
            an = rows[:, 1:4]@n
            kernel = np.abs(rows[:, 10:13]@b)
            probes[d] = {'admitted': False, 'role': 'constant-drive linear-response probe',
                         'max_1s_ratio': sig(np.convolve(an, np.ones(200)/200, 'valid').max()/required),
                         'final_50s_ratio': sig(an[-10000:].mean()/required),
                         'ba_norm_end_mps2': sig(float(np.linalg.norm(rows[-1, 4:7]))),
                         'field_axis_rotation_max_rad': sig(kernel.max()),
                         'flipped_about_field': bool(kernel.max() > FLIP_RAD)}
        record['dc_probes'] = probes

        dc, toff, ph = CONTINUATION
        prefix = float(dc)*smooth(t, 200, 10)*(1-smooth(t, toff, .5))
        tail = smooth(t, toff, .5)*2.6*_tri(float(F_SYNC)*t+float(F(ph)))
        rows, _ = run(exe, tmp, 'cont', np.outer(prefix+tail, n), bv, trec=190)
        tt = rows[:, 0]*STEP
        an = rows[:, 1:4]@n
        i0 = int(np.argmin(abs(tt-toff)))
        below = np.nonzero(an[i0:] < .9*required)[0]
        record['continuation'] = {
            'admitted_prefix': False, 'switch_time_s': toff, 'dc_prefix_mps2': float(dc),
            'ratio_at_switch': sig(an[i0]/required),
            'time_to_leave_0p9_s': sig(below[0]*STEP) if len(below) else None,
            'final_100s_ratio': sig(an[-20000:].mean()/required)}

    worst = max(r['window_16s_max_abs_ratio'] for r in admitted.values())
    record['worst_admitted_window_ratio'] = worst
    record['worst_admitted_long_ratio'] = max(abs(r['long_mean_ratio']) for r in admitted.values())
    record['all_admitted'] = all(r['marine_motion']['admitted'] for r in admitted.values())
    record['persistent_compatibility_found'] = worst >= 1
    return record


def verify_diagnostic(record):
    """Internal consistency of a committed record; no replay required."""
    if record.get('qualification') != 'OU3_COMPAT_ORBIT_SOURCE_FEASIBILITY_V1':
        raise ValueError('wrong qualification')
    if record['source_uniform_certificate'] or record['theorem_closed']:
        raise ValueError('diagnostic cannot promote')
    for key, path in (('driver_sha256', DRIVER), ('header_sha256', REPO/HEADER),
                      ('fusion_sha256', REPO/FUSION)):
        if record[key] != hashlib.sha256(path.read_bytes()).hexdigest():
            raise ValueError(f'{key}: source changed since the committed audit')
    rows = record['admitted']
    expected = len(PHASE_SCAN)*len(PHASES)+len(TUNER_SCAN)+len(COMBINED)
    if len(rows) != expected:
        raise ValueError('admitted case set changed')
    for name, row in rows.items():
        a1, f1, ph, a2, f2, _ = row['arguments']
        if row['marine_motion'] != admissibility(a1, f1, a2, f2, '0', '1', ph):
            raise ValueError(f'{name}: admissibility record differs')
        if not row['marine_motion']['admitted']:
            raise ValueError(f'{name}: history not admitted')
        if row['a21_active_step'] < 0 or row['a21_active_step']*STEP >= T0:
            raise ValueError(f'{name}: A21 not active before the input onset')
        if row['accelerometer_bias_body_mps2'] > BIAS_MAX*(1+1e-9):
            raise ValueError(f'{name}: bias above B_a')
    if record['worst_admitted_window_ratio'] != max(r['window_16s_max_abs_ratio'] for r in rows.values()):
        raise ValueError('worst window ratio is not the case maximum')
    if record['persistent_compatibility_found'] != (record['worst_admitted_window_ratio'] >= 1):
        raise ValueError('compatibility flag inconsistent')
    if record['all_admitted'] is not True:
        raise ValueError('admitted flag inconsistent')
    for probe in record['dc_probes'].values():
        if probe['admitted']:
            raise ValueError('DC probes are not admitted histories')
    if record['continuation']['admitted_prefix']:
        raise ValueError('continuation prefix is not admitted')
    return True


def _close(a, b):
    if isinstance(a, dict):
        return isinstance(b, dict) and a.keys() == b.keys() and all(_close(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return isinstance(b, list) and len(a) == len(b) and all(map(_close, a, b))
    if isinstance(a, float) and isinstance(b, (int, float)) and not isinstance(b, bool):
        return abs(a-b) <= 2e-2*max(abs(a), abs(b))+2e-3
    return a == b


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--eigen', default=EIGEN)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--expect', type=Path)
    args = parser.parse_args()
    record = diagnostic(args.eigen)
    verify_diagnostic(record)
    text = json.dumps(record, indent=2, sort_keys=True)+'\n'
    args.output.write_text(text)
    if args.expect and not _close(json.loads(args.expect.read_text()), record):
        raise SystemExit('compatibility-orbit audit differs from the committed record')
    print(text)


if __name__ == '__main__':
    main()
