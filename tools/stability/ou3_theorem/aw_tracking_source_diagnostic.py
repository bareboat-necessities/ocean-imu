"""Carried AW tracking audit on admitted jerk-limited MOVING histories.

Non-promoting finite diagnostic. A temporary copy of the shipping header gets
read-only taps (the AW mean immediately before each accelerometer
correction, and the prediction/reset/row event stream of the world-frame
array); an untapped control compiled from the unchanged header must end in
the identical state. Every history is checked against MARINE MOTION with
exact rational envelope bounds before it is run. The float replay audits one
execution; it is neither an enclosure nor an all-time service certificate.

It answers two questions of `docs/ou3-world-frame-rows.md`, section 8:
the pointwise physical AW tracking premise of Corollary A (refuted on these
admitted histories) and the signed window means used by Corollary A*
(`aw_tracking.py`), which carry a large margin on the same executions.
It also rebuilds the literal world-frame six-column array (Lemma W, with the
actual resets) on a 100-s word and compares it with the injection-free array
of Theorem G0 and with G0's floor from the word's own nominal premises.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile

from .ag_readout_source_diagnostic import HEADER, REPO

DRIVER = REPO/'tools/stability/aw_tracking_source.cpp'
TAP_ANCHOR = '''    last_acc_diag_.accepted = false;

    // Reject invalid corrections'''
TAP = '''    last_acc_diag_.accepted = false;
    trace_aw_pre = xext.template segment<3>(OFF_AW).template cast<double>();
    if (aw_event_rec) std::printf("A %.17g %.17g %.17g\\n", trace_aw_pre.x(), trace_aw_pre.y(), trace_aw_pre.z());

    // Reject invalid corrections'''
EVENT_TAPS = (
    ('''    last_gyr_bias_corrected = gyr - gyro_bias;
''', '''    last_gyr_bias_corrected = gyr - gyro_bias;
    if (aw_event_rec) std::printf("P %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\\n",
        (double)qref.w(), (double)qref.x(), (double)qref.y(), (double)qref.z(),
        (double)last_gyr_bias_corrected.x(), (double)last_gyr_bias_corrected.y(),
        (double)last_gyr_bias_corrected.z(), (double)Ts);
'''),
    ('''    ocean_imu::kalman::ou_detail::apply_left_error_reset<T, NX>(Pext, dtheta_injected);
}''', '''    if (aw_event_rec) std::printf("R %.17g %.17g %.17g %.17g %.17g %.17g %.17g\\n",
        (double)dtheta_injected.x(), (double)dtheta_injected.y(), (double)dtheta_injected.z(),
        (double)qref.w(), (double)qref.x(), (double)qref.y(), (double)qref.z());
    ocean_imu::kalman::ou_detail::apply_left_error_reset<T, NX>(Pext, dtheta_injected);
}'''),
)
MAG_SPAN = ('::measurement_update_mag_only(', '::accelerometer_measurement_func(')
MAG_ANCHOR = '    xext.noalias() += K * r;'
FIELD = (F(21), F(0), F(72))          # uT, |B|=75, horizontal fraction 7/25
PI_UPPER, PI_LOWER = F(355, 113), F(333, 106)
LIMITS = {'A': F('8.8'), 'V': F('5.5'), 'P': F('8.1'), 'J': F(100)}
# (name, A1, f1, A2, f2, A3, f3, phase): horizontal triangle (phase in
# cycles) and cosine along e_x with a C2 onset at T0 over RAMP, vertical
# cosine present from construction. The last profile phase-locks the triangle
# to the literal 21-sample AW covariance sync cycle (sync when elapsed
# time_ exceeds ADAPT_EVERY_SECS=.1 s in float32 steps of .005 s).
PROFILES = (
    ('horizontal-triangle-at-jerk-limit', '8.7', '2.8', '0', '1', '0', '1', '0'),
    ('horizontal-triangle-with-vertical-swell', '5', '4.75', '0', '1', '3', '0.1', '0'),
    ('triangle-at-adaptation-cadence', '2.5', '9.5', '0', '1', '2.8', '0.1', '0'),
    ('collinear-witness-motion', '0', '1', '-2.636', '1', '0.7688', '1', '0'),
    ('large-swell-with-chop', '4', '5.9', '0', '1', '4', '0.2', '0'),
    ('sync-locked-rectification', '2.6', '9.523809523809524', '0', '1', '0', '1', '0.4'),
)
T0, RAMP, T_END = 160, 20, 460
WINDOW_S, SLIDE_S, STEP = 16, F(1, 2), F(1, 200)
G = F('9.80665')
CORR_A_PHYSICAL_THRESHOLD = F(112383, 100000)   # signed-mean form of Corollary A
CORR_A_STAR_THRESHOLD = G/5                      # nominal mean, sigma_w=1/5


def admissibility(a1, f1, a2, f2, a3, f3, phase='0', ramp=RAMP):
    """Exact rational MARINE MOTION envelope of the supplied history.

    p_x=env d, env the quintic smoothstep over `ramp` (derivative maxima
    15/(8T), 10/(sqrt3 T^2)<29/(5T^2), 60/T^3); d'' a triangle of amplitude
    a1 (|d'|<=a1/(8f), |d|<=a1/(48 f^2)) plus a cosine; p_z a cosine.
    """
    a1, f1, a2, f2, a3, f3 = (abs(F(x)) for x in (a1, f1, a2, f2, a3, f3))
    t = F(ramp)
    e1, e2, e3 = F(15, 8)/t, F(29, 5)/t**2, F(60)/t**3
    w2l, w2u = 2*PI_LOWER*f2, 2*PI_UPPER*f2
    w3l, w3u = 2*PI_LOWER*f3, 2*PI_UPPER*f3
    disp = a1/(48*f1**2)+a2/w2l**2
    vel = a1/(8*f1)+a2/w2l
    acc_amp = a1+a2
    jerk_amp = 4*a1*f1+a2*w2u
    acc_x = acc_amp+2*e1*vel+e2*disp
    jerk_x = jerk_amp+3*e1*acc_amp+3*e2*vel+e3*disp
    vel_x = vel+e1*disp
    bounds = {'A': (acc_x, a3), 'J': (jerk_x, a3*w3u), 'V': (vel_x, a3/w3l),
              'P': (disp, a3/w3l**2)}
    squared = {k: x*x+z*z for k, (x, z) in bounds.items()}
    return {'squared_upper': {k: str(v) for k, v in squared.items()},
            'admitted': all(squared[k] <= LIMITS[k]**2 for k in LIMITS),
            'jerk_margin_squared': str(LIMITS['J']**2-squared['J']),
            'displacement_dc': '0', 'bounded_primitive': True}


def _compile(tmp):
    header = (REPO/HEADER).read_text()
    if header.count(TAP_ANCHOR) != 1:
        raise ValueError('shipping accelerometer anchor changed')
    header = header.replace(TAP_ANCHOR, TAP)
    for old, new in EVENT_TAPS:
        if header.count(old) != 1:
            raise ValueError('shipping event anchor changed')
        header = header.replace(old, new)
    a, b = header.index(MAG_SPAN[0]), header.index(MAG_SPAN[1])
    if header[a:b].count(MAG_ANCHOR) != 1:
        raise ValueError('shipping magnetic anchor changed')
    header = header[:a]+header[a:b].replace(
        MAG_ANCHOR, '    if (aw_event_rec) std::printf("M\\n");\n'+MAG_ANCHOR)+header[b:]
    tapped = Path(tmp)/'observed/kalman_ou_iii'
    tapped.mkdir(parents=True)
    (tapped/'Kalman3D_Wave_OU_III.h').write_text(header)
    out = {}
    for name, extra in (('observed', ['-DAW_TRACKING_TAP', '-I'+str(Path(tmp)/'observed')]),
                        ('control', [])):
        exe = Path(tmp)/(name+'.bin')
        subprocess.run(['g++', '-O2', '-std=c++20', *extra, '-I'+str(REPO/'src'),
                        '-I'+EIGEN, str(DRIVER), '-o', str(exe)], check=True)
        out[name] = exe
    return out


EIGEN = '/usr/include/eigen3'


def run_profile(exes, profile):
    name, *args = profile
    argv = [*args[:6], str(T0), str(RAMP), str(T_END), str(T0+RAMP), args[6]]
    observed = subprocess.run([str(exes['observed']), *argv], check=True,
                              capture_output=True, text=True).stdout.splitlines()
    control = subprocess.run([str(exes['control']), *argv], check=True,
                             capture_output=True, text=True).stdout.splitlines()
    if observed[-1] != control[-1]:
        raise ValueError(f'{name}: observer changed the carried execution')
    rows = [[float(x) for x in line.split()] for line in observed[:-1] if line[0].isdigit()]
    events = [line for line in observed[:-1] if not line[0].isdigit()]
    active = int(observed[-1].split()[1])
    return name, args, rows, active, events


def analyze(rows, active_step, t_start=T0+RAMP):
    """Pointwise and signed trapezoidal 16-s window metrics (float64)."""
    b = [float(x)/75 for x in FIELD]
    n = int(WINDOW_S/STEP)
    slide = int(SLIDE_S/STEP)
    if rows[0][0]*float(STEP) < t_start-1e-9 or active_step*float(STEP) >= t_start:
        raise ValueError('record must start after the onset ramp and inside A21')
    err = [(r[3]-r[1], r[4], r[5]-r[2]) for r in rows]
    hat = [(r[3], r[4], r[5]) for r in rows]
    sup_err = max(math.sqrt(sum(c*c for c in e)) for e in err)
    g = float(G)
    force = [math.sqrt(h[0]**2+h[1]**2+(h[2]-g)**2)/g for h in hat]
    worst_err = worst_hat = worst_perp = worst_l1 = 0.0
    for s in range(0, len(rows)-n, slide):
        me, mh = [0.0]*3, [0.0]*3
        for i in range(n+1):
            w = (.5 if i in (0, n) else 1.0)/n
            for c in range(3):
                me[c] += w*err[s+i][c]
                mh[c] += w*hat[s+i][c]
        perp = (mh[1]*b[2]-mh[2]*b[1], mh[2]*b[0]-mh[0]*b[2], mh[0]*b[1]-mh[1]*b[0])
        worst_err = max(worst_err, math.sqrt(sum(c*c for c in me)))
        worst_hat = max(worst_hat, math.sqrt(sum(c*c for c in mh)))
        worst_perp = max(worst_perp, math.sqrt(sum(c*c for c in perp)))
        worst_l1 = max(worst_l1, sum((.5 if i in (0, n) else 1.0)/n*force[s+i] for i in range(n+1)))
    fmt = lambda x: float(f'{x:.6g}')
    return {
        'A21_active_step': active_step,
        'recorded_accelerometer_rows': len(rows),
        'pointwise_aw_error_sup_mps2': fmt(sup_err),
        'pointwise_ratio_to_corollary_A': fmt(sup_err/float(CORR_A_PHYSICAL_THRESHOLD)),
        'signed_16s_mean_error_max_mps2': fmt(worst_err),
        'signed_mean_ratio_to_1_12383': fmt(worst_err/float(CORR_A_PHYSICAL_THRESHOLD)),
        'nominal_16s_mean_max_mps2': fmt(worst_hat),
        'nominal_transverse_mean_max_mps2': fmt(worst_perp),
        'corollary_A_star_ratio': fmt(worst_perp/float(CORR_A_STAR_THRESHOLD)),
        'nominal_force_L1_max': fmt(worst_l1),
        'attitude_error_max_rad': fmt(max(r[6] for r in rows)),
        'tau_aw_range_s': [fmt(min(r[7] for r in rows)), fmt(max(r[7] for r in rows))],
    }


def literal_array(events, word=(2.0, 102.0), windows=((2, 18), (82, 98))):
    """Literal vs injection-free world six-column arrays on one 100-s word.

    Lemma W with W=I: a prediction adds Gamma=R_pre' Rs' Bs to B~ (Rodrigues
    Rs=Exp(-[w]h), Bs=int_0^h Exp(-[w]s)ds); a reset multiplies A~ and B~ by
    N=R_post'(I+[d]/2)R_pre with R_pre=Exp([d])'R_post. Rows are
    [a_hat-g e_z]x[A~,B~] (accelerometer, before its reset) and [B_w]x[A~,B~]
    (applied magnetic, before its reset). Times are counted from the first
    recorded prediction. Float64; non-promoting.
    """
    import numpy as np
    g = float(G)

    def skew(v):
        return np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])

    def rot(w, x, y, z):
        n = w*w+x*x+y*y+z*z
        return np.array([[w*w+x*x-y*y-z*z, 2*(x*y-w*z), 2*(x*z+w*y)],
                         [2*(x*y+w*z), w*w-x*x+y*y-z*z, 2*(y*z-w*x)],
                         [2*(x*z-w*y), 2*(y*z+w*x), w*w-x*x-y*y+z*z]])/n

    def expm(v):
        th = np.linalg.norm(v)
        k = skew(v/th) if th > 1e-14 else skew(v)
        return np.eye(3)+(np.sin(th)*k+(1-np.cos(th))*k@k if th > 1e-14 else k)

    def bint(w, h):
        om = np.linalg.norm(w)
        W = skew(w)
        if om*h < 1e-8:
            return h*np.eye(3)-h*h/2*W+h**3/6*W@W
        return h*np.eye(3)-(1-np.cos(om*h))/om**2*W+(om*h-np.sin(om*h))/om**3*W@W

    bw = next(np.array([float(x) for x in e.split()[1:]]) for e in events if e.startswith('B '))
    b = bw/np.linalg.norm(bw)
    a, bb, b0 = np.eye(3), np.zeros((3, 3)), np.zeros((3, 3))
    t, omega, drift = 0.0, 0.0, 0.0
    gram_lit, gram_free = np.zeros((6, 6)), np.zeros((6, 6))
    win = {0: [], 1: []}
    for e in events:
        f = e.split()
        if f[0] == 'P':
            q, w, h = [float(x) for x in f[1:5]], np.array([float(x) for x in f[5:8]]), float(f[8])
            gam = rot(*q).T@expm(-w*h).T@bint(w, h)
            bb, b0, t = bb+gam, b0+gam, t+h
            if word[0] <= t <= word[1]:
                omega = max(omega, float(np.linalg.norm(w)))
        elif f[0] == 'R':
            d = np.array([float(x) for x in f[1:4]])
            post = rot(*[float(x) for x in f[4:8]])
            n = post.T@(np.eye(3)+skew(d)/2)@(expm(d).T@post)
            a, bb = n@a, n@bb
            if word[0] <= t <= word[1]:
                drift = max(drift, float(np.linalg.norm(a-np.eye(3), 2)))
        elif f[0] in 'AM' and word[0] <= t <= word[1]:
            if f[0] == 'A':
                ahat = np.array([float(x) for x in f[1:4]])
                v = ahat-np.array([0.0, 0.0, g])
                rel = t-word[0]
                for i, (lo, hi) in enumerate(windows):
                    if lo <= rel <= hi:
                        win[i].append(ahat)
            else:
                v = bw
            for gram, z in ((gram_lit, np.hstack([a, bb])), (gram_free, np.hstack([np.eye(3), b0]))):
                row = skew(v)@z
                gram += row.T@row
    floor = lambda m: float(math.sqrt(max(0.0, float(np.linalg.eigvalsh(m).min()))))
    mperp = max(float(np.linalg.norm(np.cross(np.mean(win[i], axis=0), b))) for i in (0, 1))
    l1 = max(float(np.mean([np.linalg.norm(x-np.array([0.0, 0.0, g])) for x in win[i]]))/g for i in (0, 1))
    from .aggregate_floor import optimize_eps
    up = lambda x: F(x).limit_denominator(10**6)+F(1, 10**6)
    best = optimize_eps(up(omega), 1, window=16, separation=64, sigma_w=F(1, 5),
                        m_perp=up(mperp), force_l1=up(l1),
                        rows_per_window=min(len(win[0]), len(win[1])), start_offset=2,
                        field_min=F(20))
    g0 = math.sqrt(float(best[1]['six_column_floor_squared'])) if best else 0.0
    fmt = lambda x: float(f'{x:.6g}')
    return {'word_s': [word[0], word[1]], 'literal_sigma_min': fmt(floor(gram_lit)),
            'injection_free_sigma_min': fmt(floor(gram_free)),
            'A_tilde_minus_I_max': fmt(drift), 'nominal_rate_max_rad_s': fmt(omega),
            'word_nominal_transverse_mean_mps2': fmt(mperp), 'word_nominal_force_L1': fmt(l1),
            'G0_floor_service_gap_1s': fmt(g0)}


def diagnostic():
    with tempfile.TemporaryDirectory() as tmp:
        exes = _compile(tmp)
        profiles = {}
        for profile in PROFILES:
            envelope = admissibility(*profile[1:])
            if not envelope['admitted']:
                raise ValueError(f'{profile[0]} violates MARINE MOTION')
            name, args, rows, active, events = run_profile(exes, profile)
            profiles[name] = {'arguments': args, 'marine_motion': envelope,
                              **analyze(rows, active), 'literal_array': literal_array(events)}
    worst = lambda key: max(p[key] for p in profiles.values())
    return {
        'qualification': 'OU3_AW_TRACKING_SOURCE_FEASIBILITY_V1',
        'role': 'non-promoting carried falsification/feasibility audit',
        'driver_sha256': hashlib.sha256(DRIVER.read_bytes()).hexdigest(),
        'header_sha256': hashlib.sha256((REPO/HEADER).read_bytes()).hexdigest(),
        'field_world_uT': [str(x) for x in FIELD],
        'onset_s': T0, 'ramp_s': RAMP, 'end_s': T_END, 'window_s': WINDOW_S,
        'profiles': profiles,
        'pointwise_premise_refuted_on_admitted_history': worst('pointwise_ratio_to_corollary_A') > 1,
        'worst_pointwise_ratio': worst('pointwise_ratio_to_corollary_A'),
        'worst_signed_mean_ratio': worst('signed_mean_ratio_to_1_12383'),
        'worst_corollary_A_star_ratio': worst('corollary_A_star_ratio'),
        'worst_nominal_transverse_mean_mps2': worst('nominal_transverse_mean_max_mps2'),
        'worst_nominal_force_L1': worst('nominal_force_L1_max'),
        'literal_over_injection_free_min': min(
            p['literal_array']['literal_sigma_min']/p['literal_array']['injection_free_sigma_min']
            for p in profiles.values()),
        'G0_floor_below_literal_and_injection_free': all(
            p['literal_array']['G0_floor_service_gap_1s'] <= min(
                p['literal_array']['literal_sigma_min'], p['literal_array']['injection_free_sigma_min'])
            for p in profiles.values()),
        'source_uniform_certificate': False,
        'magnetic_service_certified': False,
        'theorem_closed': False,
    }


def verify_diagnostic(record):
    """Internal consistency of a committed record; no replay required."""
    if record.get('qualification') != 'OU3_AW_TRACKING_SOURCE_FEASIBILITY_V1':
        raise ValueError('wrong qualification')
    if record['source_uniform_certificate'] or record['theorem_closed']:
        raise ValueError('diagnostic cannot promote')
    if record['driver_sha256'] != hashlib.sha256(DRIVER.read_bytes()).hexdigest():
        raise ValueError('AW tracking driver changed')
    if record['header_sha256'] != hashlib.sha256((REPO/HEADER).read_bytes()).hexdigest():
        raise ValueError('shipping header changed since the committed audit')
    names = [p[0] for p in PROFILES]
    if sorted(record['profiles']) != sorted(names):
        raise ValueError('profile set changed')
    for profile in PROFILES:
        row = record['profiles'][profile[0]]
        if row['arguments'] != list(profile[1:]) or row['marine_motion'] != admissibility(*profile[1:]):
            raise ValueError(f'{profile[0]}: admissibility record differs')
        for key, threshold in (('pointwise_aw_error_sup_mps2', CORR_A_PHYSICAL_THRESHOLD),
                               ('signed_16s_mean_error_max_mps2', CORR_A_PHYSICAL_THRESHOLD),
                               ('nominal_transverse_mean_max_mps2', CORR_A_STAR_THRESHOLD)):
            ratio = {'pointwise_aw_error_sup_mps2': 'pointwise_ratio_to_corollary_A',
                     'signed_16s_mean_error_max_mps2': 'signed_mean_ratio_to_1_12383',
                     'nominal_transverse_mean_max_mps2': 'corollary_A_star_ratio'}[key]
            if abs(row[key]/float(threshold)-row[ratio]) > 1e-5*max(1, row[ratio]):
                raise ValueError(f'{profile[0]}: {ratio} inconsistent')
        if row['nominal_transverse_mean_max_mps2'] > row['nominal_16s_mean_max_mps2']*(1+1e-6):
            raise ValueError('transverse mean exceeds the mean')
    worst = lambda key: max(p[key] for p in record['profiles'].values())
    for key in ('worst_pointwise_ratio', 'worst_signed_mean_ratio', 'worst_corollary_A_star_ratio'):
        source = {'worst_pointwise_ratio': 'pointwise_ratio_to_corollary_A',
                  'worst_signed_mean_ratio': 'signed_mean_ratio_to_1_12383',
                  'worst_corollary_A_star_ratio': 'corollary_A_star_ratio'}[key]
        if record[key] != worst(source):
            raise ValueError(f'{key} is not the profile maximum')
    for key, source in (('worst_nominal_transverse_mean_mps2', 'nominal_transverse_mean_max_mps2'),
                        ('worst_nominal_force_L1', 'nominal_force_L1_max')):
        if record[key] != worst(source):
            raise ValueError(f'{key} is not the profile maximum')
    arrays = [p['literal_array'] for p in record['profiles'].values()]
    if record['G0_floor_below_literal_and_injection_free'] != all(
            a['G0_floor_service_gap_1s'] <= min(a['literal_sigma_min'], a['injection_free_sigma_min'])
            for a in arrays):
        raise ValueError('G0 floor comparison flag inconsistent')
    if abs(record['literal_over_injection_free_min']-min(
            a['literal_sigma_min']/a['injection_free_sigma_min'] for a in arrays)) > 1e-12:
        raise ValueError('literal/injection-free ratio inconsistent')
    if record['pointwise_premise_refuted_on_admitted_history'] != (record['worst_pointwise_ratio'] > 1):
        raise ValueError('refutation flag inconsistent')
    return True


def main():
    global EIGEN
    parser = argparse.ArgumentParser()
    parser.add_argument('--eigen', default=EIGEN)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--expect', type=Path)
    args = parser.parse_args()
    EIGEN = args.eigen
    record = diagnostic()
    verify_diagnostic(record)
    text = json.dumps(record, indent=2, sort_keys=True)+'\n'
    args.output.write_text(text)
    if args.expect and json.loads(args.expect.read_text()) != record:
        raise SystemExit('carried AW tracking audit differs from the committed record')
    print(text)


if __name__ == '__main__':
    main()
