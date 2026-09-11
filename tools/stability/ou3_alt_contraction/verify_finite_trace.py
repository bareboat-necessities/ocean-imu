"""Implementation-only checks of finite physical events in literal shipping words.

No source qualification, rho search, covariance-consistency premise, or
finite-precision certificate is produced.  Local discrepancies compare host
binary32 observations to real formulas evaluated in binary64.  Tolerances are
regression tolerances, never disturbance bounds for the theorem.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
import json
from pathlib import Path
import struct

import numpy as np

BASE = 549
DT = float(np.float32(1/200))
GROUPS = ('attitude', 'gyro_bias', 'velocity', 'position', 'S', 'aw', 'bias_error', 'beta_true')
# Declared before evaluating the regression. These are not theorem bounds.
REGRESSION_RTOL = 4e-5
REGRESSION_ATOL = 4e-6


def qmul(a, b):
    aw, av, bw, bv = a[0], a[1:], b[0], b[1:]
    return np.r_[aw*bw-av@bv, aw*bv+bw*av+np.cross(av, bv)]


def qconj(q):
    return q*np.array([1., -1., -1., -1.])


def rotation(q):
    # Projective form, valid for any nonzero quaternion. Host normalization
    # error is recorded separately from this real-arithmetic rotation.
    w, v = q[0], q[1:]
    return np.eye(3)+2*(w*skew(v)+skew(v)@skew(v))/(q@q)


def skew(v):
    x, y, z = v
    return np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])


def step_quaternion(d):
    u = d@d
    if u < 1e-4:
        w, k = 1-u/8+u*u/384, .5-u/48+u*u/3840
    else:
        theta = np.sqrt(u)
        w, k = np.cos(theta/2), np.sin(theta/2)/theta
    return np.r_[w, k*d]


@lru_cache(maxsize=8192)
def physical(t: float, live: float):
    """Same analytical reference as the host harness; not a theorem source cover."""
    amp = np.array([.12, .07, .16])
    w = 2*np.pi*np.array([.17, .23, .13])
    phase = np.array([.21, -.35, .6])
    angle = w*t+phase
    p, v, a = amp*np.sin(angle), amp*w*np.cos(angle), -amp*w*w*np.sin(angle)
    S = amp/w*(np.cos(w*live+phase)-np.cos(angle))
    r, b, y = .045*np.sin(.6*t), .03*np.cos(.83*t), .07*np.sin(.37*t)
    qr = np.array([np.cos(r/2), np.sin(r/2), 0., 0.])
    qb = np.array([np.cos(b/2), 0., np.sin(b/2), 0.])
    qy = np.array([np.cos(y/2), 0., 0., np.sin(y/2)])
    q = qconj(qmul(qmul(qy, qb), qr))
    beta = np.array([.025, -.018, .012])
    bg = np.array([.0001+.00002*np.sin(.011*t), -.0002+.00001*np.cos(.009*t), .0001])
    return q, bg, v, p, S, a, beta


@dataclass
class Event:
    raw: np.ndarray

    def __post_init__(self):
        if self.raw.ndim != 1 or len(self.raw) < BASE or not np.isfinite(self.raw).all():
            raise ValueError('invalid event record')

    @property
    def kind(self): return int(self.raw[1])
    @property
    def active(self): return bool(self.raw[3])
    @property
    def x(self): return self.raw[16:37]
    @property
    def quat(self): return self.raw[37:41]
    @property
    def P(self): return self.raw[41:482].reshape(21, 21)
    @property
    def Racc(self): return self.raw[482:491].reshape(3, 3)
    @property
    def Rmag(self): return self.raw[491:500].reshape(3, 3)
    @property
    def RS(self): return self.raw[500:509].reshape(3, 3)
    @property
    def truth(self): return physical(float(self.raw[4]), float(self.raw[5]))
    @property
    def z(self):
        q, bg, v, p, S, a, beta = self.truth
        error = qmul(q, qconj(self.quat))
        if abs(error[0]) < 1e-10: raise ValueError('Cayley pole in regression')
        return np.r_[2*error[1:]/error[0], bg-self.x[3:6], v-self.x[6:9],
                     p-self.x[9:12], S-self.x[12:15], a-self.x[15:18],
                     beta-self.x[18:21], beta]
    @property
    def operands(self):
        if self.kind not in (12, 22, 32) or len(self.raw) != BASE+141:
            raise ValueError('innovation operands missing')
        a = self.raw[BASE:]
        return a[:3], a[3:66].reshape(21, 3), a[66:75].reshape(3, 3), a[75:138].reshape(21, 3), a[138:141]


def read_samples(path: Path):
    """Stream bounded-size records; malformed/truncated files fail closed."""
    with path.open('rb') as stream:
        previous = None
        while header := stream.read(8):
            if len(header) != 8: raise ValueError('truncated sample header')
            tick, count = struct.unpack('<iI', header)
            if previous is not None and tick != previous+1:
                raise ValueError('nonconsecutive physical samples')
            previous = tick
            if count > 256: raise ValueError('unreasonable event count')
            events = []
            for _ in range(count):
                length = stream.read(4)
                if len(length) != 4: raise ValueError('truncated event size')
                size, = struct.unpack('<I', length)
                if not BASE <= size <= 4096: raise ValueError('invalid event size')
                data = stream.read(8*size)
                if len(data) != 8*size: raise ValueError('truncated event payload')
                event = Event(np.frombuffer(data, dtype='<f8').copy())
                if int(event.raw[0]) != tick: raise ValueError('event from another sample')
                events.append(event)
            yield tick, events


class Audit:
    def __init__(self):
        self.errors = defaultdict(lambda: {'max_abs': 0., 'max_scaled': 0., 'checks': 0})
        self.counts = Counter()
        self.prefixes = 0
        self.min_cayley_scalar = 1.
        self.min_injection_denominator = float('inf')
        self.S_nonzero = 0
        self.latent_ba = 0
        self.local_attitude_slots_zero = True
        self.forcing_max = defaultdict(float)

    def check(self, label, actual, expected):
        a, b = np.asarray(actual), np.asarray(expected)
        if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
            raise AssertionError(f'{label}: invalid values or dimensions')
        absolute = float(np.max(np.abs(a-b), initial=0))
        # One global scale per operand: tiny off-diagonal roundoff must not be
        # treated as a relative error of infinity. Retain absolute error too.
        scale = max(1., float(np.max(np.abs(a), initial=0)), float(np.max(np.abs(b), initial=0)))
        entry = self.errors[label]
        entry['max_abs'] = max(entry['max_abs'], absolute)
        entry['max_scaled'] = max(entry['max_scaled'], absolute/scale)
        entry['checks'] += 1
        if absolute > REGRESSION_ATOL+REGRESSION_RTOL*scale:
            raise AssertionError(f'{label}: local implementation discrepancy {absolute}, scale {scale}')

    def identity_state(self, before, after, label):
        self.check(label+'.mean', after.x, before.x)
        self.check(label+'.attitude', after.quat, before.quat)

    def prediction(self, before: Event, after: Event):
        self.counts['prediction'] += 1
        h = after.raw[6]
        if len(after.raw) != BASE+360: raise ValueError('prediction operands missing')
        a = after.raw[BASE:]
        Fa, Qa, Fl, Ql = a[:36].reshape(6, 6), a[36:72].reshape(6, 6), a[72:216].reshape(12, 12), a[216:360].reshape(12, 12)
        phi = after.raw[11]
        self.check('prediction.mean.linear', after.x[6:18], Fl@before.x[6:18])
        self.check('prediction.mean.gyro_bias', after.x[3:6], before.x[3:6])
        self.check('prediction.mean.acc_bias', after.x[18:21], phi*before.x[18:21])
        qnom = step_quaternion(-after.raw[530:533]*h)
        qpred = qmul(qnom, before.quat); qpred /= np.linalg.norm(qpred)
        self.check('prediction.quaternion', rotation(after.quat), rotation(qpred))
        F = np.zeros((21, 21)); Q = np.zeros((21, 21))
        F[:6, :6], F[6:18, 6:18], F[18:21, 18:21] = Fa, Fl, phi*np.eye(3)
        Q[:6, :6], Q[6:18, 6:18] = Qa, Ql
        if after.active:
            tau = max(.001, after.raw[542])
            Q[18:21, 18:21] = after.raw[533:542].reshape(3, 3)*(-.5*tau*np.expm1(-2*h/tau))
        self.check('prediction.covariance', after.P, F@before.P@F.T+Q)
        z = before.z.copy(); tb, ta = before.truth, after.truth
        qphys = qmul(ta[0], qconj(tb[0]))
        errorq = qmul(qmul(qphys, np.r_[2., z[:3]]), qconj(qnom))
        z[:3] = 2*errorq[1:]/errorq[0]
        z[3:6] += ta[1]-tb[1]
        defect = np.r_[ta[2], ta[3], ta[4], ta[5]]-Fl@np.r_[tb[2], tb[3], tb[4], tb[5]]
        z[6:18] = Fl@before.z[6:18]+defect
        z[18:21] = phi*before.z[18:21]+ta[6]-phi*tb[6]
        z[21:24] = ta[6]
        self.check('prediction.finite_joint24', after.z, z)
        self.forcing_max['physical_translation_defect_linf'] = max(self.forcing_max['physical_translation_defect_linf'], float(np.max(abs(defect))))
        anchorq = qmul(qphys, qconj(qnom))
        self.forcing_max['physical_attitude_anchor_norm'] = max(self.forcing_max['physical_attitude_anchor_norm'], float(np.linalg.norm(2*anchorq[1:]/anchorq[0])))
        return z-after.z

    def measurement(self, begin: Event, innovation: Event, mean: Event, joseph: Event, end: Event):
        kind = {12: 'accelerometer', 22: 'magnetometer', 32: 'S_zero'}[innovation.kind]
        self.counts[kind] += 1
        if np.any(begin.x[:3] != 0): self.local_attitude_slots_zero = False
        r, N, Sigma, K, measured = innovation.operands
        H = np.zeros((3, 21)); Hbar = np.zeros((3, 24))
        c, Rhat = innovation.z[:3], rotation(innovation.quat)
        C = skew(c); D = 4+c@c; U = 4*np.eye(3)+2*C+np.outer(c, c)
        E = np.eye(3)+(4*C+2*C@C)/D
        if kind == 'S_zero':
            H[:, 12:15] = np.eye(3); Hbar[:, 12:15] = np.eye(3)
            noise = innovation.RS; nu = -innovation.truth[4]
            expected_r = -innovation.x[12:15]
            if np.linalg.norm(nu) > 1e-8: self.S_nonzero += 1
        else:
            if kind == 'accelerometer':
                f = Rhat@(innovation.x[15:18]-np.array([0., 0., innovation.raw[543]]))
                H[:, 15:18], H[:, 18:21] = Rhat, np.eye(3)
                Hbar[:, 15:18], Hbar[:, 18:21] = E@Rhat, np.eye(3)
                temp = innovation.raw[544:547]*(innovation.raw[547]-innovation.raw[548])
                expected_r = measured-f-innovation.x[18:21]-temp
                nu = measured-(E@Rhat@(innovation.truth[5]-np.array([0., 0., innovation.raw[543]]))+innovation.truth[6]+temp)
                noise = innovation.Racc
            else:
                f = Rhat@innovation.raw[527:530]
                noise = innovation.Rmag; expected_r = measured-f; nu = measured-E@f
            H[:, :3] = -skew(f); Hbar[:, :3] = -U@skew(f)/D
        Hg = H.copy()
        if not innovation.active: Hg[:, 18:21] = 0
        expected_N = innovation.P@Hg.T
        if not innovation.active: expected_N[18:21] = 0
        expected_Sigma = H@innovation.P@H.T+noise
        self.check(kind+'.residual', r, expected_r)
        self.check(kind+'.finite_residual_secant', r, Hbar@innovation.z+nu)
        self.check(kind+'.masked_numerator', N, expected_N)
        self.check(kind+'.innovation', Sigma, expected_Sigma)
        self.check(kind+'.inverse_free_gain', K@Sigma, N)
        q = np.linalg.solve(Sigma, r)
        self.check(kind+'.inverse_free_correction', N@q, K@r)
        self.check(kind+'.mean', mean.x, innovation.x+K@r)
        J = innovation.P-K@N.T-N@K.T+K@Sigma@K.T
        self.check(kind+'.Joseph', joseph.P, J)
        if kind == 'accelerometer' and not innovation.active:
            B0 = innovation.P[18:21, 18:21]
            self.check('H18.latent_BA_innovation', Sigma,
                       H[:, :18]@innovation.P[:18, :18]@H[:, :18].T+noise+B0)
            self.check('H18.held_zero_cross_blocks', innovation.P[:18, 18:21], np.zeros((18, 3)))
            if np.linalg.norm(B0) > 0: self.latent_ba += 1
        # Finite physical correction, including the exact Cayley denominator.
        z = innovation.z.copy(); delta = N@q; d = delta[:3]
        w, k = step_quaternion(d)[0], None
        u = d@d
        k = .5-u/48+u*u/3840 if u < 1e-4 else np.sin(np.sqrt(u)/2)/np.sqrt(u)
        W = 2*w+k*c@d
        self.min_injection_denominator = min(self.min_injection_denominator, abs(float(W)))
        z[:3] -= (k/W)*U@d
        z[3:21] -= delta[3:21]
        radius = innovation.raw[10]
        bhat = z[21:24]-z[18:21]
        norm = np.linalg.norm(bhat)
        if radius > 0 and norm > radius: bhat *= radius/norm
        z[18:21] = z[21:24]-bhat
        self.check(kind+'.finite_joint24', end.z, z)
        return z-end.z

    def reset(self, before, injected, transported, cleared, after):
        self.counts['reset'] += 1
        d = before.x[:3]; corr = step_quaternion(d)
        self.counts['reset_polynomial' if d@d < 1e-4 else 'reset_axis_angle'] += 1
        q = qmul(corr, before.quat); q /= np.linalg.norm(q)
        self.check('reset.quaternion', rotation(injected.quat), rotation(q))
        G = np.eye(21); G[:3, :3] += .5*skew(d)
        self.check('reset.covariance', transported.P, G@injected.P@G.T)
        self.check('reset.clear', cleared.x[:3], np.zeros(3))
        self.check('reset.nonattitude', cleared.x[3:], before.x[3:])
        self.check('reset.physical_truth', after.z[21:24], before.z[21:24])

    def projection(self, before, after):
        self.counts['projection'] += 1
        b, radius = before.x[18:21].copy(), before.raw[10]
        if radius > 0 and np.linalg.norm(b) > radius:
            b *= radius/np.linalg.norm(b); self.counts['projection_active'] += 1
        self.check('projection.mean', after.x[18:21], b)
        self.check('projection.covariance_identity', after.P, before.P)
        expected = before.truth[6]-b
        self.check('projection.same_beta_error', after.z[18:21], expected)

    def edge(self, before, after):
        self.counts['mode_setter'] += 1
        self.identity_state(before, after, 'mode')
        P = before.P.copy()
        if before.active != after.active:
            if after.active:
                self.counts['H18_to_A21'] += 1
                # Held covariance invariant makes the enable floor a fixed
                # point in this regression. Do not infer this for arbitrary P.
            else:
                self.counts['A21_to_H18'] += 1
                P[18:21, :18] = 0; P[:18, 18:21] = 0
        self.check('mode.covariance', after.P, P)

    def report(self):
        return {'operation_counts': dict(self.counts), 'observed_core_prefixes': self.prefixes,
                'local_correspondence_errors': dict(self.errors),
                'nonzero_physical_S_events': self.S_nonzero,
                'held_latent_BA_innovations_checked': self.latent_ba,
                'zero_local_attitude_slots_at_measurement_entry': self.local_attitude_slots_zero,
                'minimum_observed_injection_denominator_abs': self.min_injection_denominator,
                'observed_forcing': dict(self.forcing_max),
                'source_uniform_coverage': False, 'domain_retention_certified': False}


def verify_word(path: Path):
    audit = Audit(); n = 0; origin = None; previous_end = None; first_error = None
    for tick, events in read_samples(path):
        n += 1
        if not events or events[0].kind != 1: raise AssertionError('word sample lacks prediction entry')
        if origin is None: origin = events[0].raw[5]
        if any(e.raw[5] != origin for e in events): raise AssertionError('Live S origin was restarted')
        if any(int(e.raw[0]) != tick for e in events): raise AssertionError('source tick detached')
        if first_error is None: first_error = events[0].z.copy()
        if previous_end is not None:
            audit.identity_state(previous_end, events[0], 'sample_predecessor')
            audit.check('sample_predecessor.covariance', events[0].P, previous_end.P)
            audit.check('sample_predecessor.physical_error', events[0].z, previous_end.z)
        audit.prefixes += len(events)
        prediction_before = None; measurements = {}; resets = []; projections = []; modes = []
        last = None
        for e in events:
            if e.kind == 1: prediction_before = e
            elif e.kind == 2:
                if prediction_before is None: raise AssertionError('prediction lacks predecessor')
                audit.prediction(prediction_before, e)
            elif e.kind == 3:
                if last.kind != 2: raise AssertionError('floor not after OU prediction')
                audit.counts['aw_floor_check'] += 1
                audit.identity_state(last, e, 'floor')
                target = last.raw[518:527].reshape(3, 3)
                expected = last.P.copy()
                if bool(last.raw[13]):
                    delta = target-last.P[15:18, 15:18]; delta = (delta+delta.T)/2
                    vals, vectors = np.linalg.eigh(delta)
                    expected[15:18, 15:18] += (vectors*np.maximum(vals, 0))@vectors.T
                    audit.counts['aw_floor_applied'] += 1
                audit.check('floor.covariance', e.P, expected)
            elif e.kind == 4:
                audit.identity_state(last, e, 'symmetrize')
                audit.check('symmetrize.covariance', e.P, (last.P+last.P.T)/2)
            elif e.kind in (10, 20, 30):
                measurements[e.kind] = [e]
            elif e.kind in (12, 13, 14, 22, 23, 24, 32, 33, 34):
                measurements[(e.kind//10)*10].append(e)
            elif e.kind in (11, 21, 31):
                data = measurements.pop(e.kind-1)
                if len(data) == 1:
                    audit.counts[{11:'accelerometer_rejected',21:'magnetometer_rejected',31:'S_zero_rejected'}[e.kind]] += 1
                    audit.identity_state(data[0], e, 'rejected')
                    audit.check('rejected.covariance', e.P, data[0].P)
                elif len(data) == 4: audit.measurement(*data, e)
                else: raise AssertionError('incomplete measurement graph')
            elif e.kind == 40: resets.append([e])
            elif e.kind in (42, 43, 44): resets[-1].append(e)
            elif e.kind == 41:
                data = resets.pop()
                if len(data) != 4: raise AssertionError('reset nonfinite/early-return branch not covered by regression')
                audit.reset(*data, e)
            elif e.kind == 50: projections.append(e)
            elif e.kind == 51: audit.projection(projections.pop(), e)
            elif e.kind == 60: modes.append(e)
            elif e.kind == 61: audit.edge(modes.pop(), e)
            else: raise AssertionError(f'unhandled observed operation {e.kind}')
            last = e
        if measurements or resets or projections or modes: raise AssertionError('unclosed event scope')
        # The scheduler runs after prediction/floor/symmetrization and before
        # the accelerometer. Check its actual tolerance and service credit.
        before_schedule = next(e for e in events if e.kind == 4)
        after_schedule = next(e for e in events if e.kind in (10, 30))
        f = np.float32
        period, elapsed, dt = f(before_schedule.raw[8]), f(before_schedule.raw[9]), f(before_schedule.raw[6])
        total = f(elapsed+dt)
        tol = f(f(16)*f(np.finfo(np.float32).eps)*max(f(1), period))
        due = not (f(total+tol) < period)
        elapsed_next = (f(np.fmod(total, period)) if total >= period else f(0)) if due else total
        if due != any(e.kind == 30 for e in events): raise AssertionError('scheduled S event detached')
        audit.check('scheduler.elapsed', after_schedule.raw[9], float(elapsed_next))
        audit.counts['S_due' if due else 'S_not_due'] += 1
        previous_end = events[-1]
    if n != 600: raise AssertionError(f'expected literal 600-step word, got {n}')
    if audit.counts['prediction'] != 600 or audit.counts['accelerometer'] != 600:
        raise AssertionError('missing mandatory per-sample operations')
    if not audit.local_attitude_slots_zero: raise AssertionError('pending attitude scratch at measurement entry')
    report = audit.report(); report.update(steps=n, live_origin_s=float(origin), duration_s=600*DT,
        first_error_by_group={g: first_error[3*i:3*i+3].tolist() for i,g in enumerate(GROUPS)},
        unobserved_or_unqualified_branches=['innovation repair/rejection', 'nonfinite reset/projection',
            'active radial projection', 'tilt watchdog reset', 'frontend/tuner/guard universal graph'],
        complete_literal_runtime_prefix_coverage=False)
    return report


def verify(directory: Path):
    words = {name: verify_word(directory/(name+'.bin')) for name in ('H18', 'A21', 'H18-A21')}
    if words['H18-A21']['operation_counts'].get('H18_to_A21') != 1:
        raise AssertionError('actual H18->A21 transition missing')
    if any(words['H18']['first_error_by_group']['S']):
        raise AssertionError('fresh one-time-centered S error is not zero')
    if not words['H18']['held_latent_BA_innovations_checked']:
        raise AssertionError('H18 latent BA not exercised')
    return {'qualification': 'FINITE_PHYSICAL_EVENT_CORRESPONDENCE_REGRESSION_ONLY',
            'words': words, 'tolerances': {'absolute': REGRESSION_ATOL, 'relative': REGRESSION_RTOL},
            'tolerances_are_theorem_disturbance_bounds': False,
            'complete_source_uniform_word': False, 'storage_feasibility_attempted': False,
            'ALT_LIVE_PASS': False, 'ALT_STARTUP_PASS': False, 'ALT_END_TO_END_PASS': False}


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.directory)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({k: v['operation_counts'] for k, v in result['words'].items()}, indent=2))
