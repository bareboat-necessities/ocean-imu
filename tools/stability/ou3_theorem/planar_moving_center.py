"""Finite conditional moving-center transport on a source-bound planar word.

The center uses the observed F,Q,H,R,G and AW target/event chronology. Its
Riccati gain is recomputed from its own covariance, retaining linked P/H/S/K.
This is an explicitly conditional reference map, NOT another shipping filter
or a future center certificate: H and G still require the actual finite orbit.
No physical coordinate, scheduler phase or generated target is independently
varied. The diagnostic isolates accumulated numerical/operation defects from
the large drift incurred by returning to one fixed covariance center.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from .planar_service_audit import correction_record, reset_matrix, symmetric
from .planar_service_cell import aw_floor
from .planar_service_stream import expand, records


def relative_frobenius(value, metric):
    L = np.linalg.cholesky(symmetric(metric))
    W = np.linalg.solve(L, np.eye(21))
    return float(np.linalg.norm(W @ value @ W.T, 'fro'))


def riccati_center(P, H, R):
    """Conditional real-algebra reference; all operands come from one word."""
    S = symmetric(H @ P @ H.T + R)
    np.linalg.cholesky(S)
    K = np.linalg.solve(S, H @ P).T
    A = np.eye(21) - K @ H
    return symmetric(A @ P @ A.T + K @ R @ K.T)


def physical_arc_coordinates(P, sample):
    """Finite actual-P metric on the CENTRAL PHYSICAL chart, not nominal error.

    A physical beta parameterizes the entire curved arc. No independent
    sin/cos/BA boxes are used. The transverse maximum on [-theta,theta]
    is attained at an endpoint: after maximizing over sign, its square is
    A a(t)^2+B b(t)^2+2|C|a(t)b(t), with a=t-sin(t), b=1-cos(t)
    nonnegative increasing on [0,theta]. Gauge extrema are endpoints only
    after checking a strictly positive derivative lower comparison.
    """
    g = 9.80665; theta = 2*math.atan(1/200)
    pitch = .02*math.sin(math.pi*(sample*.005)/10)
    direction = np.array([math.cos(pitch), 0., math.sin(pitch)])
    normal = np.array([-math.sin(pitch), 0., math.cos(pitch)])
    r = np.zeros(21); r[:3] = -direction; r[19] = g
    v = np.zeros(21); v[19] = g
    w = np.zeros(21); w[18:21] = g*normal
    L = np.linalg.cholesky(symmetric(P))
    wr, wv, ww = np.linalg.solve(L, np.column_stack((r, v, w))).T
    scale = np.linalg.norm(wr); u = wr/scale
    av, aw = float(u @ wv), float(u @ ww)
    derivative_lower = float(scale-theta**2*abs(av)/2-theta*abs(aw))
    if derivative_lower <= 0:
        raise ValueError('physical gauge extrema require a different bound')
    transverse = []; gauge = []
    for beta in (-theta, theta):
        eta = (math.sin(beta)-beta)*wv+(math.cos(beta)-1)*ww
        gauge.append(abs(beta*scale+float(u @ eta)))
        transverse.append(float(np.linalg.norm(eta-u*(u @ eta))))
    return {'central_physical_chart_gauge_amplitude': max(gauge),
            'central_physical_chart_transverse_curvature': max(transverse),
            'gauge_coordinate_derivative_lower': derivative_lower}


def diagnose(path, word_samples=4000):
    if word_samples < 1:
        raise ValueError('positive word length required')
    center = root = None
    words = []
    phases = {}
    last_sample = None
    root_sample = None
    reset_expected = False
    reset_applied = False
    sample_open = False
    aw_age = None
    counts = {'prediction': 0, 'acc': 0, 'mag': 0, 'S': 0, 'AW': 0}
    prefix_defect = 0.
    arc_extrema = {}
    for kind, k, a in records(path):
        if kind == 1:
            if sample_open or reset_expected:
                raise ValueError('missing sample or reset in moving-center word')
            if last_sample is not None and k != last_sample + 1:
                raise ValueError('sample gap in moving-center word')
            sample_open = True
            P, F, Q = (expand(a[i:i+225]) for i in (0, 225, 450))
            if center is None:
                center = symmetric(P)
                root = center.copy()
                root_sample = k-1
                phases[root_sample] = {'S_elapsed': float(a[909]),
                                      'S_period': float(a[910]),
                                      'AW_samples_since_observed_sync': None}
            counts['prediction'] += 1
            if aw_age is not None:
                aw_age += 1
            center = symmetric(F @ center @ F.T + Q)
        elif kind in (2, 3, 4):
            if center is None or not sample_open or reset_expected:
                raise ValueError('unpaired correction in moving-center word')
            _, H, R, _, _, _, _ = correction_record(a)
            center = riccati_center(center, H, R)
            counts[{2: 'acc', 3: 'mag', 4: 'S'}[kind]] += 1
            reset_expected = True
        elif kind == 5:
            if not reset_expected or reset_applied:
                raise ValueError('reset without correction')
            G = reset_matrix(a[225:228])
            center = symmetric(G @ center @ G.T)
            reset_applied = True
        elif kind == 8:
            if not reset_expected or not reset_applied:
                raise ValueError('post-reset without correction')
            reset_expected = False
            reset_applied = False
        elif kind == 6:
            if center is None or not sample_open or reset_expected:
                raise ValueError('invalid AW operation placement')
            if bool(a[234]):
                center = aw_floor(center, a[225:234].reshape(3, 3), (15, 16, 17))
                counts['AW'] += 1
                aw_age = 0
        elif kind == 7:
            if center is None or not sample_open or reset_expected:
                raise ValueError('incomplete sample in moving-center word')
            sample_open = False
            P = symmetric(expand(a[:225]))
            defect = relative_frobenius(P-center, P)
            prefix_defect = max(prefix_defect, defect)
            for field, value in physical_arc_coordinates(P, k).items():
                if field not in arc_extrema or (value < arc_extrema[field]['value'] if 'lower' in field else value > arc_extrema[field]['value']):
                    arc_extrema[field] = {'value': value, 'sample': k}
            last_sample = k
            phases[k] = {'S_elapsed': float(a[281]), 'S_period': float(a[282]),
                         'AW_samples_since_observed_sync': aw_age}
            if (k-root_sample) % word_samples == 0:
                # The center is inherited, NEVER reseeded to observed P.
                drift = P-root
                L = np.linalg.cholesky(root)
                W = np.linalg.solve(L, np.eye(21))
                D = W @ drift @ W.T
                words.append({'end_sample': k, 'inherited_center_defect': defect,
                              'fixed_initial_center_drift': float(np.linalg.norm(D, 'fro')),
                              'BA_y_variance_over_initial': float(P[19, 19]/root[19, 19]),
                              'AW_variances_over_initial': (P.diagonal()[15:18]/root.diagonal()[15:18]).tolist(),
                              'phase': phases[k], 'operation_counts_so_far': dict(counts)})
    if sample_open or reset_expected or center is None or not words or last_sample != words[-1]['end_sample']:
        raise ValueError('incomplete moving-center word')
    # Directly report the two inherited roots' phase mismatch, not a phase box.
    return {'qualification': 'OU3_PLANAR_CONDITIONAL_MOVING_CENTER_V1',
            'result_type': 'FINITE DIAGNOSTIC ONLY',
            'stream_sha256': hashlib.sha256(Path(path).read_bytes()).hexdigest(),
            'driver_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'sample_range': [root_sample, last_sample], 'word_samples': word_samples,
            'initial_phase': phases[root_sample], 'words': words,
            'maximum_sample_relative_center_defect': prefix_defect,
            'finite_actual_P_physical_arc_coordinates': arc_extrema,
            'central_physical_chart_is_actual_nominal_error_chart': False,
            'physical_arc_metric_amplitude_uniform_in_time': False,
            'center_reseeded_at_word_boundaries': False,
            'center_gain_recomputed_from_linked_covariance': True,
            'center_uses_observed_future_H_and_G': True,
            'center_is_autonomous_future_construction': False,
            'uniform_center_defect_bound_verified': False,
            'nonlinear_mean_covariance_cell_verified': False,
            'joint_cell_forward_invariant': False,
            'all_time_magnetic_service_verified': False, 'theorem_closed': False,
            'structures_preserved': 'all 21 covariance coordinates, inherited center, actual ordered prediction/correction/reset/AW events and joint S/AW chronology',
            'relaxations_introduced': 'conditional real-algebra covariance reference driven by finite actual F,Q,H,R,G/targets; no independent coefficient ranges or future extension',
            'entry_in_radius_inequality': 'center defect contributes to q_P only after a uniform same-history future coefficient and arithmetic enclosure; physical C_Q alpha remains separate'}


def verify_report(out, stream_hash):
    if out.get('qualification') != 'OU3_PLANAR_CONDITIONAL_MOVING_CENTER_V1' or out.get('result_type') != 'FINITE DIAGNOSTIC ONLY':
        raise ValueError('moving center must remain a finite diagnostic')
    if out.get('stream_sha256') != stream_hash or out.get('driver_sha256') != hashlib.sha256(Path(__file__).read_bytes()).hexdigest():
        raise ValueError('moving center source/stream changed')
    for key in ('center_reseeded_at_word_boundaries', 'center_is_autonomous_future_construction',
                'central_physical_chart_is_actual_nominal_error_chart', 'physical_arc_metric_amplitude_uniform_in_time',
                'uniform_center_defect_bound_verified', 'nonlinear_mean_covariance_cell_verified',
                'joint_cell_forward_invariant', 'all_time_magnetic_service_verified', 'theorem_closed'):
        if out.get(key) is not False:
            raise ValueError('moving-center promotion: '+key)
    if out.get('center_uses_observed_future_H_and_G') is not True:
        raise ValueError('missing conditional coefficient qualification')
    return True


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stream', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    out = diagnose(args.stream)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True)+'\n')
    print(json.dumps(out, sort_keys=True))


if __name__ == '__main__':
    main()
