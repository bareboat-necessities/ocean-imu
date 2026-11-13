#!/usr/bin/env python3
"""Conditional finite-horizon information, with all motion/gyro nuisances.

Uses estimated states and the exported existing CARMA process, never truth
columns. This is a local information diagnostic, not a global nonlinear
observability theorem. Sensor-only deliberately frees acceleration; known,
uncertain and free magnetic-reference cases are separate assumptions.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import linalg, sparse
from scipy.sparse.linalg import lsmr
from live_bias_window import whiten


def skew(v):
    x, y, z = v
    return np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])


def projected_information(J, targets):
    """A' (I - B B^dagger) A; no prior or ridge on eliminated variables."""
    nuisance = np.setdiff1d(np.arange(J.shape[1]), targets)
    A = J[:, targets].toarray()
    B = J[:, nuisance].tocsc()
    scale = np.sqrt(np.asarray(B.power(2).sum(axis=0)).ravel())
    keep = scale > 1e-15
    B = B[:, keep] @ sparse.diags(1 / scale[keep])
    residual = []
    solves = []
    for column in A.T:
        fit = lsmr(B, column, atol=1e-11, btol=1e-11, conlim=1e12, maxiter=20000)
        residual.append(column - B @ fit[0])
        solves.append(dict(status=int(fit[1]), iterations=int(fit[2]), condition=float(fit[6])))
    E = np.array(residual).T
    S = E.T @ E
    return dict(eigenvalues=np.linalg.eigvalsh(S).tolist(), matrix=S.tolist(),
                solver=solves, converged=all(s['status'] in (0, 1, 2, 4, 5) for s in solves))


def analyze(prefix, start=600., length=60., prior_scale=1., field_sd=0., sensor_dt=.005):
    d = pd.read_csv(str(prefix) + '_trace.csv')
    d = d[(d['filter'] == 'QSM') & (d.t >= start) & (d.t < start + length)].copy()
    if len(d) < 20:
        raise ValueError('at least 20 one-second nodes required')
    model = json.loads(Path(str(prefix) + '_model.json').read_text())
    n = len(d)
    h = float(np.diff(d.t).mean())
    if np.max(abs(np.diff(d.t) - h)) > 1e-6:
        raise ValueError('nonuniform nodes not supported by this diagnostic')
    A = np.kron(np.eye(3), model['A6'])
    W = np.array(model['W18'])
    W = (W + W.T) / 2
    P = linalg.solve_continuous_lyapunov(A, -W)
    scale = np.sqrt(np.maximum(P.diagonal(), 1e-15))
    As = (A * scale[None, :]) / scale[:, None]
    Ws = W / scale[None, :] / scale[:, None]
    F = linalg.expm(As * h)
    E = linalg.expm(np.block([[As, Ws], [np.zeros_like(As), -As.T]]) * h)
    Q = E[:18, 18:] @ F.T
    Lq, q_receipt = whiten(Q * prior_scale)
    Lp, p_receipt = whiten(P / scale[None, :] / scale[:, None] * prior_scale)
    ps, bs, gs, fs = .01, .05, .001, 2.
    bg, ba, field = 21*n, 21*n+3, 21*n+6
    uncertain = field_sd != 0
    nvar = 21*n + 6 + (3 if uncertain else 0)
    field_rows = 3 if np.isfinite(field_sd) and field_sd > 0 else 0
    rows = 6*n + 3*(n-1) + 18*(n-1) + 18 + field_rows
    J = sparse.lil_matrix((rows, nvar))
    rotations = d[[f'r{i}{j}' for i in range(3) for j in range(3)]].to_numpy().reshape(n, 3, 3)
    acceleration = d[['ax', 'ay', 'az']].to_numpy()
    B = np.array(model['field'])[0]
    g = np.array([0., 0., 9.80665])
    # Actual single packet noise; one accel and one magnetic packet per node.
    for k, R in enumerate(rotations):
        acc = slice(3*k, 3*k+3)
        mag = slice(3*n+3*k, 3*n+3*k+3)
        pose = slice(3*k, 3*k+3)
        J[acc, pose] = R.T @ skew(acceleration[k]-g) * ps / model['acc_sigma']
        for ax in range(3):
            J[acc, 3*n+18*k+6*ax+2] = (R.T[:, ax] * scale[6*ax+2] / model['acc_sigma']).reshape(3, 1)
        J[acc, ba:ba+3] = np.eye(3) * bs / model['acc_sigma']
        J[mag, pose] = R.T @ skew(B) * ps / model['mag_sigma']
        if uncertain:
            J[mag, field:field+3] = R.T * fs / model['mag_sigma']
    offset = 6*n
    # Endpoint integration of the nominal R over each ONE-SECOND information
    # node. This is an explicit diagnostic approximation, not a change to the
    # live 200 Hz measured-endpoint integrator. The correlated-endpoint noise
    # upper bound uses the full sigma*sqrt(h*sensor_dt), not sqrt(1/2).
    sdgyro = model['gyro_sigma'] * np.sqrt(h * sensor_dt)
    for k in range(n-1):
        row = slice(offset+3*k, offset+3*k+3)
        J[row, 3*k:3*k+3] = -np.eye(3)*ps/sdgyro
        J[row, 3*k+3:3*k+6] = np.eye(3)*ps/sdgyro
        K = .5*h*(rotations[k]+rotations[k+1])
        J[row, bg:bg+3] = K*gs/sdgyro
    process_start = offset+3*(n-1)
    for k in range(n-1):
        row = slice(process_start+18*k, process_start+18*k+18)
        J[row, 3*n+18*k:3*n+18*k+18] = -Lq @ F
        J[row, 3*n+18*k+18:3*n+18*k+36] = Lq
    offset = process_start + 18*(n-1)
    J[offset:offset+18, 3*n:3*n+18] = Lp
    if field_rows:
        J[-3:, field:field+3] = np.eye(3)*fs/field_sd
    J = J.tocsc()
    targets = np.r_[0, 1, np.arange(ba, ba+3)]
    sensor_rows = np.r_[np.arange(process_start), np.arange(rows-field_rows, rows)]
    return dict(nodes=n, node_dt=h, start=start, prior_scale=prior_scale,
                field_prior_sd=None if np.isinf(field_sd) else field_sd,
                field_assumption='free' if np.isinf(field_sd) else 'known' if field_sd == 0 else 'Gaussian prior',
                target_scales=[ps, ps, bs, bs, bs], Q=q_receipt, stationary=p_receipt,
                sensor_only=projected_information(J[sensor_rows], targets),
                carma=projected_information(J, targets),
                bias_only=projected_information(J, np.arange(ba, ba+3)),
                assumptions=['frozen exported CARMA parameters over window',
                             'constant slow biases in window', 'one-second pose/gyro quadrature',
                             'no accelerometer-bias prior', 'estimated trajectory linearization'])


def structural_checks():
    """Lie-derivative local rank for stationary, zero-motion, constant biases."""
    B = np.array([20., 0., 40.])
    g = np.array([0., 0., 9.80665])
    H = np.zeros((6, 9))
    H[:3, :3] = skew(-g)
    H[:3, 6:] = np.eye(3)
    H[3:, :3] = skew(B)
    F = np.zeros((9, 9)); F[:3, 3:6] = -np.eye(3)
    O = np.vstack([H, H@F, H@F@F])
    phi = B / np.linalg.norm(B)
    null = np.r_[phi, np.zeros(3), -skew(-g)@phi]
    return dict(stationary_constant_bias_rank=int(np.linalg.matrix_rank(O, tol=1e-9)),
                dimension=9, magnetic_axis_gauge_residual=float(np.linalg.norm(O@null)),
                includes_full_three_axis_field=True, yaw_change_required=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prefix', type=Path)
    parser.add_argument('--start', type=float, default=600)
    parser.add_argument('--field-sd', type=float, default=0)
    parser.add_argument('--prior-scale', type=float, default=1)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.prefix, args.start, prior_scale=args.prior_scale, field_sd=args.field_sd)
    result['structural'] = structural_checks()
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({key: result[key] for key in ('nodes', 'field_assumption', 'structural')}))
