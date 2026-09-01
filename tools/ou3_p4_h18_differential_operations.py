#!/usr/bin/env python3
"""Shared H=18 differential operation algebra for OU-III P4 routes.

PR #449 and PR #450 use different Lyapunov/contraction arguments, but they must
not maintain different derivatives of the shipping filter.  This module is the
common finite-angle differential layer:

* exact Cayley attitude prediction with the deployed quaternion correction;
* exact H-mode integrated-OU linear-chain state transport;
* finite-angle accelerometer and magnetometer residual maps;
* the S=0 pseudo-measurement residual;
* canonical local H matrices used only to compute the shipping Joseph gain;
* accepted measurement correction followed immediately by the deployed
  quaternion/Cayley reset.

All state derivatives are outward interval AD derivatives.  The module contains
no theorem-promotion logic, no scalar contraction test, and no trajectory fit.
#450 composes these maps into a whole-word generalized Jacobian.  #449 uses the
same maps to transport/pull back directional Joseph forms before scalarization.
"""
from __future__ import annotations

import math
from typing import Sequence

from ou3_interval import Interval
import ou3_interval_ad as AD
import ou3_p5_full_h_prefix_cells as H

N = 18


def I(x: float) -> Interval:
    """Return a point interval."""
    return Interval.point(float(x))


def box_symmetric(a: float) -> Interval:
    """Return an outward symmetric interval [-|a|,|a|]."""
    a = abs(float(a))
    return Interval(math.nextafter(-a, -math.inf), math.nextafter(a, math.inf))


def ad_constant_interval(x: Interval, n: int = N) -> AD.AD:
    """Lift an interval as an AD constant."""
    return AD.constant(x, n)


def ad_matvec_interval(A, x: Sequence[AD.AD], n: int = N) -> list[AD.AD]:
    """Multiply an interval matrix by an AD vector without losing derivatives."""
    out = []
    for row in A:
        y = AD.constant(0.0, n)
        for a, b in zip(row, x):
            y = y + AD.constant(a, n) * b
        out.append(y)
    return out


def state_values(z: Sequence[AD.AD]) -> list[Interval]:
    """Extract state value intervals from an AD state."""
    return [x.val for x in z]


def prediction(z: Sequence[AD.AD], F, Rstep, domain: dict, h: float) -> list[AD.AD]:
    """Finite-angle H-mode prediction with exact 18-state cross derivatives."""
    if len(z) != N:
        raise ValueError("H differential prediction requires 18 states")
    c = list(z[:3])
    bg = list(z[3:6])
    transported = ad_matvec_interval(Rstep, c)
    Bstep = [[F[i][3 + j] for j in range(3)] for i in range(3)]
    db = ad_matvec_interval(Bstep, bg)
    wdist = float(
        domain["startup"][
            "effective_deterministic_gyro_transport_disturbance_upper_rad_s"
        ]
    )
    extra = box_symmetric(h * wdist)
    # Deterministic source disturbance enlarges the state value but is not an
    # independent filter-state coordinate, so it contributes no state derivative.
    db = [AD.AD(x.val + extra, x.der) for x in db]
    cp = AD.deployed_correct_cayley(transported, db)

    out = list(z)
    out[:3] = cp
    # H mode: gyro bias is held; v,p,S,a_w use the exact integrated-OU chain.
    for i in range(3, N):
        y = AD.constant(0.0, N)
        for j in range(3, N):
            y = y + AD.constant(F[i][j], N) * z[j]
        out[i] = y
    return out


def accelerometer_residual(
    z: Sequence[AD.AD], force: Sequence[Interval]
) -> list[AD.AD]:
    """Exact finite-angle accelerometer residual (R(c)-I)f + R(c) delta_a_w."""
    R = AD.rotation_from_cayley(z[:3])
    f = [ad_constant_interval(x) for x in force]
    aw = list(z[15:18])
    Rf = AD.matvec(R, f)
    Raw = AD.matvec(R, aw)
    return [Rf[i] - f[i] + Raw[i] for i in range(3)]


def magnetometer_residual(
    z: Sequence[AD.AD], mag: Sequence[Interval]
) -> list[AD.AD]:
    """Exact finite-angle magnetometer residual (R(c)-I)m."""
    R = AD.rotation_from_cayley(z[:3])
    m = [ad_constant_interval(x) for x in mag]
    Rm = AD.matvec(R, m)
    return [Rm[i] - m[i] for i in range(3)]


def S_residual(z: Sequence[AD.AD]) -> list[AD.AD]:
    """Exact S=0 pseudo-measurement residual."""
    return [-z[12 + i] for i in range(3)]


def zero_matrix(rows: int, cols: int):
    return [[I(0.0) for _ in range(cols)] for _ in range(rows)]


def H_acc_canonical(force: Sequence[Interval]):
    """Local shipping H_a=[-[f]_x, I_aw] used for the Joseph gain."""
    fx, fy, fz = force
    Hm = zero_matrix(3, N)
    Hm[0][1] = fz
    Hm[0][2] = -fy
    Hm[1][0] = -fz
    Hm[1][2] = fx
    Hm[2][0] = fy
    Hm[2][1] = -fx
    for i in range(3):
        Hm[i][15 + i] = I(1.0)
    return Hm


def H_mag_canonical(mag: Sequence[Interval]):
    """Local shipping H_m=-[m]_x used for the Joseph gain."""
    mx, my, mz = mag
    Hm = zero_matrix(3, N)
    Hm[0][1] = mz
    Hm[0][2] = -my
    Hm[1][0] = -mz
    Hm[1][2] = mx
    Hm[2][0] = my
    Hm[2][1] = -mx
    return Hm


def accepted_update(Pm, z: Sequence[AD.AD], Hm, Rm, residual: Sequence[AD.AD]):
    """Apply the shipping interval Joseph gain and immediate quaternion reset.

    The covariance cell is computed by the existing validated full-matrix
    backend.  The state correction uses the exact nonlinear residual AD map;
    attitude injection/reset is the deployed quaternion/Cayley composition.
    """
    cell = H._measurement_cell(Pm, Hm, Rm, [x.val for x in residual])
    K = cell["K"]
    dx = ad_matvec_interval(K, residual)
    d = [-x for x in dx[:3]]
    cp = AD.deployed_correct_cayley(z[:3], d)
    out = list(z)
    out[:3] = cp
    for i in range(3, N):
        out[i] = z[i] - dx[i]
    return cell["P_accepted"], out, cell


def residual_jacobian(residual: Sequence[AD.AD]):
    """Return the outward full-state Jacobian of one exact residual map."""
    return AD.jacobian(residual)


def state_jacobian(z: Sequence[AD.AD]):
    """Return the outward 18x18 Jacobian of one composed state map."""
    J = AD.jacobian(z)
    if len(J) != N or any(len(row) != N for row in J):
        raise RuntimeError("shared H differential state Jacobian is not 18x18")
    return J
