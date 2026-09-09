#!/usr/bin/env python3
"""Uniform Euclidean sector for the shipping radial accelerometer-bias clamp.

This is an analytical P4 ingredient, not a replay fit.  With true physical
bias beta, pre-projection corrected error e and projection radius R, shipping
projection is represented in error coordinates by

    F_R(e,beta) = beta - Pi_R(beta-e),

where Pi_R is Euclidean projection onto the closed radius-R estimate ball.
For every pair (e1,beta1),(e2,beta2), including points on different sides of
the saturation boundary,

    ||F_R(e1,beta1)-F_R(e2,beta2)||^2
       <= ||e1-e2||^2 + ||beta1-beta2||^2.                 (PSECTOR)

Proof.  Away from the boundary the Jacobian J of Pi_R is symmetric and
0<=J<=I.  Inside the ball J=I.  Outside, for x=r*u,

    J = (R/r) (I-u*u^T),

whose eigenvalues are 0,R/r,R/r.  Therefore

    dF = J de + (I-J) dbeta,
    [J,I-J][J,I-J]^T = J^2+(I-J)^2 <= I.

At ||x||=R every Clarke generalized Jacobian is a convex combination of the
one-sided limits and is again symmetric with spectrum in [0,1].  F_R is
locally Lipschitz.  Integrating its a.e. directional derivative along the line
segment between two points gives PSECTOR.  Thus no fixed projection multiplier
s and no saturation-pattern enumeration is needed.

PSECTOR is exactly the joint recurrence structure needed here: beta is a
physical/source state, not an independent bias-error input.  If beta is fixed,
projection is nonexpansive in e alone.  If beta varies under BIAS1, its actual
recurrence budget enters once through ||Delta beta||^2.  The estimate ball is
also invariant by definition: ||Pi_R(beta-e)||<=R.

This lemma closes the radial projection *operator* sector in exact real
arithmetic.  It does not enclose floating-point rounding, source-dependent
Kalman coefficients, the physical BIAS1 driver family, or the word-entry set,
and therefore does not promote P4/P5 by itself.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


def project_ball(value, radius):
    value = np.asarray(value, dtype=float)
    if value.ndim != 1 or not np.isfinite(value).all():
        raise ValueError("finite vector required")
    if not math.isfinite(radius) or radius < 0:
        raise ValueError("finite nonnegative radius required")
    if radius == 0:
        return np.zeros_like(value)
    scale = float(np.max(np.abs(value), initial=0.))
    if scale == 0:
        return value.copy()
    direction = value/scale
    length = float(np.linalg.norm(direction))
    if scale <= radius/length:
        return value.copy()
    return direction*(radius/length)


def bias_error_projection(error_pre, beta, radius):
    error_pre = np.asarray(error_pre, dtype=float)
    beta = np.asarray(beta, dtype=float)
    if error_pre.shape != beta.shape:
        raise ValueError("error and true bias must have the same shape")
    return beta-project_ball(beta-error_pre, radius)


def projection_jacobian(estimate_pre, radius):
    """One classical Jacobian away from ||estimate||=R; boundary is Clarke."""
    x = np.asarray(estimate_pre, dtype=float)
    if x.ndim != 1 or not np.isfinite(x).all():
        raise ValueError("finite estimate vector required")
    if not math.isfinite(radius) or radius < 0:
        raise ValueError("finite nonnegative radius required")
    n = len(x)
    length = float(np.linalg.norm(x))
    if radius == 0:
        return np.zeros((n, n))
    if length < radius:
        return np.eye(n)
    if length > radius:
        u = x/length
        return (radius/length)*(np.eye(n)-np.outer(u, u))
    raise ValueError("classical Jacobian is set-valued on projection boundary")


def stacked_sector_ratio(jacobian):
    """Squared norm of [J,I-J], which is <=1 for every 0<=J<=I."""
    j = np.asarray(jacobian, dtype=float)
    if j.ndim != 2 or j.shape[0] != j.shape[1]:
        raise ValueError("square Jacobian required")
    stack = np.column_stack((j, np.eye(len(j))-j))
    return float(np.linalg.norm(stack, 2)**2)


def build_report(radius):
    if not math.isfinite(radius) or radius <= 0:
        raise ValueError("positive finite radius required")
    # Exact formulas are the certificate.  Representative matrices merely pin
    # the implementation of all three regimes; they are not a sampled proof.
    inside = projection_jacobian(np.array([0.25*radius, 0., 0.]), radius)
    outside = projection_jacobian(np.array([2*radius, 0., 0.]), radius)
    tangent_limit = np.diag([0., 1., 1.])
    boundary_mid = .5*np.eye(3)+.5*tangent_limit
    ratios = {
        "inside": stacked_sector_ratio(inside),
        "outside": stacked_sector_ratio(outside),
        "boundary_clarke_midpoint": stacked_sector_ratio(boundary_mid),
    }
    return {
        "lemma": "OU3_RADIAL_BIAS_PROJECTION_JOINT_SECTOR_V1",
        "map": "F_R(e,beta)=beta-Pi_R(beta-e)",
        "sector": "||Delta F||^2 <= ||Delta e||^2 + ||Delta beta||^2",
        "projection_radius": radius,
        "exact_real_operator_sector_closed": True,
        "saturated_branch_included_analytically": True,
        "unsaturated_branch_included_analytically": True,
        "boundary_clarke_branch_included_analytically": True,
        "fixed_multiplier_required": False,
        "saturation_pattern_enumeration_required": False,
        "estimate_ball_invariant_exact_real": True,
        "representative_stacked_sector_ratios": ratios,
        "floating_point_rounding_enclosed": False,
        "source_dependent_coefficients_enclosed": False,
        "word_entry_set_qualified": False,
        "P4_MOTION_PASS": False,
        "P4_PASS": False,
        "P5_MAY_START": False,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--radius", type=float, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    report = build_report(args.radius)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+"\n")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
