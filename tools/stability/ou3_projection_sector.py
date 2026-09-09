#!/usr/bin/env python3
"""Uniform Euclidean sector for the shipping radial accelerometer-bias clamp.

This is an analytical P4 ingredient, not a replay fit. With true physical bias
beta, pre-projection corrected error e and projection radius R, shipping
projection is represented in error coordinates by

    F_R(e,beta) = beta - Pi_R(beta-e).

For every pair (e1,beta1),(e2,beta2), including points on different sides of
the saturation boundary,

    ||F_R(e1,beta1)-F_R(e2,beta2)||^2
       <= ||e1-e2||^2 + ||beta1-beta2||^2.

Away from the boundary the Euclidean projection Jacobian J is symmetric with
0<=J<=I.  At the boundary every Clarke generalized Jacobian is a convex
combination of the one-sided limits and has the same property.  Therefore

    dF = J de + (I-J) dbeta,
    [J,I-J][J,I-J]^T = J^2+(I-J)^2 <= I.

The proof is analytical. The helpers below are dependency-free regression
checks only; the retained stability package intentionally does not require
numpy merely to import its theorem modules.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def _vec(value):
    out=[float(x) for x in value]
    if not out or any(not math.isfinite(x) for x in out):
        raise ValueError("finite nonempty vector required")
    return out


def _norm(x):
    return math.sqrt(sum(v*v for v in x))


def project_ball(value, radius):
    value=_vec(value)
    if not math.isfinite(radius) or radius < 0:
        raise ValueError("finite nonnegative radius required")
    if radius == 0:
        return [0.0 for _ in value]
    scale=max(abs(x) for x in value)
    if scale == 0:
        return list(value)
    direction=[x/scale for x in value]
    length=_norm(direction)
    if scale <= radius/length:
        return list(value)
    factor=radius/(scale*length)
    return [x*factor for x in value]


def bias_error_projection(error_pre, beta, radius):
    error_pre=_vec(error_pre)
    beta=_vec(beta)
    if len(error_pre) != len(beta):
        raise ValueError("error and true bias must have the same shape")
    estimate=[b-e for b,e in zip(beta,error_pre)]
    projected=project_ball(estimate,radius)
    return [b-p for b,p in zip(beta,projected)]


def _eye(n):
    return [[1.0 if i==j else 0.0 for j in range(n)] for i in range(n)]


def projection_jacobian(estimate_pre, radius):
    """One classical Jacobian away from ||estimate||=R; boundary is Clarke."""
    x=_vec(estimate_pre)
    if not math.isfinite(radius) or radius < 0:
        raise ValueError("finite nonnegative radius required")
    n=len(x)
    length=_norm(x)
    if radius == 0:
        return [[0.0]*n for _ in range(n)]
    if length < radius:
        return _eye(n)
    if length > radius:
        u=[v/length for v in x]
        a=radius/length
        return [[a*((1.0 if i==j else 0.0)-u[i]*u[j]) for j in range(n)] for i in range(n)]
    raise ValueError("classical Jacobian is set-valued on projection boundary")


def _det3(a):
    return (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
            -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
            +a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))


def _psd3(a, tol=2e-14):
    """Regression-only PSD test by the 3x3 principal-minor criterion."""
    if len(a)!=3 or any(len(row)!=3 for row in a):
        raise ValueError("projection regression matrix must be 3x3")
    for i in range(3):
        for j in range(3):
            if not math.isfinite(float(a[i][j])):
                return False
            if abs(float(a[i][j])-float(a[j][i])) > tol:
                return False
    scale=max(1.0,max(abs(float(x)) for row in a for x in row))
    eps=tol*scale*scale*scale
    if any(float(a[i][i]) < -eps for i in range(3)):
        return False
    for i,j in ((0,1),(0,2),(1,2)):
        minor=float(a[i][i])*float(a[j][j])-float(a[i][j])*float(a[j][i])
        if minor < -eps:
            return False
    return _det3(a) >= -eps


def stacked_sector_ratio(jacobian):
    """Certified analytical upper for ||[J,I-J]||_2^2.

    The theorem proves this is at most one whenever J and I-J are PSD.  This
    helper checks that precondition for the 3D regression matrix by principal
    minors.  Unlike Gershgorin, that check does not falsely reject valid
    off-axis rank-two radial projection Jacobians.
    """
    j=[[float(x) for x in row] for row in jacobian]
    if len(j)!=3 or any(len(row)!=3 for row in j):
        raise ValueError("projection Jacobian must be 3x3")
    I=_eye(3)
    imj=[[I[i][k]-j[i][k] for k in range(3)] for i in range(3)]
    if not _psd3(j) or not _psd3(imj):
        raise ValueError("representative Jacobian lost 0<=J<=I")
    return 1.0


def build_report(radius):
    if not math.isfinite(radius) or radius <= 0:
        raise ValueError("positive finite radius required")
    inside=projection_jacobian([0.25*radius,0.0,0.0],radius)
    outside=projection_jacobian([2.0*radius,0.0,0.0],radius)
    tangent=[[0.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]
    I=_eye(3)
    boundary_mid=[[0.5*I[i][j]+0.5*tangent[i][j] for j in range(3)] for i in range(3)]
    ratios={
        "inside":stacked_sector_ratio(inside),
        "outside":stacked_sector_ratio(outside),
        "boundary_clarke_midpoint":stacked_sector_ratio(boundary_mid),
    }
    return {
        "lemma":"OU3_RADIAL_BIAS_PROJECTION_JOINT_SECTOR_V1",
        "map":"F_R(e,beta)=beta-Pi_R(beta-e)",
        "sector":"||Delta F||^2 <= ||Delta e||^2 + ||Delta beta||^2",
        "projection_radius":radius,
        "exact_real_operator_sector_closed":True,
        "saturated_branch_included_analytically":True,
        "unsaturated_branch_included_analytically":True,
        "boundary_clarke_branch_included_analytically":True,
        "fixed_multiplier_required":False,
        "saturation_pattern_enumeration_required":False,
        "estimate_ball_invariant_exact_real":True,
        "representative_stacked_sector_ratios":ratios,
        "floating_point_rounding_enclosed":False,
        "source_dependent_coefficients_enclosed":False,
        "word_entry_set_qualified":False,
        "P4_MOTION_PASS":False,
        "P4_PASS":False,
        "P5_MAY_START":False,
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--radius",type=float,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    report=build_report(args.radius)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps(report,sort_keys=True))


if __name__=="__main__":
    main()
