#!/usr/bin/env python3
"""Innovation inverse enclosure retaining S = Q + R with Q >= 0, R > 0 provenance.

This diagnostic closes only the *invertibility/enclosure* defect created when an
SPD innovation family is replaced by an independent entrywise interval box.
It does not by itself close the Joseph/K dependency, P4, or P5.

For a physical innovation

    S = H P H^T + R = Q + R,

with P positive semidefinite and R uniformly positive definite, Q is PSD and
therefore S >= R > 0 in Loewner order.  Consequently det(S) >= det(R) > 0.
For a 3x3 innovation the adjugate formula then provides an outward inverse
box without admitting singular matrices that exist only in the rectangular
entrywise hull.  Hadamard's inequality supplies det(S) <= prod_i S_ii.

The ordinary validated interval inverse remains the preferred path whenever it
closes: this preserves its tighter dependency information.  PSD-plus-R is used
only after that path rejects an entrywise hull because a pivot contains zero.
The fallback cofactor numerators are still evaluated outwardly from the ordinary
S box; only the determinant denominator uses retained covariance provenance.
This is deliberately conservative and does not close the final same-history
S/K/P/H/R dependency.
"""
from __future__ import annotations

import json
import math
from typing import Sequence

from ou3_interval import Interval, IntervalMatrix, symmetric_positive_definite_ldlt
from ou3_interval_linear_algebra import (
    IntervalPivotError,
    matrix_inverse_gauss_jordan,
    matrix_symmetric_hull,
)

QUALIFICATION = "OU3_INNOVATION_PSD_PLUS_R_INVERSE_ENCLOSURE_V1"


def _shape(A: Sequence[Sequence[Interval]]) -> tuple[int, int]:
    r = len(A)
    c = len(A[0]) if r else 0
    if any(len(row) != c for row in A):
        raise ValueError("ragged interval matrix")
    return r, c


def _prod(xs: Sequence[float]) -> float:
    y = 1.0
    for x in xs:
        y *= float(x)
    return y


def _uniform_R_det_lower(R: Sequence[Sequence[Interval]]) -> tuple[float, list[Interval]]:
    Rbox = matrix_symmetric_hull(R)
    ok, pivots = symmetric_positive_definite_ldlt(Rbox)
    if not ok or len(pivots) != len(Rbox):
        raise ValueError("R must have a validated uniformly SPD interval LDLT")
    lows = [p.lo for p in pivots]
    if any((not math.isfinite(x)) or x <= 0.0 for x in lows):
        raise ValueError("R LDLT does not provide a positive determinant lower bound")
    # det(R)=prod D_i for every exact LDLT.  The outward pivot intervals enclose
    # each exact D_i, so the product of positive lower endpoints is a hard bound.
    return math.nextafter(_prod(lows), -math.inf), pivots


def innovation_inverse_psd_plus_R_3x3(
    S: Sequence[Sequence[Interval]],
    R: Sequence[Sequence[Interval]],
) -> tuple[IntervalMatrix, dict]:
    """Enclose S^-1 for S=Q+R, preferring the tighter generic proof.

    If the ordinary validated inverse rejects the rectangular hull, the caller
    must establish Q=H P H^T with actual P PSD on the same history.  Only then
    is the determinant lower bound from S>=R consumed by the fallback.
    """
    if _shape(S) != (3, 3) or _shape(R) != (3, 3):
        raise ValueError("3x3 innovation and R required")
    S = matrix_symmetric_hull(S)
    try:
        inv = matrix_inverse_gauss_jordan(S)
        return inv, {
            "inverse_mode": "GENERIC_VALIDATED_INTERVAL",
            "generic_validated_inverse_succeeded": True,
            "PSD_plus_R_fallback_used": False,
            "singular_entrywise_hull_members_are_not_admitted": True,
            "same_history_K_dependency_closed_here": False,
            "P4_PASS": False,
            "P5_MAY_START": False,
        }
    except IntervalPivotError as exc:
        generic_failure = str(exc)

    det_lo, r_pivots = _uniform_R_det_lower(R)
    diag_hi = [S[i][i].hi for i in range(3)]
    if any((not math.isfinite(x)) or x <= 0.0 for x in diag_hi):
        raise ValueError("S diagonal upper bounds must be finite positive")
    det_hi = math.nextafter(_prod(diag_hi), math.inf)
    if not (math.isfinite(det_hi) and det_hi >= det_lo > 0.0):
        raise ValueError("invalid determinant enclosure")
    det = Interval.outward_bounds(det_lo, det_hi)

    s00, s01, s02 = S[0][0], S[0][1], S[0][2]
    s11, s12, s22 = S[1][1], S[1][2], S[2][2]
    c00 = s11 * s22 - s12 * s12
    c11 = s00 * s22 - s02 * s02
    c22 = s00 * s11 - s01 * s01
    c01 = s02 * s12 - s01 * s22
    c02 = s01 * s12 - s02 * s11
    c12 = s01 * s02 - s00 * s12

    # Every actual S is SPD, hence its principal 2x2 minors are nonnegative.
    # Intersect only those three diagonal cofactors with that structural fact.
    def principal(c: Interval) -> Interval:
        if c.hi < 0.0:
            raise ValueError("cofactor hull contradicts retained SPD provenance")
        return Interval.outward_bounds(max(0.0, c.lo), c.hi)

    c00, c11, c22 = principal(c00), principal(c11), principal(c22)
    adj = [
        [c00, c01, c02],
        [c01, c11, c12],
        [c02, c12, c22],
    ]
    inv = [[adj[i][j] / det for j in range(3)] for i in range(3)]
    return matrix_symmetric_hull(inv), {
        "inverse_mode": "PSD_PLUS_R_3X3_FALLBACK",
        "generic_validated_inverse_succeeded": False,
        "generic_failure": generic_failure,
        "PSD_plus_R_fallback_used": True,
        "retained_identity": "S=H P H^T+R",
        "required_P_property": "P>=0 on the same physical/source history",
        "derived_Q_property": "H P H^T>=0",
        "uniform_R_SPD_validated": True,
        "R_ldlt_pivots": [p.as_list() for p in r_pivots],
        "det_S_lower_from_Loewner_monotonicity": det_lo,
        "det_S_upper_from_Hadamard": det_hi,
        "singular_entrywise_hull_members_are_not_admitted": True,
        "same_history_K_dependency_closed_here": False,
        "P4_PASS": False,
        "P5_MAY_START": False,
    }


def _point_inverse(A: list[list[float]]) -> list[list[float]]:
    a,b,c=A[0]; d,e,f=A[1]; g,h,i=A[2]
    det=a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)
    return [
        [(e*i-f*h)/det,(c*h-b*i)/det,(b*f-c*e)/det],
        [(f*g-d*i)/det,(a*i-c*g)/det,(c*d-a*f)/det],
        [(d*h-e*g)/det,(b*g-a*h)/det,(a*e-b*d)/det],
    ]


def build() -> dict:
    # A deliberately wide entrywise hull for R + q q^T, |q_i|<=1.  The box
    # contains singular/indefinite matrices although every retained member is SPD.
    r = 0.04
    R = [[Interval.point(r if i == j else 0.0) for j in range(3)] for i in range(3)]
    S = [[None for _ in range(3)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            if i == j:
                S[i][j] = Interval.outward_bounds(r, r + 1.0)
            else:
                S[i][j] = Interval.outward_bounds(-1.0, 1.0)
    generic_failed = False
    try:
        matrix_inverse_gauss_jordan(S)
    except IntervalPivotError:
        generic_failed = True
    inv, meta = innovation_inverse_psd_plus_R_3x3(S, R)

    contained = True
    for q in ((0.0,0.0,0.0),(1.0,1.0,1.0),(1.0,-1.0,0.5),(-0.3,0.7,-1.0)):
        A = [[(r if i == j else 0.0) + q[i]*q[j] for j in range(3)] for i in range(3)]
        Ai = _point_inverse(A)
        for i in range(3):
            for j in range(3):
                contained &= inv[i][j].lo <= Ai[i][j] <= inv[i][j].hi
    return {
        "qualification": QUALIFICATION,
        "generic_rectangular_inverse_rejected": generic_failed,
        "provenance_aware_inverse_finite": all(math.isfinite(x.lo) and math.isfinite(x.hi) for row in inv for x in row),
        "point_regression_inverses_contained": bool(contained),
        "inverse_box": [[x.as_list() for x in row] for row in inv],
        "proof": meta,
        "diagnostic_only": True,
        "P4_PASS": False,
        "P5_MAY_START": False,
    }


def validate(d: dict) -> list[str]:
    f=[]
    if not d.get("generic_rectangular_inverse_rejected"): f.append("pathological rectangular hull did not reproduce pivot obstruction")
    if not d.get("provenance_aware_inverse_finite"): f.append("provenance-aware inverse was not finite")
    if not d.get("point_regression_inverses_contained"): f.append("point inverse escaped enclosure")
    p=d.get("proof",{})
    if p.get("PSD_plus_R_fallback_used") is not True: f.append("pathological hull did not consume PSD-plus-R fallback")
    if p.get("uniform_R_SPD_validated") is not True: f.append("R SPD validation missing")
    if p.get("singular_entrywise_hull_members_are_not_admitted") is not True: f.append("singular hull classification missing")
    if p.get("same_history_K_dependency_closed_here") is not False: f.append("diagnostic overclaims K dependency closure")
    if d.get("P4_PASS") or d.get("P5_MAY_START"): f.append("diagnostic promoted P4/P5")
    return f


def main() -> int:
    d=build(); f=validate(d); d["validation_pass"]=not f; d["validation_failures"]=f
    print(json.dumps(d, indent=2, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
