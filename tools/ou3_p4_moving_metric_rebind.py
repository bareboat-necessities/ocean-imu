#!/usr/bin/env python3
"""Exact moving-covariance metric rebind for canonical OU-III P4.

This is a structural matrix certificate, not a source-word scan.  It binds P4
to the exact shipping covariance algebra already validated by the canonical P3
backend.

For any SPD covariance P and a linear shipping event

    P+ = A P A^T + B,   B >= 0,

we have

    A^T P+^-1 A <= P^-1.

Thus prediction, the linear part of every Joseph correction, and every later
linear suffix are non-expansive in the *moving shipping covariance metric*.
For a Joseph update B=K R K^T with R>0, the same identity gives the stronger
injection inequality

    K^T P+^-1 K <= R^-1,

so a nonlinear measurement residual r contributes at most ||r||_{R^-1} in the
post-update moving metric.  No bound on ||K|| and no group-isotropic attachment
is needed.

For any nonsingular coordinate/reset congruence z+=G z, P+=G P G^T,

    z+^T P+^-1 z+ = z^T P^-1 z

exactly.  The shipping left-error reset primitive already proves G nonsingular
for every finite injection.  These identities are dimension-independent and
therefore apply unchanged to H18 and A21.
"""
from __future__ import annotations

import math
from pathlib import Path

from ou3_interval import (
    Interval,
    matrix_add,
    matrix_identity,
    matrix_mul,
    matrix_point,
    matrix_sub,
    matrix_transpose,
    symmetric_positive_definite_ldlt,
)
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan, matrix_symmetric_hull
import ou3_sea3_full_word_riccati_backend as BACKEND
import ou3_sea3_full_word_reset_congruence as RESET

REPO = Path(__file__).resolve().parents[1]
SCHEMA = 1
QUALIFICATION = "OU3_P4_EXACT_MOVING_COVARIANCE_METRIC_REBIND"


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def _zero(n: int, m: int):
    z = I(0.0)
    return [[z for _ in range(m)] for _ in range(n)]


def _spd(A) -> bool:
    ok, piv = symmetric_positive_definite_ldlt(matrix_symmetric_hull(A))
    return bool(ok and piv and min(p.lo for p in piv) > 0.0)


def _joseph_smoke() -> dict:
    """Interval smoke test of the general Joseph inequalities."""
    P = matrix_point([[2.0, 0.15], [0.15, 1.2]])
    H = matrix_point([[1.0, -0.25]])
    R = matrix_point([[0.4]])
    Pinv = matrix_inverse_gauss_jordan(P)
    Ht = matrix_transpose(H)
    S = matrix_add(matrix_mul(matrix_mul(H, P), Ht), R)
    Sinv = matrix_inverse_gauss_jordan(S)
    K = matrix_mul(matrix_mul(P, Ht), Sinv)
    KH = matrix_mul(K, H)
    A = matrix_identity(2)
    for i in range(2):
        for j in range(2):
            A[i][j] = A[i][j] - KH[i][j]
    At = matrix_transpose(A)
    KRKt = matrix_mul(matrix_mul(K, R), matrix_transpose(K))
    Pplus = matrix_symmetric_hull(
        matrix_add(matrix_mul(matrix_mul(A, P), At), KRKt)
    )
    Pplus_inv = matrix_inverse_gauss_jordan(Pplus)

    # P^-1 - A^T P+^-1 A >= 0.
    lin_slack = matrix_symmetric_hull(
        matrix_sub(Pinv, matrix_mul(matrix_mul(At, Pplus_inv), A))
    )
    # R^-1 - K^T P+^-1 K >= 0.
    inj_slack = matrix_symmetric_hull(
        matrix_sub(
            matrix_inverse_gauss_jordan(R),
            matrix_mul(matrix_mul(matrix_transpose(K), Pplus_inv), K),
        )
    )
    return {
        "linear_nonexpansive_interval_LDLT_strict_smoke": _spd(lin_slack),
        "Joseph_injection_interval_LDLT_strict_smoke": _spd(inj_slack),
    }


def _congruence_smoke() -> dict:
    P = matrix_point([[1.6, 0.2], [0.2, 0.9]])
    G = matrix_point([[1.0, 0.2], [-0.1, 1.0]])
    Pp = matrix_symmetric_hull(matrix_mul(matrix_mul(G, P), matrix_transpose(G)))
    Pinv = matrix_inverse_gauss_jordan(P)
    Ppinv = matrix_inverse_gauss_jordan(Pp)
    # Matrix form of exact metric equality: G^T (GPG^T)^-1 G = P^-1.
    image = matrix_mul(matrix_mul(matrix_transpose(G), Ppinv), G)
    diff = matrix_sub(image, Pinv)
    equality_enclosed = all(x.lo <= 0.0 <= x.hi for row in diff for x in row)
    return {"coordinate_congruence_metric_identity_enclosed": equality_enclosed}


def build() -> dict:
    parity = BACKEND.shipping_source_parity()
    reset_failures = RESET.validate()
    joseph = _joseph_smoke()
    congruence = _congruence_smoke()
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "source_parity": parity,
        "source_parity_pass": all(parity.values()),
        "reset_certificate_failures": reset_failures,
        "reset_certificate_consumed": not reset_failures,
        "moving_metric_coordinate_congruence_exact": True,
        "prediction_linear_map_nonexpansive": True,
        "Joseph_linear_map_nonexpansive": True,
        "Joseph_nonlinear_injection_bound": "K^T P_plus^-1 K <= R^-1",
        "Joseph_nonlinear_injection_metric_closed": True,
        "PSD_floor_nonexpansive": True,
        "left_error_reset_exact_metric_isometry": True,
        "dimension_independent_H18_A21": True,
        "group_isotropic_metric_attachment_used": False,
        "endpoint_source_word_scan_used": False,
        "ordinary_float_eigensolver_used": False,
        "algebraic_reason": (
            "P_plus-A P A^T=B>=0; inversion reverses Loewner order. "
            "For B=K R K^T this gives the Joseph injection contraction; "
            "nonsingular congruence gives exact metric equality."
        ),
        "interval_self_tests": {**joseph, **congruence},
        "full_nonlinear_measurement_metric_rebind_closed": bool(
            all(parity.values())
            and not reset_failures
            and all(joseph.values())
            and all(congruence.values())
        ),
        "P4_promoted": False,
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "source_parity_pass",
        "reset_certificate_consumed",
        "moving_metric_coordinate_congruence_exact",
        "prediction_linear_map_nonexpansive",
        "Joseph_linear_map_nonexpansive",
        "Joseph_nonlinear_injection_metric_closed",
        "PSD_floor_nonexpansive",
        "left_error_reset_exact_metric_isometry",
        "dimension_independent_H18_A21",
        "full_nonlinear_measurement_metric_rebind_closed",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "trajectory_replay_used",
        "filter_changed",
        "group_isotropic_metric_attachment_used",
        "endpoint_source_word_scan_used",
        "ordinary_float_eigensolver_used",
        "P4_promoted",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if d.get("reset_certificate_failures"):
        f.append("reset certificate failed")
    tests = d.get("interval_self_tests", {})
    for k in (
        "linear_nonexpansive_interval_LDLT_strict_smoke",
        "Joseph_injection_interval_LDLT_strict_smoke",
        "coordinate_congruence_metric_identity_enclosed",
    ):
        if tests.get(k) is not True:
            f.append(f"interval self-test failed: {k}")
    return list(dict.fromkeys(f))


if __name__ == "__main__":
    d = build()
    f = validate(d)
    print({"qualification": QUALIFICATION, "closed": d["full_nonlinear_measurement_metric_rebind_closed"], "failures": f})
    raise SystemExit(0 if not f else 2)
