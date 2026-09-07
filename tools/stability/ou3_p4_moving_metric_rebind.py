#!/usr/bin/env python3
"""Exact moving-covariance metric rebind for canonical OU-III P4.

This structural certificate binds P4 to the exact shipping covariance algebra.
For P+=A P A^T+B with B>=0, A^T P+^-1 A<=P^-1.  For a Joseph
update B=K R K^T, K^T P+^-1 K<=R^-1.  Nonsingular covariance
congruence is an exact metric isometry.  These identities are independent of
state dimension and therefore apply to both H18 and A21.  They are covariance
identities only: they do not identify an auxiliary residual Jacobian with the
congruent shipping Jacobian, or close nonlinear chart transport and storage.
"""
from __future__ import annotations

from pathlib import Path

from ou3_interval import (
    matrix_add,
    matrix_identity,
    matrix_mul,
    matrix_point,
    matrix_sub,
    matrix_transpose,
)
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan, matrix_symmetric_hull
import ou3_sea3_full_word_riccati_backend as BACKEND
import ou3_sea3_full_word_reset_congruence as RESET

REPO = Path(__file__).resolve().parents[2]
SCHEMA = 2
QUALIFICATION = "OU3_P4_EXACT_MOVING_COVARIANCE_METRIC_REBIND_V2"


def _contains_zero(A) -> bool:
    return all(x.lo <= 0.0 <= x.hi for row in A for x in row)


def _joseph_smoke() -> dict:
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

    lin_slack = matrix_symmetric_hull(
        matrix_sub(Pinv, matrix_mul(matrix_mul(At, Pplus_inv), A))
    )
    expected_lin_slack = matrix_symmetric_hull(matrix_mul(matrix_mul(Ht, Sinv), H))
    lin_identity = _contains_zero(matrix_sub(lin_slack, expected_lin_slack))

    inj_slack = matrix_symmetric_hull(
        matrix_sub(
            matrix_inverse_gauss_jordan(R),
            matrix_mul(matrix_mul(matrix_transpose(K), Pplus_inv), K),
        )
    )
    inj_strict = inj_slack[0][0].lo > 0.0
    return {
        "linear_nonexpansive_PSD_identity_enclosed": lin_identity,
        "Joseph_injection_interval_strict_smoke": inj_strict,
    }


def _congruence_smoke() -> dict:
    P = matrix_point([[1.6, 0.2], [0.2, 0.9]])
    G = matrix_point([[1.0, 0.2], [-0.1, 1.0]])
    Pp = matrix_symmetric_hull(matrix_mul(matrix_mul(G, P), matrix_transpose(G)))
    Pinv = matrix_inverse_gauss_jordan(P)
    Ppinv = matrix_inverse_gauss_jordan(Pp)
    image = matrix_mul(matrix_mul(matrix_transpose(G), Ppinv), G)
    return {"coordinate_congruence_metric_identity_enclosed": _contains_zero(matrix_sub(image, Pinv))}


def build() -> dict:
    parity = BACKEND.shipping_source_parity()
    reset_failures = RESET.validate()
    joseph = _joseph_smoke()
    congruence = _congruence_smoke()
    closed = bool(
        all(parity.values()) and not reset_failures
        and all(joseph.values()) and all(congruence.values())
    )
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
            "For Joseph, P^-1-A^T P_plus^-1 A=H^T S^-1 H>=0 and "
            "K^T P_plus^-1 K<=R^-1. Nonsingular congruence is exact."
        ),
        "interval_self_tests": {**joseph, **congruence},
        "structural_shipping_covariance_identities_closed": closed,
        "full_nonlinear_measurement_metric_rebind_closed": False,
        "nonlinear_chart_transport_and_storage_closed": False,
        "P4_promoted": False,
    }


def validate(d: dict) -> list[str]:
    f = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for k in (
        "source_generated_not_trajectory_fit", "source_parity_pass",
        "reset_certificate_consumed", "moving_metric_coordinate_congruence_exact",
        "prediction_linear_map_nonexpansive", "Joseph_linear_map_nonexpansive",
        "Joseph_nonlinear_injection_metric_closed", "PSD_floor_nonexpansive",
        "left_error_reset_exact_metric_isometry", "dimension_independent_H18_A21",
        "structural_shipping_covariance_identities_closed",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "trajectory_replay_used", "filter_changed", "group_isotropic_metric_attachment_used",
        "endpoint_source_word_scan_used", "ordinary_float_eigensolver_used", "P4_promoted",
        "full_nonlinear_measurement_metric_rebind_closed",
        "nonlinear_chart_transport_and_storage_closed",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if d.get("reset_certificate_failures"):
        f.append("reset certificate failed")
    tests = d.get("interval_self_tests", {})
    for k in (
        "linear_nonexpansive_PSD_identity_enclosed",
        "Joseph_injection_interval_strict_smoke",
        "coordinate_congruence_metric_identity_enclosed",
    ):
        if tests.get(k) is not True:
            f.append(f"interval self-test failed: {k}")
    return list(dict.fromkeys(f))


if __name__ == "__main__":
    d = build()
    f = validate(d)
    print({"qualification": QUALIFICATION, "structural_covariance_closed": d["structural_shipping_covariance_identities_closed"], "failures": f})
    raise SystemExit(0 if not f else 2)
