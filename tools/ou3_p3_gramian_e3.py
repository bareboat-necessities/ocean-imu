#!/usr/bin/env python3
"""Sharp spectral conversion for the OU-III P3 translation Gramian.

The retained LTV process proof already supplies a rigorous lower bound on the
four-dimensional normalized Gramian determinant.  Converting that determinant
to a smallest eigenvalue with ``det(G)/trace(G)^3`` is valid, but catastrophically
loose for the integrator chain because the four directions carry very different
powers of the horizon.

For a positive semidefinite 4x4 matrix with ordered eigenvalues
``lambda_1 <= ... <= lambda_4``, let ``e3`` be the third elementary symmetric
polynomial.  Then

    det(G) / e3(G) = 1 / sum_i 1/lambda_i <= lambda_1.

Moreover ``e3(G)`` is the sum of the four 3x3 principal minors.  Hadamard's
inequality bounds each of those minors by the product of its three diagonal
entries.  For the damped [v,p,S,a_w] endpoint-response Gramian, nonnegative OU
damping can only reduce each response magnitude relative to the undamped chain,
so the unit-process diagonal entries obey

    G_vv <= H^3/3,
    G_pp <= H^5/20,
    G_SS <= H^7/252,
    G_aa <= H.

Combining these diagonal uppers with the already-validated determinant lower
therefore gives a rigorous generalized smallest-eigenvalue lower bound without
introducing a new source assumption or a numerical eigenvalue computation.
"""
from __future__ import annotations

import itertools
import math

from ou3_interval import Interval
import ou3_translational_uco_ucc as TRANS

SCHEMA = 1


def point(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def sharpen(
    upper: list[float],
    horizon_s: float,
    unit_gramian_det_lower: float,
    q_c_min_lower: float,
) -> dict:
    """Convert a determinant lower bound to a sharp generalized eigenvalue bound.

    ``upper`` is the positive diagonal covariance comparator in [v,p,S,a_w]
    order.  ``unit_gramian_det_lower`` is the determinant lower bound before
    multiplying by process intensity; ``q_c_min_lower`` is the source process
    intensity lower bound.
    """
    if len(upper) != 4:
        raise ValueError("four translation covariance uppers required")
    U = [float(x) for x in upper]
    if any(not (math.isfinite(x) and x > 0.0) for x in U):
        raise ValueError("finite positive covariance uppers required")
    H = float(horizon_s)
    det_unit = float(unit_gramian_det_lower)
    qc = float(q_c_min_lower)
    if not (math.isfinite(H) and H > 0.0):
        raise ValueError("positive finite horizon required")
    if not (math.isfinite(det_unit) and det_unit > 0.0):
        raise ValueError("strict unit-Gramian determinant lower required")
    if not (math.isfinite(qc) and qc > 0.0):
        raise ValueError("strict process-intensity lower required")

    Hiv = point(H)
    H3 = TRANS._pow_nonnegative(Hiv, 3)
    H5 = TRANS._pow_nonnegative(Hiv, 5)
    H7 = TRANS._pow_nonnegative(Hiv, 7)

    # Undamped endpoint-response diagonal upper bounds.  OU damping is
    # nonnegative throughout the certified source box, so every damped impulse
    # response magnitude is no larger than the corresponding integrator-chain
    # polynomial response.
    gramian_diag_upper = [
        (H3 / point(3.0)).hi,
        (H5 / point(20.0)).hi,
        (H7 / point(252.0)).hi,
        Hiv.hi,
    ]
    normalized_diag_upper = [
        (point(g) / point(u)).hi for g, u in zip(gramian_diag_upper, U)
    ]

    # e3(G) is the sum of all 3x3 principal minors.  G is PSD; Hadamard gives
    # det(G[J,J]) <= prod_{j in J} G_jj for every 3-index subset J.
    e3 = point(0.0)
    for comb in itertools.combinations(range(4), 3):
        term = point(1.0)
        for j in comb:
            term = term * point(normalized_diag_upper[j])
        e3 = e3 + term
    e3_upper = e3.hi

    product_U = point(1.0)
    for u in U:
        product_U = product_U * point(u)
    det_normalized_lower = (point(det_unit) / product_U).lo
    if not (det_normalized_lower > 0.0 and math.isfinite(e3_upper) and e3_upper > 0.0):
        raise RuntimeError("e3 spectral conversion lost strict positivity")

    lambda_lower = (point(det_normalized_lower) / point(e3_upper)).lo
    rho = (point(qc) * point(lambda_lower)).lo
    if not (lambda_lower > 0.0 and rho > 0.0):
        raise RuntimeError("e3 generalized Gramian lower is not strict")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_GRAMIAN_E3_GENERALIZED_EIGENVALUE_LOWER",
        "horizon_s": H,
        "Sigma_translation_diagonal_upper": U,
        "unit_gramian_det_lower": det_unit,
        "Sigma_normalized_gramian_det_lower": det_normalized_lower,
        "unit_gramian_diagonal_upper": gramian_diag_upper,
        "Sigma_normalized_gramian_diagonal_upper": normalized_diag_upper,
        "Sigma_normalized_gramian_e3_upper": e3_upper,
        "e3_identity": "det(G)/e3(G)=1/sum_i(1/lambda_i)<=lambda_min(G)",
        "principal_minor_upper_route": "PSD_HADAMARD_DIAGONAL_PRODUCT",
        "nonnegative_damping_response_upper_route": "UNDAMPED_INTEGRATOR_CHAIN",
        "Sigma_normalized_unit_gramian_lambda_min_lower": lambda_lower,
        "q_c_min_lower": qc,
        "relative_process_floor_lower": rho,
        "validated_interval_arithmetic": True,
        "numerical_eigendecomposition_used": False,
    }


def sharpen_probe(probe: dict, upper: list[float]) -> dict:
    """Return a copy of one retained LTV probe with the e3 spectral bound."""
    d = sharpen(
        upper,
        float(probe["horizon_s"]),
        float(probe["unit_gramian_det_lower"]),
        float(probe["q_c_min_lower"]),
    )
    out = dict(probe)
    old = float(probe.get("relative_process_floor_lower", 0.0))
    out.update(d)
    out["trace3_relative_process_floor_lower_baseline"] = old
    out["improvement_over_trace3"] = math.inf if old <= 0.0 else float(out["relative_process_floor_lower"]) / old
    return out


def validate(d: dict) -> list[str]:
    f = []
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P3_GRAMIAN_E3_GENERALIZED_EIGENVALUE_LOWER":
        f.append("wrong qualification")
    for key in ("validated_interval_arithmetic",):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    if d.get("numerical_eigendecomposition_used") is not False:
        f.append("numerical eigendecomposition was used")
    for key in (
        "Sigma_normalized_gramian_det_lower",
        "Sigma_normalized_gramian_e3_upper",
        "Sigma_normalized_unit_gramian_lambda_min_lower",
        "q_c_min_lower",
        "relative_process_floor_lower",
    ):
        x = d.get(key)
        if not isinstance(x, (int, float)) or not (math.isfinite(float(x)) and float(x) > 0.0):
            f.append(f"{key} is not strict finite positive")
    return list(dict.fromkeys(f))
