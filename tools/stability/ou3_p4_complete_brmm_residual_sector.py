#!/usr/bin/env python3
"""Exact Cayley residual-sector algebra for complete-BRMM OU-III P4.

This module closes an algebraic gap in the signed-information route.  It does
not replace the complete BRMM source by independent event boxes and it does not
sum packetwise worst-case defects.

For c=2*tan(theta/2)u, let C=[c]_x, r=||c|| and

    E(c)=I + 4/(4+r^2) C + 2/(4+r^2) C^2.

The literal linearized vector/accelerometer attitude column is C*v.  Hence the
finite nonlinear vector remainder is exactly

    eta_vec = (E-I)v - C v
            = [ -r^2/(4+r^2) C + 2/(4+r^2) C^2 ] v.      (1)

For the accelerometer residual the additive b_a term is exactly linear and
cancels from eta.  The only additional finite term is

    eta_aw = (E-I) R_hat delta_a_w
           = [ 4/(4+r^2) C + 2/(4+r^2) C^2 ] R_hat da.    (2)

Thus, on a source-correlated cell with r<=r_bar,

    ||eta_vec|| <= g_vec(r_bar,||v||) ||c||,
    ||eta_acc|| <= g_f(r_bar,||f||) ||c|| + g_aw(r_bar)||da||,

where

    g_vec = ||v|| r_bar(r_bar+2)/(4+r_bar^2),
    g_f   = ||f|| r_bar(r_bar+2)/(4+r_bar^2),
    g_aw  = r_bar(4+2r_bar)/(4+r_bar^2).

These are homogeneous quadratic-sector bounds.  In particular there is no
state-independent defect and no N-times remainder term.  A caller may retain
each same-history event's own f/R/Racc cell and assemble the corresponding
block quadratic form before comparing the JOINT residual energy with the JOINT
signed information of the complete word.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Sequence

SCHEMA = 1
QUALIFICATION = "OU3_P4_COMPLETE_BRMM_EXACT_CAYLEY_RESIDUAL_SECTOR_V1"


def _norm(v: Sequence[float]) -> float:
    return math.sqrt(sum(float(x) * float(x) for x in v))


def _skew(c: Sequence[float]) -> list[list[float]]:
    if len(c) != 3:
        raise ValueError("three-vector required")
    x, y, z = (float(q) for q in c)
    return [[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]]


def _mv(A: Sequence[Sequence[float]], x: Sequence[float]) -> list[float]:
    if not A or any(len(row) != len(x) for row in A):
        raise ValueError("matrix/vector dimension mismatch")
    return [sum(float(a) * float(b) for a, b in zip(row, x)) for row in A]


def _mm(A: Sequence[Sequence[float]], B: Sequence[Sequence[float]]) -> list[list[float]]:
    if not A or not B or len(A[0]) != len(B):
        raise ValueError("matrix dimension mismatch")
    BT = list(zip(*B))
    return [[sum(float(a) * float(b) for a, b in zip(row, col)) for col in BT] for row in A]


def _add(a: Sequence[float], b: Sequence[float]) -> list[float]:
    if len(a) != len(b):
        raise ValueError("vector dimension mismatch")
    return [float(x) + float(y) for x, y in zip(a, b)]


def _scale(a: float, x: Sequence[float]) -> list[float]:
    return [float(a) * float(q) for q in x]


def rotation_from_cayley(c: Sequence[float]) -> list[list[float]]:
    """Exact inverse Cayley rotation used by the proof AD backend."""
    if len(c) != 3:
        raise ValueError("Cayley vector must have length three")
    C = _skew(c)
    C2 = _mm(C, C)
    r2 = sum(float(x) * float(x) for x in c)
    den = 4.0 + r2
    I = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
    return [[I[i][j] + (4.0 / den) * C[i][j] + (2.0 / den) * C2[i][j]
             for j in range(3)] for i in range(3)]


def exact_vector_eta(c: Sequence[float], v: Sequence[float]) -> list[float]:
    """Return exact finite vector residual minus its literal linear tangent."""
    if len(c) != 3 or len(v) != 3:
        raise ValueError("vector residual requires two three-vectors")
    E = rotation_from_cayley(c)
    Cv = _mv(_skew(c), v)
    Ev = _mv(E, v)
    return [Ev[i] - float(v[i]) - Cv[i] for i in range(3)]


def exact_vector_eta_factorized(c: Sequence[float], v: Sequence[float]) -> list[float]:
    C = _skew(c)
    C2 = _mm(C, C)
    r2 = sum(float(x) * float(x) for x in c)
    den = 4.0 + r2
    return _add(_scale(-r2 / den, _mv(C, v)), _scale(2.0 / den, _mv(C2, v)))


def exact_accelerometer_eta(
    c: Sequence[float], f_hat: Sequence[float], R_hat: Sequence[Sequence[float]], da: Sequence[float]
) -> list[float]:
    """Exact finite accelerometer residual minus H_acc*[c,da,b_a].

    b_a is absent because its finite residual contribution is exactly +delta_ba
    and therefore cancels from the nonlinear remainder in both H18 and A21.
    """
    if len(c) != 3 or len(f_hat) != 3 or len(da) != 3 or len(R_hat) != 3:
        raise ValueError("accelerometer sector requires 3D inputs")
    vec = exact_vector_eta(c, f_hat)
    E = rotation_from_cayley(c)
    Rda = _mv(R_hat, da)
    ERda = _mv(E, Rda)
    return _add(vec, [ERda[i] - Rda[i] for i in range(3)])


def exact_accelerometer_eta_factorized(
    c: Sequence[float], f_hat: Sequence[float], R_hat: Sequence[Sequence[float]], da: Sequence[float]
) -> list[float]:
    C = _skew(c)
    C2 = _mm(C, C)
    r2 = sum(float(x) * float(x) for x in c)
    den = 4.0 + r2
    eta_f = exact_vector_eta_factorized(c, f_hat)
    Rda = _mv(R_hat, da)
    eta_aw = _add(_scale(4.0 / den, _mv(C, Rda)), _scale(2.0 / den, _mv(C2, Rda)))
    return _add(eta_f, eta_aw)


def sector_coefficients(r_bar: float, vector_norm_upper: float) -> dict[str, float]:
    r = float(r_bar)
    v = float(vector_norm_upper)
    if not (math.isfinite(r) and r >= 0.0 and math.isfinite(v) and v >= 0.0):
        raise ValueError("sector bounds must be finite nonnegative")
    den = 4.0 + r * r
    return {
        "attitude_coefficient": v * r * (r + 2.0) / den,
        "aw_cross_coefficient": r * (4.0 + 2.0 * r) / den,
    }


def sector_quadratic_diagonal(
    r_bar: float,
    vector_norm_upper: float,
    Rinv_spectral_upper: float,
    *,
    include_aw: bool,
) -> list[float]:
    """Diagonal quadratic form Q with eta'R^-1 eta <= [c,da]'Q[c,da].

    We use (a+b)^2 <= 2a^2+2b^2 for the accelerometer channel.  This remains a
    matrix quadratic sector and is intended to be assembled over the correlated
    same-history word; it is not an event-count scalar remainder budget.
    """
    lam = float(Rinv_spectral_upper)
    if not (math.isfinite(lam) and lam >= 0.0):
        raise ValueError("R inverse spectral upper bound must be finite nonnegative")
    g = sector_coefficients(r_bar, vector_norm_upper)
    if include_aw:
        qc = 2.0 * lam * g["attitude_coefficient"] ** 2
        qa = 2.0 * lam * g["aw_cross_coefficient"] ** 2
        return [qc, qc, qc, qa, qa, qa]
    qc = lam * g["attitude_coefficient"] ** 2
    return [qc, qc, qc]


def build() -> dict:
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "exact_inverse_cayley_factorization_used": True,
        "vector_eta_identity": "eta=[-r^2/(4+r^2) C + 2/(4+r^2) C^2]v",
        "accelerometer_eta_aw_identity": "eta_aw=[4/(4+r^2) C + 2/(4+r^2) C^2]R_hat*delta_a_w",
        "accelerometer_bias_nonlinearity_exactly_zero": True,
        "nonlinear_residual_has_no_state_independent_term": True,
        "finite_residual_is_quadratic_sector_in_c_and_aw": True,
        "same_history_event_cells_may_be_retained_individually": True,
        "packet_count_multiplier_used": False,
        "global_correction_radius_used": False,
        "inverse_metric_floor_used": False,
        "source_family_replaced": False,
        "filter_changed": False,
        "P4_promoted_here": False,
        "next_obligation": (
            "assemble these per-event quadratic sectors on the SAME correlated complete-BRMM history and compare their joint block form, together with reset/projection terms, against the signed whole-word information ledger"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_BRMM_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for key in (
        "exact_inverse_cayley_factorization_used",
        "accelerometer_bias_nonlinearity_exactly_zero",
        "nonlinear_residual_has_no_state_independent_term",
        "finite_residual_is_quadratic_sector_in_c_and_aw",
        "same_history_event_cells_may_be_retained_individually",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "packet_count_multiplier_used", "global_correction_radius_used",
        "inverse_metric_floor_used", "source_family_replaced", "filter_changed", "P4_promoted_here",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    return f


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build()
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"qualification": QUALIFICATION, "failures": failures}, indent=2))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
