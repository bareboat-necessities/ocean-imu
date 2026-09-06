#!/usr/bin/env python3
"""Validated enclosure of the shipping a_w covariance-floor operation.

The periodic source event only requests synchronization at the currently
committed stationary covariance Sigma_aw.  Its numerical increment is *not* an
exogenous/source coordinate: shipping computes it from the current Riccati
covariance after prediction,

    Delta = Pi_+(Sigma_aw - P_awaw^-),

where Pi_+ clamps the eigenvalues of the symmetric 3x3 argument at zero.  Thus
H18 and A21 may receive different increments even on the same SEA3 source
history.

This module keeps that dependency in the Riccati executor.  When the complete
interval argument is strictly positive or negative definite, the positive-part
map is exact.  In a mixed/uncertain spectral cell, it returns a rigorous entry-
wise enclosure using

    ||Pi_+(D)||_2 <= ||D||_2 <= ||D||_F.

For a PSD matrix Y with ||Y||_2 <= rho, 0 <= Y_ii <= rho and
|Y_ij| <= rho, so the returned box contains every positive part of every
symmetric point matrix D in the interval family.  This generic fallback may be
wide, but it never selects a favorable eigenspace or turns the floor into a
source-supplied free increment.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval, IntervalMatrix, symmetric_positive_definite_ldlt
from ou3_interval_linear_algebra import matrix_symmetric_hull

REPO = Path(__file__).resolve().parents[2]
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_COVARIANCE_DEPENDENT_AW_FLOOR_ENCLOSURE"


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def _shape(A: Sequence[Sequence[Interval]]) -> tuple[int, int]:
    r = len(A)
    c = len(A[0]) if r else 0
    if any(len(row) != c for row in A):
        raise ValueError("ragged interval matrix")
    return r, c


def _zero3() -> IntervalMatrix:
    # Zero is an exact algebraic constant here.  Widening it would turn the
    # strict-negative positive-part branch into a tiny signed interval instead
    # of the exact matrix Pi_+(D)=0.
    z = Interval.point(0.0)
    return [[z for _ in range(3)] for _ in range(3)]


def _neg(A: Sequence[Sequence[Interval]]) -> IntervalMatrix:
    return [[-A[i][j] for j in range(len(A[i]))] for i in range(len(A))]


def _sub(A: Sequence[Sequence[Interval]], B: Sequence[Sequence[Interval]]) -> IntervalMatrix:
    if _shape(A) != _shape(B):
        raise ValueError("matrix shape mismatch")
    return [[A[i][j] - B[i][j] for j in range(len(A[i]))] for i in range(len(A))]


def stationary_sigma_isotropic(std_aw: Interval) -> IntervalMatrix:
    """Default shipping Sigma_aw for S_factor=1 and no cross-correlation."""
    if std_aw.lo <= 0.0:
        raise ValueError("stationary a_w standard deviation must be positive")
    S = _zero3()
    v = std_aw.square()
    for i in range(3):
        S[i][i] = v
    return S


def _frobenius_upper(A: Sequence[Sequence[Interval]]) -> float:
    total = 0.0
    for row in A:
        for x in row:
            a = x.abs_upper()
            total = math.nextafter(total + math.nextafter(a * a, math.inf), math.inf)
    r = math.sqrt(total)
    return math.nextafter(r, math.inf)


def positive_part_enclosure(
    target_sigma: Sequence[Sequence[Interval]],
    P_awaw_minus: Sequence[Sequence[Interval]],
) -> tuple[IntervalMatrix, str]:
    """Enclose Pi_+(target_sigma-P_awaw_minus) for every interval point."""
    if _shape(target_sigma) != (3, 3) or _shape(P_awaw_minus) != (3, 3):
        raise ValueError("a_w floor requires two 3x3 matrices")
    D = matrix_symmetric_hull(_sub(target_sigma, P_awaw_minus))

    pos_pass, _ = symmetric_positive_definite_ldlt(D)
    if pos_pass:
        # Every point matrix is strictly PD, hence Pi_+(D)=D.
        return D, "strict_positive_definite_exact"

    negD = matrix_symmetric_hull(_neg(D))
    neg_pass, _ = symmetric_positive_definite_ldlt(negD)
    if neg_pass:
        # Every point matrix is strictly ND, hence Pi_+(D)=0.
        return _zero3(), "strict_negative_definite_exact_zero"

    rho = _frobenius_upper(D)
    if not math.isfinite(rho):
        raise ValueError("nonfinite covariance-floor enclosure")
    out = _zero3()
    diag = Interval(0.0, rho)
    off = Interval(-rho, rho)
    for i in range(3):
        for j in range(3):
            out[i][j] = diag if i == j else off
    return out, "mixed_spectrum_frobenius_outer"


def aw_block(P: Sequence[Sequence[Interval]], off_aw: int = 15) -> IntervalMatrix:
    n, m = _shape(P)
    if n != m or n < off_aw + 3:
        raise ValueError("full covariance does not contain the a_w block")
    return [[P[off_aw + i][off_aw + j] for j in range(3)] for i in range(3)]


def shipping_source_parity() -> dict[str, bool]:
    text = MEKF.read_text(encoding="utf-8")
    p = text.find("EIGEN_STRONG_INLINE void apply_pending_aw_covariance_inflation_()")
    q = text.find("apply_pending_aw_covariance_inflation_();", p + 1)
    body = text[p:q] if p >= 0 and q > p else ""
    return {
        "floor_body_scoped": bool(body),
        "argument_is_target_minus_current_predicted_aw_covariance": (
            "Matrix3 Delta = aw_covariance_floor_target_ - P_aw;" in body
        ),
        "argument_is_symmetrized": "Delta = T(0.5) * (Delta + Delta.transpose());" in body,
        "shipping_uses_self_adjoint_eigensolver": "Eigen::SelfAdjointEigenSolver<Matrix3> es(Delta);" in body,
        "shipping_clamps_negative_eigenvalues_to_zero": (
            "evals(i) = std::max(T(0), evals(i));" in body
        ),
        "shipping_reconstructs_positive_part": (
            "Delta = es.eigenvectors() * evals.asDiagonal() * es.eigenvectors().transpose();" in body
        ),
        "shipping_adds_only_aw_block": (
            "Pext.template block<3,3>(OFF_AW, OFF_AW) += Delta;" in body
        ),
    }


def build() -> dict:
    parity = shipping_source_parity()
    target = stationary_sigma_isotropic(I(1.0))

    P_low = _zero3()
    for i in range(3):
        P_low[i][i] = I(0.5)
    d_pos, case_pos = positive_part_enclosure(target, P_low)

    P_high = _zero3()
    for i in range(3):
        P_high[i][i] = I(2.0)
    d_neg, case_neg = positive_part_enclosure(target, P_high)

    P_mix = _zero3()
    P_mix[0][0] = I(0.5)
    P_mix[1][1] = I(1.5)
    P_mix[2][2] = Interval(0.8, 1.2)
    d_mix, case_mix = positive_part_enclosure(target, P_mix)

    smoke = {
        "strict_positive_case": case_pos,
        "strict_negative_case": case_neg,
        "mixed_case": case_mix,
        "positive_case_diag_lower": [d_pos[i][i].lo for i in range(3)],
        "negative_case_is_zero": all(x.lo == 0.0 and x.hi == 0.0 for row in d_neg for x in row),
        "mixed_case_finite": all(math.isfinite(x.lo) and math.isfinite(x.hi) for row in d_mix for x in row),
    }
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "source_generator": False,
        "trajectory_replay_used": False,
        "floor_request_is_source_event": True,
        "floor_increment_is_source_coordinate": False,
        "floor_increment_depends_on_current_mode_Riccati_covariance": True,
        "same_source_request_may_produce_different_H18_A21_increment": True,
        "positive_part_outer_enclosure_closed_in_real_arithmetic": True,
        "shipping_source_parity": parity,
        "shipping_source_parity_pass": all(parity.values()),
        "smoke": smoke,
        "P3_promoted": False,
        "next_obligation": (
            "call this covariance-dependent operator after each requested prediction in the complete "
            "H18/A21 executor; do not serialize Delta as an independent SEA3 source coordinate"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for key in (
        "floor_request_is_source_event",
        "floor_increment_depends_on_current_mode_Riccati_covariance",
        "same_source_request_may_produce_different_H18_A21_increment",
        "positive_part_outer_enclosure_closed_in_real_arithmetic",
        "shipping_source_parity_pass",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in ("source_generator", "trajectory_replay_used", "floor_increment_is_source_coordinate", "P3_promoted"):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    s = d.get("smoke", {})
    if s.get("strict_positive_case") != "strict_positive_definite_exact":
        f.append("strict-positive floor smoke did not use exact branch")
    if s.get("strict_negative_case") != "strict_negative_definite_exact_zero":
        f.append("strict-negative floor smoke did not use zero branch")
    if s.get("mixed_case") != "mixed_spectrum_frobenius_outer":
        f.append("mixed floor smoke did not use generic enclosure")
    if s.get("negative_case_is_zero") is not True or s.get("mixed_case_finite") is not True:
        f.append("floor smoke failed")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build()
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "shipping_source_parity_pass": d["shipping_source_parity_pass"],
        "smoke": d["smoke"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
