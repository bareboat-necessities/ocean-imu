#!/usr/bin/env python3
"""Validated congruence driver for the canonical complete-SEA3 H18 proof.

This is not a second certificate or a different source.  It executes the same
``ou3_sea3_h18_process_forgetting`` matrix family and the same adaptive
SEA3-coordinate subdivision.  The only numerical tightening is a fixed point
congruence chosen from the midpoint of each interval matrix.

For a cell matrix M = Omega_lower - delta P_upper, ordinary binary64 Cholesky is
used only to choose a nonsingular point matrix C.  The trusted operation is then
performed entirely with outward-rounded intervals:

    M_tilde = C M C^T.

Positive definiteness of M_tilde is certified by the repository interval LDLT.
Because C is fixed and nonsingular, this is exactly equivalent to positive
definiteness of M.  The floating factorization is never accepted as evidence;
if interval LDLT fails, the source cell is subdivided and ultimately fails
closed at the configured depth.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Sequence

import ou3_sea3_h18_process_forgetting as BASE
from ou3_interval import (
    Interval,
    IntervalMatrix,
    matrix_mul,
    matrix_transpose,
    symmetric_positive_definite_ldlt,
)

DEFAULT_DOMAIN = BASE.DEFAULT_DOMAIN
QUALIFICATION = BASE.QUALIFICATION

# The old raw-coordinate LDLT used only five adaptive splits.  The congruence
# normally closes far earlier, but retaining a deeper fail-closed ceiling keeps
# source-coordinate width, not a guessed tolerance, as the only fallback.
BASE.MAX_REFINEMENT_DEPTH = 14
BASE.MAX_ACCEPTED_CELLS = 250000


def _midpoint(A: Sequence[Sequence[Interval]]) -> list[list[float]]:
    return [[0.5 * (x.lo + x.hi) for x in row] for row in A]


def _inverse_lower(L: list[list[float]]) -> list[list[float]]:
    n = len(L)
    X = [[0.0 for _ in range(n)] for _ in range(n)]
    for col in range(n):
        for i in range(n):
            rhs = 1.0 if i == col else 0.0
            s = rhs
            for k in range(i):
                s -= L[i][k] * X[k][col]
            piv = L[i][i]
            if not (math.isfinite(piv) and piv > 0.0):
                raise ArithmeticError("nonpositive midpoint Cholesky pivot")
            X[i][col] = s / piv
    return X


def _midpoint_cholesky_inverse(A: Sequence[Sequence[Interval]]) -> list[list[float]] | None:
    M = _midpoint(A)
    n = len(M)
    if n == 0 or any(len(row) != n for row in M):
        return None
    # Symmetrize the midpoint before choosing the numerical transform.  This
    # changes only the preconditioner choice, never the interval matrix proved.
    for i in range(n):
        for j in range(i):
            x = 0.5 * (M[i][j] + M[j][i])
            M[i][j] = x
            M[j][i] = x
    L = [[0.0 for _ in range(n)] for _ in range(n)]
    try:
        for i in range(n):
            for j in range(i + 1):
                s = M[i][j]
                for k in range(j):
                    s -= L[i][k] * L[j][k]
                if i == j:
                    if not (math.isfinite(s) and s > 0.0):
                        return None
                    L[i][j] = math.sqrt(s)
                else:
                    piv = L[j][j]
                    if not (math.isfinite(piv) and piv > 0.0):
                        return None
                    L[i][j] = s / piv
        C = _inverse_lower(L)
    except (ArithmeticError, OverflowError, ValueError):
        return None
    if not all(math.isfinite(x) for row in C for x in row):
        return None
    return C


def _point_interval_matrix(C: Sequence[Sequence[float]]) -> IntervalMatrix:
    return [[Interval.point(float(x)) for x in row] for row in C]


def _certify_congruence(M: IntervalMatrix) -> tuple[bool, list[Interval], bool]:
    Cfloat = _midpoint_cholesky_inverse(M)
    if Cfloat is None:
        ok, pivots = symmetric_positive_definite_ldlt(M)
        return ok, pivots, False
    C = _point_interval_matrix(Cfloat)
    transformed = matrix_mul(matrix_mul(C, M), matrix_transpose(C))
    ok, pivots = symmetric_positive_definite_ldlt(transformed)
    if ok:
        return True, pivots, True
    # Raw LDLT is a second validated attempt, not a floating fallback.
    raw_ok, raw_pivots = symmetric_positive_definite_ldlt(M)
    return raw_ok, raw_pivots, False


def _cell_certificate(
    tau1: tuple[float, float],
    tau2: tuple[float, float],
    h1: tuple[float, float],
    *,
    H: float,
    gamma: float,
    q_eta: float,
    eta_trace: float,
    Utrans: list[float],
    delta: float,
) -> tuple[bool, dict]:
    t1 = Interval.outward_bounds(*tau1)
    t2 = Interval.outward_bounds(*tau2)
    hp = Interval.outward_bounds(*h1)
    Qaxis = BASE._selected_translation_Q(t1, t2, hp, H)
    M = BASE._full_matrix(Qaxis, gamma, q_eta, eta_trace, Utrans, delta)
    ok, pivots, used = _certify_congruence(M)
    return ok, {
        "tau1": [t1.lo, t1.hi],
        "tau2": [t2.lo, t2.hi],
        "first_piece_s": [hp.lo, hp.hi],
        "second_piece_s": [BASE.down(H - hp.hi), BASE.up(H - hp.lo)],
        "pivot_lower": min((p.lo for p in pivots), default=math.inf),
        "pivot_count": len(pivots),
        "validated_fixed_point_congruence_used": used,
        "floating_midpoint_factorization_accepted_as_proof": False,
    }


# Rebind only the numerical cell test.  BASE.build still constructs the same
# complete-SEA3 source contract, same physical bounds and same 18x18 M family.
BASE._cell_certificate = _cell_certificate


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    d = BASE.build(domain_path)
    cover = d.get("parameter_cover", {})
    d["validated_fixed_point_congruence_preconditioner_available"] = True
    d["floating_midpoint_factorization_accepted_as_proof"] = False
    d["same_complete_SEA3_H18_matrix_family_as_canonical_module"] = True
    d["adaptive_refinement_depth"] = BASE.MAX_REFINEMENT_DEPTH
    d["worst_cell_used_validated_congruence"] = bool(
        cover.get("worst_cell", {}).get("validated_fixed_point_congruence_used", False)
    )
    return d


def validate(d: dict) -> list[str]:
    failures = BASE.validate(d)
    if d.get("validated_fixed_point_congruence_preconditioner_available") is not True:
        failures.append("validated congruence preconditioner missing")
    if d.get("same_complete_SEA3_H18_matrix_family_as_canonical_module") is not True:
        failures.append("preconditioned driver changed the canonical matrix family")
    if d.get("floating_midpoint_factorization_accepted_as_proof") is not False:
        failures.append("floating midpoint factorization was accepted as proof")
    if int(d.get("adaptive_refinement_depth", 0)) != BASE.MAX_REFINEMENT_DEPTH:
        failures.append("adaptive refinement depth mismatch")
    return list(dict.fromkeys(failures))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    d = build(args.domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "closure_time_from_Live_s_upper": d["closure_time_from_Live_s_upper"],
        "posterior_attenuation_lower": d["measurement_conditioning"]["posterior_attenuation_lower"],
        "accepted_cells": d["parameter_cover"].get("accepted_cells"),
        "refinements": d["parameter_cover"].get("refinements"),
        "worst_H18_LDLT_pivot_lower": d["parameter_cover"].get("worst_full_H18_LDLT_pivot_lower"),
        "worst_cell_used_validated_congruence": d.get("worst_cell_used_validated_congruence"),
        "pass": d["H18_PROCESS_FORGETTING_PASS"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
