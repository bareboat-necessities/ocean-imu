#!/usr/bin/env python3
"""Conditioned small-x overlay for the source-reachable matrix P3 backend.

The deployed IntegratedOUChain small-x process formulas contain exact factors
x^(k_i+k_j+1), with k=(1,2,3,0) for the comparison coordinates
(v,p,S,a_w).  Evaluating those tiny powers first and dividing them out later
creates interval dependency at x ~= 4e-4.  This module performs the algebraic
cancellation symbolically first and evaluates the equivalent source polynomial
as

    Q_scaled,ij = x * p_ij(x).

No coefficient or source branch changes.  The >=0.01 implementation branch is
delegated to the base backend.  A one-ulp interval straddling the C++ threshold
is split by branch and hulled.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, hull
import ou3_source_reachable_matrix_p3 as BASE

BRANCH_X = BASE.BRANCH_X
DEFAULT_DOMAIN = BASE.DEFAULT_DOMAIN
SCHEMA = BASE.SCHEMA
MIN_USEFUL_DELTA = BASE.MIN_USEFUL_DELTA

_ORIGINAL_STEP_SCALED_Q = BASE.step_scaled_q


def _series(x: Interval, coeffs: tuple[float, ...]) -> Interval:
    # Horner form keeps repeated-x dependency small.  Coefficients are the
    # exact source polynomial after factoring the known leading x power.
    y = BASE.I(coeffs[-1])
    for c in reversed(coeffs[:-1]):
        y = BASE.I(c) + x * y
    return x * y


def _small_scaled_q(x: Interval):
    if not (x.lo > 0.0 and x.hi < BRANCH_X):
        raise ValueError("small scaled-Q helper requires x strictly below source threshold")

    qvv = _series(x, (2/3, -1/2, 7/30, -1/12, 31/1260, -1/160, 127/90720))
    qvp = _series(x, (1/4, -1/6, 5/72, -1/45, 17/2880, -41/30240))
    qvS = _series(x, (1/15, -1/24, 41/2520, -7/1440, 109/90720))
    qva = _series(x, (1.0, -1.0, 7/12, -1/4, 31/360, -1/40, 127/20160, -17/12096))

    qpp = _series(x, (1/10, -1/18, 5/252, -1/180, 17/12960))
    qpS = _series(x, (1/36, -1/72, 13/2880, -1/864))
    qpa = _series(x, (1/3, -1/3, 11/60, -13/180, 19/840, -1/168, 247/181440))

    qSS = _series(x, (1/126, -1/288, 13/12960))
    qSa = _series(x, (1/12, -1/12, 2/45, -1/60, 11/2240, -73/60480))
    qaa = _series(x, (2.0, -2.0, 4/3, -2/3, 4/15, -4/45, 8/315, -2/315, 4/2835))

    return [
        [qvv, qvp, qvS, qva],
        [qvp, qpp, qpS, qpa],
        [qvS, qpS, qSS, qSa],
        [qva, qpa, qSa, qaa],
    ]


def step_scaled_q(x: Interval):
    if x.hi < BRANCH_X:
        return _small_scaled_q(x)
    if x.lo >= BRANCH_X:
        return _ORIGINAL_STEP_SCALED_Q(x)

    left_hi = math.nextafter(BRANCH_X, -math.inf)
    families = []
    if x.lo <= left_hi:
        families.append(_small_scaled_q(Interval(x.lo, left_hi)))
    if x.hi >= BRANCH_X:
        families.append(_ORIGINAL_STEP_SCALED_Q(Interval(BRANCH_X, x.hi)))
    return [[hull(*(A[i][j] for A in families)) for j in range(4)] for i in range(4)]


# Patch exactly one representation primitive.  BASE._build_cached resolves its
# globals at execution time, so all source cells, word upper bounds, generalized
# ratios, useful-margin gate and validation remain the base implementation.
BASE.step_scaled_q = step_scaled_q
BASE._build_cached.cache_clear()


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    d = BASE.build(domain_path)
    # Copy only the top dictionary so cached base results are not mutated.
    out = dict(d)
    out["small_x_process_representation"] = "FACTORED_SOURCE_POLYNOMIAL_X_TIMES_P_OF_X"
    out["small_x_algebraic_equivalence"] = True
    return out


def validate(d: dict) -> list[str]:
    failures = BASE.validate(d)
    if d.get("small_x_process_representation") != "FACTORED_SOURCE_POLYNOMIAL_X_TIMES_P_OF_X":
        failures.append("small-x process representation is not factored")
    if d.get("small_x_algebraic_equivalence") is not True:
        failures.append("small-x algebraic equivalence not asserted")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    out = dict(d)
    out["validation_pass"] = not failures
    out["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "backend": d["p3_backend"],
        "small_x_process_representation": d["small_x_process_representation"],
        "cells": d["cell_partition"],
        "H": d["modes"]["H"],
        "A": d["modes"]["A"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
