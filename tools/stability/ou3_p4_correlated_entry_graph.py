#!/usr/bin/env python3
"""Exact algebra for the fresh-wrapper OU-III linear-state entry graph.

This is an entry/source identity, NOT a P4 contraction certificate. Physical
BRMM truth and its admissible continuation are not replaced by free roots.
The shipping gain must still come from the same reachable P/H/R. Identities
valid for every gain do not license independent gain boxes in a stability
proof. All arithmetic below is exact rational arithmetic.

Let e = x_true - x_hat for x=(v,p,S,a_w), in the unchanged physical gauge.
At fresh-wrapper first Live entry x_hat=0, hence [e;x_true]=G x_true,
G=[I;I]. In particular r_S=-S_hat=[H,-H][e;x_true]=0 at that instant.
This does NOT say S_true=0, and does not say the first later S event has
zero residual: intervening accelerometer/magnetometer corrections matter.

For prediction, define d=x_true_next-F*x_true from the SAME physical
history, not as an independent selectable forcing. Then the exact joint
map is diag(F,F)[e;x_true]+G*d. A linear-state measurement correction with
residual H*(e-x_true) has map [[I-KH,KH],[0,I]]. These maps retain the
source/error correlation algebraically. Quaternion reset, bias projection,
other sensor residuals and binary32 errors remain separate obligations.
"""
from __future__ import annotations

from fractions import Fraction
import argparse
import json
from pathlib import Path
from typing import Iterable, Sequence

Q = Fraction
Matrix = tuple[tuple[Fraction, ...], ...]


def matrix(rows: Iterable[Iterable[int | Fraction]]) -> Matrix:
    result = []
    for row in rows:
        converted = []
        for value in row:
            if isinstance(value, bool) or not isinstance(value, (int, Fraction)):
                raise TypeError("certificate matrices require int or Fraction; no floats")
            converted.append(Fraction(value))
        result.append(tuple(converted))
    if not result or not result[0] or any(len(r) != len(result[0]) for r in result):
        raise ValueError("matrix must be nonempty and rectangular")
    return tuple(result)


def identity(n: int) -> Matrix:
    if isinstance(n, bool) or not isinstance(n, int) or n <= 0:
        raise ValueError("positive integer dimension required")
    return matrix([[int(i == j) for j in range(n)] for i in range(n)])


def transpose(a: Matrix) -> Matrix:
    return tuple(zip(*a))


def multiply(a: Matrix, b: Matrix) -> Matrix:
    if len(a[0]) != len(b):
        raise ValueError("incompatible matrix dimensions")
    out = [[Q(0) for _ in b[0]] for _ in a]
    for i, row in enumerate(a):
        for k, value in enumerate(row):
            if value:
                for j, other in enumerate(b[k]):
                    if other:
                        out[i][j] += value * other
    return matrix(out)


def fresh_entry_lift(n: int = 12) -> Matrix:
    """Map the SAME physical root to error and truth; not two independent balls."""
    unit = identity(n)
    return unit + unit


def innovation_rows(n: int = 12, indices: Sequence[int] = (6, 7, 8)) -> Matrix:
    """Rows selecting e_S-S_true in the [linear error; linear truth] order."""
    identity(n)  # validate dimension
    if not indices or len(set(indices)) != len(indices):
        raise ValueError("nonempty distinct innovation indices required")
    if any(isinstance(i, bool) or not isinstance(i, int) or not 0 <= i < n for i in indices):
        raise ValueError("innovation index out of range")
    return matrix([[int(j == i) - int(j == n + i) for j in range(2*n)] for i in indices])


def finite_linear_measurement_map(k: Matrix, indices: Sequence[int]) -> Matrix:
    """Exact additive LINEAR-state rows of a finite S=0 correction.

    This does not linearize the attitude injection or radial projection.
    k is the linear-state gain row block, not an independently admitted gain.
    """
    n = len(k)
    c = innovation_rows(n, indices)
    if len(k[0]) != len(indices):
        raise ValueError("gain/innovation dimensions differ")
    kc = multiply(k, c)
    out = [list(row) for row in identity(2*n)]
    for i in range(n):
        for j in range(2*n):
            out[i][j] -= kc[i][j]
    return matrix(out)


def joint_prediction_map(f: Matrix) -> Matrix:
    n = len(f)
    if len(f[0]) != n:
        raise ValueError("prediction map must be square")
    out = [[Q(0) for _ in range(2*n)] for _ in range(2*n)]
    for i in range(n):
        for j in range(n):
            out[i][j] = out[n+i][n+j] = f[i][j]
    return matrix(out)


def restrict_joint_quadratic(storage: Matrix, lift: Matrix) -> Matrix:
    """Exact pullback G^T M G; does not assert coercivity or contraction.

    The caller must retain the physical admission constraints on the root of
    G. This restriction is valid only where the entry graph is established;
    it must not be imposed again at later Live words or at A21 release.
    """
    if len(storage) != len(storage[0]) or len(storage) != len(lift):
        raise ValueError("storage/lift dimensions differ")
    if storage != transpose(storage):
        raise ValueError("storage must be symmetric")
    return multiply(transpose(lift), multiply(storage, lift))


def certificate() -> dict:
    """Check exact coefficients of the gain-affine entry identity.

    Checking zero and every gain basis coefficient is exhaustive polynomial
    coefficient verification, not sampled physical histories. It establishes
    D(K)G=G for ALL K because D(K) is affine in K, and establishes no bounds on
    K, P, physical roots, storage, or nonlinear trajectories.
    """
    n = 12
    selected = (6, 7, 8)
    g = fresh_entry_lift(n)
    c = innovation_rows(n, selected)
    zero = matrix([[0]*n for _ in selected])
    checks = {"entry_residual_coefficients_zero": multiply(c, g) == zero}
    zero_gain = matrix([[0]*3 for _ in range(n)])
    checks["measurement_constant_coefficient"] = multiply(
        finite_linear_measurement_map(zero_gain, selected), g) == g
    verified = 0
    for i in range(n):
        for j in range(3):
            basis = [[0]*3 for _ in range(n)]
            basis[i][j] = 1
            if multiply(finite_linear_measurement_map(matrix(basis), selected), g) != g:
                raise AssertionError("nonzero measurement polynomial coefficient")
            verified += 1
    checks["all_gain_coefficients_zero"] = verified == n*3
    # The prediction identity diag(F,F)G=GF is linear in each entry of F.
    prediction_coefficients = 0
    for i in range(n):
        for j in range(n):
            basis = [[0]*n for _ in range(n)]
            basis[i][j] = 1
            f = matrix(basis)
            if multiply(joint_prediction_map(f), g) != multiply(g, f):
                raise AssertionError("nonzero prediction polynomial coefficient")
            prediction_coefficients += 1
    checks["all_prediction_coefficients_zero"] = prediction_coefficients == n*n
    checks["same_driver_innovation_coefficients_zero"] = multiply(c, g) == zero
    return {
        "claim": "exact fresh-entry physical-source/error graph identities only",
        "linear_state_order": ["v", "p", "S", "a_w"],
        "entry_graph": "e_linear=x_true_linear; x_hat_linear=0",
        "physical_truth_gauge": "unchanged; no assertion p_true=S_true=0",
        "gain_polynomial_coefficients_checked": verified,
        "prediction_polynomial_coefficients_checked": prediction_coefficients,
        "checks": checks,
        "identity_certificate_pass": all(checks.values()),
        "scope_exclusions": [
            "generic goLive after a previously driven kernel",
            "subsequent Live word entrances and H18-to-A21 release",
            "residuals after intervening sensor corrections",
            "physical BRMM root/continuation qualification",
            "uniform reachable Riccati and adaptation cover",
            "nonlinear endpoint/prefix contraction and retention",
            "finite precision enclosure and finite-time startup capture"
        ]
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = certificate()
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    if not result["identity_certificate_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
