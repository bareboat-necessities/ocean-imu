"""H18 complement/input storage diagnostic, never a finite-error certificate.

For the exact finite-error relation x+ = A x + B e_b + r, r retains the
reference forcing AND the nonlinear/model/arithmetic remainder. A measured
central-difference A is only a local surrogate; it does not bound r.

Given certified A^T W+ A <= alpha W, alpha < rho < 1, Young's inequality
would give V+ <= rho V + rho/(rho-alpha) ||B e_b+r||_W+^2. This module
measures alpha, every-prefix operator gain, and an isolated linear bias gain.
It never asserts the missing finite-error, same-history or domain premises.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from decimal import Decimal, getcontext, localcontext
from pathlib import Path

try:
    from . import h18_superword as base
except ImportError:
    import h18_superword as base

D = Decimal
N = 18
Matrix = list[list[Decimal]]


def transpose(a: Matrix) -> Matrix:
    return [list(row) for row in zip(*a)]


def multiply(a: Matrix, b: Matrix) -> Matrix:
    bt = transpose(b)
    return [[base.dot(row, col) for col in bt] for row in a]


def cholesky(p: Matrix) -> Matrix:
    lower, diagonal = base.ldl_factor(p)
    return [[value * diagonal[j].sqrt() for j, value in enumerate(row)]
            for row in lower]


def forward_solve(lower: Matrix, rhs: Matrix) -> Matrix:
    out = [row[:] for row in rhs]
    for i in range(len(lower)):
        for j in range(len(rhs[0])):
            out[i][j] = (out[i][j] - sum(
                (lower[i][k] * out[k][j] for k in range(i)), D(0))) / lower[i][i]
    return out


def symmetric_max(a: Matrix, *, sweeps: int = 60) -> tuple[Decimal, Decimal]:
    """Cyclic Jacobi maximum and off-diagonal residual (not interval bounds).

    Unlike a stopped power iteration this diagonalizes the whole symmetric
    problem, including clustered dominant eigenvalues and invariant subspaces.
    The residual is numerical only: arithmetic here is not outward rounded.
    """
    n = len(a)
    if n < 1 or any(len(row) != n for row in a):
        raise base.SuperwordExportError("nonempty square matrix required")
    if any(not v.is_finite() for row in a for v in row):
        raise base.SuperwordExportError("nonfinite eigensystem")
    if any(a[i][j] != a[j][i] for i in range(n) for j in range(n)):
        raise base.SuperwordExportError("symmetric eigensystem required")
    work = [row[:] for row in a]
    tolerance = D(10) ** (-getcontext().prec + 12)
    for _ in range(sweeps):
        residual = max(sum((abs(work[i][j]) for j in range(n) if j != i), D(0))
                       for i in range(n))
        largest = max(work[i][i] for i in range(n))
        if residual <= tolerance * max(D(1), abs(largest)):
            return largest, residual
        for p in range(n - 1):
            for q in range(p + 1, n):
                apq = work[p][q]
                if apq == 0:
                    continue
                tau = (work[q][q] - work[p][p]) / (2 * apq)
                t = (D(1) if tau >= 0 else D(-1)) / (abs(tau) + (1 + tau*tau).sqrt())
                c = 1 / (1 + t*t).sqrt()
                s = t*c
                work[p][p] -= t*apq
                work[q][q] += t*apq
                work[p][q] = work[q][p] = D(0)
                for k in range(n):
                    if k in (p, q):
                        continue
                    kp, kq = work[k][p], work[k][q]
                    work[k][p] = work[p][k] = c*kp - s*kq
                    work[k][q] = work[q][k] = s*kp + c*kq
    raise base.SuperwordExportError("Jacobi spectrum did not converge")


def gram(a: Matrix) -> Matrix:
    return multiply(transpose(a), a)


def norm_squared(a: Matrix) -> Decimal:
    value, residual = symmetric_max(gram(a))
    return max(D(0), value + residual)


def young_multiplier(alpha: Decimal, rho: Decimal) -> Decimal:
    """Conditional exact scalar algebra; input alpha still needs certification."""
    if not alpha.is_finite() or not rho.is_finite() or not 0 <= alpha < rho < 1:
        raise base.SuperwordExportError("require finite 0 <= alpha < rho < 1")
    return rho / (rho - alpha)


def linear_bias_gain(m: Matrix, b: Matrix, rho: Decimal) -> Decimal:
    """Schur-complement gain for ||M y+B u||^2 <= rho ||y||^2+gamma ||u||^2.

    This gain covers the isolated linear bias input, not the unknown remainder.
    Q=rho I-M^T M must be positive definite. Gamma is the largest eigenvalue
    of B^T B+(M^T B)^T Q^-1 (M^T B). No explicit inverse is formed.
    """
    if not rho.is_finite() or not 0 < rho < 1:
        raise base.SuperwordExportError("require finite 0 < rho < 1")
    g = gram(m)
    q = [[(rho if i == j else D(0)) - value for j, value in enumerate(row)]
         for i, row in enumerate(g)]
    factor = base.ldl_factor(q)
    cross = multiply(transpose(m), b)
    solved = transpose([base.ldl_solve(factor, col) for col in transpose(cross)])
    correction = multiply(transpose(cross), solved)
    bb = gram(b)
    # LDL solves are rounded; symmetrize only their numerical quadratic form.
    matrix = [[bb[i][j] + (correction[i][j]+correction[j][i])/2
               for j in range(len(bb))] for i in range(len(bb))]
    value, residual = symmetric_max(matrix)
    return max(D(0), value + residual)


def evaluate(export: dict, *, precision: int = 60) -> dict:
    base.validate(export)
    if precision < 32:
        raise base.SuperwordExportError("at least 32 decimal digits required")
    with localcontext() as context:
        context.prec = precision
        cov = [base.symmetric_from_upper(v, 21) for v in export["covariance_upper"]]
        differences = [base.square_from_rows(v, 21) for v in export["error_difference"]]
        root_inverse = base.inverse(differences[0])
        maps = [multiply(v, root_inverse) for v in differences]
        # The split is valid only with a literal held bias and decoupled metric.
        for p, a in zip(cov, maps):
            if any(p[i][j] != 0 for i in range(N) for j in range(N, 21)):
                raise base.SuperwordExportError("held-bias covariance is not decoupled")
            if any(p[i][j] != cov[0][i][j] for i in range(N, 21) for j in range(N, 21)):
                raise base.SuperwordExportError("held-bias covariance changed")
            if any(abs(a[i][j] - (1 if i == j else 0)) > D("1e-25")
                   for i in range(N, 21) for j in range(21)):
                raise base.SuperwordExportError("held-bias identity/zero lower-left block absent")
        root = cholesky([row[:N] for row in cov[0][:N]])
        profile = []
        endpoint = None
        for k, (p, a) in enumerate(zip(cov, maps)):
            end = cholesky([row[:N] for row in p[:N]])
            m = forward_solve(end, multiply([row[:N] for row in a[:N]], root))
            b = forward_solve(end, [row[N:] for row in a[:N]])
            profile.append({"prefix": k, "operator_ratio": float(norm_squared(m)),
                            "bias_input_gain": float(norm_squared(b))})
            endpoint = m, b, end
        m, b, end = endpoint
        alpha = norm_squared(m)
        frozen = forward_solve(root, multiply([row[:N] for row in maps[-1][:N]], root))
        service = json.loads(base.CONSTANTS.read_text(encoding="utf-8"))["magnetic_service"]
        information, _ = base.applied_magnetic_information(
            export, [base.SuperwordMap(v, root_inverse) for v in differences], service)
        coarse = multiply(base.square_from_rows(export["endpoint_difference_coarse"], 21),
                          base.inverse(base.square_from_rows(export["root_difference_coarse"], 21)))
        mc = forward_solve(end, multiply([row[:N] for row in coarse[:N]], root))
        alpha_coarse = norm_squared(mc)
        gap = norm_squared([[x-y for x, y in zip(r, s)] for r, s in zip(m, mc)]).sqrt()
        budget = 1 - alpha.sqrt()
        rho = (1 + alpha)/2 if alpha < 1 else None
        peak = max(profile, key=lambda row: row["operator_ratio"])
        report = {
            "qualification": "OU3_H18_COMPLEMENT_INPUT_DIAGNOSTIC_V1",
            "role": "local central-difference surrogate, not a finite-error inequality",
            "evaluation_precision_digits": precision,
            "complement_dimension": N,
            "endpoint_operator_ratio": float(alpha),
            "coarse_endpoint_operator_ratio": float(alpha_coarse),
            "frozen_root_metric_endpoint_ratio": float(norm_squared(frozen)),
            "magnetic_information_min_eigenvalue": float(information),
            "magnetic_information_meets_floor": information >= base.exact(service["mu_M"]),
            "magnetic_information_phase": "pre_correction",
            "measured_contraction_margin": float(1-alpha),
            "whitened_map_scale_discrepancy": float(gap),
            "admissible_operator_perturbation_budget": float(budget),
            "scale_discrepancy_fits_margin": bool(budget > 0 and gap < budget),
            "discrepancy_is_a_rigorous_error_bound": False,
            "prefix_operator_ratio_max": peak["operator_ratio"],
            "prefix_of_maximum": peak["prefix"],
            "every_prefix_all_directions_evaluated": True,
            "prefix_profile": profile,
            "diagnostic_rho": float(rho) if rho is not None else None,
            "linear_bias_supply_gain": float(linear_bias_gain(m, b, rho)) if rho is not None else None,
            "total_supply_young_multiplier": float(young_multiplier(alpha, rho)) if rho is not None else None,
            "finite_error_remainder_bounded": False,
            "reference_forcing_bounded_in_storage": False,
            "source_uniform": False,
            "certificate_complete": False,
            "obligation_discharged": False,
        }
        return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--export", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--precision", type=int, default=60)
    args = parser.parse_args()
    data = args.export.read_bytes()
    report = evaluate(json.loads(data), precision=args.precision)
    report["input_sha256"] = hashlib.sha256(data).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in report.items() if k != "prefix_profile"}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
