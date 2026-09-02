#!/usr/bin/env python3
"""Cell-affine Joseph dissipation bound for finite OU-III P4 subdivision cells.

A complete-word P4 proof cannot treat every interval subdivision cell as if it
were centred at the zero-error state. For a convex entry cell X with centre
x_c, an exact nonlinear residual r(x) and a certified interval Jacobian J of
the *actual deterministic residual map from the chosen word-entry state*, the
componentwise mean-value theorem gives

    r(x) = A x + e,
    A in J,
    e = r(x_c) - A x_c.

For any epsilon in (0,1), Young's inequality in the positive innovation metric
gives

    r^T S^-1 r
      >= (1-epsilon) x^T A^T S^-1 A x
         - (1-epsilon)/epsilon * e^T S^-1 e.

Together with the exact finite-angle Joseph identity

    Delta W = r^T S^-1 r - eta^T R^-1 eta,

this yields one finite-cell inequality

    Delta W >= mu W_entry - beta.

The important promotion guard is explicit. In OU-III, accepted attitude
corrections change covariance through the reset congruence. Later P,H,S^-1 and
the correction therefore depend on the entry error. A Jacobian obtained by
propagating only the error AD variables while treating interval covariance as a
frozen parameter is useful as a design enclosure, but is not by itself the
Jacobian of that deterministic complete-word map. Callers must assert that all
state-dependent covariance/gain/reset paths have either been included in the
certified Jacobian or rigorously moved into an affine remainder. The helper
fails closed otherwise.

No K interval matrix is formed, no per-correction sector-contraction premise is
introduced, and P3's tiny positive information margin is not reinterpreted as a
physical basin.
"""
from __future__ import annotations

import argparse
import json
import math
from typing import Sequence

from ou3_interval import Interval
import ou3_interval_ad as AD
import ou3_p4_joint_word_dissipation_design as D

SCHEMA = 2


def I(x: float) -> Interval:
    return Interval.point(float(x))


def _mat_vec(A, x: Sequence[Interval]) -> list[Interval]:
    if not A or len(A[0]) != len(x):
        raise ValueError("matrix/vector dimension mismatch")
    out = []
    for row in A:
        z = I(0.0)
        for a, b in zip(row, x):
            z = z + a * b
        out.append(z)
    return out


def _quad_upper(x: Sequence[Interval], M) -> float:
    y = _mat_vec(M, x)
    z = I(0.0)
    for a, b in zip(x, y):
        z = z + a * b
    if not math.isfinite(z.hi):
        raise RuntimeError("quadratic upper bound is not finite")
    if z.hi <= 0.0:
        return 0.0
    return math.nextafter(float(z.hi), math.inf)


def _eta_Rinv_upper(eta: Sequence[Interval], R) -> float:
    if len(R) != len(eta) or any(len(row) != len(eta) for row in R):
        raise ValueError("eta/R dimension mismatch")
    total = 0.0
    for i, e in enumerate(eta):
        for j in range(len(eta)):
            if i != j and not (R[i][j].lo == 0.0 and R[i][j].hi == 0.0):
                raise RuntimeError("cell-affine eta penalty currently requires diagonal R")
        if not R[i][i].lo > 0.0:
            raise RuntimeError("measurement R lost a positive diagonal floor")
        a = e.abs_upper()
        if a == 0.0:
            continue
        term = math.nextafter(a * a / R[i][i].lo, math.inf)
        total = math.nextafter(total + term, math.inf)
    return total


def _scale_form(A, scalar: float):
    q = Interval.outward_bounds(float(scalar), float(scalar))
    return [[q * x for x in row] for row in A]


def _affine_intercept(center_residual: Sequence[Interval], J,
                      center_state: Sequence[float]) -> list[Interval]:
    if len(center_residual) != len(J):
        raise ValueError("center residual/Jacobian row mismatch")
    if J and len(J[0]) != len(center_state):
        raise ValueError("center state/Jacobian column mismatch")
    xc = [I(float(x)) for x in center_state]
    Jxc = _mat_vec(J, xc)
    return [a - b for a, b in zip(center_residual, Jxc)]


def bound_cell(residual, eta, *, center_residual: Sequence[Interval],
               center_state: Sequence[float], Sinv, R, entry_metric_upper,
               complete_word_jacobian_covariance_paths_certified: bool,
               epsilon: float = 0.5) -> dict:
    """Return a finite-cell ``Delta W >= mu W - beta`` certificate.

    The boolean covariance-path premise is intentionally mandatory. It may be
    true only if the supplied residual Jacobian encloses every dependence of the
    actual word map on its entry state, including state-dependent covariance
    reset/gain paths, or if those omitted paths have been covered by a separate
    rigorous remainder already included in ``residual``/``eta``. A caller may
    not silently promote the current frozen-covariance design AD by setting the
    flag true.
    """
    if complete_word_jacobian_covariance_paths_certified is not True:
        raise RuntimeError(
            "finite-cell promotion requires a Jacobian/remainder that includes "
            "state-dependent covariance, gain and reset paths"
        )
    eps = float(epsilon)
    if not (0.0 < eps < 1.0):
        raise ValueError("epsilon must lie strictly in (0,1)")
    J = AD.jacobian(residual)
    if not J:
        raise ValueError("empty residual Jacobian")
    if len(center_state) != len(J[0]):
        raise ValueError("entry-state dimension mismatch")
    if len(entry_metric_upper) != len(center_state):
        raise ValueError("entry metric dimension mismatch")

    raw_form = D._form(J, Sinv)
    retained = math.nextafter(1.0 - eps, -math.inf)
    directional_form = _scale_form(raw_form, retained)
    mu = float(D._generalized_margin(directional_form, entry_metric_upper))

    ebox = _affine_intercept(center_residual, J, center_state)
    eSinve = _quad_upper(ebox, Sinv)
    young_coeff = math.nextafter(retained / eps, math.inf)
    beta_affine = 0.0 if eSinve == 0.0 else math.nextafter(young_coeff * eSinve, math.inf)

    eta_box = AD.values(eta)
    beta_eta = _eta_Rinv_upper(eta_box, R)
    beta = 0.0 if beta_affine == 0.0 and beta_eta == 0.0 else math.nextafter(beta_affine + beta_eta, math.inf)
    rho = math.nextafter(1.0 - mu, math.inf) if mu > 0.0 else 1.0

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_FINITE_CELL_AFFINE_JOSEPH_DISSIPATION",
        "identity": "DeltaW=r^T S^-1 r-eta^T R^-1 eta >= mu W_entry-beta",
        "mean_value_form": "r(x)=A x+e, A in certified full-map Jacobian, e=r(x_c)-A x_c",
        "complete_word_jacobian_covariance_paths_certified": True,
        "frozen_covariance_error_only_AD_sufficient_for_promotion": False,
        "young_epsilon": eps,
        "young_retained_information_factor_lower": retained,
        "directional_generalized_margin_lower": mu,
        "rho_homogeneous_upper": rho,
        "affine_intercept_box": [[x.lo, x.hi] for x in ebox],
        "affine_intercept_information_upper": eSinve,
        "beta_affine_intercept_upper": beta_affine,
        "beta_nonlinear_eta_upper": beta_eta,
        "beta_total_upper": beta,
        "K_interval_matrix_materialized": False,
        "per_correction_sector_contraction_required": False,
        "P3_relative_margin_used_as_physical_basin": False,
        "finite_cell_not_forced_through_origin": True,
        "P4_PROMOTED_HERE": False,
        "P5_PROMOTED_HERE": False,
    }


def _self_test() -> dict:
    X = Interval.outward_bounds(1.0, 2.0)
    x = AD.independent(X, 0, 1)
    residual = [x + AD.constant(0.1, 1)]
    eta = [AD.constant(0.0, 1)]
    failures = []

    # The production guard must reject an unqualified Jacobian.
    rejected = False
    try:
        bound_cell(
            residual, eta,
            center_residual=[I(1.6)], center_state=[1.5],
            Sinv=[[I(1.0)]], R=[[I(1.0)]], entry_metric_upper=[[I(1.0)]],
            complete_word_jacobian_covariance_paths_certified=False,
            epsilon=0.5,
        )
    except RuntimeError:
        rejected = True
    if not rejected:
        failures.append("incomplete covariance-path Jacobian was not rejected")

    out = bound_cell(
        residual,
        eta,
        center_residual=[I(1.6)],
        center_state=[1.5],
        Sinv=[[I(1.0)]],
        R=[[I(1.0)]],
        entry_metric_upper=[[I(1.0)]],
        complete_word_jacobian_covariance_paths_certified=True,
        epsilon=0.5,
    )
    mu = float(out["directional_generalized_margin_lower"])
    beta = float(out["beta_total_upper"])
    if not (0.0 < mu < 0.50000000000001):
        failures.append("scalar generalized margin is not the expected strict half-information bound")
    if not (0.0 <= beta < 0.011):
        failures.append("scalar affine beta is outside expected range")
    for xv in (1.0, 1.5, 2.0):
        lhs = (xv + 0.1) ** 2
        rhs = mu * xv * xv - beta
        if lhs < rhs:
            failures.append(f"scalar inequality failed at x={xv}")
    if out.get("complete_word_jacobian_covariance_paths_certified") is not True:
        failures.append("full-map Jacobian premise missing")
    if out.get("frozen_covariance_error_only_AD_sufficient_for_promotion") is not False:
        failures.append("frozen-covariance AD was promoted")
    for key in (
        "K_interval_matrix_materialized",
        "per_correction_sector_contraction_required",
        "P3_relative_margin_used_as_physical_basin",
        "P4_PROMOTED_HERE",
        "P5_PROMOTED_HERE",
    ):
        if out.get(key) is not False:
            failures.append(f"{key} is not false")
    if out.get("finite_cell_not_forced_through_origin") is not True:
        failures.append("finite-cell affine flag missing")
    return {"pass": not failures, "failures": failures, "guard_rejected_unqualified": rejected, "bound": out}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if not args.self_test:
        ap.error("this proof primitive currently exposes only --self-test as a standalone command")
    d = _self_test()
    print(json.dumps(d, indent=2, sort_keys=True))
    return 0 if d["pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
