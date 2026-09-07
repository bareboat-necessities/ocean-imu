#!/usr/bin/env python3
"""Outward full-state Jacobians of complete-SEA3 nonlinear Joseph events.

This module differentiates the physical true-minus-estimated error map used by
P4.  It does not accept an independently selected Kalman gain as a theorem
input.  For each source event it derives

    S = H P H^T + R,
    K = P H^T S^-1

from the SAME source-correlated covariance/geometry/noise cell used by that
shipping operation.  A due S=0 event additionally requires the provenance tag
``ACTUAL_APPLIED_SPECTRALMSE_RS``: target R_S, an independent R_S interval, or
an arbitrary K is not an admissible substitute.

State order is H18=[c,b_g,v,p,S,a_w] or
A21=[c,b_g,v,p,S,a_w,b_a], with exact Cayley attitude c.  An accepted
measurement computes its exact nonlinear residual y(z), d=K y, subtracts the
additive correction, and applies the deployed left quaternion injection.  The
attitude error updates by exact right composition with the inverse correction
quaternion.  First derivatives are enclosed by ``ou3_interval_ad``.

Exact residuals:

* S=0: y=delta S;
* magnetometer/vector: y=(E(c)-I)m_hat;
* accelerometer:
    y=(E(c)-I)f_hat + E(c) R_hat delta_a_w + delta_b_a,
  with b_a absent in H18.

At z=0 these reduce exactly to the literal shipping H_S, H_mag and H_acc, so
the smooth event Jacobian is I-KH with K from that same P/H/R cell.

Shipping projects the A21 residual accelerometer-bias estimate onto the ball of
radius 0.5 m/s^2 after every state injection.  P4 may not assume that projection
is inactive: the declared nominal bound is 0.45 but the finite physical error
cell can cross the boundary.  For the Euclidean ball projection Pi_R(x),

    D Pi_R(x) = I,                                  ||x|| < R,
              = (R/r)(I - x x^T/r^2),              ||x|| > R,

and at ||x||=R the Clarke generalized Jacobian is the convex hull of the two
one-sided limits.  The routines below enclose that generalized Jacobian
outwardly and compose it with the smooth Joseph map.  An A21 theorem event must
therefore supply the same-source true residual-bias cell; omitting it fails
closed rather than silently selecting the inactive branch.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Sequence

from ou3_interval import (
    Interval,
    hull,
    matrix_add,
    matrix_identity,
    matrix_mul,
    matrix_transpose,
)
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan
import ou3_interval_ad as AD
import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_full_normal_live_word as WORD

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 3
QUALIFICATION = "OU3_P4_COMPLETE_SEA3_NONLINEAR_JOSEPH_DIFFERENTIAL_EVENTS_V3"
ACTUAL_RS_PROVENANCE = "ACTUAL_APPLIED_SPECTRALMSE_RS"

OFF_S = 12
OFF_AW = 15
OFF_BA = 18


def _dim(mode: str) -> int:
    if mode == "H":
        return 18
    if mode == "A":
        return 21
    raise ValueError("mode must be H or A")


def _shape(A) -> tuple[int, int]:
    r = len(A)
    c = len(A[0]) if r else 0
    if any(len(row) != c for row in A):
        raise ValueError("ragged matrix")
    return r, c


def _const(x, n: int) -> AD.AD:
    return AD.constant(x, n)


def _const_matrix(A: Sequence[Sequence[Interval]], n: int):
    return [[_const(A[i][j], n) for j in range(len(A[i]))] for i in range(len(A))]


def _state_ad(state: Sequence[Interval]):
    n = len(state)
    if n not in (18, 21) or any(not isinstance(x, Interval) for x in state):
        raise ValueError("state must be an H18/A21 interval vector")
    return AD.independent_vector(state, n=n, offset=0)


def _matvec_const(K: Sequence[Sequence[Interval]], y: Sequence[AD.AD]):
    n = y[0].n
    if _shape(K)[1] != len(y):
        raise ValueError("K/residual dimension mismatch")
    out = []
    for row in K:
        s = _const(0.0, n)
        for k, yy in zip(row, y):
            s = s + _const(k, n) * yy
        out.append(s)
    return out


def source_joseph_gain(P, H, R):
    """Outward K,S from one SAME source P/H/R operation cell."""
    n, n2 = _shape(P)
    if n == 0 or n != n2 or _shape(H) != (3, n) or _shape(R) != (3, 3):
        raise ValueError("source Joseph P/H/R dimensions do not match")
    PCt = matrix_mul(P, matrix_transpose(H))
    S = matrix_add(matrix_mul(H, PCt), R)
    Sinv = matrix_inverse_gauss_jordan(S)
    K = matrix_mul(PCt, Sinv)
    return K, S


def _apply_physical_correction(z, K, residual):
    n = len(z)
    if _shape(K) != (n, 3) or len(residual) != 3:
        raise ValueError("Joseph correction must use full n x 3 K")
    d = _matvec_const(K, residual)
    # Shipping injects +d into the estimate; true-minus-estimated errors
    # therefore subtract d.  Attitude uses exact E+ = E Q(d)^-1.
    cplus = AD.deployed_correct_cayley_right(z[:3], [-d[0], -d[1], -d[2]])
    out = list(cplus)
    out.extend(z[i] - d[i] for i in range(3, n))
    return out


def _sqrt_nonnegative_interval(x: Interval) -> Interval:
    """Outward square-root enclosure for a known nonnegative algebraic quantity.

    ``_norm_interval`` calls this only for a sum of squares.  Generic outward
    interval addition can widen the stored lower endpoint a few ulps below zero
    even though the exact algebraic quantity is nonnegative (notably at the
    zero-state A21 projection parity check).  Clip only that known algebraic
    lower bound.  An interval lying wholly below zero still fails closed.
    """
    if x.hi < 0.0:
        raise ValueError("nonnegative square-root enclosure lies wholly below zero")
    xlo = max(0.0, x.lo)
    lo = 0.0 if xlo == 0.0 else math.nextafter(math.sqrt(xlo), -math.inf)
    hi = math.nextafter(math.sqrt(max(0.0, x.hi)), math.inf)
    return Interval(max(0.0, lo), hi)


def _norm_interval(x: Sequence[Interval]) -> Interval:
    if len(x) != 3 or any(not isinstance(v, Interval) for v in x):
        raise ValueError("projection vector must be three outward intervals")
    r2 = Interval.point(0.0)
    for v in x:
        r2 = r2 + v.square()
    return _sqrt_nonnegative_interval(r2)


def ball_projection_enclosure(x: Sequence[Interval], radius: float) -> dict:
    """Value and Clarke-Jacobian enclosure of Euclidean ball projection.

    ``x`` is the same-source pre-projection estimated residual bias.  The
    returned Jacobian is with respect to x.  If the interval box straddles the
    projection boundary, the result hulls the inactive identity with the active
    radial projection derivative instead of choosing either branch.
    """
    R = float(radius)
    if not (math.isfinite(R) and R > 0.0):
        raise ValueError("projection radius must be finite positive")
    r = _norm_interval(x)
    I3 = matrix_identity(3)
    if r.hi < R:
        return {
            "value": list(x),
            "J": I3,
            "branch": "inactive",
            "norm": r,
        }

    # On the active subset r>=R.  In a boundary-straddling box, clipping the
    # denominator lower endpoint to R is valid because only the active subset
    # is represented by this formula; the identity branch is hulled below.
    r_active = Interval(max(R, r.lo), max(R, r.hi))
    Ri = Interval.point(R)
    alpha = Ri / r_active
    beta = Ri / (r_active * r_active * r_active)
    active_J = []
    for i in range(3):
        row = []
        for j in range(3):
            base = alpha if i == j else Interval.point(0.0)
            row.append(base - beta * x[i] * x[j])
        active_J.append(row)
    active_value = [alpha * v for v in x]

    if r.lo > R:
        return {
            "value": active_value,
            "J": active_J,
            "branch": "active",
            "norm": r,
        }

    return {
        "value": [hull(x[i], active_value[i]) for i in range(3)],
        "J": [[hull(I3[i][j], active_J[i][j]) for j in range(3)] for i in range(3)],
        "branch": "clarke_hull",
        "norm": r,
    }


def _compose_A21_bias_projection(out, bias_true, radius: float):
    """Compose shipping b_a estimate projection with smooth physical error map."""
    if len(out) != 21 or bias_true is None or len(bias_true) != 3:
        raise RuntimeError(
            "A21 projection hybrid requires the same-source true residual-bias cell"
        )
    if any(not isinstance(x, Interval) for x in bias_true):
        raise TypeError("same-source true residual-bias cell must use outward intervals")

    # Smooth correction gives u = b_true - b_hat_corr.  Therefore the estimate
    # presented to project_acc_bias_ is x_hat=b_true-u.  True bias is a frozen
    # source coordinate for this measurement event.
    n = 21
    pre_error = out[OFF_BA:OFF_BA + 3]
    xhat_ad = [_const(bias_true[i], n) - pre_error[i] for i in range(3)]
    xhat = AD.values(xhat_ad)
    proj = ball_projection_enclosure(xhat, radius)

    smooth_J = AD.jacobian(out)
    projected_ba_rows = matrix_mul(proj["J"], smooth_J[OFF_BA:OFF_BA + 3])
    J = [list(row) for row in smooth_J]
    J[OFF_BA:OFF_BA + 3] = projected_ba_rows

    state_out = AD.values(out)
    state_out[OFF_BA:OFF_BA + 3] = [
        bias_true[i] - proj["value"][i] for i in range(3)
    ]
    return J, state_out, proj


def residual_S(z):
    return [z[OFF_S + i] for i in range(3)]


def residual_magnetometer(z, m_body: Sequence[Interval]):
    if len(m_body) != 3:
        raise ValueError("magnetometer source vector must have length three")
    n = len(z)
    E = AD.rotation_from_cayley(z[:3])
    m = [_const(x, n) for x in m_body]
    Em = AD.matvec(E, m)
    return [Em[i] - m[i] for i in range(3)]


def residual_accelerometer(z, f_hat, R_hat):
    if len(f_hat) != 3 or _shape(R_hat) != (3, 3):
        raise ValueError("accelerometer source geometry must be 3-vector plus 3x3 rotation")
    n = len(z)
    E = AD.rotation_from_cayley(z[:3])
    f = [_const(x, n) for x in f_hat]
    R = _const_matrix(R_hat, n)
    da = z[OFF_AW:OFF_AW + 3]
    Rda = AD.matvec(R, da)
    Ef = AD.matvec(E, f)
    ERda = AD.matvec(E, Rda)
    y = [Ef[i] - f[i] + ERda[i] for i in range(3)]
    if n == 21:
        y = [y[i] + z[OFF_BA + i] for i in range(3)]
    return y


def _event_H(mode: str, kind: str, *, f_hat=None, R_hat=None, m_body=None):
    if kind == "S_zero":
        return WORD.H_S_zero(mode)
    if kind in ("magnetometer", "vector"):
        if m_body is None:
            raise ValueError("magnetometer event requires source vector")
        return WORD.H_magnetometer(mode, m_body)
    if kind == "accelerometer":
        if f_hat is None or R_hat is None:
            raise ValueError("accelerometer event requires source force/rotation")
        return WORD.H_accelerometer(mode, f_hat, R_hat)
    raise ValueError("unsupported nonlinear Joseph event kind")


def source_joseph_event(
    mode: str,
    state: Sequence[Interval],
    P,
    R,
    kind: str,
    *,
    f_hat=None,
    R_hat=None,
    m_body=None,
    R_provenance: str | None = None,
    bias_true: Sequence[Interval] | None = None,
    bias_projection_limit: float | None = None,
):
    """Build H/S/K and the outward physical-state Jacobian from one source cell."""
    n = _dim(mode)
    if len(state) != n or _shape(P) != (n, n) or _shape(R) != (3, 3):
        raise ValueError("state/P/R dimension mismatch")
    if kind == "S_zero" and R_provenance != ACTUAL_RS_PROVENANCE:
        raise ValueError("S event R must be the actual applied SpectralMSE R_S from this source cell")
    if mode == "A" and (
        bias_true is None
        or bias_projection_limit is None
        or not (math.isfinite(float(bias_projection_limit)) and float(bias_projection_limit) > 0.0)
    ):
        raise RuntimeError(
            "A21 projection hybrid requires same-source true bias and shipping projection limit"
        )

    H = _event_H(mode, kind, f_hat=f_hat, R_hat=R_hat, m_body=m_body)
    K, S = source_joseph_gain(P, H, R)
    z = _state_ad(state)
    if kind == "S_zero":
        y = residual_S(z)
    elif kind in ("magnetometer", "vector"):
        y = residual_magnetometer(z, m_body)
    else:
        y = residual_accelerometer(z, f_hat, R_hat)
    out = _apply_physical_correction(z, K, y)

    if mode == "A":
        J, state_out, projection = _compose_A21_bias_projection(
            out, bias_true, float(bias_projection_limit)
        )
        projection_branch = projection["branch"]
        projection_J = projection["J"]
        projection_input_norm = projection["norm"]
    else:
        J = AD.jacobian(out)
        state_out = AD.values(out)
        projection_branch = "not_applicable"
        projection_J = None
        projection_input_norm = None

    return {
        "H": H,
        "S": S,
        "K": K,
        "J_state": J,
        "state_out": state_out,
        "R_provenance": R_provenance,
        "same_P_H_R_cell": True,
        "bias_projection_branch": projection_branch,
        "bias_projection_J": projection_J,
        "bias_projection_input_norm": projection_input_norm,
    }


def phi_frame_event_jacobian(Jz, T_in, T_out):
    """Conjugate a physical-state event Jacobian into the full-rank Phi frame."""
    n1, n0 = _shape(Jz)
    if _shape(T_in) != (n0, n0) or _shape(T_out) != (n1, n1):
        raise ValueError("Phi/event differential dimensions do not match")
    Tin_inv = matrix_inverse_gauss_jordan(T_in)
    return matrix_mul(matrix_mul(T_out, Jz), Tin_inv)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    source = COMPLETE.build(path)
    literal = WORD.build(path)
    failures = (
        [f"source: {x}" for x in COMPLETE.validate(source)]
        + [f"literal: {x}" for x in WORD.validate(literal)]
    )
    if failures:
        raise RuntimeError(f"differential Joseph prerequisites failed: {failures}")
    live = json.loads(path.read_text(encoding="utf-8"))["normal_live"]
    projection_limit = float(live["active_accelerometer_bias_projection_limit_mps2"])
    state_cap = float(live["active_accelerometer_bias_state_norm_upper_mps2"])
    projection_margin = projection_limit - state_cap

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "same_complete_SEA3_event_source_required": True,
        "same_P_H_R_cell_derives_S_and_K": True,
        "independent_K_input_allowed_for_theorem": False,
        "actual_applied_RS_provenance_token": ACTUAL_RS_PROVENANCE,
        "actual_applied_RS_required_for_S_event": True,
        "H18_dimension": 18,
        "A21_dimension": 21,
        "exact_Cayley_attitude_state": True,
        "deployed_quaternion_branch_differentiated_outward": True,
        "S_residual_exact": "y=delta_S",
        "magnetometer_residual_exact": "y=(E(c)-I)m_hat",
        "accelerometer_residual_exact": "y=(E(c)-I)f_hat+E(c)R_hat*delta_a_w+delta_b_a",
        "zero_error_tangent_matches_literal_shipping_H": True,
        "physical_correction_sign": "additive_error_plus=additive_error-K*y",
        "physical_attitude_update": "E_plus=E*Q(K_theta*y)^-1",
        "Phi_frame_conjugation_available": True,
        "A21_bias_projection_limit_mps2": projection_limit,
        "A21_bias_state_norm_upper_mps2": state_cap,
        "A21_bias_projection_nominal_margin_mps2": projection_margin,
        "A21_bias_projection_generalized_Jacobian_available": True,
        "A21_bias_projection_same_source_true_bias_required": True,
        "A21_bias_projection_inactive_source_uniformly_proved_here": False,
        "projection_hybrid_silently_ignored": False,
        "packet_count_remainder_budget_used": False,
        "state_elimination_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "trajectory_replay_used": False,
        "source_family_replaced": False,
        "source_uniform_finite_angle_event_Jacobians_closed": False,
        "P4_promoted_here": False,
        "next_obligation": (
            "feed same-source P/H/R/true-bias cells into these AD events, including the exact A21 projection generalized Jacobian, and compose the resulting matrices through the literal differential word"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for key in (
        "same_complete_SEA3_event_source_required", "same_P_H_R_cell_derives_S_and_K",
        "actual_applied_RS_required_for_S_event", "exact_Cayley_attitude_state",
        "deployed_quaternion_branch_differentiated_outward",
        "zero_error_tangent_matches_literal_shipping_H", "Phi_frame_conjugation_available",
        "A21_bias_projection_generalized_Jacobian_available",
        "A21_bias_projection_same_source_true_bias_required",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "independent_K_input_allowed_for_theorem", "A21_bias_projection_inactive_source_uniformly_proved_here",
        "projection_hybrid_silently_ignored", "packet_count_remainder_budget_used", "state_elimination_used",
        "filter_changed", "declared_domain_changed", "trajectory_replay_used", "source_family_replaced",
        "source_uniform_finite_angle_event_Jacobians_closed", "P4_promoted_here",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("actual_applied_RS_provenance_token") != ACTUAL_RS_PROVENANCE:
        f.append("actual R_S provenance token changed")
    if d.get("H18_dimension") != 18 or d.get("A21_dimension") != 21:
        f.append("full-state dimensions changed")
    if not float(d.get("A21_bias_projection_limit_mps2", 0.0)) > float(
        d.get("A21_bias_state_norm_upper_mps2", 1.0)
    ):
        f.append("declared A21 bias projection has no nominal interior margin")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "source": d["canonical_source"],
        "same_cell_gain": d["same_P_H_R_cell_derives_S_and_K"],
        "actual_RS_S": d["actual_applied_RS_required_for_S_event"],
        "projection_generalized_J": d["A21_bias_projection_generalized_Jacobian_available"],
        "event_Jacobians_closed": d["source_uniform_finite_angle_event_Jacobians_closed"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())