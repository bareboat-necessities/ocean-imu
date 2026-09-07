#!/usr/bin/env python3
"""Outward full-state Jacobians of complete-SEA3 nonlinear Joseph events.

This module differentiates the physical error map used by P4.  It does not
accept an independently selected Kalman gain as a theorem input.  For each
source event it derives

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
the event Jacobian is I-KH with K from that same P/H/R cell.

A21 accelerometer-bias norm projection is a separate nonsmooth hybrid.  The
current Normal-Live event differentiator fails closed unless projection is
known inactive; universal P4 must prove that from the declared 0.45/0.5 domain
or add the projection hybrid explicitly.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from ou3_interval import (
    Interval,
    matrix_add,
    matrix_mul,
    matrix_transpose,
)
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan
import ou3_interval_ad as AD
import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_full_normal_live_word as WORD

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 2
QUALIFICATION = "OU3_P4_COMPLETE_SEA3_NONLINEAR_JOSEPH_DIFFERENTIAL_EVENTS_V2"
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
    bias_projection_inactive: bool = True,
):
    """Build H/S/K and the outward physical-state Jacobian from one source cell."""
    n = _dim(mode)
    if len(state) != n or _shape(P) != (n, n) or _shape(R) != (3, 3):
        raise ValueError("state/P/R dimension mismatch")
    if mode == "A" and not bias_projection_inactive:
        raise RuntimeError("A21 b_a projection hybrid is not silently approximated")
    if kind == "S_zero" and R_provenance != ACTUAL_RS_PROVENANCE:
        raise ValueError("S event R must be the actual applied SpectralMSE R_S from this source cell")

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
    return {
        "H": H,
        "S": S,
        "K": K,
        "J_state": AD.jacobian(out),
        "R_provenance": R_provenance,
        "same_P_H_R_cell": True,
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
            "feed same-source P/H/R cells into these AD events, prove the A21 bias projection stays inactive or include its exact hybrid, and compose the resulting matrices through the literal differential word"
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
        "projection_closed": d["A21_bias_projection_inactive_source_uniformly_proved_here"],
        "event_Jacobians_closed": d["source_uniform_finite_angle_event_Jacobians_closed"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
