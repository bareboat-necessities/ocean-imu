#!/usr/bin/env python3
"""Signed complete-SEA3 information ledger for finite-state OU-III P4.

This module composes algebra that was already individually certified but had
not yet been assembled in the form needed by the P4 master inequality.  It does
not create another source, another estimator, a packet-count remainder budget,
or a replacement Lyapunov function.

For one accepted Joseph operation on the SAME shipping P,H,R cell, write the
exact physical residual as

    y = H e + eta,
    S = H P H' + R,
    K = P H' S^-1,
    t = e - K y.

With the Joseph posterior P_J and J=P^-1, J_J=P_J^-1, the exact identity is

    t' J_J t - e' J e
      = - y' S^-1 y + eta' R^-1 eta.                 (1)

Shipping then injects the attitude correction and applies the left covariance
reset P_r=G P_J G'.  If the exact physical post-injection error is

    e_r = G t + rho,
    b = G^-1 rho,

then the exact reset-congruence identity gives

    e_r' P_r^-1 e_r - t' J_J t
      = 2 t' J_J b + b' J_J b.                       (2)

Combining (1)--(2) produces one exact signed Joseph/reset ledger row:

    Delta V
      = - I_y + E_eta + X_reset + E_reset,

    I_y     = y' S^-1 y,
    E_eta   = eta' R^-1 eta,
    X_reset = 2 t' J_J b,
    E_reset = b' J_J b.

This is the useful whole-word organization because:

* an S=0 event has the exact linear residual y=H_S e, hence eta=0;
  every due shipping S update therefore contributes the favorable term
  -y'S^-1y with its ACTUAL applied anisotropic SpectralMSE R_S and no nonlinear
  residual charge;
* accepted accelerometer/vector updates keep their exact nonlinear eta terms;
* the finite Cayley reset mismatch rho stays explicit rather than being erased
  by a tangent-only congruence argument;
* shipping prediction with P-=F P F'+Q and e-=F e cannot increase the moving
  P^-1 energy, and a PSD covariance-floor event with unchanged physical state
  cannot increase it either;
* consequently no N-times standalone remainder estimate is required.  P4 may
  compare the JOINT sum of nonlinear eta/reset costs with the JOINT signed
  information recovered by the same complete SEA3 word.

The exact whole-word endpoint transport remains an independent cross-check of
this same physical map.  P4 is still fail-closed until a source-uniform
same-history bound proves that the joint nonlinear/reset ledger fits inside the
complete-word information decrease for both H18 and A21, and until every
prefix gain/domain condition is certified.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Sequence

import ou3_p4_complete_sea3_differential_events as EVENTS
import ou3_p4_complete_word_endpoint_transport as ENDPOINT
import ou3_p4_exact_reset_transport as RESET
import ou3_sea3_riccati_metric_p3 as P3

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_P4_COMPLETE_SEA3_SIGNED_INFORMATION_LEDGER_V1"


def _shape(A: Sequence[Sequence[Any]]) -> tuple[int, int]:
    rows = len(A)
    cols = len(A[0]) if rows else 0
    if any(len(row) != cols for row in A):
        raise ValueError("ragged matrix")
    return rows, cols


def _mv(A: Sequence[Sequence[Any]], x: Sequence[Any]) -> list[Any]:
    rows, cols = _shape(A)
    if cols != len(x) or cols == 0:
        raise ValueError("matrix/vector dimension mismatch")
    out = []
    for i in range(rows):
        s = A[i][0] * x[0]
        for j in range(1, cols):
            s = s + A[i][j] * x[j]
        out.append(s)
    return out


def _add(a: Sequence[Any], b: Sequence[Any]) -> list[Any]:
    if len(a) != len(b):
        raise ValueError("vector dimension mismatch")
    return [x + y for x, y in zip(a, b)]


def _sub(a: Sequence[Any], b: Sequence[Any]) -> list[Any]:
    if len(a) != len(b):
        raise ValueError("vector dimension mismatch")
    return [x - y for x, y in zip(a, b)]


def _dot(a: Sequence[Any], b: Sequence[Any]) -> Any:
    if len(a) != len(b) or not a:
        raise ValueError("dot-product dimension mismatch")
    s = a[0] * b[0]
    for i in range(1, len(a)):
        s = s + a[i] * b[i]
    return s


def _quad(A: Sequence[Sequence[Any]], x: Sequence[Any]) -> Any:
    return _dot(x, _mv(A, x))


def joseph_signed_energy_terms(
    J_before: Sequence[Sequence[Any]],
    J_joseph: Sequence[Sequence[Any]],
    H: Sequence[Sequence[Any]],
    R_inverse: Sequence[Sequence[Any]],
    S_inverse: Sequence[Sequence[Any]],
    K: Sequence[Sequence[Any]],
    error_before: Sequence[Any],
    nonlinear_eta: Sequence[Any],
) -> dict[str, Any]:
    """Evaluate exact identity (1) in caller-supplied exact/validated arithmetic.

    This is an algebra-regression helper, not a source generator.  The caller
    is responsible for supplying one consistent Joseph tuple P,H,R,S,K and the
    corresponding inverses.  It intentionally does not reconstruct or alter a
    Riccati history.
    """
    n = len(error_before)
    if n not in (18, 21):
        raise ValueError("Joseph signed-energy helper requires H18 or A21")
    if _shape(J_before) != (n, n) or _shape(J_joseph) != (n, n):
        raise ValueError("full-state precision dimension mismatch")
    if _shape(H) != (3, n) or _shape(K) != (n, 3):
        raise ValueError("Joseph H/K dimension mismatch")
    if _shape(R_inverse) != (3, 3) or _shape(S_inverse) != (3, 3):
        raise ValueError("Joseph R/S inverse dimension mismatch")
    if len(nonlinear_eta) != 3:
        raise ValueError("nonlinear residual remainder must have length three")

    linear_residual = _mv(H, error_before)
    residual = _add(linear_residual, nonlinear_eta)
    correction = _mv(K, residual)
    tangent_posterior = _sub(error_before, correction)
    before = _quad(J_before, error_before)
    after = _quad(J_joseph, tangent_posterior)
    information = _quad(S_inverse, residual)
    nonlinear_cost = _quad(R_inverse, nonlinear_eta)
    return {
        "linear_residual": linear_residual,
        "physical_residual": residual,
        "correction": correction,
        "tangent_posterior": tangent_posterior,
        "energy_before": before,
        "energy_after_joseph": after,
        "measurement_information": information,
        "nonlinear_residual_energy": nonlinear_cost,
        "direct_delta": after - before,
        "signed_delta": -information + nonlinear_cost,
        "identity_residual": (after - before) - (-information + nonlinear_cost),
    }


def reset_signed_energy_terms(
    J_joseph: Sequence[Sequence[Any]],
    G_inverse: Sequence[Sequence[Any]],
    tangent_posterior: Sequence[Any],
    reset_defect: Sequence[Any],
) -> dict[str, Any]:
    """Evaluate exact reset-congruence contribution (2)."""
    n = len(tangent_posterior)
    if n not in (18, 21) or len(reset_defect) != n:
        raise ValueError("reset signed-energy helper requires full H18/A21 vectors")
    if _shape(J_joseph) != (n, n) or _shape(G_inverse) != (n, n):
        raise ValueError("reset precision/G inverse dimension mismatch")
    b = _mv(G_inverse, reset_defect)
    reset_cross = 2 * _dot(tangent_posterior, _mv(J_joseph, b))
    reset_energy = _quad(J_joseph, b)
    direct = _quad(J_joseph, _add(tangent_posterior, b)) - _quad(
        J_joseph, tangent_posterior
    )
    return {
        "reset_defect_in_joseph_frame": b,
        "reset_cross_term": reset_cross,
        "reset_defect_energy": reset_energy,
        "direct_delta": direct,
        "signed_delta": reset_cross + reset_energy,
        "identity_residual": direct - (reset_cross + reset_energy),
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    p3 = P3.build(path)
    endpoint = ENDPOINT.build(path)
    reset = RESET.build(path)
    failures = (
        [f"P3: {x}" for x in P3.validate(p3)]
        + [f"endpoint: {x}" for x in ENDPOINT.validate(endpoint)]
        + [f"reset: {x}" for x in RESET.validate(reset)]
    )
    if failures:
        raise RuntimeError(f"signed-information prerequisites failed: {failures}")
    if p3.get("P3_CONDITIONAL_SEA3_PASS") is not True:
        raise RuntimeError("signed-information ledger requires frozen conditional P3")

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "P3_frozen_not_modified": True,
        "P3_delta_consumed": 1.0e-18,
        "same_complete_SEA3_history_required": True,
        "same_shipping_P_H_R_S_K_cell_required": True,
        "joseph_exact_signed_identity": (
            "V_J-V=-y^T*S^-1*y+eta^T*R^-1*eta"
        ),
        "reset_exact_signed_identity": (
            "V_reset-V_J=2*t^T*J_J*(G^-1*rho)+(G^-1*rho)^T*J_J*(G^-1*rho)"
        ),
        "combined_joseph_reset_identity": (
            "DeltaV=-I_y+E_eta+X_reset+E_reset"
        ),
        "S_zero_residual_is_exactly_linear": True,
        "S_zero_nonlinear_eta_exactly_zero": True,
        "every_due_S_update_contributes_negative_information_with_actual_RS": True,
        "actual_RS_provenance_token": EVENTS.ACTUAL_RS_PROVENANCE,
        "actual_RS_retained_in_complete_word_suffixes": bool(
            endpoint["actual_RS_regularization_enters_every_applicable_suffix"]
        ),
        "all_due_S_updates_retained": bool(endpoint["all_due_S_updates_retained"]),
        "all_valid_accelerometer_updates_retained": bool(
            endpoint["all_valid_accelerometer_updates_retained"]
        ),
        "prediction_moving_information_energy_nonincrease_by_Q_PSD": True,
        "aw_covariance_floor_moving_information_energy_nonincrease_by_PSD_order": True,
        "finite_reset_defect_remains_explicit": True,
        "accelerometer_nonlinear_eta_remains_explicit": True,
        "vector_nonlinear_eta_remains_explicit": True,
        "joint_complete_word_signed_information_composition_available": True,
        "packetwise_norm_sum_used": False,
        "packet_count_multiplier_used": False,
        "standalone_eta_Rinv_budget_used": False,
        "global_correction_radius_used": False,
        "inverse_metric_floor_used": False,
        "state_elimination_used": False,
        "trajectory_replay_used": False,
        "source_family_replaced": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "source_uniform_joint_eta_reset_domination_closed": False,
        "source_uniform_endpoint_dissipation_closed_here": False,
        "source_uniform_all_prefix_gain_closed_here": False,
        "P4_promoted_here": False,
        "P5_may_start_here": False,
        "next_obligation": (
            "bound the JOINT sum of accelerometer/vector eta energy and finite reset cross/defect energy over the SAME complete SEA3 word against the JOINT signed information recovered by that word; credit every actual-R_S S event directly, and do not replace the sum by an event-count worst-case bound"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    if float(d.get("P3_delta_consumed", 0.0)) != 1.0e-18:
        f.append("frozen P3 delta changed")
    for key in (
        "P3_frozen_not_modified", "same_complete_SEA3_history_required",
        "same_shipping_P_H_R_S_K_cell_required", "S_zero_residual_is_exactly_linear",
        "S_zero_nonlinear_eta_exactly_zero",
        "every_due_S_update_contributes_negative_information_with_actual_RS",
        "actual_RS_retained_in_complete_word_suffixes", "all_due_S_updates_retained",
        "all_valid_accelerometer_updates_retained",
        "prediction_moving_information_energy_nonincrease_by_Q_PSD",
        "aw_covariance_floor_moving_information_energy_nonincrease_by_PSD_order",
        "finite_reset_defect_remains_explicit", "accelerometer_nonlinear_eta_remains_explicit",
        "vector_nonlinear_eta_remains_explicit",
        "joint_complete_word_signed_information_composition_available",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "packetwise_norm_sum_used", "packet_count_multiplier_used",
        "standalone_eta_Rinv_budget_used", "global_correction_radius_used",
        "inverse_metric_floor_used", "state_elimination_used", "trajectory_replay_used",
        "source_family_replaced", "filter_changed", "declared_domain_changed",
        "source_uniform_joint_eta_reset_domination_closed",
        "source_uniform_endpoint_dissipation_closed_here",
        "source_uniform_all_prefix_gain_closed_here", "P4_promoted_here",
        "P5_may_start_here",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("actual_RS_provenance_token") != EVENTS.ACTUAL_RS_PROVENANCE:
        f.append("actual R_S provenance token changed")
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
        "canonical_source": d["canonical_source"],
        "actual_RS_negative_information": d["every_due_S_update_contributes_negative_information_with_actual_RS"],
        "joint_signed_composition": d["joint_complete_word_signed_information_composition_available"],
        "joint_domination_closed": d["source_uniform_joint_eta_reset_domination_closed"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
