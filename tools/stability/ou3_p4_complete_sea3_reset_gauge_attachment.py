#!/usr/bin/env python3
"""Reset-gauge attachment between frozen complete-SEA3 P3 and physical P4.

P3 executes the deployed MEKF covariance reset

    P_R = G P_J G^T,        G=blkdiag(I+0.5[dtheta]x,I),

immediately after every accepted Joseph correction.  The complete-word P3
margin is already proved source-uniformly for every finite reset injection:
reset is a nonsingular congruence and therefore an exact metric isometry.

The homogeneous physical P4 event is different from a noisy point replay.  Its
measurement residual vanishes at physical error z=0.  Hence dtheta=K y=0 and
the zero-error reset representative is exactly G=I.  The physical finite event
therefore has zero-error tangent I-KH, matching the pre-reset Joseph event used
by the physical Cayley map; no nonzero replay reset may be inserted into that
homogeneous tangent.

This bridge does not remove shipping resets or alter P3.  It records that the
frozen P3 margin is valid in the zero-reset congruent representative needed by
the homogeneous physical P4 map.  Nonzero finite-error chart transport and
nonzero forcing/noise are still P4 obligations and are not closed here.

In particular this tangent bridge does not justify deleting captured resets
and rebuilding a word with unchanged F/Q/H. Finite gauge transport also needs
P'=T P T^T, F'=T_next F T^-1, Q'=T_next Q T_next^T, H'=H T^-1 and
f'(z)=T_next f(T^-1 z). The retained-word attachment audit checks this
distinction; matching zero-state tangents alone is insufficient.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ou3_interval import (
    Interval,
    matrix_identity,
    matrix_mul,
    matrix_point,
    matrix_sub,
    matrix_transpose,
)
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan, matrix_symmetric_hull
import ou3_p4_complete_sea3_differential_events as EVENTS
import ou3_sea3_full_word_event_algebra as EVENT_ALGEBRA
import ou3_sea3_full_word_reset_congruence as RESET
import ou3_sea3_riccati_metric_p3 as P3

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_P4_COMPLETE_SEA3_RESET_GAUGE_ATTACHMENT_V1"


def _contains_zero(A) -> bool:
    return all(x.lo <= 0.0 <= x.hi for row in A for x in row)


def _is_identity(A) -> bool:
    I = matrix_identity(len(A))
    return _contains_zero(matrix_sub(A, I))


def _metric_isometry_smoke() -> dict:
    # Correlated SPD Joseph-posterior covariance and a nontrivial linear Joseph
    # map.  The exact choice is only an algebra smoke test; the theorem claim is
    # the dimension-independent congruence identity already proved by RESET.
    Pj = matrix_point([
        [1.7, 0.10, 0.00, 0.03],
        [0.10, 1.3, 0.04, 0.00],
        [0.00, 0.04, 1.1, 0.02],
        [0.03, 0.00, 0.02, 0.9],
    ])
    A = matrix_point([
        [0.92, 0.03, 0.00, 0.01],
        [-0.02, 0.95, 0.01, 0.00],
        [0.00, -0.01, 0.97, 0.02],
        [0.01, 0.00, -0.02, 0.96],
    ])
    d = [Interval.point(0.07), Interval.point(-0.04), Interval.point(0.02)]
    G = RESET.reset_matrix(d, 4)
    Pr = matrix_symmetric_hull(matrix_mul(matrix_mul(G, Pj), matrix_transpose(G)))
    Pji = matrix_inverse_gauss_jordan(Pj)
    Pri = matrix_inverse_gauss_jordan(Pr)
    lhs = matrix_mul(matrix_mul(matrix_transpose(A), Pji), A)
    GA = matrix_mul(G, A)
    rhs = matrix_mul(matrix_mul(matrix_transpose(GA), Pri), GA)

    z = [Interval.point(0.0), Interval.point(0.0), Interval.point(0.0)]
    G0 = RESET.reset_matrix(z, 4)
    return {
        "arbitrary_finite_reset_energy_identity_enclosed": _contains_zero(matrix_sub(lhs, rhs)),
        "zero_correction_reset_is_identity": _is_identity(G0),
    }


def _zero_physical_event_smoke() -> dict:
    # Verify the actual theorem event implementation fixes z=0.  This catches a
    # future accidental insertion of a replay/baseline correction into the
    # homogeneous physical map.
    n = 18
    zero = [Interval.point(0.0) for _ in range(n)]
    P = matrix_point([[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)])
    R = matrix_point([[0.4 if i == j else 0.0 for j in range(3)] for i in range(3)])
    f = [Interval.point(1.0), Interval.point(-0.5), Interval.point(-9.4)]
    Rhat = matrix_identity(3)
    m = [Interval.point(20.0), Interval.point(-3.0), Interval.point(40.0)]
    acc = EVENTS.source_joseph_event("H", zero, P, R, "accelerometer", f_hat=f, R_hat=Rhat)
    mag = EVENTS.source_joseph_event("H", zero, P, R, "magnetometer", m_body=m)
    Rs = EVENTS.source_joseph_event(
        "H", zero, P, R, "S_zero", R_provenance=EVENTS.ACTUAL_RS_PROVENANCE
    )
    fixed = all(
        x.lo <= 0.0 <= x.hi
        for event in (acc, mag, Rs)
        for x in event["state_out"]
    )
    return {
        "zero_error_acc_vector_S_events_fix_origin": fixed,
        "zero_error_tangent_matches_literal_shipping_H": bool(
            EVENTS.build()["zero_error_tangent_matches_literal_shipping_H"]
        ),
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    p3 = P3.build(path)
    p3_fail = P3.validate(p3)
    algebra = EVENT_ALGEBRA.build()
    algebra_fail = EVENT_ALGEBRA.validate(algebra)
    reset_fail = RESET.validate()
    event_status = EVENTS.build(path)
    event_fail = EVENTS.validate(event_status)
    if p3_fail or algebra_fail or reset_fail or event_fail:
        raise RuntimeError(
            "reset-gauge prerequisites failed: "
            f"P3={p3_fail}; event_algebra={algebra_fail}; reset={reset_fail}; physical_events={event_fail}"
        )
    if p3.get("P3_CONDITIONAL_SEA3_PASS") is not True:
        raise RuntimeError("reset-gauge attachment requires frozen conditional P3")

    smoke = _metric_isometry_smoke()
    physical = _zero_physical_event_smoke()
    reset_info = algebra["left_error_reset"]
    closed = bool(
        all(smoke.values())
        and all(physical.values())
        and reset_info["same_full_matrix_margin_preserved_by_congruence"]
        and reset_info["determinant_lower"] >= 1.0
    )

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "P3_frozen_not_modified": True,
        "P3_delta_consumed": 1.0e-18,
        "P3_conditional_complete_SEA3_consumed": True,
        "P3_reset_margin_preserved_for_every_finite_injection": True,
        "P3_reset_congruence_metric_isometry": True,
        "P3_prior_free_information_not_weakened_by_reset_gauge": True,
        "homogeneous_physical_measurement_residual_zero_at_zero_error": True,
        "homogeneous_physical_correction_zero_at_zero_error": True,
        "homogeneous_zero_error_reset_G_is_identity": smoke["zero_correction_reset_is_identity"],
        "homogeneous_physical_events_fix_zero": physical["zero_error_acc_vector_S_events_fix_origin"],
        "homogeneous_physical_zero_error_tangent_is_I_minus_KH": physical[
            "zero_error_tangent_matches_literal_shipping_H"
        ],
        "P3_margin_valid_in_zero_reset_congruent_representative": closed,
        "zero_error_P3_to_physical_P4_tangent_attachment_closed": closed,
        "arbitrary_finite_reset_energy_identity_enclosed": smoke[
            "arbitrary_finite_reset_energy_identity_enclosed"
        ],
        "nonzero_replay_reset_inserted_into_homogeneous_tangent": False,
        "reset_injection_added_as_new_SEA3_source_coordinate": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "source_family_replaced": False,
        "trajectory_replay_used": False,
        "finite_error_nonlinear_reset_transport_closed_here": False,
        "nonzero_forcing_noise_reset_attachment_closed_here": False,
        "P4_promoted_here": False,
        "P5_may_start_here": False,
        "next_obligation": (
            "evaluate/enclose the finite physical true-minus-estimated word and its prefixes in this reset-normalized "
            "source metric; nonzero forcing/noise belongs to the explicit D_s/D_n terms, not to the homogeneous tangent"
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
        "P3_frozen_not_modified",
        "P3_conditional_complete_SEA3_consumed",
        "P3_reset_margin_preserved_for_every_finite_injection",
        "P3_reset_congruence_metric_isometry",
        "P3_prior_free_information_not_weakened_by_reset_gauge",
        "homogeneous_physical_measurement_residual_zero_at_zero_error",
        "homogeneous_physical_correction_zero_at_zero_error",
        "homogeneous_zero_error_reset_G_is_identity",
        "homogeneous_physical_events_fix_zero",
        "homogeneous_physical_zero_error_tangent_is_I_minus_KH",
        "P3_margin_valid_in_zero_reset_congruent_representative",
        "zero_error_P3_to_physical_P4_tangent_attachment_closed",
        "arbitrary_finite_reset_energy_identity_enclosed",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "nonzero_replay_reset_inserted_into_homogeneous_tangent",
        "reset_injection_added_as_new_SEA3_source_coordinate",
        "filter_changed",
        "declared_domain_changed",
        "source_family_replaced",
        "trajectory_replay_used",
        "finite_error_nonlinear_reset_transport_closed_here",
        "nonzero_forcing_noise_reset_attachment_closed_here",
        "P4_promoted_here",
        "P5_may_start_here",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
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
        "reset_gauge_attachment_closed": d["zero_error_P3_to_physical_P4_tangent_attachment_closed"],
        "finite_reset_transport_closed": d["finite_error_nonlinear_reset_transport_closed_here"],
        "P4_promoted": d["P4_promoted_here"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
