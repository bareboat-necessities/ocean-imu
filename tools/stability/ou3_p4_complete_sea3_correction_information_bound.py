#!/usr/bin/env python3
"""Exact per-operation correction inequality; global P4 radii remain OPEN.

For a shipping Kalman operation with P>0, R>0, S=HPH^T+R, K=PH^T S^-1,
let J=y^T S^-1 y and d=Ky. Cauchy--Schwarz in the S metric gives

    d d^T <= J K S K^T = J P H^T S^-1 H P <= J P.

Hence every block selector E satisfies

    ||E d||^2 <= lambda_max(E P E^T) J.

For y=H phi, H^T S^-1 H<=P^-1 also gives J<=phi^T P^-1 phi.
These dimension-independent identities retain all cross terms and apply to
H18 and A21. They do not certify a source-uniform covariance ceiling, a
posterior inverse-metric floor, a finite-angle energy ball, or a complete word.

The previous numerical construction was not a closed reset-radius proof:
(1) prediction-only propagation omitted the immediate nonorthogonal resets;
(2) the H18 ceiling was reused in A21 without a nonlinear hybrid transport;
(3) the cell minimizing a scalar contraction ratio need not minimize the
    posterior process floor; and
(4) an attitude marginal lower is not a lower on its conditional covariance
    (the relevant inverse-metric block is a Schur-complement inverse).

No such numerical bounds are emitted here. In addition, the auxiliary H0 of
the aw-free residual rewrite differs from the congruent shipping H_u=H T_E^T.
A formal Joseph identity for a trial phi storage is not yet a bound on the
physical nonlinear error. Resolve that attachment and the same-history full
SEA3 transport, including every actual-R_S S event, before deriving radii.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import ou3_p4_moving_metric_rebind as REBIND
import ou3_sea3_complete_source as COMPLETE

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 3
QUALIFICATION = "OU3_P4_COMPLETE_SEA3_OPERATION_CORRECTION_IDENTITY_ONLY_V3"


def exact_rational_self_test() -> bool:
    """A non-promoting cross-coupled check of the general algebra above."""
    # P=[[2,1],[1,2]], H=[1,-1], R=1, y=3, S=3, K=[1/3,-1/3].
    # d=[1,-1], J=3. J*P-dd^T=[[5,4],[4,5]] is strictly PSD.
    F = Fraction
    innovation_variance = F(3)
    y = F(3)
    gain = (F(1, 3), F(-1, 3))
    d = tuple(k * y for k in gain)
    j = y * y / innovation_variance
    slack = ((2*j-d[0]*d[0], j-d[0]*d[1]),
             (j-d[1]*d[0], 2*j-d[1]*d[1]))
    return bool(slack[0][0] > 0 and
                slack[0][0]*slack[1][1]-slack[0][1]*slack[1][0] > 0)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    complete = COMPLETE.build(path)
    rebind = REBIND.build()
    failures = COMPLETE.validate(complete) + REBIND.validate(rebind)
    if failures:
        raise RuntimeError(f"correction-identity prerequisites failed: {failures}")
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": complete["canonical_P3_source"],
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "source_family_replaced": False,
        "P3_frozen_not_modified": True,
        "complete_SEA3_word_retained": True,
        "all_valid_accelerometer_updates_remain_in_complete_word": True,
        "all_due_S_updates_and_actual_RS_remain_in_complete_word": True,
        "same_shipping_P_H_R_K_S_cell_required": True,
        "full_matrix_correction_inequality": "d d^T <= J K S K^T <= J P",
        "measurement_information_operator_inequality": "H^T S^-1 H <= P^-1",
        "block_correction_information_inequality": (
            "||E K y||^2 <= lambda_max(E P E^T) y^T S^-1 y"
        ),
        "measurement_energy_inequality_requires_y_equals_H_phi": True,
        "operation_algebra_identity_closed": True,
        "dimension_independent_H18_A21": True,
        "exact_rational_self_test_pass": exact_rational_self_test(),
        "candidate_angles_deg": list(domain["certificate_search"][
            "p4_complete_word_full_attitude_candidate_deg"]),
        "candidate_metric_energy_balls_derived": False,
        "source_uniform_covariance_ceiling_certified_here": False,
        "posterior_inverse_metric_floor_certified_here": False,
        "shipping_Joseph_binding_closed": False,
        "reset_transport_correction_radius_source_closed": False,
        "complete_word_nonlinear_dissipation_closed_here": False,
        "P4_promoted_here": False,
        "old_scalar_Riccati_tube_margin_consumed": False,
        "packet_count_multiplier_used": False,
        "independent_global_correction_radius_assumed": False,
        "open_obligations": [
            "Bind auxiliary residual H0 to actual congruent shipping H_u and gain.",
            "Transport the full source-correlated nonlinear coordinate through prediction, Joseph, reset and H-to-A.",
            "Prove full-matrix/conditional-covariance bounds at each required operation; do not reuse one worst scalar cell.",
            "Dominate the signed complete-word defect with every actual-R_S S update retained.",
        ],
    }


def validate(d: dict) -> list[str]:
    failures = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        failures.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        failures.append("canonical source changed")
    for key in (
        "source_generated_not_trajectory_fit", "P3_frozen_not_modified",
        "complete_SEA3_word_retained", "all_valid_accelerometer_updates_remain_in_complete_word",
        "all_due_S_updates_and_actual_RS_remain_in_complete_word",
        "same_shipping_P_H_R_K_S_cell_required", "operation_algebra_identity_closed",
        "measurement_energy_inequality_requires_y_equals_H_phi",
        "dimension_independent_H18_A21", "exact_rational_self_test_pass",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_changed",
        "source_family_replaced", "candidate_metric_energy_balls_derived",
        "source_uniform_covariance_ceiling_certified_here",
        "posterior_inverse_metric_floor_certified_here", "shipping_Joseph_binding_closed",
        "reset_transport_correction_radius_source_closed", "complete_word_nonlinear_dissipation_closed_here",
        "P4_promoted_here", "old_scalar_Riccati_tube_margin_consumed",
        "packet_count_multiplier_used", "independent_global_correction_radius_assumed",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} must remain false")
    if d.get("candidate_angles_deg") != [30.0, 25.0, 20.0, 15.0]:
        failures.append("candidate cells changed")
    if not d.get("open_obligations"):
        failures.append("open nonlinear transport obligations were hidden")
    return failures


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = build(args.domain)
    failures = validate(result)
    result.update(validation_pass=not failures, validation_failures=failures)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"operation_identity_closed": result["operation_algebra_identity_closed"],
                      "reset_radius_source_closed": result["reset_transport_correction_radius_source_closed"],
                      "validation_failures": failures}, indent=2))
    return 2 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
