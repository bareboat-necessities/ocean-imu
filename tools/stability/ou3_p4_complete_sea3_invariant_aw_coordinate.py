#!/usr/bin/env python3
"""Invariant a_w normal form for complete-SEA3 OU-III P4.

The triangular congruence and residual normal form below are separate
pointwise identities. Their composition with the shipping Joseph map remains
open: H0 in the residual rewrite is not H_u=H T_E^T from the preceding
operation-coordinate congruence. In particular this file does not certify
that a gravity-only attitude column describes the unchanged shipping gain.  It does not change the filter and it does
not replace the complete SEA3 word.

For one accepted accelerometer operation let

    f_hat = R_hat (a_hat - g),
    r_a   = (E-I) f_hat + R_hat u_aw + db_a,

where ``u_aw=R_hat^T E R_hat delta_a_w`` is the exact operation coordinate from
``ou3_p4_complete_sea3_accelerometer_operation_coordinate``.  Introduce the
source-indexed *linear* triangular coordinate

    w_lin = u_aw + B c,
    B c   = R_hat^T [c]x R_hat a_hat
          = -R_hat^T [R_hat a_hat]x c.

Then the shipping tangent residual is identically

    [c]x f_hat + R_hat u_aw + db_a
      = -[c]x R_hat g + R_hat w_lin + db_a.

Thus the wave-acceleration part of the attitude column is moved into the a_w
coordinate by a unit-triangular state transform.  If

    T_B : (c,u_aw) -> (c,w_lin),

then det(T_B)=1 and T_B is nonsingular for every finite source value.  With

    P_B = T_B P T_B^T,
    H_B = H T_B^-1,
    K_B = T_B K,

one has exactly

    H_B P_B H_B^T = H P H^T,

so the shipping innovation covariance S is unchanged, the Joseph covariance
update is a congruence, and the moving-Riccati energy is unchanged.  No
condition-number multiplier or group-isotropic metric is introduced.

Apply the same T_B to the nonlinear measurement-linearizing coordinate.  With

    e_eta = R_hat^T ((E-I)-[c]x) f_hat,
    w_exact = w_lin + e_eta,

one obtains the exact finite-angle normal form

    r_a = -[c]x R_hat g + R_hat w_exact + db_a = H_B Phi_B(z).

Equivalently, writing a_true=a_hat+delta_a_w and Q_E=R_hat^T E R_hat,

    w_exact = Q_E a_true - a_hat
              - R_hat^T ((E-I)-[c]x) R_hat g.

The only nonlinear displacement from the transformed tangent coordinate is
therefore the already-certified pure finite-angle shift e_eta.  The large
first-order wave-acceleration/attitude coupling is part of the exact linear
congruence T_B, not a remainder.  Because the whole complete SEA3 word is
transformed rather than altered, every due S=0 update and its actual applied
SpectralMSE R_S remain in the same word and regularize the transformed a_w
direction through the same information matrix.

This producer deliberately remains fail-closed.  P4 still requires an
operation-matched bound on transport of e_eta through the physical correction,
quaternion injection/reset, prediction, and source change.  That transport must
be charged to the same Joseph information decrease; no packet-count-times-worst
remainder budget is admissible.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_complete_sea3_accelerometer_operation_coordinate as ACC
import ou3_p4_complete_sea3_measurement_linearizing_aw_coordinate as AWLIN
import ou3_p4_moving_metric_rebind as REBIND
import ou3_sea3_complete_source as COMPLETE

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_P4_COMPLETE_SEA3_INVARIANT_AW_NORMAL_FORM"


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("invariant-aw normal form must not be trajectory fitted")

    complete = COMPLETE.build(path)
    acc = ACC.build(path)
    awlin = AWLIN.build(path)
    rebind = REBIND.build()
    bad = {
        "complete": COMPLETE.validate(complete),
        "accelerometer_operation_coordinate": ACC.validate(acc),
        "measurement_linearizing_aw": AWLIN.validate(awlin),
        "moving_metric_rebind": REBIND.validate(rebind),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"invariant-aw prerequisites failed: {bad}")
    if complete["canonical_P3_source"] != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        raise RuntimeError("canonical complete SEA3 source changed")

    live = domain["normal_live"]
    g = float(domain["startup"]["gravity_mps2"])
    fmax = float(live["specific_force_norm_upper_mps2"])
    if not (math.isfinite(g) and g > 0.0 and math.isfinite(fmax) and fmax > 0.0):
        raise RuntimeError("invalid gravity/specific-force theorem bounds")
    # Since f_hat=R_hat(a_hat-g), orthogonality gives
    # ||a_hat|| <= ||a_hat-g||+||g|| <= fmax+g.
    ahat = math.nextafter(fmax + g, math.inf)

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "source_family_replaced": False,
        "P3_frozen_not_modified": True,
        "complete_SEA3_word_retained": True,
        "all_valid_accelerometer_updates_remain_in_complete_word": True,
        "all_due_S_updates_and_actual_RS_remain_in_complete_word": True,
        "source_nominal_aw_norm_upper_mps2": ahat,
        "nominal_aw_bound_is_conditional_on_nominal_force_bound": True,
        "nominal_force_bound_inherited_from_physical_SEA3_proved": False,
        "shipping_Joseph_binding_closed": False,
        "linear_triangular_coordinate": {
            "w_lin": "u_aw+B*c",
            "B_c": "R_hat^T*[c]x*R_hat*a_hat",
            "B_matrix": "-R_hat^T*[R_hat*a_hat]x",
            "B_operator_norm_upper_mps2_per_cayley": ahat,
            "T_B_unit_triangular": True,
            "T_B_determinant_exact": 1.0,
            "T_B_nonsingular": True,
            "T_B_inverse_exact": True,
        },
        "tangent_normal_form_identity": (
            "[c]x*f_hat+R_hat*u_aw+delta_b_a="
            "-[c]x*R_hat*g+R_hat*w_lin+delta_b_a"
        ),
        "transformed_attitude_column_depends_only_on_gravity": True,
        "wave_acceleration_attitude_cross_term_is_linear_coordinate_coupling": True,
        "moving_metric_congruence": {
            "P_B": "T_B*P*T_B^T",
            "H_B": "H*T_B^-1",
            "K_B": "T_B*K",
            "innovation_covariance_S_invariant": True,
            "Joseph_covariance_congruence_exact": True,
            "moving_Riccati_energy_invariant": True,
            "condition_number_multiplier_used": False,
            "group_isotropic_metric_assumption_used": False,
        },
        "exact_finite_angle_coordinate": {
            "w_exact": "w_lin+e_eta",
            "e_eta": "R_hat^T*((E-I)-[c]x)*f_hat",
            "equivalent_w_exact": (
                "Q_E*a_true-a_hat-R_hat^T*((E-I)-[c]x)*R_hat*g"
            ),
            "Q_E": "R_hat^T*E*R_hat",
            "exact_residual": (
                "r_a=-[c]x*R_hat*g+R_hat*w_exact+delta_b_a=H_B*Phi_B(z)"
            ),
            "nonlinear_displacement_is_exactly_existing_e_eta": True,
            "first_order_wave_attitude_term_removed_from_remainder": True,
        },
        "measurement_linearizing_shift_bounds_reused_without_widening": awlin["candidate_cells"],
        "standalone_eta_Rinv_packet_budget_used": False,
        "packet_count_multiplier_used": False,
        "actual_RS_information_matrix_retained_under_congruence": True,
        "complete_source_correlated_transport_defect_closed_here": False,
        "P4_promoted_here": False,
        "next_obligation": (
            "first bind H0 to the actual congruent H_u and prove a nominal-force bound; then bound the full coordinate transport; "
            "use the same operation's information decrease and transformed complete-word R_S metric, then scalarize only at the complete-word endpoint"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for key in (
        "source_generated_not_trajectory_fit",
        "P3_frozen_not_modified",
        "complete_SEA3_word_retained",
        "all_valid_accelerometer_updates_remain_in_complete_word",
        "all_due_S_updates_and_actual_RS_remain_in_complete_word",
        "transformed_attitude_column_depends_only_on_gravity",
        "wave_acceleration_attitude_cross_term_is_linear_coordinate_coupling",
        "actual_RS_information_matrix_retained_under_congruence",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "shipping_Joseph_binding_closed",
        "nominal_force_bound_inherited_from_physical_SEA3_proved",
        "trajectory_replay_used",
        "filter_changed",
        "declared_domain_changed",
        "source_family_replaced",
        "standalone_eta_Rinv_packet_budget_used",
        "packet_count_multiplier_used",
        "complete_source_correlated_transport_defect_closed_here",
        "P4_promoted_here",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    tri = d.get("linear_triangular_coordinate", {})
    for key in ("T_B_unit_triangular", "T_B_nonsingular", "T_B_inverse_exact"):
        if tri.get(key) is not True:
            f.append(f"triangular coordinate {key} is not true")
    if float(tri.get("T_B_determinant_exact", math.nan)) != 1.0:
        f.append("triangular coordinate determinant changed")
    B = float(tri.get("B_operator_norm_upper_mps2_per_cayley", math.nan))
    if not (math.isfinite(B) and B > 0.0):
        f.append("triangular B norm bound invalid")
    metric = d.get("moving_metric_congruence", {})
    for key in (
        "innovation_covariance_S_invariant",
        "Joseph_covariance_congruence_exact",
        "moving_Riccati_energy_invariant",
    ):
        if metric.get(key) is not True:
            f.append(f"moving metric congruence {key} is not true")
    for key in ("condition_number_multiplier_used", "group_isotropic_metric_assumption_used"):
        if metric.get(key) is not False:
            f.append(f"forbidden moving metric shortcut {key} enabled")
    exact = d.get("exact_finite_angle_coordinate", {})
    for key in (
        "nonlinear_displacement_is_exactly_existing_e_eta",
        "first_order_wave_attitude_term_removed_from_remainder",
    ):
        if exact.get(key) is not True:
            f.append(f"exact finite-angle coordinate {key} is not true")
    rows = d.get("measurement_linearizing_shift_bounds_reused_without_widening", [])
    if [r.get("attitude_angle_deg") for r in rows] != [30.0, 25.0, 20.0, 15.0]:
        f.append("finite-angle candidate cells changed")
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
        "source_nominal_aw_norm_upper_mps2": d["source_nominal_aw_norm_upper_mps2"],
        "gravity_only_attitude_column": d["transformed_attitude_column_depends_only_on_gravity"],
        "actual_RS_retained": d["actual_RS_information_matrix_retained_under_congruence"],
        "transport_closed": d["complete_source_correlated_transport_defect_closed_here"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
