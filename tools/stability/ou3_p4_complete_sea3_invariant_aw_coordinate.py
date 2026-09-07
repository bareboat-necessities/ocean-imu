#!/usr/bin/env python3
"""Invariant a_w normal form for complete-SEA3 OU-III P4.

This lemma is a linear, source-indexed congruence on top of the exact
measurement-linearizing coordinate.  It does not change the filter and it does
not replace the complete SEA3 word.

For one accepted accelerometer operation let

    f_hat = R_hat (a_hat - g),
    r_a   = (E-I) f_hat + R_hat u_aw + db_a,

where ``u_aw=R_hat^T E R_hat delta_a_w`` is the exact operation coordinate from
``ou3_p4_complete_sea3_accelerometer_operation_coordinate``.  Introduce the
source-indexed *linear* triangular coordinate

    w_lin = delta_a_w + B c,
    B c   = R_hat^T [c]x R_hat a_hat
          = -R_hat^T [R_hat a_hat]x c.

Then the shipping tangent residual is identically

    [c]x f_hat + R_hat delta_a_w + db_a
      = -[c]x R_hat g + R_hat w_lin + db_a.

Thus the wave-acceleration part of the attitude column is moved into the a_w
coordinate by a unit-triangular state transform.  If

    T_B : (c,delta_a_w) -> (c,w_lin),

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
    epsilon_aw = (Q_E-I)delta_a_w + e_eta,
    w_exact = w_lin + epsilon_aw,

one obtains the exact finite-angle normal form

    r_a = -[c]x R_hat g + R_hat w_exact + db_a = H_B Phi_B(z).

Equivalently, writing a_true=a_hat+delta_a_w and Q_E=R_hat^T E R_hat,

    w_exact = Q_E a_true - a_hat
              - R_hat^T ((E-I)-[c]x) R_hat g.

The only nonlinear displacement from the transformed tangent coordinate is
the full shift epsilon_aw, including the mixed wave-error term. The pure
finite-angle e_eta bounds alone do not bound this displacement.  The large
first-order wave-acceleration/attitude coupling is part of the exact linear
congruence T_B, not a remainder.  Because the whole complete SEA3 word is
transformed rather than altered, every due S=0 update and its actual applied
SpectralMSE R_S remain in the same word and regularize the transformed a_w
direction through the same information matrix.

The nominal-force ceiling used below is NOT inferred from physical SEA3 alone.
The declared Normal-Live source gives the true specific-force bound

    ||a_true-g|| <= 13.80665 m/s^2,

and the already-declared P4/startup perturbation domain gives

    ||delta_a_w|| <= 2.941995 m/s^2.

Since a_hat=a_true-delta_a_w,

    ||a_hat-g|| <= 13.80665+2.941995 = 16.748645 m/s^2,
    ||a_hat||   <= 16.748645+g       = 26.555295 m/s^2.

Those are theorem-domain consequences, not replay extrema and not a new source
assumption.  The pure-e_eta candidate rows inherited from the earlier helper
were computed with the physical 13.80665 ceiling; this module outwardly widens
them by the exact nominal/physical force ratio before they are reused.

This producer deliberately remains fail-closed.  P4 still requires an
operation-matched bound on transport of full epsilon_aw through the physical
correction, quaternion injection/reset, prediction, and source change.  That
transport must be charged to the same Joseph information decrease; no
packet-count-times-worst remainder budget is admissible.
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
SCHEMA = 3
QUALIFICATION = "OU3_P4_COMPLETE_SEA3_INVARIANT_AW_NORMAL_FORM_V3"


def _up_add(a: float, b: float) -> float:
    return math.nextafter(float(a) + float(b), math.inf)


def _widen_candidate_rows(rows: list[dict], physical_force: float, nominal_force: float) -> list[dict]:
    if not (math.isfinite(physical_force) and physical_force > 0.0):
        raise RuntimeError("invalid physical force ceiling")
    if not (math.isfinite(nominal_force) and nominal_force >= physical_force):
        raise RuntimeError("invalid nominal force ceiling")
    ratio = math.nextafter(nominal_force / physical_force, math.inf)
    out: list[dict] = []
    for row in rows:
        r = dict(row)
        for key in (
            "e_eta_norm_upper_mps2",
            "e_eta_local_lipschitz_upper_mps2_per_cayley",
        ):
            value = float(r[key])
            if not (math.isfinite(value) and value > 0.0):
                raise RuntimeError(f"invalid inherited candidate bound {key}")
            r[key] = math.nextafter(value * ratio, math.inf)
        r["force_upper_mps2_used"] = nominal_force
        r["inherited_physical_force_upper_mps2"] = physical_force
        r["nominal_force_widening_ratio"] = ratio
        r["widened_for_declared_delta_aw_domain"] = True
        out.append(r)
    return out


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
    handoff = domain["startup"]["physical_handoff_coordinate_bounds"]
    g = float(domain["startup"]["gravity_mps2"])
    true_force = float(live["specific_force_norm_upper_mps2"])
    delta_aw = float(handoff["latent_acceleration_error_norm_upper_mps2"])
    if not all(math.isfinite(x) and x > 0.0 for x in (g, true_force, delta_aw)):
        raise RuntimeError("invalid gravity/true-force/delta-aw theorem bounds")

    nominal_force = _up_add(true_force, delta_aw)
    nominal_aw = _up_add(nominal_force, g)
    rows = _widen_candidate_rows(awlin["candidate_cells"], true_force, nominal_force)

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
        "nominal_force_derivation": {
            "true_specific_force_upper_mps2": true_force,
            "declared_delta_aw_error_upper_mps2": delta_aw,
            "nominal_specific_force_upper_mps2": nominal_force,
            "gravity_mps2": g,
            "nominal_aw_norm_upper_mps2": nominal_aw,
            "identity": "a_hat=a_true-delta_a_w",
            "triangle_inequality": "||a_hat-g||<=||a_true-g||+||delta_a_w||",
            "uses_complete_SEA3_physical_source_bound": True,
            "uses_already_declared_P4_error_domain": True,
            "new_source_assumption_added": False,
            "replay_extrema_used": False,
        },
        "source_nominal_specific_force_upper_mps2": nominal_force,
        "source_nominal_aw_norm_upper_mps2": nominal_aw,
        "nominal_force_bound_from_complete_SEA3_plus_declared_error_domain_proved": True,
        "nominal_force_bound_from_physical_SEA3_alone_claimed": False,
        "shipping_Joseph_binding_closed": False,
        "shipping_Joseph_binding_scope": "SOURCE_UNIFORM_NONLINEAR_TRANSPORT",
        "linear_triangular_coordinate": {
            "w_lin": "delta_a_w+B*c",
            "B_c": "R_hat^T*[c]x*R_hat*a_hat",
            "B_matrix": "-R_hat^T*[R_hat*a_hat]x",
            "B_operator_norm_upper_mps2_per_cayley": nominal_aw,
            "T_B_unit_triangular": True,
            "T_B_determinant_exact": 1.0,
            "T_B_nonsingular": True,
            "T_B_inverse_exact": True,
        },
        "tangent_normal_form_identity": (
            "[c]x*f_hat+R_hat*delta_a_w+delta_b_a="
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
            "w_exact": "w_lin+epsilon_aw",
            "epsilon_aw": "(Q_E-I)*delta_a_w+e_eta",
            "e_eta": "R_hat^T*((E-I)-[c]x)*f_hat",
            "equivalent_w_exact": (
                "Q_E*a_true-a_hat-R_hat^T*((E-I)-[c]x)*R_hat*g"
            ),
            "Q_E": "R_hat^T*E*R_hat",
            "exact_residual": (
                "r_a=-[c]x*R_hat*g+R_hat*w_exact+delta_b_a=H_B*Phi_B(z)"
            ),
            "nonlinear_displacement_is_full_aw_shift": True,
            "first_order_wave_attitude_term_removed_from_remainder": True,
        },
        "measurement_linearizing_shift_bounds_reused_without_widening": False,
        "measurement_linearizing_shift_bounds_widened_for_declared_delta_aw_domain": rows,
        "reused_candidate_bounds_cover_pure_e_eta_only": True,
        "full_nonlinear_coordinate_displacement_bounded_here": False,
        "standalone_eta_Rinv_packet_budget_used": False,
        "packet_count_multiplier_used": False,
        "actual_RS_information_matrix_retained_under_congruence": True,
        "complete_source_correlated_transport_defect_closed_here": False,
        "P4_promoted_here": False,
        "next_obligation": (
            "bound the full epsilon_aw change, including mixed wave error, under the same Joseph correction/reset and prediction; "
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
        "nominal_force_bound_from_complete_SEA3_plus_declared_error_domain_proved",
        "transformed_attitude_column_depends_only_on_gravity",
        "wave_acceleration_attitude_cross_term_is_linear_coordinate_coupling",
        "actual_RS_information_matrix_retained_under_congruence",
        "reused_candidate_bounds_cover_pure_e_eta_only",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "nominal_force_bound_from_physical_SEA3_alone_claimed",
        "shipping_Joseph_binding_closed",
        "trajectory_replay_used",
        "filter_changed",
        "declared_domain_changed",
        "source_family_replaced",
        "standalone_eta_Rinv_packet_budget_used",
        "packet_count_multiplier_used",
        "complete_source_correlated_transport_defect_closed_here",
        "full_nonlinear_coordinate_displacement_bounded_here",
        "P4_promoted_here",
        "measurement_linearizing_shift_bounds_reused_without_widening",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")

    force = d.get("nominal_force_derivation", {})
    for key in (
        "uses_complete_SEA3_physical_source_bound",
        "uses_already_declared_P4_error_domain",
    ):
        if force.get(key) is not True:
            f.append(f"nominal-force derivation {key} is not true")
    for key in ("new_source_assumption_added", "replay_extrema_used"):
        if force.get(key) is not False:
            f.append(f"nominal-force derivation {key} is not false")
    true_force = float(force.get("true_specific_force_upper_mps2", math.nan))
    delta_aw = float(force.get("declared_delta_aw_error_upper_mps2", math.nan))
    nominal_force = float(force.get("nominal_specific_force_upper_mps2", math.nan))
    g = float(force.get("gravity_mps2", math.nan))
    nominal_aw = float(force.get("nominal_aw_norm_upper_mps2", math.nan))
    if not all(math.isfinite(x) and x > 0.0 for x in (true_force, delta_aw, nominal_force, g, nominal_aw)):
        f.append("nominal-force derivation contains invalid values")
    else:
        if nominal_force < true_force + delta_aw:
            f.append("nominal specific-force ceiling dropped declared delta-aw error")
        if nominal_aw < nominal_force + g:
            f.append("nominal aw ceiling dropped gravity")

    tri = d.get("linear_triangular_coordinate", {})
    for key in ("T_B_unit_triangular", "T_B_nonsingular", "T_B_inverse_exact"):
        if tri.get(key) is not True:
            f.append(f"triangular coordinate {key} is not true")
    if float(tri.get("T_B_determinant_exact", math.nan)) != 1.0:
        f.append("triangular coordinate determinant changed")
    B = float(tri.get("B_operator_norm_upper_mps2_per_cayley", math.nan))
    if not (math.isfinite(B) and B >= nominal_aw):
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
        "nonlinear_displacement_is_full_aw_shift",
        "first_order_wave_attitude_term_removed_from_remainder",
    ):
        if exact.get(key) is not True:
            f.append(f"exact finite-angle coordinate {key} is not true")

    rows = d.get("measurement_linearizing_shift_bounds_widened_for_declared_delta_aw_domain", [])
    if [r.get("attitude_angle_deg") for r in rows] != [30.0, 25.0, 20.0, 15.0]:
        f.append("finite-angle candidate cells changed")
    for row in rows:
        if row.get("widened_for_declared_delta_aw_domain") is not True:
            f.append("finite-angle row not widened for declared delta-aw domain")
        if float(row.get("force_upper_mps2_used", math.nan)) < nominal_force:
            f.append("finite-angle row uses insufficient nominal force ceiling")
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
        "true_specific_force_upper_mps2": d["nominal_force_derivation"]["true_specific_force_upper_mps2"],
        "declared_delta_aw_error_upper_mps2": d["nominal_force_derivation"]["declared_delta_aw_error_upper_mps2"],
        "nominal_specific_force_upper_mps2": d["source_nominal_specific_force_upper_mps2"],
        "source_nominal_aw_norm_upper_mps2": d["source_nominal_aw_norm_upper_mps2"],
        "gravity_only_attitude_column": d["transformed_attitude_column_depends_only_on_gravity"],
        "actual_RS_retained": d["actual_RS_information_matrix_retained_under_congruence"],
        "transport_closed": d["complete_source_correlated_transport_defect_closed_here"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
