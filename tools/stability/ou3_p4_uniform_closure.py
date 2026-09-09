#!/usr/bin/env python3
"""Fail-closed source-uniform P4 closure audit.

This audit separates four logically different obligations that were previously
collapsed into one permissive scalar-radius test:

1. conditional admission of the declared hard word-entry error set;
2. conditional admission of the BIAS1 root/driver family and canonical P3
   execution premises;
3. source-uniform same-history Kalman/reset coefficient and signed-ledger
   domination with consecutive compatible storage; and
4. shipping finite-precision enclosure.

The first two can be closed as theorem-domain assumptions without using a
covariance ellipsoid as an error bound.  The latter two remain fail-closed until
existing outward same-history machinery actually certifies them.  In
particular, an arbitrarily tiny positive entry box may never promote P4, and a
binary32 helper for selected arithmetic may never stand in for a full shipping
rounding enclosure.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_brmm_p3_premises as PREMISES
import ou3_brmm_riccati_metric_p3 as P3
import ou3_mems_bias_contract as BIAS
import ou3_p4_projection_sector as PROJ
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_complete_brmm_differential_word as WORD
import ou3_p4_complete_brmm_finite_map_mean_value as FINITE
import ou3_p4_complete_brmm_signed_information_ledger as SIGNED
import ou3_p4_complete_brmm_joint_sector_master as JOINT
import ou3_p4_strong_linear_margin as STRONG

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
DEFAULT_CLOSURE = REPO / "tools" / "stability" / "ou3_p4_closure_domain.json"
SCHEMA = 2
QUALIFICATION = "OU3_P4_SOURCE_UNIFORM_CLOSURE_AUDIT_V2"


def _finite_nonnegative(x, label: str) -> float:
    y = float(x)
    if not math.isfinite(y) or y < 0.0:
        raise RuntimeError(f"{label} must be finite nonnegative")
    return y


def _declared_hard_entry(domain: dict) -> dict:
    startup = domain["startup"]["physical_handoff_coordinate_bounds"]
    entrance = domain["initial_filter_entrance"]
    hs = float(entrance["position"]["significant_wave_height_Hs_upper_m"])
    theta = math.radians(float(entrance["attitude"]["full_attitude_error_upper_deg"]))
    cayley = 2.0 * math.tan(theta / 2.0)
    position_component = float(entrance["position"]["component_abs_error_upper_Hs_factor"]) * hs
    return {
        "qualification": "DECLARED_HARD_P4_WORD_ENTRY_SET",
        "source": "ou3_proof_operating_domain.json",
        "attitude_cayley_norm_upper": cayley,
        "gyro_bias_error_norm_upper_rad_s": float(startup["gyro_bias_error_norm_upper_rad_s"]),
        "velocity_error_norm_upper_mps": float(startup["velocity_error_norm_upper_mps"]),
        "position_component_abs_error_upper_m": position_component,
        "integral_displacement_error_norm_upper_m_s": float(startup["integral_displacement_error_norm_upper_m_s"]),
        "latent_acceleration_error_norm_upper_mps2": float(startup["latent_acceleration_error_norm_upper_mps2"]),
        "accelerometer_bias_error_norm_upper_mps2": float(startup["accelerometer_bias_error_norm_upper_mps2"]),
        "covariance_ellipsoid_used_for_membership": False,
        "replay_fit_used_for_membership": False,
        "declared_theorem_entry_assumption": True,
        "physical_reachability_from_arbitrary_startup_proved_here": False,
    }


def _hard_entry_matches_closure_contract(entry: dict, closure: dict) -> bool:
    search = closure["hard_entry_search"]
    b = search["base_coordinate_radii"]
    checks = (
        math.isclose(entry["attitude_cayley_norm_upper"], float(b["attitude_cayley_norm"]), rel_tol=0.0, abs_tol=2e-16),
        entry["gyro_bias_error_norm_upper_rad_s"] == float(b["gyro_bias_norm_rad_s"]),
        entry["velocity_error_norm_upper_mps"] == float(b["velocity_norm_mps"]),
        entry["position_component_abs_error_upper_m"] <= float(b["position_norm_m"]),
        entry["integral_displacement_error_norm_upper_m_s"] == float(b["integral_displacement_norm_m_s"]),
        entry["latent_acceleration_error_norm_upper_mps2"] == float(b["latent_acceleration_norm_mps2"]),
        entry["accelerometer_bias_error_norm_upper_mps2"] == float(b["accelerometer_bias_error_norm_mps2"]),
        float(search["minimum_certified_scale"]) >= 0.01,
        min(map(float, search["candidate_scale_factors"])) >= float(search["minimum_certified_scale"]),
    )
    return all(checks)


def _bias1_conditional_admission(closure: dict) -> dict:
    f = closure["BIAS1_family"]
    root = _finite_nonnegative(f["root_component_abs_upper_mps2"], "BIAS1 root bound")
    amp = _finite_nonnegative(f["sinusoid_component_amplitude_abs_upper_mps2"], "BIAS1 sinusoid bound")
    det = _finite_nonnegative(f["deterministic_mismatch_component_abs_upper_mps2"], "BIAS1 deterministic bound")
    ddet = _finite_nonnegative(f["deterministic_mismatch_derivative_component_abs_upper_mps3"], "BIAS1 derivative bound")
    tau_lo, tau_hi = map(float, f["tau_true_s"])
    per_lo, per_hi = map(float, f["sinusoid_period_s"])
    phase_lo, phase_hi = map(float, f["phase_rad"])
    structurally_closed = bool(
        root > 0.0 and amp >= 0.0 and det >= 0.0 and ddet >= 0.0
        and 0.0 < tau_lo <= tau_hi and 0.0 < per_lo <= per_hi
        and phase_lo <= -math.pi and phase_hi >= math.pi
        and f["one_root_one_parameter_history_required"]
        and f["independent_per_sample_bias_slots_forbidden"]
        and f["driver_increment_definition"] == "w_i=beta(t_i)-phi_true*beta(t_{i-1})"
    )
    probe = {
        "root": [0.08, -0.05, 0.03],
        "amplitude": [0.015, 0.010, -0.008],
        "tau_s": 1200.0,
        "period_s": 600.0,
    }
    probe_member = bool(
        max(map(abs, probe["root"])) <= root
        and max(map(abs, probe["amplitude"])) <= amp
        and tau_lo <= probe["tau_s"] <= tau_hi
        and per_lo <= probe["period_s"] <= per_hi
    )
    return {
        "family": f,
        "conditional_theorem_source_family_well_posed": structurally_closed,
        "existing_nonzero_driver_probe_member": probe_member,
        "one_root_and_driver_history_required": True,
        "independent_per_sample_bias_slots_used": False,
        "conditional_BIAS1_SOURCE_ADMISSION_PASS": structurally_closed and probe_member,
        "assembled_sensor_BIAS0_deployment_qualification_pass": False,
        "scope": "conditional theorem family; not an assembled-sensor hardware qualification",
    }


def _p3_conditional_admission(p3: dict, premises: dict, closure: dict) -> dict:
    declared = premises["execution_premises"]
    bool_leaves = []
    for value in declared.values():
        if isinstance(value, bool):
            bool_leaves.append(value)
    c = closure["P3_execution_admission"]
    conditional = bool(
        p3["P3_CONDITIONAL_BRMM_PASS"]
        and float(p3["useful_gate"]) == 1.0e-18
        and float(premises["delta"]) == 1.0e-18
        and all(bool_leaves)
        and premises["projection_leaves_covariance_unchanged_source_parity"]
        and c["canonical_source"] == p3["canonical_source"]
        and c["scope_is_BRMM_family_itself"]
        and c["all_due_S_updates_with_actual_applied_RS_required"]
        and c["all_valid_accelerometer_updates_required"]
        and c["projection_covariance_identity_required"]
        and c["same_history_frontend_tuner_covariance_required"]
        and c["no_replay_or_finite_harmonic_membership_substitute"]
    )
    return {
        "conditional_execution_premises_admitted": conditional,
        "canonical_P3_delta": 1.0e-18,
        "runtime_Live_flag_is_not_admission": not premises["runtime_Live_flag_implies_all_premises"],
        "BRMM_alone_is_not_execution_admission": not premises["BRMM_alone_implies_vector_PE"],
        "physical_execution_admission_proved": premises["physical_execution_admission_proved_here"],
        "scope": premises["scope"],
    }


def build(domain_path: Path = DEFAULT_DOMAIN, closure_path: Path = DEFAULT_CLOSURE) -> dict:
    domain_path = Path(domain_path).resolve()
    closure_path = Path(closure_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    closure = json.loads(closure_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False or closure.get("trajectory_fit") is not False:
        raise RuntimeError("P4 closure audit cannot consume a trajectory-fitted domain")

    p3 = P3.build(domain_path)
    premises = PREMISES.build(domain_path)
    bias = BIAS.build(domain_path)
    projection = PROJ.build()
    events = EVENTS.build(domain_path)
    word = WORD.build(domain_path)
    finite = FINITE.build(domain_path)
    signed = SIGNED.build(domain_path)
    joint = JOINT.build(p3_contract=p3, signed_contract=signed)
    strong = STRONG.build(domain_path)

    prereq = {
        "P3": P3.validate(p3),
        "premises": PREMISES.validate(premises, domain_path),
        "bias": BIAS.validate(bias),
        "projection": PROJ.validate(projection),
        "events": EVENTS.validate(events),
        "word": WORD.validate(word),
        "finite": FINITE.validate(finite),
        "signed": SIGNED.validate(signed),
        "joint": JOINT.validate(joint),
        "strong": STRONG.validate(strong),
    }
    bad = {k: v for k, v in prereq.items() if v}
    if bad:
        raise RuntimeError(f"uniform closure prerequisites failed: {bad}")

    entry = _declared_hard_entry(domain)
    entry_declared = _hard_entry_matches_closure_contract(entry, closure)
    bias1 = _bias1_conditional_admission(closure)
    p3_admission = _p3_conditional_admission(p3, premises, closure)

    coefficient_family = bool(
        events["source_uniform_finite_angle_event_Jacobians_closed"]
        and word["source_uniform_complete_word_Jacobian_enclosed"]
        and finite["source_uniform_complete_word_generalized_Jacobian_enclosed"]
        and joint["source_uniform_same_history_joint_sector_closed"]
        and joint["source_uniform_full_augmented_LDLT_closed"]
    )
    consecutive_storage = bool(
        coefficient_family
        and finite["source_uniform_endpoint_finite_map_closed"]
        and finite["source_uniform_all_prefix_gains_closed"]
        and signed["source_uniform_joint_eta_reset_domination_closed"]
        and signed["source_uniform_endpoint_dissipation_closed_here"]
        and signed["source_uniform_all_prefix_gain_closed_here"]
        and joint["source_uniform_prefix_joint_sector_closed"]
    )

    fp_cfg = closure["finite_precision"]
    # The repository has an outward binary32 helper, but it explicitly states
    # that it does not cover compiler reassociation/FMA and it is not wired over
    # the complete Eigen/Kalman/reset execution.  Keep this false until a full
    # event-by-event shipping arithmetic enclosure exists.
    finite_precision = False

    strong_summary = {
        "H18_relative_margin_lower": float(strong["H18"]["strong_relative_margin_lower"]),
        "A21_relative_margin_lower": float(strong["A21"]["strong_relative_margin_lower"]),
        "H18_worst_LDLT_pivot_lower": float(strong["H18"]["worst_interval_LDLT_pivot_lower_at_strong_margin"]),
        "A21_first_active_ba_margin_lower": float(strong["A21"]["first_active_ba_margin_lower"]),
        "scalar_margin_architecture_rejected_for_finite_nonlinear_closure": True,
    }

    blockers = []
    if not coefficient_family:
        blockers.append("same-history source-uniform finite-angle Kalman/reset coefficient and joint-sector enclosure is not closed")
    if not consecutive_storage:
        blockers.append("source-uniform endpoint/every-prefix compatible-storage inequality is not closed")
    if not finite_precision:
        blockers.append("complete shipping binary32/Eigen Kalman-reset arithmetic is not outwardly enclosed")
    if not entry_declared:
        blockers.append("declared hard entry set no longer matches the theorem operating domain")
    if not bias1["conditional_BIAS1_SOURCE_ADMISSION_PASS"]:
        blockers.append("conditional BIAS1 driver family is not admitted")
    if not p3_admission["conditional_execution_premises_admitted"]:
        blockers.append("canonical conditional P3 execution premises are not admitted")

    p4_motion = bool(
        entry_declared
        and bias1["conditional_BIAS1_SOURCE_ADMISSION_PASS"]
        and p3_admission["conditional_execution_premises_admitted"]
        and projection["global_joint_sector_closed"]
        and coefficient_family
        and consecutive_storage
        and finite_precision
    )

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": p3["canonical_source"],
        "filter_changed": False,
        "declared_main_operating_domain_shrunk": False,
        "quality_gates_changed": False,
        "P3_delta_consumed": 1.0e-18,
        "hard_entry_error_set": entry,
        "qualified_hard_entry_error_set_as_declared_theorem_assumption": entry_declared,
        "hard_entry_reachability_from_arbitrary_startup_claimed": False,
        "shipping_covariance_used_as_entry_membership_test": False,
        "BIAS1_admission": bias1,
        "P3_execution_admission": p3_admission,
        "projection_global_nonlinear_sector_consumed": bool(projection["global_joint_sector_closed"]),
        "projection_finite_precision_closed": bool(projection["finite_precision_closed_here"]),
        "strong_linear_margin_diagnostic": strong_summary,
        "source_uniform_Kalman_reset_coefficient_family_enclosed": coefficient_family,
        "consecutive_compatible_storage_inequality_closed": consecutive_storage,
        "finite_precision_contract": fp_cfg,
        "outward_binary32_helper_available_but_not_complete_shipping_enclosure": True,
        "finite_precision_enclosure_closed": finite_precision,
        "P4_MOTION_PASS": p4_motion,
        "P4_PASS": p4_motion,
        "P5_MAY_START": p4_motion,
        "remaining_blockers": blockers,
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_BRMM_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    if float(d.get("P3_delta_consumed", 0.0)) != 1.0e-18:
        f.append("canonical P3 delta changed")
    for key in (
        "qualified_hard_entry_error_set_as_declared_theorem_assumption",
        "projection_global_nonlinear_sector_consumed",
        "outward_binary32_helper_available_but_not_complete_shipping_enclosure",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "filter_changed", "declared_main_operating_domain_shrunk", "quality_gates_changed",
        "hard_entry_reachability_from_arbitrary_startup_claimed",
        "shipping_covariance_used_as_entry_membership_test",
        "projection_finite_precision_closed",
        "source_uniform_Kalman_reset_coefficient_family_enclosed",
        "consecutive_compatible_storage_inequality_closed",
        "finite_precision_enclosure_closed", "P4_MOTION_PASS", "P4_PASS", "P5_MAY_START",
    ):
        if d.get(key) is not False:
            f.append(f"{key} must remain false until its certificate closes")
    if d.get("BIAS1_admission", {}).get("conditional_BIAS1_SOURCE_ADMISSION_PASS") is not True:
        f.append("conditional BIAS1 source family not admitted")
    if d.get("P3_execution_admission", {}).get("conditional_execution_premises_admitted") is not True:
        f.append("conditional P3 execution premises not admitted")
    if not d.get("remaining_blockers"):
        f.append("fail-closed audit unexpectedly has no blockers")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--closure-domain", type=Path, default=DEFAULT_CLOSURE)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain, args.closure_domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "hard_entry_declared": d["qualified_hard_entry_error_set_as_declared_theorem_assumption"],
        "BIAS1_conditional_admission": d["BIAS1_admission"]["conditional_BIAS1_SOURCE_ADMISSION_PASS"],
        "P3_conditional_execution_admission": d["P3_execution_admission"]["conditional_execution_premises_admitted"],
        "coefficient_family_closed": d["source_uniform_Kalman_reset_coefficient_family_enclosed"],
        "compatible_storage_closed": d["consecutive_compatible_storage_inequality_closed"],
        "finite_precision_closed": d["finite_precision_enclosure_closed"],
        "P4_PASS": d["P4_PASS"],
        "blockers": d["remaining_blockers"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
