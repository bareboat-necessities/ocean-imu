#!/usr/bin/env python3
"""Fail-closed source-uniform P4 closure audit.

This audit separates four logically different obligations:

1. the full declared hard word-entry error set, independent of covariance;
2. conditional BIAS1 and canonical P3 execution-premise admission;
3. source-uniform same-history Kalman/reset coefficients plus consecutive
   compatible-storage/signed-information domination; and
4. finite precision.

Exact-real and binary32 bias projection are now separately enclosed.  The
complete Eigen/Kalman/reset arithmetic remains fail-closed until an event-level
shipping arithmetic enclosure is available.  No fractional entry-domain shrink
or point replay can promote P4.
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
import ou3_p4_projection_binary32_enclosure as PROJFP
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_complete_brmm_differential_word as WORD
import ou3_p4_complete_brmm_finite_map_mean_value as FINITE
import ou3_p4_complete_brmm_signed_information_ledger as SIGNED
import ou3_p4_complete_brmm_joint_sector_master as JOINT
import ou3_p4_strong_linear_margin as STRONG

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
DEFAULT_CLOSURE = REPO / "tools" / "stability" / "ou3_p4_closure_domain.json"
SCHEMA = 3
QUALIFICATION = "OU3_P4_SOURCE_UNIFORM_CLOSURE_AUDIT_V3"


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
    # Position is declared componentwise |dp_i|<=4.25 m.  The closure-domain
    # base stores the implied Euclidean norm sqrt(3)*4.25.
    position_norm = math.sqrt(3.0) * entry["position_component_abs_error_upper_m"]
    return bool(
        math.isclose(entry["attitude_cayley_norm_upper"], float(b["attitude_cayley_norm"]), rel_tol=0.0, abs_tol=2e-16)
        and entry["gyro_bias_error_norm_upper_rad_s"] == float(b["gyro_bias_norm_rad_s"])
        and entry["velocity_error_norm_upper_mps"] == float(b["velocity_norm_mps"])
        and math.isclose(position_norm, float(b["position_norm_m"]), rel_tol=0.0, abs_tol=2e-15)
        and entry["integral_displacement_error_norm_upper_m_s"] == float(b["integral_displacement_norm_m_s"])
        and entry["latent_acceleration_error_norm_upper_mps2"] == float(b["latent_acceleration_norm_mps2"])
        and entry["accelerometer_bias_error_norm_upper_mps2"] == float(b["accelerometer_bias_error_norm_mps2"])
        and float(search["minimum_certified_scale"]) == 1.0
        and list(map(float, search["candidate_scale_factors"])) == [1.0]
        and search["requires_full_declared_scale"] is True
    )


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
    probe = {"root": [0.08,-0.05,0.03], "amplitude": [0.015,0.010,-0.008], "tau_s":1200.0, "period_s":600.0}
    probe_member = bool(
        max(map(abs,probe["root"])) <= root and max(map(abs,probe["amplitude"])) <= amp
        and tau_lo <= probe["tau_s"] <= tau_hi and per_lo <= probe["period_s"] <= per_hi
    )
    return {
        "family": f,
        "conditional_theorem_source_family_well_posed": structurally_closed,
        "existing_nonzero_driver_probe_member": probe_member,
        "one_root_and_driver_history_required": True,
        "independent_per_sample_bias_slots_used": False,
        "conditional_BIAS1_SOURCE_ADMISSION_PASS": structurally_closed and probe_member,
        "assembled_sensor_BIAS0_deployment_qualification_pass": False,
        "scope": "conditional theorem family; not assembled-sensor hardware qualification",
    }


def _p3_conditional_admission(p3: dict, premises: dict, closure: dict) -> dict:
    x=premises["execution_premises"]
    bool_keys=(
        "same_complete_physical_frontend_tuner_geometry_history",
        "declared_Normal_Live_branch_not_just_runtime_Live_flag",
        "lever_arm_disabled_and_vibration_guard_dormant_transparent",
        "accepted_accelerometer_at_every_valid_IMU_sample", "no_accelerometer_rejection_in_word",
        "PE_geometry_attached_to_actual_measurement_jacobians",
        "committed_OU_parameters_obey_shipping_clamp_and_commit_invariants",
        "every_due_S_with_actual_anisotropic_RS_and_full_cross_covariance",
        "shipping_configured_R_full_Q_and_PSD_aw_floors",
        "successful_Joseph_updates_and_immediate_covariance_resets",
        "no_hard_attitude_rewrite_inside_same_mode_word",
        "source_generated_Live_seed_and_shipping_H_to_A_release",
        "premises_continue_after_every_projection_and_hybrid_entry",
    )
    manifest=all(x.get(k) is True for k in bool_keys)
    c=closure["P3_execution_admission"]
    conditional=bool(
        manifest and p3["P3_CONDITIONAL_BRMM_PASS"] and float(p3["useful_gate"])==1e-18
        and float(premises["delta"])==1e-18
        and premises["projection_leaves_covariance_unchanged_source_parity"]
        and c["canonical_source"]==p3["canonical_source"] and c["scope_is_BRMM_family_itself"]
        and c["all_due_S_updates_with_actual_applied_RS_required"]
        and c["all_valid_accelerometer_updates_required"] and c["projection_covariance_identity_required"]
        and c["same_history_frontend_tuner_covariance_required"]
        and c["no_replay_or_finite_harmonic_membership_substitute"])
    return {
        "conditional_execution_premises_admitted":conditional,
        "canonical_P3_delta":1e-18,
        "runtime_Live_flag_is_not_admission":not premises["runtime_Live_flag_implies_all_premises"],
        "BRMM_alone_is_not_execution_admission":not premises["BRMM_alone_implies_vector_PE"],
        "physical_execution_admission_proved":premises["physical_execution_admission_proved_here"],
        "scope":premises["scope"],
    }


def build(domain_path: Path=DEFAULT_DOMAIN, closure_path: Path=DEFAULT_CLOSURE) -> dict:
    domain_path=Path(domain_path).resolve(); closure_path=Path(closure_path).resolve()
    domain=json.loads(domain_path.read_text()); closure=json.loads(closure_path.read_text())
    if domain.get("trajectory_fit") is not False or closure.get("trajectory_fit") is not False:
        raise RuntimeError("P4 closure audit cannot consume a trajectory-fitted domain")

    p3=P3.build(domain_path); premises=PREMISES.build(domain_path); bias=BIAS.build(domain_path)
    projection=PROJ.build(); projection_fp=PROJFP.build(); events=EVENTS.build(domain_path)
    word=WORD.build(domain_path); finite=FINITE.build(domain_path); signed=SIGNED.build(domain_path)
    joint=JOINT.build(p3_contract=p3,signed_contract=signed); strong=STRONG.build(domain_path)
    prereq={
        "P3":P3.validate(p3), "premises":PREMISES.validate(premises,domain_path), "bias":BIAS.validate(bias),
        "projection":PROJ.validate(projection), "projection_fp":PROJFP.validate(projection_fp),
        "events":EVENTS.validate(events), "word":WORD.validate(word), "finite":FINITE.validate(finite),
        "signed":SIGNED.validate(signed), "joint":JOINT.validate(joint), "strong":STRONG.validate(strong)}
    bad={k:v for k,v in prereq.items() if v}
    if bad: raise RuntimeError(f"uniform closure prerequisites failed: {bad}")

    entry=_declared_hard_entry(domain); entry_declared=_hard_entry_matches_closure_contract(entry,closure)
    bias1=_bias1_conditional_admission(closure); p3_admission=_p3_conditional_admission(p3,premises,closure)

    # These are the real word-level obligations. Norm bounds on individual K/G
    # coefficients are useful prerequisites but are not substituted for this
    # same-history generalized-Jacobian/joint-sector closure.
    coefficient_family=bool(
        events["source_uniform_finite_angle_event_Jacobians_closed"]
        and word["source_uniform_complete_word_Jacobian_enclosed"]
        and finite["source_uniform_complete_word_generalized_Jacobian_enclosed"]
        and joint["source_uniform_same_history_joint_sector_closed"]
        and joint["source_uniform_full_augmented_LDLT_closed"])
    consecutive_storage=bool(
        coefficient_family and finite["source_uniform_endpoint_finite_map_closed"]
        and finite["source_uniform_all_prefix_gains_closed"]
        and signed["source_uniform_joint_eta_reset_domination_closed"]
        and signed["source_uniform_endpoint_dissipation_closed_here"]
        and signed["source_uniform_all_prefix_gain_closed_here"]
        and joint["source_uniform_prefix_joint_sector_closed"])

    fp_cfg=closure["finite_precision"]
    projection_fp_closed=bool(projection_fp["projection_finite_precision_enclosure_closed"])
    # Full filter arithmetic remains false by explicit contract, despite the
    # now-closed projection subroutine enclosure.
    finite_precision=bool(
        projection_fp_closed and fp_cfg["complete_Eigen_expression_evaluation_order_enclosed"]
        and fp_cfg["complete_shipping_Kalman_reset_roundoff_enclosed"])

    strong_summary={
        "H18_relative_margin_lower":float(strong["H18"]["strong_relative_margin_lower"]),
        "A21_relative_margin_lower":float(strong["A21"]["strong_relative_margin_lower"]),
        "H18_worst_LDLT_pivot_lower":float(strong["H18"]["worst_interval_LDLT_pivot_lower_at_strong_margin"]),
        "A21_first_active_ba_margin_lower":float(strong["A21"]["first_active_ba_margin_lower"]),
        "scalar_margin_architecture_rejected_for_finite_nonlinear_closure":True}

    blockers=[]
    if not coefficient_family: blockers.append("same-history source-uniform finite-angle Kalman/reset generalized-Jacobian and joint-sector enclosure is not closed")
    if not consecutive_storage: blockers.append("source-uniform endpoint/every-prefix compatible-storage signed-information inequality is not closed")
    if not finite_precision: blockers.append("projection binary32 is enclosed, but complete shipping Eigen/Kalman/reset arithmetic is not")
    if not entry_declared: blockers.append("full declared hard entry set no longer matches theorem operating domain")
    if not bias1["conditional_BIAS1_SOURCE_ADMISSION_PASS"]: blockers.append("conditional BIAS1 driver family is not admitted")
    if not p3_admission["conditional_execution_premises_admitted"]: blockers.append("canonical conditional P3 execution premises are not admitted")

    p4_motion=bool(entry_declared and bias1["conditional_BIAS1_SOURCE_ADMISSION_PASS"]
        and p3_admission["conditional_execution_premises_admitted"] and projection["global_joint_sector_closed"]
        and projection_fp_closed and coefficient_family and consecutive_storage and finite_precision)
    return {
        "schema":SCHEMA,"qualification":QUALIFICATION,"canonical_source":p3["canonical_source"],
        "filter_changed":False,"declared_main_operating_domain_shrunk":False,"quality_gates_changed":False,
        "P3_delta_consumed":1e-18,"hard_entry_error_set":entry,
        "qualified_hard_entry_error_set_as_declared_theorem_assumption":entry_declared,
        "hard_entry_reachability_from_arbitrary_startup_claimed":False,
        "shipping_covariance_used_as_entry_membership_test":False,
        "BIAS1_admission":bias1,"P3_execution_admission":p3_admission,
        "projection_global_nonlinear_sector_consumed":bool(projection["global_joint_sector_closed"]),
        "projection_finite_precision_contract":projection_fp,
        "projection_finite_precision_closed":projection_fp_closed,
        "strong_linear_margin_diagnostic":strong_summary,
        "source_uniform_Kalman_reset_coefficient_family_enclosed":coefficient_family,
        "consecutive_compatible_storage_inequality_closed":consecutive_storage,
        "finite_precision_contract":fp_cfg,
        "complete_filter_finite_precision_enclosure_closed":finite_precision,
        "finite_precision_enclosure_closed":finite_precision,
        "P4_MOTION_PASS":p4_motion,"P4_PASS":p4_motion,"P5_MAY_START":p4_motion,
        "remaining_blockers":blockers}


def validate(d: dict) -> list[str]:
    f=[]
    if d.get("schema")!=SCHEMA or d.get("qualification")!=QUALIFICATION:f.append("schema/qualification mismatch")
    if d.get("canonical_source")!="COMPLETE_BRMM_NORMAL_LIVE_WORD":f.append("canonical source changed")
    if float(d.get("P3_delta_consumed",0))!=1e-18:f.append("canonical P3 delta changed")
    for k in ("qualified_hard_entry_error_set_as_declared_theorem_assumption","projection_global_nonlinear_sector_consumed","projection_finite_precision_closed"):
        if d.get(k) is not True:f.append(k+" is not true")
    for k in ("filter_changed","declared_main_operating_domain_shrunk","quality_gates_changed",
              "hard_entry_reachability_from_arbitrary_startup_claimed","shipping_covariance_used_as_entry_membership_test",
              "source_uniform_Kalman_reset_coefficient_family_enclosed","consecutive_compatible_storage_inequality_closed",
              "complete_filter_finite_precision_enclosure_closed","finite_precision_enclosure_closed",
              "P4_MOTION_PASS","P4_PASS","P5_MAY_START"):
        if d.get(k) is not False:f.append(k+" must remain false until its certificate closes")
    if d.get("BIAS1_admission",{}).get("conditional_BIAS1_SOURCE_ADMISSION_PASS") is not True:f.append("conditional BIAS1 source family not admitted")
    if d.get("P3_execution_admission",{}).get("conditional_execution_premises_admitted") is not True:f.append("conditional P3 execution premises not admitted")
    if not d.get("remaining_blockers"):f.append("fail-closed audit unexpectedly has no blockers")
    return list(dict.fromkeys(f))


def main() -> int:
    ap=argparse.ArgumentParser();ap.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN);ap.add_argument("--closure-domain",type=Path,default=DEFAULT_CLOSURE);ap.add_argument("--output",type=Path,required=True);a=ap.parse_args()
    d=build(a.domain,a.closure_domain);f=validate(d);d["validation_pass"]=not f;d["validation_failures"]=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"hard_entry":d["qualified_hard_entry_error_set_as_declared_theorem_assumption"],"BIAS1":d["BIAS1_admission"]["conditional_BIAS1_SOURCE_ADMISSION_PASS"],"P3_execution":d["P3_execution_admission"]["conditional_execution_premises_admitted"],"projection_fp":d["projection_finite_precision_closed"],"coefficient_family":d["source_uniform_Kalman_reset_coefficient_family_enclosed"],"compatible_storage":d["consecutive_compatible_storage_inequality_closed"],"full_fp":d["finite_precision_enclosure_closed"],"P4_PASS":d["P4_PASS"],"blockers":d["remaining_blockers"],"failures":f},indent=2,sort_keys=True))
    return int(bool(f))
if __name__=="__main__":raise SystemExit(main())
