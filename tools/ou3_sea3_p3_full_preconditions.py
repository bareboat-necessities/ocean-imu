#!/usr/bin/env python3
"""Complete Normal-Live source/precondition contract for canonical OU-III P3.

P3 may promote only one same-mode H18/A21 word that carries the complete
measurement-only front end, adaptive candidate/commit source, scheduler,
measurement covariances, process matrices, covariance-floor events, vector PE,
and declared Normal-Live physical bounds.  Source defaults are followed through
the actual wrapper override chain; an inner default is never substituted for a
value overwritten by the deployed outer Config.

The final numerical object is the joint (P,Psi,Omega) Riccati word.  Reduced
information/covariance products are useful diagnostics/components only and may
not promote P3.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import ou3_full_process_ucc as PROCESS
import ou3_p3_pseudo_scheduler_progress_certificate as SCHED
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_sea3_physical_admissibility as PHYSICAL
import ou3_sea3_rs_word_information as RSWORD
import ou3_sea3_windowed_vector_pe as WINDOWPE
import ou3_source_domain_contract as SOURCE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
WAVE_PERIOD = REPO / "src" / "tuner" / "WavePeriodEstimator.h"
VERTICAL = REPO / "src" / "tuner" / "VerticalAccelComplementary.h"
BANDPASS = REPO / "src" / "tuner" / "AdaptiveWaveBandPass.h"
AUTOTUNER = REPO / "src" / "tuner" / "SeaStateAutoTuner.h"
PAPER = REPO / "doc" / "kalman_ou_iii" / "w3d-iss-stability.tex-part"
SCHEMA = 4
QUALIFICATION = "OU3_SEA3_P3_COMPLETE_NORMAL_LIVE_PRECONDITION_CONTRACT"
USEFUL_GATE = 1.0e-18


def _vec3_default(text: str, name: str) -> list[float]:
    m = re.search(
        rf"Eigen::Vector3f\s+{re.escape(name)}\s*=\s*Eigen::Vector3f\(\s*"
        r"([0-9.eE+-]+)f\s*,\s*([0-9.eE+-]+)f\s*,\s*([0-9.eE+-]+)f\s*\)",
        text,
    )
    if not m:
        raise RuntimeError(f"cannot extract configured {name}")
    out = [float(m.group(i)) for i in range(1, 4)]
    if any(not (math.isfinite(x) and x > 0.0) for x in out):
        raise RuntimeError(f"configured {name} is not strictly positive")
    return out


def _config_scalar_default(text: str, name: str) -> float:
    m = re.search(rf"\bfloat\s+{re.escape(name)}\s*=\s*([0-9.eE+-]+)f\s*;", text)
    if not m:
        raise RuntimeError(f"cannot extract outer Config default {name}")
    value = float(m.group(1))
    if not math.isfinite(value):
        raise RuntimeError(f"outer Config default {name} is non-finite")
    return value


def _all(text: str, markers: tuple[str, ...]) -> bool:
    return all(marker in text for marker in markers)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("canonical P3 domain may not be trajectory fitted")
    runtime = domain["configured_runtime"]
    startup = domain["startup"]
    live = domain["normal_live"]

    dynamic = DYNAMIC.build(path)
    sched = SCHED.build(path)
    process = PROCESS.build()
    physical = PHYSICAL.build(path)
    rsword = RSWORD.build(path)
    pe = WINDOWPE.build(path)
    source = SOURCE.build(SOURCE.DEFAULT_HEADER.resolve())
    prereq = {
        "dynamic": DYNAMIC.validate(dynamic),
        "scheduler": SCHED.validate(sched),
        "process": PROCESS.validate(process),
        "physical": PHYSICAL.validate(physical),
        "four_S": RSWORD.validate(rsword),
        "windowed_PE": WINDOWPE.validate(pe),
    }
    bad = {k: v for k, v in prereq.items() if v}
    if bad:
        raise RuntimeError(f"complete P3 prerequisite validation failed: {bad}")

    wrapper = WRAPPER.read_text(encoding="utf-8")
    mekf = MEKF.read_text(encoding="utf-8")
    wave_period = WAVE_PERIOD.read_text(encoding="utf-8")
    vertical = VERTICAL.read_text(encoding="utf-8")
    bandpass = BANDPASS.read_text(encoding="utf-8")
    autotuner = AUTOTUNER.read_text(encoding="utf-8")
    paper = PAPER.read_text(encoding="utf-8")

    sigma_a = _vec3_default(wrapper, "sigma_a")
    sigma_m = _vec3_default(wrapper, "sigma_m")
    Racc = [x * x for x in sigma_a]
    Rmag = [x * x for x in sigma_m]
    inner_warmup = float(source["timing_constants_s"]["online_tune_warmup"])
    outer_warmup = _config_scalar_default(wrapper, "online_tune_warmup_sec")
    domain_warmup = float(startup["online_tune_warmup_sec"])

    floor_pos = mekf.find("apply_pending_aw_covariance_inflation_();")
    s_pos = mekf.find("applyIntegralZeroPseudoMeas();", floor_pos)
    source_parity = {
        "candidate_tau_each_valid_sample": "tune_.tau_applied   += alpha" in wrapper,
        "candidate_sigma_each_valid_sample": "tune_.sigma_applied += alpha" in wrapper,
        "candidate_RS_each_valid_sample": "tune_.RS_applied    += alpha_RS" in wrapper,
        "next_sample_schedule_commit": "void apply_pending_online_tune_()" in wrapper,
        "pseudo_cadence_from_committed_tau": "apply_pseudo_update_cadence_();" in wrapper,
        "outer_config_forwards_warmup": (
            "impl_.setOnlineTuneWarmupSec(cfg_.online_tune_warmup_sec);" in wrapper
        ),
        "inner_warmup_setter_present": "void setOnlineTuneWarmupSec(float warmup_sec)" in wrapper,
        "periodic_aw_floor_tick": "periodic_aw_cov_sync_tick_();" in wrapper,
        "aw_floor_default_is_PSD_increment": "synchronize_aw_covariance_to_stationary();" in wrapper,
        "aw_floor_applied_inside_prediction": floor_pos >= 0,
        "aw_floor_adds_PSD_increment": "Pext.template block<3,3>(OFF_AW, OFF_AW) += Delta;" in mekf,
        "S_update_runs_after_floor_in_prediction": floor_pos >= 0 and s_pos > floor_pos,
        "accelerometer_attitude_jacobian": "const Matrix3 J_att = -skew_symmetric_matrix(f_cog_b);" in mekf,
        "accelerometer_aw_jacobian": "const Matrix3 J_aw  =  R_wb();" in mekf,
        "accelerometer_bias_jacobian_active": "PCt.noalias() += P_all_ba; // J_ba = I" in mekf,
        "full_Joseph_covariance_update": "joseph_update3_(K, S_mat, PCt);" in mekf,
        "Racc_from_configured_std": "Racc(sigma_a.array().square().matrix().asDiagonal())" in mekf,
        "Rmag_from_configured_std": "Rmag(sigma_m.array().square().matrix().asDiagonal())" in mekf,
        "runtime_Racc_setter_from_std": "Racc = sigma_acc.array().square().matrix().asDiagonal();" in mekf,
        "runtime_Rmag_setter_from_std": "Rmag = sigma_mag.array().square().matrix().asDiagonal();" in mekf,
        "dormant_guard_restores_nominal_Racc": _all(
            wrapper, ("if (!(excess > 0.0f))", "mekf_->set_Racc_std(base);", "racc_effective_ = base;")
        ),
        "magnetic_reference_is_explicit_source_state": _all(
            wrapper, ("void setMagWorldRef_", "mag_world_ref_uT_", "mag_world_ref_valid_")
        ) and "void set_mag_world_ref(const Vector3& B_world)" in mekf,
        "acc_factorization_failure_is_explicit_rejection": _all(
            mekf, ("if (!safe_ldlt3_(S_mat, ldlt, Racc.norm()))", "last_acc_diag_.accepted = false;")
        ),
        "mag_factorization_failure_is_explicit_rejection": _all(
            mekf, ("if (!safe_ldlt3_(S_mat, ldlt, Rmag.norm()))", "last_mag_diag_.accepted = false;")
        ),
        "H_A_switch_explicit": "void set_acc_bias_updates_enabled(bool en)" in mekf,
        "pseudo_scheduler_progress_state_explicit": "pseudo_update_elapsed_s_" in mekf,
        "aw_floor_pending_state_explicit": "aw_covariance_floor_pending_" in mekf,
    }
    source_failures = [name for name, ok in source_parity.items() if not ok]

    front_end_parity = {
        "vertical_measurement_only_observer": _all(
            vertical, ("ahrs_", "elapsed_sec_", "up_ms2_", "verticalAccelUpMs2()")
        ),
        "wave_period_full_moment_state": _all(
            wave_period,
            (
                "high_pass_1_", "high_pass_2_", "velocity_", "elevation_",
                "velocity_mean_", "velocity_sq_", "elevation_mean_", "elevation_sq_",
                "weight_", "raw_period_sec_", "log_period_sec_", "usable_period_",
            ),
        ),
        "wave_period_one_way_usable_latch": "if (usable_period_) return;" in wave_period,
        "adaptive_bandpass_state": _all(
            bandpass, ("lowpass_low_", "band_", "p00_", "p01_", "p11_")
        ),
        "acceleration_moment_state": _all(
            autotuner, ("A_mean", "A_sq", "frequency_hz", "tau_var_sec")
        ),
        "same_weight_first_second_moments": _all(
            autotuner, ("A_mean.update(accel, alpha_var);", "A_sq.update(accel * accel, alpha_var);")
        ),
        "no_second_tuner_frequency_smoother": _all(
            autotuner, ("frequency smoothing no longer occurs here", "frequency_hz = f_eff;")
        ),
    }
    front_end_failures = [name for name, ok in front_end_parity.items() if not ok]

    paper_parity = {
        "source_reachable_schedule": "strict source-reachable family" in paper,
        "four_S_translation_observability": "Spread-selected $S$ observability" in paper,
        "windowed_eta6_PE": "\\alpha_6\\mat I_6\\preceq" in paper,
        "windowed_eta9_PE": "\\alpha_9\\mat I_9\\preceq" in paper,
        "asynchronous_PE": "finite-window asynchronous conditions" in paper,
        "mag_rejections_allowed": "rejected samples, and short outages are admissible" in paper,
        "finite_bias_correlation_route": "finite residual-bias correlation time" in paper,
        "positive_measurement_covariances": "positive accepted measurement covariances" in paper,
        "no_hard_attitude_rewrite": "no hard attitude rewrite inside a normal-Live interval" in paper,
        "full_state_block_composition": "block-lower-triangular observation operator" in paper,
    }
    paper_failures = [name for name, ok in paper_parity.items() if not ok]

    dynamic_parity = dynamic["source_parity"]
    rs_parity = rsword["source_parity"]
    joint_schedule = (
        dynamic_parity["tau_sigma_candidates_smoothed_each_valid_sample"]
        and dynamic_parity["RS_candidate_smoothed_each_valid_sample"]
        and dynamic_parity["active_schedule_commit_is_next_sample_predictable"]
        and dynamic_parity["pseudo_cadence_is_same_tau_lipschitz_image"]
        and rs_parity["spectral_mse_uses_realized_pseudo_period"]
        and rs_parity["applied_tau_and_RS_have_distinct_emas"]
    )
    axis_factors = [float(x) for x in rsword["R_S_axis_std_factors"]]
    per_axis_rs = len(axis_factors) == 3 and min(axis_factors) > 0.0

    pe_horizon = float(pe["spread_occurrence_selection"]["word_horizon_s"])
    s_horizon = float(rsword["word_horizon_s_upper"])
    horizon = max(pe_horizon, s_horizon)
    dt = float(runtime["imu_dt_s"])
    samples = int(math.ceil(horizon / dt)) + 1
    hybrid = list(source["hybrid_obligations"])

    startup_runtime = {
        "inner_filter_default_warmup_s": inner_warmup,
        "outer_config_default_warmup_s": outer_warmup,
        "declared_configured_warmup_s": domain_warmup,
        "outer_overrides_inner": (
            source_parity["outer_config_forwards_warmup"]
            and source_parity["inner_warmup_setter_present"]
            and outer_warmup != inner_warmup
        ),
        "effective_configured_value_is_outer_default": outer_warmup,
    }
    measurement_runtime = {
        "accelerometer_std_mps2": sigma_a,
        "accelerometer_variance_diag": Racc,
        "magnetometer_std_uT": sigma_m,
        "magnetometer_variance_diag": Rmag,
        "configured_defaults_source_bound": True,
        "external_constructor_override_is_outside_this_configured_certificate": True,
        "dormant_vibration_guard_preserves_accelerometer_covariance": source_parity[
            "dormant_guard_restores_nominal_Racc"
        ],
    }
    front_end_manifest = {
        "vertical_leveling": ["Mahony state", "elapsed_sec", "vertical_accel_up"],
        "wave_period": [
            "high_pass_1", "high_pass_2", "velocity", "elevation", "velocity_mean",
            "velocity_sq", "elevation_mean", "elevation_sq", "weight", "raw_period",
            "log_period", "usable_period_latch",
        ],
        "wave_band": ["lowpass_low", "band", "p00", "p01", "p11"],
        "acceleration_statistics": ["A_mean(value,weight)", "A_sq(value,weight)", "frequency", "variance_horizon"],
        "candidate_tuner": list(dynamic["adaptive_state"]),
        "active_schedule": ["tau", "sigma_aw", "R_S", "pseudo_period", "commit_phase"],
        "scheduler": ["pseudo_update_elapsed"],
        "covariance_floor": ["aw_covariance_floor_pending", "aw_covariance_floor_target"],
        "magnetic_gauge": ["mag_world_ref", "mag_world_ref_valid", "reference refinement state"],
    }

    mandatory = {
        "configured_runtime_dt_matches_source": math.isclose(
            float(source["configured_runtime_assumption"]["imu_dt_s"]), dt,
            rel_tol=0.0, abs_tol=0.0,
        ),
        "outer_warmup_override_chain_retained": startup_runtime["outer_overrides_inner"],
        "configured_warmup_matches_effective_outer_default": math.isclose(
            domain_warmup, outer_warmup, rel_tol=0.0, abs_tol=0.0
        ),
        "runtime_zero_lever_arm": runtime["imu_lever_arm_enabled"] is False,
        "runtime_dormant_transparent_vibration_guard": (
            runtime["accelerometer_vibration_guard_proof_branch"] == "dormant_transparent"
        ),
        "configured_measurement_covariances_positive": (
            min(Racc) > 0.0 and min(Rmag) > 0.0
            and source_parity["Racc_from_configured_std"]
            and source_parity["Rmag_from_configured_std"]
        ),
        "all_valid_accelerometer_updates_accepted": (
            live["accelerometer_update_required_each_valid_imu_sample_after_live_entry"] is True
            and live["accelerometer_rejection_in_normal_live_scope"] is False
        ),
        "acc_factorization_failure_excluded_by_theorem_domain": (
            source_parity["acc_factorization_failure_is_explicit_rejection"]
            and live["accelerometer_rejection_in_normal_live_scope"] is False
        ),
        "joint_SEA3_adaptive_state": dynamic["P2_DYNAMIC_SOURCE_CERTIFICATE"] == "PASS",
        "complete_front_end_generator_state_retained": not front_end_failures,
        "source_commit_phase_retained": "pending_commit_progress" in dynamic["adaptive_state"],
        "pseudo_scheduler_progress_and_recurrence": sched["scheduler_recurrence_certificate"] is True,
        "applied_RS_tau_sigma_TS_joint_source_retained": joint_schedule,
        "per_axis_RS_factors_retained": per_axis_rs,
        "measurement_covariances_retained_in_updates": all(
            source_parity[k] for k in (
                "Racc_from_configured_std", "Rmag_from_configured_std",
                "runtime_Racc_setter_from_std", "runtime_Rmag_setter_from_std",
            )
        ),
        "magnetic_reference_path_retained_not_frozen": (
            source_parity["magnetic_reference_is_explicit_source_state"]
            and "magnetic_gauge" in source["discrete_source_branches"]
        ),
        "no_hard_attitude_rewrite_inside_same_mode_word": (
            paper_parity["no_hard_attitude_rewrite"]
            and "tilt_reset" in hybrid and "tilt_relock" in hybrid
            and live.get("hard_attitude_rewrite_inside_word") is False
        ),
        "dimension_changing_H_to_A_transition_separate": (
            "held_to_active" in hybrid and source_parity["H_A_switch_explicit"]
        ),
        "four_S_information_retained": (
            rsword["P3_RS_TRANSLATION_OBSERVATION_GEOMETRY_CLOSED"] is True
            and rsword["P3_RS_TRANSLATION_INFORMATION_MATRIX_CLOSED"] is True
        ),
        "full_process_UCC_retained": process["full_process_ucc_pass"] is True,
        "aw_covariance_floor_events_retained": all(
            source_parity[k] for k in (
                "periodic_aw_floor_tick", "aw_floor_default_is_PSD_increment",
                "aw_floor_applied_inside_prediction", "aw_floor_adds_PSD_increment",
                "S_update_runs_after_floor_in_prediction", "aw_floor_pending_state_explicit",
            )
        ),
        "full_accelerometer_cross_block_information_retained": all(
            source_parity[k] for k in (
                "accelerometer_attitude_jacobian", "accelerometer_aw_jacobian",
                "accelerometer_bias_jacobian_active", "full_Joseph_covariance_update",
            )
        ),
        "windowed_asynchronous_eta6_information_retained": pe["pass"] is True,
        "A_mode_finite_bias_correlation_retained": pe["A_mode_bias_route"][
            "uses_eta6_plus_finite_bias_correlation"
        ] is True,
        "specific_force_bounds_retained": (
            float(live["specific_force_norm_lower_mps2"]) > 0.0
            and float(live["specific_force_norm_upper_mps2"])
            >= float(live["specific_force_norm_lower_mps2"])
        ),
        "magnetic_norm_bounds_retained": (
            float(live["magnetic_vector_norm_lower_uT"]) > 0.0
            and float(live["magnetic_vector_norm_upper_uT"])
            >= float(live["magnetic_vector_norm_lower_uT"])
        ),
        "vector_sine_separation_retained": float(live["vector_sine_separation_lower"]) > 0.0,
        "body_rate_bound_retained": float(live["body_rate_norm_upper_deg_s"]) > 0.0,
        "A_mode_bias_domain_retained": (
            float(live["active_accelerometer_bias_state_norm_upper_mps2"])
            < float(live["active_accelerometer_bias_projection_limit_mps2"])
        ),
        "SEA3_height_period_partition_coupling_retained": physical["three_partition_contract"][
            "independent_H_r_and_T_p_rectangular_extrema_forbidden"
        ] is True,
        "stochastic_failure_budget_retained_as_forcing_not_pruning": (
            0.0 < float(domain["stochastic"]["finite_horizon_failure_probability_budget"]) < 1.0
        ),
    }

    numeric = {
        "common_word_horizon_s": horizon,
        "common_word_samples_upper": samples,
        "H_dimension": 18,
        "A_dimension": 21,
        "one_common_event_word_required_for_H_and_A": True,
        "same_joint_source_path_feeds_F_Q_TS_RS": True,
        "same_front_end_state_path_generates_all_tuner_targets": True,
        "same_event_word_contains_accel_S_PE_and_aw_floor": True,
        "same_runtime_measurement_covariances_used_in_Riccati_updates": True,
        "every_valid_accelerometer_update_must_be_applied": True,
        "every_due_S_update_must_be_applied": True,
        "actual_applied_per_axis_RS_required": True,
        "windowed_PE_information_must_be_accumulated_not_ODR_substituted": True,
        "aw_floor_must_be_added_as_actual_PSD_event_not_marginal_Loewner_shortcut": True,
        "A_mode_bias_must_use_finite_tau_or_eta9_information": True,
        "full_18x18_and_21x21_matrix_comparison_required": True,
        "joint_P_Psi_Omega_propagation_required": True,
        "exact_decomposition_identity_required": "P_k = Psi_k P_0 Psi_k^T + Omega_k",
        "prediction_recursion_required": "P-=F P F^T+Q; Psi-=F Psi; Omega-=F Omega F^T+Q",
        "joseph_measurement_recursion_required": (
            "A=I-KH; P+=A P- A^T+K R K^T; Psi+=A Psi-; Omega+=A Omega- A^T+K R K^T"
        ),
        "aw_floor_recursion_required": (
            "P+=P+E_aw Delta E_aw^T; Psi unchanged; Omega+=Omega+E_aw Delta E_aw^T"
        ),
        "required_final_inequality": "Omega_W - delta * P_W >= 0 on full H18/A21 coordinates",
        "moving_metric_equivalence": (
            "Omega_W >= delta P_W iff Psi_W^T P_W^-1 Psi_W <= (1-delta) P_0^-1"
        ),
        "D_W_L_W_split_for_final_gate_forbidden": True,
        "zero_start_Riccati_concavity_replacement_forbidden": True,
        "blockwise_minimum_ratio_for_final_gate_forbidden": True,
        "determinant_trace_scalarization_for_final_gate_forbidden": True,
        "scalar_information_beta_for_final_gate_forbidden": True,
        "independent_tau_sigma_RS_TS_extrema_product_forbidden": True,
        "hardware_magnetometer_ODR_as_PE_recurrence_forbidden": True,
        "old_P2_graph_or_predecessor_enumeration_forbidden": True,
        "hard_attitude_rewrite_inside_word_forbidden": True,
        "front_end_state_freezing_to_replay_value_forbidden": True,
        "stochastic_noise_realization_used_as_homogeneous_pruning": False,
        "useful_gate": USEFUL_GATE,
    }

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "trajectory_fit": False,
        "canonical_architecture": "SEA3_FULL_NORMAL_LIVE_RICCATI_WORD",
        "mandatory_preconditions": mandatory,
        "all_current_machine_checkable_preconditions_present": all(mandatory.values()),
        "source_parity": source_parity,
        "source_parity_failures": source_failures,
        "front_end_state_parity": front_end_parity,
        "front_end_state_parity_failures": front_end_failures,
        "paper_parity": paper_parity,
        "paper_parity_failures": paper_failures,
        "startup_runtime": startup_runtime,
        "measurement_runtime": measurement_runtime,
        "front_end_state_manifest": front_end_manifest,
        "hybrid_obligations": hybrid,
        "R_S_is_credited_but_not_used_as_a_substitute_for_other_preconditions": True,
        "process_UCC_is_credited_but_not_used_as_a_substitute_for_measurement_preconditions": True,
        "windowed_PE_is_credited_but_not_used_as_a_substitute_for_translation_preconditions": True,
        "D_W_L_W_split_is_not_canonical_final_gate": True,
        "word": {
            "horizon_s": horizon,
            "samples_upper": samples,
            "windowed_PE": pe,
            "four_S_subcertificate": {
                "minimum_word_horizon_s": s_horizon,
                "scheduler_gap_s": sched["certified_uniform_max_gap_s"],
                "newton_coordinate_information": rsword["newton_coordinate_information"],
            },
            "aw_covariance_floor_gap_s_upper": dynamic[
                "validated_rate_and_jump_bounds"
            ]["active_commit_gap_s_upper"],
            "full_process_modes": process["modes"],
            "startup_runtime": startup_runtime,
            "measurement_runtime": measurement_runtime,
            "front_end_state_manifest": front_end_manifest,
            "same_mode_hybrid_exclusions": [
                "held_to_active", "tilt_reset", "tilt_relock", "cooldown_reentry"
            ],
            "magnetic_reference_treatment": (
                "retain source-generated reference/vector path in H_k; gauge variation is ISS forcing, not a frozen replay constant"
            ),
            "stochastic_treatment": (
                "configured R matrices stay in covariance recursion; realization noise belongs to ISS forcing"
            ),
        },
        "final_numeric_contract": numeric,
        "physical_SEA3_scope": {
            "parameter_height_period_partition_coupling_consumed": True,
            "global_finite_window_realization_left_inclusion_closed": False,
            "unqualified_RAO_coupling_used_as_hard_pruning": False,
            "canonical_P3_is_conditional_on_admitted_Normal_Live_SEA3_word": True,
        },
        "P3_promoted": False,
        "P4_may_consume": False,
        "next_obligation": (
            "propagate one exact joint (P,Psi,Omega) H18/A21 word using this complete source/event state; no reduced replacement may promote P3"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("trajectory_fit") is not False:
        f.append("complete P3 contract became trajectory fitted")
    if d.get("canonical_architecture") != "SEA3_FULL_NORMAL_LIVE_RICCATI_WORD":
        f.append("canonical architecture changed")
    for name, ok in d.get("mandatory_preconditions", {}).items():
        if ok is not True:
            f.append(f"mandatory precondition missing: {name}")
    if d.get("all_current_machine_checkable_preconditions_present") is not True:
        f.append("not all machine-checkable preconditions are present")
    f.extend(f"source parity failed: {x}" for x in d.get("source_parity_failures", []))
    f.extend(f"front-end parity failed: {x}" for x in d.get("front_end_state_parity_failures", []))
    f.extend(f"paper parity failed: {x}" for x in d.get("paper_parity_failures", []))

    startup_runtime = d.get("startup_runtime", {})
    if float(startup_runtime.get("inner_filter_default_warmup_s", math.nan)) != 5.0:
        f.append("inner filter warmup default changed")
    if float(startup_runtime.get("outer_config_default_warmup_s", math.nan)) != 10.0:
        f.append("outer configured warmup default changed")
    if startup_runtime.get("outer_overrides_inner") is not True:
        f.append("outer warmup override edge was lost")

    c = d.get("final_numeric_contract", {})
    for key in (
        "one_common_event_word_required_for_H_and_A",
        "same_joint_source_path_feeds_F_Q_TS_RS",
        "same_front_end_state_path_generates_all_tuner_targets",
        "same_event_word_contains_accel_S_PE_and_aw_floor",
        "same_runtime_measurement_covariances_used_in_Riccati_updates",
        "every_valid_accelerometer_update_must_be_applied",
        "every_due_S_update_must_be_applied",
        "actual_applied_per_axis_RS_required",
        "windowed_PE_information_must_be_accumulated_not_ODR_substituted",
        "aw_floor_must_be_added_as_actual_PSD_event_not_marginal_Loewner_shortcut",
        "A_mode_bias_must_use_finite_tau_or_eta9_information",
        "full_18x18_and_21x21_matrix_comparison_required",
        "joint_P_Psi_Omega_propagation_required",
        "D_W_L_W_split_for_final_gate_forbidden",
        "zero_start_Riccati_concavity_replacement_forbidden",
        "blockwise_minimum_ratio_for_final_gate_forbidden",
        "determinant_trace_scalarization_for_final_gate_forbidden",
        "scalar_information_beta_for_final_gate_forbidden",
        "independent_tau_sigma_RS_TS_extrema_product_forbidden",
        "hardware_magnetometer_ODR_as_PE_recurrence_forbidden",
        "old_P2_graph_or_predecessor_enumeration_forbidden",
        "hard_attitude_rewrite_inside_word_forbidden",
        "front_end_state_freezing_to_replay_value_forbidden",
    ):
        if c.get(key) is not True:
            f.append(f"final numeric contract lost {key}")
    if c.get("stochastic_noise_realization_used_as_homogeneous_pruning") is not False:
        f.append("stochastic realization entered homogeneous pruning")
    if c.get("required_final_inequality") != "Omega_W - delta * P_W >= 0 on full H18/A21 coordinates":
        f.append("canonical full-word inequality changed")
    if float(c.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("useful gate changed")
    if int(c.get("H_dimension", 0)) != 18 or int(c.get("A_dimension", 0)) != 21:
        f.append("H/A dimensions changed")

    meas = d.get("measurement_runtime", {})
    if meas.get("accelerometer_std_mps2") != [0.2, 0.2, 0.2]:
        f.append("configured accelerometer std changed")
    if meas.get("magnetometer_std_uT") != [0.3, 0.3, 0.3]:
        f.append("configured magnetometer std changed")

    phys = d.get("physical_SEA3_scope", {})
    if phys.get("global_finite_window_realization_left_inclusion_closed") is not False:
        f.append("global SEA3 left inclusion falsely promoted")
    if phys.get("unqualified_RAO_coupling_used_as_hard_pruning") is not False:
        f.append("unqualified RAO pruning entered P3")
    if d.get("P3_promoted") is not False or d.get("P4_may_consume") is not False:
        f.append("precondition contract promoted downstream theorem")
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
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "architecture": d["canonical_architecture"],
        "word_horizon_s": d["word"]["horizon_s"],
        "startup_runtime": d["startup_runtime"],
        "all_preconditions": d["all_current_machine_checkable_preconditions_present"],
        "measurement_runtime": d["measurement_runtime"],
        "final_inequality": d["final_numeric_contract"]["required_final_inequality"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
