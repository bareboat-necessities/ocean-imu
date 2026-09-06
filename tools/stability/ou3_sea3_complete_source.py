#!/usr/bin/env python3
"""Complete SEA3 source family for canonical OU-III P3.

There is exactly one promotable source family: complete admitted SEA3 Normal
Live.  No point word, tuner rectangle, independent sea/RAO product, selected
four-S word, predecessor graph, Gaussian good-event source, or independent
tau/sigma/R_S/T_S box is an alternate source.

A canonical source word carries, simultaneously,

* the compact three-partition JONSWAP/PM directional sea state;
* exact partition-energy coupling H_s^2 = sum_r H_r^2;
* per-partition H_r/T_p,r peak-steepness admissibility;
* gamma, mean direction and directional spreading and the compact SEA3
  transition relation;
* a phase-continuous SEA3 oscillator/shaping state x^s and the same admitted
  translational and rotational vessel-response state;
* the resulting pathwise CoG acceleration/body-rate sequence satisfying the
  declared Normal-Live SEA3 conditions;
* the exact shipping measurement-only vertical, WavePeriodEstimator, band-pass,
  moment, tuner-EMA, staged-commit and pseudo-scheduler states;
* the resulting committed tau, sigma_aw and R_S; and
* every Normal-Live accelerometer/magnetometer/vector-PE event required by the
  theorem.

The Kalman coordinates tau, sigma_aw, R_S, T_S, F_k and Q_k are therefore
DERIVED coordinates of the SEA3 word.  They are never independently selected.
At every due S=0 update the full-state Joseph recursion uses the actual committed
SpectralMSE R_S from that same source word, including the deployed per-axis
standard-deviation factors [0.72, 0.72, 1].

The finite-horizon stochastic calculation retained below is diagnostic/corollary
material only.  It cannot generate, prune, or promote the homogeneous P3 source
family; configured Racc/Rmag and hard SEA3 source conditions stay in the word.

This file defines the conditional theorem source contract.  It does not claim
that the complete 3 s phase-continuous family has already been interval-
materialized.  P3 remains fail-closed until that SEA3 family itself executes the
literal H18/A21 word and the final full-matrix gate closes.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, down, up
import ou3_validated_transcendentals as VT
import ou3_sea3_physical_admissibility as PHYSICAL
import ou3_sea3_p1_compatibility as P1COMPAT
import ou3_sea3_directional_response_family as RESPONSE
import ou3_sea3_wave_period_spectral_identity as PERIOD_ID
import ou3_sea3_spectral_moment_bridge as MOMENT
import ou3_sea3_wave_period_frontend as FRONTEND
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_p3_pseudo_scheduler_progress_certificate as SCHED

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
DEFAULT_RESPONSE_DOMAIN = REPO / "tools" / "ou3_sea3_directional_response_domain.json"
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
SCHEMA = 3
QUALIFICATION = "OU3_SEA3_COMPLETE_NORMAL_LIVE_SOURCE_V3"
DIM = 3


def _exp_minus_integer(t: int) -> Interval:
    if not isinstance(t, int) or isinstance(t, bool) or t < 0:
        raise ValueError("tail exponent must be a nonnegative integer")
    base = VT.exp_point(-0.5)
    out = Interval.point(1.0)
    for _ in range(2 * t):
        out = out * base
    return out


def _select_tail_exponent(samples: int, budget: float) -> dict:
    """Non-promoting stochastic diagnostic used only for a later corollary."""
    if not isinstance(samples, int) or isinstance(samples, bool) or samples <= 0:
        raise ValueError("samples must be positive")
    blo = down(float(budget))
    for t in range(1, 257):
        e = _exp_minus_integer(t)
        prob = up(float(2 * DIM * samples) * e.hi)
        if prob <= blo:
            return {
                "integer_tail_exponent": t,
                "validated_exp_minus_t": e.as_list(),
                "failure_probability_upper": prob,
                "allocated_budget_lower": blo,
            }
    raise RuntimeError("finite-horizon SEA3 concentration exponent exceeded 256")


def _trace_threshold(cap: float, t: int) -> float:
    return down(down(float(cap) * float(cap)) / float(2 * DIM * t))


def _source_rs_parity() -> dict[str, bool]:
    text = WRAPPER.read_text(encoding="utf-8")
    return {
        "deployed_law_is_SpectralMSE": (
            "RSAdaptationLaw rs_law_ = RSAdaptationLaw::SpectralMSE;" in text
        ),
        "SpectralMSE_target_uses_realized_TS": (
            "const float TS = pseudo_update_period_for_(tau);" in text
            and "return rs_mse_coeff_ * rs_qeff_pow_" in text
            and "/ std::sqrt(TS);" in text
        ),
        "SpectralMSE_skips_extra_information_rate_scale": (
            "if (rs_law_ != RSAdaptationLaw::Cubic) return 1.0f;" in text
        ),
        "applied_RS_is_sent_to_filter": (
            "const float RSbase = std::min(std::max(tune_.RS_applied, min_R_S_), max_R_S_);" in text
            and "const float RSb = RSbase * pseudo_update_information_rate_scale_();" in text
            and "mekf_->set_RS_noise(Eigen::Vector3f(" in text
        ),
        "horizontal_RS_factors_are_0p72": (
            "float R_S_x_factor_ = 0.72f;" in text
            and "float R_S_y_factor_ = 0.72f;" in text
        ),
        "pseudo_period_is_committed_tau_function": (
            "const float requested = pseudo_update_tau_ratio_ * tau;" in text
            and "mekf_->set_pseudo_update_period_s(period);" in text
        ),
        "candidate_RS_has_own_EMA": "tune_.RS_applied    += alpha_RS" in text,
        "tau_sigma_share_common_EMA": (
            "tune_.tau_applied   += alpha" in text
            and "tune_.sigma_applied += alpha" in text
        ),
        "staged_commit_applies_ou_then_RS": (
            "void apply_pending_online_tune_()" in text
            and "apply_ou_tune_(false);" in text
            and "apply_RS_tune_();" in text
        ),
    }


def build(
    domain_path: Path = DEFAULT_DOMAIN,
    response_domain_path: Path = DEFAULT_RESPONSE_DOMAIN,
) -> dict:
    domain_path = Path(domain_path).resolve()
    response_domain_path = Path(response_domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("complete SEA3 source cannot be trajectory fitted")

    physical = PHYSICAL.build(domain_path)
    p1compat = P1COMPAT.build(domain_path, response_domain_path)
    response = RESPONSE.directional_response_enclosure(REPO, response_domain_path)
    period_id = PERIOD_ID.build()
    moment = MOMENT.build()
    frontend = FRONTEND.build(REPO)
    dynamic = DYNAMIC.build(domain_path)
    scheduler = SCHED.build(domain_path)

    prereq = {
        "physical": PHYSICAL.validate(physical),
        "p1_compatibility": P1COMPAT.validate(p1compat),
        "response": RESPONSE.validate(response),
        "period_identity": PERIOD_ID.validate(period_id),
        "spectral_moment": MOMENT.validate(moment),
        "frontend": FRONTEND.validate(frontend),
        "dynamic_source_parity": DYNAMIC.validate(dynamic),
        "scheduler": SCHED.validate(scheduler),
    }
    bad = {k: v for k, v in prereq.items() if v}
    if bad:
        raise RuntimeError(f"complete SEA3 prerequisite failure: {bad}")

    dt_cpp = float(dynamic["validated_rate_and_jump_bounds"]["dt_s"])
    samples = int(math.ceil(3.0 / dt_cpp))
    if samples * dt_cpp < 3.0:
        samples += 1

    # Retained only for a later finite-horizon stochastic forcing corollary.
    # It must not generate or prune the homogeneous SEA3 Riccati source family.
    total_budget = float(domain["stochastic"]["finite_horizon_failure_probability_budget"])
    per_event_budget = down(total_budget / 2.0)
    acc_tail = _select_tail_exponent(samples, per_event_budget)
    rate_tail = _select_tail_exponent(samples, per_event_budget)
    live = domain["normal_live"]
    acc_trace_threshold = _trace_threshold(
        float(live["non_gravitational_cog_acceleration_norm_upper_mps2"]),
        acc_tail["integer_tail_exponent"],
    )
    rate_trace_threshold = _trace_threshold(
        float(live["body_rate_norm_upper_deg_s"]),
        rate_tail["integer_tail_exponent"],
    )
    combined_failure = up(
        float(acc_tail["failure_probability_upper"])
        + float(rate_tail["failure_probability_upper"])
    )

    rs_parity = _source_rs_parity()
    rs_failures = [k for k, ok in rs_parity.items() if not ok]
    sea = moment["sea_family"]
    directional = json.loads(response_domain_path.read_text(encoding="utf-8"))[
        "directional_spectrum_contract"
    ]

    source_coordinates = {
        "surface_partitions": [
            "H_r", "T_p_r", "gamma_r", "mean_direction_r", "spreading_r"
        ],
        "surface_constraints": [
            "0 <= number_of_active_partitions <= 3",
            "H_s^2 = sum_r H_r^2",
            "2*pi*H_r/(g*T_p_r^2) <= S_p,max(T_p_r)",
            "gamma_r in declared JONSWAP interval",
            "directional density nonnegative and normalized",
            "lambda_{k+1} in compact R_lambda(lambda_k)",
        ],
        "phase_continuous_sea_state": [
            "x^s oscillator/shaping state",
            "lambda compact SEA3 parameter state",
        ],
        "translational_response": ["G", "f_c", "p", "complex h(f,theta)"],
        "rotational_response": ["K_rot", "f_c_rot", "q_rot", "complex r(f,theta)"],
        "finite_window_response_state": [
            "post-RAO CoG acceleration process from the same x^s/lambda word",
            "post-RAO body-rate process from the same x^s/lambda word",
        ],
        "front_end_state": [
            "vertical observer state",
            "WavePeriodEstimator two high-pass states",
            "WavePeriodEstimator velocity/elevation states",
            "WavePeriodEstimator moment/EWMA/log-period/usable-latch states",
            "adaptive wave-band state",
            "acceleration first/second moment state",
        ],
        "adaptive_state": [
            "tau_candidate", "sigma_candidate", "R_S_candidate",
            "tau_committed", "sigma_committed", "R_S_committed",
            "commit timer/pending flag", "pseudo scheduler elapsed",
        ],
        "vector_geometry_state": [
            "attitude/rotation geometry", "gravity direction",
            "magnetic reference and accepted magnetic packet geometry",
        ],
    }

    no_fallback = {
        "point_source_generator": False,
        "independent_tau_sigma_RS_TS_generator": False,
        "dynamic_tuner_rectangle_generator": False,
        "independent_sea_x_RAO_generator": False,
        "independent_partition_height_period_generator": False,
        "selected_four_S_generator": False,
        "one_step_process_strictness_generator": False,
        "source_history_graph_generator": False,
        "predecessor_path_generator": False,
        "arbitrary_P0_generator": False,
        "trajectory_replay_generator": False,
        "gaussian_good_event_source_generator": False,
        "spectral_moment_only_source_generator": False,
        "arbitrary_bounded_input_source_generator": False,
    }

    source_contract_ready = (
        not rs_failures
        and physical["SEA3_parameter_domain_compact"] is True
        and physical["compact_transition_relation_is_theorem_domain"] is True
        and bool(p1compat["coupled_SEA3_domain_required"])
        and all(frontend["source_parity"].values())
        and all(v is False for v in no_fallback.values())
    )

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_P3_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "source_coordinates": source_coordinates,
        "no_fallback_generators": no_fallback,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_filter_domain_shrunk": False,
        "theorem_conditional_on_admitted_complete_SEA3_word": True,
        "global_physical_deployment_left_inclusion_closed_here": False,
        "retired_P2_stack_consumed": False,
        "word_horizon_s": 3.0,
        "word_samples": samples,
        "SEA3_surface_family": {
            "modes_max": int(sea["m_max"]),
            "parameter_domain_compact": bool(physical["SEA3_parameter_domain_compact"]),
            "compact_transition_relation_is_theorem_domain": bool(
                physical["compact_transition_relation_is_theorem_domain"]
            ),
            "partition_spectrum": sea["partition_frequency_shape"],
            "gamma_interval": sea["declared_gamma_interval"],
            "directional_density": directional["directional_density"],
            "mean_direction_rad": directional["mean_direction_rad"],
            "spreading_parameter": directional["spreading_parameter"],
            "partition_energy_coordinate": sea["partition_energy_coordinate"],
            "total_energy_coupling": sea["total_energy_coupling"],
            "partition_peak_steepness_constraint": physical["three_partition_contract"][
                "active_partition_constraint"
            ],
            "total_Hs_upper_m": float(physical["repository_total_Hs_upper_m"]),
            "independent_H_T_extrema_forbidden": True,
            "independent_partition_height_maxima_forbidden": True,
        },
        "SEA3_dynamic_realization": {
            "phase_continuous": True,
            "sea_parameter_state": "lambda=(H_r,T_p_r,gamma_r,beta_r,s_r)_{r=1..3}",
            "shaping_state": "x^s",
            "compact_transition_relation": "lambda_{k+1} in R_lambda(lambda_k)",
            "augmented_source_state": "zeta=(x^s,lambda,z^t,q)",
            "same_realization_drives_translation_rotation_frontend_tuner_geometry": True,
            "hard_pathwise_acceleration_and_body_rate_conditions_retained": True,
            "finite_window_family_materialized": False,
            "probabilistic_event_may_substitute_for_realization": False,
            "arbitrary_bounded_input_may_substitute_for_realization": False,
        },
        "SEA3_translational_response_family": response,
        "SEA3_response_couplings": {
            "independent_sea_x_RAO_cartesian_product_forbidden": bool(
                p1compat["coupled_SEA3_domain_required"]
            ),
            "pathwise_non_gravitational_cog_acceleration_norm_upper_mps2": float(
                live["non_gravitational_cog_acceleration_norm_upper_mps2"]
            ),
            "pathwise_body_rate_norm_upper_deg_s": float(
                live["body_rate_norm_upper_deg_s"]
            ),
            "same_H_s_partition_energy_enters_translation_and_rotation_conditions": True,
            "only_same_phase_continuous_SEA3_realization_may_generate_P3_words": True,
            "moment_or_probability_bound_may_not_generate_P3_word": True,
        },
        "stochastic_forcing_corollary": {
            "role_in_P3": domain["stochastic"]["role_in_P3"],
            "used_to_generate_P3_source_words": False,
            "used_to_prune_homogeneous_P3_family": False,
            "configured_Racc_Rmag_remain_in_every_covariance_update": True,
            "centered_Gaussian_response_diagnostic": True,
            "samples": samples,
            "total_failure_budget": total_budget,
            "acceleration_tail": acc_tail,
            "body_rate_tail": rate_tail,
            "acceleration_trace_threshold_diagnostic_m2_s4": acc_trace_threshold,
            "body_rate_trace_threshold_diagnostic_deg2_s2": rate_trace_threshold,
            "combined_failure_probability_upper": combined_failure,
            "combined_within_budget": combined_failure <= down(total_budget),
        },
        "SEA3_period_and_frontend": {
            "steady_response_weighted_period_identity": period_id[
                "continuous_time_steady_state_identity"
            ],
            "surface_multimodal_identity": moment["analytical_lemmas"][
                "multimodal_zero_crossing_identity"
            ],
            "surface_Tz_substituted_for_tuner_Tz": False,
            "exact_discrete_frontend_source_parity_pass": all(
                frontend["source_parity"].values()
            ),
            "finite_frontend_state_is_part_of_every_source_word": True,
            "frontend_state_may_be_frozen_to_replay_value": False,
        },
        "derived_adaptive_source": {
            "primitive_independent_tau_sigma_RS_TS": False,
            "same_SEA3_frontend_path_generates_tau_sigma_RS_targets": True,
            "same_candidate_snapshot_commits_tau_sigma_RS": True,
            "T_S_is_function_of_same_committed_tau": True,
            "Q_uses_same_committed_tau_sigma": True,
            "source_recurrence_rate_and_commit_bounds": dynamic[
                "validated_rate_and_jump_bounds"
            ],
            "rate_bounds_are_constraints_on_SEA3_derived_path_not_a_word_generator": True,
        },
        "R_S_regularizer": {
            "source_parity": rs_parity,
            "source_parity_failures": rs_failures,
            "deployed_law": "SpectralMSE",
            "actual_applied_R_S_required_at_every_due_S_update": True,
            "axis_std_factors": [0.72, 0.72, 1.0],
            "extra_information_rate_rescale": 1.0,
            "pseudo_scheduler_recurrence_certificate": bool(
                scheduler["scheduler_recurrence_certificate"]
            ),
            "certified_uniform_max_gap_s": scheduler["certified_uniform_max_gap_s"],
            "all_due_S_updates_remain_in_full_word": True,
            "full_P_column_S_cross_covariance_action_required": True,
            "R_S_may_not_be_replaced_by_process_strictness": True,
            "selected_four_S_subset_may_not_replace_full_scheduler_word": True,
        },
        "Normal_Live_nonsea_conditions": {
            "all_valid_accelerometer_updates_required": bool(
                live["accelerometer_update_required_each_valid_imu_sample_after_live_entry"]
            ),
            "accelerometer_rejection_in_scope": bool(
                live["accelerometer_rejection_in_normal_live_scope"]
            ),
            "magnetic_norm_uT": [
                float(live["magnetic_vector_norm_lower_uT"]),
                float(live["magnetic_vector_norm_upper_uT"]),
            ],
            "vector_sine_separation_lower": float(live["vector_sine_separation_lower"]),
            "vector_PE_recurrence_window_s": float(live["vector_pe_recurrence_window_s"]),
            "hard_attitude_rewrite_inside_word": bool(
                live["hard_attitude_rewrite_inside_word"]
            ),
        },
        "P3_source_contract_ready": source_contract_ready,
        "P3_source_family_materialized": False,
        "P3_promoted": False,
        "next_obligation": (
            "materialize the phase-continuous compact SEA3 3 s family itself, propagate the same "
            "x^s/lambda/front-end/tuner/scheduler/vector state through every shipping event, and "
            "execute the full H18/A21 Joseph/Riccati word"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_P3_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical P3 source is not complete SEA3")
    for key in (
        "theorem_conditional_on_admitted_complete_SEA3_word",
        "P3_source_contract_ready",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_filter_domain_shrunk",
        "global_physical_deployment_left_inclusion_closed_here", "retired_P2_stack_consumed",
        "P3_source_family_materialized", "P3_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    fallback = d.get("no_fallback_generators", {})
    if not fallback or any(v is not False for v in fallback.values()):
        f.append("a fallback source generator is still enabled")
    if int(d.get("word_samples", 0)) < 600:
        f.append("complete SEA3 word does not cover 3 s at 200 Hz")
    sea = d.get("SEA3_surface_family", {})
    if sea.get("modes_max") != 3 or sea.get("gamma_interval") != [1.0, 7.0]:
        f.append("SEA3 surface family changed")
    for key in (
        "parameter_domain_compact",
        "compact_transition_relation_is_theorem_domain",
        "independent_H_T_extrema_forbidden",
        "independent_partition_height_maxima_forbidden",
    ):
        if sea.get(key) is not True:
            f.append(f"SEA3 surface family lost {key}")
    realization = d.get("SEA3_dynamic_realization", {})
    for key in (
        "phase_continuous",
        "same_realization_drives_translation_rotation_frontend_tuner_geometry",
        "hard_pathwise_acceleration_and_body_rate_conditions_retained",
    ):
        if realization.get(key) is not True:
            f.append(f"SEA3 dynamic realization lost {key}")
    for key in (
        "finite_window_family_materialized",
        "probabilistic_event_may_substitute_for_realization",
        "arbitrary_bounded_input_may_substitute_for_realization",
    ):
        if realization.get(key) is not False:
            f.append(f"SEA3 dynamic realization open/forbidden flag {key} changed")
    coupled = d.get("SEA3_response_couplings", {})
    for key in (
        "independent_sea_x_RAO_cartesian_product_forbidden",
        "same_H_s_partition_energy_enters_translation_and_rotation_conditions",
        "only_same_phase_continuous_SEA3_realization_may_generate_P3_words",
        "moment_or_probability_bound_may_not_generate_P3_word",
    ):
        if coupled.get(key) is not True:
            f.append(f"SEA3 response coupling lost {key}")
    stochastic = d.get("stochastic_forcing_corollary", {})
    if stochastic.get("used_to_generate_P3_source_words") is not False:
        f.append("stochastic good event re-entered as P3 source generator")
    if stochastic.get("used_to_prune_homogeneous_P3_family") is not False:
        f.append("stochastic good event re-entered homogeneous P3 pruning")
    if stochastic.get("configured_Racc_Rmag_remain_in_every_covariance_update") is not True:
        f.append("stochastic corollary changed configured measurement covariance")
    adapt = d.get("derived_adaptive_source", {})
    for key in (
        "same_SEA3_frontend_path_generates_tau_sigma_RS_targets",
        "same_candidate_snapshot_commits_tau_sigma_RS",
        "T_S_is_function_of_same_committed_tau",
        "Q_uses_same_committed_tau_sigma",
        "rate_bounds_are_constraints_on_SEA3_derived_path_not_a_word_generator",
    ):
        if adapt.get(key) is not True:
            f.append(f"derived SEA3 adaptive path lost {key}")
    if adapt.get("primitive_independent_tau_sigma_RS_TS") is not False:
        f.append("tuner coordinates became primitive perturbations")
    rs = d.get("R_S_regularizer", {})
    if rs.get("source_parity_failures"):
        f.extend(f"R_S source parity failed: {x}" for x in rs["source_parity_failures"])
    for key in (
        "actual_applied_R_S_required_at_every_due_S_update",
        "pseudo_scheduler_recurrence_certificate",
        "all_due_S_updates_remain_in_full_word",
        "full_P_column_S_cross_covariance_action_required",
        "R_S_may_not_be_replaced_by_process_strictness",
        "selected_four_S_subset_may_not_replace_full_scheduler_word",
    ):
        if rs.get(key) is not True:
            f.append(f"R_S regularizer lost {key}")
    if rs.get("axis_std_factors") != [0.72, 0.72, 1.0]:
        f.append("R_S axis factors changed")
    if float(rs.get("extra_information_rate_rescale", math.nan)) != 1.0:
        f.append("SpectralMSE R_S incorrectly received cadence rescale")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--response-domain", type=Path, default=DEFAULT_RESPONSE_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain, args.response_domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "source": d["canonical_P3_source"],
        "samples": d["word_samples"],
        "SEA3_dynamic_realization": d["SEA3_dynamic_realization"],
        "no_fallback_generators": d["no_fallback_generators"],
        "response_couplings": d["SEA3_response_couplings"],
        "stochastic_forcing_corollary": d["stochastic_forcing_corollary"],
        "R_S": d["R_S_regularizer"],
        "P3_source_contract_ready": d["P3_source_contract_ready"],
        "P3_source_family_materialized": d["P3_source_family_materialized"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
