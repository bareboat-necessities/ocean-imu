#!/usr/bin/env python3
"""Build the continuous-source promotion contract for adaptive OU-III.

Schema 4 binds P3 and P4 to one source-correlated information geometry.  P4
uses the exact Cayley coordinate c(R)=2 tan(theta/2)u and

    W_g=s_m [c(R);xi]^T Sigma_KF(g)^-1[c(R);xi],

with one positive scalar s_m shared by every source node of a fixed-dimensional
mode.  The scalar only normalizes reported W levels; it changes neither physical
level sets nor generalized contraction ratios.  Full attitude--linear cross
terms remain in the metric.  The retired block-diagonal a_R/P_xi metric is not
a promotion fallback.

Executed replay values are diagnostic/falsification anchors only; they do not
choose the source-word horizon or the P4 metric.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_numerical_certificate as BASE
import ou3_source_domain_contract as SOURCE_DOMAIN

SCHEMA = 4


def mode_contract(info_mode: dict, completion_mode: dict, mode: str) -> dict:
    first=dict(info_mode.get("selected") or {})
    strongest=dict(info_mode.get("strongest_executed_margin") or first)
    executed={
        "first_strict_horizon_s":first.get("horizon_s"),
        "first_strict_lambda_worst":first.get("lambda_worst_information"),
        "strongest_executed_horizon_s":strongest.get("horizon_s"),
        "lambda_worst_information":strongest.get("lambda_worst_information"),
        "relative_Riccati_injection_margin_worst":strongest.get("relative_Riccati_injection_margin_worst",strongest.get("omega_relative_lambda_min_worst")),
        "Sigma_endpoint_lambda_min":strongest.get("Sigma_endpoint_lambda_min"),
        "Sigma_endpoint_lambda_max":strongest.get("Sigma_endpoint_lambda_max"),
        "replay_asymptotic_floor_b_star":completion_mode.get("asymptotic_floor_b_star_replay",completion_mode.get("invariant_level_b_replay")),
        "replay_finite_capture_level_b_eta":completion_mode.get("finite_capture_level_b_eta_replay"),
        "qualification":"DIAGNOSTIC_ONLY_NOT_A_THEOREM_HORIZON_OR_METRIC_SELECTOR",
    }
    return {
        "mode":mode,
        "required_path_metric":"CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC",
        "required_group_coordinate":"c(R)=2*tan(theta/2)*u=4*e_R/(1+tr(R))",
        "required_group_energy":"W_g=s_mode*[c(R);xi]^T Sigma_KF(g)^-1[c(R);xi]",
        "mode_global_positive_scale_required":True,
        "same_scale_on_every_source_node_in_mode":True,
        "endpoint_metric_source_correlation_required":True,
        "full_attitude_linear_cross_terms_retained":True,
        "executed_reference_only":executed,
        "required_continuous_bounds":{
            "source_complete":True,"outward_rounded":True,
            "pe_recurrence_window_s":"finite positive deployment theorem hypothesis with a certified accepted two-packet vector-PE event in every recurrence window",
            "translation_primary_route":"four spread-selected S updates observing complete [v,p,S,a_w] chain",
            "word_endpoint_relative_Riccati_injection_margin_lower":"> 0",
            "Sigma_lambda_min_lower":"> 0","Sigma_lambda_max_upper":"finite",
            "prefix_information_gain_upper":"finite positive",
            "one_sample_decrease_required":False,"joint_source_reachability_required":True,
        },
        "required_nonlinear_bounds":{
            "theta_star":"0 < theta_star < pi","certified_level_W":"> 0",
            "endpoint_relative_W_decrease_lower":"> 0","mu_W_lower":"> 0",
            "all_word_prefixes_safe":True,"accepted_correction_uses_source_series_branch":True,
            "metric_lift":"mode-global positive scalar times exact Cayley lift of matching source Kalman information metric",
            "block_diagonal_metric_fallback":False,
        },
    }


def build_contract(info: dict, completion: dict) -> dict:
    modes={"H":mode_contract(info.get("held",{}),completion.get("held",{}),"H"),
           "A":mode_contract(info.get("active",{}),completion.get("active",{}),"A")}
    return {
        "schema":SCHEMA,
        "claim":"OU3_INFORMATION_METRIC_DEPLOYMENT_PROMOTION_CONTRACT",
        "upstream_executed_information_certificate":"PASS" if info.get("status")=="PASS" else "FAIL",
        "upstream_executed_replay_funnel":"PASS" if completion.get("status")=="PASS_EXECUTED_REPLAY" else "FAIL",
        "linear_identity":"1-lambda_information=lambda_min(Sigma1^-1/2 Omega Sigma1^-1/2), Omega=Sigma1-Phi Sigma0 Phi^T",
        "metric_policy":{
            "P3_conditioning_coordinate_invariant":True,
            "P4_local_quadratic_is_positive_scalar_multiple_of_P3_information_metric":True,
            "P4_mode_global_scale_same_on_every_source_node":True,
            "P4_mode_global_scale_changes_no_physical_level_set_or_contraction_ratio":True,
            "P4_exact_Cayley_group_lift_required":True,
            "P4_endpoint_metric_matches_source_covariance":True,
            "P4_full_attitude_linear_cross_terms_retained":True,
            "node_dependent_metrics_required":True,
            "block_diagonal_group_metric_fallback_allowed":False,
            "common_Euclidean_or_common_quadratic_fallback_allowed":False,
        },
        "modes":modes,
        "source_domain_requirement":{"producer":"tools/ou3_source_domain_contract.py","claim":"OU3_SOURCE_COMPLETE_IMPLEMENTATION_DOMAIN_CONTRACT","validation":"regenerate from current implementation and declared theorem operating domain"},
        "source_word_language_requirement":{
            "producer":"tools/ou3_source_word_theorem_contract.py",
            "claim":"OU3_CONDITIONAL_SOURCE_COMPLETE_NORMAL_LIVE_WORD_LANGUAGE",
            "required_before_word_enclosure":[
                "finite positive PE recurrence window supplied as a deployment theorem hypothesis",
                "every recurrence window contains an accepted consecutive magnetic packet pair",
                "accepted accelerometer vectors are present at both vector times",
                "four spread-selected S updates cover the complete [v,p,S,a_w] translation chain",
                "all jointly source-reachable accepted/rejected branches between required PE events remain covered",
                "H and A same-mode words remain fixed-dimensional",
                "dimension-changing/reference/reset events remain separate hybrid obligations"],
            "anti_shortcut":"single favorable vector pair, three-S-only promotion, independent cell-extrema multiplication, replay-observed branch pattern, or repeated one-sample decrease cannot promote the theorem"},
        "hybrid_requirements":{
            "required_kinds":list(SOURCE_DOMAIN.HYBRID_OBLIGATIONS),
            "primitive_bounds_per_jump":["source_complete","outward_rounded","source_level_W_upper","jump_gain_upper","additive_W_upper","destination_level_W","destination_mode"],
            "held_to_active_extra":["source_dimension=18","destination_dimension=21","dimension_change_handled_by_embedding=true","new_coordinate_W_upper"],
            "periodic_aw_covariance_sync_extra":["proof_mode=PSD_NONEXPANSIVE","P_plus=P_minus+E_a Delta_plus E_a^T with Delta_plus>=0","inverse-covariance information energy nonexpansive"],
            "tilt_reset_rule":"discarded pre-reset tilt energy is not charged in multiplicative jump gain",
            "cooldown_rule":"use products over reachable cooldown word tilings, not powers of a global worst-word factor"},
        "stochastic_requirements":{
            "source_noise":"source_noise_certificate.json; standardized pre-gate Gaussian primitive covariance <= I",
            "validated_sensitivity_bounds":["source_complete","outward_rounded","localization_prefix_safe","localization_radius_standardized","word_samples_upper","finite_horizon_words","funnel_level_a","W0_upper","L_X_upper","G_bar_upper","c_zw_upper","r_star_upper","c_ww_upper","g_W_upper","h_W_upper"],
            "verifier_derives":["s2","s4","nu1","nu_W","lambda_s","sigma_s^2","Gaussian localization t_star","Freedman excursion probability","total finite-horizon failure probability"],
            "markov_union_fallback_allowed":False},
        "validated_backend_requirements":["validated arithmetic","outward rounding","source-generated enclosure, not trajectory fitting","complete continuous source coverage","joint source reachability","finite recurring vector-PE hypothesis explicit","exact deployed quaternion injection and source operation order","current-source binding"],
        "promotion_rule":"PASS requires complete-word P3 endpoint bounds, normalized exact Cayley lift of matching source information geometry, positive direct P4 gap and mu_W, prefix safety, hybrid inward margins, Gaussian/Freedman concentration, finite capture, and exact current-source binding; no older fallback route may promote the theorem",
    }


def main() -> int:
    ap=argparse.ArgumentParser(); ap.add_argument("--certificate-dir",type=Path,default=BASE.DEFAULT_OUT); args=ap.parse_args(); cert=args.certificate_dir.resolve()
    info=json.loads((cert/"information_certificate.json").read_text()); completion=json.loads((cert/"information_completion.json").read_text()); contract=build_contract(info,completion)
    (cert/"information_enclosure_contract.json").write_text(json.dumps(contract,indent=2,sort_keys=True),encoding="utf-8")
    print(json.dumps({"schema":contract["schema"],"linear":contract["upstream_executed_information_certificate"],"replay":contract["upstream_executed_replay_funnel"],"P4_metric":contract["metric_policy"]},indent=2,sort_keys=True)); return 0


if __name__=="__main__": raise SystemExit(main())
