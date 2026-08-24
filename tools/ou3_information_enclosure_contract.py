#!/usr/bin/env python3
"""Build the continuous-source promotion contract for adaptive OU-III.

Executed replay values remain sanity anchors.  Deployment promotion requires
validated source-word endpoint bounds in the paper's group-compatible node
metrics; the Kalman inverse covariance may be used to construct/condition those
bounds but is not itself required to be the nonlinear Lyapunov metric.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_numerical_certificate as BASE
import ou3_source_domain_contract as SOURCE_DOMAIN

SCHEMA = 3


def mode_contract(info_mode: dict, completion_mode: dict, mode: str) -> dict:
    first = dict(info_mode.get("selected") or {})
    strongest = dict(info_mode.get("strongest_executed_margin") or first)
    horizon = strongest.get("horizon_s", first.get("horizon_s"))
    executed = {
        "first_strict_horizon_s": first.get("horizon_s"),
        "first_strict_lambda_worst": first.get("lambda_worst_information"),
        "recommended_horizon_s": horizon,
        "lambda_worst_information": strongest.get("lambda_worst_information"),
        "relative_Riccati_injection_margin_worst": strongest.get(
            "relative_Riccati_injection_margin_worst",
            strongest.get("omega_relative_lambda_min_worst"),
        ),
        "Sigma_endpoint_lambda_min": strongest.get("Sigma_endpoint_lambda_min"),
        "Sigma_endpoint_lambda_max": strongest.get("Sigma_endpoint_lambda_max"),
        "replay_asymptotic_floor_b_star": completion_mode.get(
            "asymptotic_floor_b_star_replay", completion_mode.get("invariant_level_b_replay")
        ),
        "replay_finite_capture_level_b_eta": completion_mode.get("finite_capture_level_b_eta_replay"),
    }
    return {
        "mode": mode,
        "linear_information_metric_role": (
            "Sigma_KF^-1 and equivalent congruences may certify Riccati/information bounds but are not automatically the nonlinear Lyapunov metric"
        ),
        "required_path_metric": "Pbar_i=blkdiag((a_R_i/2) I3, P_xi_i), a_R_i>0, P_xi_i>>0",
        "recommended_word_horizon_s": horizon,
        "executed_reference_only": executed,
        "required_continuous_bounds": {
            "source_complete": True,
            "outward_rounded": True,
            "pe_recurrence_window_s": (
                "finite positive deployment theorem hypothesis; every such normal-Live window contains a certified accepted two-packet vector-PE event"
            ),
            "translation_primary_route": "four spread-selected S updates observing complete [v,p,S,a_w] chain",
            "word_endpoint_relative_Riccati_injection_margin_lower": "> 0",
            "Sigma_lambda_min_lower": "> 0",
            "Sigma_lambda_max_upper": "finite",
            "prefix_information_gain_upper": "finite positive",
            "one_sample_decrease_required": False,
            "joint_source_reachability_required": True,
        },
        "required_nonlinear_bounds": {
            "theta_star": "0 < theta_star < pi",
            "endpoint_W_ratio_upper": "0 <= ratio < 1; verifier derives mu_W=1-ratio",
            "certified_level_W": "> 0",
            "all_word_prefixes_safe": True,
            "metric_lift": "W_i=a_R_i(1-cos(theta))+xi^T P_xi_i xi",
            "attitude_linear_cross_terms_in_Pbar": False,
        },
    }


def build_contract(info: dict, completion: dict) -> dict:
    modes = {
        "H": mode_contract(info.get("held", {}), completion.get("held", {}), "H"),
        "A": mode_contract(info.get("active", {}), completion.get("active", {}), "A"),
    }
    return {
        "schema": SCHEMA,
        "claim": "OU3_INFORMATION_METRIC_DEPLOYMENT_PROMOTION_CONTRACT",
        "upstream_executed_information_certificate": "PASS" if info.get("status") == "PASS" else "FAIL",
        "upstream_executed_replay_funnel": "PASS" if completion.get("status") == "PASS_EXECUTED_REPLAY" else "FAIL",
        "linear_identity": (
            "1-lambda_information=lambda_min(Sigma1^-1/2 Omega Sigma1^-1/2), Omega=Sigma1-Phi Sigma0 Phi^T"
        ),
        "metric_policy": {
            "P3_conditioning_coordinate_invariant": True,
            "P3_metric_need_not_equal_P4_metric": True,
            "P4_group_compatible_node_metric_required": True,
            "node_dependent_metrics_allowed": True,
            "common_Euclidean_or_common_quadratic_fallback_allowed": False,
        },
        "modes": modes,
        "source_domain_requirement": {
            "producer": "tools/ou3_source_domain_contract.py",
            "claim": "OU3_SOURCE_COMPLETE_IMPLEMENTATION_DOMAIN_CONTRACT",
            "validation": "final deployment gate regenerates the source-domain contract from the current implementation",
        },
        "source_word_language_requirement": {
            "producer": "tools/ou3_source_word_theorem_contract.py",
            "claim": "OU3_CONDITIONAL_SOURCE_COMPLETE_NORMAL_LIVE_WORD_LANGUAGE",
            "required_before_word_enclosure": [
                "finite positive PE recurrence window supplied as a deployment theorem hypothesis",
                "every recurrence window contains an accepted consecutive magnetic packet pair",
                "accepted accelerometer vectors are present at both vector times",
                "four spread-selected S updates cover the complete [v,p,S,a_w] translation chain",
                "all jointly source-reachable accepted/rejected branches between required PE events remain covered",
                "H and A same-mode words remain fixed-dimensional",
                "dimension-changing/reference/reset events remain separate hybrid obligations",
            ],
            "anti_shortcut": (
                "a single favorable vector pair, three-S integrator-only route without its independent a_w hypothesis, independently multiplied cell extrema, or a replay-observed acceptance pattern is not a source-complete word certificate"
            ),
        },
        "hybrid_requirements": {
            "required_kinds": list(SOURCE_DOMAIN.HYBRID_OBLIGATIONS),
            "primitive_bounds_per_jump": [
                "source_complete", "outward_rounded", "source_level_W_upper",
                "jump_gain_upper", "additive_W_upper", "destination_level_W", "destination_mode",
            ],
            "held_to_active_extra": [
                "source_dimension=18", "destination_dimension=21",
                "dimension_change_handled_by_embedding=true", "new_coordinate_W_upper",
            ],
            "periodic_aw_covariance_sync_extra": [
                "proof_mode=PSD_NONEXPANSIVE",
                "P_plus=P_minus+E_a Delta_plus E_a^T with Delta_plus>=0",
                "inverse-covariance information energy nonexpansive",
            ],
            "tilt_reset_rule": (
                "discarded pre-reset tilt energy is not charged in the multiplicative jump gain; only continuous coordinates plus the analytic reset additive term enter"
            ),
            "cooldown_rule": "use products over reachable cooldown word tilings, not powers of a global worst-word factor",
        },
        "stochastic_requirements": {
            "source_noise": "source_noise_certificate.json; standardized pre-gate Gaussian primitive covariance <= I",
            "validated_sensitivity_bounds": [
                "source_complete", "outward_rounded", "localization_prefix_safe",
                "localization_radius_standardized", "word_samples_upper", "finite_horizon_words",
                "funnel_level_a", "W0_upper", "L_X_upper", "G_bar_upper", "c_zw_upper",
                "r_star_upper", "c_ww_upper", "g_W_upper", "h_W_upper",
            ],
            "verifier_derives": [
                "s2", "s4", "nu1", "nu_W", "lambda_s", "sigma_s^2",
                "Gaussian localization t_star", "Freedman excursion probability",
                "total finite-horizon failure probability",
            ],
            "markov_union_fallback_allowed": False,
        },
        "validated_backend_requirements": [
            "validated arithmetic", "outward rounding",
            "source-generated enclosure, not trajectory fitting",
            "complete continuous source coverage relative to the declared theorem operating envelope",
            "joint source reachability across scheduled parameters and branch state",
            "finite recurring vector-PE hypothesis explicit rather than inferred from selected packets",
            "source-domain artifact exactly matches current implementation-generated contract",
        ],
        "promotion_rule": (
            "deployment PASS requires strict validated complete-word endpoint bounds, group-compatible node metrics, prefix safety, hybrid inward margins, Gaussian/Freedman stochastic concentration, finite capture, and exact current-source binding; no older scalar, common-Euclidean, or repeated-one-step fallback path may promote the theorem"
        ),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--certificate-dir", type=Path, default=BASE.DEFAULT_OUT)
    args = ap.parse_args()
    cert = args.certificate_dir.resolve()
    info = json.loads((cert / "information_certificate.json").read_text())
    completion = json.loads((cert / "information_completion.json").read_text())
    contract = build_contract(info, completion)
    path = cert / "information_enclosure_contract.json"
    path.write_text(json.dumps(contract, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "schema": contract["schema"],
        "linear": contract["upstream_executed_information_certificate"],
        "replay": contract["upstream_executed_replay_funnel"],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
