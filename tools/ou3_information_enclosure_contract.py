#!/usr/bin/env python3
"""Build the continuous-source promotion contract for adaptive OU-III.

The executed information certificate supplies sanity anchors only. This tool
states the primitive outward-rounded quantities a source-complete validated
backend must prove. Final nonlinear, hybrid and stochastic margins are
recomputed by ``ou3_validate_enclosure.py`` and ``ou3_deployment_gate.py`` and
are never accepted as asserted PASS values.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_numerical_certificate as BASE
import ou3_source_domain_contract as SOURCE_DOMAIN

SCHEMA = 2


def mode_contract(info_mode: dict, completion_mode: dict, mode: str) -> dict:
    first = dict(info_mode.get("selected") or {})
    strongest = dict(info_mode.get("strongest_executed_margin") or first)
    horizon = strongest.get("horizon_s")
    if horizon is None:
        horizon = first.get("horizon_s")

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
        "Sigma_endpoint_condition_bound": strongest.get("Sigma_endpoint_condition_bound"),
        "information_identity_residual_max": strongest.get("information_identity_residual_max"),
        "replay_asymptotic_floor_b_star": completion_mode.get(
            "asymptotic_floor_b_star_replay",
            completion_mode.get("invariant_level_b_replay"),
        ),
        "replay_finite_capture_level_b_eta": completion_mode.get(
            "finite_capture_level_b_eta_replay"
        ),
    }

    return {
        "mode": mode,
        "metric": "M(g)=Sigma_KF(g)^(-1)",
        "recommended_word_horizon_s": horizon,
        "executed_reference_only": executed,
        "required_continuous_bounds": {
            "source_complete": True,
            "outward_rounded": True,
            "relative_Riccati_injection_margin_lower": "> 0",
            "Sigma_lambda_min_lower": "> 0",
            "Sigma_lambda_max_upper": "finite",
            "prefix_information_gain_upper": "finite positive",
        },
        "required_nonlinear_bounds": {
            "theta_star": "0 < theta_star < pi",
            "endpoint_W_ratio_upper": "0 <= ratio < 1; verifier derives mu_W=1-ratio",
            "certified_level_W": "> 0",
            "all_word_prefixes_safe": True,
            "metric_lift": "zeta^T Sigma_KF(g)^(-1) zeta with zeta=[Log(R_e); xi]",
        },
    }


def build_contract(info: dict, completion: dict) -> dict:
    info_pass = info.get("status") == "PASS"
    replay_pass = completion.get("status") == "PASS_EXECUTED_REPLAY"
    modes = {
        "H": mode_contract(info.get("held", {}), completion.get("held", {}), "H"),
        "A": mode_contract(info.get("active", {}), completion.get("active", {}), "A"),
    }
    return {
        "schema": SCHEMA,
        "claim": "OU3_INFORMATION_METRIC_DEPLOYMENT_PROMOTION_CONTRACT",
        "upstream_executed_information_certificate": "PASS" if info_pass else "FAIL",
        "upstream_executed_replay_funnel": "PASS" if replay_pass else "FAIL",
        "metric": "source-varying inverse estimator covariance",
        "linear_identity": (
            "1-lambda_information=lambda_min(Sigma1^-1/2 Omega Sigma1^-1/2), "
            "Omega=Sigma1-Phi Sigma0 Phi^T"
        ),
        "modes": modes,
        "source_domain_requirement": {
            "producer": "tools/ou3_source_domain_contract.py",
            "claim": "OU3_SOURCE_COMPLETE_IMPLEMENTATION_DOMAIN_CONTRACT",
            "validation": (
                "final deployment gate regenerates the source-domain contract from the current "
                "implementation and rejects stale, truncated, or self-asserted domain artifacts"
            ),
        },
        "hybrid_requirements": {
            "required_kinds": list(SOURCE_DOMAIN.HYBRID_OBLIGATIONS),
            "primitive_bounds_per_jump": [
                "source_complete", "outward_rounded", "source_level_W_upper",
                "jump_gain_upper", "additive_W_upper", "destination_level_W",
                "destination_mode",
            ],
            "held_to_active_extra": [
                "source_dimension=18", "destination_dimension=21",
                "dimension_change_handled_by_embedding=true",
                "new_coordinate_W_upper",
            ],
            "periodic_aw_covariance_sync_extra": [
                "proof_mode=PSD_NONEXPANSIVE",
                "P_plus=P_minus+E_a Delta_plus E_a^T with Delta_plus>=0",
                "inverse-covariance information energy nonexpansive",
            ],
            "verifier_formula": (
                "post_W_upper=jump_gain_upper*source_level_W_upper+"
                "additive_W_upper+new_coordinate_W_upper; "
                "inward_margin=destination_level_W-post_W_upper>0"
            ),
        },
        "stochastic_requirements": {
            "source_noise": "source_noise_certificate.json; standardized pre-gate Gaussian primitive covariance <= I",
            "validated_sensitivity_bounds": [
                "source_complete", "outward_rounded", "localization_prefix_safe",
                "localization_radius_standardized", "word_samples_upper",
                "finite_horizon_words", "funnel_level_a", "W0_upper",
                "L_X_upper", "G_bar_upper", "c_zw_upper", "r_star_upper",
                "c_ww_upper", "g_W_upper", "h_W_upper",
            ],
            "verifier_derives": [
                "s2", "s4", "nu1", "nu_W", "lambda_s", "sigma_s^2",
                "b_W=a", "v_W=a^2/4", "Gaussian localization t_star",
                "Freedman excursion probability", "total finite-horizon failure probability",
            ],
        },
        "validated_backend_requirements": [
            "validated arithmetic",
            "outward rounding",
            "source-generated enclosure, not trajectory fitting",
            "complete continuous source coverage",
            "source-domain artifact exactly matches current implementation-generated contract",
        ],
        "promotion_rule": (
            "executed replay values are sanity anchors only; deployment PASS requires "
            "strict validated continuous-source linear and nonlinear bounds, recomputed "
            "hybrid inward margins, recomputed stochastic concentration, independently "
            "recomputed finite capture, and exact binding to the current implementation domain"
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
        "H_horizon_s": contract["modes"]["H"]["recommended_word_horizon_s"],
        "A_horizon_s": contract["modes"]["A"]["recommended_word_horizon_s"],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
