#!/usr/bin/env python3
"""Build the continuous-source promotion contract for adaptive OU-III.

The executed certificate in ``ou3_information_certificate.py`` uses the actual
Kalman covariance as a source-varying information metric,

    M(g) = Sigma_KF(g)^(-1).

For every deterministic source word,

    Sigma_1 = Phi Sigma_0 Phi^T + Omega,

and therefore

    Phi^T Sigma_1^-1 Phi <= (1-eta) Sigma_0^-1

whenever

    Sigma_1^-1/2 Omega Sigma_1^-1/2 >= eta I,  eta > 0.

This tool does not invent a validated lower bound eta.  It turns the successful
eight-replay result into a deterministic, machine-readable contract specifying
what a continuous-source interval/Taylor-model backend must prove.  Executed
values are retained only as sanity anchors; they are never promoted by this
stage.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_numerical_certificate as BASE

SCHEMA = 2


def _finite_positive(value) -> bool:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(x) and x > 0.0


def mode_contract(info_mode: dict, completion_mode: dict, mode: str) -> dict:
    """Return one deterministic continuous-source obligation set.

    Prefer the tested horizon with the strongest executed information margin.
    A longer word is allowed by the source-word theorem and gives the validated
    backend more room to prove a robust positive Riccati injection.  The first
    strict horizon is kept as a diagnostic anchor, not as a requirement.
    """
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
            "mu_W_lower": "> 0",
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
        "hybrid_requirements": [
            "startup_handoff_into_certified_destination_sublevel",
            "held_to_active_jump_into_certified_destination_sublevel",
            "magnetic_regauge_jump_into_certified_destination_sublevel",
            "tilt_reset_jump_into_certified_destination_sublevel",
            "cooldown_reentry_into_certified_destination_sublevel",
        ],
        "stochastic_requirements": [
            "source_uniform_Sigma_bar_norm_upper",
            "source_uniform_b_W_upper",
            "source_uniform_v_W_upper",
            "finite_horizon_failure_probability_upper",
        ],
        "validated_backend_requirements": [
            "validated arithmetic",
            "outward rounding",
            "source-generated enclosure, not trajectory fitting",
            "complete continuous source coverage",
        ],
        "promotion_rule": (
            "executed replay values are sanity anchors only; deployment PASS requires "
            "strict validated continuous-source linear, nonlinear, hybrid and stochastic bounds"
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
    path.write_text(json.dumps(contract, indent=2, sort_keys=True))
    print(json.dumps({
        "schema": contract["schema"],
        "linear": contract["upstream_executed_information_certificate"],
        "replay": contract["upstream_executed_replay_funnel"],
        "H_horizon_s": contract["modes"]["H"]["recommended_word_horizon_s"],
        "A_horizon_s": contract["modes"]["A"]["recommended_word_horizon_s"],
    }, indent=2, sort_keys=True))
    # An upstream mathematical FAIL is scientific output, not a contract-tool crash.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
