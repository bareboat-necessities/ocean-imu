#!/usr/bin/env python3
"""BIAS0/1/2 premises for the unchanged complete-BRMM proof.

These are conditional model/graph obligations, not device qualification or a
new P4 certificate. True residual bias and corrected estimation error are
different coordinates. A stationary OU covariance never supplies a hard cap.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval
import ou3_validated_transcendentals as VT
import ou3_full_process_ucc as PROCESS

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools/stability/ou3_proof_operating_domain.json"
QUALIFICATION = "OU3_MEMS_BIAS_PHYSICAL_AND_WORD_PRECONDITIONS_V1"
DEPLOYMENT_BLOCKER = "BIAS0 assembled-sensor qualification and true residual-bias envelope remain open"


def homogeneous_bias_at(
    root: Sequence[Interval], tau_s: Interval, elapsed_s: Interval,
) -> tuple[Interval, Interval, Interval]:
    """Marginal enclosure derived from ONE retained root/tau, never a new root.

    Callers must retain the root and factor dependency in their word graph.
    Evaluating these marginals does not materialize a joint BRMM source cover.
    """
    if len(root) != 3 or any(not isinstance(x, Interval) for x in root):
        raise TypeError("true residual-bias root must have three outward intervals")
    if not isinstance(tau_s, Interval) or not isinstance(elapsed_s, Interval):
        raise TypeError("bias time constant and elapsed time must be intervals")
    if not all(math.isfinite(x) for x in (tau_s.lo, tau_s.hi, elapsed_s.lo, elapsed_s.hi)):
        raise ValueError("bias time coordinates must be finite")
    if tau_s.lo <= 0.0 or elapsed_s.lo < 0.0:
        raise ValueError("positive finite tau and nonnegative elapsed time required")
    if any(not math.isfinite(v) for x in root for v in (x.lo, x.hi)):
        raise ValueError("bias root must be finite")
    if elapsed_s.lo == elapsed_s.hi == 0.0:
        return tuple(root)
    phi = VT.exp_interval(-(elapsed_s / tau_s))
    return tuple(phi * x for x in root)


def ou_band_fraction(tau_s: float, omega_lo: float, omega_hi: float) -> float:
    """Two-sided +/- wave-band variance fraction; numerical corollary only."""
    if not all(math.isfinite(x) for x in (tau_s, omega_lo, omega_hi)):
        raise ValueError("finite OU band parameters required")
    if tau_s <= 0.0 or not 0.0 <= omega_lo <= omega_hi:
        raise ValueError("positive tau and ordered nonnegative angular frequencies required")
    return 2.0 / math.pi * (
        math.atan(omega_hi * tau_s) - math.atan(omega_lo * tau_s)
    )


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain = json.loads(Path(domain_path).read_text(encoding="utf-8"))
    live = domain["normal_live"]
    c = PROCESS._constants()
    return {
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "BIAS0": {
            "physical_model": "b_total=b0+K_T*dT+K_strain*strain+beta+d_nonGM",
            "shipping_centered_truth": "b_true=beta+d_det; retain all uncompensated deterministic terms",
            "assembled_sensor_qualification_closed": False,
            "qualified_tau_interval_s": None,
            "qualified_stationary_sigma_mps2": None,
            "qualified_true_residual_root_bound_mps2": None,
            "qualified_deterministic_mismatch_bound_mps2": None,
            "filter_tau_s": c["accel_bias_tau_s"],
            "filter_Q_ba_density": c["accel_bias_process_variance_density"],
            "filter_projection_radius_mps2": live["active_accelerometer_bias_projection_limit_mps2"],
            "filter_prior_is_physical_qualification": False,
            "estimated_state_cap_is_true_bias_cap": False,
        },
        "BIAS1": {
            "homogeneous_true_bias_graph": "beta_j=exp(-(t_j-t_0)/tau_true)*beta_0",
            "one_root_and_parameter_history_required": True,
            "true_bias_frozen_when_estimator_held": False,
            "H_to_A_reseeds_physical_bias": False,
            "matched_tau_required_for_homogeneous_prediction": True,
            "tau_mismatch_is_explicit_forcing": True,
            "measurement_bias_corrections_retained": True,
            "projection_uses_shipping_centered_truth": True,
            "bias_error_is_free_GM_history": False,
            "Kalman_corrections_charged_to_exogenous_ISS": False,
            "process_Q_removed_from_Riccati_word": False,
            "source_uniform_projection_compatibility_closed": False,
        },
        "BIAS2": {
            "sector": "Pi_sep=C_y^T W C_y-mu_sep X on the exact full-state same-history graph",
            "uniform_separation_constant_lower": None,
            "source_uniform_separation_closed": False,
            "requires_full_corrected_bias_error_history": True,
            "positive_separation_follows_from_BRMM_or_GM_alone": False,
            "each_consumed_sector_must_be_valid_on_its_prefix": True,
            "prefix_contraction_required": False,
        },
        "stochastic_corollary": {
            "PSD": "2*sigma_b^2*tau_b/(1+(omega*tau_b)^2)",
            "band_fraction": "2/pi*(atan(omega_hi*tau_b)-atan(omega_lo*tau_b))",
            "stationary_model_required": True,
            "Gaussian_forcing_has_hard_pathwise_cap": False,
            "PSD_used_to_prune_homogeneous_BRMM": False,
        },
        "P3_delta": 1.0e-18,
        "P3_numerical_matrices_changed": False,
        "P4_bias_preconditions_closed": False,
        "P5_capture_must_retain_bias_source_and_forcing": True,
    }


def validate(d: dict) -> list[str]:
    """Do not let descriptive model premises masquerade as numerical closure."""
    expected = build()
    failures = []
    for key, value in expected.items():
        if isinstance(value, dict):
            actual = d.get(key, {})
            if not isinstance(actual, dict):
                failures.append(f"{key} contract missing")
                continue
            for subkey, required in value.items():
                # Projection radius is domain-specific and checked against the
                # existing positive deployed radius by the source consumers.
                if subkey == "filter_projection_radius_mps2":
                    v = actual.get(subkey)
                    if not isinstance(v, (int, float)) or not math.isfinite(v) or v <= 0:
                        failures.append("BIAS0 invalid projection radius")
                elif actual.get(subkey) != required or subkey not in actual:
                    failures.append(f"{key}.{subkey} changed or falsely promoted")
        elif d.get(key) != value:
            failures.append(f"{key} changed or missing")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    d.update(validation_pass=not failures, validation_failures=failures)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(d, sort_keys=True))
    return int(bool(failures))


if __name__ == "__main__":
    raise SystemExit(main())
