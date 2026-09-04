#!/usr/bin/env python3
"""Validated SEA0 subcertificate for the WavePeriodEstimator front-end map.

The physical SEA3 spectrum is parameterized by peak period, while shipping
WavePeriodEstimator estimates a zero-crossing/moment period after two shared
high-pass stages and two leaky integrations.  This module retains two pieces of
the source-to-tuner bridge in one producer:

* the validated steady single-frequency transfer identity of the implemented
  front end; and
* the exact startup/source-causality split between the fixed tuning-frequency
  prior and the first finite WavePeriodEstimator value.

For a frozen sample interval h, leak lambda, d=exp(-lambda*h) and
q=(1-d)/lambda, the implemented recurrences imply

    H_v/H_eta = (1 - d z^-1)/q.

After the source code's subtraction of lambda^2, a sinusoid of angular
frequency omega therefore gives

    omega_hat/omega
      = [lambda*h*exp(-lambda*h/2)/(1-exp(-lambda*h))]
        * sinc(omega*h/2).

The two high-pass stages cancel from this variance ratio.  The numerical bound
below uses the repository's validated transcendental layer and outward-rounded
interval arithmetic; ordinary libm sin/exp are not used in the enclosure.

The startup result is deliberately structural rather than replay-derived.
Shipping update order calls update_tuner(..., tuner_frequency_hz_()) before the
current sample updates WavePeriodEstimator.  Until getFrequencyHz() is finite,
tuner_frequency_hz_() returns the fixed TUNE_FREQ_PRIOR_HZ; once a finite
positive estimator value exists, the following valid sample can consume it.
This takeover does not wait for WavePeriodEstimator::isReady().  Conversely the
filter's TunerReady stage is based on SeaStateAutoTuner readiness and can occur
while the period estimator is still inside its mandatory 6/lambda integrator
settling interval.

This is not the full SEA0 T_p -> tuner-T_z certificate.  It does not enclose a
multimodal response spectrum, finite EW moment transients after estimator
takeover, the log-period EMA over a changing sea, or target-platform
floating-point/libm implementation error.  Those remain explicit downstream
obligations.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

from ou3_interval import Interval
import ou3_validated_transcendentals as VT


SCHEMA_VERSION = "OU3_SEA3_WAVE_PERIOD_FRONTEND_V2"
REPO = Path(__file__).resolve().parents[1]
DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
FILTER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
ESTIMATOR = REPO / "src" / "tuner" / "WavePeriodEstimator.h"
TUNER = REPO / "src" / "tuner" / "SeaStateAutoTuner.h"
LIMITS = REPO / "src" / "tuner" / "SeaStateAdaptationLimits.h"

# Deliberately loose decimal enclosure of pi.  The endpoint decimal literals
# straddle mathematical pi by far more than one binary64 ulp; outward_bounds
# then widens once more.  This avoids treating math.pi as an exact real.
PI = Interval.outward_bounds(3.141592653589793, 3.141592653589794)
ONE = Interval.point(1.0)
TWO = Interval.point(2.0)
SIX = Interval.point(6.0)


def _decimal_interval(value: float) -> Interval:
    # JSON/source decimal constants are theorem values; one binary64 step on
    # either side safely contains the represented decimal real.
    return Interval.outward_bounds(float(value), float(value))


def _const_float(text: str, name: str) -> float:
    match = re.search(rf"\b{name}\s*=\s*([^;]+);", text)
    if not match:
        raise ValueError(f"could not find {name}")
    expr = match.group(1).replace("f", "").replace("F", "").strip()
    if "/" in expr:
        pieces = [part.strip() for part in expr.split("/")]
        if len(pieces) != 2:
            raise ValueError(f"unsupported source expression for {name}: {expr}")
        return float(pieces[0]) / float(pieces[1])
    return float(expr)


def _config_default_float(text: str, name: str) -> float:
    match = re.search(rf"\bfloat\s+{name}\s*=\s*([0-9.eE+-]+)f?\s*;", text)
    if not match:
        raise ValueError(f"could not find Config default {name}")
    return float(match.group(1))


def _default_high_pass_hz(text: str) -> float:
    match = re.search(
        r"explicit\s+WavePeriodEstimator\(float\s+high_pass_hz\s*=\s*([0-9.]+)f",
        text,
    )
    if not match:
        raise ValueError("could not find WavePeriodEstimator high-pass default")
    return float(match.group(1))


def _source_recurrence_parity(text: str) -> dict[str, bool]:
    required = {
        "exact_decay": "const float decay = std::exp(-lambda_ * dt_sec);",
        "exact_gain": (
            "const float gain = (lambda_ > 1e-9f) ? "
            "((1.0f - decay) / lambda_) : dt_sec;"
        ),
        "velocity_recurrence": (
            "velocity_ = decay * velocity_ + gain * stage2;"
        ),
        "elevation_recurrence": (
            "elevation_ = decay * elevation_ + gain * velocity_;"
        ),
        "variance_ratio": "const float ratio_sq = velocity_var / elevation_var;",
        "leak_subtraction": "const float omega_sq = ratio_sq - lambda_ * lambda_;",
        "integrator_settle_gate": "const float settle_sec = 6.0f / lambda_;",
        "integrator_settle_return": "if (elapsed_sec_ < settle_sec) return;",
        "period_ready_gate": (
            "return std::isfinite(log_period_sec_) && weight_ > 0.5f;"
        ),
    }
    compact = " ".join(text.split())
    return {
        key: " ".join(fragment.split()) in compact
        for key, fragment in required.items()
    }


def _source_startup_parity(filter_text: str, tuner_text: str) -> dict[str, bool]:
    compact_filter = " ".join(filter_text.split())
    compact_tuner = " ".join(tuner_text.split())

    selector_start = filter_text.find("float tuner_frequency_hz_() const")
    selector_end = filter_text.find("void resetTrackingState_()", selector_start)
    selector = (
        filter_text[selector_start:selector_end]
        if selector_start >= 0 and selector_end > selector_start
        else ""
    )

    stage_start = filter_text.find("void update_tuner(")
    stage_end = filter_text.find("void adapt_mekf(", stage_start)
    stage = (
        filter_text[stage_start:stage_end]
        if stage_start >= 0 and stage_end > stage_start
        else ""
    )

    tuner_call = "update_tuner(dt, a_vert_measurement, tuner_frequency_hz_());"
    period_call = "wave_period_.update(dt, wave_period_input_ms2_(direction_accel));"
    tuner_pos = filter_text.find(tuner_call)
    period_pos = filter_text.find(period_call)

    required_filter = {
        "selector_reads_wave_period_frequency": (
            "const float wave_hz = wave_period_.getFrequencyHz();"
        ),
        "selector_accepts_first_finite_positive_frequency": (
            "if (std::isfinite(wave_hz) && wave_hz > 0.0f) return wave_hz;"
        ),
        "selector_falls_back_to_fixed_prior": "return tune_freq_prior_hz_;",
        "outer_config_forwards_tuner_warmup": (
            "impl_.setOnlineTuneWarmupSec(cfg_.online_tune_warmup_sec);"
        ),
        "tunerwarm_checks_tuner_frequency_ready": (
            "if (!tuner_.isFreqReady()) return;"
        ),
        "tunerwarm_checks_tuner_ready": "if (tuner_.isReady()) {",
    }
    out = {
        key: " ".join(fragment.split()) in compact_filter
        for key, fragment in required_filter.items()
    }
    out["selector_does_not_wait_for_wave_period_isReady"] = (
        bool(selector) and "isReady" not in selector
    )
    out["tunerwarm_does_not_wait_for_wave_period_isReady"] = (
        bool(stage) and "wave_period_.isReady" not in stage
    )
    out["tuner_update_precedes_current_sample_wave_period_update"] = (
        tuner_pos >= 0 and period_pos >= 0 and tuner_pos < period_pos
    )
    out["tuner_frequency_input_is_stored_immediately"] = (
        "frequency_hz = f_eff;" in compact_tuner
    )
    out["tuner_ready_is_own_frequency_and_variance_gate"] = (
        "inline bool isReady() const { return isVarReady() && isFreqReady(); }"
        in compact_tuner
    )
    out["debiased_variance_ready_threshold_is_1e_6"] = (
        "inline bool isReady() const { return weight > 1e-6f; }" in compact_tuner
    )
    return out


def build(repo_root: Path = REPO) -> dict[str, Any]:
    repo_root = Path(repo_root)
    domain = json.loads(
        (repo_root / "tools" / "ou3_proof_operating_domain.json").read_text(
            encoding="utf-8"
        )
    )
    filter_text = (
        repo_root / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
    ).read_text(encoding="utf-8")
    estimator_text = (
        repo_root / "src" / "tuner" / "WavePeriodEstimator.h"
    ).read_text(encoding="utf-8")
    tuner_text = (
        repo_root / "src" / "tuner" / "SeaStateAutoTuner.h"
    ).read_text(encoding="utf-8")
    limits_text = (
        repo_root / "src" / "tuner" / "SeaStateAdaptationLimits.h"
    ).read_text(encoding="utf-8")

    dt_value = float(domain["configured_runtime"]["imu_dt_s"])
    hp_value = _default_high_pass_hz(estimator_text)
    f_min_value = _const_float(filter_text, "MIN_TUNE_FREQ_HZ")
    f_max_value = _const_float(filter_text, "MAX_TUNE_FREQ_HZ")
    prior_value = _const_float(filter_text, "TUNE_FREQ_PRIOR_HZ")
    warmup_value = _config_default_float(filter_text, "online_tune_warmup_sec")
    horizon_max_value = _const_float(limits_text, "kDynamicEmaHorizonMaxSec")

    dt = _decimal_interval(dt_value)
    hp = _decimal_interval(hp_value)
    freq = Interval.outward_bounds(f_min_value, f_max_value)
    prior = _decimal_interval(prior_value)
    warmup = _decimal_interval(warmup_value)
    horizon_max = _decimal_interval(horizon_max_value)

    lam = TWO * PI * hp
    x = lam * dt
    half_x = x / TWO

    # c_leak = x*exp(-x/2)/(1-exp(-x)).  Arguments are ~6e-4 and remain well
    # inside the validated exponential layer's audited |x|<=1/2 domain.
    exp_neg_half = VT.exp_interval(-half_x)
    exp_neg_x = VT.exp_interval(-x)
    leak_factor = x * exp_neg_half / (ONE - exp_neg_x)

    # omega*h/2 = pi*f*h.  On 0.03..1.2 Hz and h=5 ms this is <0.019 rad,
    # inside the validated sinc monotonicity range by two orders of magnitude.
    phase_half = PI * freq * dt
    sample_factor = VT.sinc_interval(phase_half)
    omega_hat_over_omega = leak_factor * sample_factor
    period_hat_over_period = omega_hat_over_omega.reciprocal()

    # No WavePeriodEstimator moment can be accepted before 6/lambda.  This is a
    # lower bound on first possible finite period, not a claim that a valid
    # variance ratio must exist exactly at that instant.
    estimator_settle = SIX / lam

    # SeaStateAutoTuner's first debiased moment update starts from zero weight.
    # The dynamic horizon is never above the source max, so this is a uniform
    # lower bound on that first weight increment for a valid sample.
    first_tuner_weight_lower = ONE - VT.exp_interval(-(dt / horizon_max))

    recurrence_parity = _source_recurrence_parity(estimator_text)
    startup_parity = _source_startup_parity(filter_text, tuner_text)
    parity = {**recurrence_parity, **startup_parity}

    startup_domain = domain.get("startup", {})
    domain_parity = {
        "configured_warmup_matches_source_default": (
            float(startup_domain.get("online_tune_warmup_sec", float("nan")))
            == warmup_value
        ),
        "configured_prior_matches_source_default": (
            float(startup_domain.get("tune_frequency_prior_hz", float("nan")))
            == prior_value
        ),
        "configured_high_pass_matches_source_default": (
            float(startup_domain.get("wave_period_estimator_high_pass_hz", float("nan")))
            == hp_value
        ),
        "domain_does_not_equate_tuner_and_period_readiness": (
            startup_domain.get("tuner_ready_requires_wave_period_estimator_ready")
            is False
        ),
        "domain_admits_live_before_first_period": (
            startup_domain.get("live_entry_may_precede_wave_period_estimator_first_valid_period")
            is True
        ),
    }

    return {
        "schema_version": SCHEMA_VERSION,
        "qualification": "OU3_SEA0_WAVE_PERIOD_FRONTEND_SUBCERTIFICATE",
        "proof_status": "validated_frontend_and_startup_source_causality_subcertificate",
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_operating_domain_shrunk": False,
        "SEA0_full_certificate_promoted": False,
        "P2_promoted_from_this_artifact": False,
        "P3_promoted_from_this_artifact": False,
        "source_parity": parity,
        "operating_domain_parity": domain_parity,
        "declared_inputs": {
            "imu_dt_s": dt.as_list(),
            "wave_period_high_pass_hz": hp.as_list(),
            "screened_tuning_frequency_hz": freq.as_list(),
            "fixed_tuning_frequency_prior_hz": prior.as_list(),
            "configured_online_tune_warmup_s": warmup.as_list(),
            "frequency_screen_role": (
                "committed tuning-channel range; not a physical Tp domain"
            ),
        },
        "exact_transfer_identity": {
            "H_v_over_H_eta": "(1 - d*z^-1)/q",
            "d": "exp(-lambda*h)",
            "q": "(1-d)/lambda",
            "omega_hat_over_omega": (
                "[lambda*h*exp(-lambda*h/2)/(1-exp(-lambda*h))] "
                "* sinc(omega*h/2)"
            ),
            "two_shared_high_pass_stages_cancel_from_ratio": True,
            "source_leak_square_subtraction_used": True,
        },
        "validated_intervals": {
            "lambda_per_s": lam.as_list(),
            "lambda_h": x.as_list(),
            "half_sample_phase_rad": phase_half.as_list(),
            "leak_discretization_factor": leak_factor.as_list(),
            "sampled_sinc_factor": sample_factor.as_list(),
            "omega_hat_over_omega": omega_hat_over_omega.as_list(),
            "period_hat_over_period": period_hat_over_period.as_list(),
            "wave_period_integrator_settle_lower_bound_s": estimator_settle.as_list(),
            "first_tuner_debiased_weight_increment_lower": (
                first_tuner_weight_lower.as_list()
            ),
            "validated_transcendentals_used": True,
            "ordinary_libm_transcendentals_used_for_enclosure": False,
        },
        "startup_source_language": {
            "frequency_modes": [
                "fixed_prior_until_first_finite_positive_wave_period_frequency",
                "wave_period_estimator_after_takeover",
            ],
            "tuner_ready_requires_wave_period_estimator_ready": False,
            "live_entry_may_precede_wave_period_estimator_first_valid_period": True,
            "wave_period_takeover_waits_for_isReady": False,
            "tuner_consumes_previous_sample_wave_period_state": True,
            "current_sample_wave_period_update_occurs_after_tuner_update": True,
            "first_newly_finite_wave_period_can_affect_tuner_no_earlier_than_next_valid_sample": True,
            "first_valid_tuner_update_can_satisfy_debiased_variance_ready_threshold": (
                first_tuner_weight_lower.lo > 1.0e-6
            ),
            "wave_period_integrator_settle_is_lower_bound_not_readiness_time": True,
            "required_P2_consequence": (
                "SEA3 source language must carry prior-frequency and estimator-frequency "
                "branches plus the one-sample takeover edge; TunerReady cannot be used "
                "as a proxy for settled sea-period state"
            ),
        },
        "interpretation": {
            "discrete_frontend_warping_is_current_limiter": False,
            "startup_prior_to_estimator_takeover_is_now_source_certified": True,
            "multimodal_response_moment_enclosure_still_required": True,
            "finite_EW_moment_transient_still_required": True,
            "canonical_log_period_EMA_still_required": True,
            "target_float_libm_rounding_still_required": True,
            "surface_Tz_or_sinusoid_period_may_replace_tuner_Tz": False,
        },
        "next_obligation": (
            "certify the directional vessel/IMU response-weighted finite-band SEA3 "
            "moments and propagate them through the retained two-mode prior/estimator "
            "source language, including finite EW moment/log-period dynamics, before "
            "using sea physics to prune P2"
        ),
    }


def validate(payload: dict[str, Any]) -> list[str]:
    failures: list[str] = []
    if payload.get("schema_version") != SCHEMA_VERSION:
        failures.append("schema mismatch")
    if payload.get("trajectory_replay_used") is not False:
        failures.append("front-end certificate must be replay free")
    if payload.get("SEA0_full_certificate_promoted") is not False:
        failures.append("front-end lemma must not promote full SEA0")
    parity = payload.get("source_parity", {})
    missing = [name for name, matched in parity.items() if matched is not True]
    if missing:
        failures.append(f"WavePeriodEstimator/tuner source semantics changed: {missing}")
    domain_parity = payload.get("operating_domain_parity", {})
    domain_missing = [
        name for name, matched in domain_parity.items() if matched is not True
    ]
    if domain_missing:
        failures.append(f"startup operating-domain parity changed: {domain_missing}")
    exact = payload.get("exact_transfer_identity", {})
    if exact.get("two_shared_high_pass_stages_cancel_from_ratio") is not True:
        failures.append("high-pass cancellation identity missing")
    intervals = payload.get("validated_intervals", {})
    ratio = intervals.get("omega_hat_over_omega", [])
    period = intervals.get("period_hat_over_period", [])
    settle = intervals.get("wave_period_integrator_settle_lower_bound_s", [])
    first_weight = intervals.get("first_tuner_debiased_weight_increment_lower", [])
    if len(ratio) != 2 or not (0.99 < float(ratio[0]) <= float(ratio[1]) <= 1.001):
        failures.append("unexpected frequency-warping enclosure")
    if len(period) != 2 or not (0.999 < float(period[0]) <= float(period[1]) < 1.001):
        failures.append("unexpected period-warping enclosure")
    if len(settle) != 2 or not (47.0 < float(settle[0]) <= float(settle[1]) < 49.0):
        failures.append("unexpected WavePeriodEstimator settle lower bound")
    if len(first_weight) != 2 or not (float(first_weight[0]) > 1.0e-6):
        failures.append("first tuner moment update no longer clears readiness threshold")
    if intervals.get("validated_transcendentals_used") is not True:
        failures.append("validated transcendental layer not used")
    startup = payload.get("startup_source_language", {})
    if startup.get("tuner_ready_requires_wave_period_estimator_ready") is not False:
        failures.append("TunerReady was incorrectly tied to WavePeriodEstimator::isReady")
    if startup.get("live_entry_may_precede_wave_period_estimator_first_valid_period") is not True:
        failures.append("startup prior branch was incorrectly removed")
    if startup.get("wave_period_takeover_waits_for_isReady") is not False:
        failures.append("wave-period takeover must use first finite value, not isReady")
    if startup.get("tuner_consumes_previous_sample_wave_period_state") is not True:
        failures.append("one-sample tuner/period causal order was lost")
    if startup.get(
        "first_valid_tuner_update_can_satisfy_debiased_variance_ready_threshold"
    ) is not True:
        failures.append("SeaStateAutoTuner readiness relation changed")
    interpretation = payload.get("interpretation", {})
    if interpretation.get("surface_Tz_or_sinusoid_period_may_replace_tuner_Tz") is not False:
        failures.append("single-frequency lemma must not replace tuner Tz")
    return failures


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--repo-root", type=Path, default=REPO)
    args = parser.parse_args()

    payload = build(args.repo_root)
    failures = validate(payload)
    payload["validation_pass"] = not failures
    payload["validation_failures"] = failures
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.write_text(text, encoding="utf-8")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
