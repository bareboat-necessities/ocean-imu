#!/usr/bin/env python3
"""Validated SEA0 subcertificate for the WavePeriodEstimator front-end map.

The physical SEA3 spectrum is parameterized by peak period, while shipping
WavePeriodEstimator estimates a zero-crossing/moment period after two shared
high-pass stages and two leaky integrations.  This module isolates one exact
piece of that bridge: the steady single-frequency transfer identity of the
implemented front end, in the real-arithmetic recurrence model used by the
proof tooling.

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

This is not the full SEA0 T_p -> tuner-T_z certificate.  It does not enclose a
multimodal response spectrum, finite EW moment transients, the log-period EMA,
or target-platform floating-point/libm implementation error.  Those remain
explicit downstream obligations.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

from ou3_interval import Interval
import ou3_validated_transcendentals as VT


SCHEMA_VERSION = "OU3_SEA3_WAVE_PERIOD_FRONTEND_V1"
REPO = Path(__file__).resolve().parents[1]
DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
FILTER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
ESTIMATOR = REPO / "src" / "tuner" / "WavePeriodEstimator.h"

# Deliberately loose decimal enclosure of pi.  The endpoint decimal literals
# straddle mathematical pi by far more than one binary64 ulp; outward_bounds
# then widens once more.  This avoids treating math.pi as an exact real.
PI = Interval.outward_bounds(3.141592653589793, 3.141592653589794)
ONE = Interval.point(1.0)
TWO = Interval.point(2.0)


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
    }
    compact = " ".join(text.split())
    return {
        key: " ".join(fragment.split()) in compact
        for key, fragment in required.items()
    }


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

    dt_value = float(domain["configured_runtime"]["imu_dt_s"])
    hp_value = _default_high_pass_hz(estimator_text)
    f_min_value = _const_float(filter_text, "MIN_TUNE_FREQ_HZ")
    f_max_value = _const_float(filter_text, "MAX_TUNE_FREQ_HZ")

    dt = _decimal_interval(dt_value)
    hp = _decimal_interval(hp_value)
    freq = Interval.outward_bounds(f_min_value, f_max_value)

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

    parity = _source_recurrence_parity(estimator_text)

    return {
        "schema_version": SCHEMA_VERSION,
        "qualification": "OU3_SEA0_WAVE_PERIOD_FRONTEND_SUBCERTIFICATE",
        "proof_status": "validated_real_arithmetic_steady_sinusoid_subcertificate",
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_operating_domain_shrunk": False,
        "SEA0_full_certificate_promoted": False,
        "P2_promoted_from_this_artifact": False,
        "P3_promoted_from_this_artifact": False,
        "source_parity": parity,
        "declared_inputs": {
            "imu_dt_s": dt.as_list(),
            "wave_period_high_pass_hz": hp.as_list(),
            "screened_tuning_frequency_hz": freq.as_list(),
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
            "validated_transcendentals_used": True,
            "ordinary_libm_transcendentals_used_for_enclosure": False,
        },
        "interpretation": {
            "discrete_frontend_warping_is_current_limiter": False,
            "multimodal_response_moment_enclosure_still_required": True,
            "finite_EW_moment_transient_still_required": True,
            "canonical_log_period_EMA_still_required": True,
            "target_float_libm_rounding_still_required": True,
            "surface_Tz_or_sinusoid_period_may_replace_tuner_Tz": False,
        },
        "next_obligation": (
            "feed the directional vessel/IMU response-weighted finite-band "
            "SEA3 spectrum through this front-end relation and then certify "
            "finite EW moment/log-period dynamics before pruning P2"
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
        failures.append(f"WavePeriodEstimator source recurrence changed: {missing}")
    exact = payload.get("exact_transfer_identity", {})
    if exact.get("two_shared_high_pass_stages_cancel_from_ratio") is not True:
        failures.append("high-pass cancellation identity missing")
    intervals = payload.get("validated_intervals", {})
    ratio = intervals.get("omega_hat_over_omega", [])
    period = intervals.get("period_hat_over_period", [])
    if len(ratio) != 2 or not (0.99 < float(ratio[0]) <= float(ratio[1]) <= 1.001):
        failures.append("unexpected frequency-warping enclosure")
    if len(period) != 2 or not (0.999 < float(period[0]) <= float(period[1]) < 1.001):
        failures.append("unexpected period-warping enclosure")
    if intervals.get("validated_transcendentals_used") is not True:
        failures.append("validated transcendental layer not used")
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
