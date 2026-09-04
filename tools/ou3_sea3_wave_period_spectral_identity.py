#!/usr/bin/env python3
"""Analytical steady-spectrum identity for the OU-III WavePeriodEstimator.

This is the reusable response-independent result from the superseded parallel
SEA0 study. It does not choose a vessel RAO and it does not numerically screen
a sea. It states what the continuous-time steady-state counterpart of the
shipping estimator computes for any nonnegative input specific-force spectrum
S_a with finite weighted numerator and a finite, strictly positive weighted
denominator.

With two shared high-pass stages H_hp=s/(s+lambda) and two leaky integrations,

    H_v(s)   = s^2/(s+lambda)^3,
    H_eta(s) = s^2/(s+lambda)^4.

Hence

    |H_v(iw)|^2   = (w^2+lambda^2) W(w),
    |H_eta(iw)|^2 = W(w),
    W(w) = w^4/(w^2+lambda^2)^4.

For S_a(w)>=0 with

    0 < int W(w) S_a(w) dw < infinity,
    int w^2 W(w) S_a(w) dw < infinity,

we have

    sigma_v^2/sigma_eta^2 - lambda^2
      = [int w^2 W(w) S_a(w) dw] / [int W(w) S_a(w) dw].

So the leak subtraction is exactly a response-weighted mean-square frequency
identity in the continuous-time steady-state model; it is not a narrow-band or
single-sinusoid approximation. The exact discrete finite-EWMA estimator still
requires its own transient enclosure before this lemma can prune P2.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ESTIMATOR = REPO / "src" / "tuner" / "WavePeriodEstimator.h"
SCHEMA = 2
QUALIFICATION = "OU3_SEA3_WAVE_PERIOD_STEADY_SPECTRAL_IDENTITY"

SOURCE_PARITY_KEYS = (
    "lambda_constructor",
    "decay_assignment",
    "gain_assignment",
    "stage1_assignment",
    "accel_previous_update_after_stage1",
    "stage2_assignment",
    "stage1_previous_update_after_stage2",
    "stage1_state_update_after_stage2",
    "stage2_state_update_after_stage2",
    "velocity_leaky_integration",
    "elevation_uses_updated_velocity",
    "variance_ratio",
    "leak_square_subtraction",
)


def source_parity(estimator_path: Path = ESTIMATOR) -> dict[str, bool]:
    """Check the exact source relationships used to derive the transfer identity."""
    text = Path(estimator_path).read_text(encoding="utf-8")
    compact = " ".join(text.split())

    fragments = {
        "lambda_constructor": (
            "lambda_(2.0f * 3.14159265358979323846f * "
            "std::max(1e-4f, high_pass_hz))"
        ),
        "decay_assignment": "const float decay = std::exp(-lambda_ * dt_sec);",
        "gain_assignment": (
            "const float gain = (lambda_ > 1e-9f) ? "
            "((1.0f - decay) / lambda_) : dt_sec;"
        ),
        "stage1_assignment": (
            "const float stage1 = decay * "
            "(high_pass_1_ + vertical_accel_ms2 - accel_prev_);"
        ),
        "accel_previous_update_after_stage1": "accel_prev_ = vertical_accel_ms2;",
        "stage2_assignment": (
            "const float stage2 = decay * "
            "(high_pass_2_ + stage1 - high_pass_1_prev_);"
        ),
        "stage1_previous_update_after_stage2": "high_pass_1_prev_ = stage1;",
        "stage1_state_update_after_stage2": "high_pass_1_ = stage1;",
        "stage2_state_update_after_stage2": "high_pass_2_ = stage2;",
        "velocity_leaky_integration": "velocity_ = decay * velocity_ + gain * stage2;",
        "elevation_uses_updated_velocity": "elevation_ = decay * elevation_ + gain * velocity_;",
        "variance_ratio": "const float ratio_sq = velocity_var / elevation_var;",
        "leak_square_subtraction": "const float omega_sq = ratio_sq - lambda_ * lambda_;",
    }
    present = {
        key: " ".join(fragment.split()) in compact
        for key, fragment in fragments.items()
    }

    # The transfer relation also depends on update order, not only presence.
    ordered_fragments = [
        fragments["decay_assignment"],
        fragments["gain_assignment"],
        fragments["stage1_assignment"],
        fragments["accel_previous_update_after_stage1"],
        fragments["stage2_assignment"],
        fragments["stage1_previous_update_after_stage2"],
        fragments["stage1_state_update_after_stage2"],
        fragments["stage2_state_update_after_stage2"],
        fragments["velocity_leaky_integration"],
        fragments["elevation_uses_updated_velocity"],
    ]
    positions = [compact.find(" ".join(fragment.split())) for fragment in ordered_fragments]
    order_ok = all(p >= 0 for p in positions) and all(
        positions[i] < positions[i + 1] for i in range(len(positions) - 1)
    )
    if not order_ok:
        # Keep the exact key set stable; one failed order relation invalidates
        # the recurrence keys rather than creating an unvalidated extra flag.
        for key in (
            "decay_assignment", "gain_assignment", "stage1_assignment",
            "accel_previous_update_after_stage1", "stage2_assignment",
            "stage1_previous_update_after_stage2", "stage1_state_update_after_stage2",
            "stage2_state_update_after_stage2", "velocity_leaky_integration",
            "elevation_uses_updated_velocity",
        ):
            present[key] = False
    return present


def build(estimator_path: Path = ESTIMATOR) -> dict:
    """Build the analytical identity contract without trajectory replay or RAO selection."""
    parity = source_parity(estimator_path)
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "trajectory_replay_used": False,
        "single_frequency_assumption_used": False,
        "single_RAO_used": False,
        "finite_RAO_grid_used": False,
        "source_parity": parity,
        "continuous_time_steady_state_identity": {
            "H_v": "s^2/(s+lambda)^3",
            "H_eta": "s^2/(s+lambda)^4",
            "W": "omega^4/(omega^2+lambda^2)^4",
            "velocity_weight_relation": "|H_v|^2=(omega^2+lambda^2)W",
            "elevation_weight_relation": "|H_eta|^2=W",
            "leak_subtracted_variance_ratio": (
                "sigma_v^2/sigma_eta^2-lambda^2 = "
                "int omega^2 W S_a / int W S_a"
            ),
            "input_spectrum_nonnegative": True,
            "weighted_denominator_finite_and_strictly_positive_required": True,
            "weighted_second_moment_finite_required": True,
            "holds_for_any_input_spectrum_satisfying_these_preconditions": True,
            "narrow_band_approximation": False,
            "single_sinusoid_approximation": False,
        },
        "promotion": {
            "SEA0_full_certificate_promoted": False,
            "P2_pruning_promoted": False,
            "finite_EWMA_transient_enclosed": False,
            "discrete_estimator_identified_with_continuous_steady_state": False,
        },
        "next_obligation": (
            "propagate the coupled JONSWAP-sea/RAO response spectrum through the exact discrete finite-EWMA/log-period estimator state before using the identity for source-history pruning"
        ),
    }


def validate(d: dict) -> list[str]:
    """Fail closed if the identity is narrowed, promoted, or loses source parity."""
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("trajectory_replay_used") is not False:
        f.append("spectral identity must be replay free")
    for key in (
        "single_frequency_assumption_used", "single_RAO_used", "finite_RAO_grid_used"
    ):
        if d.get(key) is not False:
            f.append(f"spectral identity overclaim guard changed: {key}")

    parity = d.get("source_parity", {})
    if set(parity) != set(SOURCE_PARITY_KEYS) or not SOURCE_PARITY_KEYS:
        f.append("shipping WavePeriodEstimator source-parity key set changed")
    elif not all(parity.get(key) is True for key in SOURCE_PARITY_KEYS):
        f.append("shipping WavePeriodEstimator source parity failed")

    ident = d.get("continuous_time_steady_state_identity", {})
    for key in (
        "input_spectrum_nonnegative",
        "weighted_denominator_finite_and_strictly_positive_required",
        "weighted_second_moment_finite_required",
        "holds_for_any_input_spectrum_satisfying_these_preconditions",
    ):
        if ident.get(key) is not True:
            f.append(f"spectral-identity precondition/quantifier changed: {key}")
    if ident.get("narrow_band_approximation") is not False:
        f.append("exact identity was weakened to narrow-band approximation")
    if ident.get("single_sinusoid_approximation") is not False:
        f.append("exact identity was weakened to a single-sinusoid approximation")

    p = d.get("promotion", {})
    for key in (
        "SEA0_full_certificate_promoted", "P2_pruning_promoted",
        "finite_EWMA_transient_enclosed",
        "discrete_estimator_identified_with_continuous_steady_state",
    ):
        if p.get(key) is not False:
            f.append(f"steady spectral identity overpromoted: {key}")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--estimator", type=Path, default=ESTIMATOR)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.estimator)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "arbitrary_spectrum_identity": d["continuous_time_steady_state_identity"]["holds_for_any_input_spectrum_satisfying_these_preconditions"],
        "positive_weighted_denominator_required": d["continuous_time_steady_state_identity"]["weighted_denominator_finite_and_strictly_positive_required"],
        "source_parity": d["source_parity"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
