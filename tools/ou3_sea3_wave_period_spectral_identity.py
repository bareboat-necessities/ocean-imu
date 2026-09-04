#!/usr/bin/env python3
"""Analytical steady-spectrum identity for the OU-III WavePeriodEstimator.

This is the reusable response-independent result from the superseded parallel
SEA0 study.  It does not choose a vessel RAO and it does not numerically screen
a sea.  It states what the continuous-time steady-state counterpart of the
shipping estimator computes for *any* nonnegative input specific-force
spectrum S_a with finite weighted moments.

With two shared high-pass stages H_hp=s/(s+lambda) and two leaky integrations,

    H_v(s)   = s^2/(s+lambda)^3,
    H_eta(s) = s^2/(s+lambda)^4.

Hence

    |H_v(iw)|^2   = (w^2+lambda^2) W(w),
    |H_eta(iw)|^2 = W(w),
    W(w) = w^4/(w^2+lambda^2)^4.

For arbitrary S_a(w)>=0,

    sigma_v^2/sigma_eta^2 - lambda^2
      = [int w^2 W(w) S_a(w) dw] / [int W(w) S_a(w) dw].

So the leak subtraction is exactly a response-weighted mean-square frequency
identity in the continuous-time steady-state model; it is not a narrow-band or
single-sinusoid approximation.  The exact discrete finite-EWMA estimator still
requires its own transient enclosure before this lemma can prune P2.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ESTIMATOR = REPO / "src" / "tuner" / "WavePeriodEstimator.h"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_WAVE_PERIOD_STEADY_SPECTRAL_IDENTITY"


def source_parity(estimator_path: Path = ESTIMATOR) -> dict[str, bool]:
    """Check that shipping source still implements the required recurrence/ratio form."""
    text = Path(estimator_path).read_text(encoding="utf-8")
    compact = " ".join(text.split())
    required = {
        "velocity_leaky_integration": "velocity_ = decay * velocity_ + gain * stage2;",
        "elevation_leaky_integration": "elevation_ = decay * elevation_ + gain * velocity_;",
        "variance_ratio": "const float ratio_sq = velocity_var / elevation_var;",
        "leak_square_subtraction": "const float omega_sq = ratio_sq - lambda_ * lambda_;",
    }
    return {
        key: " ".join(fragment.split()) in compact
        for key, fragment in required.items()
    }


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
            "holds_for_any_nonnegative_input_spectrum_with_finite_weighted_moments": True,
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
    if not all(d.get("source_parity", {}).values()):
        f.append("shipping WavePeriodEstimator source parity failed")
    ident = d.get("continuous_time_steady_state_identity", {})
    if ident.get("holds_for_any_nonnegative_input_spectrum_with_finite_weighted_moments") is not True:
        f.append("arbitrary-spectrum quantifier disappeared")
    if ident.get("narrow_band_approximation") is not False:
        f.append("exact identity was weakened to narrow-band approximation")
    p = d.get("promotion", {})
    if p.get("SEA0_full_certificate_promoted") is not False:
        f.append("steady spectral identity falsely promoted SEA0")
    if p.get("P2_pruning_promoted") is not False:
        f.append("steady spectral identity falsely promoted P2 pruning")
    if p.get("finite_EWMA_transient_enclosed") is not False:
        f.append("finite estimator transient falsely claimed closed")
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
        "arbitrary_spectrum_identity": d["continuous_time_steady_state_identity"]["holds_for_any_nonnegative_input_spectrum_with_finite_weighted_moments"],
        "source_parity": d["source_parity"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
