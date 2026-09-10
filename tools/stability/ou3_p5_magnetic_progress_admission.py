#!/usr/bin/env python3
"""Minimal finite magnetic-progress admission for the OU-III P5 theorem.

The shipping timeout path can enter Live without north.  Absolute-yaw capture
therefore cannot be proved from BRMM motion recurrence alone: MagAutoTuner
advances ``accepted_window_sec_`` only on accepted attempts, and rejected
attempts can make accepted-packet timing arbitrarily sparse in accepted time.

This module adds exactly the observability progress coordinate needed by the
shipping implementation and nothing about sea amplitude, vessel motion, tuner
frequency, tau, sigma or R_S.  It is a theorem-source admission layer, not a
filter change and not a claim that every physical deployment satisfies the
condition.

For ``with_mag=true`` an admitted history must, after the first eligible
magnetic attempt and no later than a declared finite horizon, make the literal
shipping MagAutoTuner state satisfy

    accepted_count_ >= min_samples
    accepted_window_sec_ >= min_window_sec

and the accumulated mean must satisfy the same finite/nondegenerate horizontal
checks used by ``tryFinalize_``.  With the default unweighted acquisition this
is sufficient for ``ready_=true`` because min_effective_weight is zero and hard
iron is optional/off by default.  ``with_mag=false`` needs no such premise.

The default horizon is deliberately the existing 150 s shipping startup ceiling.
It is an explicit theorem-admission constant, not inferred from packet ODR.  A
caller may widen it to any larger finite value; the capture-time bound widens by
the same amount.
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
MAG = REPO / "src" / "tuner" / "MagAutoTuner.h"
SCHEMA = 1
QUALIFICATION = "OU3_P5_MAGNETIC_PROGRESS_ADMISSION_V1"


def _cfg_number(text: str, name: str) -> float:
    m = re.search(rf"(?:float|int)\s+{re.escape(name)}\s*=\s*([0-9.eE+\-/ ]+)f?\s*;", text)
    if not m:
        raise RuntimeError(f"cannot extract shipping Config::{name}")
    raw = m.group(1).strip()
    # All values consumed here are plain literals in the shipping Config.
    if "/" in raw:
        a, b = raw.split("/", 1)
        return float(a.strip()) / float(b.strip())
    return float(raw)


def build(progress_horizon_s: float | None = None) -> dict:
    wrapper = WRAPPER.read_text(encoding="utf-8")
    mag = MAG.read_text(encoding="utf-8")

    startup_timeout = _cfg_number(wrapper, "proxy_startup_timeout_sec")
    min_samples = int(_cfg_number(wrapper, "mag_min_samples"))
    min_window_s = _cfg_number(wrapper, "mag_min_window_sec")
    min_effective_weight = _cfg_number(wrapper, "mag_min_effective_weight")
    estimate_hard_iron_default = "bool  mag_estimate_hard_iron = false;" in wrapper
    quality_weighting_default = "bool  mag_enable_quality_weighting = false;" in wrapper

    horizon = startup_timeout if progress_horizon_s is None else float(progress_horizon_s)
    if not (math.isfinite(horizon) and horizon > 0.0):
        raise ValueError("finite positive magnetic progress horizon required")

    wrapper_parity = {
        "timeout_can_enter_live_without_north": (
            "const bool ready_by_timeout" in wrapper
            and "mag_gravity_aligned_branch_;" in wrapper
        ),
        "postlive_first_north_keeps_accumulating": (
            "if (!mag_ref_set_)" in wrapper
            and "mag_auto_tuner_.addSampleWithTiltQuatDt(" in wrapper
            and "if (stage_ != Stage::Live)" in wrapper
            and "impl_.mekf().set_quaternion_boat(q_new);" in wrapper
        ),
        "north_lock_state_set_after_success": (
            "mag_ref_set_ = true;" in wrapper
            and "mag_north_lock_time_sec_ = t_;" in wrapper
        ),
    }

    tuner_parity = {
        "count_gate": "accepted_count_ < std::max(1, cfg_.min_samples)" in mag,
        "window_gate": "accepted_window_sec_ < cfg_.min_window_sec" in mag,
        "accepted_window_only_advances_on_accepted_sample": (
            "accepted_window_sec_ += dt_use;" in mag
            and mag.index("accepted_window_sec_ += dt_use;") > mag.index("if (w < cfg_.min_sample_weight)")
        ),
        "accepted_count_advances_after_window": (
            "++accepted_count_;" in mag
            and mag.index("++accepted_count_;") > mag.index("accepted_window_sec_ += dt_use;")
        ),
        "mean_must_be_finite_nonzero": (
            "if (!mean.allFinite())" in mag
            and "mean_norm > cfg_.mag_norm_min" in mag
        ),
        "horizontal_fraction_gate": "horiz_frac < cfg_.min_horizontal_fraction" in mag,
        "ready_from_finite_reference": (
            "ready_ =" in mag
            and "mag_world_ref_.allFinite()" in mag
            and "mag_world_ref_.norm() > cfg_.mag_norm_min" in mag
        ),
    }

    default_path_sufficiency = bool(
        min_effective_weight == 0.0
        and estimate_hard_iron_default
        and quality_weighting_default
        and all(wrapper_parity.values())
        and all(tuner_parity.values())
    )

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "filter_changed": False,
        "quality_gates_changed": False,
        "BRMM_motion_caps_changed": False,
        "BIAS_family_changed": False,
        "canonical_motion_source_remains_COMPLETE_BRMM": True,
        "scope": "additional observability admission only when with_mag=true and absolute/full-heading capture is claimed",
        "shipping": {
            "startup_timeout_s": startup_timeout,
            "min_accepted_samples": min_samples,
            "min_accepted_window_s": min_window_s,
            "min_effective_weight_default": min_effective_weight,
            "quality_weighting_default_off": quality_weighting_default,
            "hard_iron_estimation_default_off": estimate_hard_iron_default,
            "wrapper_parity": wrapper_parity,
            "MagAutoTuner_parity": tuner_parity,
        },
        "admission": {
            "with_mag_false_requires_progress": False,
            "with_mag_true_progress_horizon_s": horizon,
            "clock_origin": "first eligible magnetic attempt at or after startup; progress may finish before or after goLive",
            "same_history_required": True,
            "accepted_count_lower": min_samples,
            "accepted_window_sec_lower": min_window_s,
            "accepted_mean_finite": True,
            "accepted_mean_norm_strictly_above_mag_norm_min": True,
            "accepted_mean_horizontal_fraction_at_least_shipping_min": True,
            "packet_recurrence_alone_is_not_substitute": True,
            "no_heading_maneuver_required_when_hard_iron_estimation_is_off": True,
        },
        "DEFAULT_SHIPPING_PATH_PROGRESS_IMPLIES_MAGAUTOTUNER_READY": default_path_sufficiency,
        "FINITE_NORTH_EVENT_UNDER_ADMISSION_CLOSED": default_path_sufficiency,
        "finite_north_event_time_upper_from_progress_origin_s": horizon if default_path_sufficiency else None,
        "global_physical_deployment_left_inclusion_claimed": False,
        "P4_PASS": False,
        "P5_PASS": False,
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for k in (
        "canonical_motion_source_remains_COMPLETE_BRMM",
        "DEFAULT_SHIPPING_PATH_PROGRESS_IMPLIES_MAGAUTOTUNER_READY",
        "FINITE_NORTH_EVENT_UNDER_ADMISSION_CLOSED",
    ):
        if d.get(k) is not True:
            f.append(k + " not true")
    for k in (
        "filter_changed", "quality_gates_changed", "BRMM_motion_caps_changed",
        "BIAS_family_changed", "global_physical_deployment_left_inclusion_claimed",
        "P4_PASS", "P5_PASS",
    ):
        if d.get(k) is not False:
            f.append(k + " not false")
    s = d.get("shipping", {})
    if s.get("min_accepted_samples") != 128:
        f.append("shipping min accepted samples changed")
    if s.get("min_accepted_window_s") != 15.0:
        f.append("shipping min accepted window changed")
    if s.get("startup_timeout_s") != 150.0:
        f.append("shipping startup timeout changed")
    if not all(s.get("wrapper_parity", {}).values()):
        f.append("wrapper magnetic progress parity failed")
    if not all(s.get("MagAutoTuner_parity", {}).values()):
        f.append("MagAutoTuner progress parity failed")
    a = d.get("admission", {})
    if a.get("with_mag_false_requires_progress") is not False:
        f.append("no-mag branch incorrectly requires magnetic progress")
    if a.get("packet_recurrence_alone_is_not_substitute") is not True:
        f.append("packet recurrence incorrectly substituted for accepted progress")
    if not (float(a.get("with_mag_true_progress_horizon_s", 0.0)) > 0.0):
        f.append("magnetic progress horizon is not positive finite")
    return f


if __name__ == "__main__":
    d = build()
    failures = validate(d)
    print(json.dumps({**d, "validation_pass": not failures, "validation_failures": failures}, indent=2, sort_keys=True))
    raise SystemExit(1 if failures else 0)
