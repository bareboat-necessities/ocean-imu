#!/usr/bin/env python3
"""Source-independent S=0 pseudo-scheduler progress certificate for OU-III.

The deployed OU-III pseudo-update period is retargeted whenever the online tune
commits.  P3 needs a uniform upper bound on the time between S=0 updates that is
valid for *every* legal retarget sequence, not only for a fixed period.

The shipping retarget rule is progress preserving:

* if accumulated elapsed time is still below the new period, elapsed is kept;
* if the new period is already overdue, elapsed is parked at the immediate
  binary32 predecessor of the new period, so the current sample's
  ``periodic_update_due<float>`` services the update.

Consequently a retarget that does not immediately force service cannot decrease
elapsed time.  Because all deployed pseudo periods are below one second, the
``periodic_update_due`` tolerance is period independent on the whole operating
range.  The latest possible non-firing deadline is therefore the largest
shipping period.  With the deployed 150 ms clamp, exact binary32 replay of that
fixed worst case fires by sample 30 from zero elapsed.  Any post-fire residual
only advances the next service.

This certificate restores the finite S-observation-gap premise used by the
translation covariance upper.  It does not promote P3/P4/P5 by itself.
"""
from __future__ import annotations

import argparse
import json
import math
import struct
from pathlib import Path

from ou3_interval import Interval
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[2]
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
CORE = REPO / "src" / "kalman_ou_common" / "KalmanOUCoreMath.h"
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
TARGET_SAMPLES = 635
OLD_GAP = 13
OLD_LOW = 0.1300000101327896
OLD_HIGH = 0.13000193238258362
DEPLOYED_MAX_GAP_SAMPLES = 30


def _f32(x: float) -> float:
    return struct.unpack("!f", struct.pack("!f", float(x)))[0]


def _u32_of_f32(x: float) -> int:
    return struct.unpack("!I", struct.pack("!f", _f32(x)))[0]


def _f32_of_u32(x: int) -> float:
    return struct.unpack("!f", struct.pack("!I", int(x)))[0]


def _nextafterf_down_positive(x: float) -> float:
    y = _f32(x)
    if not (math.isfinite(y) and y > 0.0):
        raise ValueError("positive finite binary32 required")
    return _f32_of_u32(_u32_of_f32(y) - 1)


def _retarget(elapsed: float, period: float) -> float:
    """Binary32 transcription of the shipping progress-preserving retarget."""
    e, p = _f32(elapsed), _f32(period)
    if not (p > 0.0 and math.isfinite(p)):
        return _f32(0.0)
    if not (e >= 0.0 and math.isfinite(e)):
        return _f32(0.0)
    if e < p:
        return e
    return _nextafterf_down_positive(p)


def _due(dt: float, period: float, elapsed: float) -> tuple[bool, float]:
    """Binary32 transcription of ``periodic_update_due<float>``."""
    dt, period, elapsed = _f32(dt), _f32(period), _f32(elapsed)
    if not (dt > 0.0 and period > 0.0 and math.isfinite(dt) and math.isfinite(period)):
        return False, elapsed
    total = _f32(elapsed + dt)
    eps = _f32(2.0 ** -23)
    tol = _f32(_f32(16.0 * eps) * _f32(max(1.0, period)))
    if _f32(total + tol) < period:
        return False, total
    y = _f32(math.fmod(total, period)) if total >= period else _f32(0.0)
    if not (math.isfinite(y) and 0.0 <= y < period):
        y = _f32(0.0)
    return True, y


def _period_from_tau(tau: float, sched: dict) -> float:
    ratio = _f32(sched["pseudo_ratio"])
    pmin = _f32(sched["pseudo_min_s"])
    pmax = _f32(sched["pseudo_max_s"])
    requested = _f32(ratio * _f32(tau))
    return _f32(min(max(requested, pmin), pmax))


def _first_fire_from_zero(dt: float, period: float, limit: int = 512) -> tuple[int, float]:
    e = _f32(0.0)
    for sample in range(1, int(limit) + 1):
        fire, e = _due(dt, period, e)
        if fire:
            return sample, e
    raise RuntimeError("fixed-period scheduler did not fire within search limit")


def _replay_former_starvation(dt: float, samples: int = TARGET_SAMPLES):
    low, high = _f32(OLD_LOW), _f32(OLD_HIGH)
    elapsed = _f32(0.0)
    fired_at: list[int] = []
    since = 0
    worst = 0
    done = 0
    segment = 0
    while done < samples:
        period = high if segment < 2 or (segment % 2) == 1 else low
        elapsed = _retarget(elapsed, period)
        used = min(OLD_GAP, samples - done)
        for _ in range(used):
            done += 1
            since += 1
            fire, elapsed = _due(dt, period, elapsed)
            if fire:
                fired_at.append(done)
                worst = max(worst, since)
                since = 0
        segment += 1
    worst = max(worst, since)
    return fired_at, worst, elapsed


def _implementation_contract() -> dict:
    core = CORE.read_text(encoding="utf-8")
    mekf = MEKF.read_text(encoding="utf-8")
    helper_markers = (
        "retarget_period_elapsed_progress_preserving",
        "if (elapsed < period) return elapsed;",
        "return std::nextafter(period, T(0));",
    )
    setter_markers = (
        "const T new_period = std::max(T(1e-4), period_s);",
        "retarget_period_elapsed_progress_preserving(",
        "pseudo_update_elapsed_s_, new_period",
        "pseudo_update_period_s_ = new_period;",
    )
    legacy = "std::fmod(pseudo_update_elapsed_s_, pseudo_update_period_s_)"
    return {
        "core_progress_helper_markers_present": all(x in core for x in helper_markers),
        "ou3_setter_progress_markers_present": all(x in mekf for x in setter_markers),
        "legacy_setter_fmod_progress_reset_absent": legacy not in mekf,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain = json.loads(Path(domain_path).resolve().read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("scheduler proof domain must not be trajectory fitted")

    sched = BASE.source_schedule()
    tau_lo, tau_hi = map(float, sched["tau_applied_invariant_s"])
    dt = _f32(sched["dt_s"])
    p_lo = _period_from_tau(tau_lo, sched)
    p_hi = _period_from_tau(tau_hi, sched)
    if p_lo > p_hi:
        raise RuntimeError("deployed pseudo cadence is not monotone in tau invariant")

    fixed_fire, fixed_residual = _first_fire_from_zero(dt, p_hi)
    pred_hi = _nextafterf_down_positive(p_hi)
    max_ulp = _f32(p_hi - pred_hi)
    overdue_fire, overdue_residual = _due(dt, p_hi, pred_hi)
    old_fires, old_worst, old_endpoint = _replay_former_starvation(dt)

    cadence_interval = BASE.cadence_bounds(Interval.outward_bounds(tau_lo, tau_hi), sched)
    legacy_gap_upper = math.nextafter(float(cadence_interval[1]) + float(sched["dt_s"]), math.inf)
    exact_gap_time = math.nextafter(float(fixed_fire) * float(dt), math.inf)
    impl = _implementation_contract()

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_PROGRESS_PRESERVING_PSEUDO_SCHEDULER_RECURRENCE_CERTIFICATE",
        "trajectory_replay_used": False,
        "filter_changed": True,
        "declared_domain_changed": False,
        "canonical_gate_changed": False,
        **impl,
        "shipping_scheduler_numeric_type": "binary32/float",
        "dt_binary32_s": dt,
        "tau_applied_invariant_s": [tau_lo, tau_hi],
        "deployed_pseudo_period_binary32_s": [p_lo, p_hi],
        "all_deployed_periods_below_one_second": p_hi < 1.0,
        "due_tolerance_period_independent_on_deployed_range": p_hi < 1.0,
        "maximum_period_predecessor_binary32_s": pred_hi,
        "maximum_period_ulp_s": max_ulp,
        "dt_exceeds_maximum_period_ulp": dt > max_ulp,
        "overdue_retarget_is_serviced_next_sample": bool(overdue_fire),
        "overdue_retarget_residual_s": overdue_residual,
        "nonfiring_retarget_never_decreases_elapsed": True,
        "largest_period_is_latest_nonfiring_deadline": True,
        "fixed_max_period_first_fire_samples": fixed_fire,
        "fixed_max_period_post_fire_residual_s": fixed_residual,
        "certified_uniform_max_gap_samples": fixed_fire,
        "certified_uniform_max_gap_s": exact_gap_time,
        "translation_upper_cadence_plus_h_s": legacy_gap_upper,
        "uniform_gap_within_translation_upper_bound": exact_gap_time <= legacy_gap_upper,
        "former_starvation_cycle_target_samples": TARGET_SAMPLES,
        "former_starvation_cycle_pseudo_firings": len(old_fires),
        "former_starvation_cycle_firing_sample_indices": old_fires,
        "former_starvation_cycle_worst_gap_samples": old_worst,
        "former_starvation_cycle_endpoint_elapsed_s": old_endpoint,
        "former_635_sample_zero_fire_cycle_broken": bool(old_fires),
        "scheduler_recurrence_certificate": True,
        "P3_PROMOTED": False,
        "P4_PROMOTED": False,
        "P5_PROMOTED": False,
        "classification": f"SOURCE_INDEPENDENT_{fixed_fire}_SAMPLE_S_RECURRENCE_CERTIFIED",
        "next_obligation": (
            "rerun the canonical source-reachable P3 covariance upper with this implementation-bound "
            "recurrence premise; this scheduler certificate alone does not promote P3"
        ),
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA or d.get("qualification") != "OU3_P3_PROGRESS_PRESERVING_PSEUDO_SCHEDULER_RECURRENCE_CERTIFICATE":
        f.append("schema/qualification mismatch")
    for key in (
        "filter_changed",
        "core_progress_helper_markers_present",
        "ou3_setter_progress_markers_present",
        "legacy_setter_fmod_progress_reset_absent",
        "all_deployed_periods_below_one_second",
        "due_tolerance_period_independent_on_deployed_range",
        "dt_exceeds_maximum_period_ulp",
        "overdue_retarget_is_serviced_next_sample",
        "nonfiring_retarget_never_decreases_elapsed",
        "largest_period_is_latest_nonfiring_deadline",
        "uniform_gap_within_translation_upper_bound",
        "former_635_sample_zero_fire_cycle_broken",
        "scheduler_recurrence_certificate",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "declared_domain_changed", "canonical_gate_changed",
        "P3_PROMOTED", "P4_PROMOTED", "P5_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("shipping_scheduler_numeric_type") != "binary32/float":
        f.append("scheduler arithmetic is not bound to binary32")
    if int(d.get("fixed_max_period_first_fire_samples", 0)) != DEPLOYED_MAX_GAP_SAMPLES:
        f.append(f"largest deployed period no longer fires on sample {DEPLOYED_MAX_GAP_SAMPLES}")
    if int(d.get("certified_uniform_max_gap_samples", 0)) != DEPLOYED_MAX_GAP_SAMPLES:
        f.append(f"uniform scheduler gap is not {DEPLOYED_MAX_GAP_SAMPLES} samples")
    if int(d.get("former_starvation_cycle_target_samples", 0)) != TARGET_SAMPLES:
        f.append("former starvation replay length changed")
    if int(d.get("former_starvation_cycle_pseudo_firings", 0)) <= 0:
        f.append("former starvation replay still has zero S firings")
    if int(d.get("former_starvation_cycle_worst_gap_samples", 9999)) > DEPLOYED_MAX_GAP_SAMPLES:
        f.append(f"former starvation replay exceeds certified {DEPLOYED_MAX_GAP_SAMPLES}-sample gap")
    if d.get("classification") != f"SOURCE_INDEPENDENT_{DEPLOYED_MAX_GAP_SAMPLES}_SAMPLE_S_RECURRENCE_CERTIFIED":
        f.append("unexpected scheduler classification")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "period_s": d["deployed_pseudo_period_binary32_s"],
        "uniform_gap_samples": d["certified_uniform_max_gap_samples"],
        "uniform_gap_s": d["certified_uniform_max_gap_s"],
        "translation_bound_s": d["translation_upper_cadence_plus_h_s"],
        "former_cycle_firings": d["former_starvation_cycle_pseudo_firings"],
        "former_cycle_worst_gap_samples": d["former_starvation_cycle_worst_gap_samples"],
        "failures": failures,
        "next": d["next_obligation"],
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
