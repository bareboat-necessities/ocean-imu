#!/usr/bin/env python3
"""Non-promoting diagnostic for tau-EMA realizability of the P3 timer reset cycle.

The coarse P2 quotient already admits a 635-sample word with no S=0 pseudo
measurement because changing the tau-scaled pseudo period executes

    elapsed <- fmod(elapsed, new_period)

before the next pseudo-update test.  That result alone does not tell us whether
the within-cell tau values can be connected by the shipping tau EMA.

This diagnostic keeps the exact certified 13-sample source clock and transcribes
only the shipping float tau channel:

    tau += alpha * (tau_target - tau)
    alpha = 1 - exp(-dt / clamp(0.40 * clamp(T_sea,0.5,6),0.05,30))

with tau_target = 0.5/f_tune for the configured tau coefficient.  The float
exponential is evaluated through the platform libm ``expf`` entry point rather
than Python's binary64 ``math.exp``.  This matches the precision class of the
C++ ``std::exp(float)`` overload used by the deployed source and avoids making
an exact-image claim from a binary64 transcendental followed by a float cast.

It searches binary32 tuning frequencies whose 13-sample EMA images alternate
between two adjacent applied-tau values around the scheduler reset period.  The
high period is chosen only one scheduler-tolerance margin above the low reset
period, so this is a much sharper test than jumping between distant values
inside one P2 cell.

A positive result means the starvation mechanism is not removed merely by
remembering continuous tau EMA motion inside the current cell.  It is still not
a claim that the full WavePeriodEstimator can synthesize the required frequency
sequence from every admissible physical acceleration record; that upstream
regularity is the next obligation.
"""
from __future__ import annotations

import argparse
import ctypes
import ctypes.util
import json
import math
import re
import struct
from pathlib import Path

import ou3_p4_source_path_reachability as PATH
import ou3_source_reachable_matrix_p3 as BASE
import ou3_p3_pseudo_scheduler_starvation_witness as STARVE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
AUTO_TUNER = REPO / "src" / "tuner" / "SeaStateAutoTuner.h"
SCHEMA = 1
TARGET_SAMPLES = 635
GAP = 13

_LIBM_NAME = ctypes.util.find_library("m")
_LIBM = ctypes.CDLL(_LIBM_NAME) if _LIBM_NAME else ctypes.CDLL(None)
_EXPF = _LIBM.expf
_EXPF.argtypes = [ctypes.c_float]
_EXPF.restype = ctypes.c_float


def _f32(x: float) -> float:
    return struct.unpack("!f", struct.pack("!f", float(x)))[0]


def _bits(x: float) -> int:
    return struct.unpack("!I", struct.pack("!f", _f32(x)))[0]


def _from_bits(i: int) -> float:
    return struct.unpack("!f", struct.pack("!I", int(i)))[0]


def _next_f32(x: float) -> float:
    if not (math.isfinite(x) and x >= 0.0):
        raise ValueError("positive finite float required")
    return _from_bits(_bits(x) + 1)


def _expf(x: float) -> float:
    """Evaluate the float exponential through libm expf."""
    y = float(_EXPF(ctypes.c_float(_f32(x))))
    return _f32(y)


def _literal_member(text: str, name: str) -> float:
    m = re.search(rf"float\s+{re.escape(name)}\s*=\s*([0-9.eE+-]+)f", text)
    if not m:
        raise RuntimeError(f"cannot extract {name}")
    return _f32(float(m.group(1)))


def _effective_frequency_bounds(c: dict) -> tuple[float, float]:
    text = AUTO_TUNER.read_text(encoding="utf-8")
    tuner_lo = _literal_member(text, "f_min_hz")
    tuner_hi = _literal_member(text, "f_max_hz")
    return _f32(max(tuner_lo, c["min_freq"])), _f32(min(tuner_hi, c["max_freq"]))


def _tau_target_from_frequency(freq: float, c: dict) -> float:
    f = _f32(freq)
    lo, hi = _effective_frequency_bounds(c)
    f = _f32(min(max(f, lo), hi))
    raw = _f32(_f32(c["tau_coeff"]) * _f32(_f32(0.5) / f))
    return _f32(min(max(raw, _f32(c["min_tau"])), _f32(c["max_tau"])))


def _tau_ema_step(tau: float, tau_target: float, c: dict) -> tuple[float, float, float]:
    """Shipping float tau update; returns (tau_next, alpha, horizon)."""
    x = _f32(tau)
    target = _f32(tau_target)
    if _f32(c["tau_coeff"]) != _f32(1.0):
        raise RuntimeError("diagnostic currently pins deployed tau_coeff=1")

    # sea_time_sec = 0.5/f_tune = tau_target for tau_coeff=1.
    safe = _f32(min(max(target, _f32(c["time_scale_min"])), _f32(c["time_scale_max"])))
    horizon = _f32(_f32(c["adapt_tau_sea_periods"]) * safe)
    horizon = _f32(min(max(horizon, _f32(c["horizon_min"])), _f32(c["horizon_max"])))
    dt = _f32(c["dt"])
    arg = _f32(-_f32(dt / horizon))
    alpha = _f32(_f32(1.0) - _expf(arg))
    nxt = _f32(x + _f32(alpha * _f32(target - x)))
    return nxt, alpha, horizon


def _tau_ema_samples(tau: float, tau_target: float, c: dict, samples: int = GAP) -> float:
    x = _f32(tau)
    for _ in range(int(samples)):
        x, _, _ = _tau_ema_step(x, tau_target, c)
    return x


def _find_exact_target(tau_start: float, tau_end: float, c: dict) -> float:
    """Find a positive binary32 tau target whose GAP-step image is tau_end."""
    lo = _f32(max(c["min_tau"], c["tau_coeff"] * 0.5 / _effective_frequency_bounds(c)[1]))
    hi = _f32(min(c["max_tau"], c["tau_coeff"] * 0.5 / _effective_frequency_bounds(c)[0]))
    ilo, ihi = _bits(lo), _bits(hi)
    desired = _f32(tau_end)
    if not ilo < ihi:
        raise RuntimeError("empty legal tau-target interval")

    # The scalar EMA image is monotone in the constant positive target on this
    # configured branch. Find the first float whose image is >= desired.
    while ilo + 1 < ihi:
        mid = (ilo + ihi) // 2
        target = _from_bits(mid)
        image = _tau_ema_samples(tau_start, target, c)
        if image < desired:
            ilo = mid
        else:
            ihi = mid

    for i in range(max(_bits(lo), ihi - 8), min(_bits(hi), ihi + 8) + 1):
        target = _from_bits(i)
        if _tau_ema_samples(tau_start, target, c) == desired:
            return target
    raise RuntimeError(f"no exact binary32 tau target maps {tau_start} -> {tau_end} in {GAP} samples")


def _find_frequency_for_target(tau_target: float, c: dict) -> float:
    lo, hi = _effective_frequency_bounds(c)
    approx = _f32(_f32(0.5) / _f32(tau_target))
    ib = _bits(approx)
    ilo, ihi = _bits(lo), _bits(hi)
    for radius in range(0, 65):
        for i in (ib - radius, ib + radius):
            if i < ilo or i > ihi:
                continue
            f = _from_bits(i)
            if _tau_target_from_frequency(f, c) == _f32(tau_target):
                return f
    raise RuntimeError("could not realize tau target through a legal binary32 tuning frequency")


def _scheduler_tolerance(period: float) -> float:
    eps = _f32(2.0 ** -23)
    return _f32(_f32(16.0 * eps) * _f32(max(1.0, _f32(period))))


def _minimal_high_tau(low_period: float, ratio: float, pmin: float, pmax: float) -> tuple[float, float]:
    threshold = _f32(_f32(low_period) + _scheduler_tolerance(low_period))
    tau = _f32(_f32(low_period) / _f32(ratio))
    for _ in range(200000):
        period = STARVE._period_from_tau(tau, ratio, pmin, pmax)
        if period > threshold:
            return tau, period
        tau = _next_f32(tau)
    raise RuntimeError("could not find adjacent high tau above scheduler tolerance")


def _simulate(c: dict, sched: dict, tau_low: float, tau_high: float,
              target_up: float, target_down: float):
    h = _f32(c["dt"])
    ratio = _f32(sched["pseudo_ratio"])
    pmin = _f32(sched["pseudo_min_s"])
    pmax = _f32(sched["pseudo_max_s"])
    low = STARVE._period_from_tau(tau_low, ratio, pmin, pmax)
    high = STARVE._period_from_tau(tau_high, ratio, pmin, pmax)

    # Start immediately after a pseudo firing with the high applied tau. The
    # first H segment holds tau_high; the second steers to tau_low. Thereafter
    # L steers up and H steers down, so every L commit sees exactly the timer
    # value it reduces to zero.
    elapsed = _f32(0.0)
    wrapper_tau = _f32(tau_high)
    samples = 0
    segment = 0
    fires: list[int] = []
    trace = []
    while samples < TARGET_SAMPLES:
        if segment < 2:
            committed_tau = _f32(tau_high)
            target = _f32(tau_high if segment == 0 else target_down)
            expected_next = _f32(tau_high if segment == 0 else tau_low)
            tag = "H"
        elif segment % 2 == 0:
            committed_tau = _f32(tau_low)
            target = _f32(target_up)
            expected_next = _f32(tau_high)
            tag = "L"
        else:
            committed_tau = _f32(tau_high)
            target = _f32(target_down)
            expected_next = _f32(tau_low)
            tag = "H"

        if wrapper_tau != committed_tau:
            raise RuntimeError(f"tau boundary mismatch at segment {segment}: {wrapper_tau} != {committed_tau}")
        period = STARVE._period_from_tau(committed_tau, ratio, pmin, pmax)
        before = elapsed
        elapsed = STARVE._set_period(elapsed, period)
        after_setter = elapsed
        used = min(GAP, TARGET_SAMPLES - samples)
        for _ in range(used):
            fire, elapsed = STARVE._due(h, period, elapsed)
            wrapper_tau, _, _ = _tau_ema_step(wrapper_tau, target, c)
            samples += 1
            if fire:
                fires.append(samples)
        if used == GAP and wrapper_tau != expected_next:
            raise RuntimeError(
                f"tau EMA lost exact cycle at segment {segment}: {wrapper_tau} != {expected_next}"
            )
        trace.append({
            "segment": segment,
            "period_tag": tag,
            "committed_tau_s": committed_tau,
            "period_s": period,
            "tau_target_s": target,
            "elapsed_before_setter_s": before,
            "elapsed_after_fmod_s": after_setter,
            "elapsed_after_segment_s": elapsed,
            "tau_after_segment_s": wrapper_tau,
            "samples": used,
        })
        segment += 1
    return fires, trace, elapsed, wrapper_tau, low, high


def build(domain_path: Path = DEFAULT_DOMAIN, scheduler_witness: Path | None = None) -> dict:
    path = Path(domain_path).resolve()
    c = PATH._constants()
    sched = BASE.source_schedule()
    h = _f32(c["dt"])
    ratio = _f32(sched["pseudo_ratio"])
    pmin = _f32(sched["pseudo_min_s"])
    pmax = _f32(sched["pseudo_max_s"])

    if scheduler_witness is None:
        coarse = STARVE.build(path)
        coarse_fail = STARVE.validate(coarse)
        if coarse_fail:
            raise RuntimeError(f"coarse scheduler witness prerequisite failed: {coarse_fail}")
    else:
        coarse = json.loads(Path(scheduler_witness).read_text(encoding="utf-8"))
        if not coarse.get("validation_pass", False):
            raise RuntimeError("scheduler witness artifact did not validate")

    if int(coarse["exact_gap_samples"]) != GAP or int(coarse["pseudo_firings_in_target_word"]) != 0:
        raise RuntimeError("scheduler witness no longer supplies the gap-13 no-fire prerequisite")

    low = _f32(coarse["period_low_binary32_s"])
    tau_low = _f32(coarse["tau_low_binary32_s"])
    if STARVE._period_from_tau(tau_low, ratio, pmin, pmax) != low:
        raise RuntimeError("coarse witness low period no longer round-trips")

    tau_high, high = _minimal_high_tau(low, ratio, pmin, pmax)
    target_up = _find_exact_target(tau_low, tau_high, c)
    target_down = _find_exact_target(tau_high, tau_low, c)
    freq_up = _find_frequency_for_target(target_up, c)
    freq_down = _find_frequency_for_target(target_down, c)
    freq_hold = _find_frequency_for_target(tau_high, c)

    fires, trace, endpoint, tau_endpoint, low2, high2 = _simulate(
        c, sched, tau_low, tau_high, target_up, target_down
    )
    if low2 != low or high2 != high:
        raise RuntimeError("scheduler period reconstruction drifted")

    tau_cell = tuple(map(float, coarse["source_tau_cell_s"]))
    period_cell = tuple(map(float, coarse["source_pseudo_period_cell_s"]))
    same_cell = (
        tau_cell[0] <= tau_low <= tau_cell[1]
        and tau_cell[0] <= tau_high <= tau_cell[1]
        and period_cell[0] <= low <= period_cell[1]
        and period_cell[0] <= high <= period_cell[1]
    )
    flo, fhi = _effective_frequency_bounds(c)
    legal_freqs = all(flo <= f <= fhi for f in (freq_up, freq_down, freq_hold))

    _, alpha_up, horizon_up = _tau_ema_step(tau_low, target_up, c)
    _, alpha_down, horizon_down = _tau_ema_step(tau_high, target_down, c)
    rounded64_up = _f32(_f32(1.0) - _f32(math.exp(_f32(-_f32(h / horizon_up)))))
    rounded64_down = _f32(_f32(1.0) - _f32(math.exp(_f32(-_f32(h / horizon_down)))))

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_TAU_EMA_PSEUDO_SCHEDULER_CYCLE_DIAGNOSTIC",
        "diagnostic_only": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "scheduler_witness_consumed": True,
        "scheduler_source_node": int(coarse["source_node"]),
        "exact_stage_gap_samples": GAP,
        "exact_gap_self_edge": bool(coarse["exact_gap_self_edge"]),
        "dt_binary32_s": h,
        "tau_cell_s": list(tau_cell),
        "period_cell_s": list(period_cell),
        "tau_low_binary32_s": tau_low,
        "tau_high_binary32_s": tau_high,
        "tau_separation_s": _f32(tau_high - tau_low),
        "period_low_binary32_s": low,
        "period_high_binary32_s": high,
        "period_separation_s": _f32(high - low),
        "scheduler_tolerance_s": _scheduler_tolerance(low),
        "high_period_exceeds_low_plus_scheduler_tolerance": high > _f32(low + _scheduler_tolerance(low)),
        "both_tau_period_values_inside_same_P2_cell": same_cell,
        "tau_target_up_binary32_s": target_up,
        "tau_target_down_binary32_s": target_down,
        "frequency_for_up_target_hz": freq_up,
        "frequency_for_down_target_hz": freq_down,
        "frequency_for_hold_target_hz": freq_hold,
        "effective_tuning_frequency_bounds_hz": [flo, fhi],
        "all_required_frequencies_inside_shipping_bounds": legal_freqs,
        "frequency_roundtrip_to_up_target_exact": _tau_target_from_frequency(freq_up, c) == target_up,
        "frequency_roundtrip_to_down_target_exact": _tau_target_from_frequency(freq_down, c) == target_down,
        "frequency_roundtrip_to_hold_target_exact": _tau_target_from_frequency(freq_hold, c) == tau_high,
        "gap_step_tau_low_to_high_exact": _tau_ema_samples(tau_low, target_up, c) == tau_high,
        "gap_step_tau_high_to_low_exact": _tau_ema_samples(tau_high, target_down, c) == tau_low,
        "gap_step_tau_high_hold_exact": _tau_ema_samples(tau_high, tau_high, c) == tau_high,
        "tau_ema_dynamic_horizon_transcribed": True,
        "tau_ema_float_exp_backend": "libm_expf",
        "tau_ema_float_exp_backend_is_float": True,
        "python_binary64_exp_used_for_cycle_alpha": False,
        "alpha_up_binary32": alpha_up,
        "alpha_down_binary32": alpha_down,
        "rounded_binary64_alpha_up_binary32": rounded64_up,
        "rounded_binary64_alpha_down_binary32": rounded64_down,
        "pseudo_period_fmod_and_due_transcribed": True,
        "target_samples": TARGET_SAMPLES,
        "pseudo_firings_in_target_word": len(fires),
        "pseudo_firing_sample_indices": fires,
        "endpoint_elapsed_s": endpoint,
        "endpoint_tau_s": tau_endpoint,
        "boundary_trace": trace,
        "continuous_tau_channel_cycle_found": bool(not fires and same_cell and legal_freqs),
        "full_WavePeriodEstimator_realizability_proved": False,
        "deployed_filter_starvation_claimed": False,
        "interval_certificate": False,
        "P3_PROMOTED": False,
        "P4_PROMOTED": False,
        "P5_PROMOTED": False,
        "classification": (
            "TAU_EMA_INPUT_FAMILY_ADMITS_635_SAMPLE_NO_S_FIRING_CYCLE"
            if not fires and same_cell and legal_freqs
            else "TAU_EMA_CYCLE_NOT_FOUND"
        ),
        "next_obligation": (
            "the starvation mechanism survives the shipping float tau EMA under legal tuning-frequency inputs; "
            "do not repair P3 by tau-cell subdivision alone. Either derive an upstream WavePeriodEstimator/source-regularity theorem "
            "that excludes this frequency sequence, or formulate P3 with explicit scheduler/tuner information recurrence."
        ),
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA or d.get("qualification") != "OU3_P3_TAU_EMA_PSEUDO_SCHEDULER_CYCLE_DIAGNOSTIC":
        f.append("schema/qualification mismatch")
    for key in (
        "diagnostic_only", "scheduler_witness_consumed", "exact_gap_self_edge",
        "high_period_exceeds_low_plus_scheduler_tolerance",
        "both_tau_period_values_inside_same_P2_cell",
        "all_required_frequencies_inside_shipping_bounds",
        "frequency_roundtrip_to_up_target_exact", "frequency_roundtrip_to_down_target_exact",
        "frequency_roundtrip_to_hold_target_exact",
        "gap_step_tau_low_to_high_exact", "gap_step_tau_high_to_low_exact",
        "gap_step_tau_high_hold_exact", "tau_ema_dynamic_horizon_transcribed",
        "tau_ema_float_exp_backend_is_float",
        "pseudo_period_fmod_and_due_transcribed", "continuous_tau_channel_cycle_found",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_changed",
        "python_binary64_exp_used_for_cycle_alpha",
        "full_WavePeriodEstimator_realizability_proved", "deployed_filter_starvation_claimed",
        "interval_certificate", "P3_PROMOTED", "P4_PROMOTED", "P5_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("tau_ema_float_exp_backend") != "libm_expf":
        f.append("tau EMA float exponential backend changed")
    if int(d.get("exact_stage_gap_samples", 0)) != GAP:
        f.append("stage gap changed")
    if int(d.get("target_samples", 0)) != TARGET_SAMPLES:
        f.append("target word changed")
    if int(d.get("pseudo_firings_in_target_word", -1)) != 0:
        f.append("tau-EMA cycle contains an S firing")
    if d.get("classification") != "TAU_EMA_INPUT_FAMILY_ADMITS_635_SAMPLE_NO_S_FIRING_CYCLE":
        f.append("tau-EMA no-fire classification not established")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--scheduler-witness", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain, a.scheduler_witness)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "classification": d["classification"],
        "tau_low_high_s": [d["tau_low_binary32_s"], d["tau_high_binary32_s"]],
        "period_low_high_s": [d["period_low_binary32_s"], d["period_high_binary32_s"]],
        "tau_targets_s": [d["tau_target_down_binary32_s"], d["tau_target_up_binary32_s"]],
        "target_frequencies_hz": [d["frequency_for_down_target_hz"], d["frequency_for_up_target_hz"]],
        "alpha_binary32": [d["alpha_down_binary32"], d["alpha_up_binary32"]],
        "rounded_binary64_alpha_binary32": [d["rounded_binary64_alpha_down_binary32"], d["rounded_binary64_alpha_up_binary32"]],
        "pseudo_firings": d["pseudo_firings_in_target_word"],
        "validation_failures": vf,
        "next": d["next_obligation"],
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
