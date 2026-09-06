#!/usr/bin/env python3
"""Validated SEA3-owned OU-III tuner/scheduler state transition.

This module advances the shipping front-end/adaptation state of one already
admitted SEA3 realization.  It is not a source generator: ``a_vertical`` and
``f_wave_previous_wpe`` must come from the same SEA3 private-Mahony/WPE path.
The current measurement can move the candidate only; a staged candidate is
committed at the beginning of the following physical sample.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path
import re

from ou3_interval import Interval
import ou3_source_domain_contract as SOURCE
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[2]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
TUNER = REPO / "src" / "tuner" / "SeaStateAutoTuner.h"
BAND = REPO / "src" / "tuner" / "AdaptiveWaveBandPass.h"
LIMITS = REPO / "src" / "tuner" / "SeaStateAdaptationLimits.h"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_TUNER_SCHEDULER_INTERVAL_STEP"

PI = Interval.outward_bounds(3.141592653589793, 3.141592653589794)
ONE = Interval.point(1.0)
TWO = Interval.point(2.0)


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def ipow(x: Interval, n: int) -> Interval:
    if n < 0:
        raise ValueError("nonnegative integer power required")
    out = Interval.point(1.0)
    for _ in range(n):
        out = out * x
    return out


def clamp_interval(x: Interval, lo: float, hi: float) -> Interval:
    if not (math.isfinite(lo) and math.isfinite(hi) and lo <= hi):
        raise ValueError("finite ordered clamp required")
    return Interval(max(lo, min(hi, x.lo)), max(lo, min(hi, x.hi)))


def max_interval(x: Interval, c: float) -> Interval:
    return Interval(max(x.lo, c), max(x.hi, c))


def min_interval(x: Interval, c: float) -> Interval:
    return Interval(min(x.lo, c), min(x.hi, c))


def _wide_exp(x: Interval) -> Interval:
    if not (math.isfinite(x.lo) and math.isfinite(x.hi)):
        raise ValueError("finite exponential interval required")
    scale = 1
    magnitude = max(abs(x.lo), abs(x.hi))
    while magnitude / scale > VT.MAX_ABS_ARGUMENT:
        scale *= 2
    out = VT.exp_interval(x / I(float(scale)))
    while scale > 1:
        out = out.square()
        scale //= 2
    return out


def _point_power(x: float, n: int) -> Interval:
    return ipow(Interval.point(x), n)


def rational_power_positive(x: Interval, p: int, q: int) -> Interval:
    """Validate x**(p/q); libm is used only to seed an ulp search."""
    if x.lo <= 0.0 or p <= 0 or q <= 0:
        raise ValueError("positive interval and exponents required")
    xp = ipow(x, p)
    lo = math.pow(x.lo, p / q)
    for _ in range(4096):
        if _point_power(lo, q).hi <= xp.lo:
            break
        lo = math.nextafter(lo, -math.inf)
    else:
        raise RuntimeError("could not certify rational-power lower endpoint")
    hi = math.pow(x.hi, p / q)
    for _ in range(4096):
        if _point_power(hi, q).lo >= xp.hi:
            break
        hi = math.nextafter(hi, math.inf)
    else:
        raise RuntimeError("could not certify rational-power upper endpoint")
    return Interval(lo, hi)


def sqrt_interval(x: Interval) -> Interval:
    return rational_power_positive(x, 1, 2)


@dataclass(frozen=True)
class BandState:
    lowpass_low: Interval
    band: Interval
    p00: Interval
    p01: Interval
    p11: Interval
    ready: bool = True


@dataclass(frozen=True)
class MomentState:
    mean_value: Interval
    mean_weight: Interval
    sq_value: Interval
    sq_weight: Interval
    frequency_hz: Interval


@dataclass(frozen=True)
class CandidateState:
    tau: Interval
    sigma: Interval
    rs: Interval


@dataclass(frozen=True)
class ActiveSchedule:
    tau: Interval
    sigma: Interval
    rs_base: Interval
    pseudo_period: Interval


@dataclass(frozen=True)
class SchedulerState:
    since_last_commit_stage_s: Interval
    pending_commit: bool


@dataclass(frozen=True)
class TunerState:
    band: BandState
    moments: MomentState
    candidate: CandidateState
    active: ActiveSchedule
    scheduler: SchedulerState


@dataclass(frozen=True)
class Constants:
    dt: float
    tune_freq_min: float
    tune_freq_max: float
    tau_min: float
    tau_max: float
    sigma_max: float
    rs_min: float
    rs_max: float
    tau_coeff: float
    sigma_coeff: float
    adapt_tau_sea_periods: float
    adapt_rs_mult: float
    adapt_every_s: float
    band_low_ratio: float
    band_high_ratio: float
    band_min_hz: float
    band_max_hz: float
    acc_noise_floor_sigma: float
    pseudo_ratio: float
    pseudo_min_s: float
    pseudo_max_s: float
    rs_mse_coeff: float
    rs_accel_noise_density: float
    rs_x_factor: float
    rs_y_factor: float
    dynamic_scale_min_s: float
    dynamic_scale_max_s: float
    dynamic_horizon_min_s: float
    dynamic_horizon_max_s: float


def _member(text: str, name: str) -> float:
    m = re.search(rf"\b{re.escape(name)}\s*=\s*([0-9.eE+-]+)f\b", text)
    if not m:
        raise RuntimeError(f"cannot extract shipping member {name}")
    return float(m.group(1))


def constants() -> Constants:
    w = WRAPPER.read_text(encoding="utf-8")
    l = LIMITS.read_text(encoding="utf-8")
    names = (
        "FREQ_SMOOTHER_DT", "MIN_TUNE_FREQ_HZ", "MAX_TUNE_FREQ_HZ",
        "MIN_TAU_S", "MAX_TAU_S", "MAX_SIGMA_A", "MIN_R_S", "MAX_R_S",
        "ADAPT_TAU_SEA_PERIODS", "ADAPT_RS_MULT", "ADAPT_EVERY_SECS",
        "SIGMA_BAND_LOW_RATIO_DEFAULT", "SIGMA_BAND_HIGH_RATIO_DEFAULT",
        "SIGMA_BAND_MIN_HZ_DEFAULT", "SIGMA_BAND_MAX_HZ_DEFAULT",
        "ACC_NOISE_FLOOR_SIGMA_DEFAULT", "PSEUDO_UPDATE_TAU_RATIO_DEFAULT",
        "PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT", "PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT",
        "R_S_MSE_COEFF_DEFAULT", "R_S_ACCEL_NOISE_DENSITY_DEFAULT",
    )
    c = {name: float(SOURCE.parse_const(w, name)) for name in names}
    return Constants(
        dt=c["FREQ_SMOOTHER_DT"],
        tune_freq_min=c["MIN_TUNE_FREQ_HZ"], tune_freq_max=c["MAX_TUNE_FREQ_HZ"],
        tau_min=c["MIN_TAU_S"], tau_max=c["MAX_TAU_S"], sigma_max=c["MAX_SIGMA_A"],
        rs_min=c["MIN_R_S"], rs_max=c["MAX_R_S"],
        tau_coeff=_member(w, "tau_coeff_"), sigma_coeff=_member(w, "sigma_coeff_"),
        adapt_tau_sea_periods=c["ADAPT_TAU_SEA_PERIODS"],
        adapt_rs_mult=c["ADAPT_RS_MULT"], adapt_every_s=c["ADAPT_EVERY_SECS"],
        band_low_ratio=c["SIGMA_BAND_LOW_RATIO_DEFAULT"],
        band_high_ratio=c["SIGMA_BAND_HIGH_RATIO_DEFAULT"],
        band_min_hz=c["SIGMA_BAND_MIN_HZ_DEFAULT"], band_max_hz=c["SIGMA_BAND_MAX_HZ_DEFAULT"],
        acc_noise_floor_sigma=c["ACC_NOISE_FLOOR_SIGMA_DEFAULT"],
        pseudo_ratio=c["PSEUDO_UPDATE_TAU_RATIO_DEFAULT"],
        pseudo_min_s=c["PSEUDO_UPDATE_PERIOD_MIN_S_DEFAULT"],
        pseudo_max_s=c["PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT"],
        rs_mse_coeff=c["R_S_MSE_COEFF_DEFAULT"],
        rs_accel_noise_density=c["R_S_ACCEL_NOISE_DENSITY_DEFAULT"],
        rs_x_factor=_member(w, "R_S_x_factor_"), rs_y_factor=_member(w, "R_S_y_factor_"),
        dynamic_scale_min_s=float(SOURCE.parse_const(l, "kDynamicEmaTimeScaleMinSec")),
        dynamic_scale_max_s=float(SOURCE.parse_const(l, "kDynamicEmaTimeScaleMaxSec")),
        dynamic_horizon_min_s=float(SOURCE.parse_const(l, "kDynamicEmaHorizonMinSec")),
        dynamic_horizon_max_s=float(SOURCE.parse_const(l, "kDynamicEmaHorizonMaxSec")),
    )


def pseudo_period(tau: Interval, c: Constants) -> Interval:
    return clamp_interval(I(c.pseudo_ratio) * tau, c.pseudo_min_s, c.pseudo_max_s)


def active_rs_std_xyz(active: ActiveSchedule, c: Constants) -> list[Interval]:
    rs = clamp_interval(active.rs_base, c.rs_min, c.rs_max)
    return [I(c.rs_x_factor) * rs, I(c.rs_y_factor) * rs, rs]


def band_step(state: BandState, x: Interval, f_ref: Interval, c: Constants) -> BandState:
    if f_ref.lo <= 0.0:
        raise ValueError("positive SEA3 wave frequency required")
    upper_limit = min(c.band_max_hz, 0.45 / c.dt)
    low = clamp_interval(I(c.band_low_ratio) * f_ref, c.band_min_hz, upper_limit / 1.05)
    high = min_interval(I(c.band_high_ratio) * f_ref, upper_limit)
    high = Interval(max(high.lo, 1.05 * low.lo), max(high.hi, 1.05 * low.hi))
    high = min_interval(high, upper_limit)
    if not high.lo > low.hi:
        raise RuntimeError("band-corner interval lost high>low correlation; split SEA3 frequency cell")

    k = TWO * PI * I(c.dt)
    alpha_l = ONE - VT.exp_interval(-(k * low))
    alpha_h = ONE - VT.exp_interval(-(k * high))
    ql, qh = ONE - alpha_l, ONE - alpha_h
    low_new = ql * state.lowpass_low + alpha_l * x
    band_new = qh * state.band + alpha_h * ql * (x - state.lowpass_low)

    a00, a10, a11 = ql, -(alpha_h * ql), qh
    b0, b1 = alpha_l, alpha_h * ql
    p00 = a00.square() * state.p00 + b0.square()
    p01 = a00 * (a10 * state.p00 + a11 * state.p01) + b0 * b1
    p11 = (
        a10.square() * state.p00 + I(2.0) * a10 * a11 * state.p01
        + a11.square() * state.p11 + b1.square()
    )
    return BandState(low_new, band_new, max_interval(p00, 0.0), p01, max_interval(p11, 0.0), True)


def moment_step(state: MomentState, accel_band: Interval, f_wave: Interval, c: Constants) -> MomentState:
    f = clamp_interval(f_wave, 0.05, 5.0)
    sea_time = clamp_interval(I(0.5) / f, c.dynamic_scale_min_s, c.dynamic_scale_max_s)
    requested = clamp_interval(I(4.0) * (I(2.0) * sea_time), 0.3, 60.0)
    horizon = clamp_interval(requested, max(c.dynamic_horizon_min_s, c.dt), c.dynamic_horizon_max_s)
    alpha = ONE - VT.exp_interval(-(I(c.dt) / horizon))
    oma = ONE - alpha
    return MomentState(
        oma * state.mean_value + alpha * accel_band,
        oma * state.mean_weight + alpha,
        oma * state.sq_value + alpha * accel_band.square(),
        oma * state.sq_weight + alpha,
        f,
    )


def acceleration_variance(state: MomentState) -> Interval:
    if not (state.mean_weight.lo > 1e-6 and state.sq_weight.lo > 1e-6):
        raise RuntimeError("Live tuner moment weights must exceed readiness threshold")
    mean = state.mean_value / state.mean_weight
    return max_interval(state.sq_value / state.sq_weight - mean.square(), 0.0)


def spectral_mse_rs_target(tau: Interval, sigma: Interval, c: Constants) -> Interval:
    ts = pseudo_period(tau, c)
    sigma_ab = max_interval(sigma / I(c.sigma_coeff), 1e-6)
    u67 = rational_power_positive(sigma_ab * tau.square().square(), 6, 7)
    q14 = rational_power_positive(I(2.0 * c.rs_accel_noise_density), 1, 14)
    return I(c.rs_mse_coeff) * q14 * u67 / sqrt_interval(ts)


def targets(state: TunerState, c: Constants) -> CandidateState:
    f = clamp_interval(state.moments.frequency_hz, c.tune_freq_min, c.tune_freq_max)
    tau = clamp_interval(I(c.tau_coeff * 0.5) / f, c.tau_min, c.tau_max)
    noise_sigma = I(c.acc_noise_floor_sigma) * sqrt_interval(max_interval(state.band.p11, 1e-30))
    var_wave = max_interval(acceleration_variance(state.moments) - noise_sigma.square(), 0.0)
    sigma = min_interval(I(c.sigma_coeff) * sqrt_interval(max_interval(var_wave, 1e-6)), c.sigma_max)
    rs = clamp_interval(spectral_mse_rs_target(tau, sigma, c), c.rs_min, c.rs_max)
    return CandidateState(tau, sigma, rs)


def candidate_ema(old: CandidateState, target: CandidateState, f_tune: Interval, c: Constants) -> CandidateState:
    sea_time = clamp_interval(
        I(0.5) / clamp_interval(f_tune, c.tune_freq_min, c.tune_freq_max),
        c.dynamic_scale_min_s, c.dynamic_scale_max_s,
    )
    common_h = clamp_interval(
        I(c.adapt_tau_sea_periods) * sea_time,
        max(c.dynamic_horizon_min_s, c.dt), c.dynamic_horizon_max_s,
    )
    alpha = ONE - VT.exp_interval(-(I(c.dt) / common_h))

    rs_scale = clamp_interval(target.tau, c.dynamic_scale_min_s, c.dynamic_scale_max_s)
    rs_h = clamp_interval(
        I(c.adapt_rs_mult) * rs_scale,
        max(c.dynamic_horizon_min_s, c.dt), c.dynamic_horizon_max_s,
    )
    alpha_rs = ONE - VT.exp_interval(-(I(c.dt) / rs_h))
    return CandidateState(
        old.tau + alpha * (target.tau - old.tau),
        old.sigma + alpha * (target.sigma - old.sigma),
        old.rs + alpha_rs * (target.rs - old.rs),
    )


def commit_if_pending(state: TunerState, c: Constants) -> TunerState:
    if not state.scheduler.pending_commit:
        return state
    active = ActiveSchedule(
        state.candidate.tau, state.candidate.sigma, state.candidate.rs,
        pseudo_period(state.candidate.tau, c),
    )
    return TunerState(
        state.band, state.moments, state.candidate, active,
        SchedulerState(state.scheduler.since_last_commit_stage_s, False),
    )


def advance_after_measurement(
    state: TunerState,
    *,
    a_vertical: Interval,
    f_wave_previous_wpe: Interval,
    c: Constants | None = None,
) -> list[TunerState]:
    """Advance candidate after the current Kalman update, splitting timer boundary."""
    c = c or constants()
    if f_wave_previous_wpe.lo <= 0.0:
        raise ValueError("same-word previous WPE frequency must stay positive")
    band = band_step(state.band, a_vertical, state.moments.frequency_hz, c)
    moments = moment_step(state.moments, band.band, f_wave_previous_wpe, c)
    provisional = TunerState(band, moments, state.candidate, state.active, state.scheduler)
    candidate = candidate_ema(state.candidate, targets(provisional, c), moments.frequency_hz, c)
    elapsed = state.scheduler.since_last_commit_stage_s + I(c.dt)
    base = TunerState(band, moments, candidate, state.active, SchedulerState(elapsed, False))
    if elapsed.lo > c.adapt_every_s:
        return [TunerState(band, moments, candidate, state.active, SchedulerState(I(0.0), True))]
    if elapsed.hi <= c.adapt_every_s:
        return [base]
    return [
        base,
        TunerState(band, moments, candidate, state.active, SchedulerState(I(0.0), True)),
    ]


def build() -> dict:
    c = constants()
    w = WRAPPER.read_text(encoding="utf-8")
    t = TUNER.read_text(encoding="utf-8")
    b = BAND.read_text(encoding="utf-8")
    parity = {
        "band_state_equations": "band_ = q_high * band_prev + alpha_high * high_passed;" in b,
        "band_noise_covariance_state": "p11_new" in b and "whiteNoiseVarianceGain" in b,
        "debiased_moment_state": (
            "A_mean.update(accel, alpha_var);" in t and "A_sq.update(accel * accel, alpha_var);" in t
        ),
        "tau_target_shipping_formula": "float tau_raw = tau_coeff_ * 0.5f / f_tune;" in w,
        "spectral_mse_deployed": bool(re.search(
            r"\bRSAdaptationLaw\s+rs_law_\s*=\s*RSAdaptationLaw::SpectralMSE\s*;", w
        )),
        "spectral_mse_realized_TS": "const float TS = pseudo_update_period_for_(tau);" in w,
        "candidate_each_sample": (
            "tune_.tau_applied   += alpha" in w and "tune_.RS_applied    += alpha_RS" in w
        ),
        "next_sample_commit": (
            "void apply_pending_online_tune_()" in w and "online_tune_apply_pending_ = true;" in w
        ),
        "strict_stage_timer": bool(re.search(
            r"time_\s*-\s*last_adapt_time_sec_\s*>\s*adapt_every_secs_", w
        )),
        "deployed_rs_slew_disabled": (
            bool(re.search(r"\bADAPT_RS_SLEW_LOG\s*=\s*0\.0f\b", w))
            and bool(re.search(r"\badapt_RS_slew_log_\s*=\s*ADAPT_RS_SLEW_LOG\s*;", w))
        ),
        "spectral_mse_no_extra_cadence_scale": bool(re.search(
            r"if\s*\(rs_law_\s*!=\s*RSAdaptationLaw::Cubic\)\s*return\s+1\.0f\s*;", w
        )),
    }

    f = I(0.2)
    active = ActiveSchedule(I(1.1), I(0.5), I(2.0), pseudo_period(I(1.1), c))
    state = TunerState(
        BandState(I(0.0), I(0.1), I(0.2), I(0.0), I(0.1), True),
        MomentState(I(0.0), I(0.5), I(0.1), I(0.5), f),
        CandidateState(I(1.1), I(0.5), I(2.0)),
        active,
        SchedulerState(I(0.05), False),
    )
    successors = advance_after_measurement(state, a_vertical=I(0.2), f_wave_previous_wpe=f, c=c)
    smoke = {
        "successors": len(successors),
        "target_rs_positive": all(x.candidate.rs.lo > 0.0 for x in successors),
        "active_schedule_unchanged_until_next_sample": all(x.active == active for x in successors),
        "actual_rs_std_xyz": [x.as_list() for x in active_rs_std_xyz(active, c)],
    }
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "source_generator": False,
        "trajectory_replay_used": False,
        "independent_frequency_or_acceleration_box_promotable": False,
        "requires_same_SEA3_Mahony_and_WPE_state": True,
        "same_state_generates_tau_sigma_RS_TS": True,
        "candidate_commit_is_next_sample_predictable": True,
        "scheduler_boundary_is_split_not_selected": True,
        "spectral_MSE_actual_RS_retained": True,
        "shipping_source_parity": parity,
        "shipping_source_parity_pass": all(parity.values()),
        "smoke": smoke,
        "tuner_scheduler_step_closed": all(parity.values()) and smoke["target_rs_positive"],
        "Mahony_step_closed_here": False,
        "WavePeriodEstimator_step_closed_here": False,
        "complete_SEA3_family_materialized_here": False,
        "P3_promoted": False,
        "next_obligation": (
            "connect this recurrence to the same validated SEA3 Mahony/WPE and phase-continuous sea realization"
        ),
    }


def validate(d: dict) -> list[str]:
    failures: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        failures.append("schema/qualification mismatch")
    for key in (
        "requires_same_SEA3_Mahony_and_WPE_state", "same_state_generates_tau_sigma_RS_TS",
        "candidate_commit_is_next_sample_predictable", "scheduler_boundary_is_split_not_selected",
        "spectral_MSE_actual_RS_retained", "shipping_source_parity_pass", "tuner_scheduler_step_closed",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in (
        "source_generator", "trajectory_replay_used", "independent_frequency_or_acceleration_box_promotable",
        "Mahony_step_closed_here", "WavePeriodEstimator_step_closed_here",
        "complete_SEA3_family_materialized_here", "P3_promoted",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    if not all(d.get("shipping_source_parity", {}).values()):
        failures.append("shipping source parity failed")
    if d.get("smoke", {}).get("active_schedule_unchanged_until_next_sample") is not True:
        failures.append("current measurement leaked into current active schedule")
    return list(dict.fromkeys(failures))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build()
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "qualification": d["qualification"],
        "parity": d["shipping_source_parity"],
        "parity_pass": d["shipping_source_parity_pass"],
        "smoke": d["smoke"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
