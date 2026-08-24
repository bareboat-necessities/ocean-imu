#!/usr/bin/env python3
"""Rigorous source-uniform translational UCO/UCC constants for OU-III.

This closes the two translational existence steps that the manuscript previously
left to compactness:

* the four-firing S=0 observation operator is given an explicit positive
  singular-value/information lower bound; and
* the integrated-OU process Gramian on (v,p,S,a_w) is given an explicit positive
  eigenvalue lower bound on every configured 200 Hz Live sample.

The bounds are deliberately conservative.  They use only implementation/source
limits and validated scalar exponential arithmetic; no replay extrema or fitted
trajectory statistics enter the result.

For the process Gramian, let k(t)=exp(F t)G for one [v,p,S,a] axis.  In the
basis [1,t,t^2,exp(-lambda t)], the coefficient transform has determinant
1/(2 lambda^3).  The third divided difference of exp(-lambda t) has magnitude
at least lambda^3 exp(-lambda h)/6, so for four ordered points

    |det[k(t0),...,k(t3)]| >= V(t0,...,t3) exp(-lambda h)/12.

Andreief's identity and four separated subintervals of width h/7 then give

    det Gram >= (2025/144) (h/7)^16 exp(-2 h/tau_min).

The response norms satisfy |a|<=1, |v|<=t, |p|<=t^2/2, |S|<=t^3/6, hence
trace Gram <= h(1+h^2+h^4/4+h^6/36).  For a 4x4 PSD Gramian,
lambda_min >= det/trace^3.  Multiplying by the minimum OU driving intensity
2 sigma_aw^2/tau gives the covariance lower bound.

The S-observation bound implements the determinant/singular-value argument in
w3d-iss-stability.tex-part.  It is aligned to a pseudo-measurement firing; the
source scheduler gives Delta_min>=h and Delta_max<=T_S,max+h.  Four firings fit
in 3 Delta_max.  The effective S-measurement standard deviation is bounded by
the base R_S clamp times sqrt(T0/T_S,min), which covers the Cubic cadence
renormalization as well as the deployed SpectralMSE law.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_source_domain_contract as SOURCE
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
DEFAULT_HEADER = SOURCE.DEFAULT_HEADER
SCHEMA = 1


def _point(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def _pow_nonnegative(x: Interval, n: int) -> Interval:
    if x.lo < 0.0 or n < 0:
        raise ValueError("nonnegative integer power required")
    y = Interval.point(1.0)
    a = x
    k = n
    while k:
        if k & 1:
            y = y * a
        k >>= 1
        if k:
            a = a * a
    return y


def _exp_negative_wide(x: Interval) -> Interval:
    """Enclose exp(-X), X>=0, by power-of-two range reduction.

    The trusted Taylor core is only audited for |argument|<=1/2.  Division by a
    power of two and repeated interval squaring extends that core without a
    libm transcendental call.
    """
    if x.lo < 0.0 or not math.isfinite(x.hi):
        raise ValueError("wide negative exponential requires finite X>=0")
    scale = 1
    while x.hi / scale > VT.MAX_ABS_ARGUMENT:
        scale *= 2
    z = x / Interval.point(float(scale))
    y = VT.exp_interval(-z)
    s = scale
    while s > 1:
        y = y.square()
        s //= 2
    return y


def build(header: Path = DEFAULT_HEADER.resolve()) -> dict:
    header = header.resolve()
    source = SOURCE.build(header)
    box = source["validated_parameter_box"]
    cp = box["continuous_parameters"]
    runtime = source["configured_runtime_assumption"]

    h = Interval(*runtime["imu_dt_outward_interval_s"])
    tau = Interval(*cp["tau_aw_s"])
    sigma = Interval(*cp["sigma_aw_mps2"])
    rs_base = Interval(*cp["R_S_base"])
    ts = Interval(*cp["pseudo_update_period_s"])

    if h.lo <= 0.0 or tau.lo <= 0.0 or sigma.lo <= 0.0 or ts.lo <= 0.0:
        raise RuntimeError("source proof box lost a positive lower endpoint")

    # ---------- Process uniform complete controllability ----------
    # Use the worst decay at h_max/tau_min.
    x_decay = Interval(h.lo, h.hi) / Interval(tau.lo, tau.lo)
    e_decay = _exp_negative_wide(x_decay)
    e2_lower = _pow_nonnegative(Interval(e_decay.lo, e_decay.lo), 2).lo

    w = h / _point(7.0)
    w16 = _pow_nonnegative(w, 16)
    det_gram_lower = (
        _point(2025.0 / 144.0) * w16 * Interval(e2_lower, e2_lower)
    ).lo

    h2 = h.square()
    h4 = h2.square()
    h6 = h4 * h2
    trace_term = (
        _point(1.0)
        + h2
        + h4 / _point(4.0)
        + h6 / _point(36.0)
    )
    trace_gram_upper = (h * trace_term).hi
    gram_lambda_min_lower = math.nextafter(
        det_gram_lower / (trace_gram_upper ** 3), -math.inf
    )
    # q_c = 2 sigma^2/tau; lower endpoint uses sigma_min,tau_max.
    qc_lower = (
        _point(2.0)
        * Interval(sigma.lo, sigma.lo).square()
        / Interval(tau.hi, tau.hi)
    ).lo
    q_axis_lambda_min_lower = math.nextafter(
        qc_lower * gram_lambda_min_lower, -math.inf
    )

    # ---------- S=0 uniform complete observability ----------
    # periodic_update_due can overshoot the requested period by less than one
    # configured sample.  A word aligned to one firing therefore sees four
    # firings within 3*(T_S,max+h_max).
    delta_min = h.lo
    delta_max = math.nextafter(ts.hi + h.hi, math.inf)
    window = math.nextafter(3.0 * delta_max, math.inf)

    # determinant lower d = Delta_min^6/12 * exp(-T/tau_min)
    decay_obs = _exp_negative_wide(
        Interval.outward_bounds(window / tau.lo, window / tau.lo)
    )
    delta6 = _pow_nonnegative(Interval(delta_min, delta_min), 6).lo
    det_obs_lower = math.nextafter(
        (delta6 / 12.0) * decay_obs.lo, -math.inf
    )

    # Each row is [t^2/2,t,1,A3(t)], A3(t)<=t^3/6.  Four rows give the
    # Frobenius upper bound 2*row_norm_max.
    T = Interval.outward_bounds(0.0, window)
    row_norm2 = (
        _point(1.0)
        + T.square()
        + _pow_nonnegative(T, 4) / _point(4.0)
        + _pow_nonnegative(T, 6) / _point(36.0)
    )
    B_frob_upper = math.nextafter(2.0 * math.sqrt(row_norm2.hi), math.inf)
    obs_sigma_min_lower = math.nextafter(
        det_obs_lower / (B_frob_upper ** 3), -math.inf
    )

    # Covers the largest possible cadence normalization of the Cubic law.
    nominal_T0 = SOURCE.parse_const(header.read_text(encoding="utf-8"),
                                    "PSEUDO_UPDATE_PERIOD_NOMINAL_S")
    cadence_scale_upper = math.nextafter(math.sqrt(nominal_T0 / ts.lo), math.inf)
    rs_filter_std_upper = math.nextafter(rs_base.hi * cadence_scale_upper, math.inf)
    s_info_lambda_min_lower = math.nextafter(
        (obs_sigma_min_lower ** 2) / (rs_filter_std_upper ** 2), -math.inf
    )

    process_pass = all(math.isfinite(v) and v > 0.0 for v in (
        det_gram_lower, gram_lambda_min_lower, qc_lower, q_axis_lambda_min_lower
    ))
    observability_pass = all(math.isfinite(v) and v > 0.0 for v in (
        det_obs_lower, obs_sigma_min_lower, s_info_lambda_min_lower
    ))

    return {
        "schema": SCHEMA,
        "qualification": "VALIDATED_TRANSLATIONAL_UCO_UCC_CONFIGURED_RUNTIME",
        "source_generated_not_trajectory_fit": True,
        "validated_arithmetic": True,
        "outward_rounded": True,
        "configured_runtime": runtime,
        "process_ucc": {
            "state_order": ["v", "p", "S", "a_w"],
            "h_s": h.as_list(),
            "tau_s": tau.as_list(),
            "sigma_aw_mps2": sigma.as_list(),
            "exp_minus_h_over_tau_min_lower": e_decay.lo,
            "andreief_subinterval_width_lower_s": w.lo,
            "unit_gramian_det_lower": det_gram_lower,
            "unit_gramian_trace_upper": trace_gram_upper,
            "unit_gramian_lambda_min_lower": gram_lambda_min_lower,
            "ou_driving_intensity_lower": qc_lower,
            "Q_axis_lambda_min_lower": q_axis_lambda_min_lower,
            "pass": process_pass,
        },
        "S_observation_uco": {
            "aligned_firing_count": 4,
            "pseudo_gap_min_s": delta_min,
            "pseudo_gap_max_s": delta_max,
            "aligned_window_s": window,
            "exp_minus_window_over_tau_min_lower": decay_obs.lo,
            "observation_det_lower": det_obs_lower,
            "observation_frobenius_upper": B_frob_upper,
            "observation_sigma_min_lower": obs_sigma_min_lower,
            "R_S_filter_std_upper": rs_filter_std_upper,
            "information_gramian_lambda_min_lower": s_info_lambda_min_lower,
            "pass": observability_pass,
        },
        "translation_source_complete": bool(process_pass and observability_pass),
        "continuous_word_enclosed": False,
        "nonlinear_word_enclosed": False,
        "theorem_promotion": "NOT_ESTABLISHED",
        "next_obligation": (
            "combine these strict translational UCO/UCC constants with a source-uniform "
            "vector-packet attitude/gyro-bias information certificate, then propagate the "
            "complete Riccati/correction word in validated matrix arithmetic"
        ),
    }


def validate(payload: dict) -> list[str]:
    failures: list[str] = []
    if payload.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for flag in ("source_generated_not_trajectory_fit", "validated_arithmetic",
                 "outward_rounded", "translation_source_complete"):
        if payload.get(flag) is not True:
            failures.append(f"{flag} is not true")
    for section, keys in {
        "process_ucc": (
            "unit_gramian_det_lower", "unit_gramian_lambda_min_lower",
            "ou_driving_intensity_lower", "Q_axis_lambda_min_lower",
        ),
        "S_observation_uco": (
            "observation_det_lower", "observation_sigma_min_lower",
            "information_gramian_lambda_min_lower",
        ),
    }.items():
        row = payload.get(section, {})
        if row.get("pass") is not True:
            failures.append(f"{section} did not pass")
        for key in keys:
            value = row.get(key)
            if not isinstance(value, (int, float)) or not math.isfinite(float(value)) or not float(value) > 0.0:
                failures.append(f"{section}.{key} is not finite positive")
    if payload.get("continuous_word_enclosed") is not False:
        failures.append("translation stage must not assert full word enclosure")
    if payload.get("theorem_promotion") != "NOT_ESTABLISHED":
        failures.append("translation stage must not promote theorem")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--header", type=Path, default=DEFAULT_HEADER)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.header.resolve())
    failures = validate(out)
    out["validation_pass"] = not failures
    out["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "qualification": out["qualification"],
        "translation_source_complete": out["translation_source_complete"],
        "Q_axis_lambda_min_lower": out["process_ucc"]["Q_axis_lambda_min_lower"],
        "S_information_lambda_min_lower": out["S_observation_uco"]["information_gramian_lambda_min_lower"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
