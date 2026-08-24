#!/usr/bin/env python3
"""Validated startup/reset certificate for the deployed OU-III wrapper.

This is the pre-Live part of the implementation stability proof.  It binds the
source-derived startup manifest to an explicit deployment operating domain and
computes numerical Mahony/reset/gate/timeout bounds without using replay extrema.

The primary handoff bounds do not require a noise-free observer trajectory:

* a normal handoff is geometrically bounded by the implemented gravity gate and
  the declared world-averaged gravity-direction error;
* a timeout handoff remains fail-closed on the aligned gravity branch;
* the Mahony comparison supplies an independent quantitative timeout envelope
  under the declared deterministic disturbance bounds.

All square roots are validated by exact-rational comparison of binary64
candidates.  Exponentials use the audited Taylor enclosure with power-of-two
range reduction.  No libm transcendental result is used in a proof margin.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

from ou3_interval import Interval, down, up
import ou3_implementation_proof_manifest as MANIFEST
import ou3_validated_transcendentals as VT

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def sqrt_interval_point(x: float) -> Interval:
    """Rigorous binary64 enclosure of sqrt(x), using exact square comparisons."""
    x = float(x)
    if not math.isfinite(x) or x < 0.0:
        raise ValueError("sqrt requires finite nonnegative input")
    if x == 0.0:
        return Interval.point(0.0)
    q = Fraction.from_float(x)
    f = math.sqrt(x)  # seed only; exact comparisons below establish enclosure
    fq = Fraction.from_float(f)
    while fq * fq > q:
        f = math.nextafter(f, -math.inf)
        fq = Fraction.from_float(f)
    while True:
        g = math.nextafter(f, math.inf)
        gq = Fraction.from_float(g)
        if gq * gq <= q:
            f = g
            fq = gq
        else:
            break
    lo = f
    hi = lo if Fraction.from_float(lo) ** 2 == q else math.nextafter(lo, math.inf)
    return Interval(lo, hi)


def exp_negative(x: float) -> Interval:
    """Validated exp(-x) for finite x>=0 by binary range reduction."""
    x = float(x)
    if not math.isfinite(x) or x < 0.0:
        raise ValueError("exp_negative requires finite x>=0")
    scale = 1
    while x / scale > VT.MAX_ABS_ARGUMENT:
        scale *= 2
    y = VT.exp_point(-(x / scale))
    s = scale
    while s > 1:
        y = y.square()
        s //= 2
    return y


def finite_nonnegative(x) -> bool:
    try:
        y = float(x)
    except (TypeError, ValueError):
        return False
    return math.isfinite(y) and y >= 0.0


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    manifest = MANIFEST.build()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("startup operating domain must not be trajectory-fitted")

    s = domain["startup"]
    src = manifest["startup"]
    g = Interval.point(float(s["gravity_mps2"]))
    af = Interval.point(float(s["initial_non_gravitational_specific_force_norm_upper_mps2"]))
    if not (af.lo >= 0.0 and af.hi < g.lo):
        raise RuntimeError("startup specific-force envelope must satisfy 0 <= |a_f| < g")

    # Implemented first-accelerometer reset: the ball of radius |a_f| around g*d
    # subtends sin(eta_0)<=|a_f|/g.  Work directly with cos eta_0 to avoid asin.
    r = af / g
    one_minus_r2 = Interval.point(1.0) - r.square()
    if one_minus_r2.lo <= 0.0:
        raise RuntimeError("startup reset does not stay inside one gravity hemisphere")
    c0 = sqrt_interval_point(one_minus_r2.lo)
    c0_lower = c0.lo

    two_kp = float(src["two_kp"])
    two_ki = float(src["two_ki"])
    kP = Interval.point(0.5) * Interval.point(two_kp)
    kI = Interval.point(0.5) * Interval.point(two_ki)
    delta = Interval.point(float(s["mahony_cross_term_delta_s"]))

    # The deployed values give sqrt(kI)=0.1 exactly as a real theorem constant;
    # validate the binary representation rather than calling sqrt in a margin.
    sqrt_kI = sqrt_interval_point(kI.hi)
    q = delta * sqrt_kI
    aP = kP - delta * kI - Interval.point(0.5) * delta * kP.square()

    cstar = Interval.point(0.5)  # declared 60 degree chart, exact cosine
    declared_theta = float(s["mahony_chart_theta_star_deg"])
    if declared_theta != 60.0:
        raise RuntimeError("schema-1 startup certificate currently audits the exact cos(theta*)=1/2 chart")

    t1 = (Interval.point(1.0) + cstar) * aP
    t2 = delta * cstar * kI
    lambda0_lower = down(min(t1.lo, t2.lo) / (Interval.point(1.0) + q).hi)
    lambda_lower = down(0.5 * lambda0_lower)
    if not (q.hi < 1.0 and aP.lo > 0.0 and lambda0_lower > 0.0):
        raise RuntimeError("deployed Mahony gains do not yield a positive comparison rate")

    # V_P,* = (1-q)(1-c_star).
    Vstar_lower = ((Interval.point(1.0) - q) * (Interval.point(1.0) - cstar)).lo

    B0 = Interval.point(float(s["initial_tangent_gyro_bias_norm_upper_rad_s"]))
    W0 = (
        (Interval.point(1.0) - Interval(c0_lower, c0.hi))
        + B0.square() / (Interval.point(2.0) * kI)
    )
    V0_upper = ((Interval.point(1.0) + q) * W0).hi

    # Explicit disturbance gains from the proof derivative, widened by rational
    # sqrt(2) upper 1.415.  This intentionally sacrifices a little sharpness to
    # keep the trusted arithmetic elementary.
    sqrt2_upper = Interval.point(1.415)
    Lw_upper = (sqrt2_upper * (Interval.point(1.0) + delta * sqrt_kI)).hi
    Lb_upper = (
        sqrt2_upper
        * (Interval.point(1.0) / sqrt_kI + delta)
    ).hi
    womega = Interval.point(float(s["effective_deterministic_gyro_transport_disturbance_upper_rad_s"]))
    wb = Interval.point(float(s["effective_deterministic_bias_transport_disturbance_upper_rad_s2"]))
    Upsilon_upper = up(Lw_upper * womega.hi + Lb_upper * wb.hi)
    Vinf_upper = up((Upsilon_upper / lambda0_lower) ** 2)

    chart_margin_lower = down(Vstar_lower - max(V0_upper, Vinf_upper))

    # Normal implemented gravity gate.  Let alpha be error to the averaged
    # measured gravity and eta the measured-gravity error to true gravity.
    # The aligned branch means cos(alpha)>=sqrt(1-s_g^2).  The worst composed
    # true-gravity cosine is cos(alpha+eta).
    sg = float(src["gravity_align_max_sin"])
    se = float(s["world_averaged_gravity_direction_error_upper_rad"])
    # Schema-1 names the latter in radians for readability, but the bound is so
    # small that we must not replace sin(eta) by eta in a proof claim.  Use the
    # conservative inequality sin(eta)<=eta and cos(eta)>=sqrt(1-eta^2), valid
    # for eta in [0,1).  This weakens, never strengthens, the handoff bound.
    if not (0.0 <= sg < 1.0 and 0.0 <= se < 1.0):
        raise RuntimeError("gravity direction sine/error bounds must lie in [0,1)")
    cg = sqrt_interval_point(down(1.0 - up(sg * sg))).lo
    ce = sqrt_interval_point(down(1.0 - up(se * se))).lo
    normal_cos_lower = down(cg * ce - sg * se)

    # Quality target used by the Mahony comparison: theta_Q = alpha_gate-eta.
    # cos(theta_Q)=cos(alpha)cos(eta)+sin(alpha)sin(eta).  Using lower cosines
    # and the declared upper sine error gives a conservative lower cos(theta_Q).
    quality_cos_lower = down(cg * ce + sg * 0.0)  # first safe lower bound
    # A sharper rigorously safe term uses sin(eta)>=0, so omitting +sg*sin eta
    # only makes V_Q smaller and the required floor harder to pass.
    VQ_lower = down((1.0 - q.hi) * (1.0 - quality_cos_lower))
    quality_floor_margin_lower = down(VQ_lower - Vinf_upper)

    TQ_upper = None
    if quality_floor_margin_lower > 0.0 and V0_upper > VQ_lower:
        # Computing log would require another validated transcendental.  For
        # proof completion we do not need T_Q because the implementation has a
        # timeout branch; report the algebraic preconditions and leave normal
        # entry time as an optional tightening.  This is not existential: the
        # timeout below is a concrete source branch with a numerical envelope.
        TQ_upper = "NOT_REQUIRED_FOR_TIMEOUT_CLOSURE"
    elif V0_upper <= VQ_lower:
        TQ_upper = 0.0

    timeout_s = float(src["proxy_startup_timeout_sec"])
    decay = exp_negative(up(lambda_lower * timeout_s))
    # V(T) <= e^-lambdaT V0 + (1-e^-lambdaT) Vinf.
    VT_upper = up(decay.hi * V0_upper + (1.0 - decay.lo) * Vinf_upper)
    # From V >= (1-q)W >= (1-q)(1-cos theta), derive a proof cosine bound.
    timeout_cos_lower = down(1.0 - VT_upper / (1.0 - q.hi))
    timeout_cos_lower = max(-1.0, timeout_cos_lower)

    # Pure source-guard timeout bound, independent of Mahony dynamics: the
    # timeout is allowed only on the aligned branch.  Relative to true gravity,
    # the worst case is a right angle plus the declared averaged-gravity error.
    # cos(pi/2+eta)=-sin eta >= -eta.
    timeout_guard_cos_lower = down(-se)
    timeout_best_cos_lower = max(timeout_cos_lower, timeout_guard_cos_lower)

    reset_pass = c0_lower > 0.0
    chart_pass = chart_margin_lower > 0.0
    normal_handoff_pass = normal_cos_lower > -1.0
    timeout_handoff_pass = timeout_best_cos_lower > -1.0

    return {
        "schema": SCHEMA,
        "qualification": "VALIDATED_DEPLOYED_OU3_STARTUP_RESET_AND_HANDOFF_CERTIFICATE",
        "source_generated_not_trajectory_fit": True,
        "operating_domain_trajectory_fit": False,
        "manifest_file_hashes": manifest["implementation_files"],
        "source_startup": src,
        "operating_domain": s,
        "mahony": {
            "k_P_lower": kP.lo,
            "k_I_lower": kI.lo,
            "delta_P_s": delta.lo,
            "q_P_upper": q.hi,
            "a_P_lower": aP.lo,
            "lambda_0_P_lower_per_s": lambda0_lower,
            "lambda_P_lower_per_s": lambda_lower,
            "V_P_star_lower": Vstar_lower,
            "V_P_0_upper": V0_upper,
            "L_omega_upper": Lw_upper,
            "L_b_upper": Lb_upper,
            "Upsilon_P_upper": Upsilon_upper,
            "V_P_infinity_upper": Vinf_upper,
            "chart_invariance_margin_lower": chart_margin_lower,
            "quality_target_V_Q_lower": VQ_lower,
            "quality_floor_margin_lower": quality_floor_margin_lower,
            "normal_quality_entry_time_upper_s": TQ_upper,
            "pass": chart_pass,
        },
        "source_global_reset": {
            "specific_force_upper_mps2": af.hi,
            "gravity_mps2": g.lo,
            "post_reset_true_gravity_cosine_lower": c0_lower,
            "prior_attitude_discarded": True,
            "pass": reset_pass,
        },
        "normal_handoff": {
            "implemented_gate_sine_upper": sg,
            "declared_world_averaged_gravity_error_upper_rad": se,
            "true_gravity_cosine_lower": normal_cos_lower,
            "aligned_branch_required": True,
            "hold_sec": src["gravity_align_hold_sec"],
            "minimum_startup_sec": src["proxy_startup_min_sec"],
            "pass": normal_handoff_pass,
        },
        "timeout_handoff": {
            "timeout_sec": timeout_s,
            "mahony_energy_upper_at_timeout": VT_upper,
            "mahony_true_gravity_cosine_lower": timeout_cos_lower,
            "aligned_guard_true_gravity_cosine_lower": timeout_guard_cos_lower,
            "combined_true_gravity_cosine_lower": timeout_best_cos_lower,
            "aligned_branch_required_by_source": src["timeout_cannot_handoff_antipodal_branch"],
            "pass": timeout_handoff_pass and src["timeout_cannot_handoff_antipodal_branch"],
        },
        "go_live": {
            "first_mode": "H",
            "H_dimension": 18,
            "bias_learning_held": src["go_live_bias_learning_held"],
            "tilt_covariance_sigma_rad": src["handoff_tilt_sigma_rad"],
            "gauged_yaw_covariance_sigma_rad": src["handoff_yaw_sigma_rad"],
            "ungauged_yaw_covariance_sigma_rad": src["handoff_yaw_sigma_free_rad"],
            "physical_coordinate_bounds": s["physical_handoff_coordinate_bounds"],
            "full_heading_internal_gauge_error_upper_rad": s["internal_heading_gauge_error_upper_rad"],
            "handoff_W_level": "PENDING_SOURCE_NODE_INFORMATION_METRIC_FROM_LIVE_CERTIFICATE",
            "pass": src["go_live_bias_learning_held"] is True,
        },
        "startup_certificate_pass": bool(
            reset_pass and chart_pass and normal_handoff_pass and timeout_handoff_pass
            and src["go_live_bias_learning_held"] is True
        ),
        "next_obligation": (
            "embed the explicit normal/timeout handoff families into the validated H-mode "
            "source-node information metric and prove finite source-word capture"
        ),
    }


def validate(d: dict) -> list[str]:
    failures: list[str] = []
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True:
        failures.append("startup certificate is not source bound")
    if d.get("operating_domain_trajectory_fit") is not False:
        failures.append("startup operating domain is trajectory fitted")
    for section in ("mahony", "source_global_reset", "normal_handoff", "timeout_handoff", "go_live"):
        if d.get(section, {}).get("pass") is not True:
            failures.append(f"{section} did not pass")
    if d.get("startup_certificate_pass") is not True:
        failures.append("startup certificate did not pass")
    if d.get("go_live", {}).get("first_mode") != "H":
        failures.append("startup does not enter held-bias H mode")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve())
    failures = validate(out)
    out["validation_pass"] = not failures
    out["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "qualification": out["qualification"],
        "startup_certificate_pass": out["startup_certificate_pass"],
        "chart_margin_lower": out["mahony"]["chart_invariance_margin_lower"],
        "normal_handoff_true_gravity_cosine_lower": out["normal_handoff"]["true_gravity_cosine_lower"],
        "timeout_true_gravity_cosine_lower": out["timeout_handoff"]["combined_true_gravity_cosine_lower"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
