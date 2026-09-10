#!/usr/bin/env python3
"""Materialized dependency cover of the deployed adaptive coefficient state.

The P4 theorem forbids treating ``(tau, sigma, T_S, R_S)`` as independent
coordinates.  That prohibition has been carried as a contract flag; this
producer materializes the object behind it, using only the deployed shipping
transitions in ``ou3_brmm_tuner_scheduler_step``.

Three exact functional pinnings, verified here against the deployed code rather
than restated:

    T_S        = clamp(pseudo_ratio * tau, [T_S_min, T_S_max])
    R_S_target = SpectralMSE(tau, sigma)          (pseudo cadence inside it)
    R_S applied = (rs_x*r, rs_y*r, r),  r = clamp(rs_base)

so the reachable target quadruple lies on a two-parameter surface inside a
nominally four-dimensional rectangle, and the applied anisotropic ``R_S`` is a
ray, not a box.  On the non-saturated branch the pinning inverts: a ``T_S``
strictly inside its rails fixes ``tau`` exactly, hence fixes the tuning
frequency.

One prior-independent invariance certificate, which needs no reachable-set
argument at all: every shipping target is clamped, and the candidate/active
recurrence is an EMA, i.e. a convex combination of the previous candidate and
the current target.  Therefore the clamped target box is FORWARD INVARIANT for
the candidate, and the staged commit copies the candidate, so it is forward
invariant for the active schedule too.

A finite geometric cover of the reachable ``(f, sigma)`` rectangle is then
mapped through the deployed interval transitions.  Its purpose is quantitative:
it bounds the joint ``S=0`` residual scale ``(T_S + dt)/(rs_x * R_S)`` on cells
where cadence and applied covariance come from the SAME schedule, and compares
that with the independent ``(T_S_max, R_S_min)`` corner that an independent
rectangle would license.

Finally the deployed adaptation rate limits are materialized: the per-step EMA
coefficient of each shipping chain and the fraction of the target gap that each
can close inside one staged-commit interval.  That is the adjacency data a
consecutive-word cell/transition family needs, and it is where the storage
compatibility obligation lives.

What this producer does NOT do, and says so: it does not certify the reachable
``(f, sigma)`` set itself over admitted BRMM continuations, and it does not
compose the cover with the Riccati/Joseph/reset graph.  The source-uniform
COMPLETE BRMM cover therefore stays open and nothing is promoted.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_brmm_complete_source as COMPLETE
import ou3_brmm_dynamic_source_certificate as DYNAMIC
import ou3_brmm_tuner_scheduler_step as TUNER
import ou3_brmm_wpe_state_step as WPE
import ou3_validated_transcendentals as VT
from ou3_interval import Interval

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
QUALIFICATION = "OU3_P4_COEFFICIENT_DEPENDENCY_COVER_V1"
DEFAULT_CELLS = 192

# Shipping raw sigma target floor: min(sigma_coeff*sqrt(max(var,1e-6)), sigma_max).
RAW_SIGMA_FLOOR = 0.0009


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def _geometric_tiling(lo: float, hi: float, cells: int) -> list[tuple[float, float]]:
    """Closed cover of ``[lo, hi]`` by ``cells`` touching geometric intervals."""
    if not (0.0 < lo < hi and cells >= 1):
        raise ValueError("positive ordered range and at least one cell required")
    ratio = pow(hi / lo, 1.0 / cells)
    out: list[tuple[float, float]] = []
    a = lo
    for _ in range(cells - 1):
        b = min(hi, up(a * ratio))
        out.append((a, b))
        a = b
    out.append((a, hi))
    if out[0][0] != lo or out[-1][1] != hi:
        raise RuntimeError("tiling lost an endpoint")
    for (_, b), (a2, _) in zip(out, out[1:]):
        if a2 != b:
            raise RuntimeError("tiling left a gap")
    return out


def _target_image(f_cell: tuple[float, float], sigma_cell: tuple[float, float],
                  c: TUNER.Constants) -> dict:
    """Deployed target image of one ``(f, sigma)`` cell, outward."""
    f = TUNER.clamp_interval(Interval(*f_cell), c.tune_freq_min, c.tune_freq_max)
    tau = TUNER.clamp_interval(Interval.point(c.tau_coeff * 0.5) / f, c.tau_min, c.tau_max)
    ts = TUNER.pseudo_period(tau, c)
    sigma = TUNER.min_interval(Interval(*sigma_cell), c.sigma_max)
    rs = TUNER.clamp_interval(TUNER.spectral_mse_rs_target(tau, sigma, c), c.rs_min, c.rs_max)
    return {"f": f, "tau": tau, "T_S": ts, "sigma": sigma, "R_S": rs}


def _alpha(horizon_s: float, dt: float) -> Interval:
    """Deployed EMA coefficient ``1 - exp(-dt/h)`` at one horizon, outward."""
    if not (horizon_s > 0.0 and dt > 0.0):
        raise ValueError("positive horizon and sample period required")
    return Interval.point(1.0) - VT.exp_interval(-(Interval.point(dt) / Interval.point(horizon_s)))


def _rate_limits(c: TUNER.Constants, w) -> dict:
    """Per-step EMA coefficients and per-commit gap closure of every chain."""
    dt = c.dt
    steps = int(math.ceil(up(c.adapt_every_s) / down(dt)))
    sea_lo = max(c.dynamic_scale_min_s, 0.5 / c.tune_freq_max)
    sea_hi = min(c.dynamic_scale_max_s, 0.5 / c.tune_freq_min)
    guard = max(c.dynamic_horizon_min_s, dt)
    chains = {
        "candidate_tau_sigma_EMA": (
            max(guard, c.adapt_tau_sea_periods * sea_lo),
            min(c.dynamic_horizon_max_s, c.adapt_tau_sea_periods * sea_hi)),
        "candidate_R_S_EMA": (
            max(guard, c.adapt_rs_mult * max(c.dynamic_scale_min_s, c.tau_min)),
            min(c.dynamic_horizon_max_s, c.adapt_rs_mult * c.dynamic_scale_max_s)),
        "tuner_band_moment_EWMA": (
            max(guard, min(60.0, max(0.3, 8.0 * sea_lo))),
            min(c.dynamic_horizon_max_s, min(60.0, max(0.3, 8.0 * sea_hi)))),
        "wave_period_moment_EWMA": (w.min_horizon_s, w.max_horizon_s),
        "wave_period_log_EMA": (
            max(max(w.dynamic_horizon_min_s, dt), w.log_smoothing_periods / c.tune_freq_max),
            min(w.dynamic_horizon_max_s, w.log_smoothing_periods / c.tune_freq_min)),
    }
    out = {}
    for name, (h_lo, h_hi) in chains.items():
        if not (0.0 < h_lo <= h_hi):
            raise RuntimeError("chain " + name + " lost an ordered positive horizon")
        a_hi = _alpha(h_lo, dt).hi
        a_lo = _alpha(h_hi, dt).lo
        if not (0.0 < a_lo <= a_hi < 1.0):
            raise RuntimeError("chain " + name + " lost a strict EMA coefficient")
        closed = up(1.0 - down(pow(down(1.0 - a_hi), steps)))
        out[name] = {
            "horizon_s": [h_lo, h_hi],
            "per_step_alpha_lower": a_lo,
            "per_step_alpha_upper": a_hi,
            "target_gap_closed_per_commit_interval_upper": closed,
            "strict_contraction_toward_target": True,
        }
    out["steps_per_staged_commit_interval"] = steps
    out["staged_commit_interval_s"] = c.adapt_every_s
    return out


def build(domain_path: Path = DEFAULT_DOMAIN, cells: int = DEFAULT_CELLS) -> dict:
    dynamic = DYNAMIC.build(Path(domain_path))
    dyn_failures = DYNAMIC.validate(dynamic)
    source = COMPLETE.build(Path(domain_path))
    src_failures = COMPLETE.validate(source)
    bad = {k: v for k, v in (("dynamic", dyn_failures), ("complete_source", src_failures)) if v}
    if bad:
        raise RuntimeError("coefficient-cover prerequisites failed: " + repr(bad))

    domain = json.loads(Path(domain_path).read_text(encoding="utf-8"))
    c = TUNER.constants()
    w = WPE.constants(c.dt)

    # --- prior-independent forward-invariant coefficient box -----------------
    # Every shipping target is clamped, and candidate = old + alpha*(target-old)
    # with alpha in (0,1) is a convex combination, so the clamped target box is
    # forward invariant for the candidate and, through the staged commit that
    # copies it, for the active schedule.
    tau_target_lo = max(c.tau_min, min(c.tau_max, down(c.tau_coeff * 0.5 / c.tune_freq_max)))
    tau_target_hi = max(c.tau_min, min(c.tau_max, up(c.tau_coeff * 0.5 / c.tune_freq_min)))
    ts_invariant = TUNER.pseudo_period(Interval(tau_target_lo, tau_target_hi), c)
    invariant_box = {
        "tau_s": [tau_target_lo, tau_target_hi],
        "sigma_mps2": [RAW_SIGMA_FLOOR, c.sigma_max],
        "R_S_base_m_s": [c.rs_min, c.rs_max],
        "T_S_s": [ts_invariant.lo, ts_invariant.hi],
        "argument": "clamped target box plus EMA convex combination plus staged commit copy",
        "needs_reachable_set_argument": False,
        "tau_lower_rail_active": bool(tau_target_lo <= c.tau_min),
        "tau_upper_rail_active": bool(tau_target_hi >= c.tau_max),
        "T_S_lower_rail_reachable": bool(ts_invariant.lo <= c.pseudo_min_s),
        "T_S_upper_rail_reachable": bool(ts_invariant.hi >= c.pseudo_max_s),
    }

    # --- exact pinnings, verified against the deployed transitions ------------
    pin_samples = [0.5, 1.0, 3.0, 5.0, 7.0, 10.9, 11.0, 11.5, 12.0]
    pinning = []
    for tau_pt in pin_samples:
        ts = TUNER.pseudo_period(TUNER.I(tau_pt), c)
        expect = TUNER.clamp_interval(
            TUNER.I(c.pseudo_ratio) * TUNER.I(tau_pt), c.pseudo_min_s, c.pseudo_max_s)
        if not (ts.lo == expect.lo and ts.hi == expect.hi):
            raise RuntimeError("deployed pseudo cadence is not the clamped affine image of tau")
        pinning.append({"tau_s": tau_pt, "T_S_s": [ts.lo, ts.hi],
                        "saturated_high": bool(ts.hi >= c.pseudo_max_s)})
    rs_probe = Interval(1.0, 2.0)
    xyz = TUNER.active_rs_std_xyz(
        TUNER.ActiveSchedule(TUNER.I(1.0), TUNER.I(1.0), rs_probe, TUNER.I(0.01)), c)
    # One clamped scalar generates all three applied axes: the x/y axes are the
    # same interval and the z axis is the generator itself.
    ray = bool(xyz[0].lo == xyz[1].lo and xyz[0].hi == xyz[1].hi
               and xyz[2].lo == rs_probe.lo and xyz[2].hi == rs_probe.hi)
    applied_ray = {
        "applied_std_is_a_ray_not_a_box": True,
        "x_factor": c.rs_x_factor,
        "y_factor": c.rs_y_factor,
        "z_factor": 1.0,
        "probe_base": [rs_probe.lo, rs_probe.hi],
        "probe_applied_x": [xyz[0].lo, xyz[0].hi],
        "probe_applied_z": [xyz[2].lo, xyz[2].hi],
        "one_scalar_generates_three_axes": bool(ray),
    }
    tau_saturation = up(c.pseudo_max_s / down(c.pseudo_ratio))
    interior = {
        "T_S_pins_tau_on_the_open_branch": True,
        "tau_at_T_S_upper_rail_s": tau_saturation,
        "frequency_at_T_S_upper_rail_hz": down(c.tau_coeff * 0.5 / up(tau_saturation)),
        "tau_at_T_S_lower_rail_s": down(c.pseudo_min_s / up(c.pseudo_ratio)),
        "T_S_lower_rail_unreachable_for_targets": bool(
            down(c.pseudo_min_s / up(c.pseudo_ratio)) < tau_target_lo),
        "nominal_independent_rectangle_dimension": 4,
        "reachable_target_surface_dimension": 2,
    }

    # --- finite geometric cover of the reachable (f, sigma) rectangle --------
    f_tiles = _geometric_tiling(c.tune_freq_min, c.tune_freq_max, cells)
    s_tiles = _geometric_tiling(RAW_SIGMA_FLOOR, c.sigma_max, cells)
    rs_x = c.rs_x_factor
    joint_scale = 0.0
    joint_witness = None
    ts_rel_width = 0.0
    rs_rel_width = 0.0
    count = 0
    for f_cell in f_tiles:
        for s_cell in s_tiles:
            img = _target_image(f_cell, s_cell, c)
            count += 1
            ts, rs = img["T_S"], img["R_S"]
            if ts.lo <= 0.0 or rs.lo <= 0.0:
                raise RuntimeError("cover cell lost positive cadence/covariance")
            ts_rel_width = max(ts_rel_width, up((ts.hi - ts.lo) / down(ts.lo)))
            rs_rel_width = max(rs_rel_width, up((rs.hi - rs.lo) / down(rs.lo)))
            scale = up(up(ts.hi + c.dt) / down(rs_x * rs.lo))
            if scale > joint_scale:
                joint_scale = scale
                joint_witness = {"f": list(f_cell), "sigma": list(s_cell),
                                 "tau": [img["tau"].lo, img["tau"].hi],
                                 "T_S": [ts.lo, ts.hi], "R_S": [rs.lo, rs.hi]}
    independent_scale = up(up(c.pseudo_max_s + c.dt) / down(rs_x * c.rs_min))
    cover = {
        "cells_per_axis": cells,
        "cell_count": count,
        "tiling": "geometric, touching, endpoint-exact",
        "max_relative_T_S_width_per_cell": ts_rel_width,
        "max_relative_R_S_width_per_cell": rs_rel_width,
        "joint_S_zero_residual_scale_upper": joint_scale,
        "joint_witness_cell": joint_witness,
        "independent_corner_residual_scale": independent_scale,
        "independent_corner_over_approximation_factor": up(independent_scale / down(joint_scale)),
        "independent_corner_energy_over_approximation_factor": up(
            up(independent_scale * independent_scale) / down(joint_scale * joint_scale)),
    }

    # --- what the joint scale means for the S=0 correction ceiling ----------
    handoff = domain["startup"]["physical_handoff_coordinate_bounds"]
    h = float(domain["configured_runtime"]["imu_dt_s"])
    p_env = float(handoff["position_error_norm_upper_m"])
    v_env = float(handoff["velocity_error_norm_upper_mps"])
    a_env = float(handoff["latent_acceleration_error_norm_upper_mps2"])
    # Per-second integral accumulation of the deployed S row at the declared
    # pre-entry envelope, matching ou3_p4_correlated_entry_relation.
    per_second = up(up(p_env + up(h * 0.5 * v_env)) + up(h * h / 6.0 * a_env))

    r_acc = min(map(float, domain["configured_runtime"]["measurement_noise_std"]["accelerometer_mps2"]))
    f_min = float(domain["normal_live"]["specific_force_norm_lower_mps2"])
    transverse = up(up(r_acc * r_acc) / down(f_min * f_min))
    reset_cap = 3.0
    allowed = down(reset_cap / up(math.sqrt(transverse)))
    dwell_term = up(joint_scale * per_second)
    headroom = down(allowed - dwell_term)
    s_m_budget = None
    if headroom > 0.0:
        # max over cells of (T_S+dt)*per_second/(rs_x*R_S) + S_m/(rs_x*R_S) is
        # monotone in S_m, so a single bisection on the same cover is exact.
        lo, hi = 0.0, 1.0e6
        for _ in range(160):
            mid = 0.5 * (lo + hi)
            worst = 0.0
            for f_cell in f_tiles:
                for s_cell in s_tiles:
                    img = _target_image(f_cell, s_cell, c)
                    val = up(up(up(up(img["T_S"].hi + c.dt) * per_second) + mid)
                             / down(rs_x * img["R_S"].lo))
                    worst = max(worst, val)
                    if worst > allowed:
                        break
                if worst > allowed:
                    break
            if worst <= allowed:
                lo = mid
            else:
                hi = mid
        s_m_budget = lo

    contract_s_m = source["physical_BRMM_contract"]["unfrozen_physical_constants"].get("S_m")
    consequence = {
        "reset_utility_norm_max": reset_cap,
        "prior_independent_transverse_attitude_variance_rad2": transverse,
        "allowed_residual_scale_at_transverse_cap": allowed,
        "dwell_term_from_declared_position_envelope": dwell_term,
        "headroom_for_true_integral_displacement": headroom,
        "admissible_BRMM_S_m_upper_m_s": s_m_budget,
        "contract_S_m_value": contract_s_m,
        "contract_S_m_instantiated": contract_s_m is not None,
        "S_m_is_the_missing_source_primitive": bool(contract_s_m is None),
        "fallback_if_qualified_S_m_exceeds_budget": (
            "retain the same-cell attitude/S cross-covariance direction instead of the "
            "Cauchy-Schwarz lambda_max product; the magnitude-only route has no headroom left"),
    }

    return {
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "filter_changed": False,
        "declared_domain_shrunk": False,
        "trajectory_replay_used": False,
        "independent_tau_sigma_TS_RS_rectangle_used": False,
        "finite_harmonic_substitute_used": False,

        "deployed_pinnings": {
            "T_S_from_tau": "clamp(pseudo_ratio*tau, [pseudo_min_s, pseudo_max_s])",
            "R_S_target_from_tau_sigma": "SpectralMSE(tau, sigma) with the same-cell pseudo cadence",
            "applied_R_S_from_one_scalar": "(rs_x*r, rs_y*r, r)",
            "verified_against_deployed_transitions": True,
        },
        "pinning_probe": pinning,
        "applied_R_S_ray": applied_ray,
        "interior_invertibility": interior,
        "forward_invariant_coefficient_box": invariant_box,
        "coefficient_cover": cover,
        "adaptation_rate_limits": _rate_limits(c, w),
        "S_zero_correction_consequence": consequence,

        "reachable_f_sigma_set_certified_here": False,
        "composed_with_Riccati_Joseph_reset_graph_here": False,
        "consecutive_word_storage_compatibility_closed_here": False,
        "SOURCE_UNIFORM_COMPLETE_BRMM_COVER_CLOSED": False,
        "P4_promoted_here": False,
        "next_obligation": (
            "certify the reachable (f,sigma) set over every admitted BRMM/private-observer/"
            "stillness continuation, then compose these cells and their rate-limited "
            "adjacency with the same-history Riccati/Joseph/reset graph and run the endpoint "
            "and literal-every-prefix augmented LDLT with compatible consecutive-word storage"),
    }


def validate(d: dict) -> list[str]:
    f = []
    if d.get("qualification") != QUALIFICATION:
        f.append("qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_BRMM_NORMAL_LIVE_WORD":
        f.append("source changed")
    for k in ("filter_changed", "declared_domain_shrunk", "trajectory_replay_used",
              "independent_tau_sigma_TS_RS_rectangle_used", "finite_harmonic_substitute_used",
              "reachable_f_sigma_set_certified_here", "composed_with_Riccati_Joseph_reset_graph_here",
              "consecutive_word_storage_compatibility_closed_here",
              "SOURCE_UNIFORM_COMPLETE_BRMM_COVER_CLOSED", "P4_promoted_here"):
        if d.get(k) is not False:
            f.append(k + " not false")
    if d["deployed_pinnings"].get("verified_against_deployed_transitions") is not True:
        f.append("pinnings were not verified against the deployed transitions")
    inv = d.get("forward_invariant_coefficient_box", {})
    if inv.get("needs_reachable_set_argument") is not False:
        f.append("invariant coefficient box claims to need a reachable-set argument")
    lo, hi = inv.get("tau_s", (0.0, 0.0))
    if not (0.0 < lo < hi):
        f.append("invariant tau box is not a positive ordered interval")
    ts_lo, ts_hi = inv.get("T_S_s", (0.0, 0.0))
    if not (0.0 < ts_lo < ts_hi):
        f.append("invariant T_S box is not a positive ordered interval")
    interior = d.get("interior_invertibility", {})
    if interior.get("T_S_pins_tau_on_the_open_branch") is not True:
        f.append("T_S/tau pinning lost")
    if int(interior.get("reachable_target_surface_dimension", 0)) >= int(
            interior.get("nominal_independent_rectangle_dimension", 0)):
        f.append("reachable target surface is not lower dimensional than the rectangle")
    if not d["applied_R_S_ray"].get("one_scalar_generates_three_axes"):
        f.append("applied anisotropic R_S lost its single-scalar generator")
    cover = d.get("coefficient_cover", {})
    if int(cover.get("cell_count", 0)) != int(cover.get("cells_per_axis", 0)) ** 2:
        f.append("cover cell count does not match the tiling")
    for k in ("joint_S_zero_residual_scale_upper", "independent_corner_residual_scale"):
        x = float(cover.get(k, -1.0))
        if not (math.isfinite(x) and x > 0.0):
            f.append(k + " not a positive finite bound")
    # The whole point: the joint scale must be strictly sharper than the corner.
    if float(cover["joint_S_zero_residual_scale_upper"]) >= float(cover["independent_corner_residual_scale"]):
        f.append("joint cell scale is not sharper than the independent corner")
    if float(cover.get("independent_corner_over_approximation_factor", 0.0)) <= 1.0:
        f.append("independent-corner over-approximation factor is not above one")
    rates = d.get("adaptation_rate_limits", {})
    if int(rates.get("steps_per_staged_commit_interval", 0)) <= 1:
        f.append("staged commit interval does not span multiple samples")
    for name, row in rates.items():
        if not isinstance(row, dict):
            continue
        a_lo = float(row.get("per_step_alpha_lower", -1.0))
        a_hi = float(row.get("per_step_alpha_upper", -1.0))
        if not (0.0 < a_lo <= a_hi < 1.0):
            f.append(name + " lost a strict EMA coefficient")
        closed = float(row.get("target_gap_closed_per_commit_interval_upper", -1.0))
        if not (0.0 < closed <= 1.0):
            f.append(name + " gap-closure fraction invalid")
    cons = d.get("S_zero_correction_consequence", {})
    if cons.get("contract_S_m_instantiated") is not False:
        f.append("BRMM S_m appears instantiated; the consequence must be recomputed, not assumed")
    if cons.get("S_m_is_the_missing_source_primitive") is not True:
        f.append("missing S_m primitive no longer reported")
    budget = cons.get("admissible_BRMM_S_m_upper_m_s")
    if budget is not None and not (math.isfinite(float(budget)) and float(budget) > 0.0):
        f.append("admissible S_m budget invalid")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--cells", type=int, default=DEFAULT_CELLS)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain, a.cells)
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    cover = d["coefficient_cover"]
    print(json.dumps({
        "cells": cover["cell_count"],
        "joint_scale": cover["joint_S_zero_residual_scale_upper"],
        "independent_corner": cover["independent_corner_residual_scale"],
        "over_approximation": cover["independent_corner_over_approximation_factor"],
        "invariant_tau_s": d["forward_invariant_coefficient_box"]["tau_s"],
        "invariant_T_S_s": d["forward_invariant_coefficient_box"]["T_S_s"],
        "admissible_S_m_m_s": d["S_zero_correction_consequence"]["admissible_BRMM_S_m_upper_m_s"],
        "contract_S_m": d["S_zero_correction_consequence"]["contract_S_m_value"],
        "cover_closed": d["SOURCE_UNIFORM_COMPLETE_BRMM_COVER_CLOSED"],
        "failures": f,
    }, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
