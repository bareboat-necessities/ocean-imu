#!/usr/bin/env python3
"""Physically qualified correlated P4 entry relation for the integral state.

The declared hard entry set admits ``||e_S|| <= 300 m*s`` as an INDEPENDENT
Cartesian ball alongside independent attitude/velocity/position/latent balls.
``S`` is the integral-displacement state, so that treatment is physically wrong
in a specific, checkable way: the deployed translational factor makes ``S`` the
exact running integral of ``p``,

    S_next = S + dt*p + dt^2/2*v + phi_Sa*a_w,      0 < phi_Sa <= dt^3/6,

and the same factor drives truth and estimate, so the ERROR obeys the same
relation exactly, with no linearization remainder:

    e_S_next = e_S + dt*e_p + dt^2/2*e_v + phi_Sa*e_aw.

``e_S`` is therefore slaved to the (e_v,e_p,e_aw) history and to the deployed
repeated ``S=0`` Joseph regulation, whose error form is exact as well:

    e_S^+ = (I - K_S) e_S^- - K_S S_true,   K_S = P_SS (P_SS + R_S)^{-1}.

This module does NOT shrink 300 to a convenient constant.  It derives what the
deployed recurrence and scheduler actually permit, reports the thresholds the
integral radius has to beat, and reports explicitly which part of the
qualification is not available inside P4.

Findings, all reproduced by ``build()``:

1.  The declared 300 m*s ball has no derivation in the operating domain and is
    not conservative either: at the declared startup position envelope the free
    integral reaches ``T_handoff * 20 m*s``, above 600 m*s over the declared
    live-entry timing floor.  It is an arbitrary constant, i.e. a class-D
    entry-set modeling defect rather than a physical bound.

2.  Over a dwell window ``W`` inside the declared handoff envelope the
    accumulation is exact and small, because the deployed pseudo scheduler
    fires ``S=0`` at most one cadence plus one sample apart, not once per
    handoff interval.

3.  The deployed ``S=0`` regulation is a NON-EXPANSION of ``S_hat`` in the
    ``R_S``-weighted norm but is NOT a uniform contraction over the admitted
    cell family: at the reachable one-step process floor the weighted
    contraction factor is ``1 - mu`` with ``mu`` far below 1e-6, i.e. no useful
    anchoring.  Anchoring ``e_S`` without a dwell hypothesis therefore needs
    either a proved reachable ``P_SS`` lower bound far above the process floor
    or the P5 capture argument.  That is reported, not smuggled into P4.

Relationship to the deployed pinning.  ``ou3_p4_live_entry_reachability``
settles the entry VALUE from the startup path: the MEKF linear block is never
propagated before ``goLive``, so ``e_S`` at Live entrance is exactly
``-S_true(T)`` and its radius is the BRMM primitive ``S_m``.  This module stays
the fallback description for any admitted history that is already inside Live,
and it supplies two results the pinning does not: the chart threshold the
integral radius must beat, and the fact that the ``S=0`` regulation gives no
uniform contraction to anchor ``e_S`` with.

Nothing here promotes P4 and no declared radius is reduced: the relation is an
additional hard qualification whose consumption stays gated.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import ou3_brmm_dynamic_source_certificate as DYNAMIC
import ou3_p4_exact_reset_transport as RESET

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
RETENTION_EVIDENCE = REPO / "reports" / "results" / "rao_stability" / "entry-block-retention.json"
QUALIFICATION = "OU3_P4_CORRELATED_INTEGRAL_ENTRY_RELATION_V1"

# Shipping SpectralMSE applies the horizontal anisotropy factor to R_S.  The
# smallest applied standard deviation gives the largest residual energy.
RS_HORIZONTAL_FACTOR = 0.72

# Declared live-entry timing floor: 4/lambda with lambda = 2*pi*high_pass_hz,
# plus at least one estimated period of accumulated moment history.
LIVE_ENTRY_TIMING_FLOOR_S = 30.0

# Physically expected shipping attitude covariance trace, used only to report
# the integral radius that the same-cell Joseph correction ceiling admits once
# the attitude covariance envelope is repaired.  It is a REPORTING reference,
# never a proof premise: the correction-domain producer keeps its own envelope.
REFERENCE_ATTITUDE_TRACE = 1.0e-3


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def _accumulation_upper(window_s: float, h: float, p_env: float, v_env: float, a_env: float) -> dict:
    """Outward upper bound on ``||e_S(t0) - e_S(t0-W)||`` over a dwell window.

    The exact per-step increment is ``dt*e_p + dt^2/2*e_v + phi_Sa*e_aw`` with
    ``0 < phi_Sa <= dt^3/6``, because
    ``phi_Sa = int_0^dt (dt-s)^2/2 * exp(-s/tau) ds`` and the exponential lies
    in ``(0,1]``.  The bound is uniform in tau and deliberately avoids the
    shipping closed form ``tau^3*(x^2/2 - x + 1 - exp(-x))``, which loses every
    significant digit under plain interval arithmetic at the admitted
    ``x = dt/tau <= 0.012``.
    """
    if not (window_s > 0.0 and h > 0.0):
        raise ValueError("positive dwell window and sample period required")
    steps = int(math.ceil(up(window_s) / down(h)))
    span = up(steps * h)
    phi_sa_upper = up(h * h * h / 6.0)
    from_p = up(span * p_env)
    from_v = up(up(span * h) * 0.5 * v_env)
    from_a = up(up(steps * phi_sa_upper) * a_env)
    return {
        "window_s": window_s,
        "steps": steps,
        "covered_span_s": span,
        "phi_Sa_upper_per_step": phi_sa_upper,
        "from_position_error_m_s": from_p,
        "from_velocity_error_m_s": from_v,
        "from_latent_acceleration_error_m_s": from_a,
        "accumulation_upper_m_s": up(up(from_p + from_v) + from_a),
    }


def _s_zero_contraction(sigma_aw_lo: float, tau_hi: float, h: float, rs_std_hi: float) -> dict:
    """Uniform ``S=0`` contraction available over the admitted cell family.

    ``e_S^+ = (I-K_S) e_S^- - K_S S_true`` with ``K_S = P_SS (P_SS+R_S)^{-1}``.
    The ``R_S^{1/2}`` similarity makes ``I-K_S`` symmetric in the ``R_S^{-1}``
    weighted norm with norm ``1/(1+mu)`` and
    ``mu = lambda_min(R_S^{-1/2} P_SS R_S^{-1/2}) >= 0``, so the update never
    expands ``S_hat`` in that norm.  A uniform CONTRACTION needs a uniform
    positive lower bound on ``mu``, hence on the reachable ``P_SS``.  The only
    source-uniform lower bound available without a reachability argument is one
    prediction step of shipping process noise, ``P^- >= Q``, and

        q_SS = (2 sigma^2/tau) int_0^h phi_Sa(u)^2 du
             >= sigma^2 h^7 (1 - x/4)^2 / (126 tau)

    from ``phi_Sa(u) >= (u^3/6) (1 - u/(4 tau))``.
    """
    if not (sigma_aw_lo > 0.0 and tau_hi > 0.0 and h > 0.0 and rs_std_hi > 0.0):
        raise ValueError("positive sigma/tau/h/R_S required")
    x_hi = up(h / down(tau_hi))
    shrink = down(1.0 - up(x_hi / 4.0))
    if shrink <= 0.0:
        raise RuntimeError("admitted h/tau destroys the phi_Sa lower bound")
    q_ss = down(down(down(sigma_aw_lo * sigma_aw_lo) * down(pow(h, 7))) * down(shrink * shrink)
                / up(126.0 * tau_hi))
    mu = down(q_ss / up(rs_std_hi * rs_std_hi))
    # mu >= 0 makes 1/(1+mu) <= 1 exactly, so the outward upper is capped at 1:
    # the weighted update is a non-expansion whatever the rounding does.
    lam = min(1.0, up(1.0 / down(1.0 + mu))) if mu > 0.0 else 1.0
    return {
        "P_SS_lower_argument": "P^- >= Q over one shipping prediction step",
        "process_floor_q_SS_per_step": q_ss,
        "worst_cell_mu_lower": mu,
        "worst_cell_weighted_contraction_factor_upper": lam,
        "weighted_non_expansion_always_holds": True,
        "uniform_contraction_useful": bool(mu > 1.0e-6),
        "events_per_e_fold_lower": up(1.0 / mu) if mu > 0.0 else math.inf,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, retention_path: Path = RETENTION_EVIDENCE) -> dict:
    domain = json.loads(Path(domain_path).read_text(encoding="utf-8"))
    dynamic = DYNAMIC.build(Path(domain_path))
    failures = DYNAMIC.validate(dynamic)
    if failures:
        raise RuntimeError("correlated-entry prerequisites failed: " + repr(failures))

    h = float(domain["configured_runtime"]["imu_dt_s"])
    handoff = domain["startup"]["physical_handoff_coordinate_bounds"]
    declared_S = float(handoff["integral_displacement_error_norm_upper_m_s"])
    p_env = float(handoff["position_error_norm_upper_m"])
    v_env = float(handoff["velocity_error_norm_upper_mps"])
    a_env = float(handoff["latent_acceleration_error_norm_upper_mps2"])

    inv = dynamic["dynamic_invariant"]
    ts_hi = float(inv["pseudo_update_period_s"][1])
    tau_hi = float(inv["tau_applied_s"][1])
    sigma_lo = float(inv["sigma_aw_filter_mps2"][0])
    rs_lo, rs_hi = float(inv["R_S_applied"][0]), float(inv["R_S_applied"][1])
    rs_std_min = down(rs_lo * RS_HORIZONTAL_FACTOR)
    if not rs_std_min > 0.0:
        raise RuntimeError("applied R_S lower lost positivity")

    # The deployed scheduler cannot leave more than one cadence plus one sample
    # between consecutive S=0 Joseph events.
    gap_upper = up(ts_hi + h)
    word_s = 3.0

    windows = {
        "one_S_zero_cadence_gap": _accumulation_upper(gap_upper, h, p_env, v_env, a_env),
        "one_P4_word": _accumulation_upper(word_s, h, p_env, v_env, a_env),
    }
    contraction = _s_zero_contraction(sigma_lo, tau_hi, h, rs_hi)

    # --- threshold 1: chart retention of the retained frozen-map diagnostic ---
    retention = json.loads(Path(retention_path).read_text(encoding="utf-8"))
    retention_sha = hashlib.sha256(Path(retention_path).read_bytes()).hexdigest()
    if float(retention["declared_entry_radii"]["integral_displacement"]) != declared_S:
        raise RuntimeError("retained retention diagnostic uses a different integral radius")
    chart = {}
    for mode, m in retention["modes"].items():
        r_att = float(retention["declared_entry_radii"]["attitude"])
        chart_cap = float(m["declared_chart_cayley_norm_upper"])
        without = float(m["without_independent_integral_ball"]["attitude"])
        alone = float(m["single_ball_reach"]["attitude"]["integral_displacement"])
        # Scaling the integral radius by s scales that one term of the
        # subadditive per-prefix sum exactly by s, so every prefix satisfies
        # total <= without + s*alone by max(A+B) <= max A + max B.
        head = down(down(chart_cap / up(r_att)) - without)
        s_max = down(head / up(alone)) if (alone > 0.0 and head > 0.0) else 0.0
        chart[mode] = {
            "attitude_entry_radius": r_att,
            "declared_chart_cayley_norm_upper": chart_cap,
            "attitude_reach_without_independent_integral_ball": without,
            "attitude_reach_from_integral_ball_alone_at_declared_radius": alone,
            "admissible_integral_radius_scale_upper": s_max,
            "chart_retaining_integral_radius_upper_m_s": down(s_max * declared_S),
        }
    binding_mode = min(chart, key=lambda k: chart[k]["chart_retaining_integral_radius_upper_m_s"])
    chart_threshold = chart[binding_mode]["chart_retaining_integral_radius_upper_m_s"]

    # --- threshold 2: same-cell Joseph correction ceiling at the S=0 event ---
    # ||d_theta||^2 <= trace(P^-_theta) * ||y||^2 / rs_std_min^2 with the S=0
    # innovation y = -(e_S + S_true).  Solve for the residual that keeps the
    # ceiling inside the exact reset utility domain.
    reset_cap = float(RESET.CAYLEY_MONOTONE_NORM_MAX)
    correction = {
        "reset_utility_norm_max": reset_cap,
        "reference_attitude_covariance_trace": REFERENCE_ATTITUDE_TRACE,
        "reference_is_reporting_only_not_a_proof_premise": True,
        "S_residual_upper_m_s": down(down(reset_cap * rs_std_min)
                                     / up(math.sqrt(REFERENCE_ATTITUDE_TRACE))),
    }

    accum_gap = windows["one_S_zero_cadence_gap"]["accumulation_upper_m_s"]
    accum_word = windows["one_P4_word"]["accumulation_upper_m_s"]

    # Dwell window that meets the binding chart threshold at the declared envelope.
    per_second = up(up(p_env + up(h * 0.5 * v_env)) + up(h * h / 6.0 * a_env))
    dwell_required_s = down(chart_threshold / per_second)

    free_integral_over_live_entry = up(p_env * LIVE_ENTRY_TIMING_FLOOR_S)

    return {
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "declared_domain_shrunk": False,
        "filter_changed": False,
        "trajectory_replay_used": False,
        "covariance_confidence_ellipsoid_used": False,
        "P5_capture_assumed": False,

        "sample_period_s": h,
        "deployed_S_row": "S_next = S + dt*p + dt^2/2*v + phi_Sa*a_w",
        "deployed_S_row_coefficients": {
            "v": up(0.5 * h * h),
            "p": h,
            "S": 1.0,
            "a_w_upper": up(h * h * h / 6.0),
            "a_w_bound_argument": "phi_Sa = int_0^dt (dt-s)^2/2 * exp(-s/tau) ds in (0, dt^3/6]",
        },
        "error_obeys_same_row_exactly": True,
        "S_zero_error_update": "e_S^+ = (I-K_S) e_S^- - K_S S_true, K_S = P_SS (P_SS+R_S)^{-1}",
        "physical_minus_S_true_forcing_retained": True,

        "declared_independent_integral_radius_m_s": declared_S,
        "declared_independent_radius_has_operating_domain_derivation": False,
        "declared_independent_radius_is_conservative": bool(declared_S >= free_integral_over_live_entry),
        "declared_live_entry_timing_floor_s": LIVE_ENTRY_TIMING_FLOOR_S,
        "free_integral_reach_over_declared_live_entry_floor_m_s": free_integral_over_live_entry,

        "pre_entry_envelope": {
            "position_error_norm_upper_m": p_env,
            "velocity_error_norm_upper_mps": v_env,
            "latent_acceleration_error_norm_upper_mps2": a_env,
            "source": "operating domain startup.physical_handoff_coordinate_bounds, unchanged",
        },
        "S_zero_scheduler": {
            "pseudo_update_period_upper_s": ts_hi,
            "consecutive_event_gap_upper_s": gap_upper,
            "applied_R_S_std_interval_m_s": [rs_lo, rs_hi],
            "applied_R_S_std_lower_with_horizontal_factor": rs_std_min,
        },
        "dwell_accumulation": windows,
        "correlated_entry_relation": (
            "||e_S(t0)|| <= ||e_S(t0-W)|| + W*||e_p||_env + W*dt/2*||e_v||_env "
            "+ (W/dt)*phi_Sa_upper*||e_aw||_env, with the S=0 update non-expansive on S_hat"
        ),
        "correlated_radius_one_cadence_m_s": accum_gap,
        "correlated_radius_one_word_m_s": accum_word,

        "chart_retention_threshold": chart,
        "chart_retention_binding_mode": binding_mode,
        "chart_retaining_integral_radius_upper_m_s": chart_threshold,
        "declared_radius_meets_chart_threshold": bool(declared_S <= chart_threshold),
        "one_cadence_relation_meets_chart_threshold": bool(accum_gap <= chart_threshold),
        "one_word_relation_meets_chart_threshold": bool(accum_word <= chart_threshold),
        "chart_threshold_is_frozen_map_diagnostic_not_outward_certificate": True,
        "chart_threshold_evidence_sha256": retention_sha,

        "correction_ceiling_threshold": correction,
        "one_cadence_relation_meets_correction_ceiling": bool(accum_gap <= correction["S_residual_upper_m_s"]),
        "declared_radius_meets_correction_ceiling": bool(declared_S <= correction["S_residual_upper_m_s"]),
        "dwell_window_required_for_chart_threshold_s": dwell_required_s,
        "one_cadence_gap_meets_required_dwell": bool(gap_upper <= dwell_required_s),

        "S_zero_regulation": contraction,
        "S_zero_gives_uniform_anchor_without_dwell": bool(contraction["uniform_contraction_useful"]),
        "unconditional_entry_anchor_requires_P5_or_reachable_P_SS_lower_bound": True,
        "falsification_class": "D_entry_set_modeling_plus_E_missing_source_qualification",
        "entry_relation_consumed_by_final_gate": False,
        "entry_radii_reduced_here": False,
        "P4_promoted_here": False,
    }


def validate(d: dict) -> list[str]:
    f = []
    if d.get("qualification") != QUALIFICATION:
        f.append("qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_BRMM_NORMAL_LIVE_WORD":
        f.append("source changed")
    for k in ("error_obeys_same_row_exactly", "physical_minus_S_true_forcing_retained",
              "chart_threshold_is_frozen_map_diagnostic_not_outward_certificate",
              "unconditional_entry_anchor_requires_P5_or_reachable_P_SS_lower_bound"):
        if d.get(k) is not True:
            f.append(k + " not true")
    for k in ("declared_domain_shrunk", "filter_changed", "trajectory_replay_used",
              "covariance_confidence_ellipsoid_used", "P5_capture_assumed",
              "declared_independent_radius_has_operating_domain_derivation",
              "entry_relation_consumed_by_final_gate", "entry_radii_reduced_here",
              "P4_promoted_here"):
        if d.get(k) is not False:
            f.append(k + " not false")
    for k in ("correlated_radius_one_cadence_m_s", "correlated_radius_one_word_m_s",
              "chart_retaining_integral_radius_upper_m_s",
              "dwell_window_required_for_chart_threshold_s"):
        x = float(d.get(k, -1.0))
        if not (math.isfinite(x) and x > 0.0):
            f.append(k + " not a positive finite bound")
    # The relation must actually improve on the coordinate it reinterprets.
    if float(d["correlated_radius_one_cadence_m_s"]) >= float(d["declared_independent_integral_radius_m_s"]):
        f.append("one-cadence relation does not improve on the declared independent ball")
    # A useful uniform S=0 contraction has to be demonstrated, never assumed.
    if bool(d["S_zero_regulation"]["uniform_contraction_useful"]) != bool(d["S_zero_gives_uniform_anchor_without_dwell"]):
        f.append("S=0 anchoring claim inconsistent with its own contraction bound")
    if d["S_zero_regulation"]["weighted_non_expansion_always_holds"] is not True:
        f.append("weighted non-expansion of the S=0 update lost")
    if d["correction_ceiling_threshold"]["reference_is_reporting_only_not_a_proof_premise"] is not True:
        f.append("correction-ceiling reference trace was promoted to a premise")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--retention", type=Path, default=RETENTION_EVIDENCE)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain, a.retention)
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "declared_independent_m_s": d["declared_independent_integral_radius_m_s"],
        "declared_is_conservative": d["declared_independent_radius_is_conservative"],
        "correlated_one_cadence_m_s": d["correlated_radius_one_cadence_m_s"],
        "correlated_one_word_m_s": d["correlated_radius_one_word_m_s"],
        "chart_threshold_m_s": d["chart_retaining_integral_radius_upper_m_s"],
        "chart_binding_mode": d["chart_retention_binding_mode"],
        "correction_ceiling_threshold_m_s": d["correction_ceiling_threshold"]["S_residual_upper_m_s"],
        "one_cadence_meets_chart": d["one_cadence_relation_meets_chart_threshold"],
        "one_word_meets_chart": d["one_word_relation_meets_chart_threshold"],
        "S_zero_uniform_anchor": d["S_zero_gives_uniform_anchor_without_dwell"],
        "class": d["falsification_class"],
        "failures": f,
    }, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
