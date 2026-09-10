#!/usr/bin/env python3
"""Machine-readable A-E falsification classification of the open P4 blockers.

Saying "P4 does not close" is only useful with the reason attached.  This
producer attributes each remaining mathematical blocker of
``ou3_p4_final_closure_gate`` to one of the declared classes

  A  an actual admissible nonlinear trajectory/source counterexample;
  B  rigorous source-uniform certificate infeasibility for this theorem/domain;
  C  enclosure/conditioning/dependency-loss failure;
  D  entry-set modeling failure;
  E  missing source qualification/materialization.

Only A or a genuinely rigorous B would support the statement that P4 is not
provable on the declared physical domain.  C-E are defects of the current proof
construction and are reported as such.

The attribution is quantitative, not editorial.  Two numbers decide the
correction/reset blocker.

1.  The same-cell magnitude certificate is
    ``||d_theta||^2 <= trace(P^-_theta) * y^T R^{-1} y``.  With the retained
    endpoint-referenced covariance envelope the attitude trace upper is about
    1.2e9 rad^2.  Feeding the S=0 residual the correlated integral radius from
    ``ou3_p4_correlated_entry_relation`` instead of the independent 300 m*s
    ball moves the S=0 ceiling but leaves the ACCELEROMETER event limiting at
    about 2.4e6, so the integral entry ball is NOT the limiter of this blocker.

2.  The retained envelope carries about 1.2e9 rad^2 of attitude trace against
    the ``C^2/E_acc`` the reset domain admits, an excess of about 6.4e11.  That
    excess is a conditioning artifact of the endpoint-referenced
    Lagrange/Vandermonde inversion rather than a physical bound, which is what
    makes the blocker class C.

    RETRACTION.  An earlier revision repaired that excess with a
    prior-INDEPENDENT accelerometer posterior cap ``r/|f|^2``, from the scalar
    identity ``P^+ = P^- r/(P^- |f|^2 + r)``, and concluded that the route
    closed on the two transverse directions.  The identity holds only when
    attitude is the sole state in the residual, and the deployed accelerometer
    residual also carries the latent-acceleration and bias blocks, so a
    transverse attitude error and an ``a_w`` error produce the same residual.
    ``ou3_p4_attitude_measurement_cap`` refutes the claim on the deployed
    structure and supplies the CONDITIONAL replacement
    ``(sigma_a^2 + 2 lambda_max(P_(a_w,b_a)))/|f|^2``, about 3.35 rad^2 per
    axis, which is roughly 2800 times larger.

Under the conditional cap the magnitude route no longer closes: the ``R^{-1}``
relaxation gives 180.0 and retaining the same-cell ``S^{-1} = (HPH^T+R)^{-1}``
gives 3.387, against the exact reset utility limit 3.0.  The near miss is sharp
rather than accidental, and the module records why.  Retaining ``S^{-1}``,

    ceiling(P)^2 = 2 P E_acc r / (f^2 P + r),

is increasing in the transverse variance ``P`` with supremum
``sqrt(2 E_acc r/f^2) = 3.3876``.  That supremum is ABOVE 3.0, so the route does
not close for free; but being finite and increasing it yields an exact
threshold ``P* = C^2 r/(2 E_acc r - C^2 f^2) = 4.313e-3`` rad^2 per axis.  The
obligation closes through this route if and only if the transverse attitude
variance is bounded by ``P*``, equivalently if and only if
``lambda_max(P_(a_w,b_a)) <= 5.27e-2`` against the 56.46 currently certified.

So the one-shot measurement route is a dead end for this blocker -- one event
provably cannot separate attitude from ``a_w`` -- but it leaves a single scalar
target for the uniform observability/detectability machinery, which separates
them over a window.  The third direction, yaw about the specific force, is not
accelerometer-observed at all and needs the asynchronous magnetometer plus
gyro-bias transport, which is not established here and is reported open.

Nothing in this module promotes anything.  It consumes the fail-closed gate and
must report ``P4_MOTION_PASS`` exactly as the gate does.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_attitude_measurement_cap as CAP
import ou3_p4_correlated_entry_relation as ENTRYREL
import ou3_p4_exact_reset_transport as RESET
import ou3_p4_final_closure_gate as GATE
import ou3_p4_same_cell_correction_domain as CORR

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
QUALIFICATION = "OU3_P4_BLOCKER_FALSIFICATION_CLASSIFICATION_V1"
CLASSES = ("A", "B", "C", "D", "E")

# Blocker text is keyed exactly as the final gate emits it, so a renamed or
# retired blocker fails validation instead of silently losing its attribution.
BLOCKER_CLASSES = {
    # The three bias families are closed today, so the gate does not emit these.
    # It emits them the moment any family flag regresses, and an unclassified
    # blocker aborts this producer before it can report anything.  Classifying
    # them here keeps a regression legible instead of turning it into a crash.
    "source-uniform BIAS0 same-history physical-driver family": "E",
    "source-uniform BIAS1 same-history physical-driver family": "E",
    "source-uniform BIAS2 same-history physical-driver family": "E",
    "COMPLETE BRMM same-signal estimator/source cover with signal-derived tau and sigma": "E",
    "interdependent (tau,sigma,T_S)->R_S with literal adaptation/commit/scheduler history": "E",
    "source-uniform same-cell Joseph correction/reset domain": "C",
    "source-uniform exact-graph endpoint augmented LDLT": "E",
    "source-uniform exact-graph every-prefix augmented LDLT": "E",
    "same exact graph every-prefix hard-domain retention": "C",
}
_BIAS_FAMILY_REASON = (
    "the family lost one of its own admission, materialized joint ISS supply or "
    "family-parametric same-history graph prerequisites; it is a regression of a "
    "certificate that was closed, never an inference from another family")
BLOCKER_REASONS = {
    "source-uniform BIAS0 same-history physical-driver family": _BIAS_FAMILY_REASON,
    "source-uniform BIAS1 same-history physical-driver family": _BIAS_FAMILY_REASON,
    "source-uniform BIAS2 same-history physical-driver family": _BIAS_FAMILY_REASON,
    "COMPLETE BRMM same-signal estimator/source cover with signal-derived tau and sigma":
        "the joint estimator relation is materialized on one admitted BRMM predecessor only; "
        "the estimator-owned transition operator is not propagated over every admitted "
        "BRMM/private-observer/stillness continuation and hard-entry radial segment",
    "interdependent (tau,sigma,T_S)->R_S with literal adaptation/commit/scheduler history":
        "same missing object as the cover: the candidate/active EMA, staged commit and "
        "scheduler elapsed state exist as transitions but are not carried over the full "
        "admitted predecessor family",
    "source-uniform same-cell Joseph correction/reset domain":
        "the magnitude certificate is fed an attitude covariance envelope about 6.4e11 times "
        "larger than the deployed accelerometer event structure permits, and additionally "
        "relaxes the same-cell S^{-1} to R^{-1}; both are dependency losses, not infeasibility",
    "source-uniform exact-graph endpoint augmented LDLT": "consumer of the missing source cover",
    "source-uniform exact-graph every-prefix augmented LDLT": "consumer of the missing source cover",
    "same exact graph every-prefix hard-domain retention":
        "the retained frozen-map diagnostic leaves the declared chart only through the "
        "independent 300 m*s integral ball; with the correlated integral relation it is "
        "retained, so the open part is the outward nonlinear enclosure, not the domain",
}


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def _transverse_attitude_cap(domain: dict) -> dict:
    """Accelerometer posterior cap on transverse attitude variance.

    RETRACTION.  An earlier revision of this function asserted a PRIOR-INDEPENDENT
    cap ``r/|f|^2``, from the scalar identity ``P^+ = P^- r/(P^- |f|^2 + r)``.
    That identity holds only when attitude is the sole state in the residual.  The
    deployed ``measurement_update_acc_only`` builds ``J_att = -skew(f_cog_b)``
    alongside a latent-acceleration block and an accelerometer-bias block, so a
    transverse attitude error and an ``a_w`` error produce the same residual and
    one update cannot separate them.  ``ou3_p4_attitude_measurement_cap`` refutes
    the prior-independent claim numerically on the deployed residual structure.

    What holds instead is CONDITIONAL on the joint latent/bias covariance block,

        e^T P^+_theta,theta e <= (sigma_a^2 + 2 lambda_max(P_(a_w,b_a))) / |f|^2

    for every prior and every unit ``e`` orthogonal to ``f``, proved from the
    variational posterior identity at ``k = -(f x e)/|f|^2``.  ``lambda_max`` is
    taken from the certified BRMM covariance ceiling, so the cap is a genuine
    bound on the admitted execution rather than an assumption, but it is roughly
    three orders of magnitude larger than the retracted value.  Only the two
    directions transverse to ``f`` are bounded here; the rotation about ``f`` is
    not.
    """
    r_acc = float(CAP.build()["sigma_accelerometer"])
    f_min = float(domain["normal_live"]["specific_force_norm_lower_mps2"])
    if not (r_acc > 0.0 and f_min > 0.0):
        raise RuntimeError("configured accelerometer std / specific-force lower lost positivity")
    if domain["normal_live"]["accelerometer_update_required_each_valid_imu_sample_after_live_entry"] is not True:
        raise RuntimeError("per-sample accelerometer Joseph invariant not declared")
    capd = CAP.build()
    cap_failures = CAP.validate(capd)
    if cap_failures:
        raise RuntimeError("attitude measurement cap invalid: " + repr(cap_failures))
    if capd["prior_independent_transverse_cap_holds_for_the_deployed_filter"] is not False:
        raise RuntimeError("attitude cap module no longer records the refutation")
    per_axis = float(capd["conditional_transverse_cap_rad2"])
    return {
        "argument": (
            "e^T P^+_theta,theta e <= (sigma_a^2 + 2 lambda_max(P_(a_w,b_a)))/|f|^2 "
            "from c^T P^+ c = min_k [(c-H^T k)^T P (c-H^T k) + k^T R k] at k = -(f x e)/|f|^2"),
        "prior_independent": False,
        "retracted_prior_independent_claim": "P^+ <= r/|f|^2 for every prior",
        "retracted_prior_independent_value_rad2": float(
            capd["prior_independent_transverse_cap_claimed"]),
        "refutation_achieved_marginal_rad2": float(
            capd["refutation"]["achieved_transverse_attitude_marginal"]),
        "refutation_exceedance_factor": float(capd["refutation"]["exceedance_factor"]),
        "conditional_on_joint_latent_bias_block": True,
        "joint_latent_bias_lambda_max": float(capd["joint_latent_bias_lambda_max"]),
        "lambda_max_source": "certified BRMM covariance ceiling, mode A",
        "per_sample_accelerometer_update_declared": True,
        "accelerometer_std_lower_mps2": r_acc,
        "specific_force_norm_lower_mps2": f_min,
        "transverse_variance_upper_per_axis_rad2": per_axis,
        "transverse_trace_upper_two_axes_rad2": up(2.0 * per_axis),
        "yaw_about_specific_force_capped_here": False,
        "yaw_needs_asynchronous_magnetometer_and_gyro_bias_transport": True,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain = json.loads(Path(domain_path).read_text(encoding="utf-8"))
    gate = GATE.build()
    gate_failures = GATE.validate(gate)
    if gate_failures:
        raise RuntimeError("final gate invalid: " + repr(gate_failures))
    corr = CORR.build(Path(domain_path))
    corr_failures = CORR.validate(corr)
    if corr_failures:
        raise RuntimeError("correction-domain producer invalid: " + repr(corr_failures))
    rel = ENTRYREL.build(Path(domain_path))
    rel_failures = ENTRYREL.validate(rel)
    if rel_failures:
        raise RuntimeError("correlated entry relation invalid: " + repr(rel_failures))

    blockers = list(gate["remaining_mathematical_P4_blockers"])
    unknown = [b for b in blockers if b not in BLOCKER_CLASSES]
    if unknown:
        raise RuntimeError("unclassified P4 blocker: " + repr(unknown))

    reset_cap = float(RESET.CAYLEY_MONOTONE_NORM_MAX)
    resid = corr["residual_bounds"]
    rs_std_min = float(resid["actual_RS_std_lower_with_horizontal_factor"])
    acc_energy = float(resid["accelerometer_Rinv_energy_upper"])
    mag_energy = float(resid["magnetometer_Rinv_energy_upper"])

    # --- attribution of the correction/reset blocker ---
    correlated_S = float(rel["correlated_radius_one_cadence_m_s"])
    s_energy_correlated = up(up(correlated_S * correlated_S) / down(rs_std_min * rs_std_min))
    attribution = {}
    for mode, m in corr["modes"].items():
        patt = float(m["attitude_covariance_trace_upper"])
        events_declared = {k: float(v["same_cell_attitude_correction_norm_upper"])
                           for k, v in m["events"].items()}
        with_relation = dict(events_declared)
        with_relation["S_zero"] = up(math.sqrt(up(patt * s_energy_correlated)))
        limiting_after = max(with_relation, key=lambda k: with_relation[k])
        required_trace = down(down(reset_cap * reset_cap) / up(acc_energy))
        attribution[mode] = {
            "attitude_covariance_trace_upper_in_use": patt,
            "event_correction_ceilings_declared_entry_set": events_declared,
            "limiting_event_declared_entry_set": m["limiting_event"],
            "event_correction_ceilings_with_correlated_integral_relation": with_relation,
            "limiting_event_with_correlated_integral_relation": limiting_after,
            "integral_entry_ball_is_the_limiter": bool(limiting_after == "S_zero"),
            "reset_utility_norm_max": reset_cap,
            "attitude_trace_required_for_limiting_event": required_trace,
            "attitude_trace_excess_factor": up(patt / down(required_trace)),
            "closes_on_correlated_relation_alone": bool(with_relation[limiting_after] <= reset_cap),
        }

    cap = _transverse_attitude_cap(domain)
    tv_trace = float(cap["transverse_trace_upper_two_axes_rad2"])
    r_acc2 = down(float(cap["accelerometer_std_lower_mps2"]) ** 2)
    f_min2 = down(float(cap["specific_force_norm_lower_mps2"]) ** 2)
    # Same-cell innovation covariance on the transverse block: S = H P H^T + R,
    # so retaining S^{-1} instead of relaxing to R^{-1} discounts the charged
    # residual energy by S/R.  Note this discount grows with the cap, so the
    # product tv_trace/(S/R) saturates: the S^{-1} route cannot be rescued by a
    # looser attitude bound, which is why the conditional cap does not close
    # this obligation.
    hph = up(f_min2 * float(cap["transverse_variance_upper_per_axis_rad2"]))
    s_over_r = down(up(hph + r_acc2) / up(r_acc2))
    magnitude_route = {
        "transverse_trace_upper_rad2": tv_trace,
        "accelerometer_ceiling_with_Rinverse_relaxation": up(math.sqrt(up(tv_trace * acc_energy))),
        "same_cell_S_over_R_lower": s_over_r,
        "accelerometer_ceiling_retaining_same_cell_Sinverse":
            up(math.sqrt(up(up(tv_trace * acc_energy) / down(s_over_r)))),
        "magnetometer_ceiling_with_Rinverse_relaxation": up(math.sqrt(up(tv_trace * mag_energy))),
        "S_zero_ceiling_with_correlated_relation_and_transverse_cap":
            up(math.sqrt(up(tv_trace * s_energy_correlated))),
        "reset_utility_norm_max": reset_cap,
    }
    magnitude_route["Rinverse_relaxation_alone_closes_accelerometer"] = bool(
        magnitude_route["accelerometer_ceiling_with_Rinverse_relaxation"] <= reset_cap)
    magnitude_route["same_cell_Sinverse_closes_accelerometer"] = bool(
        magnitude_route["accelerometer_ceiling_retaining_same_cell_Sinverse"] <= reset_cap)
    magnitude_route["transverse_directions_only"] = True
    magnitude_route["third_attitude_direction_still_open"] = True
    # Under the retracted prior-independent cap the same-cell S^{-1} accounting
    # appeared to bring the accelerometer ceiling inside the reset utility
    # domain. It does not under the conditional cap, and it cannot be made to:
    # S/R grows in proportion to the cap, so tv_trace/(S/R) saturates.
    # --- the sharp threshold the one-shot route needs -------------------------
    # Retaining the same-cell S^{-1} the charged accelerometer correction is
    #
    #     ceiling(P)^2 = 2 P E_acc r / (f^2 P + r),   P = transverse variance/axis,
    #
    # strictly increasing in P with supremum sqrt(2 E_acc r / f^2).  Two facts
    # follow, and they are what actually decide this obligation.
    #
    #  (a) The supremum EXCEEDS the reset utility cap, so no attitude bound makes
    #      this route close by saturation -- the route has a genuine threshold
    #      rather than an asymptote below the cap.
    #  (b) Setting ceiling(P) = C and solving gives the exact threshold
    #      P* = C^2 r / (2 E_acc r - C^2 f^2), positive precisely because of (a).
    #      The obligation closes through this route iff the transverse attitude
    #      variance is bounded by P*.
    #
    # Propagating P* back through the conditional cap (sigma_a^2 + 2 lambda)/f^2
    # turns the open blocker into a single scalar target on the certified
    # latent/bias covariance ceiling: lambda <= (P* f^2 - sigma_a^2)/2.
    two_E_r = down(2.0 * down(acc_energy * r_acc2))
    c2 = down(reset_cap * reset_cap)
    denom = two_E_r - up(c2 * f_min2)
    route_sup = up(math.sqrt(up(two_E_r / down(f_min2))))
    threshold = {
        "ceiling_supremum_over_all_attitude_bounds": route_sup,
        "supremum_exceeds_reset_utility_cap": bool(route_sup > reset_cap),
        "saturates_below_reset_cap_for_free": bool(route_sup <= reset_cap),
    }
    if denom > 0.0:
        p_star = down(down(c2 * r_acc2) / up(denom))
        lam_star = down(down(down(p_star * f_min2) - up(r_acc2)) / 2.0)
        lam_now = float(cap["joint_latent_bias_lambda_max"])
        threshold.update({
            "required_transverse_variance_per_axis_rad2": p_star,
            "required_transverse_std_per_axis_rad": down(math.sqrt(down(p_star))),
            "attained_transverse_variance_per_axis_rad2":
                float(cap["transverse_variance_upper_per_axis_rad2"]),
            "transverse_variance_shortfall_factor":
                up(float(cap["transverse_variance_upper_per_axis_rad2"]) / down(p_star)),
            "required_joint_latent_bias_lambda_max": lam_star,
            "attained_joint_latent_bias_lambda_max": lam_now,
            "joint_latent_bias_shortfall_factor":
                up(lam_now / down(lam_star)) if lam_star > 0.0 else float("inf"),
            "requirement_is_achievable_in_principle": bool(lam_star > 0.0),
        })
    else:
        threshold["required_transverse_variance_per_axis_rad2"] = None
        threshold["requirement_is_achievable_in_principle"] = False
    magnitude_route["one_shot_route_threshold"] = threshold

    magnitude_route["one_shot_measurement_route_closes_the_correction_domain"] = bool(
        magnitude_route["same_cell_Sinverse_closes_accelerometer"]
        and magnitude_route["magnetometer_ceiling_with_Rinverse_relaxation"] <= reset_cap
        and magnitude_route["S_zero_ceiling_with_correlated_relation_and_transverse_cap"] <= reset_cap
        and not magnitude_route["third_attitude_direction_still_open"])

    per_blocker = {
        b: {
            "class": BLOCKER_CLASSES[b],
            "reason": BLOCKER_REASONS[b],
            "supports_saying_P4_is_unprovable": BLOCKER_CLASSES[b] in ("A", "B"),
        }
        for b in blockers
    }
    classes_present = sorted({v["class"] for v in per_blocker.values()})

    return {
        "qualification": QUALIFICATION,
        "canonical_source": gate["canonical_source"],
        "gate_qualification": gate["qualification"],
        "declared_classes": list(CLASSES),
        "remaining_mathematical_P4_blockers": blockers,
        "blocker_classification": per_blocker,
        "classes_present": classes_present,
        "counterexample_class_A_found": False,
        "rigorous_infeasibility_class_B_established": False,
        "P4_unprovable_on_declared_domain_supported": bool(
            any(v["class"] in ("A", "B") for v in per_blocker.values())),
        "frozen_map_diagnostic_promoted_to_lower_bound": False,
        "diagnostic_4p5788_treated_as_counterexample": False,

        "correction_domain_attribution": attribution,
        "integral_entry_ball_limits_correction_domain": bool(
            any(v["integral_entry_ball_is_the_limiter"] for v in attribution.values())),
        "conditional_transverse_attitude_cap": cap,
        "magnitude_route_after_transverse_cap": magnitude_route,

        "correlated_integral_relation": {
            "declared_independent_radius_m_s": float(rel["declared_independent_integral_radius_m_s"]),
            "correlated_one_cadence_radius_m_s": correlated_S,
            "chart_retaining_radius_upper_m_s": float(rel["chart_retaining_integral_radius_upper_m_s"]),
            "one_cadence_meets_chart_threshold": bool(rel["one_cadence_relation_meets_chart_threshold"]),
            "S_zero_gives_uniform_anchor_without_dwell": bool(rel["S_zero_gives_uniform_anchor_without_dwell"]),
            "unconditional_anchor_requires_P5_or_reachable_P_SS_lower_bound": True,
        },

        "finite_precision": {
            "conditional_mathematical_binary32_additive_ISS_closed": bool(
                gate["conditional_full_shipping_finite_precision_additive_ISS_closed"]),
            "target_toolchain_qualified": bool(gate["target_toolchain_finite_precision_qualified"]),
            "missing_toolchain_qualification_is_a_mathematical_P4_blocker": False,
            "host_arithmetic_check_is_deployment_qualification": False,
            "remaining_deployment_blockers": list(gate["remaining_deployment_blockers"]),
        },

        "P4_MOTION_PASS": bool(gate["P4_MOTION_PASS"]),
        "P4_PASS": bool(gate["P4_PASS"]),
        "P5_MAY_START": bool(gate["P5_MAY_START"]),
        "P4_promoted_here": False,
    }


def every_gate_blocker_is_classifiable(gate: dict) -> list[str]:
    """Gate blocker strings this producer could not classify.

    The gate is the only source of blocker text, so any string it can emit has
    to appear in both maps.  Callers use this to fail on a coverage gap instead
    of discovering it as a RuntimeError on the day a prerequisite regresses.
    """
    emitted = set(gate.get("remaining_mathematical_P4_blockers", ()))
    return sorted(emitted - (set(BLOCKER_CLASSES) & set(BLOCKER_REASONS)))


def validate(d: dict) -> list[str]:
    f = []
    if d.get("qualification") != QUALIFICATION:
        f.append("qualification mismatch")
    if set(BLOCKER_CLASSES) != set(BLOCKER_REASONS):
        f.append("blocker class and reason maps disagree")
    if d.get("gate_qualification") != GATE.QUALIFICATION:
        f.append("classification is not bound to the current final gate")
    if tuple(d.get("declared_classes", ())) != CLASSES:
        f.append("declared class set changed")
    blockers = d.get("remaining_mathematical_P4_blockers", [])
    table = d.get("blocker_classification", {})
    if set(table) != set(blockers):
        f.append("classification table does not cover exactly the open blockers")
    for b, row in table.items():
        if row.get("class") not in CLASSES:
            f.append(b + " has no declared class")
        if bool(row.get("supports_saying_P4_is_unprovable")) != (row.get("class") in ("A", "B")):
            f.append(b + " unprovability flag inconsistent with its class")
    for k in ("counterexample_class_A_found", "rigorous_infeasibility_class_B_established",
              "frozen_map_diagnostic_promoted_to_lower_bound",
              "diagnostic_4p5788_treated_as_counterexample", "P4_promoted_here"):
        if d.get(k) is not False:
            f.append(k + " not false")
    fp = d.get("finite_precision", {})
    for k in ("missing_toolchain_qualification_is_a_mathematical_P4_blocker",
              "host_arithmetic_check_is_deployment_qualification"):
        if fp.get(k) is not False:
            f.append("finite_precision." + k + " not false")
    if fp.get("conditional_mathematical_binary32_additive_ISS_closed") is not True:
        f.append("conditional binary32 additive ISS channel regressed")
    # The unprovability conclusion must be the exact disjunction of A/B findings.
    supported = any(row.get("class") in ("A", "B") for row in table.values())
    if bool(d.get("P4_unprovable_on_declared_domain_supported")) != supported:
        f.append("unprovability conclusion is not the exact A/B disjunction")
    # Promotion bits are mirrored from the gate and may never be raised here.
    if bool(d.get("P4_MOTION_PASS")) or bool(d.get("P4_PASS")) or bool(d.get("P5_MAY_START")):
        if not blockers:
            pass
        else:
            f.append("promotion bit set while blockers remain")
    cap = d.get("conditional_transverse_attitude_cap", {})
    if cap.get("prior_independent") is not False:
        f.append("transverse attitude cap claims prior independence again")
    if cap.get("conditional_on_joint_latent_bias_block") is not True:
        f.append("transverse attitude cap dropped its latent/bias conditioning")
    if not (float(cap.get("refutation_exceedance_factor", 0.0)) > 1.0):
        f.append("retracted prior-independent cap is no longer shown to be refuted")
    if cap.get("yaw_about_specific_force_capped_here") is not False:
        f.append("yaw direction was claimed capped by the accelerometer event")
    x = float(cap.get("transverse_variance_upper_per_axis_rad2", -1.0))
    if not (math.isfinite(x) and x > 0.0):
        f.append("transverse attitude cap is not a positive finite bound")
    route = d.get("magnitude_route_after_transverse_cap", {})
    for k in ("transverse_directions_only", "third_attitude_direction_still_open"):
        if route.get(k) is not True:
            f.append("magnitude_route." + k + " not true")
    if route.get("one_shot_measurement_route_closes_the_correction_domain") is not False:
        f.append("one-shot measurement route reported as closing the correction domain")
    th = route.get("one_shot_route_threshold", {})
    sup = float(th.get("ceiling_supremum_over_all_attitude_bounds", -1.0))
    if not (math.isfinite(sup) and sup > 0.0):
        f.append("one-shot route supremum is not a positive finite bound")
    # The threshold is only meaningful because the supremum sits above the reset
    # utility cap. If that ever inverts the route closes by saturation and this
    # whole attribution has to be redone rather than reported as-is.
    if bool(th.get("supremum_exceeds_reset_utility_cap")) == bool(
            th.get("saturates_below_reset_cap_for_free")):
        f.append("one-shot route saturation verdict is not a dichotomy")
    if th.get("supremum_exceeds_reset_utility_cap") is True:
        pstar = th.get("required_transverse_variance_per_axis_rad2")
        if pstar is None or not (float(pstar) > 0.0):
            f.append("one-shot route threshold missing while the supremum exceeds the cap")
        elif float(th.get("attained_transverse_variance_per_axis_rad2", 0.0)) <= float(pstar):
            # Attaining the threshold would close the route, which contradicts
            # the flag above; refuse to publish both.
            f.append("attitude cap meets the threshold yet the route is reported open")
        elif not (float(th.get("transverse_variance_shortfall_factor", 0.0)) > 1.0):
            f.append("threshold shortfall factor disagrees with the attained cap")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain)
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "classes_present": d["classes_present"],
        "blocker_classes": {b: r["class"] for b, r in d["blocker_classification"].items()},
        "P4_unprovable_supported": d["P4_unprovable_on_declared_domain_supported"],
        "integral_ball_limits_correction_domain": d["integral_entry_ball_limits_correction_domain"],
        "attitude_trace_excess_factor": {m: v["attitude_trace_excess_factor"]
                                         for m, v in d["correction_domain_attribution"].items()},
        "transverse_cap_rad2": d["conditional_transverse_attitude_cap"]["transverse_trace_upper_two_axes_rad2"],
        "acc_ceiling_Rinv": d["magnitude_route_after_transverse_cap"]["accelerometer_ceiling_with_Rinverse_relaxation"],
        "acc_ceiling_Sinv": d["magnitude_route_after_transverse_cap"]["accelerometer_ceiling_retaining_same_cell_Sinverse"],
        "P4_MOTION_PASS": d["P4_MOTION_PASS"],
        "failures": f,
    }, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
