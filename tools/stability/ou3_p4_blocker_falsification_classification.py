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

2.  The accelerometer Joseph event admits a prior-INDEPENDENT posterior cap on
    the two attitude directions transverse to the specific force.  For a scalar
    measurement of angle with gain ``|f|`` and noise variance ``r``,
    ``P^+ = P^- r/(P^- |f|^2 + r) <= r/|f|^2`` whatever ``P^-`` was, and the
    declared Normal-Live invariant executes that update at EVERY valid IMU
    sample.  So the reachable transverse attitude variance is at most
    ``R_acc/|f|_min^2``, about 1.19e-3 rad^2 per axis: the envelope in use
    exceeds what the deployed event structure permits by about 6.4e11.  That is
    a conditioning artifact of the endpoint-referenced Lagrange/Vandermonde
    inversion, not a physical bound, which makes the blocker class C.

The module also quantifies the one remaining lossy step of the magnitude route:
at the transverse cap the ``R^{-1}`` relaxation gives 3.39, just above the exact
reset utility limit 3.0, while retaining the same-cell ``S^{-1} = (HPH^T+R)^{-1}``
gives 2.40.  The route therefore closes on the two transverse directions once
``S^{-1}`` is kept; the third (yaw about the specific force) is not capped by
the accelerometer at all and needs the asynchronous magnetometer plus gyro-bias
transport, which is not established here and is reported open.

Nothing in this module promotes anything.  It consumes the fail-closed gate and
must report ``P4_MOTION_PASS`` exactly as the gate does.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

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
    """Prior-independent accelerometer posterior cap on transverse attitude variance.

    A scalar angle measurement with gain ``|f|`` and noise variance ``r`` gives
    ``P^+ = P^- r/(P^- |f|^2 + r)``, which is increasing in ``P^-`` with
    supremum ``r/|f|^2``.  The bound holds for every prior, so it is a
    reachability cap rather than a Riccati fixed point, and the declared
    Normal-Live invariant runs the accelerometer Joseph update at every valid
    IMU sample.  Only the two directions transverse to ``f`` are observed; the
    rotation about ``f`` is not.
    """
    noise = domain["configured_runtime"]["measurement_noise_std"]["accelerometer_mps2"]
    r_acc = min(map(float, noise))
    f_min = float(domain["normal_live"]["specific_force_norm_lower_mps2"])
    if not (r_acc > 0.0 and f_min > 0.0):
        raise RuntimeError("configured accelerometer std / specific-force lower lost positivity")
    if domain["normal_live"]["accelerometer_update_required_each_valid_imu_sample_after_live_entry"] is not True:
        raise RuntimeError("per-sample accelerometer Joseph invariant not declared")
    per_axis = up(up(r_acc * r_acc) / down(f_min * f_min))
    return {
        "argument": "P^+ = P^- r/(P^- |f|^2 + r) <= r/|f|^2 for every prior P^-",
        "prior_independent": True,
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
    # Same-cell innovation covariance on the transverse block: S = H P H^T + R.
    # At the cap H P H^T equals R exactly, so S^{-1} halves the residual energy
    # that the R^{-1} relaxation charges.
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
        "prior_independent_transverse_attitude_cap": cap,
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
    cap = d.get("prior_independent_transverse_attitude_cap", {})
    if cap.get("prior_independent") is not True:
        f.append("transverse attitude cap lost prior independence")
    if cap.get("yaw_about_specific_force_capped_here") is not False:
        f.append("yaw direction was claimed capped by the accelerometer event")
    x = float(cap.get("transverse_variance_upper_per_axis_rad2", -1.0))
    if not (math.isfinite(x) and x > 0.0):
        f.append("transverse attitude cap is not a positive finite bound")
    route = d.get("magnitude_route_after_transverse_cap", {})
    for k in ("transverse_directions_only", "third_attitude_direction_still_open"):
        if route.get(k) is not True:
            f.append("magnitude_route." + k + " not true")
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
        "transverse_cap_rad2": d["prior_independent_transverse_attitude_cap"]["transverse_trace_upper_two_axes_rad2"],
        "acc_ceiling_Rinv": d["magnitude_route_after_transverse_cap"]["accelerometer_ceiling_with_Rinverse_relaxation"],
        "acc_ceiling_Sinv": d["magnitude_route_after_transverse_cap"]["accelerometer_ceiling_retaining_same_cell_Sinverse"],
        "P4_MOTION_PASS": d["P4_MOTION_PASS"],
        "failures": f,
    }, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
