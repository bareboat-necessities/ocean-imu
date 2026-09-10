#!/usr/bin/env python3
"""Fail-closed gate for the end-to-end shipping stability theorem.

The target is not "P4 passes".  It is

    for every admitted BRMM motion, every admitted BIAS0/BIAS1/BIAS2 history and
    every admitted sensor disturbance, the shipping implementation started from
    the declared startup uncertainty reaches the certified Live basin in finite
    time through its actual Mahony/proxy startup and then stays there with
    practical ISS,

which decomposes as

    startup uncertainty --[P5 capture, finite time]--> R_P4 --[P4]--> practical ISS,

with P3 supplying the covariance/observability properties along the admitted
execution rather than being the stability statement itself.

This gate holds the composition.  It never invents a stage: every obligation is
either discharged by a named producer or listed as open with its failure class.
The three groups are

  * P5 capture, from the deployed startup state machine
    ``Cold -> TunerWarm -> TunerReady -> (external bootstrap goLive) -> Live``;
  * P4 invariance and practical ISS on the certified basin, mirrored from
    ``ou3_p4_final_closure_gate``;
  * transition coverage, so no mathematical teleport from startup into an
    assumed Live state can be counted as a proof.

The basin is deliberately handled as a MAXIMISATION.  Width certified in P4 is
width the capture stage does not have to achieve, so the gate carries the
frontier from ``ou3_p4_basin_frontier`` and the deployed Live-entry pinning from
``ou3_p5_live_entry_reachability`` rather than a single convenient radius.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_p4_basin_frontier as FRONTIER
import ou3_p4_correlated_entry_relation as ENTRYREL
import ou3_p4_final_closure_gate as P4GATE
import ou3_p5_live_entry_reachability as P5REACH

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
QUALIFICATION = "OU3_END_TO_END_SHIPPING_STABILITY_GATE_V1"

THEOREM = (
    "for all admitted BRMM motion, BIAS0/BIAS1/BIAS2 histories and sensor disturbances: "
    "declared startup uncertainty --[Mahony/proxy capture, finite time]--> certified Live "
    "basin --[practical ISS]--> bounded steady operation")

# Deployed stages, in execution order.  A proof may not skip one.
DEPLOYED_STAGES = ("Cold", "TunerWarm", "TunerReady", "goLive", "Live_held_bias_H18",
                   "magnetic_lock_or_regauge", "H18_to_A21_release", "Live_A21")

# Shipping event families that a covering proof has to carry inside Live.
DEPLOYED_LIVE_EVENTS = ("prediction", "aw_floor", "S_zero", "accelerometer",
                        "magnetometer", "bias_projection", "finite_reset")

# Stages the P4 invariance half owns rather than the capture half: once the
# estimator is inside the basin, steady A21 operation is the practical-ISS
# statement, not a capture step.
P4_COVERED_STAGES = ("Live_A21",)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain = json.loads(Path(domain_path).read_text(encoding="utf-8"))
    p4 = P4GATE.build()
    p4_failures = P4GATE.validate(p4)
    if p4_failures:
        raise RuntimeError("P4 gate invalid: " + repr(p4_failures))
    reach = P5REACH.build(Path(domain_path))
    reach_failures = P5REACH.validate(reach)
    if reach_failures:
        raise RuntimeError("live-entry reachability invalid: " + repr(reach_failures))
    frontier = FRONTIER.build()
    frontier_failures = FRONTIER.validate(frontier)
    if frontier_failures:
        raise RuntimeError("basin frontier invalid: " + repr(frontier_failures))
    rel = ENTRYREL.build(Path(domain_path))
    rel_failures = ENTRYREL.validate(rel)
    if rel_failures:
        raise RuntimeError("correlated entry relation invalid: " + repr(rel_failures))

    warmup_s = float(domain["startup"]["online_tune_warmup_sec"])
    usable_floor = domain["startup"]["wave_period_estimator_usable_time_floor_definition"]

    capture = {
        "cold_to_tuner_warm_finite_time": {
            "stage": "Cold -> TunerWarm",
            "covers_stages": ["Cold"],
            "established": True,
            "why": ("the deployed transition is a pure elapsed-time comparison against "
                    "online_tune_warmup_sec on the admitted sample cadence"),
            "bound_s": warmup_s,
            "class": None,
        },
        "tuner_warm_to_tuner_ready_finite_time": {
            "stage": "TunerWarm -> TunerReady",
            "covers_stages": ["TunerWarm"],
            "established": False,
            "why": ("needs a finite-time attainment proof of tuner freq-ready, tuner-ready "
                    "and WavePeriodEstimator usable-period over every admitted BRMM history; "
                    "the declared timing floor is a lower bound, not attainment"),
            "bound_s": None,
            "declared_timing_floor": usable_floor,
            "class": "E",
        },
        "bootstrap_tilt_and_north_acquisition": {
            "stage": "TunerReady -> goLive",
            "covers_stages": ["TunerReady"],
            "established": False,
            "why": ("goLive is called by the external bootstrap once it has tilt and north; "
                    "Mahony/proxy tilt convergence and the magnetic gauge event both need "
                    "finite-time certificates over the admitted source class"),
            "bound_s": None,
            "class": "E",
        },
        "proxy_attitude_radius_at_go_live": {
            "stage": "goLive",
            "covers_stages": ["goLive"],
            "established": False,
            "why": ("the certified Cayley radius handed to initialize_from_attitude is the "
                    "composed tilt/yaw bound of the P5 staged-capture section; it is not "
                    "instantiated against the deployed proxy"),
            "bound_s": None,
            "class": "E",
        },
        "live_entry_state_lies_in_certified_basin": {
            "stage": "goLive -> Live",
            "covers_stages": ["Live_held_bias_H18"],
            "established": False,
            "why": ("the deployed pinning makes every non-attitude entry coordinate the "
                    "negated physical truth, so the radii are the BRMM bounded-motion "
                    "primitives; V_m, P_m and S_m are declared but not instantiated"),
            "bound_s": None,
            "uninstantiated_primitives": list(reach["uninstantiated_BRMM_primitives"]),
            "class": "E",
        },
        "held_bias_to_A21_release_covered": {
            "stage": "H18 -> A21",
            "covers_stages": ["magnetic_lock_or_regauge", "H18_to_A21_release"],
            "established": False,
            "why": ("the accelerometer-bias unlock after the deployed magnetometer count is a "
                    "hybrid transition; it is outside the single-mode word certificates"),
            "bound_s": None,
            "class": "E",
        },
    }
    capture_pass = all(v["established"] for v in capture.values())

    invariance = {
        "P4_motion_practical_ISS": {
            "established": bool(p4["P4_MOTION_PASS"]),
            "open_blockers": list(p4["remaining_mathematical_P4_blockers"]),
        },
        "conditional_binary32_additive_ISS": {
            "established": bool(p4["conditional_full_shipping_finite_precision_additive_ISS_closed"]),
            "open_blockers": [],
        },
        "all_bias_families_closed": {
            "established": bool(p4["all_required_bias_families_closed"]),
            "open_blockers": [],
        },
    }
    invariance_pass = all(v["established"] for v in invariance.values())

    covered_by_capture = {stage for row in capture.values() for stage in row["covers_stages"]}
    unknown_stage = covered_by_capture - set(DEPLOYED_STAGES)
    if unknown_stage:
        raise RuntimeError("capture obligation names an undeclared stage: " + repr(sorted(unknown_stage)))
    coverage = {
        "deployed_stages_enumerated": list(DEPLOYED_STAGES),
        "deployed_live_events_enumerated": list(DEPLOYED_LIVE_EVENTS),
        "startup_to_live_teleport_permitted": False,
        "live_entry_pinning_from_deployed_code": bool(reach["reachability_argument_is_from_deployed_code"]),
        "linear_block_unpropagated_before_go_live": not bool(
            reach["mekf_linear_block_propagated_before_go_live"]),
        "attitude_linear_cross_zeroed_at_entry": bool(
            reach["attitude_linear_cross_covariance_zeroed_at_go_live"]),
        "stages_covered_by_capture_obligations": sorted(covered_by_capture),
        "stages_covered_by_invariance": sorted(P4_COVERED_STAGES),
        "every_stage_has_a_named_obligation": bool(
            covered_by_capture | set(P4_COVERED_STAGES) == set(DEPLOYED_STAGES)),
    }

    basin = {
        "objective": "maximise the certified basin, never shrink it for proof convenience",
        "declared_box_chart_overshoot_factor_H18": frontier["modes"]["H18"]["declared_box_overshoot_factor"],
        "max_uniform_scale_of_declared_box": frontier["max_uniform_scale_over_modes"],
        "max_volume_radii_H18": frontier["modes"]["H18"]["max_volume_radii"],
        "max_integral_radius_with_others_declared_m_s": frontier["max_integral_radius_with_others_declared_m_s"],
        "correlated_one_cadence_integral_radius_m_s": rel["correlated_radius_one_cadence_m_s"],
        "frontier_is_a_frozen_map_diagnostic": bool(
            frontier["frozen_map_binary64_diagnostic_not_outward_certificate"]),
        "entry_radii_reduced_for_proof_convenience": False,
    }

    end_to_end = bool(capture_pass and invariance_pass)
    return {
        "qualification": QUALIFICATION,
        "theorem": THEOREM,
        "canonical_source": p4["canonical_source"],
        "P3_delta": p4["P3_delta"],
        "filter_changed": False,
        "quality_gates_changed": False,
        "declared_domain_shrunk": False,

        "P5_capture_obligations": capture,
        "P5_CAPTURE_PASS": capture_pass,
        "P4_invariance_obligations": invariance,
        "P4_INVARIANCE_PASS": invariance_pass,
        "transition_coverage": coverage,
        "certified_basin": basin,

        "P4_MOTION_PASS": bool(p4["P4_MOTION_PASS"]),
        "P4_PASS": bool(p4["P4_PASS"]),
        "P5_MAY_START": bool(p4["P5_MAY_START"]),
        "END_TO_END_STABILITY_PASS": end_to_end,
        "END_TO_END_DEPLOYMENT_PASS": bool(
            end_to_end and p4["target_toolchain_finite_precision_qualified"]),

        "remaining_P5_capture_obligations": [k for k, v in capture.items() if not v["established"]],
        "remaining_P4_mathematical_blockers": list(p4["remaining_mathematical_P4_blockers"]),
        "remaining_deployment_blockers": list(p4["remaining_deployment_blockers"]),
    }


def validate(d: dict) -> list[str]:
    f = []
    if d.get("qualification") != QUALIFICATION:
        f.append("qualification mismatch")
    if d.get("theorem") != THEOREM:
        f.append("theorem statement changed")
    if float(d.get("P3_delta", 0.0)) != 1e-18:
        f.append("P3 delta changed")
    for k in ("filter_changed", "quality_gates_changed", "declared_domain_shrunk"):
        if d.get(k) is not False:
            f.append(k + " not false")
    cov = d.get("transition_coverage", {})
    if tuple(cov.get("deployed_stages_enumerated", ())) != DEPLOYED_STAGES:
        f.append("deployed stage list changed")
    if tuple(cov.get("deployed_live_events_enumerated", ())) != DEPLOYED_LIVE_EVENTS:
        f.append("deployed live event list changed")
    if cov.get("startup_to_live_teleport_permitted") is not False:
        f.append("startup-to-Live teleport permitted")
    for k in ("live_entry_pinning_from_deployed_code",
              "linear_block_unpropagated_before_go_live",
              "attitude_linear_cross_zeroed_at_entry",
              "every_stage_has_a_named_obligation"):
        if cov.get(k) is not True:
            f.append("transition_coverage." + k + " not true")
    capture = d.get("P5_capture_obligations", {})
    if not capture:
        f.append("no P5 capture obligations enumerated")
    for name, row in capture.items():
        if not row.get("covers_stages"):
            f.append(name + " covers no deployed stage")
        if any(s not in DEPLOYED_STAGES for s in row.get("covers_stages", ())):
            f.append(name + " names an undeclared stage")
        if row.get("established") is not True and row.get("class") not in ("A", "B", "C", "D", "E"):
            f.append(name + " is open without a declared failure class")
        if row.get("established") is True and row.get("class") is not None:
            f.append(name + " is established but still carries a failure class")
    if bool(d.get("P5_CAPTURE_PASS")) != all(bool(v["established"]) for v in capture.values()):
        f.append("P5 capture promotion is not the exact conjunction")
    inv = d.get("P4_invariance_obligations", {})
    if bool(d.get("P4_INVARIANCE_PASS")) != all(bool(v["established"]) for v in inv.values()):
        f.append("P4 invariance promotion is not the exact conjunction")
    expected = bool(d.get("P5_CAPTURE_PASS")) and bool(d.get("P4_INVARIANCE_PASS"))
    if bool(d.get("END_TO_END_STABILITY_PASS")) != expected:
        f.append("end-to-end promotion is not the exact conjunction")
    if d.get("END_TO_END_STABILITY_PASS") and not d.get("P4_PASS"):
        f.append("end-to-end passed while P4 is false")
    # END_TO_END_DEPLOYMENT_PASS is the strongest bit in the payload, so it needs
    # its own conjunction check rather than riding on the stability one.
    if d.get("END_TO_END_DEPLOYMENT_PASS") and not d.get("END_TO_END_STABILITY_PASS"):
        f.append("deployment pass set while the end-to-end theorem is open")
    if d.get("END_TO_END_DEPLOYMENT_PASS") and d.get("remaining_deployment_blockers"):
        f.append("deployment pass set while deployment blockers remain")
    # The bits mirrored from the P4 gate must stay consistent with the blocker
    # list they are mirrored alongside.
    if d.get("P4_PASS") and d.get("remaining_P4_mathematical_blockers"):
        f.append("P4 pass set while mathematical blockers remain")
    if d.get("P4_MOTION_PASS") and d.get("remaining_P4_mathematical_blockers"):
        f.append("P4 motion pass set while mathematical blockers remain")
    if d.get("P5_MAY_START") and not d.get("P4_PASS"):
        f.append("P5 may start while P4 is false")
    if bool(d.get("P4_MOTION_PASS")) != bool(
            d.get("P4_invariance_obligations", {}).get("P4_motion_practical_ISS", {}).get("established")):
        f.append("mirrored P4_MOTION_PASS disagrees with its own invariance row")
    basin = d.get("certified_basin", {})
    if basin.get("entry_radii_reduced_for_proof_convenience") is not False:
        f.append("basin was reduced for proof convenience")
    if basin.get("frontier_is_a_frozen_map_diagnostic") is not True:
        f.append("basin frontier was promoted past its diagnostic status")
    if float(basin.get("max_uniform_scale_of_declared_box", 0.0)) <= 0.0:
        f.append("basin frontier carries no positive uniform scale")
    # Open obligation lists must agree with their own tables.
    if set(d.get("remaining_P5_capture_obligations", ())) != {
            k for k, v in capture.items() if not v["established"]}:
        f.append("remaining P5 obligation list disagrees with the table")
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
        "END_TO_END_STABILITY_PASS": d["END_TO_END_STABILITY_PASS"],
        "P5_CAPTURE_PASS": d["P5_CAPTURE_PASS"],
        "P4_INVARIANCE_PASS": d["P4_INVARIANCE_PASS"],
        "remaining_P5": d["remaining_P5_capture_obligations"],
        "remaining_P4": d["remaining_P4_mathematical_blockers"],
        "basin": d["certified_basin"]["max_volume_radii_H18"],
        "failures": f,
    }, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
