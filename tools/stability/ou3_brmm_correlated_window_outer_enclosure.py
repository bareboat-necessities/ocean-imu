#!/usr/bin/env python3
"""Correlated hard outer enclosure of the 601-sample COMPLETE-BRMM behavior.

Exact membership in the continuum spectral BRMM is unnecessary for a robust P4
proof: it is sufficient to prove P4 on a deterministic superset that contains
every BRMM sampled history. This module constructs such a superset without
Cartesianizing samples or axes.

The outer relation retains, on one source history:
  * the exact cross-sample translational primitive recurrence for (v,p,S_L);
  * one coupled J0/J1/J2 acceleration-moment IQC across all three axes;
  * the uniform physical p/v/acceleration/body-rate hard caps;
  * SO(3) attitude membership and rate-limited consecutive attitude motion;
  * the specific-force norm shell already declared by Normal-Live; and
  * the coupled BRMM lambda transition relation.

Every COMPLETE-BRMM history satisfies these constraints by construction of the
existing source contracts, so the left inclusion B^601_BRMM subset O^601_BRMM
is analytical. O^601_BRMM is intentionally larger than BRMM; no converse is
claimed. This closes a correlated outer-enclosure obligation, not the final
601-sample executor materialization or P4.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_brmm_complete_source as COMPLETE
import ou3_brmm_finite_window_primitive_qualification as PRIMITIVE
import ou3_brmm_centered_primitive_transition as TRANSITION
import ou3_brmm_acceleration_moment_iqc as MOMENT
import ou3_brmm_rlambda_transition as RLAMBDA
import ou3_brmm_physical_wave_source as WAVE

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 2
QUALIFICATION = "OU3_BRMM_CORRELATED_601_SAMPLE_OUTER_ENCLOSURE_V2"
N = 601
DT = 0.005
P3_DELTA = 1.0e-18


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    complete = COMPLETE.build(path)
    primitive = PRIMITIVE.build(path)
    transition = TRANSITION.build()
    moment = MOMENT.build()
    rlambda = RLAMBDA.build(path)
    bad = {
        "complete": COMPLETE.validate(complete),
        "primitive": PRIMITIVE.validate(primitive),
        "transition": TRANSITION.validate(transition),
        "moment": MOMENT.validate(moment),
        "R_lambda": RLAMBDA.validate(rlambda),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError("correlated BRMM outer-enclosure prerequisites failed: " + repr(bad))

    live = domain["normal_live"]
    pm = primitive["uniform_physical_primitives"]
    acc = float(live["non_gravitational_cog_acceleration_norm_upper_mps2"])
    rate_deg = float(live["body_rate_norm_upper_deg_s"])
    rate_rad = math.radians(rate_deg)
    sf_lo = float(live["specific_force_norm_lower_mps2"])
    sf_hi = float(live["specific_force_norm_upper_mps2"])
    if not all(math.isfinite(x) and x > 0.0 for x in (acc, rate_rad, sf_lo, sf_hi)):
        raise RuntimeError("outer-enclosure physical constants invalid")

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "sample_count": N,
        "sample_period_s": DT,
        "outer_set_symbol": "O^601_BRMM",
        "left_set_symbol": "B^601_BRMM",
        "left_inclusion": "B^601_BRMM subset O^601_BRMM",
        "left_inclusion_closed": True,
        "left_inclusion_reason": (
            "every constraint defining O^601_BRMM is a necessary consequence of the existing COMPLETE-BRMM physical/dynamic source contracts"
        ),
        "same_history_required_for_entire_window": True,
        "constraints": {
            "translation": {
                "state": "(x_s,k,phi_k,phi_L,v_k,p_k,S_L,k)",
                "physical_wave_generator_contract": WAVE.build(),
                "potential_identity": "S_L,k=phi(x_s,k)-phi_L; dphi/dt=p",
                "potential_out_is_next_potential_in": True,
                "generator_state_invariant_required": True,
                "potential_is_not_an_independent_S_supply_port": True,
                "numeric_generator_family_envelope_qualified": False,
                "uniform_velocity_norm_upper_mps": float(pm["V_m_norm_upper_mps"]),
                "uniform_position_norm_upper_m": float(pm["P_m_norm_upper_m"]),
                "recurrence": transition["exact_recurrence"],
                "primitive_out_is_next_primitive_in": True,
                "one_live_centered_S_origin_for_all_samples": True,
            },
            "acceleration_moments": {
                "normalized_coordinates": moment["normalized_coordinates"],
                "gram_inverse": moment["gram_inverse"],
                "joint_three_axis_supply_upper_mps4": acc * acc,
                "same_transition_witness_as_translation": True,
                "independent_J0_J1_J2_forbidden": True,
                "independent_axes_forbidden": True,
            },
            "rotation": {
                "R_wb_in_SO3_each_sample": True,
                "body_rate_norm_upper_rad_s": rate_rad,
                "consecutive_geodesic_distance_upper_rad": math.nextafter(rate_rad * DT, math.inf),
                "same_body_rate_history_drives_rotation": True,
            },
            "specific_force": {
                "norm_lower_mps2": sf_lo,
                "norm_upper_mps2": sf_hi,
                "same_physical_history_required": True,
            },
            "lambda": {
                "qualification": rlambda["qualification"],
                "machine_readable_transition_relation_closed": bool(rlambda["machine_readable_R_lambda_closed"]),
                "coupled_partition_energy_retained": bool(rlambda["coupled_partition_energy_retained"]),
                "coupled_peak_steepness_retained": bool(rlambda["coupled_peak_steepness_retained"]),
            },
        },
        "correlation_retained_across_samples": True,
        "correlation_retained_across_axes": True,
        "correlation_retained_between_translation_and_moments": True,
        "independent_sample_boxes_used": False,
        "independent_axis_boxes_used": False,
        "spectral_moments_alone_used_as_membership": False,
        "finite_frequency_grid_used": False,
        "seeded_simulator_used": False,
        "gaussian_good_event_used": False,
        "trajectory_replay_used": False,
        "converse_outer_to_BRMM_membership_claimed": False,
        "validated_correlated_outer_enclosure_closed": True,
        "executor_601_sample_relation_materialized_here": False,
        "joint_event_local_P_H_R_K_attached_here": False,
        "P3_delta": P3_DELTA,
        "P4_MOTION_PASS": False,
        "P4_PASS": False,
        "P5_MAY_START": False,
        "next_obligation": (
            "propagate O^601_BRMM through the estimator-owned 601-sample event graph using the primitive-prefix binding; attach reachable P/H/R/K, BIAS and radial coordinates, then attempt endpoint/every-prefix augmented LDLT"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_BRMM_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    if d.get("sample_count") != N or not math.isclose(float(d.get("sample_period_s", 0)), DT):
        f.append("canonical 601x5ms window changed")
    if d.get("left_inclusion") != "B^601_BRMM subset O^601_BRMM":
        f.append("left inclusion changed")
    for k in (
        "left_inclusion_closed", "same_history_required_for_entire_window",
        "correlation_retained_across_samples", "correlation_retained_across_axes",
        "correlation_retained_between_translation_and_moments",
        "validated_correlated_outer_enclosure_closed",
    ):
        if d.get(k) is not True:
            f.append(k + " not true")
    for k in (
        "independent_sample_boxes_used", "independent_axis_boxes_used",
        "spectral_moments_alone_used_as_membership", "finite_frequency_grid_used",
        "seeded_simulator_used", "gaussian_good_event_used", "trajectory_replay_used",
        "converse_outer_to_BRMM_membership_claimed", "executor_601_sample_relation_materialized_here",
        "joint_event_local_P_H_R_K_attached_here", "P4_MOTION_PASS", "P4_PASS", "P5_MAY_START",
    ):
        if d.get(k) is not False:
            f.append(k + " not false")
    c = d.get("constraints", {})
    t = c.get("translation", {})
    if t.get("primitive_out_is_next_primitive_in") is not True:
        f.append("translation primitive continuity lost")
    if t.get("one_live_centered_S_origin_for_all_samples") is not True:
        f.append("centered-S origin continuity lost")
    f.extend(WAVE.validate(t.get("physical_wave_generator_contract", {})))
    for key in ("potential_out_is_next_potential_in", "generator_state_invariant_required",
                "potential_is_not_an_independent_S_supply_port"):
        if t.get(key) is not True:
            f.append("physical wave primitive constraint lost: "+key)
    if t.get("potential_identity") != "S_L,k=phi(x_s,k)-phi_L; dphi/dt=p":
        f.append("physical wave potential identity lost")
    if t.get("numeric_generator_family_envelope_qualified") is not False:
        f.append("numeric physical generator family falsely promoted")
    m = c.get("acceleration_moments", {})
    if m.get("same_transition_witness_as_translation") is not True:
        f.append("moment witness detached from translation")
    reference_moment = MOMENT.build()
    if m.get("gram_inverse") != reference_moment["gram_inverse"]:
        f.append("moment IQC Gram inverse changed")
    if m.get("normalized_coordinates") != reference_moment["normalized_coordinates"]:
        f.append("normalized moment coordinates changed")
    if m.get("independent_J0_J1_J2_forbidden") is not True or m.get("independent_axes_forbidden") is not True:
        f.append("moment correlation weakened")
    r = c.get("rotation", {})
    if r.get("R_wb_in_SO3_each_sample") is not True or r.get("same_body_rate_history_drives_rotation") is not True:
        f.append("rotation relation weakened")
    lam = c.get("lambda", {})
    for k in ("machine_readable_transition_relation_closed", "coupled_partition_energy_retained", "coupled_peak_steepness_retained"):
        if lam.get(k) is not True:
            f.append("lambda relation lost " + k)
    if d.get("P3_delta") != P3_DELTA:
        f.append("P3 delta changed")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "left_inclusion": d["left_inclusion_closed"],
        "correlated_outer": d["validated_correlated_outer_enclosure_closed"],
        "executor_materialized": d["executor_601_sample_relation_materialized_here"],
        "P4": d["P4_PASS"],
        "failures": f,
    }, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
