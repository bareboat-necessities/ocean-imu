#!/usr/bin/env python3
"""Compact sampled behavior set for the complete BRMM finite window.

B^601_BRMM is the closure of sampled joint physical outputs of admitted
complete-BRMM realizations carrying one common spectral/phase/response witness.
The set is compact in its finite-dimensional sampled projection.

A robust proof does not require an exact membership oracle if it instead uses a
validated deterministic correlated superset.  The companion
``ou3_brmm_correlated_window_outer_enclosure`` now provides O^601_BRMM with
B^601_BRMM subset O^601_BRMM, retaining cross-sample primitive recurrence,
all-axis acceleration-moment IQC, rate/SO(3) constraints and the coupled lambda
relation.  This closes the correlated-outer-enclosure route while exact BRMM
membership/separation and final provider/executor materialization remain open.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_brmm_complete_source as COMPLETE
import ou3_brmm_continuum_phase_state as PHASE
import ou3_brmm_rlambda_transition as RLAMBDA
import ou3_brmm_correlated_window_outer_enclosure as OUTER

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 2
QUALIFICATION = "OU3_BRMM_COMPACT_SAMPLED_WINDOW_BEHAVIOR_V2"
N = 601


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    complete = COMPLETE.build(path)
    phase = PHASE.build()
    rlambda = RLAMBDA.build(path)
    outer = OUTER.build(path)
    bad = {
        "complete": COMPLETE.validate(complete),
        "phase": PHASE.validate(phase),
        "R_lambda": RLAMBDA.validate(rlambda),
        "outer": OUTER.validate(outer),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"BRMM hard-window behavior prerequisites failed: {bad}")

    live = domain["normal_live"]
    acc = float(live["non_gravitational_cog_acceleration_norm_upper_mps2"])
    rate_deg = float(live["body_rate_norm_upper_deg_s"])
    rate_rad = math.radians(rate_deg)
    sf_lo = float(live["specific_force_norm_lower_mps2"])
    sf_hi = float(live["specific_force_norm_upper_mps2"])
    if not all(math.isfinite(x) and x > 0.0 for x in (acc, rate_deg, rate_rad, sf_lo, sf_hi)):
        raise RuntimeError("Normal-Live hard bounds are not positive finite")
    if sf_lo > sf_hi:
        raise RuntimeError("specific-force norm bounds are reversed")

    per_sample_scalar_coordinates = 3 + 3 + 3 + 9
    sampled_projection_dimension = N * per_sample_scalar_coordinates
    outer_closed = bool(
        outer["left_inclusion_closed"]
        and outer["validated_correlated_outer_enclosure_closed"]
        and outer["same_history_required_for_entire_window"]
        and not outer["independent_sample_boxes_used"]
        and not outer["independent_axis_boxes_used"]
    )

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "sample_count": N,
        "sampled_projection_coordinates": [
            "f_cog_body[3]", "omega_body_corrected[3]",
            "specific_force_body[3]", "R_wb[3x3]",
        ],
        "sampled_projection_dimension": sampled_projection_dimension,
        "behavior_set_symbol": "B^601_BRMM",
        "behavior_set_definition": (
            "closure of sampled joint physical outputs of admitted complete BRMM realizations carrying one common (x^s,lambda,response) witness"
        ),
        "closure_is_part_of_definition": True,
        "finite_dimensional_sample_projection": True,
        "boundedness_witness": {
            "non_gravitational_cog_acceleration_norm_upper_mps2": acc,
            "body_rate_norm_upper_deg_s": rate_deg,
            "body_rate_norm_upper_rad_s": rate_rad,
            "specific_force_norm_lower_mps2": sf_lo,
            "specific_force_norm_upper_mps2": sf_hi,
            "R_wb_in_SO3": True,
        },
        "compactness_argument": (
            "B^601_BRMM is closed by definition and bounded by Normal-Live hard caps and compact SO(3); finite-dimensional Heine-Borel gives compactness"
        ),
        "sampled_behavior_set_compact": True,
        "reference_parameter_domain_compact": True,
        "continuum_phase_coordinate_set_closed": phase["continuum_phase_coordinate_set_closed"],
        "phase_continuous_propagation_closed": phase["phase_continuous_propagation_closed"],
        "machine_readable_R_lambda_closed": rlambda["machine_readable_R_lambda_closed"],
        "membership_requires_common_BRMM_witness": True,
        "normal_live_caps_are_membership_sufficient": False,
        "arbitrary_bounded_sequence_is_member": False,
        "independent_sample_boxes_define_behavior_set": False,
        "independent_axis_boxes_define_behavior_set": False,
        "finite_frequency_grid_used": False,
        "seeded_simulator_used": False,
        "gaussian_good_event_used": False,
        "spectral_moments_alone_used_as_membership": False,
        "validated_membership_or_separation_oracle_closed": False,
        "validated_correlated_outer_enclosure_closed": outer_closed,
        "correlated_outer_enclosure_qualification": outer["qualification"],
        "correlated_outer_set_symbol": outer["outer_set_symbol"],
        "correlated_outer_left_inclusion_closed": outer["left_inclusion_closed"],
        "correlated_outer_retains_cross_sample_and_axis_dependence": bool(
            outer["correlation_retained_across_samples"]
            and outer["correlation_retained_across_axes"]
            and outer["correlation_retained_between_translation_and_moments"]
        ),
        "provider_artifact_materialized_here": False,
        "P3_promoted": False,
        "next_obligation": (
            "the correlated outer-enclosure route is closed; propagate O^601_BRMM through the 601-sample estimator/Riccati event graph and primitive-prefix binding, then materialize the provider artifact"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for key in (
        "closure_is_part_of_definition", "finite_dimensional_sample_projection",
        "sampled_behavior_set_compact", "reference_parameter_domain_compact",
        "continuum_phase_coordinate_set_closed", "phase_continuous_propagation_closed",
        "machine_readable_R_lambda_closed", "membership_requires_common_BRMM_witness",
        "validated_correlated_outer_enclosure_closed",
        "correlated_outer_left_inclusion_closed",
        "correlated_outer_retains_cross_sample_and_axis_dependence",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "normal_live_caps_are_membership_sufficient", "arbitrary_bounded_sequence_is_member",
        "independent_sample_boxes_define_behavior_set", "independent_axis_boxes_define_behavior_set",
        "finite_frequency_grid_used", "seeded_simulator_used", "gaussian_good_event_used",
        "spectral_moments_alone_used_as_membership",
        "validated_membership_or_separation_oracle_closed",
        "provider_artifact_materialized_here", "P3_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("canonical_source") != "COMPLETE_BRMM_NORMAL_LIVE_WORD":
        f.append("behavior set detached from canonical BRMM source")
    if d.get("behavior_set_symbol") != "B^601_BRMM":
        f.append("behavior set symbol changed")
    if d.get("correlated_outer_set_symbol") != "O^601_BRMM":
        f.append("correlated outer set symbol changed")
    if int(d.get("sample_count", 0)) != N:
        f.append("behavior set does not cover canonical 601-sample window")
    witness = d.get("boundedness_witness", {})
    for key in (
        "non_gravitational_cog_acceleration_norm_upper_mps2",
        "body_rate_norm_upper_deg_s", "body_rate_norm_upper_rad_s",
        "specific_force_norm_lower_mps2", "specific_force_norm_upper_mps2",
    ):
        x = float(witness.get(key, math.nan))
        if not (math.isfinite(x) and x > 0.0):
            f.append(f"invalid boundedness witness {key}")
    if witness.get("R_wb_in_SO3") is not True:
        f.append("SO(3) boundedness/compactness witness disappeared")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "behavior_set": d["behavior_set_symbol"],
        "dimension": d["sampled_projection_dimension"],
        "compact": d["sampled_behavior_set_compact"],
        "membership_oracle_closed": d["validated_membership_or_separation_oracle_closed"],
        "correlated_outer_closed": d["validated_correlated_outer_enclosure_closed"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
