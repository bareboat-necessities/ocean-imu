#!/usr/bin/env python3
"""Compact sampled behavior set for the complete SEA3 finite window.

This module does not create a bounded-input source.  It gives a precise
finite-dimensional target set for the numerical SEA0 enclosure:

    B_SEA3^N = closure { sampled joint physical outputs of admitted complete
                         SEA3 realizations on N=601 samples }.

Membership in the set still requires one common SEA3 spectral/phase/response
witness, the coupled lambda transition relation, and the existing Normal-Live
conditions.  The Normal-Live acceleration/body-rate caps are used only to prove
that the sampled projection is bounded.  They are not sufficient membership
conditions and may not be used to generate arbitrary sample sequences.

Because the sampled physical projection is finite dimensional, its closure is
closed; the existing hard Normal-Live bounds make it bounded; Heine--Borel then
gives compactness.  This closes an important topological part of the hard
finite-window constraint without inventing a spectral-amplitude constant.

What remains open is computational: a validated membership/separation or outer-
enclosure oracle that preserves the common SEA3 witness.  Until that exists,
this compact behavior set cannot be serialized as independent sample boxes and
cannot open the canonical provider/P3 gate.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_continuum_phase_state as PHASE
import ou3_sea3_rlambda_transition as RLAMBDA

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_COMPACT_SAMPLED_WINDOW_BEHAVIOR_V1"
N = 601


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    complete = COMPLETE.build(path)
    phase = PHASE.build()
    rlambda = RLAMBDA.build(path)
    bad = {
        "complete": COMPLETE.validate(complete),
        "phase": PHASE.validate(phase),
        "R_lambda": RLAMBDA.validate(rlambda),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"SEA3 hard-window behavior prerequisites failed: {bad}")

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

    # These are finite-dimensional sampled *output* coordinates.  The phase
    # field itself remains continuum and is referenced by a common witness; it
    # is not discretized into these coordinates.
    per_sample_scalar_coordinates = 3 + 3 + 3 + 9
    sampled_projection_dimension = N * per_sample_scalar_coordinates

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "sample_count": N,
        "sampled_projection_coordinates": [
            "f_cog_body[3]",
            "omega_body_corrected[3]",
            "specific_force_body[3]",
            "R_wb[3x3]",
        ],
        "sampled_projection_dimension": sampled_projection_dimension,
        "behavior_set_symbol": "B^601_SEA3",
        "behavior_set_definition": (
            "closure of sampled joint physical outputs of admitted complete SEA3 realizations carrying one common (x^s,lambda,response) witness"
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
            "B^601_SEA3 is closed by definition and bounded by the existing Normal-Live physical caps and compact SO(3) rotation coordinate; in finite-dimensional sampled output space Heine-Borel gives compactness"
        ),
        "sampled_behavior_set_compact": True,
        "SEA3_parameter_domain_compact": True,
        "continuum_phase_coordinate_set_closed": phase[
            "continuum_phase_coordinate_set_closed"
        ],
        "phase_continuous_propagation_closed": phase[
            "phase_continuous_propagation_closed"
        ],
        "machine_readable_R_lambda_closed": rlambda[
            "machine_readable_R_lambda_closed"
        ],
        "membership_requires_common_SEA3_witness": True,
        "normal_live_caps_are_membership_sufficient": False,
        "arbitrary_bounded_sequence_is_member": False,
        "independent_sample_boxes_define_behavior_set": False,
        "independent_axis_boxes_define_behavior_set": False,
        "finite_frequency_grid_used": False,
        "seeded_simulator_used": False,
        "gaussian_good_event_used": False,
        "spectral_moments_alone_used_as_membership": False,
        "validated_membership_or_separation_oracle_closed": False,
        "validated_correlated_outer_enclosure_closed": False,
        "provider_artifact_materialized_here": False,
        "P3_promoted": False,
        "next_obligation": (
            "construct a validated correlated membership/separation or outer-enclosure oracle for B^601_SEA3 that preserves the common continuum SEA3 witness; do not replace B^601_SEA3 by its per-sample norm hull"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for key in (
        "closure_is_part_of_definition",
        "finite_dimensional_sample_projection",
        "sampled_behavior_set_compact",
        "SEA3_parameter_domain_compact",
        "continuum_phase_coordinate_set_closed",
        "phase_continuous_propagation_closed",
        "machine_readable_R_lambda_closed",
        "membership_requires_common_SEA3_witness",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "normal_live_caps_are_membership_sufficient",
        "arbitrary_bounded_sequence_is_member",
        "independent_sample_boxes_define_behavior_set",
        "independent_axis_boxes_define_behavior_set",
        "finite_frequency_grid_used",
        "seeded_simulator_used",
        "gaussian_good_event_used",
        "spectral_moments_alone_used_as_membership",
        "validated_membership_or_separation_oracle_closed",
        "validated_correlated_outer_enclosure_closed",
        "provider_artifact_materialized_here",
        "P3_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("behavior set detached from canonical SEA3 source")
    if d.get("behavior_set_symbol") != "B^601_SEA3":
        f.append("behavior set symbol changed")
    if int(d.get("sample_count", 0)) != N:
        f.append("behavior set does not cover canonical 601-sample window")
    witness = d.get("boundedness_witness", {})
    for key in (
        "non_gravitational_cog_acceleration_norm_upper_mps2",
        "body_rate_norm_upper_deg_s",
        "body_rate_norm_upper_rad_s",
        "specific_force_norm_lower_mps2",
        "specific_force_norm_upper_mps2",
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
        "oracle_closed": d["validated_membership_or_separation_oracle_closed"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
