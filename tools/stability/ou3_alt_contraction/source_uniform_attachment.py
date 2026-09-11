#!/usr/bin/env python3
"""Falsifiable ALT attachment audit for the actual source-uniform shipping word.

This consumes existing COMPLETE-BRMM/JOINT event machinery read-only. It asks
which ingredients required by the ALT finite-increment joint24 word are present,
and refuses to promote a differential or single-execution cell into a paired
finite-increment relation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT / "tools/stability"), str(ROOT)]

import ou3_p4_source_uniform_estimator_event_attachment as ATTACH
import ou3_brmm_complete_source as COMPLETE

QUALIFICATION = "ALT_SOURCE_UNIFORM_FINITE_INCREMENT_ATTACHMENT_AUDIT_V1"


def build() -> dict:
    source = COMPLETE.build()
    sf = COMPLETE.validate(source)
    if sf:
        raise RuntimeError(f"COMPLETE-BRMM source invalid: {sf}")

    # Execute the event-local relation without importing the old route's startup
    # theorem. This proves the mechanics needed by ALT are callable.
    smoke = ATTACH._smoke()
    event_local_closed = bool(
        smoke["joint_successors_retained"] > 0
        and smoke["all_H_cells_bound"]
        and smoke["all_A_cells_bound"]
        and smoke["all_event_orders_equal"]
        and smoke["authoritative_next_frontend_is_joint"]
    )

    # Deliberately try the inherited theorem-facing builder. On current main it
    # reaches into the separate Mahony startup invariant. ALT Live-word
    # attachment must not silently acquire that prerequisite.
    inherited = {"builder_completed": False, "startup_coupled": False, "error": None}
    try:
        d = ATTACH.build()
    except RuntimeError as exc:
        msg = str(exc)
        inherited.update(
            error=msg,
            startup_coupled=(
                "continuous Mahony invariant invalid" in msg
                or "continuous_all_live_PI_invariant_closed" in msg
                or "initial seed angle not closed" in msg
            ),
        )
    else:
        inherited.update(
            builder_completed=True,
            relation_closed=bool(
                d.get("source_uniform_estimator_owned_event_attachment_relation_closed")
            ),
        )

    paired_fields = {
        "paired_admissible_predecessor_states": False,
        "paired_covariance_states_P0_P1": False,
        "paired_frontend_tuner_states": False,
        "paired_guard_outcomes_or_hard_guard_graph": False,
        "paired_measurement_residuals_r0_r1": False,
        "paired_gain_numerators_N0_N1": False,
        "paired_innovations_S0_S1": False,
        "nominal_solve_variable_q0_bound_to_S0_q0_eq_r0": False,
        "deltaN_and_deltaS_bound_to_same_paired_history": False,
        "thin_product_ports_uN_uS_have_hard_same_history_product_graphs": False,
    }

    return {
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "event_local_attachment_executes_without_startup_theorem": event_local_closed,
        "event_local_smoke": smoke,
        "inherited_theorem_builder": inherited,
        "ALT_must_not_inherit_startup_for_Live_word_attachment": True,
        "available_actual_shipping_relations": {
            "same_history_estimator_owned_event_cells_executable": event_local_closed,
            "correlated_P_H_R_same_cell": True,
            "actual_applied_anisotropic_R_S": True,
            "frontend_tuner_scheduler_ancestry": True,
            "H18_and_A21_literal_event_order": True,
            "A21_same_history_true_bias_coordinate": True,
            "physical_wave_payload_in_literal_event_cells": True,
        },
        "paired_finite_increment_requirements": paired_fields,
        "paired_finite_increment_relation_closed": False,
        "rank3_structure_safe_to_use": {
            "measurement_dimension": 3,
            "inverse_free_form": "S*q=r, correction=N*q",
            "thin_increment_ports": "u_S=delta_S*q0; u_N=delta_N*q0",
            "thin_master_port_dimension_joint24": 33,
            "old_vec_increment_dimension_joint24": 87,
            "storage_measurement_delta_rank_upper": 6,
            "full_state_reduction_permitted": False,
            "A21_motion_bias_cross_terms_retained": True,
        },
        "ALT_ACTUAL_SOURCE_UNIFORM_FINITE_INCREMENT_WORD_ATTACHED": False,
        "ALT_LIVE_PASS": False,
        "P4_promoted": False,
        "classification": (
            "proof-method/representation attachment failure plus unwanted startup coupling"
        ),
        "limiter": (
            "event-local same-history H/A attachment executes, but the theorem-facing "
            "inherited builder is startup-coupled and ALT still lacks a paired coefficient "
            "graph tying r0/r1, N0/N1, S0/S1, q0/q1 and thin products u_N/u_S to one "
            "paired COMPLETE-BRMM history"
        ),
        "does_not_invalidate": [
            "inverse-free finite-increment algebra",
            "joint24 bounded-neutral-supply target",
            "rank-three thin-factor arithmetic",
            "independent old P2/P3/P4/P5 route",
        ],
        "next_falsifiable_experiment": (
            "construct a Live-only paired JOINT+kernel transition that takes two estimator "
            "states over one admitted COMPLETE-BRMM source continuation, emits paired literal "
            "event cells, and enforces u_S=delta_S*q0 and u_N=delta_N*q0 as hard same-history "
            "product graphs; keep startup separate"
        ),
    }


def validate(d: dict) -> list[str]:
    f = []
    if d.get("qualification") != QUALIFICATION:
        f.append("qualification mismatch")
    if d.get("event_local_attachment_executes_without_startup_theorem") is not True:
        f.append("event-local Live attachment did not execute")
    inherited = d.get("inherited_theorem_builder", {})
    if inherited.get("builder_completed") is False and inherited.get("startup_coupled") is not True:
        f.append("inherited builder failed for an unclassified reason")
    if d.get("paired_finite_increment_relation_closed") is not False:
        f.append("paired finite-increment relation must remain fail-closed")
    if d.get("ALT_ACTUAL_SOURCE_UNIFORM_FINITE_INCREMENT_WORD_ATTACHED") is not False:
        f.append("ALT actual word promoted without paired relation")
    if d.get("ALT_LIVE_PASS") is not False or d.get("P4_promoted") is not False:
        f.append("proof gate promoted")
    r = d.get("rank3_structure_safe_to_use", {})
    if r.get("measurement_dimension") != 3 or r.get("full_state_reduction_permitted") is not False:
        f.append("rank-three structure mischaracterized")
    if not int(r.get("thin_master_port_dimension_joint24", 999)) < int(
        r.get("old_vec_increment_dimension_joint24", 0)
    ):
        f.append("thin product-port representation did not reduce lifted dimension")
    return f


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build()
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "event_local": d["event_local_attachment_executes_without_startup_theorem"],
                "startup_coupled": d["inherited_theorem_builder"].get("startup_coupled"),
                "paired": d["paired_finite_increment_relation_closed"],
                "limiter": d["limiter"],
                "failures": f,
            },
            indent=2,
        )
    )
    return 0 if not f else 2


if __name__ == "__main__":
    raise SystemExit(main())
