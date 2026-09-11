#!/usr/bin/env python3
"""Exact same-history A21+[b_true] event cocycle for all mandatory bias families.

For a retained source-selector lineage this module constructs the actual joint
24-state differential cocycle

    X = [e_A21 ; b_true]

and the suffix-propagated source map for

    s_k = [w_bias,k ; m_tau,k].

The shipping A21 prediction differential comes from the exact finite predictor.
It is lifted with the family-specific physical phi_true interval.  The SAME w
column enters both accelerometer-bias error and physical true bias; tau mismatch
enters only the error.  Every accepted S=0/accelerometer/magnetometer event is
recomputed from the selector's captured same-event P/H/R cell, including actual
applied R_S and the deployed projection Clarke Jacobian, then lifted so the
same persistent b_true coordinate passes through the measurement/reset event.

The accumulated source map is not a packet-norm shortcut.  If A_k is the joint
event differential and B_k is nonzero only on prediction events, the recurrence
is exactly

    Phi_{k+1} = A_k Phi_k + B_k s_k,

and the returned B_word is formed by repeated variation of constants.  Thus all
later Joseph/reset/projection events are inside each prediction source suffix.

This constructor works on a concrete retained selector lineage.  Universality
comes from the separate source-reachable selector-family relation; this file
does not collapse that relation to one independent 601-sample box and does not
promote P4.  The H18 word has no active accelerometer-bias coordinate and is
handled by the existing exact nonlinear H18 lineage.  H18->A21 release is a
separate inter-word hybrid because shipping guarantees at least one complete
3 s H18 word before release.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval, matrix_identity, matrix_mul
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_p4_a21_bias1_24state_event_lift as LIFT
import ou3_p4_bias_family_joint_iss_supply as SUPPLY
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_complete_brmm_differential_prediction as PREDICTION
import ou3_p4_complete_brmm_same_history_prefix_selectors as SELECTORS
import ou3_p4_same_history_nonlinear_graph_lineage as GRAPH
import ou3_p4_source_reachable_selector_family as FAMILY
import ou3_p4_source_uniform_bias_prefix_lineage as BIASPREFIX

SCHEMA = 1
QUALIFICATION = "OU3_P4_ALL_BIAS_A21_24STATE_AUGMENTED_COCYCLE_V1"
P3_DELTA = 1.0e-18


def _shape(A):
    return len(A), len(A[0]) if A else 0


def _zero(r: int, c: int):
    z = Interval.point(0.0)
    return [[z for _ in range(c)] for _ in range(r)]


def _hstack(A, B):
    ar, ac = _shape(A)
    br, bc = _shape(B)
    if ar != br:
        raise ValueError("horizontal stack row mismatch")
    return [list(A[i]) + list(B[i]) for i in range(ar)]


def _propagate_source_map(A, Bprev, Bnew=None):
    propagated = matrix_mul(A, Bprev) if _shape(Bprev)[1] else _zero(_shape(A)[0], 0)
    return propagated if Bnew is None else _hstack(propagated, Bnew)


@dataclass(frozen=True)
class JointCocycleResult:
    family: str
    endpoint_source_cell_id: str
    A_word: tuple[tuple[Interval, ...], ...]
    B_word: tuple[tuple[Interval, ...], ...]
    event_kinds: tuple[str, ...]
    prediction_source_column_ranges: tuple[tuple[int, int], ...]
    projection_branches: tuple[str, ...]
    predictions: int
    measurements: int


def materialize_A21_joint_cocycle(
    selectors: Sequence[SELECTORS.PrefixSelector],
    endpoint_source_cell_id: str,
    *,
    family: str,
    initial_error_state: Sequence[Interval],
    domain_path: Path = GRAPH.DEFAULT_DOMAIN,
) -> JointCocycleResult:
    """Build the exact 24-state event cocycle and stacked prediction supply map."""
    if family not in BIASPREFIX.FAMILIES:
        raise ValueError("family must be BIAS0/BIAS1/BIAS2")
    if len(initial_error_state) != 21 or any(not isinstance(x, Interval) for x in initial_error_state):
        raise ValueError("A21 initial error state must be Interval[21]")

    attached = FAMILY.attach_endpoint_family(selectors, endpoint_source_cell_id)
    lineage = list(attached.selector_lineage)
    cert = attached.bias_lineages[family]
    constants = KERNEL._process_constants(domain_path)
    projection_limit = GRAPH._projection_limit(domain_path)

    Aword = matrix_identity(24)
    Bword = _zero(24, 0)
    state = list(initial_error_state)
    kinds: list[str] = []
    source_ranges: list[tuple[int, int]] = []
    projection_branches: list[str] = []
    predictions = 0
    measurements = 0

    for expected, selector in enumerate(lineage, start=1):
        if selector.prefix_length != expected:
            raise RuntimeError("selector lineage skipped a global prefix")
        cells = selector.A_event_cells
        if tuple(c.kind for c in cells) != tuple(selector.A_events_this_sample):
            raise RuntimeError("A21 event cells detached from shipping event slice")
        if not cells or cells[0].kind != "prediction":
            raise RuntimeError("A21 sample does not begin with prediction")
        beta_true = cert.at(selector.source_cell_id)

        for cell in cells:
            if cell.mode != "A":
                raise RuntimeError("non-A21 event in A21 joint cocycle")
            if cell.kind == "prediction":
                pred = PREDICTION.prediction_event(
                    "A", state, selector.sample_coordinates.omega_body_corrected,
                    constants.h, selector.active_schedule.tau,
                    tau_ba=constants.accel_bias_tau_s,
                )
                Aevent, Bsource = LIFT.prediction_lift(pred["J_state"], cert.phi_true)
                start = _shape(Bword)[1]
                Bword = _propagate_source_map(Aevent, Bword, Bsource)
                source_ranges.append((start, start + 6))
                Aword = matrix_mul(Aevent, Aword)
                state = list(pred["state_out"])
                predictions += 1
                kinds.append("prediction")
                continue

            if cell.kind == "aw_floor":
                Aevent = matrix_identity(24)
                Bword = _propagate_source_map(Aevent, Bword)
                Aword = matrix_mul(Aevent, Aword)
                kinds.append("aw_floor")
                continue

            if cell.kind not in ("S_zero", "accelerometer", "magnetometer"):
                raise RuntimeError(f"unsupported A21 event {cell.kind}")
            if cell.H is None or cell.R is None or not cell.P_before:
                raise RuntimeError(f"{cell.kind} lost same-event P/H/R")
            kwargs = GRAPH._measurement_geometry(selector, cell)
            provenance = EVENTS.ACTUAL_RS_PROVENANCE if cell.kind == "S_zero" else None
            if cell.kind == "S_zero" and cell.actual_rs_from_committed_schedule is not True:
                raise RuntimeError("S=0 event lost actual applied R_S provenance")
            source_event = EVENTS.source_joseph_event(
                "A", state, cell.P_before, cell.R, cell.kind,
                R_provenance=provenance,
                bias_true=beta_true,
                bias_projection_limit=projection_limit,
                **kwargs,
            )
            if source_event["H"] != cell.H:
                raise RuntimeError(f"{cell.kind} nonlinear H detached from captured shipping H")
            Aevent = LIFT.measurement_lift(source_event)
            Bword = _propagate_source_map(Aevent, Bword)
            Aword = matrix_mul(Aevent, Aword)
            state = list(source_event["state_out"])
            measurements += 1
            kinds.append(cell.kind)
            projection_branches.append(str(source_event["bias_projection_branch"]))

    if _shape(Aword) != (24, 24):
        raise RuntimeError("joint A21 cocycle dimension drifted")
    if _shape(Bword) != (24, 6 * predictions):
        raise RuntimeError("joint source map lost one 6-column block per prediction")
    if not LIFT.source_map_same_w(
        LIFT.prediction_lift(matrix_identity(21), cert.phi_true)[1]
    ):
        raise RuntimeError("family prediction lift lost shared w column")
    return JointCocycleResult(
        family=family,
        endpoint_source_cell_id=endpoint_source_cell_id,
        A_word=tuple(tuple(x for x in row) for row in Aword),
        B_word=tuple(tuple(x for x in row) for row in Bword),
        event_kinds=tuple(kinds),
        prediction_source_column_ranges=tuple(source_ranges),
        projection_branches=tuple(projection_branches),
        predictions=predictions,
        measurements=measurements,
    )


def build() -> dict:
    selector_family = FAMILY.build()
    ff = FAMILY.validate(selector_family)
    supply = SUPPLY.build()
    sf = SUPPLY.validate(supply)
    lift = LIFT.build()
    lf = LIFT.validate(lift)
    graph = GRAPH.build()
    gf = GRAPH.validate(graph)
    bad = {"selector_family": ff, "supply": sf, "lift": lf, "graph": gf}
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError("24-state cocycle prerequisites failed: " + repr(bad))

    family_rows = {}
    for family in BIASPREFIX.FAMILIES:
        s = supply["family_supply"][family]
        family_rows[family] = {
            "phi_true_interval": s["phi_true_interval"],
            "joint_supply_norm_upper_per_prediction_mps2": s[
                "joint_supply_norm_upper_per_prediction_mps2"
            ],
            "shared_w_prediction_lift": bool(
                s["event_lift"]["same_w_column_shared_by_error_and_truth"]
            ),
            "physical_factor_carried_in_truth_block": bool(
                s["event_lift"]["physical_factor_carried_in_truth_block"]
            ),
        }
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "source_reachable_selector_family_relation_consumed": selector_family[
            "source_reachable_COMPLETE_BRMM_selector_family_relation_closed"
        ],
        "joint_state": "[e_A21;b_true]",
        "joint_dimension": 24,
        "prediction_supply": "[w_bias;m_tau] per prediction",
        "family_lifts": family_rows,
        "all_three_bias_families_use_same_24state_event_architecture": True,
        "one_persistent_true_bias_state_across_word": True,
        "same_w_enters_error_and_true_bias": True,
        "tau_mismatch_enters_error_only": True,
        "same_event_projection_Clarke_Jacobian_lifted": True,
        "actual_RS_S_events_retained": True,
        "suffix_propagated_prediction_supply_map_materializer_available": True,
        "packetwise_source_norm_substitution_used": False,
        "independent_per_event_true_bias_state_used": False,
        "independent_K_or_RS_schedule_used": False,
        "finite_source_enumeration_used": False,
        "production_all_selector_lineages_numerically_executed_here": False,
        "source_uniform_joint_augmented_LDLT_closed_here": False,
        "first_exit_retention_closed_here": False,
        "P3_delta": P3_DELTA,
        "P4_MOTION_PASS": False,
        "P4_PASS": False,
        "P5_MAY_START": False,
        "next_obligation": (
            "for the universal selector-family relation, attach the exact chord/reset/S-zero centered graph sectors "
            "to this 24-state cocycle and its suffix-propagated family supply map, then close endpoint/every-prefix LDLT"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_BRMM_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    if d.get("joint_dimension") != 24:
        f.append("joint dimension changed")
    if set(d.get("family_lifts", {})) != set(BIASPREFIX.FAMILIES):
        f.append("family lift table incomplete")
    for k in (
        "source_reachable_selector_family_relation_consumed",
        "all_three_bias_families_use_same_24state_event_architecture",
        "one_persistent_true_bias_state_across_word",
        "same_w_enters_error_and_true_bias",
        "tau_mismatch_enters_error_only",
        "same_event_projection_Clarke_Jacobian_lifted",
        "actual_RS_S_events_retained",
        "suffix_propagated_prediction_supply_map_materializer_available",
    ):
        if d.get(k) is not True:
            f.append(k + " not true")
    for family, row in d.get("family_lifts", {}).items():
        if row.get("shared_w_prediction_lift") is not True:
            f.append(f"{family} lost shared-w lift")
        if row.get("physical_factor_carried_in_truth_block") is not True:
            f.append(f"{family} lost physical factor")
    for k in (
        "packetwise_source_norm_substitution_used",
        "independent_per_event_true_bias_state_used",
        "independent_K_or_RS_schedule_used",
        "finite_source_enumeration_used",
        "production_all_selector_lineages_numerically_executed_here",
        "source_uniform_joint_augmented_LDLT_closed_here",
        "first_exit_retention_closed_here",
        "P4_MOTION_PASS", "P4_PASS", "P5_MAY_START",
    ):
        if d.get(k) is not False:
            f.append(k + " not false")
    if d.get("P3_delta") != P3_DELTA:
        f.append("P3 delta changed")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build()
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "all_bias_joint_cocycle": d["all_three_bias_families_use_same_24state_event_architecture"],
        "source_map_materializer": d["suffix_propagated_prediction_supply_map_materializer_available"],
        "P4": d["P4_PASS"],
        "failures": f,
    }, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
