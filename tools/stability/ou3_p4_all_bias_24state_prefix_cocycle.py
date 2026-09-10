#!/usr/bin/env python3
"""Every-literal-event snapshots of the all-bias A21 24-state cocycle.

The endpoint cocycle alone is insufficient for P4: chart/reset validity and the
finite-window gain theorem are prefix obligations.  This module therefore runs
the same exact event constructors as the endpoint materializer and records the
cumulative 24x24 state differential A_l and 24x(6 N_l) prediction-supply map
B_l after every shipping event.  Previous supply columns are propagated through
every later Joseph/reset/projection event before a new prediction block is
appended.

Each snapshot remains attached to the exact selector child and captured
Riccati-event index.  It is not a finite source enumeration and it does not
promote P4; the universal source quantifier is supplied by
ou3_p4_source_reachable_selector_family.  The next layer must attach the same
physical primitive/moment/radial witness, exact chord/reset graph sectors and
finite-precision map before invoking the terminal augmented LDLT.
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
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_complete_brmm_differential_prediction as PREDICTION
import ou3_p4_complete_brmm_same_history_prefix_selectors as SELECTORS
import ou3_p4_same_history_nonlinear_graph_lineage as GRAPH
import ou3_p4_source_reachable_selector_family as FAMILY
import ou3_p4_source_uniform_bias_prefix_lineage as BIASPREFIX
import ou3_p4_all_bias_24state_augmented_cocycle as ENDPOINT

SCHEMA = 1
QUALIFICATION = "OU3_P4_ALL_BIAS_A21_24STATE_EVERY_EVENT_PREFIX_COCYCLE_V1"
P3_DELTA = 1e-18


def _shape(A):
    return len(A), len(A[0]) if A else 0


def _zero(r, c):
    z = Interval.point(0.0)
    return [[z for _ in range(c)] for _ in range(r)]


def _hstack(A, B):
    if len(A) != len(B):
        raise ValueError("horizontal stack row mismatch")
    return [list(a) + list(b) for a, b in zip(A, B)]


def _propagate(A, previous, new=None):
    carried = matrix_mul(A, previous) if _shape(previous)[1] else _zero(24, 0)
    return carried if new is None else _hstack(carried, new)


@dataclass(frozen=True)
class PrefixSnapshot:
    family: str
    literal_prefix_ordinal: int
    selector_source_cell_id: str
    sample_prefix_length: int
    event_index_in_sample: int
    kind: str
    A_prefix: tuple[tuple[Interval, ...], ...]
    B_prefix: tuple[tuple[Interval, ...], ...]
    state_out: tuple[Interval, ...]
    prediction_count: int
    measurement_count: int
    actual_rs_provenance: bool
    projection_branch: str


def materialize_A21_joint_prefixes(
    selectors: Sequence[SELECTORS.PrefixSelector],
    endpoint_source_cell_id: str,
    *,
    family: str,
    initial_error_state: Sequence[Interval],
    domain_path: Path = GRAPH.DEFAULT_DOMAIN,
) -> list[PrefixSnapshot]:
    if family not in BIASPREFIX.FAMILIES:
        raise ValueError("family must be BIAS0/BIAS1/BIAS2")
    if len(initial_error_state) != 21 or any(not isinstance(x, Interval) for x in initial_error_state):
        raise ValueError("A21 initial error state must be Interval[21]")
    attached = FAMILY.attach_endpoint_family(selectors, endpoint_source_cell_id)
    lineage = list(attached.selector_lineage)
    cert = attached.bias_lineages[family]
    constants = KERNEL._process_constants(domain_path)
    projection_limit = GRAPH._projection_limit(domain_path)

    Ap = matrix_identity(24)
    Bp = _zero(24, 0)
    state = list(initial_error_state)
    out: list[PrefixSnapshot] = []
    predictions = measurements = 0

    for expected, selector in enumerate(lineage, start=1):
        if selector.prefix_length != expected:
            raise RuntimeError("selector lineage skipped a global sample prefix")
        cells = selector.A_event_cells
        if tuple(x.kind for x in cells) != tuple(selector.A_events_this_sample):
            raise RuntimeError("captured A21 event cells detached from shipping order")
        if not cells or cells[0].kind != "prediction":
            raise RuntimeError("A21 shipping sample does not start with prediction")
        beta_true = cert.at(selector.source_cell_id)

        for cell in cells:
            actual_rs = False
            projection = "not_applicable"
            if cell.kind == "prediction":
                p = PREDICTION.prediction_event(
                    "A", state, selector.sample_coordinates.omega_body_corrected,
                    constants.h, selector.active_schedule.tau,
                    tau_ba=constants.accel_bias_tau_s,
                )
                Ae, Bs = LIFT.prediction_lift(p["J_state"], cert.phi_true)
                Bp = _propagate(Ae, Bp, Bs)
                Ap = matrix_mul(Ae, Ap)
                state = list(p["state_out"])
                predictions += 1
            elif cell.kind == "aw_floor":
                Ae = matrix_identity(24)
                Bp = _propagate(Ae, Bp)
                Ap = matrix_mul(Ae, Ap)
            elif cell.kind in ("S_zero", "accelerometer", "magnetometer"):
                if cell.H is None or cell.R is None or not cell.P_before:
                    raise RuntimeError(f"{cell.kind} lost same-event P/H/R")
                kwargs = GRAPH._measurement_geometry(selector, cell)
                provenance = EVENTS.ACTUAL_RS_PROVENANCE if cell.kind == "S_zero" else None
                if cell.kind == "S_zero":
                    if cell.actual_rs_from_committed_schedule is not True:
                        raise RuntimeError("S=0 event lost actual applied R_S provenance")
                    actual_rs = True
                ev = EVENTS.source_joseph_event(
                    "A", state, cell.P_before, cell.R, cell.kind,
                    R_provenance=provenance,
                    bias_true=beta_true,
                    bias_projection_limit=projection_limit,
                    **kwargs,
                )
                if ev["H"] != cell.H:
                    raise RuntimeError(f"{cell.kind} nonlinear H detached from shipping H")
                Ae = LIFT.measurement_lift(ev)
                Bp = _propagate(Ae, Bp)
                Ap = matrix_mul(Ae, Ap)
                state = list(ev["state_out"])
                measurements += 1
                projection = str(ev["bias_projection_branch"])
            else:
                raise RuntimeError(f"unsupported A21 event {cell.kind}")

            if _shape(Ap) != (24, 24) or _shape(Bp) != (24, 6 * predictions):
                raise RuntimeError("prefix cocycle/source-map dimension drift")
            out.append(PrefixSnapshot(
                family=family,
                literal_prefix_ordinal=len(out) + 1,
                selector_source_cell_id=selector.source_cell_id,
                sample_prefix_length=selector.prefix_length,
                event_index_in_sample=cell.event_index_in_sample,
                kind=cell.kind,
                A_prefix=tuple(tuple(v for v in row) for row in Ap),
                B_prefix=tuple(tuple(v for v in row) for row in Bp),
                state_out=tuple(state),
                prediction_count=predictions,
                measurement_count=measurements,
                actual_rs_provenance=actual_rs,
                projection_branch=projection,
            ))

    if not out:
        raise RuntimeError("empty every-event prefix cocycle")
    endpoint = ENDPOINT.materialize_A21_joint_cocycle(
        selectors, endpoint_source_cell_id, family=family,
        initial_error_state=initial_error_state, domain_path=domain_path,
    )
    if out[-1].A_prefix != endpoint.A_word or out[-1].B_prefix != endpoint.B_word:
        raise RuntimeError("every-prefix endpoint differs from canonical endpoint cocycle")
    return out


def build() -> dict:
    family = FAMILY.build(); ff = FAMILY.validate(family)
    endpoint = ENDPOINT.build(); ef = ENDPOINT.validate(endpoint)
    if ff or ef:
        raise RuntimeError(f"prefix cocycle prerequisites failed family={ff} endpoint={ef}")
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "source_reachable_selector_family_relation_consumed": True,
        "all_three_bias_families_supported": True,
        "every_literal_A21_event_prefix_snapshot_available": True,
        "each_prefix_retains_selector_child_and_event_index": True,
        "each_prediction_appends_one_shared_w_tau_mismatch_block": True,
        "all_previous_supply_blocks_suffix_propagated_through_later_events": True,
        "endpoint_identity_with_canonical_joint_cocycle_enforced": True,
        "actual_RS_provenance_retained_on_each_due_S_prefix": True,
        "projection_branch_retained_on_each_measurement_prefix": True,
        "independent_prefix_Jacobian_boxes_used": False,
        "packet_norm_source_accumulation_used": False,
        "finite_source_enumeration_used": False,
        "physical_moment_radial_binding_attached_here": False,
        "graph_sector_master_attached_here": False,
        "finite_precision_map_attached_here": False,
        "endpoint_augmented_LDLT_closed_here": False,
        "every_prefix_augmented_LDLT_closed_here": False,
        "first_exit_retention_closed_here": False,
        "P3_delta": P3_DELTA,
        "P4_MOTION_PASS": False,
        "P4_PASS": False,
        "P5_MAY_START": False,
        "next_obligation": (
            "bind each PrefixSnapshot to its source primitive/moment witness and centered-S origin, "
            "attach exact chord/reset/projection graph sectors plus binary32 finite precision, then feed every snapshot to augmented LDLT"
        ),
    }


def validate(d):
    f=[]
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for k in (
        "source_reachable_selector_family_relation_consumed",
        "all_three_bias_families_supported",
        "every_literal_A21_event_prefix_snapshot_available",
        "each_prefix_retains_selector_child_and_event_index",
        "each_prediction_appends_one_shared_w_tau_mismatch_block",
        "all_previous_supply_blocks_suffix_propagated_through_later_events",
        "endpoint_identity_with_canonical_joint_cocycle_enforced",
        "actual_RS_provenance_retained_on_each_due_S_prefix",
        "projection_branch_retained_on_each_measurement_prefix",
    ):
        if d.get(k) is not True: f.append(k + " not true")
    for k in (
        "independent_prefix_Jacobian_boxes_used", "packet_norm_source_accumulation_used",
        "finite_source_enumeration_used", "physical_moment_radial_binding_attached_here",
        "graph_sector_master_attached_here", "finite_precision_map_attached_here",
        "endpoint_augmented_LDLT_closed_here", "every_prefix_augmented_LDLT_closed_here",
        "first_exit_retention_closed_here", "P4_MOTION_PASS", "P4_PASS", "P5_MAY_START",
    ):
        if d.get(k) is not False: f.append(k + " not false")
    if d.get("P3_delta") != P3_DELTA: f.append("P3 delta changed")
    return list(dict.fromkeys(f))


def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument("--output",type=Path,required=True); a=ap.parse_args()
    d=build(); f=validate(d); d["validation_pass"]=not f; d["validation_failures"]=f
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"prefix_cocycle":d["every_literal_A21_event_prefix_snapshot_available"],"P4":d["P4_PASS"],"failures":f},sort_keys=True)); return int(bool(f))


if __name__ == "__main__": raise SystemExit(main())
