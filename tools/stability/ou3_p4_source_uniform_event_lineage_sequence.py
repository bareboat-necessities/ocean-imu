#!/usr/bin/env python3
"""Literal cross-sample ancestry for synchronized source-uniform event lineages.

A trusted sample selector has a sample-level parent/child id, while the terminal
P4 every-prefix certificate requires one literal event chain.  Therefore the
first event of sample k+1 must point to the last literal event of sample k; its
*estimator* predecessor remains the prior sample selector id.  These are two
different ancestries and must not be conflated.

This module stitches already synchronized AttachedSampleLineage objects without
changing any state, P/H/R, estimator coefficient, bias, radial or geometry
coordinate.  Only the literal event predecessor label of the first event in each
noninitial sample is replaced by the preceding literal event token.
"""
from __future__ import annotations

import argparse
import copy
from dataclasses import dataclass, replace
import json
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_brmm_full_normal_live_word as WORD
import ou3_p4_joint_brmm_frontend_transition as JOINT
import ou3_p4_complete_brmm_source_cover_contract as COVER
import ou3_p4_typed_sample_source_cell_binding as BIND
import ou3_p4_source_uniform_estimator_event_attachment as ATTACH
import ou3_p4_source_reachable_selector_family as FAMILY

SCHEMA = 1
QUALIFICATION = "OU3_P4_SOURCE_UNIFORM_LITERAL_EVENT_LINEAGE_SEQUENCE_V1"
P3_DELTA = 1.0e-18


def I(x: float) -> Interval:
    return Interval.point(float(x))


@dataclass(frozen=True)
class StitchedEventLineage:
    mode: str
    samples: tuple[ATTACH.AttachedSampleLineage, ...]
    cells: tuple[COVER.SourceCoverCell, ...]


def stitch_samples(samples: Sequence[ATTACH.AttachedSampleLineage]) -> StitchedEventLineage:
    if not samples:
        raise ValueError("nonempty synchronized sample lineage required")
    mode = samples[0].mode
    if mode not in ("H", "A"):
        raise ValueError("lineage mode must be H/A")
    flat: list[COVER.SourceCoverCell] = []
    previous_selector: str | None = None
    previous_event: str | None = None
    previous_prefix = 0

    for si, sample in enumerate(samples):
        if sample.mode != mode:
            raise ValueError("mixed H/A sample lineage")
        selector = sample.selector
        if selector.prefix_length != previous_prefix + 1:
            raise ValueError("sample-prefix ancestry skipped or duplicated a prefix")
        if previous_selector is not None and selector.parent_source_cell_id != previous_selector:
            raise ValueError("sample selector parent is not previous selector child")
        if sample.image.source_token != selector.source_cell_id:
            raise ValueError("sample estimator image detached from selector child")
        if sample.image.predecessor_token != selector.parent_source_cell_id:
            raise ValueError("sample estimator predecessor detached from selector parent")
        if not sample.cells:
            raise ValueError("sample lineage contains no literal events")

        cells = list(sample.cells)
        # The per-sample constructor uses the sample-level parent as its initial
        # literal predecessor.  Replace only that label after the first sample.
        if previous_event is not None:
            cells[0] = replace(cells[0], predecessor_token=previous_event)
        for ei, cell in enumerate(cells):
            failures = COVER.validate_cell(cell, require_estimator_provenance=True)
            if failures:
                raise RuntimeError(f"sample {si} event {ei} invalid: {failures}")
            if cell.estimator_source_token != selector.source_cell_id:
                raise RuntimeError("literal event lost estimator ownership")
            expected_pred = previous_event if ei == 0 and previous_event is not None else (
                selector.parent_source_cell_id if ei == 0 else cells[ei-1].source_token
            )
            if cell.predecessor_token != expected_pred:
                raise RuntimeError("literal event predecessor chain is not contiguous")
        if BIND.validate_event_cells_against_selector(selector, cells, mode=mode):
            raise RuntimeError("cross-sample relabel changed typed selector binding")

        flat.extend(cells)
        previous_event = cells[-1].source_token
        previous_selector = selector.source_cell_id
        previous_prefix = selector.prefix_length

    # Global literal ancestry check.
    for i in range(1, len(flat)):
        if flat[i].predecessor_token != flat[i-1].source_token:
            raise RuntimeError(f"global literal ancestry broken at event {i}")
    return StitchedEventLineage(mode, tuple(samples), tuple(flat))


def _identity(n: int):
    return [[I(1.0 if i == j else 0.0) for j in range(n)] for i in range(n)]


def _sample() -> KERNEL.SampleCoordinates:
    return KERNEL.SampleCoordinates(
        gyro_measurement=KERNEL.MAHONY.Vec3(I(.01), I(-.02), I(.005)),
        omega_body_corrected=(I(.01), I(-.02), I(.005)),
        specific_force=KERNEL.MAHONY.Vec3(I(.2), I(-.1), I(-9.75)),
        f_cog_body=(I(0), I(0), I(-9.80665)),
        R_wb=_identity(3), due_S=True, aw_floor_requested=True,
        magnetometer_events_after_imu=(KERNEL.MagneticEvent((I(20), I(0), I(40))),),
    )


def _two_sample_smoke() -> dict:
    js0 = JOINT._smoke_state()
    branch0 = KERNEL.ExecutionBranch(
        frontend=copy.deepcopy(js0.frontend),
        H=WORD.initialize_word("H", _identity(18)),
        A=WORD.initialize_word("A", _identity(21)),
        source_cell_id="root",
    )
    args = dict(
        radial_scale=Interval(0.0, 1.0), true_bias=[I(0), I(0), I(0)],
        bias_projection_limit=0.4, tau_ba=I(1800.0), sample=_sample(),
    )
    p0 = ATTACH.synchronize_sample(
        branch=branch0, joint_state=js0,
        state_in_H=[I(0) for _ in range(18)], state_in_A=[I(0) for _ in range(21)],
        sample_index=0, next_cell_prefix="seq-k0", **args,
    )
    if not p0:
        raise RuntimeError("first synchronized sample emitted no pair")
    h0, a0 = p0[0]
    branch1 = KERNEL.ExecutionBranch(
        frontend=copy.deepcopy(h0.image.state.frontend),
        H=copy.deepcopy(h0.selector.H_after),
        A=copy.deepcopy(a0.selector.A_after),
        source_cell_id=h0.selector.source_cell_id,
    )
    p1 = ATTACH.synchronize_sample(
        branch=branch1, joint_state=h0.image.state,
        state_in_H=h0.state_out, state_in_A=a0.state_out,
        sample_index=1, next_cell_prefix="seq-k1", **args,
    )
    if not p1:
        raise RuntimeError("second synchronized sample emitted no pair")
    h1, a1 = p1[0]
    H = stitch_samples((h0, h1)); A = stitch_samples((a0, a1))
    return {
        "H_literal_events": len(H.cells),
        "A_literal_events": len(A.cells),
        "H_cross_sample_predecessor_is_prior_literal_event": H.samples[1].selector.parent_source_cell_id != H.cells[len(h0.cells)].predecessor_token and H.cells[len(h0.cells)].predecessor_token == H.cells[len(h0.cells)-1].source_token,
        "A_cross_sample_predecessor_is_prior_literal_event": A.samples[1].selector.parent_source_cell_id != A.cells[len(a0.cells)].predecessor_token and A.cells[len(a0.cells)].predecessor_token == A.cells[len(a0.cells)-1].source_token,
        "H_estimator_parent_remains_prior_selector": H.cells[len(h0.cells)].estimator_predecessor_token == h0.selector.source_cell_id,
        "A_estimator_parent_remains_prior_selector": A.cells[len(a0.cells)].estimator_predecessor_token == a0.selector.source_cell_id,
    }


def build() -> dict:
    family = FAMILY.build(); ff = FAMILY.validate(family)
    attach = ATTACH.build(); af = ATTACH.validate(attach)
    if ff or af:
        raise RuntimeError(f"lineage-sequence prerequisites failed family={ff} attachment={af}")
    smoke = _two_sample_smoke()
    closed = bool(
        family["source_reachable_COMPLETE_BRMM_selector_family_relation_closed"]
        and attach["source_uniform_estimator_owned_event_attachment_relation_closed"]
        and smoke["H_cross_sample_predecessor_is_prior_literal_event"]
        and smoke["A_cross_sample_predecessor_is_prior_literal_event"]
        and smoke["H_estimator_parent_remains_prior_selector"]
        and smoke["A_estimator_parent_remains_prior_selector"]
    )
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "source_uniform_estimator_event_attachment_consumed": True,
        "sample_selector_and_literal_event_ancestry_kept_distinct": True,
        "first_event_of_next_sample_points_to_previous_literal_event": closed,
        "estimator_predecessor_remains_previous_selector_child": closed,
        "state_P_H_R_coefficients_unchanged_by_stitch": True,
        "radial_bias_geometry_coordinates_unchanged_by_stitch": True,
        "global_literal_event_lineage_constructor_closed": closed,
        "finite_source_enumeration_used": False,
        "favorable_branch_selection_part_of_constructor": False,
        "smoke": smoke,
        "production_common_augmented_coordinate_attached_here": False,
        "production_every_prefix_LDLT_closed_here": False,
        "P3_delta": P3_DELTA,
        "P4_MOTION_PASS": False,
        "P4_PASS": False,
        "P5_MAY_START": False,
        "next_obligation": "attach centered-S/moment/radial, family-specific bias supplies, graph sectors and binary32 coordinates to each flattened literal prefix, then run the outward augmented LDLT",
    }


def validate(d: dict) -> list[str]:
    f=[]
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for k in (
        "source_uniform_estimator_event_attachment_consumed",
        "sample_selector_and_literal_event_ancestry_kept_distinct",
        "first_event_of_next_sample_points_to_previous_literal_event",
        "estimator_predecessor_remains_previous_selector_child",
        "state_P_H_R_coefficients_unchanged_by_stitch",
        "radial_bias_geometry_coordinates_unchanged_by_stitch",
        "global_literal_event_lineage_constructor_closed",
    ):
        if d.get(k) is not True: f.append(k + " not true")
    for k in (
        "finite_source_enumeration_used", "favorable_branch_selection_part_of_constructor",
        "production_common_augmented_coordinate_attached_here", "production_every_prefix_LDLT_closed_here",
        "P4_MOTION_PASS", "P4_PASS", "P5_MAY_START",
    ):
        if d.get(k) is not False: f.append(k + " not false")
    if d.get("P3_delta") != P3_DELTA: f.append("P3 delta changed")
    s=d.get("smoke",{})
    for k in ("H_cross_sample_predecessor_is_prior_literal_event","A_cross_sample_predecessor_is_prior_literal_event","H_estimator_parent_remains_prior_selector","A_estimator_parent_remains_prior_selector"):
        if s.get(k) is not True: f.append("smoke " + k + " not true")
    return list(dict.fromkeys(f))


def main() -> int:
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("--output",type=Path,required=True); a=p.parse_args()
    d=build(); f=validate(d); d["validation_pass"]=not f; d["validation_failures"]=f
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"lineage":d["global_literal_event_lineage_constructor_closed"],"smoke":d["smoke"],"P4":d["P4_PASS"],"failures":f},sort_keys=True)); return int(bool(f))


if __name__ == "__main__": raise SystemExit(main())
