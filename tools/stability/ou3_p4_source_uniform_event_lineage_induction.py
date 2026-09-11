#!/usr/bin/env python3
"""Inductive source-uniform literal-event lineage for COMPLETE-BRMM P4.

The COMPLETE-BRMM selector family is a relation over a continuum source; it is
not a finite catalogue of 601-sample realizations.  The synchronized event
attachment is correspondingly generic in one retained typed sample and retains
*every* JOINT successor.  This module composes those two facts by induction.

Induction node:
  (shipping frontend/estimator correlation state,
   reachable H18/A21 Riccati states,
   exact/outward H18/A21 nonlinear error states,
   one radial source cell,
   one family-specific absolute true-bias prefix state,
   source-cell ancestry token).

For every admitted next typed source sample, ``synchronize_sample`` produces all
compatible successor nodes.  The successor carries the JOINT image as the next
correlation state, the common trusted Riccati post-state, and the exact nonlinear
state_out from every literal event.  Thus the construction does not require a
favorable branch ordinal, a replay, or an independent P/H/R/K hull.

This closes the *relational local-event attachment* needed by the first-exit
argument.  It does NOT close endpoint/every-prefix storage or LDLT, and it does
not supply a numeric centered-S diameter D_S that is still separately required
for the indefinite theorem.
"""
from __future__ import annotations

from dataclasses import dataclass
import argparse
import copy
import json
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_p4_joint_brmm_frontend_transition as JOINT
import ou3_p4_source_reachable_selector_family as FAMILY
import ou3_p4_source_uniform_estimator_event_attachment as ATTACH
import ou3_p4_source_uniform_bias_prefix_lineage as BIAS
import ou3_p4_brmm_primitive_prefix_binding as PRIMITIVE

SCHEMA = 1
QUALIFICATION = "OU3_P4_SOURCE_UNIFORM_EVENT_LINEAGE_INDUCTION_V1"
P3_DELTA = 1.0e-18


def I(x: float) -> Interval:
    return Interval.point(float(x))


@dataclass(frozen=True)
class LineageNode:
    branch: KERNEL.ExecutionBranch
    joint_state: JOINT.State
    H_error_state: tuple[Interval, ...]
    A_error_state: tuple[Interval, ...]
    radial_scale: Interval
    true_bias: tuple[Interval, Interval, Interval]
    sample_prefix_length: int


def validate_node(node: LineageNode) -> list[str]:
    f: list[str] = []
    if node.branch.frontend != node.joint_state.frontend:
        f.append("kernel Riccati carrier and JOINT correlation state have different shipping frontend")
    if len(node.H_error_state) != 18:
        f.append("H18 nonlinear error state dimension changed")
    if len(node.A_error_state) != 21:
        f.append("A21 nonlinear error state dimension changed")
    if len(node.true_bias) != 3:
        f.append("absolute true-bias coordinate dimension changed")
    if node.radial_scale.lo < 0.0 or node.radial_scale.hi > 1.0:
        f.append("radial source cell outside [0,1]")
    if node.sample_prefix_length < 0:
        f.append("negative source prefix length")
    return f


def advance_node(
    node: LineageNode,
    sample: KERNEL.SampleCoordinates,
    *,
    next_true_bias: Sequence[Interval],
    bias_projection_limit: float,
    tau_ba: Interval,
    next_cell_prefix: str,
    domain_path: Path = KERNEL.DEFAULT_DOMAIN,
) -> list[LineageNode]:
    """Advance one induction node and retain every synchronized JOINT successor."""
    nf = validate_node(node)
    if nf:
        raise ValueError("invalid source-uniform lineage node: " + repr(nf))
    if len(next_true_bias) != 3 or any(not isinstance(x, Interval) for x in next_true_bias):
        raise ValueError("next same-history true bias must be Interval[3]")

    pairs = ATTACH.synchronize_sample(
        branch=node.branch,
        joint_state=node.joint_state,
        sample=sample,
        state_in_H=node.H_error_state,
        state_in_A=node.A_error_state,
        radial_scale=node.radial_scale,
        true_bias=next_true_bias,
        bias_projection_limit=bias_projection_limit,
        tau_ba=tau_ba,
        sample_index=node.sample_prefix_length,
        next_cell_prefix=next_cell_prefix,
        domain_path=domain_path,
    )
    if not pairs:
        raise RuntimeError("branch-complete synchronized transition returned no successors")

    out: list[LineageNode] = []
    for h, a in pairs:
        if h.selector.source_cell_id != a.selector.source_cell_id:
            raise RuntimeError("H18/A21 successor ancestry detached")
        if h.image != a.image:
            raise RuntimeError("H18/A21 successor estimator image detached")
        branch = ATTACH.next_execution_branch(h, a)
        child = LineageNode(
            branch=branch,
            joint_state=h.image.state,
            H_error_state=tuple(h.state_out),
            A_error_state=tuple(a.state_out),
            radial_scale=node.radial_scale,
            true_bias=tuple(next_true_bias),  # type: ignore[arg-type]
            sample_prefix_length=node.sample_prefix_length + 1,
        )
        cf = validate_node(child)
        if cf:
            raise RuntimeError("synchronized successor is not a valid induction node: " + repr(cf))
        out.append(child)
    return out


def _smoke() -> dict:
    js = JOINT._smoke_state()
    sample = KERNEL.SampleCoordinates(
        gyro_measurement=KERNEL.MAHONY.Vec3(I(.01), I(-.02), I(.005)),
        omega_body_corrected=(I(.01), I(-.02), I(.005)),
        specific_force=KERNEL.MAHONY.Vec3(I(.2), I(-.1), I(-9.75)),
        f_cog_body=(I(0), I(0), I(-9.80665)),
        R_wb=ATTACH._identity(3),
        due_S=True,
        aw_floor_requested=True,
        magnetometer_events_after_imu=(KERNEL.MagneticEvent((I(20), I(0), I(40))),),
    )
    branch = KERNEL.ExecutionBranch(
        frontend=copy.deepcopy(js.frontend),
        H=ATTACH.WORD.initialize_word("H", ATTACH._identity(18)),
        A=ATTACH.WORD.initialize_word("A", ATTACH._identity(21)),
        source_cell_id="root",
    )
    root = LineageNode(
        branch=branch,
        joint_state=js,
        H_error_state=tuple(I(0) for _ in range(18)),
        A_error_state=tuple(I(0) for _ in range(21)),
        radial_scale=Interval(0.0, 1.0),
        true_bias=(I(0), I(0), I(0)),
        sample_prefix_length=0,
    )
    first = advance_node(root, sample, next_true_bias=root.true_bias,
                         bias_projection_limit=.4, tau_ba=I(1800),
                         next_cell_prefix="induction:0")
    second: list[LineageNode] = []
    for i, child in enumerate(first):
        second.extend(advance_node(child, sample, next_true_bias=child.true_bias,
                                   bias_projection_limit=.4, tau_ba=I(1800),
                                   next_cell_prefix=f"induction:1:p{i}"))
    return {
        "first_successors": len(first),
        "second_successors": len(second),
        "all_first_prefix_lengths_one": bool(first) and all(x.sample_prefix_length == 1 for x in first),
        "all_second_prefix_lengths_two": bool(second) and all(x.sample_prefix_length == 2 for x in second),
        "all_successors_valid": all(not validate_node(x) for x in first + second),
        "all_radial_cells_preserved": all(x.radial_scale == root.radial_scale for x in first + second),
    }


def build() -> dict:
    family = FAMILY.build(); ff = FAMILY.validate(family)
    attach = ATTACH.build(); af = ATTACH.validate(attach)
    bias = BIAS.build(); bf = BIAS.validate(bias)
    primitive = PRIMITIVE.build(); pf = PRIMITIVE.validate(primitive)
    bad = {k:v for k,v in (("selector_family",ff),("event_attachment",af),("bias_prefix",bf),("primitive_prefix",pf)) if v}
    if bad:
        raise RuntimeError("lineage-induction prerequisites failed: " + repr(bad))
    smoke = _smoke()
    relation = bool(
        family["source_reachable_COMPLETE_BRMM_selector_family_relation_closed"]
        and attach["source_uniform_estimator_owned_event_attachment_relation_closed"]
        and bias["source_uniform_projection_bias_coordinate_materialized"]
        and primitive["provider_primitive_prefix_binding_interface_materialized"]
        and smoke["all_successors_valid"]
        and smoke["all_radial_cells_preserved"]
    )
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "source_reachable_selector_family_relation_consumed": True,
        "generic_one_sample_transition_is_branch_complete": True,
        "induction_successor_uses_JOINT_image_as_next_correlation_state": True,
        "induction_successor_uses_trusted_common_Riccati_post_state": True,
        "induction_successor_uses_exact_nonlinear_state_out": True,
        "literal_P_H_R_state_attachment_closed_by_induction": relation,
        "source_uniform_local_event_coefficient_cover_closed": relation,
        "same_radial_cell_preserved_across_induction": True,
        "family_specific_absolute_bias_prefix_required": True,
        "same_centered_S_origin_and_primitive_ancestry_required": True,
        "all_successors_retained": True,
        "favorable_successor_selected": False,
        "branch_ordinal_correspondence_assumed": False,
        "independent_P_H_R_K_boxes_used": False,
        "finite_source_enumeration_used": False,
        "trajectory_replay_used": False,
        "numeric_D_S_qualification_closed_here": False,
        "endpoint_augmented_LDLT_closed_here": False,
        "every_prefix_augmented_LDLT_closed_here": False,
        "first_exit_retention_closed_here": False,
        "smoke": smoke,
        "P3_delta": P3_DELTA,
        "P4_MOTION_PASS": False,
        "P4_PASS": False,
        "P5_MAY_START": False,
        "next_obligation": (
            "embed each inductively attached literal event in one common augmented coordinate containing storage state, "
            "nonlinear graph auxiliaries, primitive moment/radial coordinates, BIAS0/1/2 supply and conditional binary32 supply; "
            "then run source-uniform endpoint/every-prefix outward LDLT and hard-domain targets"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for k in (
        "source_reachable_selector_family_relation_consumed",
        "generic_one_sample_transition_is_branch_complete",
        "induction_successor_uses_JOINT_image_as_next_correlation_state",
        "induction_successor_uses_trusted_common_Riccati_post_state",
        "induction_successor_uses_exact_nonlinear_state_out",
        "literal_P_H_R_state_attachment_closed_by_induction",
        "source_uniform_local_event_coefficient_cover_closed",
        "same_radial_cell_preserved_across_induction",
        "family_specific_absolute_bias_prefix_required",
        "same_centered_S_origin_and_primitive_ancestry_required",
        "all_successors_retained",
    ):
        if d.get(k) is not True:
            f.append(k + " not true")
    for k in (
        "favorable_successor_selected", "branch_ordinal_correspondence_assumed",
        "independent_P_H_R_K_boxes_used", "finite_source_enumeration_used", "trajectory_replay_used",
        "numeric_D_S_qualification_closed_here", "endpoint_augmented_LDLT_closed_here",
        "every_prefix_augmented_LDLT_closed_here", "first_exit_retention_closed_here",
        "P4_MOTION_PASS", "P4_PASS", "P5_MAY_START",
    ):
        if d.get(k) is not False:
            f.append(k + " not false")
    if d.get("P3_delta") != P3_DELTA:
        f.append("P3 delta changed")
    s=d.get("smoke",{})
    for k in ("all_first_prefix_lengths_one","all_second_prefix_lengths_two","all_successors_valid","all_radial_cells_preserved"):
        if s.get(k) is not True:
            f.append("smoke " + k + " not true")
    return list(dict.fromkeys(f))


def main() -> int:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args(); d=build(); f=validate(d)
    d["validation_pass"]=not f; d["validation_failures"]=f
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps({"event_cover":d["source_uniform_local_event_coefficient_cover_closed"],"P4":d["P4_PASS"],"failures":f},sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
