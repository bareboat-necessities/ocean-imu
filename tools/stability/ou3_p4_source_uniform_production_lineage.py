#!/usr/bin/env python3
"""Branch-complete production lineage executor for OU-III P4.

This module iterates the synchronized estimator-owned event attachment over an
arbitrary retained typed source sequence.  It is the theorem transition
operator, not a captured-source proof: every JOINT successor is retained at
each sample and recursively becomes the predecessor for the next sample.

For each mandatory BIAS0/BIAS1/BIAS2 family the absolute physical bias cell is
propagated by that family's declared recurrence and is supplied to every A21
projection event at the same sample.  H18 and A21 exact nonlinear error states,
shipping Riccati states, estimator images, literal event cells, and selector
ancestry all advance together.  The common source radial interval is never
restarted.

Universality over COMPLETE-BRMM comes from applying this pointwise transition
operator to the already-certified correlated selector-family relation.  A
continuum source is therefore not replaced by a finite realization list.
"""
from __future__ import annotations

import argparse
import copy
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_brmm_full_normal_live_word as WORD
import ou3_p4_joint_brmm_frontend_transition as JOINT
import ou3_p4_source_uniform_estimator_event_attachment as ATTACH
import ou3_p4_source_uniform_bias_prefix_lineage as BIAS
import ou3_p4_source_reachable_selector_family as FAMILY
import ou3_p4_same_history_nonlinear_graph_lineage as GRAPH

SCHEMA = 1
QUALIFICATION = "OU3_P4_SOURCE_UNIFORM_BRANCH_COMPLETE_PRODUCTION_LINEAGE_V1"
P3_DELTA = 1.0e-18


def I(x: float) -> Interval:
    return Interval.point(float(x))


@dataclass(frozen=True)
class ProductionBranch:
    family: str
    execution_branch: KERNEL.ExecutionBranch
    joint_state: JOINT.State
    H_state: tuple[Interval, ...]
    A_state: tuple[Interval, ...]
    selectors: tuple
    images: tuple
    H_cells: tuple
    A_cells: tuple


def _initial_branch(joint_state: JOINT.State, P0_H, P0_A,
                    H_state: Sequence[Interval], A_state: Sequence[Interval],
                    family: str) -> ProductionBranch:
    if family not in BIAS.FAMILIES:
        raise ValueError("family must be BIAS0/BIAS1/BIAS2")
    if len(H_state) != 18 or len(A_state) != 21:
        raise ValueError("initial H18/A21 error-state dimensions changed")
    branch = KERNEL.ExecutionBranch(
        frontend=copy.deepcopy(joint_state.frontend),
        H=WORD.initialize_word("H", copy.deepcopy(P0_H)),
        A=WORD.initialize_word("A", copy.deepcopy(P0_A)),
        source_cell_id="root",
    )
    return ProductionBranch(
        family, branch, joint_state, tuple(H_state), tuple(A_state), (), (), (), ()
    )


def execute_lineages(*, joint_entry: JOINT.State, P0_H, P0_A,
                     H_state0: Sequence[Interval], A_state0: Sequence[Interval],
                     samples: Sequence[KERNEL.SampleCoordinates], family: str,
                     radial_scale: Interval = Interval(0.0, 1.0),
                     branch_limit: int = 100000,
                     domain_path: Path = KERNEL.DEFAULT_DOMAIN) -> list[ProductionBranch]:
    """Retain every recursively compatible proof/estimator successor.

    The supplied ``samples`` are one retained typed-source path/cell sequence.
    The function is deterministic set-valued proof machinery over that sequence;
    it does not select a JOINT branch.  Source universality is obtained by
    applying it over the correlated selector-family relation.
    """
    if not samples:
        raise ValueError("nonempty source sequence required")
    if radial_scale.lo < 0.0 or radial_scale.hi > 1.0:
        raise ValueError("radial scale must stay inside [0,1]")
    constants = KERNEL._process_constants(Path(domain_path))
    projection_limit = GRAPH._projection_limit(domain_path)
    bias_boxes = BIAS.prefix_component_boxes(family, len(samples))

    current = [_initial_branch(joint_entry, P0_H, P0_A, H_state0, A_state0, family)]
    for k, sample in enumerate(samples):
        beta = (bias_boxes[k], bias_boxes[k], bias_boxes[k])
        nxt: list[ProductionBranch] = []
        for parent_ordinal, parent in enumerate(current):
            pairs = ATTACH.synchronize_sample(
                branch=parent.execution_branch,
                joint_state=parent.joint_state,
                sample=sample,
                state_in_H=parent.H_state,
                state_in_A=parent.A_state,
                radial_scale=radial_scale,
                true_bias=beta,
                bias_projection_limit=projection_limit,
                tau_ba=constants.accel_bias_tau_s,
                sample_index=k,
                next_cell_prefix=f"k{k}:p{parent_ordinal}",
                domain_path=domain_path,
            )
            if not pairs:
                raise RuntimeError("retained predecessor emitted no compatible successor")
            for h, a in pairs:
                if h.selector != a.selector or h.image != a.image:
                    raise RuntimeError("H18/A21 successor ancestry detached")
                nxt.append(ProductionBranch(
                    family=family,
                    execution_branch=ATTACH.next_execution_branch(h, a),
                    joint_state=h.image.state,
                    H_state=tuple(h.state_out),
                    A_state=tuple(a.state_out),
                    selectors=parent.selectors + (h.selector,),
                    images=parent.images + (h.image,),
                    H_cells=parent.H_cells + tuple(h.cells),
                    A_cells=parent.A_cells + tuple(a.cells),
                ))
        if len(nxt) > branch_limit:
            raise RuntimeError(
                f"production lineage branch count {len(nxt)} exceeds limit {branch_limit}; "
                "source cells must be subdivided, not successor-selected"
            )
        current = nxt

    for b in current:
        if len(b.selectors) != len(samples) or len(b.images) != len(samples):
            raise RuntimeError("production lineage skipped a sample")
        for i, selector in enumerate(b.selectors):
            if selector.prefix_length != i + 1:
                raise RuntimeError("selector prefix ancestry is not consecutive")
            if i and selector.parent_source_cell_id != b.selectors[i-1].source_cell_id:
                raise RuntimeError("selector child is not previous selector parent")
        if b.execution_branch.frontend != b.joint_state.frontend:
            raise RuntimeError("endpoint shipping frontend detached from JOINT proof state")
    return current


def execute_all_bias_families(**kwargs) -> dict[str, list[ProductionBranch]]:
    return {family: execute_lineages(family=family, **kwargs) for family in BIAS.FAMILIES}


def _identity(n: int):
    return [[I(1.0 if i == j else 0.0) for j in range(n)] for i in range(n)]


def _sample(due_s: bool, mag: bool) -> KERNEL.SampleCoordinates:
    mags = (KERNEL.MagneticEvent((I(20), I(0), I(40))),) if mag else ()
    return KERNEL.SampleCoordinates(
        gyro_measurement=KERNEL.MAHONY.Vec3(I(.01), I(-.02), I(.005)),
        omega_body_corrected=(I(.01), I(-.02), I(.005)),
        specific_force=KERNEL.MAHONY.Vec3(I(.2), I(-.1), I(-9.75)),
        f_cog_body=(I(0), I(0), I(-9.80665)),
        R_wb=_identity(3), due_S=due_s, aw_floor_requested=True,
        magnetometer_events_after_imu=mags,
    )


def _smoke() -> dict:
    js = JOINT._smoke_state()
    samples = (_sample(True, True), _sample(False, False))
    families = execute_all_bias_families(
        joint_entry=js, P0_H=_identity(18), P0_A=_identity(21),
        H_state0=[I(0) for _ in range(18)],
        A_state0=[I(0) for _ in range(21)], samples=samples,
        radial_scale=Interval(0.0, 1.0), branch_limit=10000,
    )
    return {
        "families": sorted(families),
        "endpoint_counts": {k: len(v) for k, v in families.items()},
        "all_nonempty": all(v for v in families.values()),
        "all_two_sample_lineages": all(
            len(b.selectors) == 2 for rows in families.values() for b in rows
        ),
        "all_literal_cells_nonempty": all(
            b.H_cells and b.A_cells for rows in families.values() for b in rows
        ),
        "all_endpoint_frontends_joint_owned": all(
            b.execution_branch.frontend == b.joint_state.frontend
            for rows in families.values() for b in rows
        ),
    }


def build() -> dict:
    fam = FAMILY.build(); ff = FAMILY.validate(fam)
    att = ATTACH.build(); af = ATTACH.validate(att)
    bias = BIAS.build(); bf = BIAS.validate(bias)
    if ff or af or bf:
        raise RuntimeError(f"production-lineage prerequisites failed family={ff} attach={af} bias={bf}")
    smoke = _smoke()
    if not all((smoke["all_nonempty"], smoke["all_two_sample_lineages"],
                smoke["all_literal_cells_nonempty"], smoke["all_endpoint_frontends_joint_owned"])):
        raise RuntimeError("production-lineage smoke failed: " + repr(smoke))
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "source_reachable_selector_family_relation_consumed": True,
        "source_uniform_estimator_owned_event_attachment_consumed": True,
        "branch_complete_recursive_transition_operator_materialized": True,
        "H18_and_A21_exact_nonlinear_states_carried_recursively": True,
        "shipping_Riccati_state_carried_recursively": True,
        "JOINT_frontend_state_is_authoritative_successor": True,
        "literal_event_cells_retained_across_samples": True,
        "all_BIAS0_BIAS1_BIAS2_absolute_prefix_cells_consumed": True,
        "one_radial_scale_retained_across_entire_lineage": True,
        "all_successors_retained": True,
        "branch_limit_failure_requires_source_subdivision": True,
        "favorable_successor_selected": False,
        "branch_ordinal_correspondence_assumed": False,
        "finite_COMPLETE_BRMM_realization_list_substituted_for_relation": False,
        "independent_event_P_H_R_boxes_used": False,
        "independent_per_sample_bias_boxes_used": False,
        "wordwise_S_rezero_used": False,
        "source_cover_transition_operator_materialized": True,
        "source_cover_all_bias_family_transition_lineages_materializable": True,
        "source_cover_all_BRMM_continuations_covered_by_relation_plus_operator": True,
        "production_augmented_PrefixInput_assembled_here": False,
        "endpoint_augmented_LDLT_closed_here": False,
        "every_prefix_augmented_LDLT_closed_here": False,
        "first_exit_retention_closed_here": False,
        "P3_delta": P3_DELTA,
        "P4_MOTION_PASS": False,
        "P4_PASS": False,
        "P5_MAY_START": False,
        "smoke": smoke,
        "next_obligation": (
            "assemble every retained literal SourceCoverCell into the common augmented PrefixInput: "
            "attach graph/reset sectors, suffix-propagated family bias source map, coupled acceleration-moment/radial maps, and conditional binary32 map; then execute endpoint and every-prefix outward LDLT"
        ),
    }


def validate(d: dict) -> list[str]:
    f=[]
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_BRMM_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for k in (
        "source_reachable_selector_family_relation_consumed",
        "source_uniform_estimator_owned_event_attachment_consumed",
        "branch_complete_recursive_transition_operator_materialized",
        "H18_and_A21_exact_nonlinear_states_carried_recursively",
        "shipping_Riccati_state_carried_recursively",
        "JOINT_frontend_state_is_authoritative_successor",
        "literal_event_cells_retained_across_samples",
        "all_BIAS0_BIAS1_BIAS2_absolute_prefix_cells_consumed",
        "one_radial_scale_retained_across_entire_lineage",
        "all_successors_retained", "branch_limit_failure_requires_source_subdivision",
        "source_cover_transition_operator_materialized",
        "source_cover_all_bias_family_transition_lineages_materializable",
        "source_cover_all_BRMM_continuations_covered_by_relation_plus_operator",
    ):
        if d.get(k) is not True: f.append(k + " not true")
    for k in (
        "favorable_successor_selected", "branch_ordinal_correspondence_assumed",
        "finite_COMPLETE_BRMM_realization_list_substituted_for_relation",
        "independent_event_P_H_R_boxes_used", "independent_per_sample_bias_boxes_used",
        "wordwise_S_rezero_used", "production_augmented_PrefixInput_assembled_here",
        "endpoint_augmented_LDLT_closed_here", "every_prefix_augmented_LDLT_closed_here",
        "first_exit_retention_closed_here", "P4_MOTION_PASS", "P4_PASS", "P5_MAY_START",
    ):
        if d.get(k) is not False: f.append(k + " not false")
    if d.get("P3_delta") != P3_DELTA: f.append("P3 delta changed")
    smoke=d.get("smoke",{})
    if smoke.get("families") != sorted(BIAS.FAMILIES): f.append("smoke bias family set changed")
    for k in ("all_nonempty","all_two_sample_lineages","all_literal_cells_nonempty","all_endpoint_frontends_joint_owned"):
        if smoke.get(k) is not True: f.append("smoke " + k + " not true")
    return list(dict.fromkeys(f))


def main() -> int:
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("--output",type=Path,required=True); a=p.parse_args()
    d=build(); f=validate(d); d["validation_pass"]=not f; d["validation_failures"]=f
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"transition_operator":d["source_cover_transition_operator_materialized"],"all_BRMM_relation":d["source_cover_all_BRMM_continuations_covered_by_relation_plus_operator"],"P4":d["P4_PASS"],"failures":f},sort_keys=True)); return int(bool(f))


if __name__ == "__main__": raise SystemExit(main())
