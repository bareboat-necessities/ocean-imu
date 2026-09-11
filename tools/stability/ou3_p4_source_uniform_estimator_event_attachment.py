#!/usr/bin/env python3
"""Estimator-owned event attachment for source-uniform P4 lineages.

For one retained source branch/sample this module executes, from the SAME
predecessor state and SAME typed sample:

  * the trusted Riccati/kernel transition with passive P/H/R event capture; and
  * the correlated JOINT Mahony/stillness/WPE/tuner transition.

The trusted kernel is authoritative only for the current sample's Riccati event
slice.  Its auxiliary FRONT successor must NOT be used to choose the next proof
frontend state: FRONT does not carry the theorem's same-history stillness/
central-statistics refinement.  The JOINT image is the authoritative next
frontend state.  This is valid because the current Riccati schedule is committed
before the post-measurement tuner transition; every retained JOINT image must
therefore agree with the kernel on the current active tau/sigma and actual R_S.

All kernel frontend successors share the same post-Riccati H18/A21 state because
frontend branching occurs only after the current IMU Riccati events and before
asynchronous magnetic events, which are then applied identically.  This module
verifies that equality and uses the common post-Riccati state as the covariance
carrier for EACH retained JOINT successor.  No favorable kernel frontend child
is selected and no branch-ordinal correspondence is assumed.

For each JOINT successor, SourceCoverCell events are instantiated directly from
the trusted event-local P_before/H/R snapshots while the physical error state is
advanced through the exact deployed nonlinear prediction/Joseph/reset maps.
Only theorem ancestry labels are assigned; numerical/statistical coordinates are
unchanged.  P4 remains fail-closed until the resulting literal lineages are
embedded in the common augmented PrefixInput and endpoint/every-prefix LDLT plus
hard-domain retention close.
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
import ou3_brmm_frontend_state_step as FRONT
import ou3_brmm_full_normal_live_word as WORD
import ou3_p4_joint_brmm_frontend_transition as JOINT
import ou3_p4_complete_brmm_same_history_prefix_selectors as SELECTORS
import ou3_p4_complete_brmm_source_cover_contract as COVER
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_complete_brmm_differential_prediction as PRED
import ou3_p4_typed_sample_source_cell_binding as BIND
import ou3_p4_typed_sample_nonlinear_chain as CHAIN
import ou3_p4_source_reachable_selector_family as FAMILY

SCHEMA = 2
QUALIFICATION = "OU3_P4_SOURCE_UNIFORM_ESTIMATOR_OWNED_EVENT_ATTACHMENT_V2"
P3_DELTA = 1.0e-18


def I(x: float) -> Interval:
    return Interval.point(float(x))


def _same_interval(a: Interval, b: Interval) -> bool:
    return isinstance(a, Interval) and isinstance(b, Interval) and a.lo == b.lo and a.hi == b.hi


def _same_vec(a, b) -> bool:
    return len(a) == len(b) and all(_same_interval(x, y) for x, y in zip(a, b))


def _same_mat(A, B) -> bool:
    return len(A) == len(B) and all(_same_vec(x, y) for x, y in zip(A, B))


def _same_word(a: WORD.LiteralWordState, b: WORD.LiteralWordState) -> bool:
    return (
        a.mode == b.mode
        and _same_mat(a.riccati.P, b.riccati.P)
        and a.event_log == b.event_log
        and a.imu_samples == b.imu_samples
        and a.accel_updates == b.accel_updates
        and a.S_updates == b.S_updates
        and a.mag_updates == b.mag_updates
        and a.aw_floor_applications == b.aw_floor_applications
    )


@dataclass(frozen=True)
class AttachedSampleLineage:
    selector: SELECTORS.PrefixSelector
    image: JOINT.Image
    mode: str
    cells: tuple[COVER.SourceCoverCell, ...]
    state_out: tuple[Interval, ...]
    J_word: object


def _retoken_image(image: JOINT.Image, *, child: str, parent: str) -> JOINT.Image:
    """Change theorem ancestry labels only; retain every numerical coordinate."""
    state = replace(image.state, source_token=child)
    return replace(image, source_token=child, predecessor_token=parent, state=state)


def _common_kernel_post(children: Sequence[KERNEL.ExecutionBranch]):
    if not children:
        raise RuntimeError("trusted kernel emitted no child")
    h = children[0].H
    a = children[0].A
    if any(not _same_word(c.H, h) or not _same_word(c.A, a) for c in children[1:]):
        raise RuntimeError("kernel frontend branch unexpectedly changed current Riccati post-state")
    return h, a


def _selector_for_image(
    branch: KERNEL.ExecutionBranch,
    sample: KERNEL.SampleCoordinates,
    image: JOINT.Image,
    meta: dict,
    H_after: WORD.LiteralWordState,
    A_after: WORD.LiteralWordState,
    *,
    sample_index: int,
    successor_ordinal: int,
    child_id: str,
) -> SELECTORS.PrefixSelector:
    H_cells = tuple(copy.deepcopy(meta["H_event_cells"]))
    A_cells = tuple(copy.deepcopy(meta["A_event_cells"]))
    return SELECTORS.PrefixSelector(
        sample_index=sample_index,
        prefix_length=sample_index + 1,
        parent_branch_ordinal=0,
        successor_ordinal=successor_ordinal,
        parent_source_cell_id=branch.source_cell_id,
        source_cell_id=child_id,
        sample_coordinates=copy.deepcopy(sample),
        active_schedule=copy.deepcopy(image.active_schedule_for_current_riccati),
        actual_rs_std_xyz=copy.deepcopy(tuple(image.actual_rs_std_xyz_for_current_riccati)),
        H_before=copy.deepcopy(branch.H),
        H_after=copy.deepcopy(H_after),
        A_before=copy.deepcopy(branch.A),
        A_after=copy.deepcopy(A_after),
        H_events_this_sample=tuple(c.kind for c in H_cells),
        A_events_this_sample=tuple(c.kind for c in A_cells),
        H_event_cells=H_cells,
        A_event_cells=A_cells,
        H_floor_case=meta.get("H_floor_case"),
        A_floor_case=meta.get("A_floor_case"),
    )


def next_execution_branch(h: AttachedSampleLineage, a: AttachedSampleLineage) -> KERNEL.ExecutionBranch:
    """Carry the common Riccati state with the authoritative JOINT frontend."""
    if h.selector.source_cell_id != a.selector.source_cell_id or h.image != a.image:
        raise ValueError("H/A attachments are not the same retained JOINT successor")
    return KERNEL.ExecutionBranch(
        frontend=copy.deepcopy(h.image.state.frontend),
        H=copy.deepcopy(h.selector.H_after),
        A=copy.deepcopy(a.selector.A_after),
        source_cell_id=h.selector.source_cell_id,
    )


def _event_cells(
    selector: SELECTORS.PrefixSelector,
    image: JOINT.Image,
    *,
    mode: str,
    state_in: Sequence[Interval],
    dt_s: Interval,
    radial_scale: Interval,
    true_bias: Sequence[Interval] | None,
    bias_projection_limit: float | None,
    tau_ba: Interval | None,
) -> tuple[tuple[COVER.SourceCoverCell, ...], tuple[Interval, ...], object]:
    trusted = selector.H_event_cells if mode == "H" else selector.A_event_cells
    n = 18 if mode == "H" else 21
    if len(state_in) != n:
        raise ValueError(f"{mode} entry state must have dimension {n}")
    if mode == "A":
        if true_bias is None or len(true_bias) != 3:
            raise ValueError("A21 attachment requires one same-history true-bias coordinate")
        if not (bias_projection_limit is not None and float(bias_projection_limit) > 0.0):
            raise ValueError("A21 attachment requires shipping projection radius")
        if tau_ba is None or tau_ba.lo <= 0.0:
            raise ValueError("A21 attachment requires positive tau_ba")

    state = list(state_in)
    cells: list[COVER.SourceCoverCell] = []
    predecessor = selector.parent_source_cell_id
    pseudo = I((selector.sample_index + 1) * dt_s.hi)

    for ordinal, event in enumerate(trusted):
        token = f"{selector.source_cell_id}:e{ordinal}"
        kwargs = dict(
            mode=mode,
            sample_index=selector.sample_index,
            event_ordinal=ordinal,
            kind=event.kind,
            state=list(state),
            P=copy.deepcopy(event.P_before),
            dt_s=dt_s,
            pseudo_elapsed_s=pseudo,
            radial_scale=radial_scale,
            event_source_token=token,
            event_predecessor_token=predecessor,
            wave_primitive=selector.sample_coordinates.wave_primitive,
        )
        if mode == "A":
            kwargs["true_bias"] = list(true_bias or ())
            kwargs["bias_projection_limit"] = float(bias_projection_limit)
        if event.kind == "accelerometer":
            kwargs["R"] = copy.deepcopy(event.R)
            kwargs["f_hat"] = list(selector.sample_coordinates.f_cog_body)
            kwargs["R_hat"] = copy.deepcopy(selector.sample_coordinates.R_wb)
        elif event.kind == "magnetometer":
            kwargs["R"] = copy.deepcopy(event.R)
            mi = event.magnetic_event_index
            if mi is None or not (0 <= mi < len(selector.sample_coordinates.magnetometer_events_after_imu)):
                raise ValueError("trusted magnetic event lost typed-sample index")
            kwargs["m_body"] = list(selector.sample_coordinates.magnetometer_events_after_imu[mi].m_body)

        cell = COVER.source_cell_from_joint_image(image, **kwargs)
        cells.append(cell)
        if event.kind == "prediction":
            p = PRED.prediction_event(
                mode, state, selector.sample_coordinates.omega_body_corrected,
                dt_s, cell.tau_applied_s, tau_ba=tau_ba,
            )
            state = list(p["state_out"])
        elif event.kind == "aw_floor":
            pass
        elif event.kind in ("S_zero", "accelerometer", "magnetometer"):
            j = EVENTS.source_joseph_event(**COVER.joseph_event_kwargs(cell))
            state = list(j["state_out"])
        else:
            raise ValueError(f"unsupported sample-local event {event.kind!r}")
        predecessor = token

    failures = BIND.validate_event_cells_against_selector(selector, cells, mode=mode)
    if failures:
        raise RuntimeError("synchronized event attachment failed selector binding: " + repr(failures))
    chain = CHAIN.materialize_sample_chain(selector, cells, mode=mode, tau_ba=tau_ba)
    if tuple(chain["event_kinds"]) != tuple(c.kind for c in cells):
        raise RuntimeError("nonlinear chain changed trusted event order")
    return tuple(cells), tuple(chain["state_out"]), chain["J_word"]


def synchronize_sample(
    *,
    branch: KERNEL.ExecutionBranch,
    joint_state: JOINT.State,
    sample: KERNEL.SampleCoordinates,
    state_in_H: Sequence[Interval],
    state_in_A: Sequence[Interval],
    radial_scale: Interval,
    true_bias: Sequence[Interval],
    bias_projection_limit: float,
    tau_ba: Interval,
    sample_index: int,
    next_cell_prefix: str,
    domain_path: Path = KERNEL.DEFAULT_DOMAIN,
) -> list[tuple[AttachedSampleLineage, AttachedSampleLineage]]:
    """Attach every retained JOINT successor to the common current Riccati slice."""
    constants = KERNEL._process_constants(Path(domain_path))
    if branch.frontend != joint_state.frontend:
        raise ValueError("kernel and JOINT predecessors do not share the same shipping frontend state")

    children, meta, images = KERNEL.advance_branch_with_joint_frontend(
        branch, sample, joint_state=joint_state, constants=constants,
        next_cell_prefix=next_cell_prefix + ":kernel",
        joint_child_prefix=next_cell_prefix + ":joint",
        capture_riccati_event_cells=True,
    )
    if not images:
        raise RuntimeError("JOINT transition emitted an empty successor family")
    if meta.get("same_active_schedule_verified") is not True or meta.get("same_actual_RS_verified") is not True:
        raise RuntimeError("trusted kernel did not verify current estimator schedule/R_S")
    H_after, A_after = _common_kernel_post(children)

    committed = KERNEL.TUNER.commit_if_pending(branch.frontend.tuner, constants.tuner)
    expected_active = committed.active
    expected_rs = tuple(KERNEL.TUNER.active_rs_std_xyz(expected_active, constants.tuner))
    for image in images:
        if image.active_schedule_for_current_riccati != expected_active:
            raise RuntimeError("JOINT image current active schedule detached from trusted Riccati schedule")
        if tuple(image.actual_rs_std_xyz_for_current_riccati) != expected_rs:
            raise RuntimeError("JOINT image actual R_S detached from trusted Riccati schedule")

    out=[]; dt=I(constants.tuner.dt)
    for ordinal, raw_image in enumerate(images):
        child_id=f"{next_cell_prefix}:joint{ordinal}"
        image=_retoken_image(raw_image, child=child_id, parent=branch.source_cell_id)
        selector=_selector_for_image(
            branch, sample, image, meta, H_after, A_after,
            sample_index=sample_index, successor_ordinal=ordinal, child_id=child_id,
        )
        Hcells,Hout,HJ=_event_cells(
            selector,image,mode="H",state_in=state_in_H,dt_s=dt,radial_scale=radial_scale,
            true_bias=None,bias_projection_limit=None,tau_ba=None)
        Acells,Aout,AJ=_event_cells(
            selector,image,mode="A",state_in=state_in_A,dt_s=dt,radial_scale=radial_scale,
            true_bias=true_bias,bias_projection_limit=bias_projection_limit,tau_ba=tau_ba)
        out.append((
            AttachedSampleLineage(selector,image,"H",Hcells,Hout,HJ),
            AttachedSampleLineage(selector,image,"A",Acells,Aout,AJ),
        ))
    return out


def _identity(n: int):
    return [[I(1.0 if i == j else 0.0) for j in range(n)] for i in range(n)]


def _smoke() -> dict:
    js=JOINT._smoke_state()
    sample=KERNEL.SampleCoordinates(
        gyro_measurement=KERNEL.MAHONY.Vec3(I(.01),I(-.02),I(.005)),
        omega_body_corrected=(I(.01),I(-.02),I(.005)),
        specific_force=KERNEL.MAHONY.Vec3(I(.2),I(-.1),I(-9.75)),
        f_cog_body=(I(0),I(0),I(-9.80665)),R_wb=_identity(3),due_S=True,
        aw_floor_requested=True,magnetometer_events_after_imu=(KERNEL.MagneticEvent((I(20),I(0),I(40))),))
    branch=KERNEL.ExecutionBranch(frontend=copy.deepcopy(js.frontend),H=WORD.initialize_word("H",_identity(18)),A=WORD.initialize_word("A",_identity(21)),source_cell_id="root")
    pairs=synchronize_sample(
        branch=branch,joint_state=js,sample=sample,state_in_H=[I(0) for _ in range(18)],state_in_A=[I(0) for _ in range(21)],
        radial_scale=Interval(0.0,1.0),true_bias=[I(0),I(0),I(0)],bias_projection_limit=.4,tau_ba=I(1800),sample_index=0,next_cell_prefix="sync-smoke")
    next_branches=[next_execution_branch(h,a) for h,a in pairs]
    return {
        "joint_successors_retained":len(pairs),
        "all_H_cells_bound":bool(pairs) and all(not BIND.validate_event_cells_against_selector(h.selector,h.cells,mode="H") for h,_ in pairs),
        "all_A_cells_bound":bool(pairs) and all(not BIND.validate_event_cells_against_selector(a.selector,a.cells,mode="A") for _,a in pairs),
        "all_event_orders_equal":bool(pairs) and all(tuple(c.kind for c in h.cells)==h.selector.H_events_this_sample and tuple(c.kind for c in a.cells)==a.selector.A_events_this_sample for h,a in pairs),
        "authoritative_next_frontend_is_joint":bool(next_branches) and all(b.frontend==h.image.state.frontend for b,(h,_) in zip(next_branches,pairs)),
    }


def build() -> dict:
    family=FAMILY.build(); ff=FAMILY.validate(family)
    bind=BIND.build(); bf=BIND.validate(bind)
    chain=CHAIN.build(); cf=CHAIN.validate(chain)
    if ff or bf or cf: raise RuntimeError(f"attachment prerequisites failed family={ff} binding={bf} chain={cf}")
    smoke=_smoke()
    closed=bool(family["source_reachable_COMPLETE_BRMM_selector_family_relation_closed"] and smoke["joint_successors_retained"]>0 and smoke["all_H_cells_bound"] and smoke["all_A_cells_bound"] and smoke["all_event_orders_equal"] and smoke["authoritative_next_frontend_is_joint"])
    return {
        "schema":SCHEMA,"qualification":QUALIFICATION,"canonical_source":"COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "source_reachable_selector_family_relation_consumed":True,
        "same_predecessor_same_sample_kernel_and_joint_execution":True,
        "obsolete_measured_only_frontend_not_used_in_joint_attachment":True,
        "physical_wave_payload_retained_in_every_literal_event_cell":True,
        "physical_generator_to_joint24_forcing_attachment_closed_here":False,
        "kernel_authoritative_for_current_Riccati_slice_only":True,
        "JOINT_image_authoritative_for_next_frontend_state":True,
        "current_schedule_committed_before_post_measurement_adaptation":True,
        "every_joint_successor_retained":True,
        "all_kernel_frontend_children_share_current_post_Riccati_state":True,
        "kernel_frontend_successor_used_to_select_next_theorem_state":False,
        "kernel_joint_child_frontend_equality_required":False,
        "branch_ordinal_correspondence_assumed":False,
        "theorem_ancestry_retoken_only_numeric_coordinates_unchanged":True,
        "trusted_event_local_P_H_R_imported_without_reconstruction":True,
        "exact_nonlinear_state_succession_materialized_between_literal_events":True,
        "same_finite_map_generates_state_and_Jacobian":True,
        "actual_applied_RS_owned_by_same_joint_image":True,
        "H18_and_A21_event_attachments_constructible_for_each_retained_joint_successor":closed,
        "source_uniform_estimator_owned_event_attachment_relation_closed":closed,
        "favorable_successor_selected":False,"independent_P_H_R_K_reconstruction_used":False,"numeric_COMPLETE_BRMM_enumeration_used":False,
        "smoke":smoke,"production_augmented_PrefixInput_assembled_here":False,"production_endpoint_LDLT_closed_here":False,"production_every_prefix_LDLT_closed_here":False,"first_exit_retention_closed_here":False,
        "P3_delta":P3_DELTA,"P4_MOTION_PASS":False,"P4_PASS":False,"P5_MAY_START":False,
        "next_obligation":"stitch these estimator-owned event cells into the common augmented PrefixInput coordinate with centered-S/moment/radial, BIAS0/1/2 supply and conditional binary32 maps; then run endpoint/every-prefix outward LDLT and hard-ball first-exit targets",
    }


def validate(d:dict)->list[str]:
    f=[]
    if d.get("schema")!=SCHEMA or d.get("qualification")!=QUALIFICATION:f.append("schema/qualification mismatch")
    for k in ("source_reachable_selector_family_relation_consumed","same_predecessor_same_sample_kernel_and_joint_execution","kernel_authoritative_for_current_Riccati_slice_only","JOINT_image_authoritative_for_next_frontend_state","current_schedule_committed_before_post_measurement_adaptation","every_joint_successor_retained","all_kernel_frontend_children_share_current_post_Riccati_state","theorem_ancestry_retoken_only_numeric_coordinates_unchanged","trusted_event_local_P_H_R_imported_without_reconstruction","exact_nonlinear_state_succession_materialized_between_literal_events","same_finite_map_generates_state_and_Jacobian","actual_applied_RS_owned_by_same_joint_image","H18_and_A21_event_attachments_constructible_for_each_retained_joint_successor","source_uniform_estimator_owned_event_attachment_relation_closed"):
        if d.get(k) is not True:f.append(k+" not true")
    for k in ("kernel_frontend_successor_used_to_select_next_theorem_state","kernel_joint_child_frontend_equality_required","branch_ordinal_correspondence_assumed","favorable_successor_selected","independent_P_H_R_K_reconstruction_used","numeric_COMPLETE_BRMM_enumeration_used","production_augmented_PrefixInput_assembled_here","production_endpoint_LDLT_closed_here","production_every_prefix_LDLT_closed_here","first_exit_retention_closed_here","P4_MOTION_PASS","P4_PASS","P5_MAY_START"):
        if d.get(k) is not False:f.append(k+" not false")
    if d.get("P3_delta")!=P3_DELTA:f.append("P3 delta changed")
    s=d.get("smoke",{})
    if int(s.get("joint_successors_retained",0))<=0:f.append("smoke produced no joint successors")
    for k in ("all_H_cells_bound","all_A_cells_bound","all_event_orders_equal","authoritative_next_frontend_is_joint"):
        if s.get(k) is not True:f.append("smoke "+k+" not true")
    return list(dict.fromkeys(f))


def main()->int:
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--output",type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d["validation_pass"]=not f;d["validation_failures"]=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n");print(json.dumps({"attachment":d["source_uniform_estimator_owned_event_attachment_relation_closed"],"successors":d["smoke"]["joint_successors_retained"],"P4":d["P4_PASS"],"failures":f},sort_keys=True));return int(bool(f))


if __name__=="__main__":raise SystemExit(main())
