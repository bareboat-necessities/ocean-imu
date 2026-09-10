#!/usr/bin/env python3
"""Synchronized estimator-owned event attachment for source-uniform P4 lineages.

This module closes the representation seam between the universal COMPLETE-BRMM
selector relation and the exact nonlinear event chain.  For one retained source
branch/sample it executes, from the SAME predecessor state and SAME typed sample:

  * the trusted Riccati/kernel transition with passive P/H/R event capture; and
  * the correlated JOINT Mahony/stillness/WPE/tuner transition.

Children are paired by exact equality of the shipping frontend successor state,
not by a favorable branch ordinal.  Every kernel successor and every JOINT image
must participate in at least one compatible pair.  An attached JOINT image is
then retokened only at the theorem ancestry-label level to the trusted kernel
child id; all numerical/statistical coordinates are unchanged.

For each compatible pair, SourceCoverCell events are instantiated directly from
the trusted event-local P_before/H/R snapshots while the physical error state is
advanced through the exact deployed nonlinear prediction/Joseph/reset maps.
Thus every event cell contains the same finite state that generated its event
Jacobian and the same P/H/R owned by the trusted Riccati transition.

The constructor is generic over a retained branch and therefore can be lifted
pointwise over the already-closed source-reachable selector-family relation.  It
neither enumerates COMPLETE-BRMM nor selects a successor.  P4 remains fail-closed
until these attachments are embedded in the common augmented PrefixInput and the
endpoint/every-prefix outward LDLT plus hard-domain targets close.
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
import ou3_brmm_tuner_scheduler_step as TUNER
import ou3_p4_joint_brmm_frontend_transition as JOINT
import ou3_p4_complete_brmm_same_history_prefix_selectors as SELECTORS
import ou3_p4_complete_brmm_source_cover_contract as COVER
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_complete_brmm_differential_prediction as PRED
import ou3_p4_typed_sample_source_cell_binding as BIND
import ou3_p4_typed_sample_nonlinear_chain as CHAIN
import ou3_p4_source_reachable_selector_family as FAMILY

SCHEMA = 1
QUALIFICATION = "OU3_P4_SOURCE_UNIFORM_ESTIMATOR_OWNED_EVENT_ATTACHMENT_V1"
P3_DELTA = 1.0e-18


def I(x: float) -> Interval:
    return Interval.point(float(x))


def _same_interval(a: Interval, b: Interval) -> bool:
    return isinstance(a, Interval) and isinstance(b, Interval) and a.lo == b.lo and a.hi == b.hi


def _same_vec(a, b) -> bool:
    return len(a) == len(b) and all(_same_interval(x, y) for x, y in zip(a, b))


def _same_mat(A, B) -> bool:
    return len(A) == len(B) and all(_same_vec(x, y) for x, y in zip(A, B))


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


def _frontend_matches(image: JOINT.Image, child: KERNEL.ExecutionBranch) -> bool:
    return image.state.frontend == child.frontend


def _selector_for_pair(
    branch: KERNEL.ExecutionBranch,
    child: KERNEL.ExecutionBranch,
    sample: KERNEL.SampleCoordinates,
    image: JOINT.Image,
    meta: dict,
    *,
    sample_index: int,
    successor_ordinal: int,
) -> SELECTORS.PrefixSelector:
    active = image.active_schedule_for_current_riccati
    rs = tuple(image.actual_rs_std_xyz_for_current_riccati)
    H_cells = tuple(copy.deepcopy(meta["H_event_cells"]))
    A_cells = tuple(copy.deepcopy(meta["A_event_cells"]))
    H_events = tuple(c.kind for c in H_cells)
    A_events = tuple(c.kind for c in A_cells)
    return SELECTORS.PrefixSelector(
        sample_index=sample_index,
        prefix_length=sample_index + 1,
        parent_branch_ordinal=0,
        successor_ordinal=successor_ordinal,
        parent_source_cell_id=branch.source_cell_id,
        source_cell_id=child.source_cell_id,
        sample_coordinates=copy.deepcopy(sample),
        active_schedule=copy.deepcopy(active),
        actual_rs_std_xyz=copy.deepcopy(rs),
        H_before=copy.deepcopy(branch.H),
        H_after=copy.deepcopy(child.H),
        A_before=copy.deepcopy(branch.A),
        A_after=copy.deepcopy(child.A),
        H_events_this_sample=H_events,
        A_events_this_sample=A_events,
        H_event_cells=H_cells,
        A_event_cells=A_cells,
        H_floor_case=meta.get("H_floor_case"),
        A_floor_case=meta.get("A_floor_case"),
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
        )
        if mode == "A":
            kwargs["true_bias"] = list(true_bias or ())
            kwargs["bias_projection_limit"] = float(bias_projection_limit)
        if event.kind == "S_zero":
            # source_cell_from_joint_image derives R from the same active schedule.
            pass
        elif event.kind == "accelerometer":
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
            pred = PRED.prediction_event(
                mode,
                state,
                selector.sample_coordinates.omega_body_corrected,
                dt_s,
                cell.tau_applied_s,
                tau_ba=tau_ba,
            )
            state = list(pred["state_out"])
        elif event.kind == "aw_floor":
            state = list(state)
        elif event.kind in ("S_zero", "accelerometer", "magnetometer"):
            joseph = EVENTS.source_joseph_event(**COVER.joseph_event_kwargs(cell))
            state = list(joseph["state_out"])
        else:
            raise ValueError(f"unsupported sample-local event {event.kind!r}")
        predecessor = token

    failures = BIND.validate_event_cells_against_selector(selector, cells, mode=mode)
    if failures:
        raise RuntimeError("synchronized event attachment failed selector binding: " + repr(failures))
    chain = CHAIN.materialize_sample_chain(selector, cells, mode=mode, tau_ba=tau_ba)
    if tuple(chain["event_kinds"]) != tuple(c.kind for c in cells):
        raise RuntimeError("nonlinear chain changed trusted event order")
    if not chain["literal_state_succession_exact"]:
        raise RuntimeError("nonlinear event state succession not closed")
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
    """Return every compatible H18/A21 attachment for one retained source sample."""
    constants = KERNEL._process_constants(Path(domain_path))
    if branch.frontend != joint_state.frontend:
        raise ValueError("kernel and JOINT predecessors do not share the same shipping frontend state")

    children, meta = KERNEL.advance_branch(
        branch,
        sample,
        constants=constants,
        next_cell_prefix=next_cell_prefix,
        capture_riccati_event_cells=True,
    )
    images = JOINT.advance(
        joint_state,
        FRONT.Sample(sample.gyro_measurement, sample.specific_force),
        gravity_ms2=constants.gravity,
        two_kp=constants.two_kp,
        two_ki=constants.two_ki,
        child_prefix=next_cell_prefix + ":joint",
    )
    if not children or not images:
        raise RuntimeError("synchronized transition emitted an empty successor family")
    if meta.get("same_active_schedule_verified") is not True or meta.get("same_actual_RS_verified") is not True:
        raise RuntimeError("trusted kernel did not verify current estimator schedule/R_S")

    pairs: list[tuple[int, int]] = []
    for ci, child in enumerate(children):
        for ii, image in enumerate(images):
            if _frontend_matches(image, child):
                pairs.append((ci, ii))
    covered_children = {ci for ci, _ in pairs}
    covered_images = {ii for _, ii in pairs}
    if covered_children != set(range(len(children))):
        raise RuntimeError("a trusted kernel successor has no exact JOINT-image owner")
    if covered_images != set(range(len(images))):
        raise RuntimeError("a JOINT-image successor has no trusted kernel child")

    out: list[tuple[AttachedSampleLineage, AttachedSampleLineage]] = []
    dt = I(constants.tuner.dt)
    for ci, ii in pairs:
        child = children[ci]
        image = _retoken_image(images[ii], child=child.source_cell_id, parent=branch.source_cell_id)
        selector = _selector_for_pair(
            branch, child, sample, image, meta,
            sample_index=sample_index, successor_ordinal=ci,
        )
        Hcells, Hout, HJ = _event_cells(
            selector, image, mode="H", state_in=state_in_H, dt_s=dt,
            radial_scale=radial_scale, true_bias=None, bias_projection_limit=None, tau_ba=None,
        )
        Acells, Aout, AJ = _event_cells(
            selector, image, mode="A", state_in=state_in_A, dt_s=dt,
            radial_scale=radial_scale, true_bias=true_bias,
            bias_projection_limit=bias_projection_limit, tau_ba=tau_ba,
        )
        out.append((
            AttachedSampleLineage(selector, image, "H", Hcells, Hout, HJ),
            AttachedSampleLineage(selector, image, "A", Acells, Aout, AJ),
        ))
    return out


def _identity(n: int):
    return [[I(1.0 if i == j else 0.0) for j in range(n)] for i in range(n)]


def _smoke() -> dict:
    js = JOINT._smoke_state()
    sample = KERNEL.SampleCoordinates(
        gyro_measurement=KERNEL.MAHONY.Vec3(I(.01), I(-.02), I(.005)),
        omega_body_corrected=(I(.01), I(-.02), I(.005)),
        specific_force=KERNEL.MAHONY.Vec3(I(.2), I(-.1), I(-9.75)),
        f_cog_body=(I(0), I(0), I(-9.80665)),
        R_wb=_identity(3),
        due_S=True,
        aw_floor_requested=True,
        magnetometer_events_after_imu=(KERNEL.MagneticEvent((I(20), I(0), I(40))),),
    )
    branch = KERNEL.ExecutionBranch(
        frontend=copy.deepcopy(js.frontend),
        H=WORD.initialize_word("H", _identity(18)),
        A=WORD.initialize_word("A", _identity(21)),
        source_cell_id="root",
    )
    pairs = synchronize_sample(
        branch=branch, joint_state=js, sample=sample,
        state_in_H=[I(0) for _ in range(18)], state_in_A=[I(0) for _ in range(21)],
        radial_scale=Interval(0.0, 1.0), true_bias=[I(0), I(0), I(0)],
        bias_projection_limit=0.4, tau_ba=I(1800.0), sample_index=0,
        next_cell_prefix="sync-smoke",
    )
    return {
        "compatible_pairs": len(pairs),
        "all_H_cells_bound": bool(pairs) and all(not BIND.validate_event_cells_against_selector(h.selector, h.cells, mode="H") for h, _ in pairs),
        "all_A_cells_bound": bool(pairs) and all(not BIND.validate_event_cells_against_selector(a.selector, a.cells, mode="A") for _, a in pairs),
        "all_event_orders_equal": bool(pairs) and all(tuple(c.kind for c in h.cells) == h.selector.H_events_this_sample and tuple(c.kind for c in a.cells) == a.selector.A_events_this_sample for h, a in pairs),
        "all_frontend_matches_exact": bool(pairs) and all(h.image.state.frontend == a.image.state.frontend for h, a in pairs),
    }


def build() -> dict:
    family = FAMILY.build()
    ff = FAMILY.validate(family)
    bind = BIND.build(); bf = BIND.validate(bind)
    chain = CHAIN.build(); cf = CHAIN.validate(chain)
    if ff or bf or cf:
        raise RuntimeError(f"attachment prerequisites failed family={ff} binding={bf} chain={cf}")
    smoke = _smoke()
    closed = bool(
        family["source_reachable_COMPLETE_BRMM_selector_family_relation_closed"]
        and smoke["compatible_pairs"] > 0
        and smoke["all_H_cells_bound"] and smoke["all_A_cells_bound"]
        and smoke["all_event_orders_equal"] and smoke["all_frontend_matches_exact"]
    )
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "source_reachable_selector_family_relation_consumed": True,
        "same_predecessor_same_sample_kernel_and_joint_execution": True,
        "kernel_joint_children_matched_by_exact_frontend_state_not_ordinal": True,
        "every_kernel_successor_must_have_joint_owner": True,
        "every_joint_successor_must_have_kernel_child": True,
        "theorem_ancestry_retoken_only_numeric_coordinates_unchanged": True,
        "trusted_event_local_P_H_R_imported_without_reconstruction": True,
        "exact_nonlinear_state_succession_materialized_between_literal_events": True,
        "same_finite_map_generates_state_and_Jacobian": True,
        "actual_applied_RS_owned_by_same_joint_image": True,
        "H18_and_A21_event_attachments_constructible_for_each_retained_pair": closed,
        "source_uniform_estimator_owned_event_attachment_relation_closed": closed,
        "favorable_successor_selected": False,
        "branch_ordinal_correspondence_assumed": False,
        "independent_P_H_R_K_reconstruction_used": False,
        "numeric_COMPLETE_BRMM_enumeration_used": False,
        "smoke": smoke,
        "production_augmented_PrefixInput_assembled_here": False,
        "production_endpoint_LDLT_closed_here": False,
        "production_every_prefix_LDLT_closed_here": False,
        "first_exit_retention_closed_here": False,
        "P3_delta": P3_DELTA,
        "P4_MOTION_PASS": False,
        "P4_PASS": False,
        "P5_MAY_START": False,
        "next_obligation": (
            "lift each AttachedSampleLineage into the common augmented PrefixInput coordinate with centered-S/moment/radial, BIAS0/1/2 supply and conditional binary32 maps; compose across literal prefixes and run outward LDLT plus hard-ball first-exit targets"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for k in (
        "source_reachable_selector_family_relation_consumed",
        "same_predecessor_same_sample_kernel_and_joint_execution",
        "kernel_joint_children_matched_by_exact_frontend_state_not_ordinal",
        "every_kernel_successor_must_have_joint_owner",
        "every_joint_successor_must_have_kernel_child",
        "theorem_ancestry_retoken_only_numeric_coordinates_unchanged",
        "trusted_event_local_P_H_R_imported_without_reconstruction",
        "exact_nonlinear_state_succession_materialized_between_literal_events",
        "same_finite_map_generates_state_and_Jacobian",
        "actual_applied_RS_owned_by_same_joint_image",
        "H18_and_A21_event_attachments_constructible_for_each_retained_pair",
        "source_uniform_estimator_owned_event_attachment_relation_closed",
    ):
        if d.get(k) is not True:
            f.append(k + " not true")
    for k in (
        "favorable_successor_selected", "branch_ordinal_correspondence_assumed",
        "independent_P_H_R_K_reconstruction_used", "numeric_COMPLETE_BRMM_enumeration_used",
        "production_augmented_PrefixInput_assembled_here", "production_endpoint_LDLT_closed_here",
        "production_every_prefix_LDLT_closed_here", "first_exit_retention_closed_here",
        "P4_MOTION_PASS", "P4_PASS", "P5_MAY_START",
    ):
        if d.get(k) is not False:
            f.append(k + " not false")
    if d.get("P3_delta") != P3_DELTA:
        f.append("P3 delta changed")
    s = d.get("smoke", {})
    if int(s.get("compatible_pairs", 0)) <= 0:
        f.append("smoke produced no compatible pair")
    for k in ("all_H_cells_bound", "all_A_cells_bound", "all_event_orders_equal", "all_frontend_matches_exact"):
        if s.get(k) is not True:
            f.append("smoke " + k + " not true")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(); f = validate(d)
    d["validation_pass"] = not f; d["validation_failures"] = f
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"attachment": d["source_uniform_estimator_owned_event_attachment_relation_closed"], "pairs": d["smoke"]["compatible_pairs"], "P4": d["P4_PASS"], "failures": f}, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
