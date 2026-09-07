#!/usr/bin/env python3
"""Branch-correlated every-prefix selectors for complete-SEA3 OU-III P4.

P4 cannot be certified from an endpoint word object alone.  The finite-state
theorem requires the endpoint and every finite prefix to refer to ONE
same-history complete SEA3 execution.  In particular, the nonlinear graph
assembler must not reconstruct independent boxes for the tuner schedule,
actual applied R_S, Riccati state, source sample, or front-end branch.

This module is a non-promoting bridge from the retained trusted typed execution
kernel to the P4 joint-sector machinery.  It does not create a source family,
does not bypass the canonical hard-window provider gate, and does not change
shipping algebra.  Every transition is executed by
``ou3_sea3_complete_window_execution_kernel.advance_branch``.

For every retained front-end successor after sample k, one PrefixSelector keeps

* the exact parent and child source-cell identifiers, so ancestry is explicit;
* the provider sample coordinate object used on that transition;
* the active tau/sigma/pseudo-period schedule committed before that sample;
* the actual anisotropic per-axis R_S derived from that same schedule;
* deep snapshots of the H18/A21 literal-word state before and after the sample;
* the exact shipping event slice appended on that sample.

Those objects are sufficient to attach the next rigorous objects -- event-local
P/H/R cells, Cayley/reset residual coordinates, A21 projection graph sectors,
and suffix/prefix selectors -- without detaching them from execution history.

The bridge intentionally stops before claiming a source-uniform nonlinear graph
or P4.  The hard SEA3 provider remains open upstream.
"""
from __future__ import annotations

import argparse
import copy
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval
import ou3_p4_complete_sea3_joint_sector_master as JOINT
import ou3_sea3_complete_window_execution_kernel as KERNEL
import ou3_sea3_frontend_state_step as FRONTEND
import ou3_sea3_full_normal_live_word as WORD
import ou3_sea3_tuner_scheduler_step as TUNER

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_P4_COMPLETE_SEA3_SAME_HISTORY_PREFIX_SELECTORS_V1"


@dataclass(frozen=True)
class PrefixSelector:
    """One retained child of one exact same-history sample transition.

    Snapshot fields are deep copies.  They are proof inputs and must be treated
    as read-only by downstream assemblers.
    """

    sample_index: int
    prefix_length: int
    parent_branch_ordinal: int
    successor_ordinal: int
    parent_source_cell_id: str
    source_cell_id: str
    sample_coordinates: KERNEL.SampleCoordinates
    active_schedule: TUNER.ActiveSchedule
    actual_rs_std_xyz: tuple[Interval, Interval, Interval]
    H_before: WORD.LiteralWordState
    H_after: WORD.LiteralWordState
    A_before: WORD.LiteralWordState
    A_after: WORD.LiteralWordState
    H_events_this_sample: tuple[str, ...]
    A_events_this_sample: tuple[str, ...]
    H_floor_case: str | None
    A_floor_case: str | None


def _active_schedule_and_rs(
    branch: KERNEL.ExecutionBranch,
    constants: KERNEL.KernelConstants,
) -> tuple[TUNER.ActiveSchedule, tuple[Interval, Interval, Interval]]:
    """Read the shipping schedule that advance_branch will use this sample."""
    committed = TUNER.commit_if_pending(branch.frontend.tuner, constants.tuner)
    active = committed.active
    rs = tuple(TUNER.active_rs_std_xyz(active, constants.tuner))
    if len(rs) != 3:
        raise RuntimeError("shipping actual R_S must have exactly three axes")
    return active, rs  # type: ignore[return-value]


def execute_with_prefix_selectors(
    *,
    frontend_entry: FRONTEND.FrontEndState,
    P0_H,
    P0_A,
    samples: Sequence[KERNEL.SampleCoordinates],
    domain_path: Path = DEFAULT_DOMAIN,
    branch_limit: int = 100000,
) -> tuple[list[KERNEL.ExecutionBranch], list[PrefixSelector], dict]:
    """Execute the trusted kernel while retaining every branch-correlated prefix.

    This intentionally mirrors only the branch-retention loop.  The physical
    and Riccati transition itself is not reimplemented: ``advance_branch`` is
    the sole transition primitive.
    """
    constants = KERNEL._process_constants(domain_path)
    branches = [
        KERNEL.ExecutionBranch(
            frontend=copy.deepcopy(frontend_entry),
            H=WORD.initialize_word("H", copy.deepcopy(P0_H)),
            A=WORD.initialize_word("A", copy.deepcopy(P0_A)),
            source_cell_id="root",
        )
    ]
    selectors: list[PrefixSelector] = []
    floor_cases: dict[str, int] = {}
    max_branches = 1
    all_frontend_before_mag = True

    for k, sample in enumerate(samples):
        next_branches: list[KERNEL.ExecutionBranch] = []
        for j, branch in enumerate(branches):
            active, rs_xyz = _active_schedule_and_rs(branch, constants)
            H_before = copy.deepcopy(branch.H)
            A_before = copy.deepcopy(branch.A)

            successors, meta = KERNEL.advance_branch(
                branch,
                sample,
                constants=constants,
                next_cell_prefix=f"k{k}:parent{j}",
            )
            if meta.get("same_active_schedule_verified") is not True:
                raise RuntimeError("trusted kernel did not verify same active schedule")
            if meta.get("same_actual_RS_verified") is not True:
                raise RuntimeError("trusted kernel did not verify actual R_S provenance")
            if not successors:
                raise RuntimeError("trusted kernel transition produced no successor")

            for i, child in enumerate(successors):
                H_events = tuple(child.H.event_log[len(H_before.event_log):])
                A_events = tuple(child.A.event_log[len(A_before.event_log):])
                selectors.append(
                    PrefixSelector(
                        sample_index=k,
                        prefix_length=k + 1,
                        parent_branch_ordinal=j,
                        successor_ordinal=i,
                        parent_source_cell_id=branch.source_cell_id,
                        source_cell_id=child.source_cell_id,
                        sample_coordinates=copy.deepcopy(sample),
                        active_schedule=copy.deepcopy(active),
                        actual_rs_std_xyz=copy.deepcopy(rs_xyz),
                        H_before=copy.deepcopy(H_before),
                        H_after=copy.deepcopy(child.H),
                        A_before=copy.deepcopy(A_before),
                        A_after=copy.deepcopy(child.A),
                        H_events_this_sample=H_events,
                        A_events_this_sample=A_events,
                        H_floor_case=meta.get("H_floor_case"),
                        A_floor_case=meta.get("A_floor_case"),
                    )
                )
            next_branches.extend(successors)
            all_frontend_before_mag = all_frontend_before_mag and bool(
                meta["frontend_completed_before_async_mag"]
            )
            for key in ("H_floor_case", "A_floor_case"):
                case = meta[key]
                if case is not None:
                    floor_cases[case] = floor_cases.get(case, 0) + 1

        if len(next_branches) > branch_limit:
            raise RuntimeError(
                f"front-end branch count {len(next_branches)} exceeds validated execution limit {branch_limit}; "
                "provider must partition the source cell more tightly, not select a successor"
            )
        branches = next_branches
        max_branches = max(max_branches, len(branches))

    return branches, selectors, {
        "samples_executed": len(samples),
        "endpoint_branches": len(branches),
        "prefix_selectors": len(selectors),
        "max_branch_count": max_branches,
        "floor_cases": floor_cases,
        "same_word_executed_H18_A21": True,
        "frontend_completed_before_async_mag": all_frontend_before_mag,
        "favorable_frontend_successor_selected": False,
        "shipping_transition_reimplemented": False,
        "trusted_advance_branch_is_only_transition_primitive": True,
        "kernel_self_test_only_not_P4": True,
    }


def selector_index(selectors: Sequence[PrefixSelector]) -> dict[str, PrefixSelector]:
    """Index unique child cells for ancestry and prefix lookup."""
    out: dict[str, PrefixSelector] = {}
    for selector in selectors:
        if selector.source_cell_id == "root":
            raise ValueError("root is not a child prefix selector")
        if selector.source_cell_id in out:
            raise ValueError(f"duplicate source cell id {selector.source_cell_id}")
        out[selector.source_cell_id] = selector
    return out


def lineage_for_endpoint(
    selectors: Sequence[PrefixSelector],
    endpoint_source_cell_id: str,
) -> list[PrefixSelector]:
    """Return root-to-endpoint selector lineage for one retained endpoint."""
    by_child = selector_index(selectors)
    if endpoint_source_cell_id not in by_child:
        raise ValueError(f"unknown endpoint source cell {endpoint_source_cell_id}")
    reverse: list[PrefixSelector] = []
    cell = endpoint_source_cell_id
    while cell != "root":
        selector = by_child.get(cell)
        if selector is None:
            raise ValueError(f"broken same-history ancestry at {cell}")
        reverse.append(selector)
        cell = selector.parent_source_cell_id
    lineage = list(reversed(reverse))
    for expected, selector in enumerate(lineage, start=1):
        if selector.prefix_length != expected:
            raise ValueError("same-history lineage skipped or duplicated a prefix")
    return lineage


def validate_selector_graph(
    selectors: Sequence[PrefixSelector],
    *,
    samples_executed: int,
    endpoint_source_cell_ids: Sequence[str],
) -> list[str]:
    """Validate ancestry/event/source invariants without claiming source closure."""
    failures: list[str] = []
    try:
        by_child = selector_index(selectors)
    except ValueError as exc:
        return [str(exc)]

    for selector in selectors:
        if selector.prefix_length != selector.sample_index + 1:
            failures.append(f"{selector.source_cell_id}: prefix/sample index mismatch")
        if selector.sample_index == 0:
            if selector.parent_source_cell_id != "root":
                failures.append(f"{selector.source_cell_id}: first prefix parent is not root")
        else:
            parent = by_child.get(selector.parent_source_cell_id)
            if parent is None:
                failures.append(f"{selector.source_cell_id}: missing parent selector")
            elif parent.prefix_length + 1 != selector.prefix_length:
                failures.append(f"{selector.source_cell_id}: parent is not previous prefix")

        for mode, events in (
            ("H18", selector.H_events_this_sample),
            ("A21", selector.A_events_this_sample),
        ):
            if not events or events[0] != "prediction":
                failures.append(f"{selector.source_cell_id}: {mode} sample does not begin with prediction")
            if "accelerometer" not in events:
                failures.append(f"{selector.source_cell_id}: {mode} sample lost accelerometer")
            if selector.sample_coordinates.due_S and "S_zero" not in events:
                failures.append(f"{selector.source_cell_id}: {mode} due S event missing")
            if not selector.sample_coordinates.due_S and "S_zero" in events:
                failures.append(f"{selector.source_cell_id}: {mode} inserted non-due S event")
            mag_count = len(selector.sample_coordinates.magnetometer_events_after_imu)
            if events.count("magnetometer") != mag_count:
                failures.append(f"{selector.source_cell_id}: {mode} async magnetometer count changed")
            if mag_count and events.index("magnetometer") < events.index("accelerometer"):
                failures.append(f"{selector.source_cell_id}: {mode} magnetometer moved before accelerometer")

        if selector.H_after.imu_samples != selector.H_before.imu_samples + 1:
            failures.append(f"{selector.source_cell_id}: H18 sample counter mismatch")
        if selector.A_after.imu_samples != selector.A_before.imu_samples + 1:
            failures.append(f"{selector.source_cell_id}: A21 sample counter mismatch")
        if selector.sample_coordinates.due_S:
            if selector.H_after.S_updates != selector.H_before.S_updates + 1:
                failures.append(f"{selector.source_cell_id}: H18 due S counter mismatch")
            if selector.A_after.S_updates != selector.A_before.S_updates + 1:
                failures.append(f"{selector.source_cell_id}: A21 due S counter mismatch")

    for endpoint in endpoint_source_cell_ids:
        try:
            lineage = lineage_for_endpoint(selectors, endpoint)
        except ValueError as exc:
            failures.append(str(exc))
            continue
        if len(lineage) != samples_executed:
            failures.append(f"{endpoint}: endpoint lineage does not contain every sample")
        elif lineage and lineage[-1].source_cell_id != endpoint:
            failures.append(f"{endpoint}: lineage endpoint changed")

    return list(dict.fromkeys(failures))


def _point_sample() -> KERNEL.SampleCoordinates:
    z = KERNEL.MAHONY.I(0.0)
    return KERNEL.SampleCoordinates(
        gyro_measurement=KERNEL.MAHONY.Vec3(
            KERNEL.MAHONY.I(0.01),
            KERNEL.MAHONY.I(-0.02),
            KERNEL.MAHONY.I(0.005),
        ),
        omega_body_corrected=(
            KERNEL.PRED.I(0.01),
            KERNEL.PRED.I(-0.02),
            KERNEL.PRED.I(0.005),
        ),
        specific_force=KERNEL.MAHONY.Vec3(
            KERNEL.MAHONY.I(0.2),
            KERNEL.MAHONY.I(-0.1),
            KERNEL.MAHONY.I(-9.75),
        ),
        f_cog_body=(
            KERNEL.PRED.I(0.0),
            KERNEL.PRED.I(0.0),
            KERNEL.PRED.I(-9.80665),
        ),
        R_wb=[
            [KERNEL.PRED.I(1.0), z, z],
            [z, KERNEL.PRED.I(1.0), z],
            [z, z, KERNEL.PRED.I(1.0)],
        ],
        due_S=True,
        aw_floor_requested=True,
        magnetometer_events_after_imu=(
            KERNEL.MagneticEvent(
                (
                    KERNEL.PRED.I(20.0),
                    KERNEL.PRED.I(0.0),
                    KERNEL.PRED.I(40.0),
                )
            ),
        ),
    )


def _smoke(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    samples = [_point_sample(), _point_sample()]
    endpoints, selectors, meta = execute_with_prefix_selectors(
        frontend_entry=FRONTEND._point_state(),
        P0_H=KERNEL._diag_P(18, 2.0),
        P0_A=KERNEL._diag_P(21, 2.0),
        samples=samples,
        domain_path=domain_path,
        branch_limit=128,
    )
    endpoint_ids = [branch.source_cell_id for branch in endpoints]
    failures = validate_selector_graph(
        selectors,
        samples_executed=len(samples),
        endpoint_source_cell_ids=endpoint_ids,
    )
    lineages = [lineage_for_endpoint(selectors, endpoint) for endpoint in endpoint_ids]
    return {
        **meta,
        "selector_graph_failures": failures,
        "selector_graph_valid": not failures,
        "all_endpoint_lineages_cover_every_prefix": bool(lineages)
        and all(len(lineage) == len(samples) for lineage in lineages),
        "all_selectors_retain_actual_RS": bool(selectors)
        and all(len(s.actual_rs_std_xyz) == 3 for s in selectors),
        "all_selectors_retain_H18_A21_before_after": bool(selectors)
        and all(
            s.H_before.mode == "H"
            and s.H_after.mode == "H"
            and s.A_before.mode == "A"
            and s.A_after.mode == "A"
            for s in selectors
        ),
        "all_selectors_retain_exact_event_slices": bool(selectors)
        and all(
            "prediction" in s.H_events_this_sample
            and "accelerometer" in s.H_events_this_sample
            and "prediction" in s.A_events_this_sample
            and "accelerometer" in s.A_events_this_sample
            for s in selectors
        ),
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    kernel = KERNEL.build(domain_path)
    kf = KERNEL.validate(kernel)
    if kf:
        raise RuntimeError(f"trusted typed execution kernel invalid: {kf}")
    joint = JOINT.build()
    jf = JOINT.validate(joint)
    if jf:
        raise RuntimeError(f"joint-sector master prerequisite invalid: {jf}")
    smoke = _smoke(domain_path)
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "P3_delta_preserved": 1.0e-18,
        "source_generator": False,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "quality_gates_changed": False,
        "declared_domain_changed": False,
        "canonical_provider_gate_bypassed": False,
        "trusted_typed_kernel_consumed": True,
        "shipping_transition_reimplemented": False,
        "branch_correlated_every_prefix_selector_available": True,
        "parent_child_source_cell_ancestry_retained": True,
        "same_provider_sample_coordinates_retained_per_transition": True,
        "same_committed_active_schedule_retained_per_transition": True,
        "actual_applied_anisotropic_RS_retained_per_transition": True,
        "H18_A21_before_after_literal_states_retained_per_transition": True,
        "exact_shipping_event_slice_retained_per_transition": True,
        "joint_sector_master_consumed_without_promotion": True,
        "source_uniform_provider_family_closed_here": False,
        "event_local_same_P_H_R_cells_materialized_here": False,
        "nonlinear_residual_history_graph_materialized_here": False,
        "A21_projection_graph_attached_here": False,
        "source_uniform_endpoint_joint_sector_closed_here": False,
        "source_uniform_every_prefix_joint_sector_closed_here": False,
        "P4_promoted_here": False,
        "smoke": smoke,
        "next_obligation": (
            "consume these branch-correlated prefix selectors to materialize event-local same-P/H/R "
            "prediction/Joseph/reset/projection graph cells and the joint nonlinear residual history "
            "for endpoint plus every prefix; then run the full augmented outward LDLT on the same lineage"
        ),
    }


def validate(d: dict) -> list[str]:
    failures: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        failures.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        failures.append("canonical source changed")
    if float(d.get("P3_delta_preserved", 0.0)) != 1.0e-18:
        failures.append("frozen P3 delta changed")
    for key in (
        "trusted_typed_kernel_consumed",
        "branch_correlated_every_prefix_selector_available",
        "parent_child_source_cell_ancestry_retained",
        "same_provider_sample_coordinates_retained_per_transition",
        "same_committed_active_schedule_retained_per_transition",
        "actual_applied_anisotropic_RS_retained_per_transition",
        "H18_A21_before_after_literal_states_retained_per_transition",
        "exact_shipping_event_slice_retained_per_transition",
        "joint_sector_master_consumed_without_promotion",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in (
        "source_generator",
        "trajectory_replay_used",
        "filter_changed",
        "quality_gates_changed",
        "declared_domain_changed",
        "canonical_provider_gate_bypassed",
        "shipping_transition_reimplemented",
        "source_uniform_provider_family_closed_here",
        "event_local_same_P_H_R_cells_materialized_here",
        "nonlinear_residual_history_graph_materialized_here",
        "A21_projection_graph_attached_here",
        "source_uniform_endpoint_joint_sector_closed_here",
        "source_uniform_every_prefix_joint_sector_closed_here",
        "P4_promoted_here",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    smoke = d.get("smoke", {})
    for key in (
        "same_word_executed_H18_A21",
        "frontend_completed_before_async_mag",
        "trusted_advance_branch_is_only_transition_primitive",
        "selector_graph_valid",
        "all_endpoint_lineages_cover_every_prefix",
        "all_selectors_retain_actual_RS",
        "all_selectors_retain_H18_A21_before_after",
        "all_selectors_retain_exact_event_slices",
    ):
        if smoke.get(key) is not True:
            failures.append(f"smoke lost {key}")
    for key in (
        "favorable_frontend_successor_selected",
        "shipping_transition_reimplemented",
    ):
        if smoke.get(key) is not False:
            failures.append(f"smoke changed {key}")
    if int(smoke.get("samples_executed", 0)) != 2:
        failures.append("smoke did not retain two complete prefixes")
    if int(smoke.get("prefix_selectors", 0)) < 2:
        failures.append("smoke did not materialize every prefix selector")
    if smoke.get("selector_graph_failures") != []:
        failures.append("selector graph smoke reported failures")
    return list(dict.fromkeys(failures))


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
    print(
        json.dumps(
            {
                "prefix_selectors_available": d[
                    "branch_correlated_every_prefix_selector_available"
                ],
                "actual_RS_retained": d[
                    "actual_applied_anisotropic_RS_retained_per_transition"
                ],
                "P4_promoted_here": d["P4_promoted_here"],
                "smoke": d["smoke"],
                "failures": failures,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
