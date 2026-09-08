#!/usr/bin/env python3
"""Branch-correlated every-prefix selectors for complete-SEA3 OU-III P4.

P4 needs the endpoint and every finite prefix of ONE correlated complete-SEA3
history.  Reconstructing tuner schedules, actual R_S, Riccati covariances, or
Joseph gains from independent boxes would destroy exactly the source
correlation the theorem must retain.

This non-promoting bridge executes every transition through
``ou3_sea3_complete_window_execution_kernel.advance_branch`` and records, for
each retained child prefix:

* explicit parent/child source-cell ancestry;
* the exact provider sample coordinates;
* the committed tau/sigma/pseudo-period schedule and actual anisotropic R_S;
* H18/A21 literal-word states before and after the sample;
* the exact shipping event slice; and
* passive event-local Riccati cells captured inside that same trusted
  transition: P-before/P-after plus F/Q, floor increment, or Joseph H/R.

The bridge still does not create the hard SEA3 provider family, nonlinear
residual/projection graph sectors, or an augmented endpoint/prefix LDLT.
Therefore it cannot promote P4.
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
import ou3_sea3_live_covariance_seed as LIVE
import ou3_sea3_tuner_scheduler_step as TUNER

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 2
QUALIFICATION = "OU3_P4_COMPLETE_SEA3_SAME_HISTORY_PREFIX_SELECTORS_V2"


@dataclass(frozen=True)
class PrefixSelector:
    """One retained child of one exact same-history sample transition."""

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
    H_event_cells: tuple[KERNEL.RiccatiEventCell, ...]
    A_event_cells: tuple[KERNEL.RiccatiEventCell, ...]
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

    Only the branch-retention loop lives here.  The physical and Riccati
    transition is not reimplemented.  Event-local cells are passive snapshots
    emitted by the same ``advance_branch`` call that creates the child state.
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
                capture_riccati_event_cells=True,
            )
            if meta.get("same_active_schedule_verified") is not True:
                raise RuntimeError("trusted kernel did not verify same active schedule")
            if meta.get("same_actual_RS_verified") is not True:
                raise RuntimeError("trusted kernel did not verify actual R_S provenance")
            if meta.get("riccati_event_cells_captured") is not True:
                raise RuntimeError("trusted kernel did not capture event-local Riccati cells")
            if not successors:
                raise RuntimeError("trusted kernel transition produced no successor")

            H_cells = tuple(copy.deepcopy(meta["H_event_cells"]))
            A_cells = tuple(copy.deepcopy(meta["A_event_cells"]))
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
                        H_event_cells=copy.deepcopy(H_cells),
                        A_event_cells=copy.deepcopy(A_cells),
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
        "event_local_cells_captured_inside_same_transition": True,
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
    visited: set[str] = set()
    while cell != "root":
        if cell in visited:
            raise ValueError("cycle in same-history ancestry")
        visited.add(cell)
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


def _validate_event_cells(
    selector: PrefixSelector,
    mode: str,
    events: tuple[str, ...],
    cells: tuple[KERNEL.RiccatiEventCell, ...],
) -> list[str]:
    failures: list[str] = []
    label = "H18" if mode == "H" else "A21"
    prefix = f"{selector.source_cell_id}: {label}"

    if tuple(cell.kind for cell in cells) != events:
        failures.append(f"{prefix} event-cell order detached from shipping event slice")
    if any(cell.mode != mode for cell in cells):
        failures.append(f"{prefix} event-cell mode mismatch")
    if any(cell.event_index_in_sample != i for i, cell in enumerate(cells)):
        failures.append(f"{prefix} event-cell sample index is not contiguous")

    measurement_kinds = {"S_zero", "accelerometer", "magnetometer"}
    for cell in cells:
        if not cell.P_before or not cell.P_after:
            failures.append(f"{prefix} {cell.kind} lost covariance snapshots")
        if cell.kind == "prediction":
            if cell.F is None or cell.Q is None:
                failures.append(f"{prefix} prediction lost same-event F/Q")
        elif cell.kind == "aw_floor":
            if cell.floor_increment is None:
                failures.append(f"{prefix} floor lost covariance increment")
        elif cell.kind in measurement_kinds:
            if cell.H is None or cell.R is None:
                failures.append(f"{prefix} {cell.kind} lost same-event P/H/R")

    expected_s_R = WORD.R_S_zero(selector.actual_rs_std_xyz)
    s_cells = [cell for cell in cells if cell.kind == "S_zero"]
    expected_s_count = 1 if selector.sample_coordinates.due_S else 0
    if len(s_cells) != expected_s_count:
        failures.append(f"{prefix} captured due-S count changed")
    for cell in s_cells:
        if cell.actual_rs_from_committed_schedule is not True:
            failures.append(f"{prefix} S cell lost actual committed R_S provenance")
        if cell.R != expected_s_R:
            failures.append(f"{prefix} S cell R differs from exact actual applied R_S")
        if cell.H != WORD.H_S_zero(mode):
            failures.append(f"{prefix} S cell H differs from shipping H_S")

    acc_cells = [cell for cell in cells if cell.kind == "accelerometer"]
    if len(acc_cells) != 1:
        failures.append(f"{prefix} captured accelerometer count changed")
    elif acc_cells[0].H != WORD.H_accelerometer(
        mode,
        selector.sample_coordinates.f_cog_body,
        selector.sample_coordinates.R_wb,
    ):
        failures.append(f"{prefix} accelerometer H detached from provider geometry")

    mag_cells = [cell for cell in cells if cell.kind == "magnetometer"]
    magnetic_events = selector.sample_coordinates.magnetometer_events_after_imu
    if len(mag_cells) != len(magnetic_events):
        failures.append(f"{prefix} captured magnetometer count changed")
    else:
        for i, (cell, event) in enumerate(zip(mag_cells, magnetic_events)):
            if cell.magnetic_event_index != i:
                failures.append(f"{prefix} magnetometer event index changed")
            if cell.H != WORD.H_magnetometer(mode, event.m_body):
                failures.append(f"{prefix} magnetometer H detached from provider vector")

    if cells:
        initial_P = selector.H_before.riccati.P if mode == "H" else selector.A_before.riccati.P
        final_P = selector.H_after.riccati.P if mode == "H" else selector.A_after.riccati.P
        if cells[0].P_before != initial_P:
            failures.append(f"{prefix} first captured covariance differs from parent word")
        for previous, current in zip(cells, cells[1:]):
            if previous.P_after != current.P_before:
                failures.append(f"{prefix} covariance continuity lost before {current.kind}")
        if cells[-1].P_after != final_P:
            failures.append(f"{prefix} final captured covariance differs from child word")
    return failures


def validate_selector_graph(
    selectors: Sequence[PrefixSelector],
    *,
    samples_executed: int,
    endpoint_source_cell_ids: Sequence[str],
) -> list[str]:
    """Validate ancestry/event/source invariants without claiming source closure."""
    failures: list[str] = []
    if samples_executed <= 0 or not selectors or not endpoint_source_cell_ids:
        return ["nonempty executed prefix family and endpoints required"]
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
            elif (
                selector.H_before != parent.H_after
                or selector.A_before != parent.A_after
            ):
                failures.append(f"{selector.source_cell_id}: word state detached from parent prefix")

        for mode, events, cells in (
            ("H", selector.H_events_this_sample, selector.H_event_cells),
            ("A", selector.A_events_this_sample, selector.A_event_cells),
        ):
            label = "H18" if mode == "H" else "A21"
            if not events or events[0] != "prediction":
                failures.append(
                    f"{selector.source_cell_id}: {label} sample does not begin with prediction"
                )
            if "accelerometer" not in events:
                failures.append(f"{selector.source_cell_id}: {label} sample lost accelerometer")
            if selector.sample_coordinates.due_S and "S_zero" not in events:
                failures.append(f"{selector.source_cell_id}: {label} due S event missing")
            if not selector.sample_coordinates.due_S and "S_zero" in events:
                failures.append(f"{selector.source_cell_id}: {label} inserted non-due S event")
            mag_count = len(selector.sample_coordinates.magnetometer_events_after_imu)
            if events.count("magnetometer") != mag_count:
                failures.append(
                    f"{selector.source_cell_id}: {label} async magnetometer count changed"
                )
            if mag_count and events.index("magnetometer") < events.index("accelerometer"):
                failures.append(
                    f"{selector.source_cell_id}: {label} magnetometer moved before accelerometer"
                )
            failures.extend(_validate_event_cells(selector, mode, events, cells))

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


def _live_structured_point_covariance_fixture(
    frontend_entry: FRONTEND.FrontEndState,
    domain_path: Path,
):
    """Build a Live-structured covariance fixture, not a reachable A21 entry.

    The previous 2*I fixture was an arbitrary PSD matrix and caused enclosure
    loss on the second Joseph prefix.  This fixture instead consumes the
    shipping Live seed contract: full-heading tilt/yaw handoff covariance,
    constructor b_g/v/p/S seeds, a_w reset to the committed stationary
    covariance, and the held->active b_a diagonal release floor.  The point
    geometry is identity, so the body-down projector is the z axis and the
    attitude seed is diag(tilt^2, tilt^2, yaw^2).

    The H block uses the Live seed structure, while the A block appends the
    shipping bias-floor variance solely to exercise the full 21-state code.
    No pre-release H history or rectangular H->A event is executed here.
    Neither this fixture nor the point frontend proves SEA3 reachability.
    """
    live = LIVE.build(Path(domain_path).resolve())
    lf = LIVE.validate(live)
    if lf:
        raise RuntimeError(f"shipping Live covariance seed invalid: {lf}")

    constants = KERNEL._process_constants(Path(domain_path).resolve())
    active = TUNER.commit_if_pending(frontend_entry.tuner, constants.tuner).active
    if active.sigma.lo <= 0.0:
        raise RuntimeError("point frontend active sigma is not strictly positive")

    def zero_matrix(n: int):
        z = KERNEL.PRED.I(0.0)
        return [[z for _ in range(n)] for _ in range(n)]

    def set_diag(P, offset: int, values):
        for i, value in enumerate(values):
            P[offset + i][offset + i] = value

    tilt = KERNEL.PRED.I(float(live["full_heading_gauged_live_attitude_seed"]["tilt_std_rad"]))
    yaw = KERNEL.PRED.I(float(live["full_heading_gauged_live_attitude_seed"]["yaw_std_rad"]))
    tilt_var = tilt * tilt
    yaw_var = yaw * yaw
    bg_var = KERNEL.PRED.I(float(live["constructor"]["P_bg_variance"]))
    v_var = KERNEL.PRED.I(float(live["translation_seed"]["P_v"]))
    p_var = KERNEL.PRED.I(float(live["translation_seed"]["P_p"]))
    s_var = KERNEL.PRED.I(float(live["translation_seed"]["P_S"]))
    aw_var = active.sigma * active.sigma
    ba_std = KERNEL.PRED.I(float(live["constructor"]["sigma_ba0"]))
    ba_var = ba_std * ba_std

    P_H = zero_matrix(WORD.H_DIM)
    set_diag(P_H, WORD.OFF_TH, (tilt_var, tilt_var, yaw_var))
    set_diag(P_H, WORD.OFF_BG, (bg_var, bg_var, bg_var))
    set_diag(P_H, WORD.OFF_V, (v_var, v_var, v_var))
    set_diag(P_H, WORD.OFF_P, (p_var, p_var, p_var))
    set_diag(P_H, WORD.OFF_S, (s_var, s_var, s_var))
    set_diag(P_H, WORD.OFF_AW, (aw_var, aw_var, aw_var))

    P_A = zero_matrix(WORD.A_DIM)
    for i in range(WORD.H_DIM):
        for j in range(WORD.H_DIM):
            P_A[i][j] = P_H[i][j]
    set_diag(P_A, WORD.OFF_BA, (ba_var, ba_var, ba_var))

    return P_H, P_A, {
        "live_seed_contract_consumed": True,
        "arbitrary_2I_covariance_used": False,
        "identity_point_attitude_seed_uses_tilt_yaw_split": True,
        "aw_seed_uses_same_committed_sigma": True,
        "A21_ba_floor_variance_attached_for_fixture_only": True,
        "same_history_H_to_A_release_executed": False,
        "A21_entry_reachability_certified": False,
        "complete_SEA3_source_membership_certified": False,
    }


def _prefix_evidence(selector: PrefixSelector) -> dict:
    """Small numerical CI record from the captured cells, not a certificate."""
    result = {
        "prefix_length": selector.prefix_length,
        "parent_source_cell_id": selector.parent_source_cell_id,
        "source_cell_id": selector.source_cell_id,
    }
    for mode, word, cells in (
        ("H18", selector.H_after, selector.H_event_cells),
        ("A21", selector.A_after, selector.A_event_cells),
    ):
        P = word.riccati.P
        result[mode] = {
            "dimension": len(P),
            "events": [cell.kind for cell in cells],
            "P_after_diagonal": [[P[i][i].lo, P[i][i].hi] for i in range(len(P))],
            "P_after_max_entry_width": max(x.hi - x.lo for row in P for x in row),
            "actual_R_S_diagonals": [
                [[cell.R[i][i].lo, cell.R[i][i].hi] for i in range(3)]
                for cell in cells if cell.kind == "S_zero"
            ],
        }
    return result


def _smoke(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    samples = [_point_sample(), _point_sample()]
    frontend_entry = FRONTEND._point_state()
    P0_H, P0_A, seed_meta = _live_structured_point_covariance_fixture(
        frontend_entry, domain_path
    )
    endpoints, selectors, meta = execute_with_prefix_selectors(
        frontend_entry=frontend_entry,
        P0_H=P0_H,
        P0_A=P0_A,
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
    all_cells = [cell for s in selectors for cell in s.H_event_cells + s.A_event_cells]
    measurement_cells = [
        cell for cell in all_cells if cell.kind in ("S_zero", "accelerometer", "magnetometer")
    ]
    return {
        **meta,
        "live_structured_covariance_fixture": seed_meta,
        "evidence_scope": "TWO_SAMPLE_FIXTURE_NOT_COMPLETE_SEA3_FAMILY",
        "prefix_evidence": [_prefix_evidence(s) for s in selectors],
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
            tuple(cell.kind for cell in s.H_event_cells) == s.H_events_this_sample
            and tuple(cell.kind for cell in s.A_event_cells) == s.A_events_this_sample
            for s in selectors
        ),
        "all_measurement_event_cells_retain_same_P_H_R": bool(measurement_cells)
        and all(cell.P_before and cell.P_after and cell.H is not None and cell.R is not None for cell in measurement_cells),
        "all_due_S_cells_retain_actual_committed_RS": bool(
            [cell for cell in all_cells if cell.kind == "S_zero"]
        )
        and all(
            cell.actual_rs_from_committed_schedule
            for cell in all_cells
            if cell.kind == "S_zero"
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
        "event_local_cells_captured_inside_same_trusted_transition": True,
        "event_local_same_P_H_R_cells_materialized_on_typed_execution": True,
        "event_local_due_S_cells_use_exact_actual_committed_RS": True,
        "joint_sector_master_consumed_without_promotion": True,
        "source_uniform_provider_family_closed_here": False,
        "source_uniform_event_local_same_P_H_R_cells_closed_here": False,
        "nonlinear_residual_history_graph_materialized_here": False,
        "A21_projection_graph_attached_here": False,
        "source_uniform_endpoint_joint_sector_closed_here": False,
        "source_uniform_every_prefix_joint_sector_closed_here": False,
        "P4_promoted_here": False,
        "smoke": smoke,
        "next_obligation": (
            "feed these same-history event-local P/H/R cells plus source residual/true-bias coordinates "
            "into the exact nonlinear Joseph/reset/projection graph, assemble the full augmented master "
            "at endpoint and every prefix, and only then attempt source-uniform outward LDLT"
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
        "event_local_cells_captured_inside_same_trusted_transition",
        "event_local_same_P_H_R_cells_materialized_on_typed_execution",
        "event_local_due_S_cells_use_exact_actual_committed_RS",
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
        "source_uniform_event_local_same_P_H_R_cells_closed_here",
        "nonlinear_residual_history_graph_materialized_here",
        "A21_projection_graph_attached_here",
        "source_uniform_endpoint_joint_sector_closed_here",
        "source_uniform_every_prefix_joint_sector_closed_here",
        "P4_promoted_here",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    smoke = d.get("smoke", {})
    seed = smoke.get("live_structured_covariance_fixture", {})
    for key in (
        "live_seed_contract_consumed",
        "identity_point_attitude_seed_uses_tilt_yaw_split",
        "aw_seed_uses_same_committed_sigma",
        "A21_ba_floor_variance_attached_for_fixture_only",
    ):
        if seed.get(key) is not True:
            failures.append(f"smoke covariance fixture lost {key}")
    for key in (
        "arbitrary_2I_covariance_used",
        "same_history_H_to_A_release_executed",
        "A21_entry_reachability_certified",
        "complete_SEA3_source_membership_certified",
    ):
        if seed.get(key) is not False:
            failures.append(f"smoke covariance fixture changed {key}")
    if smoke.get("evidence_scope") != "TWO_SAMPLE_FIXTURE_NOT_COMPLETE_SEA3_FAMILY":
        failures.append("smoke was relabeled as a complete SEA3 source family")
    if len(smoke.get("prefix_evidence", [])) != smoke.get("prefix_selectors"):
        failures.append("numerical CI evidence does not cover every retained prefix")
    for key in (
        "same_word_executed_H18_A21",
        "frontend_completed_before_async_mag",
        "trusted_advance_branch_is_only_transition_primitive",
        "event_local_cells_captured_inside_same_transition",
        "selector_graph_valid",
        "all_endpoint_lineages_cover_every_prefix",
        "all_selectors_retain_actual_RS",
        "all_selectors_retain_H18_A21_before_after",
        "all_selectors_retain_exact_event_slices",
        "all_measurement_event_cells_retain_same_P_H_R",
        "all_due_S_cells_retain_actual_committed_RS",
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
                "same_P_H_R_cells": d[
                    "event_local_same_P_H_R_cells_materialized_on_typed_execution"
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
