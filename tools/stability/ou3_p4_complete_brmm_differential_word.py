#!/usr/bin/env python3
"""Literal full-state differential cocycle for complete-BRMM OU-III P4.

This is the matrix object consumed by canonical P4. It does not generate a
source word or replace shipping by independent event boxes. A caller supplies
outward event Jacobian enclosures produced from ONE same-history complete BRMM
source path. Measurement events must come from the same-cell P/H/R constructor
in ``ou3_p4_complete_brmm_differential_events``; a bare arbitrary K/J matrix
cannot be labeled into the theorem path.

Required fixed-mode IMU order is

    prediction -> optional a_w covariance floor (state Jacobian I)
               -> every due S=0 Joseph/reset
               -> required accelerometer Joseph/reset,

with asynchronous accepted vector/magnetometer Joseph/reset events inserted at
their actual source positions. The H18->A21 release is a separate hybrid event.
Its homogeneous derivative is exactly ``[I_18;0]``: the already-active H18
coordinates are continuous, while the newly activated b_a coordinate has no
derivative with respect to H18. The held b_a error is NOT set to zero; it is a
separate hybrid/ISS forcing term. The shipping b_a covariance seed floor is
likewise a separate covariance/metric event. Both facts are mandatory metadata
on the release event.

Every S event must carry ``ACTUAL_APPLIED_SPECTRALMSE_RS`` provenance from the
same source cell. If J_k is the exact/outward differential map of event k,

    J_W = J_{N-1} ... J_0.

No packet norm, packet count, scalar beta, selected-S replacement word, or
independent R_S schedule appears. ``certify_word`` sends the full cocycle to

    rho M_0 - J_W^T M_1 J_W > 0

through the full interval-LDLT routine of the Phi pullback differential metric.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import argparse
import json
from pathlib import Path
from typing import Any, Sequence

from ou3_interval import Interval, matrix_identity, matrix_mul
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_complete_brmm_phi_differential_metric as DIFF
import ou3_brmm_complete_source as COMPLETE
import ou3_brmm_full_normal_live_word as WORD
import ou3_brmm_riccati_metric_p3 as P3

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 3
QUALIFICATION = "OU3_P4_COMPLETE_BRMM_LITERAL_DIFFERENTIAL_WORD_V3"


def _shape(A: Sequence[Sequence[Any]]) -> tuple[int, int]:
    r = len(A)
    c = len(A[0]) if r else 0
    if any(len(row) != c for row in A):
        raise ValueError("ragged differential event matrix")
    return r, c


@dataclass(frozen=True)
class DifferentialEvent:
    kind: str
    J: list[list[Interval]]
    source_token: str
    R_provenance: str | None = None
    same_P_H_R_cell: bool = False
    accepted: bool = True
    held_ba_forcing_separate: bool = False
    release_covariance_floor_retained: bool = False


@dataclass
class DifferentialWord:
    start_mode: str
    source_token: str
    J_word: list[list[Interval]]
    current_dim: int
    events: list[str] = field(default_factory=list)
    predictions: int = 0
    floors: int = 0
    S_updates: int = 0
    accelerometer_updates: int = 0
    vector_updates: int = 0
    hybrid_lifts: int = 0
    held_ba_forcing_separate: bool = False
    release_covariance_floor_retained: bool = False

    @property
    def start_dim(self) -> int:
        return 18 if self.start_mode == "H" else 21


def initialize(mode: str, source_token: str) -> DifferentialWord:
    if mode not in ("H", "A"):
        raise ValueError("differential word mode must be H or A")
    if not isinstance(source_token, str) or not source_token:
        raise ValueError("same-history complete BRMM source token required")
    n = 18 if mode == "H" else 21
    return DifferentialWord(
        start_mode=mode,
        source_token=source_token,
        J_word=matrix_identity(n),
        current_dim=n,
    )


def event_from_source_joseph(kind: str, source_token: str, source_event: dict) -> DifferentialEvent:
    """Promote only the same-cell event object emitted by EVENTS.source_joseph_event."""
    if source_event.get("same_P_H_R_cell") is not True:
        raise ValueError("Joseph differential event did not derive K from the same P/H/R cell")
    return DifferentialEvent(
        kind=kind,
        J=source_event["J_state"],
        source_token=source_token,
        R_provenance=source_event.get("R_provenance"),
        same_P_H_R_cell=True,
    )


def H_to_A_release_event(source_token: str) -> DifferentialEvent:
    """Canonical homogeneous 21x18 release derivative plus required metadata.

    Held-mode b_a error is an exogenous bounded coordinate, not part of H18.
    Its finite value therefore enters the A21 state as a separate hybrid/ISS
    forcing and cannot appear as a derivative column of this homogeneous lift.
    The shipping diagonal b_a covariance floor changes only the metric/Riccati
    state and is retained separately from this physical-state Jacobian.
    """
    if not isinstance(source_token, str) or not source_token:
        raise ValueError("same-history complete BRMM source token required")
    J = [[Interval.point(0.0) for _ in range(18)] for _ in range(21)]
    for i in range(18):
        J[i][i] = Interval.point(1.0)
    return DifferentialEvent(
        kind="H_to_A",
        J=J,
        source_token=source_token,
        held_ba_forcing_separate=True,
        release_covariance_floor_retained=True,
    )


def _is_canonical_H_to_A_lift(J) -> bool:
    if _shape(J) != (21, 18):
        return False
    for i in range(21):
        for j in range(18):
            expected = Interval.point(1.0 if i == j and i < 18 else 0.0)
            if J[i][j] != expected:
                return False
    return True


def apply_event(word: DifferentialWord, event: DifferentialEvent) -> None:
    if event.source_token != word.source_token:
        raise ValueError("differential event detached from same complete BRMM source history")
    if not event.accepted:
        raise ValueError("only executed/accepted nonlinear events belong in the cocycle")
    r, c = _shape(event.J)
    kind = event.kind
    if kind == "H_to_A":
        if word.hybrid_lifts != 0 or word.current_dim != 18 or (r, c) != (21, 18):
            raise ValueError("H->A differential lift must be the unique 21x18 hybrid event")
    elif c != word.current_dim:
        raise ValueError("differential event input dimension does not match current word")

    if kind == "prediction":
        if r != c:
            raise ValueError("same-mode prediction must be square")
        word.predictions += 1
    elif kind == "aw_floor":
        if r != c:
            raise ValueError("covariance-floor state event must be square identity")
        I = matrix_identity(c)
        if any(event.J[i][j] != I[i][j] for i in range(c) for j in range(c)):
            raise ValueError("a_w covariance floor must have identity state Jacobian")
        word.floors += 1
    elif kind == "S_zero":
        if r != c:
            raise ValueError("S=0 Joseph/reset event must be same-mode square")
        if event.same_P_H_R_cell is not True:
            raise ValueError("S=0 event did not derive its gain from the same P/H/R cell")
        if event.R_provenance != EVENTS.ACTUAL_RS_PROVENANCE:
            raise ValueError("S=0 event lacks actual applied SpectralMSE R_S provenance")
        word.S_updates += 1
    elif kind == "accelerometer":
        if r != c:
            raise ValueError("accelerometer Joseph/reset event must be same-mode square")
        if event.same_P_H_R_cell is not True:
            raise ValueError("accelerometer event did not derive its gain from the same P/H/R cell")
        word.accelerometer_updates += 1
    elif kind in ("magnetometer", "vector"):
        if r != c:
            raise ValueError("vector Joseph/reset event must be same-mode square")
        if event.same_P_H_R_cell is not True:
            raise ValueError("vector event did not derive its gain from the same P/H/R cell")
        word.vector_updates += 1
    elif kind == "H_to_A":
        if not _is_canonical_H_to_A_lift(event.J):
            raise ValueError("H->A homogeneous derivative must be canonical [I18;0]")
        if event.held_ba_forcing_separate is not True:
            raise ValueError("H->A must retain held b_a error as a separate hybrid forcing")
        if event.release_covariance_floor_retained is not True:
            raise ValueError("H->A must retain the shipping b_a covariance floor separately")
        word.hybrid_lifts = 1
        word.held_ba_forcing_separate = True
        word.release_covariance_floor_retained = True
    else:
        raise ValueError(f"unknown differential event kind {kind}")

    word.J_word = matrix_mul(event.J, word.J_word)
    word.current_dim = r
    word.events.append(kind)


def certify_word(word, P0_inverse, P1_inverse, T0, T1, rho: float):
    if word.predictions <= 0 or word.accelerometer_updates <= 0:
        raise ValueError("complete differential word is missing prediction/accelerometer events")
    return DIFF.certify_differential_contraction(
        word.J_word, P0_inverse, P1_inverse, T0, T1, rho
    )


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    source = COMPLETE.build(path)
    literal = WORD.build(path)
    p3 = P3.build(path)
    diff = DIFF.build(path)
    events = EVENTS.build(path)
    failures = (
        [f"source: {x}" for x in COMPLETE.validate(source)]
        + [f"literal word: {x}" for x in WORD.validate(literal)]
        + [f"P3: {x}" for x in P3.validate(p3)]
        + [f"differential metric: {x}" for x in DIFF.validate(diff)]
        + [f"differential event: {x}" for x in EVENTS.validate(events)]
    )
    if failures:
        raise RuntimeError(f"differential word prerequisites failed: {failures}")
    event_closed = bool(events["source_uniform_finite_angle_event_Jacobians_closed"])
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "P3_frozen_not_modified": True,
        "P3_delta_consumed": 1.0e-18,
        "same_complete_BRMM_source_token_required_for_every_event": True,
        "literal_shipping_event_order_required": True,
        "prediction_required_every_valid_IMU_sample": bool(literal["every_valid_imu_sample_requires_prediction"]),
        "accelerometer_required_every_valid_IMU_sample": bool(literal["every_valid_imu_sample_requires_accelerometer_Joseph"]),
        "all_due_S_updates_required": bool(literal["S_scheduler_is_executed_not_replaced_by_selected_four"]),
        "same_P_H_R_cell_required_for_every_Joseph_event": True,
        "independent_K_input_allowed_for_theorem": False,
        "actual_applied_per_axis_RS_required_on_every_S_event": True,
        "actual_RS_provenance_token": EVENTS.ACTUAL_RS_PROVENANCE,
        "asynchronous_vector_events_retained": bool(literal["magnetometer_is_asynchronous_external_event_family"]),
        "aw_covariance_floor_state_jacobian_is_identity": True,
        "A21_bias_projection_generalized_Jacobian_available": bool(
            events["A21_bias_projection_generalized_Jacobian_available"]
        ),
        "A21_bias_projection_source_uniform_attachment_closed": False,
        "H_to_A_unique_rectangular_event_required": True,
        "H_to_A_homogeneous_lift_constructor_available": True,
        "H_to_A_homogeneous_lift": "[I18;0]",
        "H_to_A_held_ba_error_retained_as_separate_forcing": True,
        "H_to_A_covariance_floor_retained_as_separate_metric_event": True,
        "full_state_dimensions": {"H18": 18, "A21": 21},
        "differential_event_qualification": events["qualification"],
        "cocycle_definition": "J_W=J_{N-1}*...*J_0",
        "full_matrix_interval_LDLT_is_terminal_gate": True,
        "packetwise_norm_sum_used": False,
        "packet_count_multiplier_used": False,
        "selected_S_subset_used_as_word": False,
        "independent_RS_schedule_used": False,
        "state_elimination_used": False,
        "trajectory_replay_used": False,
        "source_family_replaced": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "source_uniform_event_Jacobian_enclosures_supplied": event_closed,
        "source_uniform_complete_word_Jacobian_enclosed": False,
        "P4_promoted_here": False,
        "next_obligation": (
            "close same-history source cells for prediction/Joseph/projection events, carry the explicit H->A release forcing/covariance event separately, then compose the literal cocycle; do not use independent event boxes"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_BRMM_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    if float(d.get("P3_delta_consumed", 0.0)) != 1.0e-18:
        f.append("frozen P3 delta changed")
    for key in (
        "P3_frozen_not_modified", "same_complete_BRMM_source_token_required_for_every_event",
        "literal_shipping_event_order_required", "prediction_required_every_valid_IMU_sample",
        "accelerometer_required_every_valid_IMU_sample", "all_due_S_updates_required",
        "same_P_H_R_cell_required_for_every_Joseph_event",
        "actual_applied_per_axis_RS_required_on_every_S_event", "asynchronous_vector_events_retained",
        "aw_covariance_floor_state_jacobian_is_identity",
        "A21_bias_projection_generalized_Jacobian_available",
        "H_to_A_unique_rectangular_event_required", "H_to_A_homogeneous_lift_constructor_available",
        "H_to_A_held_ba_error_retained_as_separate_forcing",
        "H_to_A_covariance_floor_retained_as_separate_metric_event",
        "full_matrix_interval_LDLT_is_terminal_gate",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "independent_K_input_allowed_for_theorem", "A21_bias_projection_source_uniform_attachment_closed",
        "packetwise_norm_sum_used", "packet_count_multiplier_used", "selected_S_subset_used_as_word",
        "independent_RS_schedule_used", "state_elimination_used", "trajectory_replay_used",
        "source_family_replaced", "filter_changed", "declared_domain_changed",
        "source_uniform_event_Jacobian_enclosures_supplied",
        "source_uniform_complete_word_Jacobian_enclosed", "P4_promoted_here",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("actual_RS_provenance_token") != EVENTS.ACTUAL_RS_PROVENANCE:
        f.append("actual R_S provenance token changed")
    if d.get("H_to_A_homogeneous_lift") != "[I18;0]":
        f.append("H->A homogeneous lift changed")
    if d.get("full_state_dimensions") != {"H18": 18, "A21": 21}:
        f.append("full-state dimensions changed")
    return list(dict.fromkeys(f))


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
    print(json.dumps({
        "source": d["canonical_source"],
        "literal_order": d["literal_shipping_event_order_required"],
        "same_cell_Joseph": d["same_P_H_R_cell_required_for_every_Joseph_event"],
        "actual_RS_each_S": d["actual_applied_per_axis_RS_required_on_every_S_event"],
        "projection_map": d["A21_bias_projection_generalized_Jacobian_available"],
        "H_to_A_lift": d["H_to_A_homogeneous_lift"],
        "event_Jacobians_closed": d["source_uniform_event_Jacobian_enclosures_supplied"],
        "P4_promoted_here": d["P4_promoted_here"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
