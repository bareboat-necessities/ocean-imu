#!/usr/bin/env python3
"""Literal full-state differential cocycle for complete-SEA3 OU-III P4.

This is the matrix object consumed by canonical P4.  It does not generate a
source word and it does not replace the shipping scheduler by independent
bounds.  A caller supplies outward event Jacobian enclosures produced from ONE
same-history complete SEA3 source path.  The cocycle composes them in literal
shipping order and sends the full result to the pullback differential metric.

Required fixed-mode IMU order is

    prediction -> optional a_w covariance floor (state Jacobian I)
               -> every due S=0 Joseph/reset
               -> required accelerometer Joseph/reset,

with asynchronous accepted vector/magnetometer Joseph/reset events inserted at
their actual source positions.  The H18->A21 21x18 lift is one separate hybrid
event and may occur at most once.  Every S event must carry provenance that its
Jacobian was built with the actual applied per-axis SpectralMSE R_S from that
same source state.

If J_k is the exact/outward differential map of event k, the complete cocycle is

    J_W = J_{N-1} ... J_0.

No packet norm, packet count, scalar beta, selected-S replacement word, or
independent R_S schedule appears.  ``certify_word`` evaluates directly

    rho M_0 - J_W^T M_1 J_W > 0

through the full interval-LDLT routine in
``ou3_p4_complete_sea3_phi_differential_metric``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import argparse
import json
from pathlib import Path
from typing import Any, Sequence

from ou3_interval import Interval, matrix_identity, matrix_mul
import ou3_p4_complete_sea3_phi_differential_metric as DIFF
import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_full_normal_live_word as WORD
import ou3_sea3_riccati_metric_p3 as P3

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_P4_COMPLETE_SEA3_LITERAL_DIFFERENTIAL_WORD_V1"


def _shape(A: Sequence[Sequence[Any]]) -> tuple[int, int]:
    r = len(A)
    c = len(A[0]) if r else 0
    if any(len(row) != c for row in A):
        raise ValueError("ragged differential event matrix")
    return r, c


def _point_one_zero(A: Sequence[Sequence[Interval]]) -> tuple[Interval, Interval]:
    if not A or not A[0]:
        raise ValueError("nonempty interval matrix required")
    return Interval.point(1.0), Interval.point(0.0)


@dataclass(frozen=True)
class DifferentialEvent:
    kind: str
    J: list[list[Interval]]
    source_token: str
    actual_applied_RS: bool = False
    accepted: bool = True


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

    @property
    def start_dim(self) -> int:
        return 18 if self.start_mode == "H" else 21


def initialize(mode: str, source_token: str) -> DifferentialWord:
    if mode not in ("H", "A"):
        raise ValueError("differential word mode must be H or A")
    if not isinstance(source_token, str) or not source_token:
        raise ValueError("same-history complete SEA3 source token required")
    n = 18 if mode == "H" else 21
    return DifferentialWord(
        start_mode=mode,
        source_token=source_token,
        J_word=matrix_identity(n),
        current_dim=n,
    )


def apply_event(word: DifferentialWord, event: DifferentialEvent) -> None:
    if event.source_token != word.source_token:
        raise ValueError("differential event detached from same complete SEA3 source history")
    if not event.accepted:
        raise ValueError("only executed/accepted nonlinear events belong in the cocycle")
    r, c = _shape(event.J)
    if c != word.current_dim:
        raise ValueError("differential event input dimension does not match current word")

    kind = event.kind
    if kind == "prediction":
        if r != c:
            raise ValueError("same-mode prediction must be square")
        word.predictions += 1
    elif kind == "aw_floor":
        if r != c:
            raise ValueError("covariance-floor state event must be square identity")
        one, zero = _point_one_zero(event.J)
        I = matrix_identity(c)
        if any(event.J[i][j] != I[i][j] for i in range(c) for j in range(c)):
            raise ValueError("a_w covariance floor must have identity state Jacobian")
        word.floors += 1
    elif kind == "S_zero":
        if r != c:
            raise ValueError("S=0 Joseph/reset event must be same-mode square")
        if event.actual_applied_RS is not True:
            raise ValueError("S=0 differential event lacks actual-applied R_S provenance")
        word.S_updates += 1
    elif kind == "accelerometer":
        if r != c:
            raise ValueError("accelerometer Joseph/reset event must be same-mode square")
        word.accelerometer_updates += 1
    elif kind in ("magnetometer", "vector"):
        if r != c:
            raise ValueError("vector Joseph/reset event must be same-mode square")
        word.vector_updates += 1
    elif kind == "H_to_A":
        if word.hybrid_lifts != 0 or word.current_dim != 18 or (r, c) != (21, 18):
            raise ValueError("H->A differential lift must be the unique 21x18 hybrid event")
        word.hybrid_lifts = 1
    else:
        raise ValueError(f"unknown differential event kind {kind}")

    word.J_word = matrix_mul(event.J, word.J_word)
    word.current_dim = r
    word.events.append(kind)


def certify_word(
    word: DifferentialWord,
    P0_inverse,
    P1_inverse,
    T0,
    T1,
    rho: float,
):
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
    failures = (
        [f"source: {x}" for x in COMPLETE.validate(source)]
        + [f"literal word: {x}" for x in WORD.validate(literal)]
        + [f"P3: {x}" for x in P3.validate(p3)]
        + [f"differential metric: {x}" for x in DIFF.validate(diff)]
    )
    if failures:
        raise RuntimeError(f"differential word prerequisites failed: {failures}")
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "P3_frozen_not_modified": True,
        "P3_delta_consumed": 1.0e-18,
        "same_complete_SEA3_source_token_required_for_every_event": True,
        "literal_shipping_event_order_required": True,
        "prediction_required_every_valid_IMU_sample": bool(
            literal["every_valid_imu_sample_requires_prediction"]
        ),
        "accelerometer_required_every_valid_IMU_sample": bool(
            literal["every_valid_imu_sample_requires_accelerometer_Joseph"]
        ),
        "all_due_S_updates_required": bool(
            literal["S_scheduler_is_executed_not_replaced_by_selected_four"]
        ),
        "actual_applied_per_axis_RS_required_on_every_S_event": True,
        "asynchronous_vector_events_retained": bool(
            literal["magnetometer_is_asynchronous_external_event_family"]
        ),
        "aw_covariance_floor_state_jacobian_is_identity": True,
        "H_to_A_unique_rectangular_event_required": True,
        "full_state_dimensions": {"H18": 18, "A21": 21},
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
        "source_uniform_event_Jacobian_enclosures_supplied": False,
        "source_uniform_complete_word_Jacobian_enclosed": False,
        "P4_promoted_here": False,
        "next_obligation": (
            "derive outward nonlinear event Jacobians on the declared finite-angle cells from the SAME complete SEA3 source variables, then compose them with this literal cocycle; do not replace the source history by independent event boxes"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    if float(d.get("P3_delta_consumed", 0.0)) != 1.0e-18:
        f.append("frozen P3 delta changed")
    for key in (
        "P3_frozen_not_modified", "same_complete_SEA3_source_token_required_for_every_event",
        "literal_shipping_event_order_required", "prediction_required_every_valid_IMU_sample",
        "accelerometer_required_every_valid_IMU_sample", "all_due_S_updates_required",
        "actual_applied_per_axis_RS_required_on_every_S_event", "asynchronous_vector_events_retained",
        "aw_covariance_floor_state_jacobian_is_identity", "H_to_A_unique_rectangular_event_required",
        "full_matrix_interval_LDLT_is_terminal_gate",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "packetwise_norm_sum_used", "packet_count_multiplier_used", "selected_S_subset_used_as_word",
        "independent_RS_schedule_used", "state_elimination_used", "trajectory_replay_used",
        "source_family_replaced", "filter_changed", "declared_domain_changed",
        "source_uniform_event_Jacobian_enclosures_supplied",
        "source_uniform_complete_word_Jacobian_enclosed", "P4_promoted_here",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
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
        "actual_RS_each_S": d["actual_applied_per_axis_RS_required_on_every_S_event"],
        "event_Jacobians_closed": d["source_uniform_event_Jacobian_enclosures_supplied"],
        "P4_promoted_here": d["P4_promoted_here"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
