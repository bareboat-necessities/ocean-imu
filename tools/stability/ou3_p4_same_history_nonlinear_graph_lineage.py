#!/usr/bin/env python3
"""Same-history nonlinear H18/A21 graph lineage for complete-SEA3 P4.

This module is the theorem-facing bridge between the branch-correlated prefix
selectors and the retained exact nonlinear prediction/Joseph/reset/projection
maps.  It deliberately does not generate a SEA3 source family.  A caller must
supply one already-connected selector lineage and, for A21, one absolute
accelerometer-bias history attached to that same lineage.

No event may provide an independent P, H, R, K, tuner schedule, R_S, body rate,
or geometry box.  Prediction coordinates come from the selector's provider
sample and committed schedule.  Joseph events consume the passive event-local
P-before/H/R cells captured inside the same trusted shipping transition.  Every
due S=0 retains ACTUAL_APPLIED_SPECTRALMSE_RS provenance.  A21 projection is
composed after every Joseph/reset event exactly as shipping does and therefore
requires the same-source absolute true residual-bias coordinate needed to
reconstruct the estimated bias presented to the deployed 0.4 m/s^2 ball
projection.

The helper point smoke uses the identically-zero physical bias history.  That is
a single homogeneous finite-tau_b point history used only to falsify/validate
this composition layer; it is not a source-uniform bias box and cannot promote
P4.  Source-uniform physical-bias materialization remains an upstream Phase-3
obligation.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Mapping, Sequence

from ou3_interval import Interval, matrix_identity
import ou3_p4_complete_sea3_differential_events as EVENTS
import ou3_p4_complete_sea3_differential_prediction as PREDICTION
import ou3_p4_complete_sea3_differential_word as DWORD
import ou3_p4_complete_sea3_same_history_prefix_selectors as SELECTORS
import ou3_sea3_complete_window_execution_kernel as KERNEL
import ou3_sea3_frontend_state_step as FRONTEND
import ou3_mems_bias_contract as BIAS

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_P4_COMPLETE_SEA3_SAME_HISTORY_NONLINEAR_GRAPH_LINEAGE_V1"
CANONICAL_SOURCE = "COMPLETE_SEA3_NORMAL_LIVE_WORD"


@dataclass(frozen=True)
class SameHistoryBiasLineage:
    """Absolute A21 physical bias cells keyed by the exact selector child id."""

    endpoint_source_cell_id: str
    bias_true_by_source_cell_id: Mapping[
        str, tuple[Interval, Interval, Interval]
    ]
    source_uniform_materialization: bool = False
    homogeneous_root: tuple[Interval, Interval, Interval] | None = None
    root_source_cell_id: str | None = None
    tau_s: Interval | None = None

    def validate_homogeneous(self, lineage, constants, projection_limit):
        """Check BIAS1 ancestry and matched dynamics, not just matching ID sets."""
        if self.source_uniform_materialization:
            raise RuntimeError("a bias point/hull lineage is not a qualified source-uniform family")
        if self.homogeneous_root is None or self.tau_s is None:
            raise RuntimeError("BIAS1 requires one true-bias root and common GM parameter")
        if self.root_source_cell_id != lineage[0].parent_source_cell_id:
            raise RuntimeError("BIAS1 true-bias root detached from selector ancestry")
        if self.tau_s != constants.accel_bias_tau_s:
            raise RuntimeError("true/filter tau mismatch requires explicit forcing, not the homogeneous map")
        norm_sq = Interval.point(0.0)
        for x in self.homogeneous_root:
            norm_sq = norm_sq + x * x
        radius_sq = Interval.point(projection_limit) * Interval.point(projection_limit)
        if norm_sq.hi > radius_sq.lo:
            raise RuntimeError("homogeneous true-bias root does not certify projection zero-error invariance")
        for selector in lineage:
            elapsed = constants.h * Interval.point(float(selector.prefix_length))
            expected = BIAS.homogeneous_bias_at(self.homogeneous_root, self.tau_s, elapsed)
            if self.at(selector.source_cell_id) != expected:
                raise RuntimeError("BIAS1 bias prefix detached from the common GM root")

    def at(self, source_cell_id: str) -> tuple[Interval, Interval, Interval]:
        try:
            value = self.bias_true_by_source_cell_id[source_cell_id]
        except KeyError as exc:
            raise RuntimeError(
                f"A21 same-history absolute bias missing at {source_cell_id}"
            ) from exc
        if len(value) != 3 or any(not isinstance(x, Interval) for x in value):
            raise TypeError("A21 same-history bias cell must be three outward intervals")
        return value


@dataclass(frozen=True)
class NonlinearEventRecord:
    source_cell_id: str
    mode: str
    kind: str
    event_index_in_sample: int
    same_P_H_R_cell: bool
    actual_rs_provenance: bool
    projection_branch: str


@dataclass(frozen=True)
class NonlinearLineageResult:
    endpoint_source_cell_id: str
    source_token: str
    mode: str
    state_out: tuple[Interval, ...]
    J_word: tuple[tuple[Interval, ...], ...]
    event_records: tuple[NonlinearEventRecord, ...]
    prefixes: int
    predictions: int
    floors: int
    S_updates: int
    accelerometer_updates: int
    vector_updates: int


def _state_point(n: int, value: float = 0.0) -> list[Interval]:
    return [Interval.point(float(value)) for _ in range(n)]


def _projection_limit(domain_path: Path) -> float:
    domain = json.loads(Path(domain_path).resolve().read_text(encoding="utf-8"))
    limit = float(domain["normal_live"]["active_accelerometer_bias_projection_limit_mps2"])
    if limit != 0.4:
        raise RuntimeError("Phase-2 graph requires deployed 0.4 m/s^2 A21 projection limit")
    return limit


def _word_event(kind: str, source_token: str, J, *, source_event=None):
    if source_event is not None:
        return DWORD.event_from_source_joseph(kind, source_token, source_event)
    return DWORD.DifferentialEvent(kind=kind, J=J, source_token=source_token)


def _measurement_geometry(selector: SELECTORS.PrefixSelector, cell: KERNEL.RiccatiEventCell):
    sample = selector.sample_coordinates
    if cell.kind == "S_zero":
        return {}
    if cell.kind == "accelerometer":
        return {"f_hat": sample.f_cog_body, "R_hat": sample.R_wb}
    if cell.kind == "magnetometer":
        i = cell.magnetic_event_index
        if i is None or i < 0 or i >= len(sample.magnetometer_events_after_imu):
            raise RuntimeError("magnetometer event cell lost same-history source index")
        return {"m_body": sample.magnetometer_events_after_imu[i].m_body}
    raise ValueError(f"unsupported measurement event {cell.kind}")


def _consume_mode_lineage(
    *,
    mode: str,
    lineage: Sequence[SELECTORS.PrefixSelector],
    initial_state: Sequence[Interval],
    source_token: str,
    constants: KERNEL.KernelConstants,
    bias_lineage: SameHistoryBiasLineage | None,
    projection_limit: float,
) -> NonlinearLineageResult:
    n = 18 if mode == "H" else 21
    if len(initial_state) != n or any(not isinstance(x, Interval) for x in initial_state):
        raise ValueError(f"{mode} initial physical error cell must have dimension {n}")
    if mode == "A" and bias_lineage is None:
        raise RuntimeError("A21 nonlinear lineage requires same-history absolute bias")
    if mode == "H" and bias_lineage is not None:
        raise ValueError("H18 lineage must not consume an A21 absolute-bias history")
    if mode == "A":
        bias_lineage.validate_homogeneous(lineage, constants, projection_limit)

    word = DWORD.initialize(mode, source_token)
    state = list(initial_state)
    records: list[NonlinearEventRecord] = []

    for expected_prefix, selector in enumerate(lineage, start=1):
        if selector.prefix_length != expected_prefix:
            raise RuntimeError("nonlinear graph lineage skipped or duplicated a selector prefix")
        cells = selector.H_event_cells if mode == "H" else selector.A_event_cells
        events = selector.H_events_this_sample if mode == "H" else selector.A_events_this_sample
        if tuple(cell.kind for cell in cells) != tuple(events):
            raise RuntimeError("event-local cells detached from selector shipping slice")
        if not cells or cells[0].kind != "prediction":
            raise RuntimeError("selector sample does not begin with shipping prediction")

        bias_true = bias_lineage.at(selector.source_cell_id) if mode == "A" else None

        for cell in cells:
            if cell.mode != mode:
                raise RuntimeError("event-local Riccati cell mode detached from lineage")
            if cell.kind == "prediction":
                pred = PREDICTION.prediction_event(
                    mode,
                    state,
                    selector.sample_coordinates.omega_body_corrected,
                    constants.h,
                    selector.active_schedule.tau,
                    tau_ba=constants.accel_bias_tau_s if mode == "A" else None,
                )
                event = _word_event("prediction", source_token, pred["J_state"])
                DWORD.apply_event(word, event)
                state = list(pred["state_out"])
                records.append(
                    NonlinearEventRecord(
                        selector.source_cell_id,
                        mode,
                        "prediction",
                        cell.event_index_in_sample,
                        False,
                        False,
                        "not_applicable",
                    )
                )
                continue

            if cell.kind == "aw_floor":
                event = _word_event("aw_floor", source_token, matrix_identity(n))
                DWORD.apply_event(word, event)
                records.append(
                    NonlinearEventRecord(
                        selector.source_cell_id,
                        mode,
                        "aw_floor",
                        cell.event_index_in_sample,
                        False,
                        False,
                        "not_applicable",
                    )
                )
                continue

            if cell.kind not in ("S_zero", "accelerometer", "magnetometer"):
                raise RuntimeError(f"unsupported shipping event in nonlinear lineage: {cell.kind}")
            if cell.H is None or cell.R is None or not cell.P_before or not cell.P_after:
                raise RuntimeError(f"{cell.kind} lost same-event P/H/R cell")

            kwargs = _measurement_geometry(selector, cell)
            provenance = EVENTS.ACTUAL_RS_PROVENANCE if cell.kind == "S_zero" else None
            if cell.kind == "S_zero" and cell.actual_rs_from_committed_schedule is not True:
                raise RuntimeError("due S event lost actual committed R_S provenance")
            source_event = EVENTS.source_joseph_event(
                mode,
                state,
                cell.P_before,
                cell.R,
                cell.kind,
                R_provenance=provenance,
                bias_true=bias_true,
                bias_projection_limit=projection_limit if mode == "A" else None,
                **kwargs,
            )
            if source_event["H"] != cell.H:
                raise RuntimeError(f"{cell.kind} nonlinear H detached from captured shipping H")
            event = _word_event(
                cell.kind,
                source_token,
                source_event["J_state"],
                source_event=source_event,
            )
            DWORD.apply_event(word, event)
            state = list(source_event["state_out"])
            records.append(
                NonlinearEventRecord(
                    selector.source_cell_id,
                    mode,
                    cell.kind,
                    cell.event_index_in_sample,
                    True,
                    provenance == EVENTS.ACTUAL_RS_PROVENANCE,
                    str(source_event["bias_projection_branch"]),
                )
            )

    if word.events != [record.kind for record in records]:
        raise RuntimeError("nonlinear differential word detached from physical event ledger")
    return NonlinearLineageResult(
        endpoint_source_cell_id=lineage[-1].source_cell_id,
        source_token=source_token,
        mode=mode,
        state_out=tuple(state),
        J_word=tuple(tuple(x for x in row) for row in word.J_word),
        event_records=tuple(records),
        prefixes=len(lineage),
        predictions=word.predictions,
        floors=word.floors,
        S_updates=word.S_updates,
        accelerometer_updates=word.accelerometer_updates,
        vector_updates=word.vector_updates,
    )


def consume_endpoint_lineage(
    selectors: Sequence[SELECTORS.PrefixSelector],
    endpoint_source_cell_id: str,
    *,
    initial_H_state: Sequence[Interval],
    initial_A_state: Sequence[Interval],
    A21_bias_lineage: SameHistoryBiasLineage,
    domain_path: Path = DEFAULT_DOMAIN,
) -> tuple[NonlinearLineageResult, NonlinearLineageResult]:
    """Compose exact nonlinear H18/A21 event maps on one selector lineage."""
    lineage = SELECTORS.lineage_for_endpoint(selectors, endpoint_source_cell_id)
    if not lineage:
        raise ValueError("endpoint nonlinear lineage is empty")
    if A21_bias_lineage.endpoint_source_cell_id != endpoint_source_cell_id:
        raise RuntimeError("A21 absolute-bias history belongs to a different endpoint lineage")
    lineage_ids = {selector.source_cell_id for selector in lineage}
    bias_ids = set(A21_bias_lineage.bias_true_by_source_cell_id)
    if bias_ids != lineage_ids:
        raise RuntimeError("A21 absolute-bias history must match every and only selector prefix")

    source_token = f"{CANONICAL_SOURCE}:{endpoint_source_cell_id}"
    constants = KERNEL._process_constants(domain_path)
    limit = _projection_limit(domain_path)
    H = _consume_mode_lineage(
        mode="H",
        lineage=lineage,
        initial_state=initial_H_state,
        source_token=source_token,
        constants=constants,
        bias_lineage=None,
        projection_limit=limit,
    )
    A = _consume_mode_lineage(
        mode="A",
        lineage=lineage,
        initial_state=initial_A_state,
        source_token=source_token,
        constants=constants,
        bias_lineage=A21_bias_lineage,
        projection_limit=limit,
    )
    return H, A


def homogeneous_bias_lineage(
    lineage: Sequence[SELECTORS.PrefixSelector],
    root: tuple[Interval, Interval, Interval],
    domain_path: Path = DEFAULT_DOMAIN,
) -> SameHistoryBiasLineage:
    """Retain one root and matched tau alongside their derived prefix hulls.

    This is the zero-forcing BIAS1 model only, conditional on BIAS0. The
    interval AD consumer is still an enclosure, not an exhaustive SEA3 cover.
    """
    if not lineage:
        raise ValueError("nonempty selector lineage required")
    constants = KERNEL._process_constants(domain_path)
    return SameHistoryBiasLineage(
        endpoint_source_cell_id=lineage[-1].source_cell_id,
        bias_true_by_source_cell_id={
            selector.source_cell_id: BIAS.homogeneous_bias_at(
                root, constants.accel_bias_tau_s,
                constants.h * Interval.point(float(selector.prefix_length)),
            ) for selector in lineage
        },
        source_uniform_materialization=False,
        homogeneous_root=root,
        root_source_cell_id=lineage[0].parent_source_cell_id,
        tau_s=constants.accel_bias_tau_s,
    )


def _zero_bias_lineage(lineage, domain_path: Path = DEFAULT_DOMAIN) -> SameHistoryBiasLineage:
    z = (Interval.point(0.0), Interval.point(0.0), Interval.point(0.0))
    return homogeneous_bias_lineage(lineage, z, domain_path)


def _smoke(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    frontend = FRONTEND._point_state()
    P0_H, P0_A, seed_meta = SELECTORS._live_structured_point_covariance_fixture(
        frontend, domain_path
    )
    sample = SELECTORS._point_sample()
    endpoints, selectors, meta = SELECTORS.execute_with_prefix_selectors(
        frontend_entry=frontend,
        P0_H=P0_H,
        P0_A=P0_A,
        samples=[sample, sample],
        domain_path=domain_path,
        branch_limit=128,
    )
    endpoint = endpoints[0].source_cell_id
    lineage = SELECTORS.lineage_for_endpoint(selectors, endpoint)
    bias = _zero_bias_lineage(lineage, domain_path)
    H, A = consume_endpoint_lineage(
        selectors,
        endpoint,
        initial_H_state=_state_point(18),
        initial_A_state=_state_point(21),
        A21_bias_lineage=bias,
        domain_path=domain_path,
    )
    expected_events = sum(len(selector.H_event_cells) for selector in lineage)
    A_measurements = [
        record for record in A.event_records
        if record.kind in ("S_zero", "accelerometer", "magnetometer")
    ]
    S_records = [record for record in H.event_records + A.event_records if record.kind == "S_zero"]
    return {
        "selector_seed": seed_meta,
        "selector_meta": meta,
        "endpoint_source_cell_id": endpoint,
        "prefixes": len(lineage),
        "H_events": len(H.event_records),
        "A_events": len(A.event_records),
        "expected_events_per_mode": expected_events,
        "H_full_word_dimension": [len(H.J_word), len(H.J_word[0])],
        "A_full_word_dimension": [len(A.J_word), len(A.J_word[0])],
        "same_source_token_H_A": H.source_token == A.source_token,
        "all_measurements_use_same_P_H_R": all(
            record.same_P_H_R_cell for record in H.event_records + A.event_records
            if record.kind in ("S_zero", "accelerometer", "magnetometer")
        ),
        "all_due_S_retain_actual_RS": bool(S_records)
        and all(record.actual_rs_provenance for record in S_records),
        "A21_projection_attached_every_measurement": bool(A_measurements)
        and all(record.projection_branch != "not_applicable" for record in A_measurements),
        "zero_bias_point_projection_inactive": bool(A_measurements)
        and all(record.projection_branch == "inactive" for record in A_measurements),
        "point_bias_history_source_uniform": bias.source_uniform_materialization,
        "point_history_only_not_P4": True,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    selectors = SELECTORS.build(domain_path)
    sf = SELECTORS.validate(selectors)
    events = EVENTS.build(domain_path)
    ef = EVENTS.validate(events)
    prediction = PREDICTION.build(domain_path)
    pf = PREDICTION.validate(prediction)
    dword = DWORD.build(domain_path)
    wf = DWORD.validate(dword)
    bad = {"selectors": sf, "events": ef, "prediction": pf, "word": wf}
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"same-history nonlinear lineage prerequisites failed: {bad}")
    smoke = _smoke(domain_path)
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": CANONICAL_SOURCE,
        "mems_bias_preconditions": BIAS.build(domain_path),
        "P3_delta_preserved": 1.0e-18,
        "branch_correlated_prefix_selector_consumed": True,
        "event_local_same_P_H_R_cells_consumed": True,
        "prediction_uses_same_selector_omega_tau_h": True,
        "actual_applied_RS_provenance_retained": True,
        "full_H18_nonlinear_lineage_composition_available": True,
        "full_A21_nonlinear_lineage_composition_available": True,
        "A21_absolute_bias_history_required": True,
        "A21_projection_generalized_Jacobian_composed_after_every_Joseph": True,
        "independent_true_bias_event_boxes_allowed": False,
        "point_zero_bias_history_is_homogeneous_and_nonpromoting": True,
        "source_uniform_absolute_bias_history_materialized_here": False,
        "source_uniform_SEA3_window_family_materialized_here": False,
        "joint_graph_sectors_assembled_here": False,
        "endpoint_augmented_LDLT_closed_here": False,
        "every_prefix_augmented_LDLT_closed_here": False,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "quality_gates_changed": False,
        "declared_domain_changed": False,
        "P4_promoted_here": False,
        "smoke": smoke,
        "next_obligation": (
            "materialize the correlated finite-window SEA3 source cells, including the same-history absolute A21 physical-bias coordinate, then lift this exact nonlinear lineage into dense joint graph sectors before attempting the canonical-point augmented master"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    f.extend(f"MEMS bias: {x}" for x in BIAS.validate(d.get("mems_bias_preconditions", {})))
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != CANONICAL_SOURCE:
        f.append("canonical source changed")
    if float(d.get("P3_delta_preserved", 0.0)) != 1.0e-18:
        f.append("frozen P3 delta changed")
    for key in (
        "branch_correlated_prefix_selector_consumed",
        "event_local_same_P_H_R_cells_consumed",
        "prediction_uses_same_selector_omega_tau_h",
        "actual_applied_RS_provenance_retained",
        "full_H18_nonlinear_lineage_composition_available",
        "full_A21_nonlinear_lineage_composition_available",
        "A21_absolute_bias_history_required",
        "A21_projection_generalized_Jacobian_composed_after_every_Joseph",
        "point_zero_bias_history_is_homogeneous_and_nonpromoting",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "independent_true_bias_event_boxes_allowed",
        "source_uniform_absolute_bias_history_materialized_here",
        "source_uniform_SEA3_window_family_materialized_here",
        "joint_graph_sectors_assembled_here",
        "endpoint_augmented_LDLT_closed_here",
        "every_prefix_augmented_LDLT_closed_here",
        "trajectory_replay_used",
        "filter_changed",
        "quality_gates_changed",
        "declared_domain_changed",
        "P4_promoted_here",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    smoke = d.get("smoke", {})
    for key in (
        "same_source_token_H_A",
        "all_measurements_use_same_P_H_R",
        "all_due_S_retain_actual_RS",
        "A21_projection_attached_every_measurement",
        "zero_bias_point_projection_inactive",
        "point_history_only_not_P4",
    ):
        if smoke.get(key) is not True:
            f.append(f"smoke lost {key}")
    if smoke.get("point_bias_history_source_uniform") is not False:
        f.append("point bias smoke incorrectly claimed source uniformity")
    if smoke.get("H_full_word_dimension") != [18, 18]:
        f.append("H18 nonlinear word dimension changed")
    if smoke.get("A_full_word_dimension") != [21, 21]:
        f.append("A21 nonlinear word dimension changed")
    if int(smoke.get("prefixes", 0)) != 2:
        f.append("nonlinear smoke did not retain both prefixes")
    if int(smoke.get("H_events", 0)) != int(smoke.get("expected_events_per_mode", -1)):
        f.append("H18 nonlinear event count changed")
    if int(smoke.get("A_events", 0)) != int(smoke.get("expected_events_per_mode", -1)):
        f.append("A21 nonlinear event count changed")
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
        "nonlinear_H18": d["full_H18_nonlinear_lineage_composition_available"],
        "nonlinear_A21": d["full_A21_nonlinear_lineage_composition_available"],
        "A21_absolute_bias_required": d["A21_absolute_bias_history_required"],
        "source_uniform": d["source_uniform_SEA3_window_family_materialized_here"],
        "P4_promoted": d["P4_promoted_here"],
        "smoke": d["smoke"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
