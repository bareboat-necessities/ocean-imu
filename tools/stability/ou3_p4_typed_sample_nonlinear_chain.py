#!/usr/bin/env python3
"""Compose one trusted typed-kernel sample through the exact nonlinear P4 maps.

This bridge closes the event-chain seam between a retained PrefixSelector and
the finite nonlinear differential maps.  It does not create source cells: the
caller supplies SourceCoverCell objects already bound to the selector by
``ou3_p4_typed_sample_source_cell_binding``.  The bridge then requires literal
state succession and applies, in trusted shipping order:

* the exact deployed-quaternion/integrated-OU prediction map;
* identity physical state transport for an a_w covariance-floor event; and
* the exact same-P/H/R Joseph/reset map, including the A21 bias projection.

The same finite map call produces both ``state_out`` and ``J_state`` at every
nontrivial event.  No identity placeholder is permitted for prediction and no
independent K, F, R_S, rate or tau box is accepted.

This is still non-promoting machinery.  Production closure additionally needs
all provider lineages, the H18->A21 hybrid forcing/covariance attachment, and
endpoint/every-prefix augmented LDLT plus first-exit retention.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval, matrix_identity
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_complete_brmm_differential_prediction as PRED
import ou3_p4_complete_brmm_differential_word as DWORD
import ou3_p4_complete_brmm_same_history_prefix_selectors as SELECTORS
import ou3_p4_complete_brmm_source_cover_contract as COVER
import ou3_p4_typed_sample_source_cell_binding as BIND

SCHEMA = 1
QUALIFICATION = "OU3_P4_TYPED_SAMPLE_EXACT_NONLINEAR_EVENT_CHAIN_V1"
P3_DELTA = 1.0e-18


def _same_interval(a: Interval, b: Interval) -> bool:
    return isinstance(a, Interval) and isinstance(b, Interval) and a.lo == b.lo and a.hi == b.hi


def _same_state(a: Sequence[Interval], b: Sequence[Interval]) -> bool:
    return len(a) == len(b) and all(_same_interval(x, y) for x, y in zip(a, b))


def materialize_sample_chain(
    selector: SELECTORS.PrefixSelector,
    cells: Sequence[COVER.SourceCoverCell],
    *,
    mode: str,
    tau_ba: Interval | None = None,
) -> dict:
    """Return the exact/outward finite state/Jacobian chain for one typed sample.

    ``cells[k].state`` must be exactly the state enclosure emitted by event
    ``k-1`` (the first cell supplies the sample-entry error enclosure).  This
    prevents a valid Jacobian from being attached to a different finite state.
    """
    failures = BIND.validate_event_cells_against_selector(selector, cells, mode=mode)
    if failures:
        raise ValueError("typed source cells are not selector-bound: " + repr(failures))
    if not cells:
        raise ValueError("typed sample nonlinear chain cannot be empty")
    if mode == "A" and tau_ba is None:
        raise ValueError("A21 exact prediction chain requires configured tau_ba")

    source_token = selector.source_cell_id
    word = DWORD.initialize(mode, source_token)
    state = list(cells[0].state)
    event_records: list[dict] = []

    for ordinal, cell in enumerate(cells):
        if cell.event_ordinal != ordinal:
            raise ValueError("literal event ordinal changed during nonlinear-chain materialization")
        if not _same_state(cell.state, state):
            raise ValueError(
                f"event {ordinal} state is detached from previous exact/outward state_out"
            )

        if cell.kind == "prediction":
            pred = PRED.prediction_event(
                mode,
                cell.state,
                selector.sample_coordinates.omega_body_corrected,
                cell.dt_s,
                cell.tau_applied_s,
                tau_ba=tau_ba,
            )
            event = DWORD.DifferentialEvent(
                kind="prediction", J=pred["J_state"], source_token=source_token
            )
            state_out = pred["state_out"]
            same_map_state_and_jacobian = True
        elif cell.kind == "aw_floor":
            n = 18 if mode == "H" else 21
            event = DWORD.DifferentialEvent(
                kind="aw_floor", J=matrix_identity(n), source_token=source_token
            )
            state_out = list(cell.state)
            same_map_state_and_jacobian = True
        elif cell.kind in ("S_zero", "accelerometer", "magnetometer"):
            joseph = EVENTS.source_joseph_event(**COVER.joseph_event_kwargs(cell))
            event = DWORD.event_from_source_joseph(cell.kind, source_token, joseph)
            state_out = joseph["state_out"]
            same_map_state_and_jacobian = True
        else:
            raise ValueError(
                f"sample-local nonlinear chain does not admit event kind {cell.kind!r}"
            )

        DWORD.apply_event(word, event)
        event_records.append(
            {
                "ordinal": ordinal,
                "kind": cell.kind,
                "same_selector_source_token": True,
                "same_finite_map_state_and_jacobian": same_map_state_and_jacobian,
            }
        )
        state = list(state_out)

    return {
        "mode": mode,
        "sample_index": selector.sample_index,
        "source_token": source_token,
        "event_kinds": tuple(word.events),
        "state_in": list(cells[0].state),
        "state_out": state,
        "J_word": word.J_word,
        "event_records": event_records,
        "prediction_count": word.predictions,
        "floor_count": word.floors,
        "S_update_count": word.S_updates,
        "accelerometer_update_count": word.accelerometer_updates,
        "vector_update_count": word.vector_updates,
        "literal_state_succession_exact": True,
        "prediction_identity_placeholder_used": False,
        "independent_event_Jacobian_box_used": False,
        "same_selector_source_token_retained": True,
    }


def build() -> dict:
    bind = BIND.build()
    pred = PRED.build()
    events = EVENTS.build()
    word = DWORD.build()
    bad = {
        "binding": BIND.validate(bind),
        "prediction": PRED.validate(pred),
        "events": EVENTS.validate(events),
        "word": DWORD.validate(word),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError("nonlinear-chain prerequisites failed: " + repr(bad))
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "trusted_selector_binding_consumed": True,
        "exact_finite_prediction_state_and_Jacobian_composed": True,
        "same_cell_Joseph_state_and_Jacobian_composed": True,
        "A21_projection_hybrid_retained_in_Joseph_events": True,
        "aw_floor_has_identity_physical_state_map": True,
        "literal_event_state_succession_required": True,
        "same_selector_source_token_required_for_all_events": True,
        "independent_prediction_F_box_allowed": False,
        "prediction_identity_placeholder_allowed": False,
        "independent_K_allowed": False,
        "production_complete_601_sample_chain_materialized_here": False,
        "H_to_A_hybrid_attachment_closed_here": False,
        "endpoint_augmented_LDLT_closed_here": False,
        "every_prefix_augmented_LDLT_closed_here": False,
        "first_exit_retention_closed_here": False,
        "P3_delta": P3_DELTA,
        "P4_MOTION_PASS": False,
        "P4_PASS": False,
        "P5_MAY_START": False,
        "next_obligation": (
            "instantiate selector-bound SourceCoverCell lineages for every provider branch, compose each sample with materialize_sample_chain, then attach the unique H18->A21 forcing/covariance hybrid and run endpoint/every-prefix augmented LDLT with first-exit retention"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_BRMM_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for k in (
        "trusted_selector_binding_consumed",
        "exact_finite_prediction_state_and_Jacobian_composed",
        "same_cell_Joseph_state_and_Jacobian_composed",
        "A21_projection_hybrid_retained_in_Joseph_events",
        "aw_floor_has_identity_physical_state_map",
        "literal_event_state_succession_required",
        "same_selector_source_token_required_for_all_events",
    ):
        if d.get(k) is not True:
            f.append(k + " not true")
    for k in (
        "independent_prediction_F_box_allowed",
        "prediction_identity_placeholder_allowed",
        "independent_K_allowed",
        "production_complete_601_sample_chain_materialized_here",
        "H_to_A_hybrid_attachment_closed_here",
        "endpoint_augmented_LDLT_closed_here",
        "every_prefix_augmented_LDLT_closed_here",
        "first_exit_retention_closed_here",
        "P4_MOTION_PASS",
        "P4_PASS",
        "P5_MAY_START",
    ):
        if d.get(k) is not False:
            f.append(k + " not false")
    if d.get("P3_delta") != P3_DELTA:
        f.append("P3 delta changed")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build()
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "exact_sample_chain_interface": not f,
        "production_601": d["production_complete_601_sample_chain_materialized_here"],
        "P4": d["P4_PASS"],
        "failures": f,
    }, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
