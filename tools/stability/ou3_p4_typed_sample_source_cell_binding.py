#!/usr/bin/env python3
"""Bind one trusted typed executor sample to its literal nonlinear P4 cells.

The trusted 601-sample kernel already emits a branch-correlated PrefixSelector
containing the exact SampleCoordinates and event-local Riccati cells produced by
one `advance_branch` call.  The nonlinear P4 source-cover uses SourceCoverCell
objects.  This module closes the representation seam between them: every
literal P4 event for a sample must use the same sample index/estimator ancestry,
and accelerometer geometry must equal that typed sample's f_cog_body/R_wb
exactly.  The event ordering must equal the trusted kernel event slice.

No P/H/R/K is reconstructed here and no favorable branch is selected.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_p4_complete_brmm_same_history_prefix_selectors as SELECTORS
import ou3_p4_complete_brmm_source_cover_contract as COVER

SCHEMA = 1
QUALIFICATION = "OU3_P4_TYPED_SAMPLE_SOURCE_CELL_BINDING_V1"
P3_DELTA = 1.0e-18


def _same_interval(a: Interval, b: Interval) -> bool:
    return isinstance(a, Interval) and isinstance(b, Interval) and a.lo == b.lo and a.hi == b.hi


def _same_vec(a, b) -> bool:
    return len(a) == len(b) and all(_same_interval(x, y) for x, y in zip(a, b))


def _same_mat(A, B) -> bool:
    return len(A) == len(B) and all(_same_vec(a, b) for a, b in zip(A, B))


def validate_event_cells_against_selector(
    selector: SELECTORS.PrefixSelector,
    cells: Sequence[COVER.SourceCoverCell],
    *,
    mode: str,
) -> list[str]:
    """Validate literal nonlinear cells against one trusted kernel sample slice."""
    f: list[str] = []
    if mode not in ("H", "A"):
        return ["mode must be H/A"]
    expected_events = selector.H_events_this_sample if mode == "H" else selector.A_events_this_sample
    expected_cells = selector.H_event_cells if mode == "H" else selector.A_event_cells
    if len(cells) != len(expected_events):
        f.append("literal source-cell count differs from trusted event slice")
        return f
    if tuple(c.kind for c in cells) != tuple(expected_events):
        f.append("literal source-cell order differs from trusted event slice")
    if tuple(ec.kind for ec in expected_cells) != tuple(expected_events):
        f.append("trusted event-cell slice internally inconsistent")

    for i, c in enumerate(cells):
        cf = COVER.validate_cell(c, require_estimator_provenance=True)
        f.extend(f"event {i}: {x}" for x in cf)
        if c.mode != mode:
            f.append(f"event {i}: mode detached from trusted selector")
        if c.sample_index != selector.sample_index:
            f.append(f"event {i}: sample index detached from trusted selector")
        if c.event_ordinal != i:
            f.append(f"event {i}: event ordinal differs from trusted event slice")
        if c.estimator_source_token != selector.source_cell_id:
            f.append(f"event {i}: estimator/source ancestry detached from trusted selector")

        trusted = expected_cells[i]
        if not _same_mat(c.P, trusted.P_before):
            f.append(f"event {i}: source-cell P is not trusted event P_before")
        if c.kind in ("S_zero", "accelerometer", "magnetometer"):
            if c.R is None or trusted.R is None or not _same_mat(c.R, trusted.R):
                f.append(f"event {i}: Joseph R detached from trusted event cell")

        if c.kind == "accelerometer":
            sample = selector.sample_coordinates
            if c.f_hat is None or not _same_vec(c.f_hat, sample.f_cog_body):
                f.append(f"event {i}: f_hat detached from typed sample f_cog_body")
            if c.R_hat is None or not _same_mat(c.R_hat, sample.R_wb):
                f.append(f"event {i}: R_hat detached from typed sample R_wb")
            if trusted.H is None:
                f.append(f"event {i}: trusted accelerometer H missing")
        if c.kind == "magnetometer":
            # The selector retains asynchronous events in their actual order.
            if trusted.magnetic_event_index is None:
                f.append(f"event {i}: trusted magnetic event index missing")
            else:
                mi = trusted.magnetic_event_index
                mags = selector.sample_coordinates.magnetometer_events_after_imu
                if not (0 <= mi < len(mags)):
                    f.append(f"event {i}: magnetic event index outside typed sample")
                elif c.m_body is None or not _same_vec(c.m_body, mags[mi].m_body):
                    f.append(f"event {i}: magnetic geometry detached from typed sample")
    return list(dict.fromkeys(f))


def build() -> dict:
    selectors = SELECTORS.build()
    sf = SELECTORS.validate(selectors)
    cover = COVER.build()
    cf = COVER.validate(cover)
    if sf or cf:
        raise RuntimeError(f"typed-sample binding prerequisites failed selectors={sf} cover={cf}")
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "trusted_prefix_selector_sample_coordinates_consumed": True,
        "trusted_event_local_P_H_R_cells_consumed": True,
        "nonlinear_source_cell_P_must_equal_trusted_event_P_before": True,
        "accelerometer_f_hat_must_equal_typed_sample_f_cog_body": True,
        "accelerometer_R_hat_must_equal_typed_sample_R_wb": True,
        "magnetic_geometry_must_equal_indexed_typed_sample_event": True,
        "literal_event_order_must_equal_trusted_kernel_slice": True,
        "estimator_source_ancestry_must_equal_trusted_selector_child": True,
        "P_H_R_K_reconstructed_independently": False,
        "favorable_frontend_successor_selected": False,
        "shipping_transition_reimplemented": False,
        "production_complete_601_sample_binding_closed_here": False,
        "production_all_H18_cells_bound_here": False,
        "production_all_A21_cells_bound_here": False,
        "P3_delta": P3_DELTA,
        "P4_MOTION_PASS": False,
        "P4_PASS": False,
        "P5_MAY_START": False,
        "next_obligation": (
            "for every retained provider/typed-kernel branch, instantiate SourceCoverCell objects directly from the selector event-local cells and require validate_event_cells_against_selector for H18 and A21 before assembling nonlinear prefix Jacobians/LDLT"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for k in (
        "trusted_prefix_selector_sample_coordinates_consumed",
        "trusted_event_local_P_H_R_cells_consumed",
        "nonlinear_source_cell_P_must_equal_trusted_event_P_before",
        "accelerometer_f_hat_must_equal_typed_sample_f_cog_body",
        "accelerometer_R_hat_must_equal_typed_sample_R_wb",
        "magnetic_geometry_must_equal_indexed_typed_sample_event",
        "literal_event_order_must_equal_trusted_kernel_slice",
        "estimator_source_ancestry_must_equal_trusted_selector_child",
    ):
        if d.get(k) is not True:
            f.append(k + " not true")
    for k in (
        "P_H_R_K_reconstructed_independently", "favorable_frontend_successor_selected",
        "shipping_transition_reimplemented", "production_complete_601_sample_binding_closed_here",
        "production_all_H18_cells_bound_here", "production_all_A21_cells_bound_here",
        "P4_MOTION_PASS", "P4_PASS", "P5_MAY_START",
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
    print(json.dumps({"binding_interface": not f, "production": d["production_complete_601_sample_binding_closed_here"], "P4": d["P4_PASS"], "failures": f}, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
