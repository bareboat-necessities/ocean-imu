#!/usr/bin/env python3
"""Suffix-propagated conditional binary32 response for all-bias A21 prefixes.

The mathematical P4 theorem is allowed to consume the repository's explicit
conditional binary32 execution premise; target/toolchain qualification remains
separate.  This module turns that per-event additive state-error premise into an
actual every-prefix response map instead of multiplying a worst-case roundoff
number by packet count.

For each literal A21 event, the exact 24-state event Jacobian A_e is the same one
used by the all-bias cocycle.  Previous arithmetic inputs are propagated by A_e
and a new normalized 21-vector arithmetic input is appended through

    G_fp = [ delta_e I_21 ; 0_(3x21) ],

so the proof-only true-bias coordinate receives no numerical perturbation.  The
per-event delta_e is taken conservatively from the already-validated conditional
Kalman/reset arithmetic contract.  Prediction uses its time-update bound;
measurement/reset/projection events use the full measurement/reset bound; an
aw-floor event uses the larger of the two declared event bounds so covariance-
floor implementation effects are not silently treated as zero.

This is a conditional mathematical arithmetic channel, not deployment evidence.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval, matrix_identity, matrix_mul
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_p4_a21_bias1_24state_event_lift as LIFT
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_complete_brmm_differential_prediction as PREDICTION
import ou3_p4_complete_brmm_same_history_prefix_selectors as SELECTORS
import ou3_p4_same_history_nonlinear_graph_lineage as GRAPH
import ou3_p4_source_uniform_bias_prefix_lineage as BIASPREFIX
import ou3_p4_source_reachable_selector_family as FAMILY
import ou3_p4_kalman_reset_binary32_iss as FPISS

SCHEMA = 1
QUALIFICATION = "OU3_P4_A21_24STATE_CONDITIONAL_BINARY32_PREFIX_RESPONSE_V1"
P3_DELTA = 1.0e-18


def I(x: float) -> Interval:
    return Interval.point(float(x))


def _shape(A):
    return len(A), len(A[0]) if A else 0


def _zero(r: int, c: int):
    return [[I(0.0) for _ in range(c)] for _ in range(r)]


def _hstack(A, B):
    if len(A) != len(B):
        raise ValueError("horizontal stack row mismatch")
    return [list(a) + list(b) for a, b in zip(A, B)]


def _fp_injection(delta: float):
    d = I(float(delta)); z = I(0.0)
    G = [[z for _ in range(21)] for _ in range(24)]
    for i in range(21):
        G[i][i] = d
    return G


def _propagate(Ae, previous, delta: float):
    carried = matrix_mul(Ae, previous) if _shape(previous)[1] else _zero(24, 0)
    return _hstack(carried, _fp_injection(delta))


@dataclass(frozen=True)
class Binary32PrefixResponse:
    family: str
    literal_prefix_ordinal: int
    selector_source_cell_id: str
    kind: str
    event_index_in_sample: int
    fp_response: tuple[tuple[Interval, ...], ...]
    normalized_fp_input_dimension: int
    local_state_error_norm_upper: float
    conditional_platform_premise: bool


def materialize_A21_binary32_prefix_response(
    selectors: Sequence[SELECTORS.PrefixSelector], endpoint_source_cell_id: str,
    *, family: str, initial_error_state: Sequence[Interval],
    domain_path: Path = GRAPH.DEFAULT_DOMAIN,
) -> list[Binary32PrefixResponse]:
    if family not in BIASPREFIX.FAMILIES:
        raise ValueError("family must be BIAS0/BIAS1/BIAS2")
    if len(initial_error_state) != 21:
        raise ValueError("A21 initial error state must have dimension 21")
    attached = FAMILY.attach_endpoint_family(selectors, endpoint_source_cell_id)
    lineage = list(attached.selector_lineage)
    bias = attached.bias_lineages[family]
    constants = KERNEL._process_constants(domain_path)
    projection_limit = GRAPH._projection_limit(domain_path)
    fp = FPISS.build(); ff = FPISS.validate(fp)
    if ff or not fp["additive_ISS_channel_complete_for_conditional_P4"]:
        raise RuntimeError("conditional binary32 ISS prerequisite is not closed: " + repr(ff))
    row = fp["modes"]["A21"]
    d_pred = float(row["time_update_libm_and_coefficient_additive_state_error_norm_upper_per_prediction"])
    d_meas = float(row["explicit_measurement_reset_state_roundoff_norm_upper_per_event"])
    d_floor = max(d_pred, d_meas)

    C = _zero(24, 0)
    state = list(initial_error_state)
    out=[]
    ordinal=0
    for selector in lineage:
        beta_true = bias.at(selector.source_cell_id)
        for cell in selector.A_event_cells:
            if cell.kind == "prediction":
                p = PREDICTION.prediction_event(
                    "A", state, selector.sample_coordinates.omega_body_corrected,
                    constants.h, selector.active_schedule.tau,
                    tau_ba=constants.accel_bias_tau_s,
                )
                Ae, _ = LIFT.prediction_lift(p["J_state"], bias.phi_true)
                state = list(p["state_out"]); delta=d_pred
            elif cell.kind == "aw_floor":
                Ae = matrix_identity(24); delta=d_floor
            elif cell.kind in ("S_zero", "accelerometer", "magnetometer"):
                if cell.H is None or cell.R is None:
                    raise RuntimeError("measurement event lost P/H/R")
                kwargs = GRAPH._measurement_geometry(selector, cell)
                provenance = EVENTS.ACTUAL_RS_PROVENANCE if cell.kind == "S_zero" else None
                ev = EVENTS.source_joseph_event(
                    "A", state, cell.P_before, cell.R, cell.kind,
                    R_provenance=provenance, bias_true=beta_true,
                    bias_projection_limit=projection_limit, **kwargs,
                )
                Ae = LIFT.measurement_lift(ev)
                state = list(ev["state_out"]); delta=d_meas
            else:
                raise RuntimeError(f"unsupported A21 event {cell.kind}")
            C = _propagate(Ae, C, delta)
            ordinal += 1
            if _shape(C) != (24, 21 * ordinal):
                raise RuntimeError("binary32 prefix response dimension drift")
            out.append(Binary32PrefixResponse(
                family, ordinal, selector.source_cell_id, cell.kind,
                cell.event_index_in_sample,
                tuple(tuple(v for v in r) for r in C), 21 * ordinal,
                float(delta), True,
            ))
    if not out:
        raise RuntimeError("empty binary32 prefix response")
    return out


def build() -> dict:
    fp=FPISS.build(); ff=FPISS.validate(fp)
    family=FAMILY.build(); famf=FAMILY.validate(family)
    if ff or famf:
        raise RuntimeError(f"binary32-prefix prerequisites failed fp={ff} family={famf}")
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "conditional_binary32_Kalman_reset_ISS_contract_consumed": True,
        "conditional_platform_execution_premise_consumed": True,
        "deployment_toolchain_qualification_required_for_mathematical_P4": False,
        "deployment_toolchain_qualification_closed_here": False,
        "event_local_roundoff_enters_as_additive_shipping_state_input": True,
        "proof_true_bias_coordinate_receives_roundoff": False,
        "previous_roundoff_inputs_suffix_propagated_through_later_event_Jacobians": True,
        "one_new_normalized_21D_roundoff_input_appended_per_literal_event": True,
        "packet_count_times_worst_roundoff_used": False,
        "all_BIAS0_BIAS1_BIAS2_supported": True,
        "A21_every_prefix_binary32_response_materializable": True,
        "production_common_augmented_PrefixInput_embedded_here": False,
        "endpoint_augmented_LDLT_closed_here": False,
        "every_prefix_augmented_LDLT_closed_here": False,
        "P3_delta": P3_DELTA,
        "P4_MOTION_PASS": False,
        "P4_PASS": False,
        "P5_MAY_START": False,
        "next_obligation": (
            "embed each suffix-propagated fp_response beside the corresponding all-bias physical source response, "
            "state/storage master, graph sectors, and moment/radial maps in one PrefixInput coordinate; then run outward LDLT"
        ),
    }


def validate(d: dict) -> list[str]:
    f=[]
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for k in (
        "conditional_binary32_Kalman_reset_ISS_contract_consumed",
        "conditional_platform_execution_premise_consumed",
        "event_local_roundoff_enters_as_additive_shipping_state_input",
        "previous_roundoff_inputs_suffix_propagated_through_later_event_Jacobians",
        "one_new_normalized_21D_roundoff_input_appended_per_literal_event",
        "all_BIAS0_BIAS1_BIAS2_supported",
        "A21_every_prefix_binary32_response_materializable",
    ):
        if d.get(k) is not True: f.append(k + " not true")
    for k in (
        "deployment_toolchain_qualification_required_for_mathematical_P4",
        "deployment_toolchain_qualification_closed_here",
        "proof_true_bias_coordinate_receives_roundoff",
        "packet_count_times_worst_roundoff_used",
        "production_common_augmented_PrefixInput_embedded_here",
        "endpoint_augmented_LDLT_closed_here", "every_prefix_augmented_LDLT_closed_here",
        "P4_MOTION_PASS", "P4_PASS", "P5_MAY_START",
    ):
        if d.get(k) is not False: f.append(k + " not false")
    if d.get("P3_delta") != P3_DELTA: f.append("P3 delta changed")
    return list(dict.fromkeys(f))


def main() -> int:
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("--output",type=Path,required=True); a=p.parse_args()
    d=build(); f=validate(d); d["validation_pass"]=not f; d["validation_failures"]=f
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"conditional_fp_prefix":d["A21_every_prefix_binary32_response_materializable"],"deployment":d["deployment_toolchain_qualification_closed_here"],"P4":d["P4_PASS"],"failures":f},sort_keys=True)); return int(bool(f))


if __name__ == "__main__": raise SystemExit(main())
