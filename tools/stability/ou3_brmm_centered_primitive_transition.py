#!/usr/bin/env python3
"""Exact same-history centered translational primitive transition for BRMM.

The theorem-facing source must not pass only instantaneous acceleration samples to
P4.  The same physical history also owns velocity, position and the one-time
Live-centered integral displacement

    S_L(t) = S(t) - S(t_L),        S_L(t_L) = 0.

For one physical interval of duration h, define the three *coupled* acceleration
moments of one vector history a(s), 0<=s<=h,

    J0 = integral a(s) ds,
    J1 = integral (h-s) a(s) ds,
    J2 = integral (h-s)^2 a(s)/2 ds.

Then the exact primitive recurrence is

    v1   = v0 + J0,
    p1   = p0 + h v0 + J1,
    S_L1 = S_L0 + h p0 + h^2 v0/2 + J2.

J0/J1/J2 are not independent disturbance ports.  They are three linear
functionals of ONE acceleration history and must retain one source-transition
witness.  Scalar norm bounds are emitted only as consequences useful for
screening; they cannot generate theorem cells.

This module is exact rational algebra.  It changes no shipping code and uses no
absolute S origin, fresh 300 m*s entrance ball, covariance membership, or
wordwise S reset.
"""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path
from typing import Sequence

import ou3_brmm_contract as BRMM
import ou3_brmm_centered_S_recurrence as SREC

SCHEMA = 1
QUALIFICATION = "OU3_BRMM_CENTERED_PRIMITIVE_TRANSITION_V1"
H = F(1, 200)  # shipping/proof source period, exactly 5 ms


@dataclass(frozen=True)
class MomentWitness:
    """Three moments owned by one source transition; never Cartesianized."""

    source_transition_witness_id: str
    J0: tuple[F, F, F]
    J1: tuple[F, F, F]
    J2: tuple[F, F, F]


@dataclass(frozen=True)
class PrimitiveState:
    primitive_id: str
    centered_S_origin_witness_id: str
    v: tuple[F, F, F]
    p: tuple[F, F, F]
    S_L: tuple[F, F, F]


def _v3(values: Sequence[F | int]) -> tuple[F, F, F]:
    if len(values) != 3:
        raise ValueError("three-vector required")
    return tuple(F(x) for x in values)  # type: ignore[return-value]


def advance(
    state: PrimitiveState,
    moments: MomentWitness,
    *,
    next_primitive_id: str,
    h: F = H,
) -> PrimitiveState:
    if not state.primitive_id or not state.centered_S_origin_witness_id:
        raise ValueError("primitive/source origin ids required")
    if not moments.source_transition_witness_id or not next_primitive_id:
        raise ValueError("source transition and next primitive ids required")
    if h <= 0:
        raise ValueError("positive interval required")
    v1 = tuple(state.v[i] + moments.J0[i] for i in range(3))
    p1 = tuple(state.p[i] + h * state.v[i] + moments.J1[i] for i in range(3))
    s1 = tuple(
        state.S_L[i] + h * state.p[i] + h * h * state.v[i] / 2 + moments.J2[i]
        for i in range(3)
    )
    return PrimitiveState(
        next_primitive_id,
        state.centered_S_origin_witness_id,
        _v3(v1),
        _v3(p1),
        _v3(s1),
    )


def _matrix() -> list[list[str]]:
    # scalar-axis map [v,p,S,J0,J1,J2] -> [v+,p+,S+]
    A = [
        [F(1), F(0), F(0), F(1), F(0), F(0)],
        [H, F(1), F(0), F(0), F(1), F(0)],
        [H * H / 2, H, F(1), F(0), F(0), F(1)],
    ]
    return [[str(x) for x in row] for row in A]


def build() -> dict:
    source = BRMM.build()
    recurrence = SREC.build()
    if BRMM.validate(source):
        raise RuntimeError("BRMM declaration failed")
    if SREC.validate(recurrence):
        raise RuntimeError("centered-S recurrence prerequisite failed")
    Amax = F(str(source["hard_caps"]["acceleration_norm_mps2"]))

    # Exact mutation-sensitive smoke: nonzero root and moments on all axes.
    root = PrimitiveState(
        "q0", "live-origin",
        _v3((F(1, 3), F(-2, 5), F(7, 11))),
        _v3((F(2, 7), F(3, 8), F(-5, 9))),
        _v3((F(0), F(0), F(0))),
    )
    m = MomentWitness(
        "tr0",
        _v3((F(1, 1000), F(-1, 2000), F(3, 4000))),
        _v3((F(1, 200000), F(1, 300000), F(-1, 250000))),
        _v3((F(1, 60000000), F(-1, 80000000), F(1, 90000000))),
    )
    child = advance(root, m, next_primitive_id="q1")

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "sample_period_exact_s": str(H),
        "state_order_per_axis": ["v", "p", "S_L"],
        "same_history_acceleration_moments": {
            "J0": "integral_0^h a(s) ds",
            "J1": "integral_0^h (h-s) a(s) ds",
            "J2": "integral_0^h (h-s)^2 a(s)/2 ds",
            "one_source_transition_witness_required": True,
            "may_be_selected_independently": False,
        },
        "exact_scalar_axis_transition_matrix": _matrix(),
        "exact_recurrence": {
            "v_next": "v + J0",
            "p_next": "p + h*v + J1",
            "S_L_next": "S_L + h*p + h^2*v/2 + J2",
        },
        "moment_norm_consequences_from_pointwise_acceleration_cap": {
            "A_max_mps2": str(Amax),
            "J0_norm_upper": str(Amax * H),
            "J1_norm_upper": str(Amax * H * H / 2),
            "J2_norm_upper": str(Amax * H * H * H / 6),
            "consequences_are_not_independent_ports": True,
        },
        "cross_word_primitive_out_is_next_word_primitive_in_required": True,
        "centered_S_origin_id_is_invariant_under_advance": (
            child.centered_S_origin_witness_id == root.centered_S_origin_witness_id
        ),
        "wordwise_S_rezero_used": False,
        "absolute_S_origin_used": False,
        "legacy_300_m_s_fresh_entry_ball_used": False,
        "covariance_membership_used": False,
        "shipping_filter_changed": False,
        "D_S_numeric_value_used": None,
        "P3_delta": 1e-18,
        "P4_PASS": False,
        "P5_MAY_START": False,
        "smoke": {
            "root": {
                "v": [str(x) for x in root.v],
                "p": [str(x) for x in root.p],
                "S_L": [str(x) for x in root.S_L],
            },
            "child": {
                "v": [str(x) for x in child.v],
                "p": [str(x) for x in child.p],
                "S_L": [str(x) for x in child.S_L],
            },
        },
        "primitive_transition_operator_materialized": True,
        "full_complete_BRMM_moment_set_materialized_here": False,
        "source_uniform_endpoint_or_prefix_closed_here": False,
        "next_obligation": (
            "bind every provider transition's v/p/S_L ingress and egress to this exact map "
            "with one coupled J0/J1/J2 source witness, then carry that primitive state into "
            "the same event-local P/H/R/K and joint24 storage lineage"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("sample_period_exact_s") != "1/200":
        f.append("sample period changed")
    moments = d.get("same_history_acceleration_moments", {})
    if moments.get("one_source_transition_witness_required") is not True:
        f.append("moment witness detached")
    if moments.get("may_be_selected_independently") is not False:
        f.append("moments became independent ports")
    if d.get("exact_scalar_axis_transition_matrix") != _matrix():
        f.append("exact primitive transition matrix changed")
    for k in (
        "cross_word_primitive_out_is_next_word_primitive_in_required",
        "centered_S_origin_id_is_invariant_under_advance",
        "primitive_transition_operator_materialized",
    ):
        if d.get(k) is not True:
            f.append(k + " not true")
    for k in (
        "wordwise_S_rezero_used",
        "absolute_S_origin_used",
        "legacy_300_m_s_fresh_entry_ball_used",
        "covariance_membership_used",
        "shipping_filter_changed",
        "full_complete_BRMM_moment_set_materialized_here",
        "source_uniform_endpoint_or_prefix_closed_here",
        "P4_PASS",
        "P5_MAY_START",
    ):
        if d.get(k) is not False:
            f.append(k + " not false")
    if d.get("D_S_numeric_value_used") is not None:
        f.append("numeric D_S injected into primitive transition")
    if d.get("P3_delta") != 1e-18:
        f.append("P3 delta changed")
    c = d.get("moment_norm_consequences_from_pointwise_acceleration_cap", {})
    if c.get("consequences_are_not_independent_ports") is not True:
        f.append("moment consequences became generator ports")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", required=True, type=Path)
    a = ap.parse_args()
    d = build()
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "operator": d["primitive_transition_operator_materialized"],
        "moments_independent": d["same_history_acceleration_moments"]["may_be_selected_independently"],
        "P4": d["P4_PASS"],
        "failures": f,
    }, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
