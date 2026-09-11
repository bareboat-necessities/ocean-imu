#!/usr/bin/env python3
"""Exact translational observation-kernel lemma for strengthened COMPLETE-BRMM.

The IMU measures acceleration, not the integration constants of velocity and
position.  The finite-window COMPLETE-BRMM source therefore admitted the exact
quiet ambiguity exposed by ``ou3_brmm_infinite_continuation``.  Once the
origin-invariant indefinite centered-S recurrence is qualified, however, the
translation-only zero-output kernel collapses exactly.

For two physical histories with identical non-gravitational acceleration after
the common Live handoff, their differences satisfy

    dv(h) = c_v,
    dp(h) = c_p + h c_v,
    dS_L(h) = h c_p + h^2 c_v / 2,

because both proof references use S_L(t_L)=0.  If each history has uniformly
bounded physical position, then dp must remain bounded and hence c_v=0.  If
each history also has uniformly bounded handoff-centered S, then dS_L must
remain bounded and hence c_p=0.  Thus no nonzero velocity/position integration
constant can remain invisible to the sensor history inside the strengthened
source family.

This is an exact uniqueness/necessary-structure result, not a P4 contraction
certificate.  It does not assign D_S, does not bound absolute S, does not use
covariance membership, and does not modify the shipping estimator.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction as F
from pathlib import Path

import ou3_brmm_contract as BRMM
import ou3_brmm_centered_S_recurrence as SREC

SCHEMA = 1
QUALIFICATION = "OU3_BRMM_TRANSLATION_OBSERVATION_KERNEL_V1"


def _poly_trim(p: tuple[F, ...]) -> tuple[F, ...]:
    q = list(p)
    while len(q) > 1 and q[-1] == 0:
        q.pop()
    return tuple(q)


def _poly_add(a: tuple[F, ...], b: tuple[F, ...]) -> tuple[F, ...]:
    n = max(len(a), len(b))
    out = [F(0) for _ in range(n)]
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] += x
    return _poly_trim(tuple(out))


def _poly_integral(a: tuple[F, ...], constant: F = F(0)) -> tuple[F, ...]:
    return _poly_trim((constant,) + tuple(a[i] / F(i + 1) for i in range(len(a))))


def _encoded(p: tuple[F, ...]) -> list[str]:
    return [str(x) for x in _poly_trim(p)]


def build() -> dict:
    source = BRMM.build()
    recurrence = SREC.build()
    if BRMM.validate(source):
        raise RuntimeError("BRMM source prerequisite failed")
    if SREC.validate(recurrence):
        raise RuntimeError("centered-S recurrence prerequisite failed")
    if recurrence["P3_delta"] != 1e-18:
        raise RuntimeError("canonical P3 changed")

    # Work one scalar axis; Cartesian axes are independent for these exact
    # additive primitive identities.  c_v and c_p are symbolic basis columns.
    # Polynomial coefficients are ordered [1,h,h^2,...].
    dv_cv = (F(1),)
    dp_cv = _poly_integral(dv_cv)             # h
    dS_cv = _poly_integral(dp_cv)             # h^2/2
    dv_cp = (F(0),)
    dp_cp = (F(1),)
    dS_cp = _poly_integral(dp_cp)             # h

    exact = (
        dp_cv == (F(0), F(1))
        and dS_cv == (F(0), F(0), F(1, 2))
        and dp_cp == (F(1),)
        and dS_cp == (F(0), F(1))
    )
    if not exact:
        raise RuntimeError("primitive polynomial identities failed")

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "same_acceleration_history_assumed_for_kernel": True,
        "common_fresh_handoff_centered_S_zero": True,
        "difference_basis": {
            "velocity_constant_c_v": {
                "delta_v": _encoded(dv_cv),
                "delta_p": _encoded(dp_cv),
                "delta_S_L": _encoded(dS_cv),
            },
            "position_constant_c_p": {
                "delta_v": _encoded(dv_cp),
                "delta_p": _encoded(dp_cp),
                "delta_S_L": _encoded(dS_cp),
            },
        },
        "bounded_position_forces_c_v_zero": True,
        "bounded_centered_S_forces_c_p_zero_after_c_v_zero": True,
        "translation_zero_output_kernel_dimension_before_indefinite_primitives_per_axis": 2,
        "translation_zero_output_kernel_dimension_after_bounded_position_per_axis": 1,
        "translation_zero_output_kernel_dimension_after_centered_S_recurrence_per_axis": 0,
        "all_three_axes_kernel_collapses": True,
        "absolute_S_origin_needed_for_kernel_collapse": False,
        "legacy_300_m_s_entry_ball_used": False,
        "position_reanchoring_used": False,
        "wordwise_S_rezero_used": False,
        "D_S_numeric_value_used": None,
        "only_finiteness_of_D_S_used": True,
        "covariance_membership_used": False,
        "shipping_filter_changed": False,
        "P3_delta": 1e-18,
        "P4_PASS": False,
        "P5_MAY_START": False,
        "result": (
            "conditional on finite uniform handoff-centered S recurrence, identical "
            "acceleration histories cannot differ by nonzero bounded translational "
            "integration constants"
        ),
        "next_obligation": (
            "enclose nonzero same-history acceleration/source variations and the "
            "shipping correction graph; this lemma removes only the exact unobservable "
            "translation kernel"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for k in (
        "same_acceleration_history_assumed_for_kernel",
        "common_fresh_handoff_centered_S_zero",
        "bounded_position_forces_c_v_zero",
        "bounded_centered_S_forces_c_p_zero_after_c_v_zero",
        "all_three_axes_kernel_collapses",
        "only_finiteness_of_D_S_used",
    ):
        if d.get(k) is not True:
            f.append(k + " not true")
    for k in (
        "absolute_S_origin_needed_for_kernel_collapse",
        "legacy_300_m_s_entry_ball_used",
        "position_reanchoring_used",
        "wordwise_S_rezero_used",
        "covariance_membership_used",
        "shipping_filter_changed",
        "P4_PASS",
        "P5_MAY_START",
    ):
        if d.get(k) is not False:
            f.append(k + " not false")
    if d.get("D_S_numeric_value_used") is not None:
        f.append("numeric D_S must not enter the kernel lemma")
    if d.get("translation_zero_output_kernel_dimension_before_indefinite_primitives_per_axis") != 2:
        f.append("wrong initial translation kernel dimension")
    if d.get("translation_zero_output_kernel_dimension_after_bounded_position_per_axis") != 1:
        f.append("bounded-position kernel dimension wrong")
    if d.get("translation_zero_output_kernel_dimension_after_centered_S_recurrence_per_axis") != 0:
        f.append("centered-S kernel did not collapse")
    if d.get("P3_delta") != 1e-18:
        f.append("P3 delta changed")
    basis = d.get("difference_basis", {})
    try:
        if basis["velocity_constant_c_v"]["delta_p"] != ["0", "1"]:
            f.append("c_v position polynomial changed")
        if basis["velocity_constant_c_v"]["delta_S_L"] != ["0", "0", "1/2"]:
            f.append("c_v S polynomial changed")
        if basis["position_constant_c_p"]["delta_p"] != ["1"]:
            f.append("c_p position polynomial changed")
        if basis["position_constant_c_p"]["delta_S_L"] != ["0", "1"]:
            f.append("c_p S polynomial changed")
    except Exception:
        f.append("difference basis missing")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build()
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "kernel_before": d["translation_zero_output_kernel_dimension_before_indefinite_primitives_per_axis"],
        "after_P_bound": d["translation_zero_output_kernel_dimension_after_bounded_position_per_axis"],
        "after_centered_S": d["translation_zero_output_kernel_dimension_after_centered_S_recurrence_per_axis"],
        "P4": d["P4_PASS"],
        "failures": f,
    }, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
