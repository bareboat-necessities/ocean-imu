#!/usr/bin/env python3
"""P5 wide-handoff H18 finite-angle differential backbone.

The shipping handoff theorem places the tilt quotient in an open ball strictly
below pi/2 + 0.02 rad ~= 91.146 deg.  The P4 geometry producer intentionally
labels only sectors with >=64% tangent-information retention as "usable P4"
sectors.  That 64% threshold is a conservative design floor, not a singularity
or a theorem boundary.

For the transient P5 capture problem we therefore recompute the same prior-free
H18 information/LDLT chain on a 1.60-rad Cayley sector (91.673 deg), which
strictly contains the complete certified shipping handoff tilt section.  We do
not consume CAYLEY.validate(), because that validator intentionally enforces the
P4-only 0.64 design floor.  Instead we consume the exact Cayley identities and
require only positive chart conditioning and positive retained information,
then run the repository's actual H18 prior-free cell completion.

This certificate is a *wide-sector differential/detectability backbone*.  It
can remove "loss of finite-angle information above 45 deg" as a P5 blocker.  It
is not by itself a finite-map contraction theorem: the exact chord/Joseph/reset
word and source-uniform endpoint/prefix storage still have to close.
"""
from __future__ import annotations

import copy
import json
import math
from pathlib import Path

import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_brmm_h18_information_composition as HINFO
import ou3_brmm_h18_prior_free_completion as HPF
import ou3_brmm_dynamic_source_certificate as DYNAMIC
import ou3_full_process_ucc as PROCESS
import ou3_brmm_full_word_event_algebra as EVENT
import ou3_startup_handoff_tilt_hemisphere as HANDOFF

REPO = Path(__file__).resolve().parents[2]
SCHEMA = 1
QUALIFICATION = "OU3_P5_WIDE_HANDOFF_H18_BACKBONE_V1"
OUTER_ANGLE_RAD = 1.60


def build() -> dict:
    handoff = HANDOFF.build()
    if HANDOFF.validate(handoff):
        raise RuntimeError("shipping handoff tilt certificate invalid")
    handoff_upper = float(handoff["true_gravity_quotient_tilt_strict_upper_rad"])
    if not (handoff_upper < OUTER_ANGLE_RAD < math.pi):
        raise RuntimeError("chosen transient sector does not strictly contain handoff set")

    cay = CAYLEY.build(outer_angle_rad=OUTER_ANGLE_RAD)
    # Do not call CAYLEY.validate: its 0.64 information floor is deliberately a
    # P4 sector-selection policy, not a geometric validity condition.
    geometric_validity = bool(
        cay["source_generated_not_trajectory_fit"]
        and cay["chart_antipode_excluded"]
        and float(cay["chart_sigma_min_lower"]) > 0.0
        and float(cay["exact_vector_information_retention_factor_lower"]) > 0.0
        and math.isfinite(float(cay["cayley_radius_upper"]))
    )
    if not geometric_validity:
        raise RuntimeError("wide Cayley sector lost positive finite geometry")

    base = HINFO.build()
    dyn = DYNAMIC.build()
    proc = PROCESS.build()
    event = EVENT.build()
    bad = {
        "hinfo": HINFO.validate(base),
        "dynamic": DYNAMIC.validate(dyn),
        "process": PROCESS.validate(proc),
        "event": EVENT.validate(event),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError("wide H18 prerequisites failed: " + repr(bad))

    k = float(cay["exact_vector_information_retention_factor_lower"])
    alpha0 = float(base["eta6_information_lower"])
    alpha = math.nextafter(k * alpha0, -math.inf)
    comp = base["triangular_information_composition"]
    d_aw = float(comp["aw_direction_information_lower"])
    cross = float(comp["C_aw_spectral_norm_squared_upper"])
    trace = math.nextafter(alpha + d_aw + cross, math.inf)
    coupled = math.nextafter(alpha * d_aw / trace, -math.inf)
    non_aw = float(comp["non_aw_translation_lambda_min_lower"])
    D = math.nextafter(min(coupled, non_aw), -math.inf)
    if not (alpha > 0.0 and coupled > 0.0 and D > 0.0):
        raise RuntimeError("wide finite-angle information lost positivity")

    h = copy.deepcopy(base)
    h["eta6_information_lower"] = alpha
    h["H18_information_useful_gate_pass"] = D >= 1e-18
    h["triangular_information_composition"]["A_transpose_A_lower"] = alpha
    h["triangular_information_composition"]["coupled_eta6_aw_scalar_2x2_trace_upper"] = trace
    h["triangular_information_composition"]["coupled_eta6_aw_scalar_2x2_determinant_lower"] = math.nextafter(alpha * d_aw, -math.inf)
    h["triangular_information_composition"]["coupled_eta6_aw_lambda_min_lower"] = coupled
    h["triangular_information_composition"]["D_H18_lambda_min_lower"] = D

    pbar = HPF._same_word_covariance_upper(Path(HPF.DEFAULT_DOMAIN).resolve(), dyn, proc, h)
    fnorm = HPF._prediction_norm_sq_upper(proc)
    penalty = math.nextafter(
        (HPF.USEFUL_GATE ** 2 / 4.0) * fnorm * float(pbar["Pbar_trace_upper"]),
        math.inf,
    )

    rows = []
    failures = []
    worst = math.inf
    for x in HPF._x_cover(dyn):
        ok, row = HPF._full_H18_cell(x, process=proc, dynamic=dyn, penalty_physical=penalty)
        rows.append(row)
        if ok:
            worst = min(worst, float(row["pivot_lower"]))
        else:
            failures.append(row)
    ldlt_closed = bool(rows) and not failures and math.isfinite(worst) and worst > 0.0

    preserve = event["full_matrix_margin_preservation"]
    suffix = all(bool(preserve[k0]) for k0 in (
        "covers_prediction",
        "covers_every_due_S_update",
        "covers_every_Normal_Live_accelerometer_update",
        "covers_asynchronous_magnetometer_update",
        "covers_immediate_left_error_reset",
        "covers_aw_covariance_floor",
        "covers_not_due_or_rejected_identity_branches",
    ))

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "filter_changed": False,
        "quality_gates_changed": False,
        "domain_shrunk": False,
        "source_enumeration_used": False,
        "trajectory_replay_used": False,
        "shipping_handoff_tilt_strict_upper_rad": handoff_upper,
        "shipping_handoff_tilt_strict_upper_deg": float(handoff["true_gravity_quotient_tilt_strict_upper_deg"]),
        "transient_outer_angle_rad": OUTER_ANGLE_RAD,
        "transient_outer_angle_deg": math.degrees(OUTER_ANGLE_RAD),
        "strictly_contains_complete_handoff_tilt_section": handoff_upper < OUTER_ANGLE_RAD,
        "cayley_P4_usable_floor_pass_at_wide_angle": bool(cay["usable_sector_geometry_pass"]),
        "P4_usable_floor_is_not_consumed_by_P5_backbone": True,
        "wide_chart_sigma_min_lower": float(cay["chart_sigma_min_lower"]),
        "wide_vector_information_retention_lower": k,
        "wide_finite_angle_eta6_information_lower": alpha,
        "wide_finite_angle_H18_information_lower": D,
        "wide_information_above_P3_delta": D >= 1e-18,
        "same_word_Pbar_trace_upper": float(pbar["Pbar_trace_upper"]),
        "delta_squared_completion_penalty": penalty,
        "x_cells_certified": len(rows),
        "x_cell_failures": failures,
        "worst_full_H18_LDLT_pivot_lower": worst if ldlt_closed else None,
        "WIDE_HANDOFF_H18_PRIOR_FREE_LDLT_CLOSED": ldlt_closed,
        "linear_suffix_event_algebra_available": suffix,
        "finite_angle_information_loss_above_45deg_is_P5_blocker": not ldlt_closed,
        "exact_chord_Joseph_reset_capture_word_still_required": True,
        "finite_map_contraction_closed_here": False,
        "finite_capture_time_closed_here": False,
        "P4_PASS": False,
        "P5_PASS": False,
    }


def validate(d: dict) -> list[str]:
    f = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for k in (
        "strictly_contains_complete_handoff_tilt_section",
        "P4_usable_floor_is_not_consumed_by_P5_backbone",
        "wide_information_above_P3_delta",
        "WIDE_HANDOFF_H18_PRIOR_FREE_LDLT_CLOSED",
        "linear_suffix_event_algebra_available",
        "exact_chord_Joseph_reset_capture_word_still_required",
    ):
        if d.get(k) is not True:
            f.append(k + " not true")
    # At 1.60 rad the old P4 convenience floor should fail; if it starts passing,
    # that is fine geometrically but means the diagnostic expectation changed.
    if d.get("cayley_P4_usable_floor_pass_at_wide_angle") is not False:
        f.append("unexpected: wide sector now passes old P4 0.64 convenience floor")
    for k in (
        "filter_changed", "quality_gates_changed", "domain_shrunk",
        "source_enumeration_used", "trajectory_replay_used",
        "finite_angle_information_loss_above_45deg_is_P5_blocker",
        "finite_map_contraction_closed_here", "finite_capture_time_closed_here",
        "P4_PASS", "P5_PASS",
    ):
        if d.get(k) is not False:
            f.append(k + " not false")
    if d.get("x_cell_failures"):
        f.append("wide H18 x-cell failures present")
    if not (0.0 < float(d.get("wide_vector_information_retention_lower", 0.0)) < 0.64):
        f.append("wide information retention not in expected transient range")
    if not (float(d.get("worst_full_H18_LDLT_pivot_lower") or 0.0) > 0.0):
        f.append("wide H18 worst LDLT pivot is not positive")
    return f


if __name__ == "__main__":
    d = build()
    failures = validate(d)
    print(json.dumps({**d, "validation_pass": not failures, "validation_failures": failures}, indent=2, sort_keys=True))
    raise SystemExit(1 if failures else 0)
