#!/usr/bin/env python3
"""V32: exact theta-x Delta-C row refinement over the V31 current-subbox lift.

V31 leaves 15 of the 64 authoritative current-Cayley subboxes open after the
V30 rowwise gain refinement.  V30 still bounded the theta-x row of
``Delta C_theta`` by V12D's full attitude-block operator bound.  In the exact
V4/V11 proof gauge the nominal sample-1 attitude covariance is

    P_theta = L_x(d) diag(a,a,Y) L_x(d)^T,

and ``L_x(d)`` has first row/column ``[1,0,0]``.  Therefore the theta-x
covariance row is exactly ``[a,0,0]`` and

    ||e_x^T Delta C_theta||
      <= dP ||H_theta|| + a dH + dP dH + dP.

The final ``dP`` term is the theta/a_w cross-block perturbation already present
in V12D.  V32 substitutes this row-specific ``Delta C_x`` into the same V30
resolvent

    ||Delta K_x|| <= (dC_x + k_parallel dS) ||S'^{-1}||,

then reruns V31 unchanged over all 64 current subboxes.  Y/Z gain perturbation,
V29 yz/radial parents, V28 signed residual source geometry, V23 exact current
residuals, V16/V15/V18 composition, source domain, estimator, six-radian limit,
and theorem-promotion state are unchanged.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_sample1_v30_current_subbox_lift_v31 as V31

DEFAULT_DOMAIN = V31.DEFAULT_DOMAIN
SCHEMA = 3200
FULL = V31.FULL
Q_TARGET = V31.Q_TARGET
WITNESS = V31.WITNESS


def _sum_up(*xs: float) -> float:
    y = 0.0
    for x in xs:
        y = FULL.up(y + float(x))
    return y


def _theta_x_deltac_upper(*, dP: float, dH: float,
                          htheta_norm: float, a_row_norm: float) -> float:
    vals = (dP, dH, htheta_norm, a_row_norm)
    if not all(math.isfinite(float(x)) and float(x) >= 0.0 for x in vals):
        raise ValueError("finite nonnegative theta-x Delta-C inputs required")
    return _sum_up(
        FULL.up(dP * htheta_norm),
        FULL.up(a_row_norm * dH),
        FULL.up(dP * dH),
        dP,
    )


def _exact_theta_x_gain_detail(path: Path, *, source_pieces: int,
                               source_cell_index: int, p_pieces: int,
                               base: dict, vr: dict) -> dict:
    V12D = V31.V23.V22.V21B.V21.V12D
    V11 = V12D.V11
    V10 = V11.V10
    dom = json.loads(path.read_text(encoding="utf-8"))
    src, phase = V10.RG._source_phase_children(source_pieces)[source_cell_index]
    if phase != "due":
        raise RuntimeError("V32 focused refinement requires first due source cell")

    first = V11.FIRST.build(path, source_pieces=source_pieces)
    fr = first["source_cells"][source_cell_index]
    p_all = Interval.outward_bounds(*map(float, fr["P_aw_variance_interval"]))
    pcells = V11.SUB.parts(p_all.lo, p_all.hi, p_pieces)
    pi = int(base["p_cell"])
    p = pcells[pi]

    hstep = float(src["dt_s"])
    g = float(dom["startup"]["gravity_mps2"])
    tilt, _yaw, eps = V10.RG._attitude_covariance_epsilon(path, hstep)
    t = Interval.outward_bounds(tilt, FULL.up(tilt + eps))
    vec = V11.VECTOR.build()
    r = FULL._R_diag(float(vec["configured_measurement_bounds"]["acc_measurement_std_mps2"]))[0][0]
    D = FULL.I(g * g) * t + p + r
    a = t * (p + r) / D
    a_row = a.abs_upper()

    dP = float(vr["total_reduced_covariance_perturbation_upper"])
    dH = float(vr["sample1_H_perturbation_upper"])
    htheta = float(vr["sample1_Htheta_operator_upper"])
    dS = float(vr["sample1_innovation_perturbation_upper"])
    inv = float(vr["actual_innovation_inverse_operator_upper"])
    kpar = float(base["Ktheta_parallel_block_upper"])
    parent_dC = float(vr["sample1_attitude_cross_covariance_perturbation_upper"])
    parent_dk = float(vr["sample1_attitude_gain_operator_perturbation_upper"])

    dC_x = _theta_x_deltac_upper(
        dP=dP, dH=dH, htheta_norm=htheta, a_row_norm=a_row)
    if dC_x > FULL.up(parent_dC):
        raise RuntimeError("exact theta-x Delta-C exceeded V12D parent")
    numerator = FULL.up(dC_x + FULL.up(kpar * dS))
    dk_x = FULL.up(numerator * inv)
    if dk_x > FULL.up(parent_dk):
        raise RuntimeError("exact theta-x gain perturbation exceeded V12D parent")

    return {
        "theta_x_nominal_attitude_covariance_row_norm_upper": a_row,
        "theta_x_DeltaC_operator_upper": dC_x,
        "V12D_full_DeltaC_operator_upper": parent_dC,
        "sample1_attitude_cross_covariance_perturbation_upper": dC_x,
        "sample1_innovation_perturbation_upper": dS,
        "actual_innovation_inverse_operator_upper": inv,
        "nominal_theta_x_gain_row_norm_upper": kpar,
        "theta_x_gain_perturbation_operator_upper": dk_x,
        "V12D_full_attitude_gain_perturbation_operator_upper": parent_dk,
        "theta_x_covariance_row_exactly_a_0_0": True,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24,
          residual_x_pieces: int = 6, parallel_pieces: int = 6,
          current_component_pieces: int = 4) -> dict:
    path = Path(domain_path).resolve()
    V12D = V31.V23.V22.V21B.V21.V12D
    V10 = V12D.V11.V10
    v12 = V12D.build(path, source_pieces=source_pieces,
                     source_cell_index=source_cell_index, p_pieces=p_pieces,
                     tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    core = V10.build(path, source_pieces=source_pieces,
                     source_cell_index=source_cell_index, p_pieces=p_pieces,
                     tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
    failures = [f"V12D: {x}" for x in V12D.validate(v12)]
    failures += [f"V10: {x}" for x in V10.validate(core)]
    try:
        vr = V31.V30._witness_row(v12)
        base = V31.V30._witness_row(core)
        detail = _exact_theta_x_gain_detail(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            base=base, vr=vr)
        old = V31.V30._theta_x_gain_perturbation_upper(vr, base)
        if float(detail["theta_x_gain_perturbation_operator_upper"]) > FULL.up(
                float(old["theta_x_gain_perturbation_operator_upper"])):
            raise RuntimeError("V32 theta-x gain refinement exceeded V30 parent")
    except Exception as exc:
        failures.append(f"V32 theta-x Delta-C refinement: {exc}")
        detail = None

    original = V31.V30._theta_x_gain_perturbation_upper
    def exact_row(vr_row: dict, base_row: dict) -> dict:
        if detail is None:
            return original(vr_row, base_row)
        return dict(detail)
    V31.V30._theta_x_gain_perturbation_upper = exact_row
    try:
        parent = V31.build(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            tangent_pieces=tangent_pieces, axial_pieces=axial_pieces,
            residual_x_pieces=residual_x_pieces,
            parallel_pieces=parallel_pieces,
            current_component_pieces=current_component_pieces)
    finally:
        V31.V30._theta_x_gain_perturbation_upper = original

    failures += [f"V31: {x}" for x in V31.validate(parent)]
    if parent.get("P5_SAMPLE1_V30_CURRENT_SUBBOX_LIFT_V31") != "PASS":
        failures.append("V31 current-subbox lift prerequisite did not pass")

    out = dict(parent)
    out.update({
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_EXACT_THETA_X_DELTAC_LIFT_V32",
        "V31_current_subbox_lift_parent_retained": True,
        "V30_theta_x_row_resolvent_retained": True,
        "V12D_full_DeltaC_parent_retained": True,
        "theta_x_nominal_covariance_row_exact_a_0_0_used": True,
        "theta_x_Ptheta_DeltaH_uses_exact_a_row_norm": True,
        "theta_x_exact_DeltaC_gain_detail": detail,
        "deployed_correction_limit_rad": 6.0,
        "deployed_correction_limit_increased": False,
        "q_target": Q_TARGET,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_EXACT_THETA_X_DELTAC_LIFT_V32": (
            "PASS" if not failures else "NOT_ESTABLISHED"),
        "next_obligation": (
            "LIFT_V32_ROW_DIRECTIONAL_REFINEMENT_INTO_FULL_SOURCE_CELL0_Q8_COVER"
            if parent.get("focused_first_witness_signed_subcell_closed_by_V30_lift") is True
            else "REFINE_FIRST_REMAINING_V32_CURRENT_SUBBOX_WITH_DIRECTIONAL_DELTA_S_OR_YZ_GAIN_ROWS"),
        "failures": list(dict.fromkeys(failures)),
    })
    return out


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_EXACT_THETA_X_DELTAC_LIFT_V32":
        f.append("qualification mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "V31_current_subbox_lift_parent_retained",
        "V30_theta_x_row_resolvent_retained",
        "V12D_full_DeltaC_parent_retained",
        "theta_x_nominal_covariance_row_exact_a_0_0_used",
        "theta_x_Ptheta_DeltaH_uses_exact_a_row_norm",
        "V22_exact_current_residual_parent_retained",
        "V28_split_gravity_signed_source_enclosure_retained",
        "V29_yz_and_radial_perturbation_parents_retained",
        "V23_current_partition_and_q_ball_projection_retained",
        "current_dependent_and_source_directional_correction_enclosures_intersected",
        "V16_axis_cone_V15_geodesic_V18_yz_support_retained",
        "all_candidate_current_subboxes_accounted",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "deployed_correction_limit_increased",
        "q8_composed_here", "q8_word_promoted_here", "whole_word_promoted_here",
        "N_H_words_set_here",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    gd = d.get("theta_x_exact_DeltaC_gain_detail") or {}
    dc = float(gd.get("theta_x_DeltaC_operator_upper", math.inf))
    dcp = float(gd.get("V12D_full_DeltaC_operator_upper", -math.inf))
    dk = float(gd.get("theta_x_gain_perturbation_operator_upper", math.inf))
    dkp = float(gd.get("V12D_full_attitude_gain_perturbation_operator_upper", -math.inf))
    if not (math.isfinite(dc) and 0.0 <= dc <= FULL.up(dcp)):
        f.append("invalid exact theta-x Delta-C refinement")
    if not (math.isfinite(dk) and 0.0 <= dk <= FULL.up(dkp)):
        f.append("invalid exact theta-x gain refinement")
    if gd.get("theta_x_covariance_row_exactly_a_0_0") is not True:
        f.append("theta-x exact covariance row flag missing")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        f.append("deployed correction limit changed")
    if float(d.get("q_target", 0.0)) != Q_TARGET:
        f.append("q target changed")
    if d.get("P5_SAMPLE1_EXACT_THETA_X_DELTAC_LIFT_V32") not in ("PASS", "NOT_ESTABLISHED"):
        f.append("invalid V32 status")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=4)
    ap.add_argument("--source-cell-index", type=int, default=0)
    ap.add_argument("--p-pieces", type=int, default=24)
    ap.add_argument("--tangent-pieces", type=int, default=24)
    ap.add_argument("--axial-pieces", type=int, default=24)
    ap.add_argument("--residual-x-pieces", type=int, default=6)
    ap.add_argument("--parallel-pieces", type=int, default=6)
    ap.add_argument("--current-component-pieces", type=int, default=4)
    ap.add_argument("--output", type=Path, required=True)
    x = ap.parse_args()
    d = build(
        x.domain, source_pieces=x.source_pieces,
        source_cell_index=x.source_cell_index, p_pieces=x.p_pieces,
        tangent_pieces=x.tangent_pieces, axial_pieces=x.axial_pieces,
        residual_x_pieces=x.residual_x_pieces,
        parallel_pieces=x.parallel_pieces,
        current_component_pieces=x.current_component_pieces)
    vf = validate(d); d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_EXACT_THETA_X_DELTAC_LIFT_V32"],
        "theta_x_detail": d.get("theta_x_exact_DeltaC_gain_detail"),
        "candidate": d.get("candidate_current_subboxes"),
        "closed": d.get("closed_current_subboxes"),
        "open": d.get("open_current_subboxes"),
        "minimum_best_q": d.get("minimum_best_q_upper"),
        "maximum_best_q": d.get("maximum_best_q_upper"),
        "first_open": d.get("first_open_current_subbox"),
        "worst_open": d.get("worst_open_current_subbox"),
        "witness_closed": d.get("focused_first_witness_signed_subcell_closed_by_V30_lift"),
        "next": d.get("next_obligation"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
