#!/usr/bin/env python3
"""V49: authoritative theta-y/z Delta-C term decomposition after V48.

V48 proves that splitting the duplicated y/z correction radius is real but does
not close the authoritative V41 survivor: the y cap shrinks while the z cap is
still pinned at the V31 parent.  V33 previously attempted rowwise theta-y/z
Delta-C bounds but correctly failed closed when its scalar row triangle bound
exceeded the already-certified V12D full-operator parent.

V49 does not promote that failed candidate.  It reconstructs the authoritative
V40/V45 witness and decomposes the V33 row candidate

    dC_i <= dP ||H_theta|| + ||P_theta[i,:]|| dH + dP dH + dP

for i=y,z into its four separately outward-rounded terms.  It also combines
those candidate rows with V34's certified first-row Delta-S refinement to show
which numerator term dominates the y/z gain-resolvent budget.  This is a
source-derived diagnostic used to choose the next component-matrix proof; V12D
remains the binding parent and no q<8, sample-1, whole-word, P5, or N_H
promotion occurs here.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_sample1_authoritative_componentwise_yz_v48 as V48
import ou3_p5_sample1_exact_theta_yz_gain_rows_lift_v33 as V33

DEFAULT_DOMAIN = V48.DEFAULT_DOMAIN
SCHEMA = 4900
WITNESS = V48.WITNESS
FULL = V48.FULL


def _sum_up(*xs: float) -> float:
    y = 0.0
    for x in xs:
        y = FULL.up(y + float(x))
    return y


def _deltac_terms(*, dP: float, dH: float, htheta: float,
                  row_norm: float) -> dict:
    vals = (dP, dH, htheta, row_norm)
    if not all(math.isfinite(float(x)) and float(x) >= 0.0 for x in vals):
        raise ValueError("finite nonnegative Delta-C term inputs required")
    projected_dP_H = FULL.up(dP * htheta)
    nominal_row_dH = FULL.up(row_norm * dH)
    mixed_dP_dH = FULL.up(dP * dH)
    theta_aw_cross = FULL.up(dP)
    total = _sum_up(projected_dP_H, nominal_row_dH, mixed_dP_dH,
                    theta_aw_cross)
    terms = {
        "projected_DeltaP_Htheta_upper": projected_dP_H,
        "nominal_Ptheta_row_DeltaH_upper": nominal_row_dH,
        "mixed_DeltaP_DeltaH_upper": mixed_dP_dH,
        "theta_aw_cross_block_parent_upper": theta_aw_cross,
        "row_DeltaC_candidate_upper": total,
    }
    dominant = max(
        (k for k in terms if k != "row_DeltaC_candidate_upper"),
        key=lambda k: terms[k])
    terms["dominant_term"] = dominant
    terms["dominant_fraction"] = 0.0 if total == 0.0 else terms[dominant] / total
    return terms


def _nominal_yz_rows(path: Path, *, source_pieces: int,
                     source_cell_index: int, p_pieces: int,
                     base: dict) -> tuple[float, float]:
    """Reconstruct V33's nominal theta-y/z covariance-row norms."""
    V12D = V48.V12D
    V11 = V12D.V11
    V10 = V11.V10
    dom = json.loads(path.read_text(encoding="utf-8"))
    src, phase = V10.RG._source_phase_children(source_pieces)[source_cell_index]
    if phase != "due":
        raise RuntimeError("V49 focused decomposition requires first due source cell")

    first = V11.FIRST.build(path, source_pieces=source_pieces)
    fr = first["source_cells"][source_cell_index]
    p_all = Interval.outward_bounds(*map(float, fr["P_aw_variance_interval"]))
    p = V11.SUB.parts(p_all.lo, p_all.hi, p_pieces)[int(base["p_cell"])]

    hstep = float(src["dt_s"])
    g = float(dom["startup"]["gravity_mps2"])
    tilt, yaw, eps = V10.RG._attitude_covariance_epsilon(path, hstep)
    t = Interval.outward_bounds(tilt, FULL.up(tilt + eps))
    Y = Interval.outward_bounds(yaw, FULL.up(yaw + eps))
    vec = V11.VECTOR.build()
    r = FULL._R_diag(float(
        vec["configured_measurement_bounds"]["acc_measurement_std_mps2"]))[0][0]
    F, Q, _ = FULL._transition_and_Q(src, dom)
    alpha = F[15][15]
    qaw = Q[15][15]

    rt = Interval.outward_bounds(
        *map(float, base["first_tangent_residual_magnitude_mps2"]))
    rz = Interval.outward_bounds(*map(float, base["first_axial_residual_mps2"]))
    d = Interval.outward_bounds(*map(float, base["first_attitude_correction_rad"]))
    D = FULL.I(g * g) * t + p + r
    fy = -(alpha * (p / D) * rt)
    fz = FULL.I(g) + alpha * (p / (p + r)) * rz
    Pn, _Hn, _Sn = V11._nominal_sample1_matrices(
        t=t, Y=Y, p=p, r=r, g=g, alpha=alpha, qaw=qaw,
        d=d, fy=fy, fz=fz)
    return V33._row_norm_upper(Pn[1][:3]), V33._row_norm_upper(Pn[2][:3])


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 4,
          source_cell_index: int = 0, p_pieces: int = 24,
          tangent_pieces: int = 24, axial_pieces: int = 24) -> dict:
    path = Path(domain_path).resolve()
    failures: list[str] = []
    try:
        _core, _v12, base, vr, row_failures = V48._build_v40_rows(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            tangent_pieces=tangent_pieces, axial_pieces=axial_pieces)
        failures += row_failures
        if (int(base.get("p_cell", -1)),
            int(base.get("tangent_residual_cell", -1)),
            int(base.get("axial_residual_cell", -1))) != WITNESS:
            failures.append("V49 did not reconstruct authoritative witness")

        py, pz = _nominal_yz_rows(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            base=base)
        ds = V48.V34._directional_delta_s_detail(
            path, source_pieces=source_pieces,
            source_cell_index=source_cell_index, p_pieces=p_pieces,
            base=base, vr=vr)

        dP = float(vr["total_reduced_covariance_perturbation_upper"])
        dH = float(vr["sample1_H_perturbation_upper"])
        htheta = float(vr["sample1_Htheta_operator_upper"])
        parent_dC = float(vr["sample1_attitude_cross_covariance_perturbation_upper"])
        parent_dK = float(vr["sample1_attitude_gain_operator_perturbation_upper"])
        inv = float(vr["actual_innovation_inverse_operator_upper"])
        dS0 = float(ds["first_measurement_row_DeltaS_intersected_upper"])
        ky = float(ds["nominal_theta_y_gain_row_norm_upper"])
        kz = float(ds["nominal_theta_z_gain_row_norm_upper"])

        y = _deltac_terms(dP=dP, dH=dH, htheta=htheta, row_norm=py)
        z = _deltac_terms(dP=dP, dH=dH, htheta=htheta, row_norm=pz)
        for row, kval in ((y, ky), (z, kz)):
            dc = float(row["row_DeltaC_candidate_upper"])
            ds_term = FULL.up(kval * dS0)
            numerator = FULL.up(dc + ds_term)
            dk_candidate = FULL.up(numerator * inv)
            parent_numerator = FULL.up(parent_dC + ds_term)
            parent_based_dk = FULL.up(parent_numerator * inv)
            row.update({
                "V12D_full_DeltaC_parent_upper": parent_dC,
                "candidate_over_parent_ratio": (math.inf if parent_dC == 0.0 else dc / parent_dC),
                "candidate_within_V12D_parent": dc <= FULL.up(parent_dC),
                "directional_nominalK_DeltaS_term_upper": ds_term,
                "gain_resolvent_candidate_numerator_upper": numerator,
                "gain_resolvent_candidate_upper": dk_candidate,
                "gain_resolvent_using_V12D_DeltaC_parent_upper": parent_based_dk,
                "V12D_full_gain_perturbation_parent_upper": parent_dK,
            })

        dominant = z["dominant_term"]
        if dominant == "theta_aw_cross_block_parent_upper":
            next_obligation = "REFINE_THETA_Z_AW_CROSS_BLOCK_COMPONENT_OF_DELTAP_ON_AUTHORITATIVE_V40_PARENT"
        elif dominant == "projected_DeltaP_Htheta_upper":
            next_obligation = "REFINE_THETA_Z_PROJECTED_DELTAP_HTHETA_COMPONENT_MATRIX_ON_AUTHORITATIVE_V40_PARENT"
        elif dominant == "nominal_Ptheta_row_DeltaH_upper":
            next_obligation = "REFINE_THETA_Z_NOMINAL_ROW_DELTAH_COMPONENT_ON_AUTHORITATIVE_V40_PARENT"
        else:
            next_obligation = "REFINE_THETA_Z_MIXED_DELTAP_DELTAH_COMPONENT_ON_AUTHORITATIVE_V40_PARENT"
    except Exception as exc:
        failures.append(f"V49 authoritative Delta-C decomposition: {exc}")
        y = z = None
        next_obligation = "REPAIR_V49_AUTHORITATIVE_DELTAC_DECOMPOSITION"

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SAMPLE1_AUTHORITATIVE_YZ_DELTAC_DECOMPOSITION_V49",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "authoritative_V40_rows_reconstructed": y is not None and z is not None,
        "V12D_full_DeltaC_parent_retained": True,
        "V34_directional_DeltaS_retained": True,
        "failed_V33_row_candidate_promoted": False,
        "theta_y_DeltaC_term_decomposition": y,
        "theta_z_DeltaC_term_decomposition": z,
        "q8_composed_here": False,
        "q8_word_promoted_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_established_here": False,
        "P5_SAMPLE1_AUTHORITATIVE_YZ_DELTAC_DECOMPOSITION_V49": (
            "PASS" if not failures else "NOT_ESTABLISHED"),
        "next_obligation": next_obligation,
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P5_SAMPLE1_AUTHORITATIVE_YZ_DELTAC_DECOMPOSITION_V49":
        f.append("qualification mismatch")
    for k in ("source_generated_not_trajectory_fit",
              "authoritative_V40_rows_reconstructed",
              "V12D_full_DeltaC_parent_retained",
              "V34_directional_DeltaS_retained"):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in ("source_replay_used", "filter_changed",
              "failed_V33_row_candidate_promoted", "q8_composed_here",
              "q8_word_promoted_here", "whole_word_promoted_here",
              "N_H_words_set_here", "P5_established_here"):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    for label in ("theta_y_DeltaC_term_decomposition",
                  "theta_z_DeltaC_term_decomposition"):
        row = d.get(label) or {}
        dc = float(row.get("row_DeltaC_candidate_upper", math.inf))
        parent = float(row.get("V12D_full_DeltaC_parent_upper", -math.inf))
        if not (math.isfinite(dc) and dc >= 0.0 and math.isfinite(parent) and parent >= 0.0):
            f.append(f"invalid {label}")
    if d.get("P5_SAMPLE1_AUTHORITATIVE_YZ_DELTAC_DECOMPOSITION_V49") not in (
            "PASS", "NOT_ESTABLISHED"):
        f.append("invalid V49 status")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=4)
    ap.add_argument("--source-cell-index", type=int, default=0)
    ap.add_argument("--p-pieces", type=int, default=24)
    ap.add_argument("--tangent-pieces", type=int, default=24)
    ap.add_argument("--axial-pieces", type=int, default=24)
    ap.add_argument("--output", type=Path, required=True)
    x = ap.parse_args()
    d = build(x.domain, source_pieces=x.source_pieces,
              source_cell_index=x.source_cell_index, p_pieces=x.p_pieces,
              tangent_pieces=x.tangent_pieces, axial_pieces=x.axial_pieces)
    vf = validate(d)
    d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_SAMPLE1_AUTHORITATIVE_YZ_DELTAC_DECOMPOSITION_V49"],
        "theta_y": d.get("theta_y_DeltaC_term_decomposition"),
        "theta_z": d.get("theta_z_DeltaC_term_decomposition"),
        "next": d.get("next_obligation"),
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
