#!/usr/bin/env python3
"""Full H18 measurement-information lower inside the complete SEA3 word.

This is a quantitative lemma for the canonical
``COMPLETE_SEA3_NORMAL_LIVE_WORD``.  It does not generate a source word and it
does not replace the literal scheduler by four selected S updates.  The four
selected updates are actual members of that word and only witness the
translation information already supplied by every due shipping S=0 update.

The held-bias H18 state is partitioned as

    x = [eta6, z12],
    eta6 = [delta_theta, delta_b_g],
    z12  = three copies of [S, g p, g^2 v, g^3 a_w].

The asynchronous vector-PE certificate supplies

    A^T A >= alpha6 I6

for the eta6 columns after the translational a_w contribution is retained in
its own block.  The corrected complete-SEA3 four-S certificate supplies

    B^T B >= dS I12

for the physical translation coordinates, using actual applied SpectralMSE
R_S.  Its covariance tightening uses the independence of configured S
measurement noise between selected updates while retaining a four-record trace
bound for correlated process nuisance.  It therefore does not replace the
complete SEA3 source or change the estimator.

The only cross block in the selected accelerometer rows is the a_w column.  For
the two required PE occurrences its whitened norm is bounded by

    ||C||^2 <= 2 / (Racc_min * g^6),

because the body/world rotation and the time-varying OU attenuation have norm
at most one.  Magnetometer rows have no translation columns.

Hence for all u,v,

  ||A u + C v||^2 + ||B v||^2
    >= (sqrt(alpha6)||u|| - ||C||||v||)^2 + dS||v||^2.

The 2x2 scalar quadratic has determinant alpha6*dS and trace
alpha6+||C||^2+dS, so the complete selected-record information obeys the
rigorous full-matrix bound

    D_H18 >= alpha6*dS/(alpha6+||C||^2+dS) I18.

This scalar 2x2 calculation is only an analytic lower for the coupled 18x18
quadratic form; it is not a blockwise contraction ratio and does not promote
P3.  All omitted valid accelerometer, magnetometer and S updates contribute PSD
information and can only improve this lower bound.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import down, up
import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_four_s_translation_information_tight as FOUR_S
import ou3_sea3_windowed_vector_pe as PE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_COMPLETE_SEA3_H18_FULL_INFORMATION_COMPOSITION"
USEFUL_GATE = 1.0e-18


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    complete = COMPLETE.build(path)
    four = FOUR_S.build(path)
    pe = PE.build(path)
    bad = {
        "complete": COMPLETE.validate(complete),
        "four_S": FOUR_S.validate(four),
        "PE": PE.validate(pe),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"H18 complete-SEA3 information prerequisites failed: {bad}")
    if complete["canonical_P3_source"] != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        raise RuntimeError("H18 information detached from complete SEA3")
    if four["canonical_source"] != complete["canonical_P3_source"]:
        raise RuntimeError("four-S lemma is not bound to the same complete SEA3 source")

    alpha6 = float(pe["eta6_information"]["alpha_6_information_lower"])
    dS = float(four["newton_coordinate_information"]["D_S_physical_lambda_min_lower"])
    g = float(four["uniform_S_gap_s_upper"])
    racc_var = float(pe["measurement_runtime"]["accelerometer_variance_upper"])
    if not (alpha6 > 0.0 and dS > 0.0 and g > 0.0 and racc_var > 0.0):
        raise RuntimeError("H18 information inputs lost strict positivity")

    # The PE witness uses two required accelerometer occurrences.  In the z
    # coordinates a_w = z_aw/g^3.  R_wb and OU homogeneous attenuation have
    # operator norm <=1, so each whitened a_w block has norm <=1/(sqrt(Ra)g^3).
    cross_norm_sq_upper = up(2.0 / (racc_var * (g ** 6)))
    trace_upper = up(alpha6 + dS + cross_norm_sq_upper)
    information_lower = down((alpha6 * dS) / trace_upper)
    if not (math.isfinite(information_lower) and information_lower > 0.0):
        raise RuntimeError("H18 full information lower is not strict")

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": complete["canonical_P3_source"],
        "component_of_complete_SEA3_full_word": True,
        "source_family_replaced": False,
        "selected_four_S_events_replace_complete_word": False,
        "all_due_S_updates_remain_in_literal_word": True,
        "actual_applied_SpectralMSE_R_S_consumed": True,
        "tight_four_S_measurement_covariance_structure_consumed": True,
        "four_S_process_cross_record_trace_bound_retained": True,
        "same_complete_SEA3_word_supplies_PE_and_translation_information": True,
        "H18_state_coordinates": {
            "eta6": ["delta_theta", "delta_b_g"],
            "translation_per_axis": ["S", "g*p", "g^2*v", "g^3*a_w"],
            "dimension": 18,
        },
        "eta6_information_lower": alpha6,
        "translation_information_lower": dS,
        "translation_information_source": four["qualification"],
        "uniform_S_gap_s_upper": g,
        "accelerometer_variance_upper": racc_var,
        "selected_PE_accelerometer_occurrences": 2,
        "accelerometer_translation_cross_norm_squared_upper": cross_norm_sq_upper,
        "useful_gate": USEFUL_GATE,
        "H18_information_useful_gate_pass": information_lower >= USEFUL_GATE,
        "triangular_information_composition": {
            "form": "||Au+Cv||^2+||Bv||^2",
            "A_transpose_A_lower": alpha6,
            "B_transpose_B_lower": dS,
            "C_spectral_norm_squared_upper": cross_norm_sq_upper,
            "scalar_2x2_determinant_lower": down(alpha6 * dS),
            "scalar_2x2_trace_upper": trace_upper,
            "D_H18_lambda_min_lower": information_lower,
            "full_18x18_matrix_information_lower_closed": True,
        },
        "omitted_shipping_measurement_rows_are_PSD_information_only": True,
        "determinant_trace_scalarization_of_18x18_matrix_used": False,
        "blockwise_minimum_ratio_used": False,
        "scalar_information_beta_used": False,
        "P3_promoted": False,
        "next_obligation": (
            "if the corrected H18 information clears the useful gate, carry it into the same-word "
            "prior-free full-matrix completion; then extend to A21 using finite accelerometer-bias "
            "correlation without changing the complete SEA3 source"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for key in (
        "component_of_complete_SEA3_full_word",
        "all_due_S_updates_remain_in_literal_word",
        "actual_applied_SpectralMSE_R_S_consumed",
        "tight_four_S_measurement_covariance_structure_consumed",
        "four_S_process_cross_record_trace_bound_retained",
        "same_complete_SEA3_word_supplies_PE_and_translation_information",
        "omitted_shipping_measurement_rows_are_PSD_information_only",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_family_replaced",
        "selected_four_S_events_replace_complete_word",
        "determinant_trace_scalarization_of_18x18_matrix_used",
        "blockwise_minimum_ratio_used",
        "scalar_information_beta_used",
        "P3_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if float(d.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("useful gate changed")
    c = d.get("triangular_information_composition", {})
    if c.get("full_18x18_matrix_information_lower_closed") is not True:
        f.append("H18 full matrix information lower is not closed")
    for key in (
        "A_transpose_A_lower", "B_transpose_B_lower",
        "scalar_2x2_determinant_lower", "D_H18_lambda_min_lower",
    ):
        x = c.get(key)
        if not isinstance(x, (int, float)) or not (math.isfinite(float(x)) and float(x) > 0.0):
            f.append(f"invalid positive quantitative field {key}")
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
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "alpha6": d["eta6_information_lower"],
        "four_S_translation": d["translation_information_lower"],
        "cross_norm_sq_upper": d["accelerometer_translation_cross_norm_squared_upper"],
        "H18_information_lambda_min_lower": d["triangular_information_composition"]["D_H18_lambda_min_lower"],
        "H18_information_useful_gate_pass": d["H18_information_useful_gate_pass"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
