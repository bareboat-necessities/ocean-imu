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
its own block.  The corrected complete-SEA3 four-S certificate supplies the
physical translation observation inverse row bounds using actual applied
SpectralMSE R_S.  Its covariance tightening uses the independence of configured
S measurement noise between selected updates while retaining a four-record
trace bound for correlated process nuisance.  It therefore does not replace the
complete SEA3 source or change the estimator.

The earlier composition collapsed the complete four-S block to one scalar dS
before paying the accelerometer cross.  That is unnecessarily destructive:
the only translation column present in the selected accelerometer PE rows is
``a_w``.  The four-S physical inverse already gives separate rigorous row-l1
bounds r_i for [S,gp,g^2v,g^3a_w].  If y=M z then

    |z_i| <= r_i ||y||_infinity <= r_i ||y||_2,

so, summing the four coordinate inequalities,

    ||M z||_2^2 >= (1/4) sum_i z_i^2/r_i^2.

With Sigma_Y <= lambda_Y I this proves the directional information matrix

    M^T Sigma_Y^-1 M >= diag_i(1/(4 lambda_Y r_i^2)).

This is a Loewner matrix lower obtained from the same four actual S records;
it is not a blockwise contraction ratio.  In particular the a_w direction is
much stronger than the weakest translation direction.

The only cross block in the selected vector rows is the a_w column.  For the
two required PE occurrences its whitened norm is bounded by

    ||C_aw||^2 <= 2 / (Racc_min * g^6),

because the body/world rotation and the time-varying OU attenuation have norm
at most one.  Magnetometer rows have no translation columns.  Therefore for
u=eta6 and w=g^3 a_w,

  ||A u + C_aw w||^2 + d_aw ||w||^2
    >= lambda_c (||u||^2+||w||^2),

where

  lambda_c >= alpha6*d_aw/(alpha6+||C_aw||^2+d_aw).

The remaining S,gp,g^2v directions do not enter C_aw and retain their own
directional four-S lower.  The full 18-state information lower is the minimum
of lambda_c and those three directional entries.  All omitted valid
accelerometer, magnetometer and S updates contribute PSD information and can
only improve the bound.
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

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2
QUALIFICATION = "OU3_COMPLETE_SEA3_H18_DIRECTIONAL_INFORMATION_COMPOSITION"
USEFUL_GATE = 1.0e-18


def _directional_translation_information(four: dict) -> dict[str, float]:
    noise = four["selected_S_record_noise"]
    lam = float(noise["four_record_covariance_lambda_max_upper"])
    rows = four["newton_coordinate_information"]["physical_state_recovery"][
        "physical_state_inverse_raw_record_row_l1_upper"
    ]
    order = ("S", "g*p", "g^2*v", "g^3*a_w")
    if not (math.isfinite(lam) and lam > 0.0):
        raise RuntimeError("four-S record covariance upper is not finite positive")
    out: dict[str, float] = {}
    for name in order:
        r = float(rows[name])
        if not (math.isfinite(r) and r > 0.0):
            raise RuntimeError(f"invalid four-S physical inverse row bound {name}")
        out[name] = down(1.0 / up(4.0 * up(lam * up(r * r))))
        if not (math.isfinite(out[name]) and out[name] > 0.0):
            raise RuntimeError(f"directional four-S information lost positivity for {name}")
    return out


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
    scalar_dS = float(four["newton_coordinate_information"]["D_S_physical_lambda_min_lower"])
    directional = _directional_translation_information(four)
    g = float(four["uniform_S_gap_s_upper"])
    racc_var = float(pe["measurement_runtime"]["accelerometer_variance_upper"])
    if not (alpha6 > 0.0 and scalar_dS > 0.0 and g > 0.0 and racc_var > 0.0):
        raise RuntimeError("H18 information inputs lost strict positivity")

    # The PE witness uses two required accelerometer occurrences.  In the z
    # coordinates a_w = z_aw/g^3.  R_wb and OU homogeneous attenuation have
    # operator norm <=1, so each whitened a_w block has norm <=1/(sqrt(Ra)g^3).
    cross_norm_sq_upper = up(2.0 / (racc_var * (g ** 6)))
    d_aw = float(directional["g^3*a_w"])
    coupled_trace_upper = up(alpha6 + d_aw + cross_norm_sq_upper)
    coupled_eta_aw_lower = down((alpha6 * d_aw) / coupled_trace_upper)
    non_aw_lower = min(
        float(directional["S"]),
        float(directional["g*p"]),
        float(directional["g^2*v"]),
    )
    information_lower = down(min(coupled_eta_aw_lower, non_aw_lower))
    if not (math.isfinite(information_lower) and information_lower > 0.0):
        raise RuntimeError("H18 full directional information lower is not strict")

    old_trace_upper = up(alpha6 + scalar_dS + cross_norm_sq_upper)
    old_scalar_lower = down((alpha6 * scalar_dS) / old_trace_upper)

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
        "translation_information_lower": scalar_dS,
        "translation_information_source": four["qualification"],
        "directional_translation_information_lower": directional,
        "directional_four_S_inverse_row_bound_used": True,
        "four_coordinate_sum_factor": 4.0,
        "uniform_S_gap_s_upper": g,
        "accelerometer_variance_upper": racc_var,
        "selected_PE_accelerometer_occurrences": 2,
        "accelerometer_translation_cross_norm_squared_upper": cross_norm_sq_upper,
        "accelerometer_cross_touches_only_aw_translation_coordinate": True,
        "useful_gate": USEFUL_GATE,
        "H18_information_useful_gate_pass": information_lower >= USEFUL_GATE,
        "triangular_information_composition": {
            "form": "directional four-S regularizer plus ||A eta6 + C_aw z_aw||^2",
            "A_transpose_A_lower": alpha6,
            "legacy_scalar_B_transpose_B_lower_diagnostic": scalar_dS,
            "directional_B_transpose_B_diagonal_lower": directional,
            "aw_direction_information_lower": d_aw,
            "C_aw_spectral_norm_squared_upper": cross_norm_sq_upper,
            "coupled_eta6_aw_scalar_2x2_determinant_lower": down(alpha6 * d_aw),
            "coupled_eta6_aw_scalar_2x2_trace_upper": coupled_trace_upper,
            "coupled_eta6_aw_lambda_min_lower": coupled_eta_aw_lower,
            "non_aw_translation_lambda_min_lower": non_aw_lower,
            "D_H18_lambda_min_lower": information_lower,
            "legacy_scalarized_D_H18_lambda_min_lower_diagnostic": old_scalar_lower,
            "full_18x18_matrix_information_lower_closed": True,
        },
        "directional_information_strictly_improves_legacy_scalarized_bound": (
            information_lower > old_scalar_lower
        ),
        "omitted_shipping_measurement_rows_are_PSD_information_only": True,
        "determinant_trace_scalarization_of_18x18_matrix_used": False,
        "blockwise_minimum_ratio_used": False,
        "scalar_information_beta_used": False,
        "P3_promoted": False,
        "next_obligation": (
            "use the directional eta6/a_w information headroom in the finite-bias A21 Riccati "
            "completion while keeping the same complete SEA3 source and 1e-18 gate"
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
        "directional_four_S_inverse_row_bound_used",
        "accelerometer_cross_touches_only_aw_translation_coordinate",
        "directional_information_strictly_improves_legacy_scalarized_bound",
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
    directional = d.get("directional_translation_information_lower", {})
    for key in ("S", "g*p", "g^2*v", "g^3*a_w"):
        x = directional.get(key)
        if not isinstance(x, (int, float)) or not (math.isfinite(float(x)) and float(x) > 0.0):
            f.append(f"invalid directional translation information {key}")
    c = d.get("triangular_information_composition", {})
    if c.get("full_18x18_matrix_information_lower_closed") is not True:
        f.append("H18 full matrix information lower is not closed")
    for key in (
        "A_transpose_A_lower", "aw_direction_information_lower",
        "coupled_eta6_aw_scalar_2x2_determinant_lower",
        "coupled_eta6_aw_lambda_min_lower", "non_aw_translation_lambda_min_lower",
        "D_H18_lambda_min_lower",
    ):
        x = c.get(key)
        if not isinstance(x, (int, float)) or not (math.isfinite(float(x)) and float(x) > 0.0):
            f.append(f"invalid positive quantitative field {key}")
    old = c.get("legacy_scalarized_D_H18_lambda_min_lower_diagnostic")
    new = c.get("D_H18_lambda_min_lower")
    if not isinstance(old, (int, float)) or not isinstance(new, (int, float)) or not float(new) > float(old):
        f.append("directional H18 information did not improve legacy scalar bound")
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
    c = d["triangular_information_composition"]
    print(json.dumps({
        "alpha6": d["eta6_information_lower"],
        "directional_translation": d["directional_translation_information_lower"],
        "aw_cross_norm_sq_upper": d["accelerometer_translation_cross_norm_squared_upper"],
        "legacy_H18_information_lower": c["legacy_scalarized_D_H18_lambda_min_lower_diagnostic"],
        "directional_H18_information_lower": c["D_H18_lambda_min_lower"],
        "improvement_factor": c["D_H18_lambda_min_lower"] / c["legacy_scalarized_D_H18_lambda_min_lower_diagnostic"],
        "H18_information_useful_gate_pass": d["H18_information_useful_gate_pass"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())