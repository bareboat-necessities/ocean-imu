#!/usr/bin/env python3
"""Non-promoting numerical diagnostic after correlated innovation repair.

This diagnostic intentionally stops before startup/capture composition.  It asks:
with the false rectangular-innovation singularity removed, what quantitative
post-Live H18 information certificate is obtained on the current COMPLETE-BRMM
word?  It prints the actual directional information lower bounds and their
ratio to the frozen P3/P4 useful gate 1e-18.

No value here promotes P3/P4/P5.  If the information lower is positive but below
1e-18, classify the result as quantitative enclosure/headroom insufficiency,
not literal rank deficiency and not an innovation singularity.
"""
from __future__ import annotations

import json
import math

import ou3_brmm_h18_information_composition as HINFO
import ou3_correlated_innovation_family as CORR

GATE = 1.0e-18


def build() -> dict:
    h = HINFO.build()
    hf = HINFO.validate(h)
    if hf:
        raise RuntimeError("H18 information diagnostic prerequisite invalid: " + repr(hf))
    c = h["triangular_information_composition"]
    lam = float(c["D_H18_lambda_min_lower"])
    directional = {k: float(v) for k, v in h["directional_translation_information_lower"].items()}
    regression = CORR.validate_point_regression()
    return {
        "qualification": "OU3_P4_POST_LIVE_CORRELATED_INNOVATION_NUMERICS_V1",
        "innovation_covariance_rectangular_singularity_removed": True,
        "same_P_H_R_used_through_S_inverse_K": bool(regression["same_P_H_R_used_for_PHt_S_Sinv_K"]),
        "H18_information_lambda_min_lower": lam,
        "H18_information_gate": GATE,
        "H18_information_to_gate_ratio": lam / GATE,
        "H18_information_strictly_positive": lam > 0.0,
        "H18_information_gate_pass": lam >= GATE,
        "eta6_information_lower": float(h["eta6_information_lower"]),
        "directional_translation_information_lower": directional,
        "aw_cross_norm_squared_upper": float(h["accelerometer_translation_cross_norm_squared_upper"]),
        "coupled_eta6_aw_lambda_min_lower": float(c["coupled_eta6_aw_lambda_min_lower"]),
        "non_aw_translation_lambda_min_lower": float(c["non_aw_translation_lambda_min_lower"]),
        "legacy_scalarized_H18_information_lower": float(c["legacy_scalarized_D_H18_lambda_min_lower_diagnostic"]),
        "failure_class_if_gate_missed": None if lam >= GATE else "C: quantitative enclosure/headroom below frozen 1e-18 gate",
        "literal_information_rank_loss_claimed": False,
        "innovation_singularity_claimed": False,
        "P3_delta": GATE,
        "P4_PASS": False,
        "P5_MAY_START": False,
    }


def validate(d: dict) -> list[str]:
    f=[]
    if d.get("innovation_covariance_rectangular_singularity_removed") is not True:
        f.append("fake innovation singularity returned")
    if d.get("same_P_H_R_used_through_S_inverse_K") is not True:
        f.append("correlated P/H/R provenance absent")
    lam=float(d.get("H18_information_lambda_min_lower", math.nan))
    if not (math.isfinite(lam) and lam > 0.0):
        f.append("H18 information lower is not strictly positive")
    if float(d.get("P3_delta", 0.0)) != GATE:
        f.append("frozen 1e-18 gate changed")
    if d.get("literal_information_rank_loss_claimed") is not False:
        f.append("positive-but-small information misclassified as rank loss")
    if d.get("innovation_singularity_claimed") is not False:
        f.append("old innovation singularity misclassified as current blocker")
    if d.get("P4_PASS") or d.get("P5_MAY_START"):
        f.append("diagnostic promoted P4/P5")
    return f


def main() -> int:
    d=build(); f=validate(d); d["validation_pass"]=not f; d["validation_failures"]=f
    print(json.dumps(d, indent=2, sort_keys=True))
    return int(bool(f))

if __name__ == "__main__":
    raise SystemExit(main())
