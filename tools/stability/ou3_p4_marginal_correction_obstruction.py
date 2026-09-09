#!/usr/bin/env python3
"""Obstruction certificate for marginal-Pbar reset-domain proofs in OU-III P4.

This is deliberately *not* a reachable-source counterexample.  It proves a
more limited but useful impossibility statement about proof architecture.

Suppose a reset-domain proof uses only

* the canonical endpoint-referenced diagonal covariance ceiling P_ii <= Pbar_i;
* the declared hard S-error ball ||delta S|| <= r_S; and
* the source-uniform applied R_S interval,

while discarding the same-history correlation that determines the actual
attitude/S covariance block.  Those marginals do not imply a finite useful
attitude-correction bound.

For one horizontal S=0 channel, select the two-dimensional principal covariance
block on [theta_x,S_x]

    P2 = [[p_theta, c], [c, p_S]],
    p_S = R_x,
    c = rho sqrt(p_theta p_S),       0 < rho < 1.

P2 is strictly positive definite because det(P2)=p_theta p_S(1-rho^2)>0.
It can be embedded in a full positive diagonal covariance while respecting all
reported Pbar diagonal ceilings.  The S=0 gain from this admissible marginal
block is

    K_theta = c/(p_S+R_x) = rho/2 sqrt(p_theta/R_x).

Taking a hard-entry residual delta S_x=r_S therefore gives

    |d_theta| = r_S rho/2 sqrt(p_theta/R_x).

The executable witness below uses p_theta=Pbar_theta/2 and rho=1/2, so it stays
strictly inside the covariance ceiling and away from a singular PSD boundary.
If this already lies outside the exact-reset utility chart, then no theorem that
uses only these detached marginals can establish reset admission.  A successful
P4 proof must retain additional reachable-history/same-cell information (for
example the actual correlated P/H/R/K event family or an equivalent stronger
reachable covariance invariant).

No claim is made that this witness covariance is reachable by the shipping
Riccati recursion.  The certificate therefore cannot falsify P4; it only rules
out the marginal-envelope proof route.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_brmm_riccati_tube_factored as TUBE
import ou3_brmm_dynamic_source_certificate as DYNAMIC
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_exact_reset_transport as RESET

QUALIFICATION = "OU3_P4_MARGINAL_PBAR_CORRECTION_OBSTRUCTION_V1"
RHO = 0.5
PBAR_FRACTION = 0.5


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def build() -> dict:
    tube = TUBE.build()
    dynamic = DYNAMIC.build()
    entry = ENTRY.build()
    bad = {
        "endpoint_covariance_envelope": TUBE.validate_covariance_ceiling(tube),
        "dynamic": DYNAMIC.validate(dynamic),
        "entry": ENTRY.validate(entry),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError("marginal correction obstruction prerequisites failed: " + repr(bad))

    rs_lo = float(dynamic["dynamic_invariant"]["R_S_applied"][0])
    rs_std_x = down(0.72 * rs_lo)
    if not rs_std_x > 0.0:
        raise RuntimeError("horizontal applied R_S lower lost positivity")
    R = down(rs_std_x * rs_std_x)
    rS = float(entry["coordinate_radii"]["integral_displacement_norm_m_s"])

    modes = {}
    for mode, key in (("H18", "H"), ("A21", "A")):
        pbar = list(map(float, tube["modes"][key]["Pbar_diagonal_variance_upper"]))
        ptheta_bar = pbar[0]
        pS_bar = pbar[12]
        ptheta = down(PBAR_FRACTION * ptheta_bar)
        pS = R
        if not (ptheta > 0.0 and pS > 0.0 and pS < pS_bar and ptheta < ptheta_bar):
            raise RuntimeError(mode + " obstruction witness does not lie strictly inside diagonal Pbar")
        c = down(RHO * math.sqrt(down(ptheta * pS)))
        det = down(ptheta * pS - c * c)
        if not det > 0.0:
            raise RuntimeError(mode + " two-state covariance witness lost strict positive definiteness")
        gain = down(c / up(pS + R))
        correction = down(rS * gain)
        modes[mode] = {
            "Pbar_theta_variance_upper": ptheta_bar,
            "Pbar_S_variance_upper": pS_bar,
            "witness_theta_variance": ptheta,
            "witness_S_variance": pS,
            "witness_theta_S_covariance": c,
            "witness_correlation_coefficient": RHO,
            "witness_2x2_determinant_lower": det,
            "witness_strictly_inside_diagonal_Pbar": True,
            "S_zero_horizontal_R_variance": R,
            "hard_S_residual_norm_used_on_one_axis": rS,
            "same_cell_S_zero_attitude_gain": gain,
            "same_cell_attitude_correction_norm": correction,
            "outside_exact_reset_utility_domain": correction > RESET.CAYLEY_MONOTONE_NORM_MAX,
        }

    obstruction = all(m["outside_exact_reset_utility_domain"] for m in modes.values())
    return {
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "statement_scope": "INSUFFICIENCY_OF_DETACHED_MARGINAL_PBAR_PLUS_HARD_ERROR_BOUNDS",
        "reachable_shipping_covariance_counterexample_claimed": False,
        "canonical_P4_falsified": False,
        "endpoint_covariance_diagonal_ceiling_consumed": True,
        "actual_applied_RS_horizontal_factor": 0.72,
        "hard_entry_S_ball_consumed": True,
        "witness_covariance_is_strict_SPD_principal_block": True,
        "witness_uses_only_marginal_information": True,
        "same_history_reachability_correlation_intentionally_not_imposed": True,
        "marginal_only_reset_domain_proof_architecture_obstructed": obstruction,
        "required_replacement": (
            "source-correlated reachable P/H/R/K event family or an equivalent stronger reachable covariance invariant"
        ),
        "rowwise_or_trace_Pbar_retuning_can_close_P4": False,
        "modes": modes,
        "P4_promoted_here": False,
    }


def validate(d: dict) -> list[str]:
    f = []
    if d.get("qualification") != QUALIFICATION:
        f.append("qualification mismatch")
    for k in (
        "endpoint_covariance_diagonal_ceiling_consumed",
        "hard_entry_S_ball_consumed",
        "witness_covariance_is_strict_SPD_principal_block",
        "witness_uses_only_marginal_information",
        "same_history_reachability_correlation_intentionally_not_imposed",
        "marginal_only_reset_domain_proof_architecture_obstructed",
    ):
        if d.get(k) is not True:
            f.append(k + " not true")
    for k in (
        "reachable_shipping_covariance_counterexample_claimed",
        "canonical_P4_falsified",
        "rowwise_or_trace_Pbar_retuning_can_close_P4",
        "P4_promoted_here",
    ):
        if d.get(k) is not False:
            f.append(k + " not false")
    if float(d.get("actual_applied_RS_horizontal_factor", 0.0)) != 0.72:
        f.append("R_S horizontal factor changed")
    for mode, m in d.get("modes", {}).items():
        if m.get("witness_strictly_inside_diagonal_Pbar") is not True:
            f.append(mode + " witness not inside Pbar")
        if not float(m.get("witness_2x2_determinant_lower", 0.0)) > 0.0:
            f.append(mode + " witness is not strictly SPD")
        if m.get("outside_exact_reset_utility_domain") is not True:
            f.append(mode + " marginal witness did not leave reset utility domain")
    return f


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build()
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "obstructed": d["marginal_only_reset_domain_proof_architecture_obstructed"],
        "H18_correction": d["modes"]["H18"]["same_cell_attitude_correction_norm"],
        "A21_correction": d["modes"]["A21"]["same_cell_attitude_correction_norm"],
        "reset_utility_max": RESET.CAYLEY_MONOTONE_NORM_MAX,
        "failures": f,
    }, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
