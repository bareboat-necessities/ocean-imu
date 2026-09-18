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
block on [theta_h,S_h]

    P2 = [[p_theta, c], [c, p_S]],
    p_S = R_h,
    c = rho sqrt(p_theta p_S),       0 < rho < 1.

P2 is strictly positive definite because det(P2)=p_theta p_S(1-rho^2)>0.
It can be embedded in a full positive diagonal covariance while respecting all
reported Pbar diagonal ceilings.  The S=0 gain from this admissible marginal
block is

    K_theta = c/(p_S+R_h) = rho/2 sqrt(p_theta/R_h).

Taking a hard-entry residual delta S_h=r_S therefore gives

    |d_theta| = r_S rho/2 sqrt(p_theta/R_h).

The executable witness below uses p_theta=Pbar_theta/2 and rho=1/2, so it stays
strictly inside the covariance ceiling and away from a singular PSD boundary.
If this already lies outside the exact-reset utility chart, then no theorem that
uses only these detached marginals can establish reset admission.  A successful
P4 proof must retain additional reachable-history/same-cell information (for
example the actual correlated P/H/R/K event family or an equivalent stronger
reachable covariance invariant).

The deployed horizontal channels are anisotropic, R_h=(rho_h r_S)^2 with rho_h
read out of the shipping header per axis.  The correction is not monotone the
way the innovation magnitudes are: it goes as 1/rho_h, so a *smaller* deployed
factor makes the obstruction easier and the axis with the *largest* factor is
the hardest one.  Both axes are therefore carried and each must leave the reset
utility chart, with the hardest axis reported as the binding one; the x branch
is no longer binding by construction.

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
import ou3_brmm_complete_source as SOURCE
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_exact_reset_transport as RESET

QUALIFICATION = "OU3_P4_MARGINAL_PBAR_CORRECTION_OBSTRUCTION_V1"
RHO = 0.5
PBAR_FRACTION = 0.5


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def deployed_horizontal_rs_factors() -> list[float]:
    """Deployed [rho_x, rho_y]; fail closed on an unreadable pair."""
    horizontal, reason = SOURCE.deployed_horizontal_rs_factors()
    if horizontal is None:
        raise RuntimeError("deployed R_S horizontal factors unreadable: " + reason)
    return horizontal


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
    horizontal = deployed_horizontal_rs_factors()
    rS = float(entry["coordinate_radii"]["integral_displacement_norm_m_s"])

    modes = {}
    for mode, key in (("H18", "H"), ("A21", "A")):
        pbar = list(map(float, tube["modes"][key]["Pbar_diagonal_variance_upper"]))
        ptheta_bar = pbar[0]
        pS_bar = pbar[12]
        ptheta = down(PBAR_FRACTION * ptheta_bar)
        axes = {}
        for name, factor in zip(("x", "y"), horizontal, strict=True):
            rs_std = down(factor * rs_lo)
            if not rs_std > 0.0:
                raise RuntimeError(mode + " horizontal applied R_S lower lost positivity on " + name)
            R = down(rs_std * rs_std)
            pS = R
            if not (ptheta > 0.0 and pS > 0.0 and pS < pS_bar and ptheta < ptheta_bar):
                raise RuntimeError(mode + " obstruction witness does not lie strictly inside diagonal Pbar on " + name)
            c = down(RHO * math.sqrt(down(ptheta * pS)))
            det = down(ptheta * pS - c * c)
            if not det > 0.0:
                raise RuntimeError(mode + " two-state covariance witness lost strict positive definiteness on " + name)
            gain = down(c / up(pS + R))
            correction = down(rS * gain)
            axes[name] = {
                "deployed_RS_std_factor": factor,
                "witness_S_variance": pS,
                "witness_theta_S_covariance": c,
                "witness_2x2_determinant_lower": det,
                "witness_strictly_inside_diagonal_Pbar": True,
                "S_zero_horizontal_R_variance": R,
                "same_cell_S_zero_attitude_gain": gain,
                "same_cell_attitude_correction_norm": correction,
                "outside_exact_reset_utility_domain": correction > RESET.CAYLEY_MONOTONE_NORM_MAX,
            }
        # The correction goes as 1/rho_h, so the smallest correction comes from
        # the largest deployed factor: that axis is the hardest and binds the
        # obstruction.  Requiring it keeps the statement true on both axes.
        binding = min(axes, key=lambda n: axes[n]["same_cell_attitude_correction_norm"])
        bound = axes[binding]
        modes[mode] = {
            "Pbar_theta_variance_upper": ptheta_bar,
            "Pbar_S_variance_upper": pS_bar,
            "witness_theta_variance": ptheta,
            "witness_correlation_coefficient": RHO,
            "hard_S_residual_norm_used_on_one_axis": rS,
            "horizontal_axes": axes,
            "binding_horizontal_axis": binding,
            "witness_S_variance": bound["witness_S_variance"],
            "witness_theta_S_covariance": bound["witness_theta_S_covariance"],
            "witness_2x2_determinant_lower": bound["witness_2x2_determinant_lower"],
            "witness_strictly_inside_diagonal_Pbar": all(
                a["witness_strictly_inside_diagonal_Pbar"] for a in axes.values()
            ),
            "S_zero_horizontal_R_variance": bound["S_zero_horizontal_R_variance"],
            "same_cell_S_zero_attitude_gain": bound["same_cell_S_zero_attitude_gain"],
            "same_cell_attitude_correction_norm": bound["same_cell_attitude_correction_norm"],
            "outside_exact_reset_utility_domain": all(
                a["outside_exact_reset_utility_domain"] for a in axes.values()
            ),
        }

    obstruction = all(m["outside_exact_reset_utility_domain"] for m in modes.values())
    return {
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "statement_scope": "INSUFFICIENCY_OF_DETACHED_MARGINAL_PBAR_PLUS_HARD_ERROR_BOUNDS",
        "reachable_shipping_covariance_counterexample_claimed": False,
        "canonical_P4_falsified": False,
        "endpoint_covariance_diagonal_ceiling_consumed": True,
        "actual_applied_RS_horizontal_factors": horizontal,
        "RS_axis_factors_read_from_deployed_source": True,
        "every_deployed_horizontal_axis_obstructed": obstruction,
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
        "RS_axis_factors_read_from_deployed_source",
        "every_deployed_horizontal_axis_obstructed",
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
    horizontal, reason = SOURCE.deployed_horizontal_rs_factors()
    if horizontal is None:
        f.append("deployed R_S horizontal factors unreadable: " + reason)
    elif [float(x) for x in d.get("actual_applied_RS_horizontal_factors", [])] != horizontal:
        f.append("R_S horizontal pair no longer matches the deployed source")
    for mode, m in d.get("modes", {}).items():
        axes = m.get("horizontal_axes", {})
        if horizontal is not None and [
            float(a.get("deployed_RS_std_factor", 0.0)) for a in axes.values()
        ] != horizontal:
            f.append(mode + " witness axes are not the deployed horizontal pair")
        for name, a in axes.items():
            if a.get("witness_strictly_inside_diagonal_Pbar") is not True:
                f.append(mode + " " + name + " witness not inside Pbar")
            if not float(a.get("witness_2x2_determinant_lower", 0.0)) > 0.0:
                f.append(mode + " " + name + " witness is not strictly SPD")
            if a.get("outside_exact_reset_utility_domain") is not True:
                f.append(mode + " " + name + " marginal witness did not leave reset utility domain")
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
