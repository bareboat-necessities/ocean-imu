#!/usr/bin/env python3
"""Global nonlinear budget inputs for the source-complete OU-III P4 word.

This producer is the bridge between the repaired linear P3 certificate and the
new 0.8-rad nonlinear sector.  It intentionally does *not* replay a fixed word
or accumulate interval-AD remainders.  Instead it converts the exact global
vector-remainder identities into diagonal homogeneous quadratic forms that a
source-complete signed-Joseph word inequality can subtract operation by
operation.

For one accepted magnetometer update,

    eta_m^T R_m^-1 eta_m
        <= s^2 m_max^2 / r_m * ||c||^2.

For one accepted accelerometer update, Young's inequality gives

    eta_a^T R_a^-1 eta_a
        <= A(eps) ||c||^2 + B(eps) ||a_w||^2,

with the exact A-mode accelerometer-bias coefficient equal to zero.  The
coefficients are emitted directly in H/A state order.  When normalized by a
diagonal covariance upper bound, the optimal Young parameter is available in
closed form: eps*=B0/A0 and the minimized max-coordinate charge is A0+B0.
That optimization is a proof algebra improvement, not a fitted parameter.

The remaining P4 obligation is now narrow and explicit: construct a rigorous
source-complete lower form for the *positive* residual-information terms
r^T S^-1 r over all accepted vector branches and show that, after subtracting
these emitted eta forms, the complete H/A word has a strict generalized margin.
Until that matrix inequality is certified, P4 is not promoted.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p3_source_uniform_certificate as P3
import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_p4_source_word_timing as TIMING
import ou3_p4_vector_remainder_sector as REM
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _positive(x, label: str) -> float:
    y = float(x)
    if not math.isfinite(y) or y <= 0.0:
        raise RuntimeError(f"{label} must be finite positive, got {x!r}")
    return y


def _penalty_diagonal(mode: str, theta_coeff: float, aw_coeff: float) -> list[float]:
    n = 18 if mode == "H" else 21
    out = [0.0] * n
    for i in range(3):
        out[i] = up(theta_coeff)
    for i in range(15, 18):
        out[i] = up(aw_coeff)
    return out


def _mode(mode: str, p3: dict, domain: dict, rem: dict, vector: dict) -> dict:
    m = p3["modes"][mode]
    sigma = list(m["Sigma_diagonal_upper"])
    n = 18 if mode == "H" else 21
    if len(sigma) != n:
        raise RuntimeError(f"{mode} P3 directional Sigma dimension mismatch")

    live = domain["normal_live"]
    vc = vector["configured_measurement_bounds"]
    fhi = _positive(live["specific_force_norm_upper_mps2"], "specific-force upper")
    mhi = _positive(live["magnetic_vector_norm_upper_uT"], "magnetic upper")
    ra = _positive(vc["acc_measurement_variance_upper"], "accelerometer variance")
    rm = _positive(vc["mag_measurement_variance_upper"], "magnetometer variance")

    mag_geom = _positive(
        rem["mag_eta_squared_over_linear_rotation_squared_upper"],
        "mag remainder coefficient",
    )
    acc_att_geom = _positive(
        rem["acc_eta_force_rotation_quadratic_coefficient_upper"],
        "accelerometer attitude remainder coefficient",
    )
    acc_aw_geom = _positive(
        rem["acc_eta_aw_quadratic_coefficient_upper"],
        "accelerometer aw remainder coefficient",
    )

    mag_theta = up(mag_geom * up(fhi * 0.0 + mhi * mhi) / rm)
    acc_theta = up(acc_att_geom * up(fhi * fhi) / ra)
    acc_aw = up(acc_aw_geom / ra)

    # Normalize the diagonal eta form by the retained P3 directional Sigma
    # upper.  This is useful for choosing Young epsilon and sizing the later
    # signed lower form.  It is not, by itself, a word contraction statement.
    ptheta = up(max(float(x) for x in sigma[0:3]))
    paw = up(max(float(x) for x in sigma[15:18]))
    mag_metric_charge = up(mag_theta * ptheta)
    acc_att_metric_charge = up(acc_theta * ptheta)
    acc_aw_metric_charge = up(acc_aw * paw)

    # For generic Young epsilon, normalized terms have the form
    # A0(1+eps), B0(1+1/eps).  Their minimax optimum is eps=B0/A0 and the
    # common upper value is A0+B0.  The current REM certificate was emitted at
    # eps=1, so also provide the geometry-independent optimized coefficients
    # directly from s^2.
    s2 = _positive(rem["sin_half_angle_squared_upper"], "sin^2 half-angle")
    A0 = up(s2 * up(fhi * fhi) * ptheta / ra)
    B0 = up(up(4.0 * s2) * paw / ra)
    eps_star = up(B0 / A0) if A0 > 0.0 else math.inf
    acc_minimax_metric_charge = up(A0 + B0)
    acc_theta_opt = up((1.0 + eps_star) * s2 * up(fhi * fhi) / ra)
    acc_aw_opt = up((1.0 + 1.0 / eps_star) * up(4.0 * s2) / ra)

    return {
        "dimension": n,
        "P3_relative_Riccati_injection_margin_lower": m[
            "relative_Riccati_injection_margin_lower"
        ],
        "P3_linear_rho_upper": up(1.0 - float(m["relative_Riccati_injection_margin_lower"])),
        "P3_Sigma_diagonal_upper": sigma,
        "mag_eta_penalty_diagonal_per_accepted_update": _penalty_diagonal(
            mode, mag_theta, 0.0
        ),
        "acc_eta_penalty_diagonal_per_accepted_update_eps1": _penalty_diagonal(
            mode, acc_theta, acc_aw
        ),
        "accelerometer_bias_eta_coefficient": 0.0,
        "accelerometer_bias_cancels_exactly": True,
        "normalized_eps1": {
            "mag_theta_charge_upper": mag_metric_charge,
            "acc_theta_charge_upper": acc_att_metric_charge,
            "acc_aw_charge_upper": acc_aw_metric_charge,
        },
        "optimized_accelerometer_young": {
            "epsilon_star": eps_star,
            "attitude_coefficient_upper": acc_theta_opt,
            "aw_coefficient_upper": acc_aw_opt,
            "minimax_metric_charge_upper": acc_minimax_metric_charge,
            "derivation": "eps*=B0/A0; min max(A0(1+eps),B0(1+1/eps))=A0+B0",
            "trajectory_fitted": False,
        },
        "positive_residual_information_lower_form_built_here": False,
        "signed_word_generalized_margin_lower": None,
        "rho_full_nonlinear_word_upper": None,
        "P4_PROMOTED": False,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    p3 = P3.build(path)
    p3f = P3.validate(p3)
    cayley = CAYLEY.build(path)
    cf = CAYLEY.validate(cayley)
    rem = REM.build(path, outer_angle_rad=float(cayley["outer_angle_rad"]), young_epsilon=1.0)
    rf = REM.validate(rem)
    timing = TIMING.build(path)
    tf = TIMING.validate(timing)
    vector = VECTOR.build()
    vf = VECTOR.validate(vector)

    failures = [f"P3: {x}" for x in p3f]
    failures += [f"Cayley: {x}" for x in cf]
    failures += [f"remainder: {x}" for x in rf]
    failures += [f"timing: {x}" for x in tf]
    failures += [f"vector: {x}" for x in vf]

    modes = {}
    if not failures:
        for mode in ("H", "A"):
            modes[mode] = _mode(mode, p3, domain, rem, vector)

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_GLOBAL_SIGNED_JOSEPH_SECTOR_INPUTS",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "outer_angle_rad": cayley.get("outer_angle_rad"),
        "declared_filter_entrance_covered": cayley.get(
            "declared_filter_entrance_covered"
        ),
        "exact_vector_information_retention_factor_lower": cayley.get(
            "exact_vector_information_retention_factor_lower"
        ),
        "S_timing_consumed_by_linear_P3": timing.get(
            "S_timing_consumed_by_linear_P3_translation_UCO"
        ),
        "S_nonlinear_eta_identically_zero": timing.get("S_nonlinear_eta_identically_zero"),
        "modes": modes,
        "global_packet_count_times_defect_used_for_promotion": False,
        "fixed_terminal_schedule_used": False,
        "interval_AD_long_prefix_used": False,
        "signed_word_information_lower_form_remaining": True,
        "P3_ESTABLISHED": not p3f,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED": False,
        "P5_FINITE_CAPTURE_ESTABLISHED": False,
        "next_obligation": (
            "build the source-complete positive residual-information lower form for all accepted "
            "vector branches, subtract the emitted homogeneous eta forms in the same H/A metric, "
            "and certify a strict generalized margin mu>0 (rho=1-mu<1) over theta<=0.8 rad"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "declared_filter_entrance_covered",
        "S_timing_consumed_by_linear_P3",
        "S_nonlinear_eta_identically_zero",
        "signed_word_information_lower_form_remaining",
        "P3_ESTABLISHED",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used",
        "filter_changed",
        "global_packet_count_times_defect_used_for_promotion",
        "fixed_terminal_schedule_used",
        "interval_AD_long_prefix_used",
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED",
        "P5_FINITE_CAPTURE_ESTABLISHED",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    for mode in ("H", "A"):
        m = d.get("modes", {}).get(mode, {})
        if m.get("accelerometer_bias_cancels_exactly") is not True:
            failures.append(f"{mode} accelerometer-bias cancellation missing")
        q = m.get("optimized_accelerometer_young", {}).get("minimax_metric_charge_upper")
        if not isinstance(q, (int, float)) or not math.isfinite(float(q)) or float(q) <= 0.0:
            failures.append(f"{mode} optimized accelerometer charge invalid")
        if m.get("positive_residual_information_lower_form_built_here") is not False:
            failures.append(f"{mode} falsely claims residual-information lower form")
        if m.get("P4_PROMOTED") is not False:
            failures.append(f"{mode} falsely promotes P4")
    return list(dict.fromkeys(failures))


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
        "P3": d["P3_ESTABLISHED"],
        "outer_angle_rad": d["outer_angle_rad"],
        "H": d.get("modes", {}).get("H"),
        "A": d.get("modes", {}).get("A"),
        "P4": d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED"],
        "next_obligation": d["next_obligation"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
