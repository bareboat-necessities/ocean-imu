#!/usr/bin/env python3
"""Same-cell Joseph S^-1 directional-weight bridge for the OU-III P4 route.

The usable P4 proof must preserve directional information until a complete H/A
word has accumulated enough rank.  In particular it must not replace an
accepted vector packet by a fictitious positive scalar full-state margin.

For one source-correlated measurement cell let

    S = H P H^T + R,

with the actual covariance family P, measurement Jacobian family H, and
measurement covariance R kept together until S is formed.  Every concrete S is
symmetric positive definite because P is a covariance and R is positive
definite.  If s_bar is a validated upper bound on ||S||_2, then

    S^-1 >= (1/s_bar) I.

For diagonal R with r_min I <= R this gives the directional Loewner bridge

    H^T S^-1 H >= c_J H^T R^-1 H,
    c_J = r_min / s_bar > 0.

The inequality is deliberately directional: rank and nullspaces of H are
preserved.  No interval Kalman-gain matrix, covariance condition number, or
per-packet scalarization is introduced.

This producer also evaluates the bridge on a source-uniform goLive covariance
envelope for H=18 and A=21.  The vector Jacobian families over-cover every
admitted orientation componentwise, so the reported factors are conservative
initial-word diagnostics.  They are not yet complete-word factors: each prefix
covariance still has to be propagated through prediction, Joseph update and
immediate reset before P4 may be promoted.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, matrix_abs_row_sum_upper
import ou3_p4_candidate_full_word as CAND
import ou3_p4_covariance_primitives as COV
import ou3_p4_effective_vector_inputs as EFFECTIVE
import ou3_p4_joint_joseph as JOSEPH
import ou3_p4_source_word_timing as TIMING
import ou3_source_reachable_matrix_p3 as P3
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def I(x: float) -> Interval:
    return Interval.point(float(x))


def box(a: float) -> Interval:
    a = abs(float(a))
    return Interval(down(-a), up(a))


def zero(rows: int, cols: int):
    return [[I(0.0) for _ in range(cols)] for _ in range(rows)]


def _r_diagonal_bounds(R) -> tuple[float, float]:
    n = len(R)
    if n == 0 or any(len(row) != n for row in R):
        raise ValueError("R must be a nonempty square matrix")
    lo = math.inf
    hi = 0.0
    for i in range(n):
        if not R[i][i].lo > 0.0:
            raise RuntimeError("measurement covariance lost positive diagonal floor")
        lo = min(lo, R[i][i].lo)
        hi = max(hi, R[i][i].hi)
        for j in range(n):
            if i != j and not (R[i][j].lo == 0.0 and R[i][j].hi == 0.0):
                raise RuntimeError("current Joseph directional bridge requires diagonal R")
    return down(lo), up(hi)


def same_cell_bridge(P, H, R) -> dict:
    """Return a rigorous scalar Loewner lower for one correlated Joseph cell.

    For symmetric positive definite S, lambda_max(S) <= ||S||_2.  Because S is
    symmetric, ||S||_2 <= ||S||_inf, and ``matrix_abs_row_sum_upper`` provides
    an outward-rounded upper bound on that norm for every matrix in the interval
    family.  Therefore ``w=1/s_bar`` is a valid uniform lower eigenvalue bound
    for S^-1.
    """
    _, S = JOSEPH.innovation(P, H, R)
    s_bar = matrix_abs_row_sum_upper(S)
    if not (math.isfinite(s_bar) and s_bar > 0.0):
        raise RuntimeError("innovation spectral upper is not finite positive")
    w = down(1.0 / s_bar)
    rlo, rhi = _r_diagonal_bounds(R)
    c = down(w * rlo)
    if not (math.isfinite(w) and w > 0.0 and math.isfinite(c) and c > 0.0):
        raise RuntimeError("Joseph directional weight lost strict positivity")
    # Since S=HPH'+R >= R for every concrete covariance family member, the true
    # attenuation never exceeds one.  Clamp only the displayed upper-safe scalar
    # against one; a value above one would indicate arithmetic/source drift.
    if c > 1.00000000000001:
        raise RuntimeError("Joseph-vs-R attenuation exceeded one")
    return {
        "innovation_abs_row_sum_spectral_upper": s_bar,
        "S_inverse_isotropic_Loewner_weight_lower": w,
        "R_diagonal_variance_lower": rlo,
        "R_diagonal_variance_upper": rhi,
        "joseph_vs_R_inverse_directional_attenuation_lower": c,
        "directional_inequality": "H^T S^-1 H >= c_J H^T R^-1 H",
        "rank_and_nullspaces_preserved": True,
        "K_interval_matrix_materialized": False,
        "condition_number_conversion_used": False,
    }


def _wide_source_cell() -> dict:
    sched = P3.source_schedule()
    return {
        "dt_s": float(sched["dt_s"]),
        "tau_s": Interval(*map(float, sched["tau_applied_invariant_s"])),
        "sigma_aw_mps2": Interval(*map(float, sched["sigma_aw_applied_safety"])),
        "R_S_filter_std": Interval(*map(float, sched["R_S_applied_invariant"])),
        "R_S_axis_std_factors": list(map(float, sched["R_S_axis_std_factors"])),
    }


def _acc_H(mode: str, fmax: float):
    n = 18 if mode == "H" else 21
    H = zero(3, n)
    # J_att=-[f]_x.  A componentwise box with |f_i|<=||f|| over-covers every
    # admitted force orientation and magnitude.
    H[0][1] = box(fmax); H[0][2] = box(fmax)
    H[1][0] = box(fmax); H[1][2] = box(fmax)
    H[2][0] = box(fmax); H[2][1] = box(fmax)
    # J_aw=R_wb is orthogonal in the source.  [-1,1] component boxes are a
    # conservative family enclosure and intentionally do not assume a gauge.
    for i in range(3):
        for j in range(3):
            H[i][15 + j] = box(1.0)
    if mode == "A":
        # Shipping active-bias branch uses J_ba=I.
        for i in range(3):
            H[i][18 + i] = I(1.0)
    return H


def _mag_H(mode: str, mmax: float):
    n = 18 if mode == "H" else 21
    H = zero(3, n)
    H[0][1] = box(mmax); H[0][2] = box(mmax)
    H[1][0] = box(mmax); H[1][2] = box(mmax)
    H[2][0] = box(mmax); H[2][1] = box(mmax)
    return H


def _S_H(mode: str):
    n = 18 if mode == "H" else 21
    H = zero(3, n)
    for i in range(3):
        H[i][12 + i] = I(1.0)
    return H


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("Joseph directional bridge must not be trajectory fitted")

    effective = EFFECTIVE.build(path)
    timing = TIMING.build(path)
    vector = VECTOR.build()
    failures = [f"effective-input: {x}" for x in EFFECTIVE.validate(effective)]
    failures += [f"timing: {x}" for x in TIMING.validate(timing)]
    failures += [f"vector: {x}" for x in VECTOR.validate(vector)]

    if effective.get("acc_eta_in_aw_measurement_range_exact") is not True:
        failures.append("accelerometer effective-input range identity missing")
    if effective.get("mag_radial_residual_gain_null_exact") is not True:
        failures.append("magnetometer radial gain-null identity missing")
    if timing.get("S_nonlinear_eta_identically_zero") is not True:
        failures.append("S=0 nonlinear eta is not identically zero")

    src = _wide_source_cell()
    live = domain["normal_live"]
    fmax = float(live["specific_force_norm_upper_mps2"])
    mmax = float(live["magnetic_vector_norm_upper_uT"])
    vc = vector["configured_measurement_bounds"]
    Racc = COV.R_diag(float(vc["acc_measurement_std_mps2"]))
    Rmag = COV.R_diag(float(vc["mag_measurement_std_uT"]))
    RS = COV.R_S(src)

    modes = {}
    for mode in ("H", "A"):
        P0 = CAND._initial_covariance(mode, src, path)
        acc = same_cell_bridge(P0, _acc_H(mode, fmax), Racc)
        mag = same_cell_bridge(P0, _mag_H(mode, mmax), Rmag)
        s0 = same_cell_bridge(P0, _S_H(mode), RS)
        minimum = min(
            float(acc["joseph_vs_R_inverse_directional_attenuation_lower"]),
            float(mag["joseph_vs_R_inverse_directional_attenuation_lower"]),
            float(s0["joseph_vs_R_inverse_directional_attenuation_lower"]),
        )
        if not minimum > 0.0:
            failures.append(f"{mode}: initial Joseph directional attenuation is not positive")
        modes[mode] = {
            "dimension": 18 if mode == "H" else 21,
            "source_uniform_goLive_covariance_envelope": True,
            "all_orientation_vector_Jacobian_overcover": True,
            "accelerometer": acc,
            "magnetometer": mag,
            "S_zero": s0,
            "minimum_initial_directional_attenuation_lower": minimum,
            "complete_word_prefix_covariances_propagated_here": False,
            "complete_word_directional_credit_accumulated_here": False,
            "P4_PROMOTED": False,
        }

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_SOURCE_CORRELATED_JOSEPH_DIRECTIONAL_WEIGHT_BRIDGE",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "same_cell_P_H_R_correlation_retained_until_innovation": True,
        "innovation_S_SPD_by_covariance_plus_positive_R": True,
        "interval_K_materialized": False,
        "covariance_condition_number_conversion_used": False,
        "per_packet_scalar_full_state_margin_claimed": False,
        "directional_rank_nullspaces_preserved": True,
        "modes": modes,
        "P4_JOSEPH_DIRECTIONAL_WEIGHT_PRIMITIVE_ESTABLISHED": not failures,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED": False,
        "P5_FINITE_CAPTURE_ESTABLISHED": False,
        "next_obligation": (
            "propagate the source-correlated H/A covariance family through every prediction, accepted/rejected measurement and immediate reset; "
            "evaluate this same-cell S^-1 bridge at each accepted operation; transport the resulting directional forms to the word endpoint and "
            "accumulate them before taking the first generalized scalar margin"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "same_cell_P_H_R_correlation_retained_until_innovation",
        "innovation_S_SPD_by_covariance_plus_positive_R",
        "directional_rank_nullspaces_preserved",
        "P4_JOSEPH_DIRECTIONAL_WEIGHT_PRIMITIVE_ESTABLISHED",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_replay_used",
        "filter_changed",
        "interval_K_materialized",
        "covariance_condition_number_conversion_used",
        "per_packet_scalar_full_state_margin_claimed",
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED",
        "P5_FINITE_CAPTURE_ESTABLISHED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    for mode, n in (("H", 18), ("A", 21)):
        m = d.get("modes", {}).get(mode, {})
        if m.get("dimension") != n:
            f.append(f"{mode}: dimension mismatch")
        if m.get("source_uniform_goLive_covariance_envelope") is not True:
            f.append(f"{mode}: source-uniform goLive covariance envelope missing")
        for op in ("accelerometer", "magnetometer", "S_zero"):
            r = m.get(op, {})
            x = r.get("joseph_vs_R_inverse_directional_attenuation_lower")
            if not isinstance(x, (int, float)) or not math.isfinite(float(x)) or not (0.0 < float(x) <= 1.00000000000001):
                f.append(f"{mode}.{op}: invalid Joseph attenuation")
            if r.get("rank_and_nullspaces_preserved") is not True:
                f.append(f"{mode}.{op}: directional rank preservation missing")
            if r.get("K_interval_matrix_materialized") is not False:
                f.append(f"{mode}.{op}: interval K was materialized")
        if m.get("complete_word_prefix_covariances_propagated_here") is not False:
            f.append(f"{mode}: primitive prematurely claims prefix propagation")
        if m.get("complete_word_directional_credit_accumulated_here") is not False:
            f.append(f"{mode}: primitive prematurely claims word accumulation")
        if m.get("P4_PROMOTED") is not False:
            f.append(f"{mode}: P4 prematurely promoted")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain.resolve())
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": "PASS" if not vf else "FAIL",
        "H_min_initial_attenuation": d["modes"]["H"]["minimum_initial_directional_attenuation_lower"],
        "A_min_initial_attenuation": d["modes"]["A"]["minimum_initial_directional_attenuation_lower"],
        "P4": d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED"],
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
