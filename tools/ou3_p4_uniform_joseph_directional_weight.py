#!/usr/bin/env python3
"""Source-uniform Joseph directional weights from the retained P3 metric.

The initial Joseph bridge is not enough for a complete word: P4 needs a lower
S^-1 weight at every recurrent source node.  P3 already supplies the required
covariance object.  Its matrix upper argument constructs a diagonal dominator
``D_m`` with

    P_m <= D_m

in Loewner order for every source-uniform H/A metric cell represented by the
certificate.  P3 also certifies homogeneous prefix information gain at most one
in those matching source metrics.

For any concrete measurement Jacobian H and positive diagonal R,

    S = H P H^T + R <= H D_m H^T + R.

A validated spectral upper ``s_bar`` on the right-hand side gives

    S^-1 >= (1/s_bar) I

and therefore the same directional comparison used by the same-cell bridge,

    H^T S^-1 H >= c_J H^T R^-1 H,
    c_J = r_min/s_bar > 0.

This producer computes conservative H/A constants using an all-orientation
vector-Jacobian over-cover and the full deployed R_S source interval.  Rank and
nullspaces remain untouched.  It does not yet add the pulled-back forms over a
complete word and therefore does not promote P4.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p3_source_uniform_certificate as P3CERT
import ou3_p4_cayley_sector_certificate as SECTOR
import ou3_p4_covariance_primitives as COV
import ou3_p4_joseph_directional_weight as J
import ou3_source_reachable_matrix_p3 as P3
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def I(x: float) -> Interval:
    return Interval.point(float(x))


def _dominator(diagonal_upper: list[float]):
    n = len(diagonal_upper)
    if n not in (18, 21):
        raise RuntimeError("P3 covariance dominator has unexpected dimension")
    out = [[I(0.0) for _ in range(n)] for _ in range(n)]
    for i, x in enumerate(diagonal_upper):
        x = float(x)
        if not (math.isfinite(x) and x > 0.0):
            raise RuntimeError("P3 covariance diagonal dominator is not finite positive")
        out[i][i] = Interval.outward_bounds(x, x)
    return out


def _wide_source_cell() -> dict:
    sched = P3.source_schedule()
    return {
        "dt_s": float(sched["dt_s"]),
        "tau_s": Interval(*map(float, sched["tau_applied_invariant_s"])),
        "sigma_aw_mps2": Interval(*map(float, sched["sigma_aw_applied_safety"])),
        "R_S_filter_std": Interval(*map(float, sched["R_S_applied_invariant"])),
        "R_S_axis_std_factors": list(map(float, sched["R_S_axis_std_factors"])),
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("uniform Joseph weights must not be trajectory fitted")

    p3 = P3CERT.build(path)
    sector = SECTOR.build(path)
    vector = VECTOR.build()
    failures = [f"P3: {x}" for x in P3CERT.validate(p3)]
    failures += [f"sector: {x}" for x in SECTOR.validate(sector)]
    failures += [f"vector: {x}" for x in VECTOR.validate(vector)]

    matrix = p3.get("matrix_certificate", {})
    upper_argument = str(matrix.get("matrix_upper_argument", ""))
    if "Loewner upper relation" not in upper_argument:
        failures.append("P3 matrix upper argument no longer certifies a Loewner dominator")
    for mode in ("H", "A"):
        if p3["modes"][mode].get("prefix_information_gain_upper") != 1.0:
            failures.append(f"{mode}: P3 prefix information gain is not one")

    live = domain["normal_live"]
    fmax = float(live["specific_force_norm_upper_mps2"])
    mmax = float(live["magnetic_vector_norm_upper_uT"])
    vc = vector["configured_measurement_bounds"]
    Racc = COV.R_diag(float(vc["acc_measurement_std_mps2"]))
    Rmag = COV.R_diag(float(vc["mag_measurement_std_uT"]))
    RS = COV.R_S(_wide_source_cell())

    alpha6_R = float(P3.vector_alpha6(live, vector))
    retention = float(sector["exact_vector_information_retention_factor_lower"])
    if not (alpha6_R > 0.0 and 0.0 < retention <= 1.0):
        failures.append("vector/sector directional prerequisites are not strictly positive")

    modes = {}
    for mode, n in (("H", 18), ("A", 21)):
        diag = list(map(float, p3["modes"][mode]["Sigma_diagonal_upper"]))
        if len(diag) != n:
            failures.append(f"{mode}: P3 Sigma diagonal upper has wrong dimension")
            continue
        D = _dominator(diag)
        acc = J.same_cell_bridge(D, J._acc_H(mode, fmax), Racc)
        mag = J.same_cell_bridge(D, J._mag_H(mode, mmax), Rmag)
        s0 = J.same_cell_bridge(D, J._S_H(mode), RS)
        c_vec = down(min(
            float(acc["joseph_vs_R_inverse_directional_attenuation_lower"]),
            float(mag["joseph_vs_R_inverse_directional_attenuation_lower"]),
        ))
        c_s = float(s0["joseph_vs_R_inverse_directional_attenuation_lower"])
        alpha_sector_R = down(alpha6_R * retention)
        alpha_sector_J = down(alpha_sector_R * c_vec)
        if not (c_vec > 0.0 and c_s > 0.0 and alpha_sector_J > 0.0):
            failures.append(f"{mode}: source-uniform Joseph directional weight lost positivity")
        modes[mode] = {
            "dimension": n,
            "P3_Sigma_diagonal_Loewner_dominator": diag,
            "P3_prefix_information_gain_upper": 1.0,
            "accelerometer": acc,
            "magnetometer": mag,
            "S_zero": s0,
            "uniform_vector_Joseph_attenuation_lower": c_vec,
            "uniform_S_Joseph_attenuation_lower": c_s,
            "P3_R_inverse_vector_gyro_bias_alpha6_lower": alpha6_R,
            "finite_angle_vector_information_retention_lower": retention,
            "uniform_finite_angle_Joseph_vector_gyro_bias_alpha6_lower": alpha_sector_J,
            "directional_rank_preserved": True,
            "complete_word_directional_sum_emitted_here": False,
            "P4_PROMOTED": False,
        }

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_SOURCE_UNIFORM_JOSEPH_DIRECTIONAL_WEIGHT_FROM_P3_METRIC",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "P3_covariance_Loewner_upper_consumed": True,
        "P3_prefix_metric_contract_consumed": True,
        "same_cell_P_H_R_correlation_replaced_by_valid_P3_Loewner_upper_only": True,
        "interval_K_materialized": False,
        "condition_number_conversion_used": False,
        "per_packet_full_state_scalarization_used": False,
        "modes": modes,
        "P4_UNIFORM_JOSEPH_DIRECTIONAL_WEIGHTS_ESTABLISHED": not failures,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED": False,
        "P5_FINITE_CAPTURE_ESTABLISHED": False,
        "next_obligation": (
            "pull the uniformly weighted vector and S directional forms through the invertible source predictions and exact reset congruences, "
            "sum them on a complete recurrent H/A word, and certify a positive full-rank generalized endpoint margin before charging the remaining Cayley reset/effective-input nonlinear defect"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "P3_covariance_Loewner_upper_consumed",
        "P3_prefix_metric_contract_consumed",
        "same_cell_P_H_R_correlation_replaced_by_valid_P3_Loewner_upper_only",
        "P4_UNIFORM_JOSEPH_DIRECTIONAL_WEIGHTS_ESTABLISHED",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_replay_used", "filter_changed", "interval_K_materialized",
        "condition_number_conversion_used", "per_packet_full_state_scalarization_used",
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED", "P5_FINITE_CAPTURE_ESTABLISHED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    for mode, n in (("H", 18), ("A", 21)):
        m = d.get("modes", {}).get(mode, {})
        if m.get("dimension") != n:
            f.append(f"{mode}: dimension mismatch")
        if len(m.get("P3_Sigma_diagonal_Loewner_dominator", [])) != n:
            f.append(f"{mode}: missing P3 covariance dominator")
        for key in (
            "uniform_vector_Joseph_attenuation_lower",
            "uniform_S_Joseph_attenuation_lower",
            "P3_R_inverse_vector_gyro_bias_alpha6_lower",
            "finite_angle_vector_information_retention_lower",
            "uniform_finite_angle_Joseph_vector_gyro_bias_alpha6_lower",
        ):
            x = m.get(key)
            if not isinstance(x, (int, float)) or not math.isfinite(float(x)) or not float(x) > 0.0:
                f.append(f"{mode}.{key} is not finite positive")
        if m.get("directional_rank_preserved") is not True:
            f.append(f"{mode}: directional rank not preserved")
        if m.get("complete_word_directional_sum_emitted_here") is not False:
            f.append(f"{mode}: prematurely emits complete word sum")
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
        "H_vector_cJ": d.get("modes", {}).get("H", {}).get("uniform_vector_Joseph_attenuation_lower"),
        "A_vector_cJ": d.get("modes", {}).get("A", {}).get("uniform_vector_Joseph_attenuation_lower"),
        "H_S_cJ": d.get("modes", {}).get("H", {}).get("uniform_S_Joseph_attenuation_lower"),
        "A_S_cJ": d.get("modes", {}).get("A", {}).get("uniform_S_Joseph_attenuation_lower"),
        "H_alpha_J_sector": d.get("modes", {}).get("H", {}).get("uniform_finite_angle_Joseph_vector_gyro_bias_alpha6_lower"),
        "A_alpha_J_sector": d.get("modes", {}).get("A", {}).get("uniform_finite_angle_Joseph_vector_gyro_bias_alpha6_lower"),
        "P4": d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED"],
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
