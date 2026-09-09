#!/usr/bin/env python3
"""Source-uniform same-cell Joseph correction-domain diagnostic for P4.

For one accepted three-vector Joseph update, with the ACTUAL same-cell
P,H,R,S,K, let d_theta be the attitude component of d=K y.  The covariance
reduction identity gives

    K S K^T = P^- - P^+ <= P^- .

Hence, using the first three rows only,

    ||d_theta||^2
      <= lambda_max(K_theta S K_theta^T) * y^T S^-1 y
      <= lambda_max(P^-_{theta,theta}) * y^T R^-1 y.

This is not a rowwise K bound and does not detach K from P/H/R.  It is a
same-cell Joseph metric consequence.  The source-uniform Riccati tube supplies
an outward upper bound on the attitude covariance block; the deterministic hard
entry set and exact finite-angle geometry supply event residual bounds.

The result is intentionally diagnostic.  If the certified correction ceilings
fit inside the exact-reset chart utility, they can be consumed as production
reset-domain certificates.  If they are too large, the production proof must
retain more of the same-cell residual/K direction instead of weakening the
bound further.
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
import ou3_p4_cayley_sector_certificate as CAYLEY

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
QUALIFICATION = "OU3_P4_SAME_CELL_JOSEPH_CORRECTION_DOMAIN_DIAGNOSTIC_V1"


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def _residual_bounds(domain: dict, entry: dict, dynamic: dict) -> dict:
    live = domain["normal_live"]
    q = float(entry["coordinate_radii"]["attitude_cayley_norm"])
    theta = 2.0 * math.atan(0.5 * q)
    rot_minus_I = up(2.0 * math.sin(0.5 * theta))
    fmax = float(live["specific_force_norm_upper_mps2"])
    mmax = float(live["magnetic_vector_norm_upper_uT"])
    aw = float(entry["coordinate_radii"]["latent_acceleration_norm_mps2"])
    ba = float(entry["coordinate_radii"].get("accelerometer_bias_error_norm_mps2", 0.0))
    S = float(entry["coordinate_radii"]["integral_displacement_norm_m_s"])

    noise = domain["configured_runtime"]["measurement_noise_std"]
    acc_std = min(map(float, noise["accelerometer_mps2"]))
    mag_std = min(map(float, noise["magnetometer_uT"]))
    if not (acc_std > 0.0 and mag_std > 0.0):
        raise RuntimeError("configured vector measurement std lost positivity")

    inv = dynamic["dynamic_invariant"]
    rs_lo = float(inv["R_S_applied"][0])
    # Shipping horizontal factors are 0.72 and vertical is 1.0.  The smallest
    # standard deviation gives the largest R^-1 energy and is therefore safe.
    rs_std_min = down(rs_lo * 0.72)
    if not rs_std_min > 0.0:
        raise RuntimeError("actual applied R_S lower lost positivity")

    acc_norm = up(rot_minus_I * fmax + aw + ba)
    mag_norm = up(rot_minus_I * mmax)
    return {
        "attitude_cayley_norm_upper": q,
        "attitude_angle_rad": theta,
        "rotation_minus_identity_norm_upper": rot_minus_I,
        "accelerometer_residual_norm_upper": acc_norm,
        "magnetometer_residual_norm_upper": mag_norm,
        "S_zero_residual_norm_upper": S,
        "accelerometer_measurement_std_lower": acc_std,
        "magnetometer_measurement_std_lower": mag_std,
        "accelerometer_Rinv_energy_upper": up(acc_norm * acc_norm / down(acc_std * acc_std)),
        "magnetometer_Rinv_energy_upper": up(mag_norm * mag_norm / down(mag_std * mag_std)),
        "S_zero_Rinv_energy_upper": up(S * S / down(rs_std_min * rs_std_min)),
        "actual_RS_std_lower_with_horizontal_factor": rs_std_min,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    tube = TUBE.build(path)
    entry = ENTRY.build()
    dynamic = DYNAMIC.build(path)
    cayley = CAYLEY.build(path)
    bad = {
        "tube": TUBE.validate(tube),
        "entry": ENTRY.validate(entry),
        "dynamic": DYNAMIC.validate(dynamic),
        "cayley": CAYLEY.validate(cayley),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError("same-cell correction-domain prerequisites failed: " + repr(bad))

    residual = _residual_bounds(domain, entry, dynamic)
    modes = {}
    for mode, key in (("H18", "H"), ("A21", "A")):
        p = list(map(float, tube["modes"][key]["Pbar_diagonal_variance_upper"]))
        patt_trace = up(sum(p[:3]))
        events = {}
        for name, energy in (
            ("S_zero", residual["S_zero_Rinv_energy_upper"]),
            ("accelerometer", residual["accelerometer_Rinv_energy_upper"]),
            ("magnetometer", residual["magnetometer_Rinv_energy_upper"]),
        ):
            delta2 = up(patt_trace * float(energy))
            delta = up(math.sqrt(max(0.0, delta2)))
            events[name] = {
                "attitude_covariance_trace_upper": patt_trace,
                "residual_Rinv_energy_upper": float(energy),
                "same_cell_attitude_correction_norm_upper": delta,
                "inside_reset_utility_domain": delta <= RESET.CAYLEY_MONOTONE_NORM_MAX,
                "inside_one_radian": delta <= 1.0,
            }
        worst_name = max(events, key=lambda n: events[n]["same_cell_attitude_correction_norm_upper"])
        worst = events[worst_name]["same_cell_attitude_correction_norm_upper"]
        modes[mode] = {
            "attitude_covariance_trace_upper": patt_trace,
            "events": events,
            "limiting_event": worst_name,
            "correction_norm_upper": worst,
            "reset_utility_domain_closed": worst <= RESET.CAYLEY_MONOTONE_NORM_MAX,
        }

    all_reset = all(m["reset_utility_domain_closed"] for m in modes.values())
    return {
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "same_cell_Joseph_covariance_identity_consumed": True,
        "identity": "K*S*K^T=Pminus-Pplus<=Pminus",
        "attitude_block_only_used": True,
        "rowwise_K_bound_used": False,
        "independent_K_box_used": False,
        "hard_entry_residual_geometry_consumed": True,
        "source_uniform_Riccati_attitude_covariance_consumed": True,
        "actual_applied_RS_lower_consumed": True,
        "residual_bounds": residual,
        "modes": modes,
        "all_events_inside_exact_reset_utility_domain": all_reset,
        "production_same_graph_correction_domain_promoted_here": False,
        "endpoint_augmented_LDLT_closed_here": False,
        "every_prefix_augmented_LDLT_closed_here": False,
        "P4_promoted_here": False,
    }


def validate(d: dict) -> list[str]:
    f = []
    if d.get("qualification") != QUALIFICATION:
        f.append("qualification mismatch")
    for k in (
        "same_cell_Joseph_covariance_identity_consumed",
        "attitude_block_only_used",
        "hard_entry_residual_geometry_consumed",
        "source_uniform_Riccati_attitude_covariance_consumed",
        "actual_applied_RS_lower_consumed",
    ):
        if d.get(k) is not True:
            f.append(k + " not true")
    for k in (
        "rowwise_K_bound_used",
        "independent_K_box_used",
        "production_same_graph_correction_domain_promoted_here",
        "endpoint_augmented_LDLT_closed_here",
        "every_prefix_augmented_LDLT_closed_here",
        "P4_promoted_here",
    ):
        if d.get(k) is not False:
            f.append(k + " not false")
    for mode, m in d.get("modes", {}).items():
        for ev, e in m.get("events", {}).items():
            x = float(e.get("same_cell_attitude_correction_norm_upper", math.nan))
            if not (math.isfinite(x) and x >= 0.0):
                f.append(mode + " " + ev + " correction bound invalid")
    return f


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain)
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "modes": d["modes"],
        "all_reset_utility": d["all_events_inside_exact_reset_utility_domain"],
        "failures": f,
    }, sort_keys=True))
    return int(bool(f))


if __name__ == "__main__":
    raise SystemExit(main())
