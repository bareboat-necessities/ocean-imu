#!/usr/bin/env python3
"""Windowed Normal-Live vector-PE certificate for canonical OU-III P3.

Canonical P3 must use the theorem's recurring asynchronous PE premise, not the
narrower historical shortcut that required two consecutive 25 Hz magnetic
packets.  The declared Normal-Live domain gives a recurrence window W_PE and
physical vector/rate bounds.  Select one accepted PE occurrence from [0,W_PE]
and another from [2 W_PE,3 W_PE].  Their separation is therefore in
[W_PE,3 W_PE] regardless of magnetometer ODR, rejected magnetic packets, or
short outages between the required occurrences.

At each required occurrence the declared specific-force and magnetic norms and
the declared sine separation give a strict attitude-information lower bound.
The two separated occurrences then distinguish attitude from gyro bias.  This
is a sufficient machine witness for the paper's finite-window condition

    W_eta6(k,N) >= alpha_6 I_6,

while retaining the paper's asynchronous semantics.  It is not claimed to be
the only way the PE condition can hold.

For A mode the paper permits the same eta6 condition together with the finite
residual accelerometer-bias correlation time.  We therefore report the exact
homogeneous bias contraction over the same word rather than inventing an eta9
pointwise packet condition.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_full_process_ucc as PROCESS
import ou3_validated_transcendentals as VT
import ou3_vector_uco_certificate as VECTOR
from ou3_interval import Interval

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
PAPER = REPO / "doc" / "kalman_ou_iii" / "w3d-iss-stability.tex-part"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_WINDOWED_ASYNCHRONOUS_VECTOR_PE"


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _exp_negative_wide(x: float) -> Interval:
    if not (math.isfinite(x) and x >= 0.0):
        raise ValueError("finite nonnegative exponential argument required")
    scale = 1
    while x / scale > VT.MAX_ABS_ARGUMENT:
        scale *= 2
    y = VT.exp_interval(Interval.outward_bounds(-x / scale, -x / scale))
    s = scale
    while s > 1:
        y = y.square()
        s //= 2
    return y


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("windowed PE certificate may not be trajectory fitted")
    live = domain["normal_live"]

    vector = VECTOR.build()
    vf = VECTOR.validate(vector)
    process = PROCESS.build()
    pf = PROCESS.validate(process)
    if vf or pf:
        raise RuntimeError(f"upstream source constants invalid: vector={vf}, process={pf}")

    paper = PAPER.read_text(encoding="utf-8")
    paper_markers = (
        "finite-window asynchronous conditions",
        "Unequal sensor rates,",
        "rejected samples, and short outages are admissible",
        "No pointwise accelerometer--magnetometer",
        "\\alpha_6\\mat I_6\\preceq",
        "finite residual-bias correlation time",
    )
    missing = [m for m in paper_markers if m not in paper]
    if missing:
        raise RuntimeError(f"paper PE contract changed: {missing}")

    W = float(live["vector_pe_recurrence_window_s"])
    fmin = float(live["specific_force_norm_lower_mps2"])
    mmin = float(live["magnetic_vector_norm_lower_uT"])
    sine = float(live["vector_sine_separation_lower"])
    rate_deg = float(live["body_rate_norm_upper_deg_s"])
    if not (W > 0.0 and fmin > 0.0 and mmin > 0.0 and 0.0 < sine < 1.0 and rate_deg > 0.0):
        raise RuntimeError("invalid declared Normal-Live PE envelope")

    # Two recurrence cells separated by one whole recurrence cell.
    # t0 in [0,W], t1 in [2W,3W].
    delta_min = down(W)
    delta_max = up(3.0 * W)
    word_horizon = delta_max

    vc = vector["configured_measurement_bounds"]
    ra = up(float(vc["acc_measurement_variance_upper"]))
    rm = up(float(vc["mag_measurement_variance_upper"]))
    if not (ra > 0.0 and rm > 0.0):
        raise RuntimeError("measurement variance upper lost positivity")

    # Stable evaluation of 1-sqrt(1-s^2).
    root = up(math.sqrt(max(0.0, 1.0 - sine * sine)))
    angular = down((sine * sine) / up(1.0 + root))
    mu_theta = down(min(fmin * fmin / ra, mmin * mmin / rm) * angular)

    omega = up(rate_deg * math.pi / 180.0)
    # Same exact-transport lower used by the existing vector certificate, but
    # now over source-guaranteed separated recurrence cells rather than hardware ODR.
    bracket = down(1.0 - up(0.5 * omega * delta_max))
    if not bracket > 0.0:
        raise RuntimeError(
            "declared PE recurrence/rate box is too wide for the current two-occurrence transport bound"
        )
    gamma_sigma = down(delta_min * bracket)
    # The scaled eta6 coordinates use a one-second bg scale in the retained theorem.
    Tbg = float(vector["operating_envelope"]["gyro_bias_time_scale_s"])
    gamma = down(gamma_sigma / Tbg)
    alpha6 = down(mu_theta / up(1.0 + up(2.0 / down(gamma * gamma))))

    tau_ba = float(process["source_constants"]["accel_bias_tau_s"])
    ba_exp = _exp_negative_wide(word_horizon / tau_ba)
    ba_contraction_upper = ba_exp.hi
    ba_contraction_gap_lower = down(1.0 - ba_contraction_upper)

    passed = all(math.isfinite(x) and x > 0.0 for x in (
        delta_min, delta_max, angular, mu_theta, bracket, gamma_sigma, gamma,
        alpha6, ba_contraction_gap_lower,
    )) and ba_contraction_upper < 1.0

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "trajectory_fit": False,
        "paper_windowed_PE_semantics_consumed": True,
        "asynchronous_magnetometer_semantics_consumed": True,
        "hardware_magnetometer_ODR_used_as_PE_recurrence": False,
        "two_consecutive_accepted_magnetic_packets_required": False,
        "rejected_magnetic_packets_between_required_occurrences_allowed": True,
        "all_valid_accelerometer_packets_required": bool(
            live["accelerometer_update_required_each_valid_imu_sample_after_live_entry"]
        ),
        "accelerometer_rejection_branch_present": bool(
            live["accelerometer_rejection_in_normal_live_scope"]
        ),
        "declared_normal_live_PE": {
            "recurrence_window_s": W,
            "specific_force_norm_lower_mps2": fmin,
            "magnetic_vector_norm_lower_uT": mmin,
            "vector_sine_separation_lower": sine,
            "body_rate_norm_upper_deg_s": rate_deg,
        },
        "spread_occurrence_selection": {
            "first_occurrence_window_s": [0.0, W],
            "second_occurrence_window_s": [2.0 * W, 3.0 * W],
            "separation_lower_s": delta_min,
            "separation_upper_s": delta_max,
            "word_horizon_s": word_horizon,
        },
        "eta6_information": {
            "acc_measurement_variance_upper": ra,
            "mag_measurement_variance_upper": rm,
            "angular_factor_lower": angular,
            "single_occurrence_attitude_information_lower": mu_theta,
            "transport_bracket_lower": bracket,
            "Gamma_g_sigma_min_lower_s": gamma_sigma,
            "alpha_6_information_lower": alpha6,
        },
        "A_mode_bias_route": {
            "uses_eta9_pointwise_packet_shortcut": False,
            "uses_eta6_plus_finite_bias_correlation": True,
            "accel_bias_tau_s": tau_ba,
            "homogeneous_bias_contraction_upper_over_word": ba_contraction_upper,
            "homogeneous_bias_contraction_gap_lower": ba_contraction_gap_lower,
        },
        "pass": passed,
        "P3_promoted": False,
        "next_obligation": (
            "compose this windowed eta6 information, the finite A-mode bias contraction, all accepted accelerometer updates, recurrent S=0 information, process UCC, covariance-floor events, and the joint SEA3 adaptive source in one H18/A21 word"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for key in (
        "paper_windowed_PE_semantics_consumed",
        "asynchronous_magnetometer_semantics_consumed",
        "rejected_magnetic_packets_between_required_occurrences_allowed",
        "all_valid_accelerometer_packets_required",
        "pass",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_fit",
        "hardware_magnetometer_ODR_used_as_PE_recurrence",
        "two_consecutive_accepted_magnetic_packets_required",
        "accelerometer_rejection_branch_present",
        "P3_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    alpha = d.get("eta6_information", {}).get("alpha_6_information_lower")
    if not isinstance(alpha, (int, float)) or not (math.isfinite(float(alpha)) and float(alpha) > 0.0):
        f.append("eta6 alpha_6 lower is not strict")
    ba = d.get("A_mode_bias_route", {})
    if ba.get("uses_eta6_plus_finite_bias_correlation") is not True:
        f.append("A-mode finite-bias-correlation route missing")
    if ba.get("uses_eta9_pointwise_packet_shortcut") is not False:
        f.append("A mode reintroduced pointwise eta9 shortcut")
    c = ba.get("homogeneous_bias_contraction_upper_over_word")
    if not isinstance(c, (int, float)) or not (0.0 < float(c) < 1.0):
        f.append("A-mode bias contraction is not strict")
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
        "word_horizon_s": d["spread_occurrence_selection"]["word_horizon_s"],
        "alpha_6_lower": d["eta6_information"]["alpha_6_information_lower"],
        "A_bias_contraction_gap_lower": d["A_mode_bias_route"]["homogeneous_bias_contraction_gap_lower"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
