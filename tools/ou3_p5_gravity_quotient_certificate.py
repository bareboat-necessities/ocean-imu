#!/usr/bin/env python3
"""Detectable gravity-only quotient prerequisite for OU-III P5.

The full yaw-only quotient cannot be strictly contracting because a body-frame
gyro-bias error parallel to gravity is exact zero dynamics when body rate is
zero: it produces only yaw, which the quotient removes.  This producer builds
the corrected deterministic quotient used by the paper.

Strict coordinates are

    (tilt_perp, b_g,perp, v, p, S, a_w),

while absolute yaw is quotiented and b_g,parallel is retained in the source
state only as a bounded detectability input.  No estimator state is removed or
changed; this is a theorem-coordinate construction.

For each accepted quotient gravity packet assume

    f_b = -g u_g + eps,   ||eps|| <= f_q < g.

Then f_min=g-f_q and the angle between -f_b and u_g has cosine at least
c_q=(g-f_q)/(g+f_q).  Hence for every tangent tilt x perpendicular to u_g,

    x' [f]_x' [f]_x x >= f_min^2 c_q^2 ||x||^2.

Two accepted packets separated by Delta in [Delta_min,Delta_max] expose the two
transverse gyro-bias coordinates.  With the same transport estimate as the
full-heading analytical certificate,

    gamma_q >= Delta_min/T_bg * (1-omega_bar Delta_max/2) > 0

and the scaled 4-state quotient Gramian has the explicit lower bound

    alpha_4q >= mu_q/(1+2/gamma_q^2) > 0,
    mu_q = f_min^2 c_q^2 / r_a,+.

The complete translational [v,p,S,a_w] part remains covered by the existing
four-S source-word qualification.  This file therefore establishes the reduced
observable/detectable quotient geometry and its strict information constant; it
does not yet claim the complete nonlinear quotient funnel word.  That later
word must retain b_g,parallel as an input term and enclose every exact prefix.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_implementation_word_language as WORDS
import ou3_p5_yaw_quotient_word_certificate as OLDQ
import ou3_translational_uco_ucc as TRANS
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
T_BG_S = 1.0


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("gravity quotient domain must not be trajectory fitted")

    old = OLDQ.build(domain_path)
    old_fail = OLDQ.validate(old)
    words = WORDS.build(domain_path)
    word_fail = WORDS.validate(words)
    trans = TRANS.build(TRANS.DEFAULT_HEADER)
    trans_fail = TRANS.validate(trans)
    vector = VECTOR.build()
    vector_fail = VECTOR.validate(vector)

    failures = []
    failures += [f"yaw-only audit: {x}" for x in old_fail]
    failures += [f"word language: {x}" for x in word_fail]
    failures += [f"translation: {x}" for x in trans_fail]
    failures += [f"vector configuration: {x}" for x in vector_fail]

    live = domain["normal_live"]
    qh = live.get("gravity_quotient", {})
    g = float(domain["startup"]["gravity_mps2"])
    fq = float(qh.get("non_gravitational_specific_force_norm_upper_mps2", math.nan))
    recurrence = float(qh.get("accepted_packet_recurrence_window_s", math.nan))
    dmin = float(qh.get("accepted_packet_separation_min_s", math.nan))
    dmax = float(qh.get("accepted_packet_separation_max_s", math.nan))
    omega = up(float(live["body_rate_norm_upper_deg_s"]) * math.pi / 180.0)
    if not (0.0 <= fq < g):
        failures.append("quotient non-gravitational force bound must be below gravity")
    if not (0.0 < dmin <= dmax <= recurrence):
        failures.append("quotient accepted-packet separation/recurrence is invalid")
    if qh.get("axial_gyro_bias_is_neutral_input") is not True:
        failures.append("quotient domain does not mark axial gyro bias as neutral input")

    acc_var_upper = float(vector["configured_measurement_bounds"]["acc_measurement_variance_upper"])
    if not (math.isfinite(acc_var_upper) and acc_var_upper > 0.0):
        failures.append("configured accelerometer variance upper is not positive")

    fmin = down(g - fq)
    c_align = down(fmin / up(g + fq)) if fmin > 0.0 else 0.0
    mu_q = down(down(fmin * fmin) * down(c_align * c_align) / acc_var_upper) if c_align > 0.0 else 0.0
    bracket = down(1.0 - up(0.5 * omega * dmax))
    gamma = down(dmin * bracket / T_BG_S) if bracket > 0.0 else 0.0
    gamma2 = down(gamma * gamma) if gamma > 0.0 else 0.0
    denom = up(1.0 + up(2.0 / gamma2)) if gamma2 > 0.0 else math.inf
    alpha4 = down(mu_q / denom) if math.isfinite(denom) and denom > 0.0 else 0.0

    strict = all(math.isfinite(x) and x > 0.0 for x in (fmin, c_align, mu_q, bracket, gamma, alpha4))
    if not strict:
        failures.append("reduced gravity quotient information lower bound is not strict")

    translation_complete = bool(trans.get("translation_source_complete"))
    if not translation_complete:
        failures.append("four-S translation qualification is not source complete")

    old_status = old.get("P5_YAW_ONLY_QUOTIENT_OBSTRUCTION_IDENTIFIED")
    if old_status != "PASS":
        failures.append("old yaw-only quotient obstruction was not independently established")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_DETECTABLE_GRAVITY_ONLY_QUOTIENT_PREREQUISITE",
        "claim": "YAW_QUOTIENT_WITH_TRANSVERSE_GYRO_BIAS_STRICT_AND_AXIAL_BIAS_AS_INPUT",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "old_yaw_only_full_bias_metric_rejected": True,
        "absolute_yaw_quotiented": True,
        "axial_gyro_bias_removed_from_strict_metric": True,
        "axial_gyro_bias_removed_from_filter_state": False,
        "axial_gyro_bias_role": "NEUTRAL_BOUNDED_DETECTABILITY_INPUT_UNTIL_EXCITATION_OR_MAGNETIC_REGAUGE",
        "detectable_coordinates": [
            "gravity_tilt_2d",
            "gyro_bias_transverse_2d",
            "v_3d",
            "p_3d",
            "S_3d",
            "a_w_3d"
        ],
        "schur_metric_policy": {
            "horizontal_yaw_minimized_by_group_action": True,
            "axial_bias_minimized_only_for_detectable_Schur_metric": True,
            "axial_bias_declared_physical_symmetry": False,
            "source_information_cross_terms_retained": True,
        },
        "gravity_packet_hypothesis": {
            "recurrence_window_s": recurrence,
            "packet_separation_s": [dmin, dmax],
            "non_gravitational_specific_force_norm_upper_mps2": fq,
            "gravity_mps2": g,
            "specific_force_norm_lower_mps2": fmin,
            "gravity_alignment_cosine_lower": c_align,
            "accepted_accelerometer_packets_per_word_lower": 2,
            "hypothesis_origin": "DECLARED_DEPLOYMENT_THEOREM_ASSUMPTION_NOT_REPLAY_FIT",
        },
        "reduced_attitude_bias_information": {
            "acc_measurement_variance_upper": acc_var_upper,
            "tilt_packet_information_mu_lower": mu_q,
            "body_rate_rad_s_upper": omega,
            "transport_bracket_lower": bracket,
            "scaled_Gamma_sigma_min_lower": gamma,
            "alpha_4_quotient_information_lower": alpha4,
            "strict": strict,
        },
        "translation_word": {
            "state_order": ["v", "p", "S", "a_w"],
            "four_S_complete_chain_required": True,
            "source_complete": translation_complete,
            "Q_axis_lambda_min_lower": trans["process_ucc"]["Q_axis_lambda_min_lower"],
        },
        "old_yaw_only_counterexample_consumed": True,
        "old_yaw_only_counterexample_status": old_status,
        "P5_GRAVITY_QUOTIENT_REDUCED_DETECTABILITY_CERTIFICATE": "PASS" if not failures else "FAIL",
        "P5_UNGAUGED_TIMEOUT_QUOTIENT_WORD_CERTIFICATE": "NOT_ESTABLISHED",
        "next_obligation": (
            "outward-enclose the complete exact gravity-only H source word in the detectable Schur metric, "
            "with b_g_parallel charged as an input and every quaternion/Joseph/S-prefix kept source correlated"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for flag in ("source_generated_not_trajectory_fit", "absolute_yaw_quotiented",
                 "axial_gyro_bias_removed_from_strict_metric", "old_yaw_only_full_bias_metric_rejected"):
        if d.get(flag) is not True:
            failures.append(f"{flag} is not true")
    if d.get("source_replay_used") is not False:
        failures.append("gravity quotient uses replay")
    if d.get("filter_changed") is not False:
        failures.append("gravity quotient changes filter")
    if d.get("axial_gyro_bias_removed_from_filter_state") is not False:
        failures.append("proof incorrectly removes axial gyro bias from estimator state")
    if d.get("axial_gyro_bias_role") != "NEUTRAL_BOUNDED_DETECTABILITY_INPUT_UNTIL_EXCITATION_OR_MAGNETIC_REGAUGE":
        failures.append("axial gyro-bias role is not the corrected neutral-input role")
    sm = d.get("schur_metric_policy", {})
    if sm.get("axial_bias_declared_physical_symmetry") is not False:
        failures.append("axial bias was incorrectly declared a gauge symmetry")
    if sm.get("source_information_cross_terms_retained") is not True:
        failures.append("detectable Schur metric discards source information cross terms")
    info = d.get("reduced_attitude_bias_information", {})
    if info.get("strict") is not True:
        failures.append("reduced attitude/bias information is not strict")
    alpha = info.get("alpha_4_quotient_information_lower")
    if not (isinstance(alpha, (int, float)) and math.isfinite(float(alpha)) and float(alpha) > 0.0):
        failures.append("alpha_4 quotient lower is not finite positive")
    tr = d.get("translation_word", {})
    if tr.get("four_S_complete_chain_required") is not True or tr.get("source_complete") is not True:
        failures.append("quotient does not retain complete four-S translation word")
    if d.get("old_yaw_only_counterexample_consumed") is not True:
        failures.append("old yaw-only zero-dynamics witness was not consumed")
    if not failures and d.get("P5_GRAVITY_QUOTIENT_REDUCED_DETECTABILITY_CERTIFICATE") != "PASS":
        failures.append("reduced gravity quotient prerequisite did not pass")
    if d.get("P5_UNGAUGED_TIMEOUT_QUOTIENT_WORD_CERTIFICATE") != "NOT_ESTABLISHED":
        failures.append("complete nonlinear quotient word was promoted prematurely")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve())
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_GRAVITY_QUOTIENT_REDUCED_DETECTABILITY_CERTIFICATE"],
        "gravity_packet": out["gravity_packet_hypothesis"],
        "information": out["reduced_attitude_bias_information"],
        "translation": out["translation_word"],
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
