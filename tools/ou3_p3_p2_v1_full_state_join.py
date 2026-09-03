#!/usr/bin/env python3
"""Join source-complete translation P3 with H/A attitude/bias blocks.

The correlated translation certificate computes a post-measurement conditional
covariance lower for T=(v,p,S,a_w) while conditioning attitude/bias selected
modes known.  To obtain a full H/A lower without charging the old whole-word
translation information loss a second time, select an independent final-step
attitude/gyro-bias process block (and accelerometer-bias process block in A
mode), condition T known, and certify its final accelerometer/magnetometer
posterior.

For the joint selected Gaussian posterior let the precision be

    J = [[J_TT, J_TB], [J_BT, J_BB]].

Because J is positive semidefinite, block Cauchy-Schwarz gives

    J <= 2 diag(J_TT, J_BB),

hence by order reversal under inversion

    J^-1 >= 1/2 diag(J_TT^-1, J_BB^-1).

The two inverse diagonal precision blocks are exactly the conditional
covariances Cov(T|B,y) and Cov(B|T,y).  Therefore any rigorous lower bounds on
those two conditional covariances join with only the universal factor 1/2.
This avoids multiplying the new translation result by the retained lifted
translation attenuation again.

The attitude/bias covariance ceilings are evaluated from each legal P2-V1
Pareto history first.  They depend on the same-history maximum pseudo cadence
(through the word horizon) and maximum q_c.  Only after evaluating each legal
history are the physical covariance results enveloped.  No Cartesian
(tau,sigma,R_S) extrema are reintroduced.

This is a canonical-P3 candidate, not an alternate definition: the unique
canonical gate still decides PASS from the frozen source/history, measurement,
H/A and 1e-18 obligations.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p2_correlation_path_memory as CORR
import ou3_p3_p2_v1_history_frontier as HIST
import ou3_p3_p2_v1_stage_phase_translation as TRANS
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
KALMAN = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
SCHEMA = 1
JOIN_FACTOR = 0.5


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def _source_contract(domain: dict) -> None:
    runtime = domain.get("configured_runtime", {})
    if runtime.get("imu_lever_arm_enabled") is not False:
        raise RuntimeError("full-state P3 requires declared zero-lever-arm branch")
    if runtime.get("accelerometer_vibration_guard_proof_branch") != "dormant_transparent":
        raise RuntimeError("full-state P3 requires dormant transparent vibration guard")
    text = KALMAN.read_text(encoding="utf-8")
    markers = (
        "const Matrix3 J_att = -skew_symmetric_matrix(f_cog_b);",
        "const Matrix3 J_aw  =  R_wb();",
        "constexpr int off_S = OFF_S;",
        "joseph_update3_(K, S_mat, PCt);",
    )
    for marker in markers:
        if marker not in text:
            raise RuntimeError(f"shipping measurement source marker changed: {marker}")


def _common_blocks(domain: dict) -> dict:
    live = domain["normal_live"]
    vector = BASE.VECTOR.build()
    process = BASE.PROCESS.build()
    for label, failures in (
        ("vector", BASE.VECTOR.validate(vector)),
        ("process", BASE.PROCESS.validate(process)),
    ):
        if failures:
            raise RuntimeError(f"{label} prerequisite failed: {failures}")
    sched = BASE.source_schedule()
    alpha6 = BASE.pos(BASE.vector_alpha6(live, vector), "declared alpha6")
    return {
        "live": live,
        "vector": vector,
        "process": process,
        "sched": sched,
        "alpha6": alpha6,
    }


def _history_word(summary: dict, domain: dict, sched: dict) -> dict:
    Tpe = BASE.pos(domain["normal_live"]["vector_pe_recurrence_window_s"], "PE recurrence")
    h = float(sched["dt_s"])
    cadence_hi = float(summary["pseudo_update_cadence_s"][1])
    qc = float(summary["q_c_upper"])
    if not all(math.isfinite(x) and x > 0.0 for x in (Tpe, h, cadence_hi, qc)):
        raise RuntimeError("same-history H/A word statistics lost positivity")
    gap = BASE.up(cadence_hi + h)
    spacing = BASE.up(max(Tpe, 2.0 * gap))
    Tobs = BASE.up(2.0 * spacing + gap)
    Tword = BASE.up(Tobs + Tpe)
    if float(summary["history_duration_s"][0]) < Tword:
        raise RuntimeError("same-history H/A summary is shorter than its covariance word")
    return {
        "Tpe_s": Tpe,
        "gap_s_upper": gap,
        "observation_horizon_s_upper": Tobs,
        "word_horizon_s_upper": Tword,
        "q_c_upper": qc,
    }


def _history_bias_upper(summary: dict, mode: str, domain: dict, blocks: dict) -> dict:
    if summary.get("all_statistics_from_one_legal_P2_history") is not True:
        raise RuntimeError("H/A covariance upper requires one legal P2 history")
    if summary.get("independent_global_source_extrema_used") is not False:
        raise RuntimeError("H/A covariance upper forbids independent global source extrema")
    live = blocks["live"]
    vector = blocks["vector"]
    process = blocks["process"]
    sched = blocks["sched"]
    alpha6 = float(blocks["alpha6"])
    word = _history_word(summary, domain, sched)

    vc = vector["configured_measurement_bounds"]
    ra = BASE.down(BASE.pos(vc["acc_measurement_std_mps2"], "acc std") ** 2)
    rm = BASE.down(BASE.pos(vc["mag_measurement_std_uT"], "mag std") ** 2)
    pc = process["source_constants"]
    qg = BASE.down(BASE.pos(pc["gyro_noise_density_rad_sqrt_s"], "gyro noise") ** 2)
    qb = BASE.down(BASE.pos(pc["gyro_bias_rw_variance_density"], "gyro bias process"))
    qba = BASE.down(BASE.pos(pc["accel_bias_process_variance_density"], "acc bias process"))
    tau_ba = BASE.pos(pc["accel_bias_tau_s"], "acc bias tau")
    pba = BASE.up(max(0.004 ** 2, qba * tau_ba / 2.0))
    fhi = BASE.pos(live["specific_force_norm_upper_mps2"], "force upper")
    mhi = BASE.pos(live["magnetic_vector_norm_upper_uT"], "mag upper")
    pair = BASE.pos(vector["operating_envelope"]["packet_gap_s"][1], "vector packet gap")
    qc = float(word["q_c_upper"])

    # This is the retained source-reachable attitude/bias covariance ceiling,
    # but with q_c and T taken from one legal history rather than a Cartesian
    # source cell/global box.
    qab = BASE.up(3.0 * (qg * pair + qb * (pair + pair ** 3 / 3.0)))
    whitened = BASE.up(
        (fhi * fhi / ra + mhi * mhi / rm) * qab
        + (3.0 * qc * pair) / ra
        + (3.0 * qba * pair) / ra
        + (6.0 * pba) / ra
    )
    u0 = BASE.up((1.0 + whitened) / alpha6)
    T = float(word["word_horizon_s_upper"])
    uab_prop = BASE.up(6.0 * (qg * T + qb * (T + T ** 3 / 3.0)))
    utheta = BASE.up(2.0 * (1.0 + T * T) * u0 + uab_prop)
    ubg = BASE.up(2.0 * u0 + uab_prop)
    return {
        "mode": mode,
        "theta_covariance_upper": utheta,
        "gyro_bias_covariance_upper": ubg,
        "accel_bias_covariance_upper": pba if mode == "A" else None,
        "word": word,
        "same_history_q_c_and_cadence_used": True,
    }


def _conditional_bias_floor(mode: str, domain: dict, blocks: dict) -> dict:
    live = blocks["live"]
    vector = blocks["vector"]
    process = blocks["process"]
    vc = vector["configured_measurement_bounds"]
    ra = BASE.down(BASE.pos(vc["acc_measurement_std_mps2"], "acc std") ** 2)
    rm = BASE.down(BASE.pos(vc["mag_measurement_std_uT"], "mag std") ** 2)
    fhi = BASE.pos(live["specific_force_norm_upper_mps2"], "force upper")
    mhi = BASE.pos(live["magnetic_vector_norm_upper_uT"], "mag upper")

    qtheta = BASE.pos(process["attitude_gyro_bias"]["theta_diagonal_lower"], "theta process")
    qbg = BASE.pos(process["attitude_gyro_bias"]["gyro_bias_diagonal_lower"], "bias process")
    cross = float(process["attitude_gyro_bias"]["cross_norm_upper"])
    rho_att = BASE.down(1.0 - cross / math.sqrt(qtheta * qbg))
    if rho_att <= 0.0:
        raise RuntimeError("attitude/gyro-bias process comparison lost positivity")
    ltheta = BASE.down(rho_att * qtheta)
    lbg = BASE.down(rho_att * qbg)
    qba = None
    if mode == "A":
        qba = BASE.pos(
            process["active_accelerometer_bias"]["Q_accel_bias_lambda_min_lower"],
            "active accel-bias discrete process",
        )

    # Selected B modes are injected on the final prediction.  Conditioning T
    # known removes translation columns.  S has no B column; only the final
    # accelerometer/magnetometer packets can attenuate these fresh modes.
    beta_theta_acc = (I(3.0) * I(ltheta) * I(fhi).square() / I(ra)).hi
    beta_theta_mag = (I(3.0) * I(ltheta) * I(mhi).square() / I(rm)).hi
    beta_ba_acc = 0.0
    if qba is not None:
        beta_ba_acc = (I(3.0) * I(qba) / I(ra)).hi
    beta = (I(beta_theta_acc) + I(beta_theta_mag) + I(beta_ba_acc)).hi
    attenuation = (I(1.0) / (I(1.0) + I(beta))).lo
    if attenuation <= 0.0:
        raise RuntimeError("conditional H/A measurement attenuation lost positivity")
    return {
        "mode": mode,
        "attitude_process_lower": ltheta,
        "gyro_bias_process_lower": lbg,
        "accel_bias_process_lower": qba,
        "beta_attitude_acc_upper": beta_theta_acc,
        "beta_attitude_mag_upper": beta_theta_mag,
        "beta_accel_bias_acc_upper": beta_ba_acc,
        "beta_total_upper": beta,
        "conditional_measurement_attenuation_lower": attenuation,
        "attitude_conditional_posterior_lower": BASE.down(attenuation * ltheta),
        "gyro_bias_conditional_posterior_lower": BASE.down(attenuation * lbg),
        "accel_bias_conditional_posterior_lower": (
            None if qba is None else BASE.down(attenuation * qba)
        ),
        "translation_columns_conditioned_known": True,
        "fresh_final_prediction_modes_only": True,
    }


def _mapped_labels(fr: dict, endpoint: int, positive_phase: bool) -> set[tuple[int, int, int, int]]:
    labels = HIST.endpoint_labels(fr, int(endpoint))
    out: set[tuple[int, int, int, int]] = set()
    rank = fr["stats"]["node_ranks"][int(endpoint)]
    for label in labels:
        q = HIST.update_label(label, rank) if positive_phase else label
        HIST.pareto_insert(out, q)
    return out


def _uniform_bias_upper(mode: str, domain: dict, blocks: dict, fr: dict) -> dict:
    max_theta = 0.0
    max_bg = 0.0
    max_ba = 0.0
    rows = 0
    limiting_theta = None
    limiting_bg = None
    for endpoint in fr["endpoint_nodes"]:
        for positive in (False, True):
            labels = _mapped_labels(fr, int(endpoint), positive)
            if not labels:
                raise RuntimeError("H/A source-history envelope lost endpoint labels")
            for label in labels:
                summary = HIST.label_summary(label, fr)
                u = _history_bias_upper(summary, mode, domain, blocks)
                rows += 1
                if u["theta_covariance_upper"] > max_theta:
                    max_theta = float(u["theta_covariance_upper"])
                    limiting_theta = {"endpoint": int(endpoint), "positive_phase": positive, "label": list(label)}
                if u["gyro_bias_covariance_upper"] > max_bg:
                    max_bg = float(u["gyro_bias_covariance_upper"])
                    limiting_bg = {"endpoint": int(endpoint), "positive_phase": positive, "label": list(label)}
                if mode == "A":
                    max_ba = max(max_ba, float(u["accel_bias_covariance_upper"]))
    if not (max_theta > 0.0 and max_bg > 0.0 and (mode == "H" or max_ba > 0.0)):
        raise RuntimeError("same-history H/A covariance envelope is not positive")
    return {
        "mode": mode,
        "same_history_rows_evaluated": rows,
        "theta_covariance_upper": max_theta,
        "gyro_bias_covariance_upper": max_bg,
        "accel_bias_covariance_upper": None if mode == "H" else max_ba,
        "limiting_theta_history": limiting_theta,
        "limiting_gyro_bias_history": limiting_bg,
        "same_history_evaluated_before_uniform_envelope": True,
        "raw_tuner_cartesian_extrema_used": False,
    }


def _mode_join(mode: str, translation: dict, domain: dict, blocks: dict, fr: dict) -> dict:
    cond = _conditional_bias_floor(mode, domain, blocks)
    upper = _uniform_bias_upper(mode, domain, blocks, fr)
    ratios = {
        "attitude": BASE.down(
            float(cond["attitude_conditional_posterior_lower"]) / float(upper["theta_covariance_upper"])
        ),
        "gyro_bias": BASE.down(
            float(cond["gyro_bias_conditional_posterior_lower"]) / float(upper["gyro_bias_covariance_upper"])
        ),
    }
    if mode == "A":
        ratios["accel_bias"] = BASE.down(
            float(cond["accel_bias_conditional_posterior_lower"])
            / float(upper["accel_bias_covariance_upper"])
        )
    bias_ratio = min(ratios.values())
    trans_ratio = float(translation["worst_translation_margin_lower"])
    prejoin = min(trans_ratio, bias_ratio)
    delta = BASE.down(JOIN_FACTOR * prejoin)
    return {
        "mode": mode,
        "dimension": 18 if mode == "H" else 21,
        "translation_conditional_relative_margin_lower": trans_ratio,
        "attitude_bias_conditional_relative_margins": ratios,
        "attitude_bias_conditional_relative_margin_lower": bias_ratio,
        "pre_precision_join_block_margin_lower": prejoin,
        "precision_block_join_factor": JOIN_FACTOR,
        "relative_Riccati_injection_margin_lower": delta,
        "useful_gate": BASE.MIN_USEFUL_DELTA,
        "useful_margin_established": delta >= BASE.MIN_USEFUL_DELTA,
        "conditional_bias_floor": cond,
        "same_history_bias_covariance_upper": upper,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, translation_candidate: dict | None = None) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("full-state P3 join must not be trajectory fitted")
    _source_contract(domain)
    if float(BASE.MIN_USEFUL_DELTA) != 1.0e-18:
        raise RuntimeError("canonical useful gate changed")

    translation = TRANS.build(path) if translation_candidate is None else translation_candidate
    tf = TRANS.validate(translation)
    if tf:
        raise RuntimeError(f"source-complete translation P3 invalid: {tf}")
    if translation.get("P2_correlation_interface_version") != CORR.INTERFACE_VERSION:
        raise RuntimeError("translation candidate lost frozen P2 V1 binding")

    fr = HIST.frontier_runtime(path)
    blocks = _common_blocks(domain)
    modes = {m: _mode_join(m, translation, domain, blocks, fr) for m in ("H", "A")}
    worst = min(float(modes[m]["relative_Riccati_injection_margin_lower"]) for m in ("H", "A"))
    passed = all(modes[m]["useful_margin_established"] for m in ("H", "A"))

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_P2_V1_SOURCE_COMPLETE_HA_PRECISION_BLOCK_JOIN",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "zero_lever_arm_branch": True,
        "dormant_transparent_vibration_guard_branch": True,
        "P2_correlation_interface_consumed": True,
        "P2_correlation_interface_version": CORR.INTERFACE_VERSION,
        "process_covariance_measurement_bounds_same_source_history": True,
        "independent_cartesian_tau_sigma_RS_extrema_used": False,
        "independent_cartesian_tau_sigma_R_S_extrema_used": False,
        "time_varying_tuner_over_word_covered": True,
        "interleaved_accelerometer_and_S_measurements_covered": True,
        "finite_clock_13_26_stage_language_covered": True,
        "frozen_clock_absorbing_hold_covered": translation.get("frozen_clock_absorbing_hold_branch_included") is True,
        "translation_full_matrix_samplewise_measurements_consumed": True,
        "translation_whole_word_lift_charged_again": False,
        "attitude_bias_fresh_final_prediction_modes_used": True,
        "conditional_precision_block_theorem_used": True,
        "precision_block_inequality": "J <= 2 diag(J_TT,J_BB) => J^-1 >= 0.5 diag(J_TT^-1,J_BB^-1)",
        "precision_block_join_factor": JOIN_FACTOR,
        "same_history_bias_upper_evaluated_before_uniform_envelope": True,
        "useful_gate": BASE.MIN_USEFUL_DELTA,
        "translation_subcertificate_established": translation.get("P3_TRANSLATION_CERTIFICATE_ESTABLISHED") is True,
        "translation_worst_margin_lower": float(translation["worst_translation_margin_lower"]),
        "modes": modes,
        "worst_H_A_relative_Riccati_injection_margin_lower": worst,
        "P3_PRODUCER_NUMERIC_PASS": passed,
        "P3_CANONICAL_PROMOTED": False,
        "P4_PROMOTED": False,
        "next_obligation": (
            "submit this exact producer artifact to the unique canonical P3 gate; only that gate may enable P4"
            if passed else
            "keep the frozen P2-V1/full-matrix/precision-block route and tighten only source-faithful numerical bounds; do not redefine P3 or relax 1e-18"
        ),
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P3_P2_V1_SOURCE_COMPLETE_HA_PRECISION_BLOCK_JOIN":
        f.append("wrong qualification")
    for key in (
        "source_generated_not_trajectory_fit", "zero_lever_arm_branch",
        "dormant_transparent_vibration_guard_branch", "P2_correlation_interface_consumed",
        "process_covariance_measurement_bounds_same_source_history",
        "time_varying_tuner_over_word_covered", "interleaved_accelerometer_and_S_measurements_covered",
        "finite_clock_13_26_stage_language_covered", "frozen_clock_absorbing_hold_covered",
        "translation_full_matrix_samplewise_measurements_consumed",
        "attitude_bias_fresh_final_prediction_modes_used", "conditional_precision_block_theorem_used",
        "same_history_bias_upper_evaluated_before_uniform_envelope",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_changed",
        "independent_cartesian_tau_sigma_RS_extrema_used",
        "independent_cartesian_tau_sigma_R_S_extrema_used",
        "translation_whole_word_lift_charged_again", "P3_CANONICAL_PROMOTED", "P4_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("P2_correlation_interface_version") != CORR.INTERFACE_VERSION:
        f.append("full-state P3 lost frozen P2 V1 binding")
    if float(d.get("precision_block_join_factor", math.nan)) != JOIN_FACTOR:
        f.append("precision-block join factor changed")
    if float(d.get("useful_gate", math.nan)) != 1.0e-18:
        f.append("full-state P3 useful gate changed")
    modes = d.get("modes", {})
    numeric_pass = True
    for mode in ("H", "A"):
        row = modes.get(mode, {})
        delta = row.get("relative_Riccati_injection_margin_lower")
        if not isinstance(delta, (int, float)) or isinstance(delta, bool) or not math.isfinite(float(delta)) or float(delta) < 0.0:
            f.append(f"{mode}: invalid full-state margin")
            numeric_pass = False
            continue
        expected = float(delta) >= BASE.MIN_USEFUL_DELTA
        if row.get("useful_margin_established") is not expected:
            f.append(f"{mode}: useful flag does not match unchanged gate")
        numeric_pass = numeric_pass and expected
    if d.get("P3_PRODUCER_NUMERIC_PASS") is not numeric_pass:
        f.append("producer numeric pass flag does not match H/A margins")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--translation-candidate", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    tc = json.loads(args.translation_candidate.read_text(encoding="utf-8")) if args.translation_candidate else None
    d = build(args.domain, tc)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "translation_margin": d["translation_worst_margin_lower"],
        "H_margin": d["modes"]["H"]["relative_Riccati_injection_margin_lower"],
        "A_margin": d["modes"]["A"]["relative_Riccati_injection_margin_lower"],
        "worst_H_A": d["worst_H_A_relative_Riccati_injection_margin_lower"],
        "producer_numeric_pass": d["P3_PRODUCER_NUMERIC_PASS"],
        "canonical_promoted_here": d["P3_CANONICAL_PROMOTED"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
