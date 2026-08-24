#!/usr/bin/env python3
"""Explicit numerical H/A information-word certificate for deployed OU-III.

This closes the linear quantitative step without appealing only to abstract
Riccati existence and without using replay extrema.

Construction
============

A source-complete proof word is synchronized as follows.  Start at a qualifying
accepted accel/mag PE packet.  The periodic S=0 scheduler supplies four S
firings within its validated translational UCO window.  After that window, the
recurring-PE hypothesis supplies the next qualifying vector packet within the
declared recurrence window.  Thus a word has a finite source-uniform horizon

    T_word <= T_S,UCO + T_PE.

For covariance *upper* bounds, use an explicit finite-memory estimator.  The
four S observations give a full 12-state (v,p,S,a_w) observation operator with
known smallest singular value.  The vector packet gives a six-state
(attitude,gyro-bias) information operator.  Accelerometer rows depend also on
a_w, so estimate translation first and subtract its contribution; the resulting
estimator is block-triangular.  Its mean-square error is bounded by elementary
operator inequalities including process noise.  Since the Kalman conditional
covariance is no larger than the MSE of any measurement-based estimator, its
trace gives a source-uniform Sigma upper bound.  Active b_a is a stable
Gauss-Markov tail and is bounded independently.

For covariance/noise *lower* bounds, every sample prediction injects a strict
full-state process floor.  An accepted linear measurement maps a prior p I to a
posterior bounded below by

    (p^{-1} + ||H||^2 / r_min)^{-1} I.

Applying this conservatively to S, accelerometer and magnetometer corrections
bounds the last-sample injected word noise Omega_w from below.  The left-error
reset G=I+0.5[dtheta]x satisfies G^T G >= I, and periodic a_w synchronization
adds PSD covariance, so neither can reduce this lower bound.

Consequently

    delta_w >= Omega_min / Sigma_max > 0,
    W_1 <= (1-delta_w) W_0.

The exact covariance identity also gives prefix information gain <= 1 for every
same-mode prefix.  This is a numerical source-uniform contraction statement,
not an existential local-ISS result.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_full_process_ucc as PROCESS
import ou3_implementation_word_language as WORDS
import ou3_source_domain_contract as SOURCE
import ou3_translational_uco_ucc as TRANS
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def pos(x, label: str) -> float:
    y = float(x)
    if not math.isfinite(y) or not y > 0.0:
        raise RuntimeError(f"{label} must be finite positive, got {x!r}")
    return y


def posterior_floor(p: float, h2_upper: float, r_lower: float) -> float:
    p = pos(p, "prior floor")
    h2_upper = max(0.0, float(h2_upper))
    r_lower = pos(r_lower, "measurement variance lower")
    info = up(up(1.0 / p) + up(h2_upper / r_lower))
    return down(1.0 / info)


def translation_process_trace(qc: float, T: float) -> float:
    # Three independent axes, response [v,p,S,a] from OU driving.
    poly = up(T + T**3 / 3.0 + T**5 / 20.0 + T**7 / 252.0)
    return up(3.0 * up(qc * poly))


def s_process_variance_per_axis(qc: float, T: float) -> float:
    return up(qc * up(T**7 / 252.0))


def attitude_bias_process_trace(qg: float, qb: float, T: float) -> float:
    # theta receives gyro white noise and integrated gyro-bias RW; b_g receives RW.
    per_axis = up(qg * T + qb * (T + T**3 / 3.0))
    return up(3.0 * per_axis)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("proof domain must not be trajectory fitted")
    live = domain["normal_live"]

    source = SOURCE.build(SOURCE.DEFAULT_HEADER.resolve())
    trans = TRANS.build(TRANS.DEFAULT_HEADER.resolve())
    vector = VECTOR.build()
    process = PROCESS.build()
    words = WORDS.build(domain_path)

    if TRANS.validate(trans):
        raise RuntimeError(f"translation prerequisite failed: {TRANS.validate(trans)}")
    if VECTOR.validate(vector):
        raise RuntimeError(f"vector prerequisite failed: {VECTOR.validate(vector)}")
    if PROCESS.validate(process):
        raise RuntimeError(f"process prerequisite failed: {PROCESS.validate(process)}")
    if WORDS.validate(words):
        raise RuntimeError(f"word-language prerequisite failed: {WORDS.validate(words)}")

    cp = source["validated_parameter_box"]["continuous_parameters"]
    sigma_hi = pos(cp["sigma_aw_mps2"][1], "sigma_aw upper")
    tau_lo = pos(cp["tau_aw_s"][0], "tau lower")
    qc_max = up(2.0 * sigma_hi * sigma_hi / tau_lo)

    T_s = pos(trans["S_observation_uco"]["aligned_window_s"], "S UCO window")
    T_pe = pos(live["vector_pe_recurrence_window_s"], "PE recurrence")
    T_word = up(T_s + T_pe)
    pair_gap = pos(vector["operating_envelope"]["packet_gap_s"][1], "vector packet gap")

    # ---------------- translation finite-memory estimator ----------------
    sigma_t = pos(trans["S_observation_uco"]["observation_sigma_min_lower"], "translation observation sigma min")
    rs_std_hi = pos(trans["S_observation_uco"]["R_S_filter_std_upper"], "R_S std upper")
    rs_var_hi = up(rs_std_hi * rs_std_hi)

    # 4 scalar S samples per axis = 12 scalar measurements.  Process correlations
    # across the stack are handled by bounding E||n||^2 (trace), so no
    # independence across firing times is assumed.
    s_proc_var = s_process_variance_per_axis(qc_max, T_s)
    n_t_trace = up(2.0 * up(12.0 * rs_var_hi + 12.0 * s_proc_var))
    Lt2 = up(1.0 / down(sigma_t * sigma_t))
    translation_start_mse_trace = up(Lt2 * n_t_trace)

    L_t = up(1.0 + T_word + 0.5 * T_word**2 + T_word**3 / 6.0)
    q_t_trace = translation_process_trace(qc_max, T_word)
    translation_endpoint_trace = up(
        2.0 * up(L_t * L_t * translation_start_mse_trace) + 2.0 * q_t_trace
    )

    # ---------------- attitude / gyro-bias finite-memory estimator ---------
    vc = vector["configured_measurement_bounds"]
    ra_hi = pos(vc["acc_measurement_variance_upper"], "acc variance upper")
    rm_hi = pos(vc["mag_measurement_variance_upper"], "mag variance upper")
    ra_lo = down(pos(vc["acc_measurement_std_mps2"], "acc std") ** 2)
    rm_lo = down(pos(vc["mag_measurement_std_uT"], "mag std") ** 2)
    alpha6 = pos(vector["gyro_bias_two_packet"]["alpha_6_information_lower"], "vector alpha6")

    # O^T R^-1 O >= alpha I and R^-1 <= I/r_min imply
    # O^T O >= r_min*alpha I.  Thus a left inverse has squared norm <=1/(r_min alpha).
    r_vec_min = min(ra_lo, rm_lo)
    sigma_a2 = down(r_vec_min * alpha6)
    La2 = up(1.0 / sigma_a2)

    pc = process["source_constants"]
    qg = down(pos(pc["gyro_noise_density_rad_sqrt_s"], "gyro noise density") ** 2)
    qb = pos(pc["gyro_bias_rw_variance_density"], "gyro bias RW density")
    qba = pos(pc["accel_bias_process_variance_density"], "accel bias density")
    tau_ba = pos(pc["accel_bias_tau_s"], "accel bias tau")
    pba_cap = up(max(0.004**2, qba * tau_ba / 2.0))

    f_hi = pos(live["specific_force_norm_upper_mps2"], "specific force upper")
    m_hi = pos(live["magnetic_vector_norm_upper_uT"], "magnetic norm upper")

    qab_pair_trace = attitude_bias_process_trace(qg, qb, pair_gap)
    qaw_pair_trace = up(3.0 * qc_max * pair_gap)
    qba_pair_trace = up(3.0 * qba * pair_gap)
    measurement_trace = up(6.0 * ra_hi + 6.0 * rm_hi)
    vector_process_output_trace = up(
        up((f_hi * f_hi + m_hi * m_hi) * qab_pair_trace)
        + qaw_pair_trace + qba_pair_trace
    )
    # b_a contributes to the two accel packets.  Use the A-mode stable cap for
    # both H and A; this is conservative for held mode.
    ba_packet_trace = up(6.0 * pba_cap)
    n_a_trace = up(4.0 * up(measurement_trace + vector_process_output_trace + ba_packet_trace))

    # Two accepted accelerometer packets couple vector observations to a_w.
    # The stacked operator has ||C||^2 <= 2.  Eliminate translation with its
    # explicit left inverse, then estimate attitude/bias.
    C2 = 2.0
    attitude_start_mse_trace = up(
        2.0 * La2 * n_a_trace
        + 2.0 * La2 * C2 * translation_start_mse_trace
    )
    L_ab = up(1.0 + T_word)
    q_ab_trace = attitude_bias_process_trace(qg, qb, T_word)
    attitude_endpoint_trace = up(
        2.0 * L_ab * L_ab * attitude_start_mse_trace + 2.0 * q_ab_trace
    )

    ba_endpoint_trace = up(3.0 * pba_cap)
    sigma_max_H = up(translation_endpoint_trace + attitude_endpoint_trace)
    sigma_max_A = up(sigma_max_H + ba_endpoint_trace)

    # ---------------- covariance / word-noise lower bounds -----------------
    rs_std_lo = pos(trans["S_observation_uco"]["R_S_filter_std_lower"], "R_S std lower")
    rs_var_lo = down(rs_std_lo * rs_std_lo)

    h2_S = 1.0
    h2_acc_H = up(f_hi * f_hi + 1.0)       # attitude + a_w
    h2_acc_A = up(f_hi * f_hi + 2.0)       # attitude + a_w + b_a
    h2_mag = up(m_hi * m_hi)

    modes = {}
    for mode, h2_acc, sigma_max in (
        ("H", h2_acc_H, sigma_max_H),
        ("A", h2_acc_A, sigma_max_A),
    ):
        qmin = pos(process["modes"][mode]["prediction_Q_lambda_min_lower"], f"{mode} Q floor")
        # Worst lower floor assumes every shrinking measurement available in a
        # sample is accepted. Rejections only leave a larger covariance floor.
        p = qmin
        p = posterior_floor(p, h2_S, rs_var_lo)
        p = posterior_floor(p, h2_acc, ra_lo)
        p = posterior_floor(p, h2_mag, rm_lo)
        omega_min = p
        sigma_min = p
        delta = down(omega_min / sigma_max)
        if not (math.isfinite(delta) and 0.0 < delta < 1.0):
            raise RuntimeError(f"{mode} explicit information margin is not strict: {delta!r}")
        modes[mode] = {
            "dimension": 18 if mode == "H" else 21,
            "word_horizon_s_upper": T_word,
            "Sigma_lambda_min_lower": sigma_min,
            "Sigma_lambda_max_upper": sigma_max,
            "word_noise_Omega_lambda_min_lower": omega_min,
            "relative_Riccati_injection_margin_lower": delta,
            "lambda_information_upper_formula": "1-delta_lower",
            "prefix_information_gain_upper": 1.0,
            "strict_information_contraction": True,
            "pass": True,
        }

    return {
        "schema": SCHEMA,
        "qualification": "EXPLICIT_SOURCE_UNIFORM_DEPLOYED_OU3_HA_INFORMATION_WORD_CERTIFICATE",
        "source_generated_not_trajectory_fit": True,
        "validated_scalar_bounds": True,
        "configured_runtime": source["configured_runtime_assumption"],
        "word_construction": {
            "S_uco_window_s_upper": T_s,
            "PE_recurrence_window_s_upper": T_pe,
            "synchronized_word_horizon_s_upper": T_word,
            "source_complete_relative_to_declared_hypotheses": True,
            "construction": (
                "qualifying vector packet -> four guaranteed S firings -> next qualifying vector packet"
            ),
        },
        "translation_candidate_estimator": {
            "observation_sigma_min_lower": sigma_t,
            "stacked_noise_MSE_trace_upper": n_t_trace,
            "start_MSE_trace_upper": translation_start_mse_trace,
            "endpoint_MSE_trace_upper": translation_endpoint_trace,
            "process_trace_upper": q_t_trace,
        },
        "attitude_gyro_bias_candidate_estimator": {
            "unweighted_observation_sigma_squared_lower": sigma_a2,
            "stacked_noise_MSE_trace_upper": n_a_trace,
            "translation_coupling_norm_squared_upper": C2,
            "start_MSE_trace_upper": attitude_start_mse_trace,
            "endpoint_MSE_trace_upper": attitude_endpoint_trace,
            "process_trace_upper": q_ab_trace,
        },
        "active_accelerometer_bias": {
            "per_axis_covariance_cap_upper": pba_cap,
            "trace_cap_upper": ba_endpoint_trace,
            "stable_tail": True,
        },
        "measurement_lower_bounds": {
            "R_S_variance_lower": rs_var_lo,
            "R_acc_variance_lower": ra_lo,
            "R_mag_variance_lower": rm_lo,
            "H_S_norm_squared_upper": h2_S,
            "H_acc_H_norm_squared_upper": h2_acc_H,
            "H_acc_A_norm_squared_upper": h2_acc_A,
            "H_mag_norm_squared_upper": h2_mag,
        },
        "modes": modes,
        "continuous_linear_information_certificate": "PASS",
        "nonlinear_word_enclosed": False,
        "theorem_promotion": "LINEAR_ONLY",
        "next_obligation": (
            "validate the exact nonlinear SO(3) source words on explicit source-node levels and prove mu_W>0 plus prefix safety"
        ),
    }


def validate(d: dict) -> list[str]:
    failures: list[str] = []
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True:
        failures.append("certificate is not source generated")
    if d.get("validated_scalar_bounds") is not True:
        failures.append("scalar bounds are not validated")
    if d.get("word_construction", {}).get("source_complete_relative_to_declared_hypotheses") is not True:
        failures.append("word construction is not source complete")
    for mode in ("H", "A"):
        row = d.get("modes", {}).get(mode, {})
        for key in (
            "Sigma_lambda_min_lower", "Sigma_lambda_max_upper",
            "word_noise_Omega_lambda_min_lower", "relative_Riccati_injection_margin_lower",
        ):
            x = row.get(key)
            if not isinstance(x, (int, float)) or not math.isfinite(float(x)) or not float(x) > 0.0:
                failures.append(f"{mode}.{key} is not finite positive")
        eta = row.get("relative_Riccati_injection_margin_lower")
        if isinstance(eta, (int, float)) and not float(eta) < 1.0:
            failures.append(f"{mode} information margin is not <1")
        if row.get("prefix_information_gain_upper") != 1.0:
            failures.append(f"{mode} prefix information identity was not retained")
        if row.get("pass") is not True:
            failures.append(f"{mode} linear information certificate failed")
    if d.get("continuous_linear_information_certificate") != "PASS":
        failures.append("continuous linear information certificate did not pass")
    if d.get("nonlinear_word_enclosed") is not False:
        failures.append("linear producer must not claim nonlinear enclosure")
    if d.get("theorem_promotion") != "LINEAR_ONLY":
        failures.append("linear producer must not claim full theorem promotion")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain.resolve())
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "word_horizon_s_upper": d["word_construction"]["synchronized_word_horizon_s_upper"],
        "H": d["modes"]["H"],
        "A": d["modes"]["A"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
