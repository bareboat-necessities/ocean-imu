#!/usr/bin/env python3
"""Explicit numerical H/A information-word certificate for deployed OU-III.

This producer supplies the quantitative linear step of the implementation proof
without replay promotion.  The construction avoids global physical-unit box
conditioning in the candidate-estimator upper bound:

* S=0 observations are selected from well-separated guaranteed firing windows;
* the three-firing (v,p,S) inverse is enclosed directly with validated interval
  Gauss-Jordan arithmetic instead of det/Frobenius singular-value relaxation;
* the full-heading information bound is evaluated on the declared physical PE
  deployment envelope rather than on near-zero generic guard values;
* the vector-UCO constant is used in its native noise-weighted information
  metric, so measurement noise is not counted twice;
* only a_w, not the complete translational state, is charged as nuisance in the
  accelerometer vector observation.

The remaining scalar covariance/noise lower bound is deliberately reported
separately.  If it is still the active conditioning bottleneck, P4 must replace
that scalar floor with a block/source-scaled relative bound rather than hide it
behind a tiny local-neighborhood claim.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, up as ia_up
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan
import ou3_full_process_ucc as PROCESS
import ou3_implementation_word_language as WORDS
import ou3_source_domain_contract as SOURCE
import ou3_translational_uco_ucc as TRANS
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 3


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
    # Three independent axes, response [v,p,S,a] from OU driving.  Replacing
    # exponential decay by one only enlarges every impulse-response component.
    poly = up(T + T**3 / 3.0 + T**5 / 20.0 + T**7 / 252.0)
    return up(3.0 * up(qc * poly))


def s_process_variance_per_axis(qc: float, T: float) -> float:
    return up(qc * up(T**7 / 252.0))


def attitude_bias_process_trace(qg: float, qb: float, T: float) -> float:
    per_axis = up(qg * T + qb * (T + T**3 / 3.0))
    return up(3.0 * per_axis)


def interval_frobenius_upper(A) -> float:
    total = 0.0
    for row in A:
        for x in row:
            a = x.abs_upper()
            total = ia_up(total + ia_up(a * a))
    return ia_up(math.sqrt(total))


def validated_integrator_inverse_bound(gap_max: float, spacing: float) -> dict:
    """Enclose B^-1 for all three source-reachable selected S firing times.

    The selected firings lie in windows
        t0 in [0,g], t1 in [s,s+g], t2 in [2s,2s+g].
    Use column order [S,p,v], so each row is [1,t,t^2/2].  Fixed-pivot interval
    elimination is safe because the first pivot is exactly one and the later
    divided-difference pivots remain strictly separated by s-g>0.
    """
    g = pos(gap_max, "pseudo gap max")
    s = pos(spacing, "S window spacing")
    if not s > g:
        raise RuntimeError("selected S windows are not strictly separated")
    windows = (
        Interval.outward_bounds(0.0, g),
        Interval.outward_bounds(s, s + g),
        Interval.outward_bounds(2.0 * s, 2.0 * s + g),
    )
    half = Interval.point(0.5)
    one = Interval.point(1.0)
    B = [[one, t, half * t.square()] for t in windows]
    Binv = matrix_inverse_gauss_jordan(B)
    inv_frob = interval_frobenius_upper(Binv)
    sigma_min = down(1.0 / inv_frob)
    if not sigma_min > 0.0:
        raise RuntimeError("validated three-firing integrator inverse lost rank")
    return {
        "time_windows_s": [t.as_list() for t in windows],
        "inverse_frobenius_upper": inv_frob,
        "observation_sigma_min_lower": sigma_min,
        "inverse_interval": [[[x.lo, x.hi] for x in row] for row in Binv],
        "backend": "OUTWARD_INTERVAL_FIXED_PIVOT_GAUSS_JORDAN",
    }


def declared_vector_information(live: dict, vector: dict) -> dict:
    """Re-evaluate the vector-UCO lemma on the declared deployment PE box."""
    base = vector["operating_envelope"]
    for name in (
        "specific_force_norm_lower_mps2",
        "magnetic_vector_norm_lower_uT",
        "vector_sine_separation_lower",
        "body_rate_norm_upper_deg_s",
    ):
        if name not in live:
            raise KeyError(name)

    f_min = pos(live["specific_force_norm_lower_mps2"], "declared force floor")
    m_min = pos(live["magnetic_vector_norm_lower_uT"], "declared magnetic floor")
    s = pos(live["vector_sine_separation_lower"], "declared vector sine floor")
    omega_deg = pos(live["body_rate_norm_upper_deg_s"], "declared body-rate ceiling")
    if not s < 1.0:
        raise RuntimeError("declared vector sine floor must be <1")

    if f_min < float(base["specific_force_norm_lower_mps2"]):
        raise RuntimeError("declared force floor is weaker than vector-UCO source floor")
    if m_min < float(base["magnetic_vector_norm_lower_uT"]):
        raise RuntimeError("declared magnetic floor is weaker than vector-UCO source floor")
    if s < float(base["vector_sine_separation_lower"]):
        raise RuntimeError("declared vector separation is weaker than vector-UCO source floor")
    if omega_deg > float(base["body_rate_norm_upper_deg_s"]):
        raise RuntimeError("declared body-rate ceiling exceeds vector-UCO source ceiling")

    vc = vector["configured_measurement_bounds"]
    ra = up(pos(vc["acc_measurement_variance_upper"], "acc variance upper"))
    rm = up(pos(vc["mag_measurement_variance_upper"], "mag variance upper"))
    af = down(f_min * f_min / ra)
    am = down(m_min * m_min / rm)
    root = up(math.sqrt(max(0.0, 1.0 - s * s)))
    angular = down((s * s) / up(1.0 + root))
    mu = down(min(af, am) * angular)

    gap_lo, gap_hi = map(float, base["packet_gap_s"])
    omega = up(omega_deg * math.pi / 180.0)
    bracket = down(1.0 - up(0.5 * omega * gap_hi))
    g_min = down(gap_lo * bracket)
    Tbg = pos(base["gyro_bias_time_scale_s"], "gyro-bias time scale")
    gamma = down(g_min / Tbg)
    alpha6 = down(mu / up(1.0 + up(2.0 / down(gamma * gamma))))
    if not alpha6 > 0.0:
        raise RuntimeError("declared vector information lower bound is not strict")
    return {
        "specific_force_norm_lower_mps2": f_min,
        "magnetic_vector_norm_lower_uT": m_min,
        "vector_sine_separation_lower": s,
        "body_rate_norm_upper_deg_s": omega_deg,
        "angular_factor_lower": angular,
        "mu_theta_lower": mu,
        "alpha_6_information_lower": alpha6,
        "packet_gap_s": [gap_lo, gap_hi],
        "source": "declared theorem operating domain; not replay fitted",
    }


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

    for label, failures in (
        ("translation", TRANS.validate(trans)),
        ("vector", VECTOR.validate(vector)),
        ("process", PROCESS.validate(process)),
        ("word-language", WORDS.validate(words)),
    ):
        if failures:
            raise RuntimeError(f"{label} prerequisite failed: {failures}")

    cp = source["validated_parameter_box"]["continuous_parameters"]
    sigma_hi = pos(cp["sigma_aw_mps2"][1], "sigma_aw upper")
    tau_lo = pos(cp["tau_aw_s"][0], "tau lower")
    qc_max = up(2.0 * sigma_hi * sigma_hi / tau_lo)

    T_pe = pos(live["vector_pe_recurrence_window_s"], "PE recurrence")
    pseudo_gap_max = pos(
        trans["S_observation_uco"]["pseudo_gap_max_s"], "pseudo firing gap upper"
    )
    spacing = up(max(T_pe, 2.0 * pseudo_gap_max))
    selected_gap = down(spacing - pseudo_gap_max)
    if not selected_gap > 0.0:
        raise RuntimeError("cannot construct separated S firing windows")
    T_obs = up(2.0 * spacing + pseudo_gap_max)
    T_word = up(T_obs + T_pe)

    # ---------------- (v,p,S) finite-memory estimator ----------------------
    inv_cert = validated_integrator_inverse_bound(pseudo_gap_max, spacing)
    sigma_int = pos(inv_cert["observation_sigma_min_lower"], "integrator sigma min")
    Lint2 = up(1.0 / down(sigma_int * sigma_int))

    rs_std_hi = pos(
        trans["S_observation_uco"]["R_S_filter_std_upper"], "R_S std upper"
    )
    rs_var_hi = up(rs_std_hi * rs_std_hi)
    phiS_hi = up(T_obs**3 / 6.0)
    aw_nuisance_trace = up(9.0 * sigma_hi * sigma_hi * phiS_hi * phiS_hi)
    s_proc_var = s_process_variance_per_axis(qc_max, T_obs)
    measurement_trace = up(9.0 * rs_var_hi)
    process_output_trace = up(9.0 * s_proc_var)
    stacked_noise_trace = up(
        measurement_trace + aw_nuisance_trace + process_output_trace
    )
    integrator_start_mse = up(Lint2 * stacked_noise_trace)

    L_integrator = up(1.0 + T_word + 0.5 * T_word * T_word)
    aw_response_energy = up(
        1.0 + T_word**2 + T_word**4 / 4.0 + T_word**6 / 36.0
    )
    aw_endpoint_trace = up(3.0 * sigma_hi * sigma_hi * aw_response_energy)
    q_t_trace = translation_process_trace(qc_max, T_word)
    translation_endpoint_trace = up(
        2.0 * L_integrator * L_integrator * integrator_start_mse
        + 2.0 * up(aw_endpoint_trace + q_t_trace)
    )
    aw_cov_trace_upper = up(3.0 * sigma_hi * sigma_hi)

    # ---------------- attitude / gyro-bias finite-memory estimator ---------
    vec_decl = declared_vector_information(live, vector)
    vc = vector["configured_measurement_bounds"]
    ra_hi = pos(vc["acc_measurement_variance_upper"], "acc variance upper")
    rm_hi = pos(vc["mag_measurement_variance_upper"], "mag variance upper")
    ra_lo = down(pos(vc["acc_measurement_std_mps2"], "acc std") ** 2)
    rm_lo = down(pos(vc["mag_measurement_std_uT"], "mag std") ** 2)
    alpha6 = pos(vec_decl["alpha_6_information_lower"], "declared vector alpha6")

    pc = process["source_constants"]
    qg = down(pos(pc["gyro_noise_density_rad_sqrt_s"], "gyro noise density") ** 2)
    qb = pos(pc["gyro_bias_rw_variance_density"], "gyro bias RW density")
    qba = pos(pc["accel_bias_process_variance_density"], "accel bias density")
    tau_ba = pos(pc["accel_bias_tau_s"], "accel bias tau")
    pba_cap = up(max(0.004**2, qba * tau_ba / 2.0))

    f_hi = pos(live["specific_force_norm_upper_mps2"], "specific force upper")
    m_hi = pos(live["magnetic_vector_norm_upper_uT"], "magnetic norm upper")
    pair_gap = pos(vec_decl["packet_gap_s"][1], "vector packet gap")
    qab_pair_trace = attitude_bias_process_trace(qg, qb, pair_gap)
    qaw_pair_trace = up(3.0 * qc_max * pair_gap)
    qba_pair_trace = up(3.0 * qba * pair_gap)

    # alpha6 is already lambda_min(O^T R^-1 O), so weighted least squares has
    # measurement-noise covariance <= alpha6^-1 I.  Additional process/nuisance
    # outputs are charged in the same whitened measurement coordinates.
    whitened_process_trace = up(
        up((f_hi * f_hi / ra_lo + m_hi * m_hi / rm_lo) * qab_pair_trace)
        + up(qaw_pair_trace / ra_lo)
        + up(qba_pair_trace / ra_lo)
        + up(6.0 * pba_cap / ra_lo)
    )
    measurement_noise_state_trace = up(6.0 / alpha6)
    nuisance_state_trace = up(whitened_process_trace / alpha6)
    attitude_start_mse = up(measurement_noise_state_trace + nuisance_state_trace)

    L_ab = up(1.0 + T_word)
    q_ab_trace = attitude_bias_process_trace(qg, qb, T_word)
    attitude_endpoint_trace = up(
        2.0 * L_ab * L_ab * attitude_start_mse + 2.0 * q_ab_trace
    )

    ba_endpoint_trace = up(3.0 * pba_cap)
    sigma_max_H = up(translation_endpoint_trace + attitude_endpoint_trace)
    sigma_max_A = up(sigma_max_H + ba_endpoint_trace)

    # ---------------- conservative covariance/noise lower bounds -----------
    rs_std_lo = pos(
        trans["S_observation_uco"]["R_S_filter_std_lower"], "R_S std lower"
    )
    rs_var_lo = down(rs_std_lo * rs_std_lo)
    h2_S = 1.0
    h2_acc_H = up(f_hi * f_hi + 1.0)
    h2_acc_A = up(f_hi * f_hi + 2.0)
    h2_mag = up(m_hi * m_hi)

    modes = {}
    for mode, h2_acc, sigma_max in (
        ("H", h2_acc_H, sigma_max_H),
        ("A", h2_acc_A, sigma_max_A),
    ):
        qmin = pos(
            process["modes"][mode]["prediction_Q_lambda_min_lower"], f"{mode} Q floor"
        )
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
            "scalar_lower_bound_is_p4_conditioning_bottleneck": delta < 1.0e-12,
            "pass": True,
        }

    return {
        "schema": SCHEMA,
        "qualification": "EXPLICIT_SOURCE_UNIFORM_DEPLOYED_OU3_HA_INFORMATION_WORD_CERTIFICATE",
        "source_generated_not_trajectory_fit": True,
        "validated_scalar_bounds": True,
        "configured_runtime": source["configured_runtime_assumption"],
        "word_construction": {
            "pseudo_gap_max_s": pseudo_gap_max,
            "selected_S_window_spacing_s": spacing,
            "selected_S_gap_lower_s": selected_gap,
            "selected_S_observation_horizon_s_upper": T_obs,
            "PE_recurrence_window_s_upper": T_pe,
            "synchronized_word_horizon_s_upper": T_word,
            "source_complete_relative_to_declared_hypotheses": True,
            "construction": (
                "qualifying vector packet -> three well-separated guaranteed S firings -> next qualifying vector packet"
            ),
        },
        "declared_vector_PE_information": vec_decl,
        "translation_candidate_estimator": {
            "state_order": ["v", "p", "S"],
            "validated_observation_inverse": inv_cert,
            "selected_observation_sigma_min_lower": sigma_int,
            "stacked_measurement_trace_upper": measurement_trace,
            "initial_aw_nuisance_trace_upper": aw_nuisance_trace,
            "process_output_trace_upper": process_output_trace,
            "start_MSE_trace_upper": integrator_start_mse,
            "endpoint_MSE_trace_upper": translation_endpoint_trace,
            "endpoint_aw_trace_component_upper": aw_endpoint_trace,
            "process_trace_upper": q_t_trace,
            "route": "validated three-firing inverse + detectable (v,p,S) + stable a_w tail",
        },
        "attitude_gyro_bias_candidate_estimator": {
            "weighted_information_lambda_min_lower": alpha6,
            "measurement_noise_state_trace_upper": measurement_noise_state_trace,
            "whitened_process_and_nuisance_trace_upper": whitened_process_trace,
            "nuisance_state_trace_upper": nuisance_state_trace,
            "aw_only_translation_coupling_trace_upper": aw_cov_trace_upper,
            "start_MSE_trace_upper": attitude_start_mse,
            "endpoint_MSE_trace_upper": attitude_endpoint_trace,
            "process_trace_upper": q_ab_trace,
            "route": "noise-weighted vector UCO; no second measurement-variance scaling",
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
            "replace the isolated one-sample scalar covariance floor with a block/source-scaled relative word-noise bound, then use that metric in the exact SO(3) word majorant"
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
    if not float(d.get("word_construction", {}).get("selected_S_gap_lower_s", 0.0)) > 0.0:
        failures.append("selected S firing separation is not strict")
    pe = d.get("declared_vector_PE_information", {})
    if not float(pe.get("alpha_6_information_lower", 0.0)) > 0.0:
        failures.append("declared vector PE information is not strict")
    inv = d.get("translation_candidate_estimator", {}).get("validated_observation_inverse", {})
    if inv.get("backend") != "OUTWARD_INTERVAL_FIXED_PIVOT_GAUSS_JORDAN":
        failures.append("integrator observation inverse is not validated")
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
        "word_construction": d["word_construction"],
        "declared_vector_PE_information": d["declared_vector_PE_information"],
        "translation_candidate_estimator": d["translation_candidate_estimator"],
        "attitude_gyro_bias_candidate_estimator": d["attitude_gyro_bias_candidate_estimator"],
        "H": d["modes"]["H"],
        "A": d["modes"]["A"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
