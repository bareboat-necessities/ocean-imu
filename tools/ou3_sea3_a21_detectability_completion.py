#!/usr/bin/env python3
"""Quantitative A21 finite-bias detectability for canonical complete SEA3.

The OU-III paper permits two active accelerometer-bias routes: full eta9 PE, or
eta6 PE together with the implemented finite residual-bias correlation time.
The canonical deployment uses the latter here.  This module deliberately does
NOT invent an eta9 packet condition.

The construction is source-bound to ``COMPLETE_SEA3_NORMAL_LIVE_WORD``:

* the H18 active coordinates use the already-certified complete-SEA3 H18 word
  contraction at delta_H = 1e-18;
* the residual accelerometer bias has exact homogeneous Gauss--Markov decay
  exp(-T/tau_b) over the same 3 s word;
* while held, b_a is decoupled, has identity homogeneous dynamics, no process
  injection, zero covariance cross terms and frozen measurement rows;
* H->A is a one-time bounded source transition: enabling only floors the b_a
  diagonal before exact GM prediction/correction resumes;
* the active A21 process covariance is strictly positive on all 21 coordinates.

For detectability, take any stabilizing complete-word observer for the H18
subsystem supplied by the certified H word.  Extend its output injection by a
zero b_a row.  Because the physical homogeneous prediction is block diagonal
between H18 and b_a, the resulting A21 comparison-observer word map is block
upper triangular,

    E_A = [[E_H, C_Hb],
           [  0, Phi_b]],

where C_Hb is finite on the compact 3 s SEA3 word.  Its diagonal blocks have
strict energy gaps delta_H and 1-Phi_b^2 respectively.  Since Phi_b is much
more contractive than E_H, the finite upper-right block changes only the UES
prefactor, not the asymptotic word rate.  Thus the detectability route retains
an asymptotic energy gap min(delta_H, 1-Phi_b^2) = delta_H.

A deliberately crude finite numerical bound on the bias->H coupling is also
reported in log10 form.  It uses only the source-uniform H covariance ceiling,
configured accelerometer covariance, bounded specific force, and the finite
number of IMU events in the word.  This bound is not used as a contraction
shortcut; it certifies that the triangular comparison observer has a finite
mode-coupling constant.

Together with the shipping A21 process UCC, this closes the paper's finite-tau_b
A21 *detectability/UES hypothesis*.  It does not by itself assert the stronger
canonical implementation-word matrix inequality Omega-delta*P >= 0.  That
bridge remains explicit and fail-closed below.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_full_process_ucc as PROCESS
import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_full_normal_live_word as WORD
import ou3_sea3_h18_prior_free_completion as H18
import ou3_sea3_live_covariance_seed as LIVE
import ou3_sea3_windowed_vector_pe as PE

DEFAULT_DOMAIN = COMPLETE.DEFAULT_DOMAIN
SCHEMA = 1
QUALIFICATION = "OU3_COMPLETE_SEA3_A21_FINITE_BIAS_DETECTABILITY"
USEFUL_GATE = 1.0e-18
HORIZON_S = 3.0


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    complete = COMPLETE.build(path)
    word = WORD.build(path)
    h18 = H18.build(path)
    live = LIVE.build(path)
    pe = PE.build(path)
    process = PROCESS.build()
    bad = {
        "complete": COMPLETE.validate(complete),
        "literal_word": WORD.validate(word),
        "H18": H18.validate(h18),
        "Live_seed": LIVE.validate(live),
        "PE": PE.validate(pe),
        "process": PROCESS.validate(process),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"A21 detectability prerequisites failed: {bad}")
    source = complete["canonical_P3_source"]
    if source != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        raise RuntimeError("A21 detectability detached from complete SEA3")
    if float(complete["word_horizon_s"]) != HORIZON_S:
        raise RuntimeError("A21 detectability no longer uses the canonical 3 s word")
    if not h18["H18_prior_free_completion_closed"]:
        raise RuntimeError("A21 detectability requires the certified H18 complete-word contraction")

    route = pe["A_mode_bias_route"]
    if route["uses_eta9_pointwise_packet_shortcut"] is not False:
        raise RuntimeError("eta9 point-packet shortcut re-entered A21")
    if route["uses_eta6_plus_finite_bias_correlation"] is not True:
        raise RuntimeError("paper finite-bias detectability route is not active")

    phi_b = float(route["homogeneous_bias_contraction_upper_over_word"])
    bias_gap = float(route["homogeneous_bias_contraction_gap_lower"])
    if not (0.0 < phi_b < 1.0 and bias_gap > 0.0):
        raise RuntimeError("active b_a homogeneous decay is not strict")
    bias_energy_gap = down(1.0 - up(phi_b * phi_b))
    if not bias_energy_gap > 0.0:
        raise RuntimeError("active b_a energy contraction is not strict")

    held = live["held_ba"]
    release = live["H_to_A_release"]
    if not all(bool(held[k]) for k in (
        "excluded_from_H18", "identity_homogeneous_dynamics",
        "no_process_injection_while_held", "cross_covariances_zero",
        "measurement_rows_frozen",
    )):
        raise RuntimeError("shipping held-bias semantics changed")
    seed_var = float(held["seed_variance"])
    floor_var = float(release["bias_diagonal_floor_variance"])
    if not (seed_var > 0.0 and floor_var > 0.0 and floor_var == seed_var):
        raise RuntimeError("H->A release floor is not the shipping b_a seed")

    pc = process["source_constants"]
    qba_density = float(pc["accel_bias_process_variance_density"])
    tau_ba = float(pc["accel_bias_tau_s"])
    qba_d = float(process["active_accelerometer_bias"]["Q_accel_bias_lambda_min_lower"])
    pba_stationary = up(qba_density * tau_ba / 2.0)
    pba_uniform_upper = up(max(seed_var, pba_stationary))
    if not (qba_density > 0.0 and tau_ba > 0.0 and qba_d > 0.0 and pba_uniform_upper > 0.0):
        raise RuntimeError("active b_a process/release bound lost positivity")

    # A very conservative finite comparison-observer cross-coupling bound.
    # This is used only to prove finiteness of the triangular UES prefactor.
    pbar_trace = float(h18["same_word_diffuse_prior_covariance_upper"]["Pbar_trace_upper"])
    domain = json.loads(path.read_text(encoding="utf-8"))
    fmax = float(domain["normal_live"]["specific_force_norm_upper_mps2"])
    racc_std = min(map(float, pe["measurement_runtime"]["accelerometer_std_mps2"]))
    rmin = down(racc_std * racc_std)
    Hacc_norm = up(math.sqrt(fmax * fmax + 1.0))
    Kacc_norm = up(pbar_trace * Hacc_norm / rmin)
    Aacc_norm = up(1.0 + Kacc_norm * Hacc_norm)
    Fstep_norm = up(math.sqrt(float(h18["prediction_F_spectral_norm_squared_upper"])))
    step_norm = up(max(1.0, Fstep_norm) * max(1.0, Aacc_norm))
    samples = int(word["imu_samples_upper"])
    if samples <= 0 or not (step_norm > 1.0 and math.isfinite(step_norm)):
        raise RuntimeError("cannot establish finite A21 comparison-observer coupling bound")
    log10_step = math.log10(step_norm)
    # Sum of at most N injected bias terms through at most N bounded H maps.
    # log10(N*K*L^N) is finite and avoids intentionally gigantic binary64 powers.
    log10_cross = up(
        math.log10(max(1, samples))
        + math.log10(max(1.0, Kacc_norm))
        + samples * log10_step
    )

    h_gap = float(h18["useful_gate"])
    if h_gap != USEFUL_GATE:
        raise RuntimeError("H18 useful gap changed")
    asymptotic_gap = min(h_gap, bias_energy_gap)
    detectability_closed = (
        asymptotic_gap >= USEFUL_GATE
        and math.isfinite(log10_cross)
        and bool(process["modes"]["A"]["pass"])
        and float(process["modes"]["A"]["prediction_Q_lambda_min_lower"]) > 0.0
    )

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": source,
        "complete_word_horizon_s": HORIZON_S,
        "component_of_complete_SEA3_full_word": True,
        "paper_active_bias_route": "ETA6_PLUS_FINITE_RESIDUAL_BIAS_CORRELATION",
        "eta9_point_packet_shortcut_used": False,
        "H18_complete_word_contraction_consumed": True,
        "H18_word_energy_gap_lower": h_gap,
        "bias_homogeneous_contraction_upper_over_word": phi_b,
        "bias_homogeneous_contraction_gap_lower": bias_gap,
        "bias_homogeneous_energy_gap_lower": bias_energy_gap,
        "A21_detectability_asymptotic_word_energy_gap_lower": asymptotic_gap,
        "A21_detectability_useful_gate": USEFUL_GATE,
        "A21_detectability_useful_gate_pass": detectability_closed,
        "held_bias_semantics": held,
        "H_to_A_release": release,
        "H_to_A_release_is_bounded_one_time_mode_jump": True,
        "H_to_A_release_requires_preceding_three_second_H_interval": False,
        "active_bias_process": {
            "tau_ba_s": tau_ba,
            "continuous_driving_variance_density": qba_density,
            "one_sample_Q_ba_lambda_min_lower": qba_d,
            "stationary_variance_upper": pba_stationary,
            "uniform_release_and_active_variance_upper": pba_uniform_upper,
            "full_A21_process_UCC_pass": bool(process["modes"]["A"]["pass"]),
        },
        "triangular_detectability_observer": {
            "form": "[[E_H,C_Hb],[0,Phi_b]]",
            "H_diagonal_block_uses_complete_SEA3_H18_certificate": True,
            "bias_diagonal_block_uses_shipping_GM_decay": True,
            "upper_right_coupling_finite_on_compact_word": True,
            "finite_coupling_log10_upper": log10_cross,
            "finite_coupling_changes_prefactor_not_asymptotic_rate": True,
            "comparison_observer_only_not_alternate_estimator": True,
            "shipping_filter_changed": False,
        },
        "A21_finite_bias_detectability_closed": detectability_closed,
        "A21_paper_UES_hypotheses_closed": detectability_closed,
        "full_21x21_Omega_minus_delta_P_LDLT_closed_here": False,
        "P3_CANONICAL_PASS": False,
        "P4_MAY_CONSUME_P3": False,
        "source_family_replaced": False,
        "trajectory_replay_used": False,
        "independent_tau_sigma_RS_source_created": False,
        "next_obligation": (
            "bridge this certified complete-SEA3 finite-tau_b detectability route to the canonical full A21 "
            "Riccati word comparison (or make the canonical gate consume the theorem-equivalent quantitative "
            "detectability/UCC certificate) without weakening the 1e-18 requirement"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for key in (
        "component_of_complete_SEA3_full_word",
        "H18_complete_word_contraction_consumed",
        "A21_detectability_useful_gate_pass",
        "H_to_A_release_is_bounded_one_time_mode_jump",
        "A21_finite_bias_detectability_closed",
        "A21_paper_UES_hypotheses_closed",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "eta9_point_packet_shortcut_used",
        "H_to_A_release_requires_preceding_three_second_H_interval",
        "full_21x21_Omega_minus_delta_P_LDLT_closed_here",
        "P3_CANONICAL_PASS", "P4_MAY_CONSUME_P3",
        "source_family_replaced", "trajectory_replay_used",
        "independent_tau_sigma_RS_source_created",
    ):
        if d.get(key) is not False:
            f.append(f"forbidden/open flag {key} changed")
    if d.get("paper_active_bias_route") != "ETA6_PLUS_FINITE_RESIDUAL_BIAS_CORRELATION":
        f.append("wrong active-bias theorem route")
    if float(d.get("A21_detectability_useful_gate", math.nan)) != USEFUL_GATE:
        f.append("A21 detectability gate changed")
    if float(d.get("A21_detectability_asymptotic_word_energy_gap_lower", 0.0)) < USEFUL_GATE:
        f.append("A21 detectability energy gap is below useful gate")
    tri = d.get("triangular_detectability_observer", {})
    for key in (
        "H_diagonal_block_uses_complete_SEA3_H18_certificate",
        "bias_diagonal_block_uses_shipping_GM_decay",
        "upper_right_coupling_finite_on_compact_word",
        "finite_coupling_changes_prefactor_not_asymptotic_rate",
        "comparison_observer_only_not_alternate_estimator",
    ):
        if tri.get(key) is not True:
            f.append(f"triangular detectability property lost: {key}")
    if tri.get("shipping_filter_changed") is not False:
        f.append("detectability comparison changed shipping filter")
    x = tri.get("finite_coupling_log10_upper")
    if not isinstance(x, (int, float)) or not math.isfinite(float(x)):
        f.append("finite bias-to-H coupling bound is not finite")
    proc = d.get("active_bias_process", {})
    for key in (
        "tau_ba_s", "continuous_driving_variance_density",
        "one_sample_Q_ba_lambda_min_lower", "uniform_release_and_active_variance_upper",
    ):
        x = proc.get(key)
        if not isinstance(x, (int, float)) or not (math.isfinite(float(x)) and float(x) > 0.0):
            f.append(f"invalid A bias process field {key}")
    if proc.get("full_A21_process_UCC_pass") is not True:
        f.append("full A21 process UCC not retained")
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
        "route": d["paper_active_bias_route"],
        "H_gap": d["H18_word_energy_gap_lower"],
        "bias_energy_gap": d["bias_homogeneous_energy_gap_lower"],
        "A21_asymptotic_gap": d["A21_detectability_asymptotic_word_energy_gap_lower"],
        "finite_coupling_log10_upper": d["triangular_detectability_observer"]["finite_coupling_log10_upper"],
        "A21_detectability_closed": d["A21_finite_bias_detectability_closed"],
        "full_A21_Riccati_matrix_closed": d["full_21x21_Omega_minus_delta_P_LDLT_closed_here"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
