#!/usr/bin/env python3
"""Source-uniform post-measurement LTV P3 candidate for deployed OU-III.

This closes the measurement-attenuation gap left by
:mod:`ou3_p3_ltv_translation_ucc_probe` without multiplying a Joseph loss once
per sample.

The translation determinant used by that producer has an equivalent finite-mode
interpretation.  Split a suffix of length H into seven equal pieces of width
w=H/7 and retain one unit-L2 boxcar process mode on pieces 0,2,4,6 for each of
the three translation axes.  These twelve Gaussian coordinates are independent.
Multilinearity of the determinant gives

    det(C C^T) >= (225/16) q0^4 w^16 exp(-2 Lambda),

which is exactly the existing ``2025/144`` determinant constant.  Thus the
premeasurement LTV lower may be read as the endpoint covariance of a selected
12-dimensional process-noise subspace, not merely as an abstract Gramian bound.

Condition every unselected process input, the initial state and all nuisance
states as known.  Extra conditioning can only reduce Gaussian conditional
covariance, so the resulting selected-mode posterior is a lower bound on the
full Kalman/Joseph covariance.  If A is the complete measurement matrix seen by
those standardized selected modes, then

    Cov(zeta | y) = (I + A^T R^-1 A)^-1
                  >= I / (1 + beta),
    beta >= lambda_max(A^T R^-1 A).

We bound beta once by the trace over the whole suffix.  Every accepted
accelerometer packet is admitted, every possible S=0 packet is admitted, and
one final magnetometer packet is admitted.  Rejected/omitted packets only make
the true posterior larger.  The deployed CoG proof branch has no lever arm,
J_aw=R_wb is orthogonal, S=0 selects S exactly, and magnetometer has no
translation Jacobian.

The same conditioning argument includes a strict last-prediction attitude /
gyro-bias process block and, in A mode, the strict accelerometer-bias process
block.  This yields a full H/A linear P3 candidate against source-uniform
covariance uppers.  Promotion is numerical and fail-closed: the certificate is
PASS only when the worst H and A generalized margins exceed the unchanged
1e-18 usefulness gate.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_p3_ltv_translation_ucc_probe as LTV
import ou3_source_reachable_matrix_p3 as BASE

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
KALMAN = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
SCHEMA = 1


def point(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def _source_contract(domain: dict) -> None:
    runtime = domain.get("configured_runtime", {})
    if runtime.get("imu_lever_arm_enabled") is not False:
        raise RuntimeError("P3 lifted measurement bound requires the declared zero-lever-arm branch")
    if runtime.get("accelerometer_vibration_guard_proof_branch") != "dormant_transparent":
        raise RuntimeError("P3 lifted measurement bound requires the dormant transparent vibration-guard branch")
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


def _source_uniform_endpoint_blocks(domain: dict) -> dict:
    """Evaluate the retained P3 covariance upper on the complete source box."""
    live = domain["normal_live"]
    vector = BASE.VECTOR.build()
    process = BASE.PROCESS.build()
    sched = BASE.source_schedule()
    h = float(sched["dt_s"])
    tau_lo, tau_hi = map(float, sched["tau_applied_invariant_s"])
    sigma_lo, sigma_hi = map(float, sched["sigma_aw_applied_safety"])
    rs_lo, rs_hi = map(float, sched["R_S_applied_invariant"])
    x = Interval.outward_bounds(BASE.down(h / tau_hi), BASE.up(h / tau_lo))
    sigma = Interval.outward_bounds(sigma_lo, sigma_hi)
    rs = Interval.outward_bounds(rs_lo, rs_hi)
    alpha6 = BASE.vector_alpha6(live, vector)

    rows = {}
    for mode in ("H", "A"):
        # rho_trans=1 isolates the retained strict attitude/gyro-bias process
        # comparison while leaving every covariance-upper formula unchanged.
        raw = BASE.mode_cell(mode, x, 1.0, sigma, rs, live, vector, process, sched, alpha6)
        rows[mode] = raw
    return {"H": rows["H"], "A": rows["A"], "vector": vector, "process": process, "schedule": sched}


def _info_attenuation(probe: dict, mode: str, blocks: dict, domain: dict) -> dict:
    """One-shot information attenuation of the selected Gaussian noise modes."""
    H = float(probe["horizon_s"])
    q0 = float(probe["q_c_min_lower"])
    if not (H > 0.0 and q0 > 0.0):
        raise RuntimeError("strict LTV selected-mode horizon/process intensity required")

    sched = blocks["schedule"]
    vector = blocks["vector"]
    raw = blocks[mode]
    h = float(sched["dt_s"])
    # +2 covers an arbitrary left/right endpoint phase.  The deployed minimum
    # pseudo period is one IMU sample, so at most one S packet can occur per
    # valid sample; treating every slot as both S and accelerometer is stronger
    # than every source branch.
    packets = int(math.ceil(H / h)) + 2
    if float(sched["pseudo_min_s"]) + 1e-15 < h:
        raise RuntimeError("pseudo update can occur more than once per IMU sample")

    vc = vector["configured_measurement_bounds"]
    ra = BASE.down(BASE.pos(vc["acc_measurement_std_mps2"], "acc std") ** 2)
    rm = BASE.down(BASE.pos(vc["mag_measurement_std_uT"], "mag std") ** 2)
    rs = Interval(*map(float, sched["R_S_applied_invariant"]))
    rS_min = BASE.rs_variance_lower(rs, sched)
    if not (ra > 0.0 and rm > 0.0 and rS_min > 0.0):
        raise RuntimeError("measurement covariance lower lost positivity")

    Hiv = point(H)
    w = Hiv / point(7.0)
    selected_energy = point(q0) * w

    # Four boxcar modes x three axes.  For accelerometer, J_aw is orthogonal
    # and the damped a_w kernel has norm <=1, so each standardized mode has
    # squared sensitivity <= q0*w.
    beta_acc_translation = (
        point(float(packets)) * point(12.0) * selected_energy / point(ra)
    ).hi

    # For S, |Phi_Sa(t,s)| <= (t-s)^3/6 <= H^3/6.  Using the strongest S axis
    # for all three axes is conservative.  12*(H^6/36) = H^6/3.
    H6 = BASE.ipow(Hiv, 6)
    beta_S_translation = (
        point(float(packets)) * selected_energy * H6 / (point(3.0) * point(rS_min))
    ).hi

    scales = [float(x) for x in raw["comparison_scale_diagonal_squared"]]
    upper = [float(x) for x in raw["Sigma_diagonal_upper"]]
    rho_att = float(raw["process_scaled_lambda_min_lower"])
    qtheta = scales[0]
    qbg = scales[3]
    ltheta = BASE.down(rho_att * qtheta)
    lbg = BASE.down(rho_att * qbg)
    utheta = upper[0]
    ubg = upper[3]

    live = domain["normal_live"]
    fhi = BASE.pos(live["specific_force_norm_upper_mps2"], "force upper")
    mhi = BASE.pos(live["magnetic_vector_norm_upper_uT"], "mag upper")

    # Last-prediction attitude modes do not exist before the final packet.
    # ||skew(v)||_F^2 = 2||v||^2; factor 3 is a deliberately simple upper.
    beta_att_acc = (point(3.0) * point(ltheta) * point(fhi).square() / point(ra)).hi
    beta_att_mag = (point(3.0) * point(ltheta) * point(mhi).square() / point(rm)).hi

    beta_ba_acc = 0.0
    ba_ratio = math.inf
    qba = None
    uba = None
    if mode == "A":
        qba = scales[-1]
        uba = upper[-1]
        beta_ba_acc = (point(3.0) * point(qba) / point(ra)).hi
        ba_ratio = BASE.down(qba / uba)

    beta_total = (
        point(beta_acc_translation)
        + point(beta_S_translation)
        + point(beta_att_acc)
        + point(beta_att_mag)
        + point(beta_ba_acc)
    ).hi
    attenuation = (point(1.0) / (point(1.0) + point(beta_total))).lo
    if not attenuation > 0.0:
        raise RuntimeError("lifted measurement attenuation lost strict positivity")

    theta_ratio = BASE.down(ltheta / utheta)
    bg_ratio = BASE.down(lbg / ubg)
    translation_ratio = float(probe["relative_process_floor_lower"])
    block_ratio = min(theta_ratio, bg_ratio, translation_ratio, ba_ratio)
    delta = BASE.down(attenuation * block_ratio)

    return {
        "mode": mode,
        "horizon_s": H,
        "measurement_slots_upper": packets,
        "selected_translation_modes": 12,
        "selected_translation_boxcar_width_s": [w.lo, w.hi],
        "selected_base_process_intensity_lower": q0,
        "acc_measurement_variance_lower": ra,
        "mag_measurement_variance_lower": rm,
        "S_measurement_variance_lower": rS_min,
        "beta_acc_translation_upper": beta_acc_translation,
        "beta_S_translation_upper": beta_S_translation,
        "beta_final_attitude_acc_upper": beta_att_acc,
        "beta_final_attitude_mag_upper": beta_att_mag,
        "beta_final_accel_bias_acc_upper": beta_ba_acc,
        "beta_total_upper": beta_total,
        "joint_selected_mode_attenuation_lower": attenuation,
        "attitude_relative_noise_floor_lower": theta_ratio,
        "gyro_bias_relative_noise_floor_lower": bg_ratio,
        "translation_relative_noise_floor_before_measurements": translation_ratio,
        "accel_bias_relative_noise_floor_lower": None if mode == "H" else ba_ratio,
        "pre_attenuation_full_block_relative_floor_lower": block_ratio,
        "relative_Riccati_injection_margin_lower": delta,
        "useful_gate": BASE.MIN_USEFUL_DELTA,
        "useful": delta >= BASE.MIN_USEFUL_DELTA,
        "last_step_attitude_process_lower": {"theta": ltheta, "gyro_bias": lbg},
        "last_step_accel_bias_process_lower": qba,
        "covariance_upper": {"theta": utheta, "gyro_bias": ubg, "accel_bias": uba},
    }


def _mode_certificate(pre: dict, mode: str, blocks: dict, domain: dict) -> dict:
    rows = []
    for endpoint in pre["source_complete_translation"]["endpoint_rows"]:
        candidates = []
        for probe in endpoint["candidates"]:
            row = _info_attenuation(probe, mode, blocks, domain)
            row["endpoint_tau_index"] = int(endpoint["endpoint_tau_index"])
            row["clock_phase_decay_exponent_upper"] = float(probe["decay_exponent_upper"])
            candidates.append(row)
        best = max(candidates, key=lambda x: x["relative_Riccati_injection_margin_lower"])
        rows.append({
            "endpoint_tau_index": int(endpoint["endpoint_tau_index"]),
            "best": best,
            "candidates": candidates,
        })
    worst = min(rows, key=lambda x: x["best"]["relative_Riccati_injection_margin_lower"])
    delta = float(worst["best"]["relative_Riccati_injection_margin_lower"])
    return {
        "mode": mode,
        "endpoint_tau_cells_scanned": len(rows),
        "endpoint_rows": rows,
        "worst_endpoint": worst,
        "relative_Riccati_injection_margin_lower": delta,
        "useful_gate": BASE.MIN_USEFUL_DELTA,
        "useful_margin_established": delta >= BASE.MIN_USEFUL_DELTA,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("LTV postmeasurement P3 must not be trajectory fitted")
    _source_contract(domain)

    # Exact-node diagnostics are intentionally omitted here; only the active
    # source-complete LTV route is consumed.
    pre = LTV.build(path, source_node_indices=())
    pf = LTV.validate(pre)
    if pf:
        raise RuntimeError(f"premeasurement LTV producer failed: {pf}")
    blocks = _source_uniform_endpoint_blocks(domain)

    modes = {m: _mode_certificate(pre, m, blocks, domain) for m in ("H", "A")}
    worst = min(float(modes[m]["relative_Riccati_injection_margin_lower"]) for m in ("H", "A"))
    passed = all(modes[m]["useful_margin_established"] for m in ("H", "A"))

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P3_SOURCE_UNIFORM_LTV_POSTMEASUREMENT_CERTIFICATE",
        "source_generated_not_trajectory_fit": True,
        "linear_P3_only": True,
        "zero_lever_arm_branch": True,
        "dormant_transparent_vibration_guard_branch": True,
        "four_boxcar_modes_per_axis": True,
        "selected_translation_modes": 12,
        "finite_mode_determinant_identity": "2025/144 = (15/4)^2 = 225/16",
        "all_accelerometer_packets_admitted": True,
        "all_possible_S_packets_admitted": True,
        "final_magnetometer_packet_admitted": True,
        "rejected_or_missing_measurements_only_increase_covariance": True,
        "nuisance_states_and_unselected_process_conditioned_known": True,
        "measurement_attenuation_applied_once_in_lifted_information_space": True,
        "per_sample_multiplicative_Joseph_loss_used": False,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "useful_gate": BASE.MIN_USEFUL_DELTA,
        "modes": modes,
        "worst_H_A_relative_Riccati_injection_margin_lower": worst,
        "P3_LINEAR_CERTIFICATE_ESTABLISHED": passed,
        "P3_PROMOTED": passed,
        "next_obligation": (
            "if P3 passes, feed these source-uniform directional covariance/noise ratios into the retained finite-angle P4 complete-word route; "
            "P3 promotion is linear only and does not itself establish nonlinear P4 dissipation"
        ) if passed else (
            "tighten only source-faithful lifted information/covariance bounds; do not relax the 1e-18 usefulness gate"
        ),
    }


def validate(d: dict) -> list[str]:
    f = []
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P3_SOURCE_UNIFORM_LTV_POSTMEASUREMENT_CERTIFICATE":
        f.append("wrong qualification")
    for key in (
        "source_generated_not_trajectory_fit", "linear_P3_only", "zero_lever_arm_branch",
        "dormant_transparent_vibration_guard_branch", "four_boxcar_modes_per_axis",
        "all_accelerometer_packets_admitted", "all_possible_S_packets_admitted",
        "final_magnetometer_packet_admitted", "rejected_or_missing_measurements_only_increase_covariance",
        "nuisance_states_and_unselected_process_conditioned_known",
        "measurement_attenuation_applied_once_in_lifted_information_space",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in ("per_sample_multiplicative_Joseph_loss_used", "trajectory_replay_used", "filter_changed"):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if int(d.get("selected_translation_modes", 0)) != 12:
        f.append("selected translation mode count changed")
    for mode in ("H", "A"):
        row = d.get("modes", {}).get(mode, {})
        delta = row.get("relative_Riccati_injection_margin_lower")
        if not isinstance(delta, (int, float)) or not (math.isfinite(float(delta)) and float(delta) > 0.0):
            f.append(f"{mode}: postmeasurement margin is not strict")
        if int(row.get("endpoint_tau_cells_scanned", 0)) != 10:
            f.append(f"{mode}: did not scan ten endpoint tau cells")
    established = all(
        float(d.get("modes", {}).get(m, {}).get("relative_Riccati_injection_margin_lower", 0.0)) >= BASE.MIN_USEFUL_DELTA
        for m in ("H", "A")
    )
    if bool(d.get("P3_LINEAR_CERTIFICATE_ESTABLISHED")) != established:
        f.append("P3 establishment flag disagrees with numerical useful gate")
    if bool(d.get("P3_PROMOTED")) != established:
        f.append("P3 promotion flag disagrees with numerical useful gate")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "P3_LINEAR_CERTIFICATE_ESTABLISHED": d["P3_LINEAR_CERTIFICATE_ESTABLISHED"],
        "worst_H_A_delta": d["worst_H_A_relative_Riccati_injection_margin_lower"],
        "H": {
            "delta": d["modes"]["H"]["relative_Riccati_injection_margin_lower"],
            "worst_tau": d["modes"]["H"]["worst_endpoint"]["endpoint_tau_index"],
            "best_horizon_s": d["modes"]["H"]["worst_endpoint"]["best"]["horizon_s"],
            "attenuation": d["modes"]["H"]["worst_endpoint"]["best"]["joint_selected_mode_attenuation_lower"],
        },
        "A": {
            "delta": d["modes"]["A"]["relative_Riccati_injection_margin_lower"],
            "worst_tau": d["modes"]["A"]["worst_endpoint"]["endpoint_tau_index"],
            "best_horizon_s": d["modes"]["A"]["worst_endpoint"]["best"]["horizon_s"],
            "attenuation": d["modes"]["A"]["worst_endpoint"]["best"]["joint_selected_mode_attenuation_lower"],
        },
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
