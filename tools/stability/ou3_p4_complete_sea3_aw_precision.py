#!/usr/bin/env python3
"""Source-uniform full-posterior a_w precision bound for complete SEA3 P4.

This module closes one narrow P4 obligation without constructing another source
family.  It consumes the same source-complete H/A process cells used by P3 and
bounds the *a_w block of the full posterior inverse*, not the inverse of an a_w
marginal covariance.

For every fixed-dimensional Normal-Live IMU packet, shipping predicts first:

    P^- = F P F^T + Q >= Q > 0.

Hence, in Loewner order,

    (P^-)^-1 <= Q^-1.

The translational process block is the exact correlated [v,p,S,a_w] integrated
OU covariance.  We invert its full 4x4 one-axis matrix cell-by-cell over the
same h/tau source interval used by the canonical P3 backend.  If D is the
physical diagonal similarity used there,

    Q_axis = D Q_scaled D,       D_aa = sigma_aw,

so the full-inverse a_w entry obeys

    e_a^T Q_axis^-1 e_a
       <= sup (Q_scaled^-1)_{aa} / sigma_aw,min^2.

The shipping operation order then makes the block bound especially simple:

* the pending a_w covariance floor is a PSD increment, so it can only decrease
  full precision;
* an S=0 Joseph update adds H_S^T R_S^-1 H_S to full information, whose a_w-a_w
  block is exactly zero;
* the accepted accelerometer update adds H_a^T R_acc^-1 H_a, and its a_w column
  is the orthogonal R_hat, so its a_w-a_w block is at most R_acc,min^-1;
* magnetometer information has zero a_w column;
* every immediate left-error reset is block diagonal with identity on a_w, so
  E_aw P_reset^-1 E_aw^T is exactly unchanged by the reset congruence.

Therefore this is a bound on E_aw P_J^-1 E_aw^T for the actual full H18/A21
posterior, with all cross-covariances retained implicitly through Q^-1.  No
marginal inverse, selected replay, packet-count remainder, independent R_S
schedule, or alternate estimator is used.

The accelerometer-noise value is the configured deployment value already bound
into canonical conditional P3.  The theorem branch requires the vibration guard
to be dormant/transparent, so the enabled vibration-aware covariance path does
not lower that covariance; if this proof branch changes, validation fails.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan
import ou3_source_reachable_matrix_p3 as P3M
import ou3_vector_uco_certificate as VECTOR
import ou3_sea3_riccati_metric_p3 as P3

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
SCHEMA = 1
QUALIFICATION = "OU3_P4_COMPLETE_SEA3_FULL_POSTERIOR_AW_PRECISION_V1"


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _inverse_aw_entry_upper(cell: Interval, depth: int = 0) -> tuple[float, int]:
    """Outward upper for (Q_scaled^-1)[a_w,a_w] on one h/tau cell."""
    try:
        inv = matrix_inverse_gauss_jordan(P3M.step_scaled_q(cell))
        value = float(inv[3][3].hi)
        if not (math.isfinite(value) and value > 0.0):
            raise RuntimeError("nonpositive/nonfinite inverse a_w entry")
        return up(value), 1
    except Exception:
        if depth >= 18:
            raise
        mid = math.sqrt(cell.lo * cell.hi)
        if not (cell.lo < mid < cell.hi):
            raise RuntimeError("cannot subdivide process inverse cell")
        left, nl = _inverse_aw_entry_upper(
            Interval.outward_bounds(cell.lo, mid), depth + 1
        )
        right, nr = _inverse_aw_entry_upper(
            Interval.outward_bounds(mid, cell.hi), depth + 1
        )
        return up(max(left, right)), nl + nr


def _source_contract(domain: dict) -> list[str]:
    failures: list[str] = []
    live = domain.get("normal_live", {})
    runtime = domain.get("configured_runtime", {})
    if runtime.get("accelerometer_vibration_guard_proof_branch") != "dormant_transparent":
        failures.append("a_w precision proof requires dormant transparent vibration-guard branch")
    if live.get("accelerometer_update_on_every_valid_imu_sample") is not True:
        failures.append("Normal-Live accelerometer update-on-every-valid-sample contract changed")
    if live.get("accelerometer_rejection_after_live_allowed") is not False:
        failures.append("Normal-Live theorem unexpectedly admits accelerometer rejection")

    w = WRAPPER.read_text(encoding="utf-8")
    m = MEKF.read_text(encoding="utf-8")
    wrapper_markers = (
        "apply_pending_online_tune_();",
        "if (drive_mekf) apply_racc_vibration_inflation_();",
        "mekf_->time_update(gyro, dt);",
        "mekf_->measurement_update_acc_only(acc_in, tempC);",
        "const Eigen::Vector3f base = racc_base_std_();",
        "const float excess = accel_guard_.excessRms();",
        "if (!(excess > 0.0f))",
    )
    mekf_markers = (
        "const Matrix3 J_aw  =  R_wb();",
        "applyIntegralZeroPseudoMeas",
        "measurement_update_mag_only",
        "apply_error_state_reset_jacobian_(dtheta);",
    )
    for marker in wrapper_markers:
        if marker not in w:
            failures.append(f"shipping wrapper precision semantic changed: {marker}")
    for marker in mekf_markers:
        if marker not in m:
            failures.append(f"shipping MEKF precision semantic changed: {marker}")
    return failures


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("a_w precision proof must not be trajectory fitted")

    p3 = P3.build(path)
    vector = VECTOR.build()
    failures = P3.validate(p3) + VECTOR.validate(vector) + _source_contract(domain)
    if failures or p3.get("P3_CONDITIONAL_SEA3_PASS") is not True:
        raise RuntimeError(f"unchanged P3/configured-vector prerequisites failed: {failures}")

    sched = P3M.source_schedule()
    h = float(sched["dt_s"])
    tau_lo, tau_hi = map(float, sched["tau_applied_invariant_s"])
    sigma_lo, sigma_hi = map(float, sched["sigma_aw_applied_safety"])
    if not (0.0 < h and 0.0 < tau_lo <= tau_hi and 0.0 < sigma_lo <= sigma_hi):
        raise RuntimeError("source schedule lost positive h/tau/sigma bounds")

    xlo = down(h / tau_hi)
    xhi = up(h / tau_lo)
    edges = P3M.geom_edges(xlo, xhi, 48)
    if xlo < P3M.BRANCH_X < xhi:
        edges = sorted(set(edges + [P3M.BRANCH_X]))

    qscaled_aw_inv_upper = 0.0
    inverse_cells = 0
    for cell in P3M.interval_cells(edges):
        bound, count = _inverse_aw_entry_upper(cell)
        qscaled_aw_inv_upper = max(qscaled_aw_inv_upper, bound)
        inverse_cells += count
    qscaled_aw_inv_upper = up(qscaled_aw_inv_upper)

    # D_aa = sigma_aw in the P3 physical similarity.  The smallest admitted
    # sigma gives the largest physical inverse a_w information.
    prediction_aw_precision_diag_upper = up(
        qscaled_aw_inv_upper / down(sigma_lo * sigma_lo)
    )

    vc = vector["configured_measurement_bounds"]
    racc_std = float(vc["acc_measurement_std_mps2"])
    racc_var_lower = down(racc_std * racc_std)
    if not (math.isfinite(racc_var_lower) and racc_var_lower > 0.0):
        raise RuntimeError("configured accelerometer covariance lost positive lower bound")
    accel_aw_information_diag_upper = up(1.0 / racc_var_lower)

    posterior_aw_precision_diag_upper = up(
        prediction_aw_precision_diag_upper + accel_aw_information_diag_upper
    )
    posterior_aw_precision_trace_upper = up(3.0 * posterior_aw_precision_diag_upper)

    modes = {}
    for mode, dim in (("H", 18), ("A", 21)):
        modes[mode] = {
            "dimension": dim,
            "prediction_full_inverse_aw_diag_upper": prediction_aw_precision_diag_upper,
            "S_zero_added_aw_information_diag_upper": 0.0,
            "accelerometer_added_aw_information_diag_upper": accel_aw_information_diag_upper,
            "magnetometer_added_aw_information_diag_upper": 0.0,
            "reset_changes_full_inverse_aw_block": False,
            "full_posterior_aw_precision_diag_upper": posterior_aw_precision_diag_upper,
            "full_posterior_aw_precision_trace_upper": posterior_aw_precision_trace_upper,
            "source_uniform_full_posterior_aw_precision_closed": True,
        }

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "P3_frozen_not_modified": True,
        "P3_conditional_complete_SEA3_consumed": True,
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "source_family_replaced": False,
        "independent_parameter_sequence_used": False,
        "independent_RS_schedule_used": False,
        "marginal_covariance_inverse_used": False,
        "full_process_cross_covariances_retained": True,
        "packet_count_multiplier_used": False,
        "configured_runtime_dt_s": h,
        "tau_applied_invariant_s": [tau_lo, tau_hi],
        "sigma_aw_applied_safety_mps2": [sigma_lo, sigma_hi],
        "x_h_over_tau": [xlo, xhi],
        "validated_process_inverse_cells": inverse_cells,
        "scaled_process_inverse_aw_diag_upper": qscaled_aw_inv_upper,
        "configured_accelerometer_std_mps2": racc_std,
        "configured_accelerometer_variance_lower": racc_var_lower,
        "proof_identity": {
            "prediction": "Pminus>=Q => Pminus^-1<=Q^-1",
            "S_zero": "P_J^-1=P^-1+H_S^T R_S^-1 H_S; E_aw H_S^T=0",
            "accelerometer": "E_aw H_a^T Racc^-1 H_a E_aw^T=R_hat^T Racc^-1 R_hat",
            "magnetometer": "E_aw H_m^T=0",
            "left_reset": "G_ext is identity on a_w, so E_aw G_ext^-T J G_ext^-1 E_aw^T=E_aw J E_aw^T",
            "aw_floor": "Pplus=Pminus+E_aw Delta_+ E_aw^T >= Pminus, so full precision cannot increase",
        },
        "modes": modes,
        "source_uniform_full_posterior_aw_precision_closed": True,
        "source_uniform_operation_covariance_ceiling_closed_here": False,
        "source_correlated_correction_radius_closed_here": False,
        "full_epsilon_aw_transport_closed_here": False,
        "complete_word_nonlinear_dissipation_closed_here": False,
        "P4_promoted_here": False,
        "next_obligation": (
            "combine this full-posterior a_w precision bound with a reset/floor-complete "
            "source-uniform correction covariance/radius bound; then transport the full epsilon_aw "
            "through the same complete word"
        ),
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for key in (
        "P3_frozen_not_modified", "P3_conditional_complete_SEA3_consumed",
        "source_generated_not_trajectory_fit", "full_process_cross_covariances_retained",
        "source_uniform_full_posterior_aw_precision_closed",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "source_family_replaced", "independent_parameter_sequence_used",
        "independent_RS_schedule_used", "marginal_covariance_inverse_used",
        "packet_count_multiplier_used", "source_uniform_operation_covariance_ceiling_closed_here",
        "source_correlated_correction_radius_closed_here", "full_epsilon_aw_transport_closed_here",
        "complete_word_nonlinear_dissipation_closed_here", "P4_promoted_here",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if not isinstance(d.get("validated_process_inverse_cells"), int) or d["validated_process_inverse_cells"] <= 0:
        f.append("no validated process inverse cells")
    for key in (
        "scaled_process_inverse_aw_diag_upper",
        "configured_accelerometer_variance_lower",
    ):
        x = d.get(key)
        if not isinstance(x, (int, float)) or not (math.isfinite(float(x)) and float(x) > 0.0):
            f.append(f"invalid positive scalar {key}")
    for mode, dim in (("H", 18), ("A", 21)):
        row = d.get("modes", {}).get(mode, {})
        if row.get("dimension") != dim:
            f.append(f"{mode} dimension changed")
        if row.get("source_uniform_full_posterior_aw_precision_closed") is not True:
            f.append(f"{mode} full posterior a_w precision not closed")
        if row.get("reset_changes_full_inverse_aw_block") is not False:
            f.append(f"{mode} reset incorrectly changes a_w precision block")
        for key in (
            "prediction_full_inverse_aw_diag_upper",
            "accelerometer_added_aw_information_diag_upper",
            "full_posterior_aw_precision_diag_upper",
            "full_posterior_aw_precision_trace_upper",
        ):
            x = row.get(key)
            if not isinstance(x, (int, float)) or not (math.isfinite(float(x)) and float(x) > 0.0):
                f.append(f"{mode}.{key} invalid")
        if row.get("S_zero_added_aw_information_diag_upper") != 0.0:
            f.append(f"{mode} S=0 falsely adds direct a_w information")
        if row.get("magnetometer_added_aw_information_diag_upper") != 0.0:
            f.append(f"{mode} magnetometer falsely adds direct a_w information")
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
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "qualification": d["qualification"],
        "process_inverse_cells": d["validated_process_inverse_cells"],
        "scaled_Qinv_aw_diag_upper": d["scaled_process_inverse_aw_diag_upper"],
        "Racc_variance_lower": d["configured_accelerometer_variance_lower"],
        "H18_aw_precision_trace_upper": d["modes"]["H"]["full_posterior_aw_precision_trace_upper"],
        "A21_aw_precision_trace_upper": d["modes"]["A"]["full_posterior_aw_precision_trace_upper"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
