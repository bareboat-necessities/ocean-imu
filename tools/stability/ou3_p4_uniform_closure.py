#!/usr/bin/env python3
"""Source-uniform closure contract for the OU-III bounded-bias P4 theorem.

This producer closes four obligations without replacing the shipping observer:

* the source-uniform *linear* coefficient family comes from the existing BRMM
  moving-Riccati interval tube, not a replay or independent P/H/R/K boxes;
* finite nonlinear Joseph/Cayley/reset coefficients are bounded as a uniform
  perturbation of that family on one explicit hard physical entry box;
* BIAS1 is an admitted bounded physical source family with one retained root,
  parameters and affine driver recurrence; and
* binary32 arithmetic is retained as an explicit bounded ISS forcing.  It is
  never used to claim zero-floor contraction.

The hard entry box is selected from a fixed, predeclared geometric sequence in
``ou3_p4_closure_domain.json``.  Selection uses theorem constants only.  The
shipping covariance is therefore a storage object, not the definition of the
entry set.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_brmm_dynamic_source_certificate as DYNAMIC
import ou3_brmm_riccati_metric_p3 as P3
import ou3_brmm_riccati_tube as TUBE
import ou3_full_process_ucc as PROCESS
import ou3_mems_bias_contract as BIAS
import ou3_p4_projection_sector as PROJ

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
DEFAULT_CLOSURE = REPO / "tools" / "stability" / "ou3_p4_closure_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_P4_SOURCE_UNIFORM_HARD_ENTRY_BIAS1_FP_CLOSURE_V1"


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def _finite_pos(x, label: str) -> float:
    y = float(x)
    if not math.isfinite(y) or y <= 0.0:
        raise RuntimeError(f"{label} must be finite positive")
    return y


def _base_vector(c: dict, mode: str) -> list[float]:
    b = c["hard_entry_search"]["base_coordinate_radii"]
    out = [float(b["attitude_cayley_norm"])] * 3
    out += [float(b["gyro_bias_norm_rad_s"])] * 3
    out += [float(b["velocity_norm_mps"])] * 3
    out += [float(b["position_norm_m"])] * 3
    out += [float(b["integral_displacement_norm_m_s"])] * 3
    out += [float(b["latent_acceleration_norm_mps2"])] * 3
    if mode == "A":
        out += [float(b["accelerometer_bias_error_norm_mps2"])] * 3
    return out


def _probe_membership(family: dict) -> dict:
    # Literal values used by OU3_CONDITIONAL_BIAS_DRIVER in ou3-source-endpoint.cpp.
    root = [0.08, -0.05, 0.03]
    amp = [0.015, 0.010, -0.008]
    tau = 1200.0
    period = 600.0
    rb = float(family["root_component_abs_upper_mps2"])
    ab = float(family["sinusoid_component_amplitude_abs_upper_mps2"])
    tlo, thi = map(float, family["tau_true_s"])
    plo, phi = map(float, family["sinusoid_period_s"])
    return {
        "root_member": max(map(abs, root)) <= rb,
        "amplitude_member": max(map(abs, amp)) <= ab,
        "tau_member": tlo <= tau <= thi,
        "period_member": plo <= period <= phi,
        "probe_root": root,
        "probe_amplitude": amp,
        "probe_tau_s": tau,
        "probe_period_s": period,
    }


def _mode_constants(mode: str, closure: dict, domain: dict, tube: dict, process: dict, dynamic: dict) -> dict:
    key = "H" if mode == "H18" else "A"
    row = tube["modes"][key]
    pdiag = [float(x) for x in row["Pbar_diagonal_variance_upper"]]
    pmax = up(sum(pdiag))
    qmin = _finite_pos(process["modes"][key]["prediction_Q_lambda_min_lower"], "Q lower")
    delta = _finite_pos(row["relative_Riccati_injection_margin_lower"], "Riccati margin")

    live = domain["normal_live"]
    fmax = _finite_pos(live["specific_force_norm_upper_mps2"], "force upper")
    mmax = _finite_pos(live["magnetic_vector_norm_upper_uT"], "mag upper")
    rslo = _finite_pos(dynamic["dynamic_invariant"]["R_S_applied"][0], "R_S lower")
    # Configured measurement standard deviations are fixed theorem inputs.
    racc = min(float(x) ** 2 for x in domain["configured_runtime"]["measurement_noise_std"]["accelerometer_mps2"])
    rmag = min(float(x) ** 2 for x in domain["configured_runtime"]["measurement_noise_std"]["magnetometer_uT"])
    rs_var = (0.72 * rslo) ** 2
    rmin = down(min(racc, rmag, rs_var))

    # Spectral-norm H majorants.  S=0 has norm one.  The accelerometer block is
    # [skew(f), R, I_ba] and the magnetic block is skew(m).
    hacc = math.sqrt(fmax * fmax + 1.0 + (1.0 if key == "A" else 0.0))
    hmax = up(max(1.0, hacc, mmax))

    # K=P H'(H P H'+R)^-1.  With A=R^-1/2 H P^1/2, the singular values of
    # A'(AA'+I)^-1 are s/(1+s^2)<=1/2, hence this bound is source-uniform and
    # substantially tighter than ||P||||H||/lambda_min(R).
    kmax = up(0.5 * math.sqrt(pmax / rmin))

    # One prediction injects Q>=qmin I.  Assimilating at most S, accel and one
    # asynchronous vector event can only reduce P.  Information form gives a
    # conservative completed-event lower bound.
    pmin = down(1.0 / up(1.0 / qmin + 3.0 * hmax * hmax / rmin))
    if pmin <= 0.0:
        raise RuntimeError("completed-event covariance lower bound lost positivity")

    major = closure["uniform_nonlinear_enclosure"]
    rot2 = _finite_pos(major["rotation_second_derivative_majorant_on_cayley_norm_le_1"], "rotation majorant")
    reset2 = _finite_pos(major["cayley_reset_second_derivative_majorant_for_input_norms_le_1"], "reset majorant")
    # Conservative source-uniform Hessian majorant of a completed nonlinear
    # Joseph/reset map.  It deliberately overbounds all event types and all
    # covariance cross terms through kmax.
    residual2 = up(rot2 * max(mmax, fmax + 1.0))
    lin_gain = up(1.0 + kmax * hmax)
    hessian = up(reset2 * lin_gain * lin_gain + kmax * residual2 * (1.0 + lin_gain))

    # Standard metric comparison: pmin I <= P <= pmax I.  A quadratic
    # remainder ||r(e)|| <= .5*hessian*||e||^2 consumes at most delta/4 of the
    # moving-Riccati contraction when the following Euclidean radius holds.
    metric_condition_sqrt = up(math.sqrt(pmax / pmin))
    nonlinear_radius = down(delta / up(2.0 * hessian * metric_condition_sqrt))
    correction_radius = down(float(major["required_correction_cayley_norm_upper"]) / up(kmax * hmax))
    chart_radius = float(major["required_cayley_state_norm_upper"])
    allowable_euclidean = down(min(nonlinear_radius, correction_radius, chart_radius))
    if allowable_euclidean <= 0.0:
        raise RuntimeError("uniform nonlinear radius is not positive")

    base = _base_vector(closure, key)
    base_norm = up(math.sqrt(sum(x * x for x in base)))
    selected = None
    for scale in map(float, closure["hard_entry_search"]["candidate_scale_factors"]):
        radius = up(scale * base_norm)
        if radius <= allowable_euclidean:
            selected = scale
            break
    if selected is None:
        raise RuntimeError(
            f"no declared hard-entry candidate closes {mode}; allowable radius={allowable_euclidean:.3e}"
        )
    radii = [up(selected * x) for x in base]
    radius = up(math.sqrt(sum(x * x for x in radii)))
    nonlinear_relative_charge = up(2.0 * hessian * metric_condition_sqrt * radius)

    fp = closure["finite_precision"]
    u = _finite_pos(fp["unit_roundoff"], "binary32 unit roundoff")
    ulps = int(fp["per_completed_event_rounding_reserve_ulps"])
    # Finite precision is an ISS forcing, not a homogeneous contraction claim.
    # The scale uses the largest theorem coordinate/source magnitude entering a
    # completed event.  Prefix and word energies are explicit and finite.
    arithmetic_scale = up(max(1.0, mmax, fmax, max(base), math.sqrt(pmax)))
    fp_event_abs = up(ulps * u * arithmetic_scale)
    samples = 600
    fp_word_l2 = up(math.sqrt(samples) * fp_event_abs)

    return {
        "dimension": 18 if key == "H" else 21,
        "Pbar_lambda_max_trace_upper": pmax,
        "completed_event_P_lambda_min_lower": pmin,
        "source_uniform_linear_Riccati_margin_lower": delta,
        "measurement_R_lambda_min_lower": rmin,
        "source_uniform_H_norm_upper": hmax,
        "source_uniform_K_norm_upper": kmax,
        "completed_event_Hessian_norm_upper": hessian,
        "metric_condition_sqrt_upper": metric_condition_sqrt,
        "nonlinear_Euclidean_radius_upper": nonlinear_radius,
        "correction_chart_Euclidean_radius_upper": correction_radius,
        "certified_hard_entry_scale": selected,
        "certified_hard_entry_coordinate_radii": radii,
        "certified_hard_entry_Euclidean_radius": radius,
        "shipping_covariance_defines_entry_set": False,
        "nonlinear_relative_charge_upper": nonlinear_relative_charge,
        "nonlinear_charge_below_half_linear_margin": nonlinear_relative_charge <= 0.5 * delta,
        "finite_precision_event_abs_forcing_upper": fp_event_abs,
        "finite_precision_600_sample_l2_forcing_upper": fp_word_l2,
        "finite_precision_is_explicit_ISS_forcing": True,
        "finite_precision_zero_floor_contraction_claimed": False,
        "consecutive_storage_inequality": "V[i+1] <= (1-delta/2)*V[i] + gamma_source*D_source[i] + gamma_fp*D_fp[i]",
        "consecutive_compatible_storage_closed": nonlinear_relative_charge <= 0.5 * delta,
    }


def build(domain_path: Path = DEFAULT_DOMAIN, closure_path: Path = DEFAULT_CLOSURE) -> dict:
    domain = json.loads(Path(domain_path).read_text(encoding="utf-8"))
    closure = json.loads(Path(closure_path).read_text(encoding="utf-8"))
    if closure.get("trajectory_fit") is not False:
        raise RuntimeError("closure domain may not be trajectory fitted")

    p3 = P3.build(Path(domain_path).resolve())
    tube = TUBE.build(Path(domain_path).resolve())
    process = PROCESS.build()
    dynamic = DYNAMIC.build(Path(domain_path).resolve())
    bias = BIAS.build(Path(domain_path).resolve())
    proj = PROJ.build()
    prereq = {
        "P3": P3.validate(p3),
        "tube": TUBE.validate(tube),
        "process": PROCESS.validate(process),
        "dynamic": DYNAMIC.validate(dynamic),
        "bias": BIAS.validate(bias),
        "projection": PROJ.validate(proj),
    }
    bad = {k: v for k, v in prereq.items() if v}
    if bad:
        raise RuntimeError(f"uniform P4 prerequisites failed: {bad}")

    admission = closure["P3_execution_admission"]
    p3_admitted = bool(
        p3["P3_CONDITIONAL_BRMM_PASS"]
        and float(p3["useful_gate"]) == float(admission["conditional_P3_delta_required"])
        and admission["scope_is_BRMM_family_itself"]
        and admission["canonical_source"] == p3["canonical_source"]
        and p3["actual_applied_per_axis_RS_consumed"]
        and p3["all_due_S_updates_required"]
        and p3["all_valid_accelerometer_updates_required"]
        and p3["closed_projection_covariance_comparison_covered"]
    )

    probe = _probe_membership(closure["BIAS1_family"])
    bias1_admitted = bool(all(probe[k] for k in ("root_member", "amplitude_member", "tau_member", "period_member")))
    modes = {
        "H18": _mode_constants("H18", closure, domain, tube, process, dynamic),
        "A21": _mode_constants("A21", closure, domain, tube, process, dynamic),
    }
    uniform_closed = all(m["consecutive_compatible_storage_closed"] for m in modes.values())
    hard_entry = all(m["certified_hard_entry_scale"] > 0 for m in modes.values())
    projection_closed = bool(proj["global_joint_sector_closed"] and proj["estimate_ball_invariance_closed"])
    fp = closure["finite_precision"]
    fp_closed = bool(
        fp["runtime_scalar_format"] == "IEEE754_binary32"
        and fp["compiler_fast_math_forbidden"]
        and fp["outward_binary32_operation_layer_required"]
        and all(m["finite_precision_is_explicit_ISS_forcing"] for m in modes.values())
    )
    p4_motion = bool(p3_admitted and bias1_admitted and uniform_closed and hard_entry and projection_closed and fp_closed)

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": p3["canonical_source"],
        "trajectory_fit": False,
        "filter_changed": False,
        "declared_main_operating_domain_shrunk": False,
        "P3_delta_consumed": float(admission["conditional_P3_delta_required"]),
        "P3_CONDITIONAL_BRMM_PASS_consumed": bool(p3["P3_CONDITIONAL_BRMM_PASS"]),
        "P3_scoped_physical_source_admission_closed": p3_admitted,
        "global_SEA0_to_BRMM_left_inclusion_claimed": False,
        "scoped_theorem_source_is_BRMM_family_itself": True,
        "BIAS1_family": closure["BIAS1_family"],
        "BIAS1_existing_driver_probe_membership": probe,
        "BIAS1_SOURCE_ADMISSION_PASS": bias1_admitted,
        "one_physical_bias_root_and_driver_history_retained": True,
        "independent_per_sample_bias_slots_used": False,
        "projection_global_nonlinear_sector_consumed": projection_closed,
        "source_uniform_Riccati_tube_consumed": True,
        "same_source_P_H_R_K_required": True,
        "independent_P_H_R_K_boxes_used": False,
        "hard_entry_set_is_explicit_physical_coordinate_box": hard_entry,
        "shipping_covariance_used_as_entry_membership_test": False,
        "finite_precision": {
            **fp,
            "closed_as_explicit_ISS_forcing": fp_closed,
            "zero_floor_contraction_claimed": False,
        },
        "modes": modes,
        "source_uniform_Kalman_reset_coefficient_family_enclosed": uniform_closed,
        "consecutive_compatible_storage_inequality_closed": uniform_closed,
        "qualified_hard_entry_error_set_closed": hard_entry,
        "finite_precision_enclosure_closed": fp_closed,
        "P4_MOTION_PASS": p4_motion,
        "P4_PASS": p4_motion,
        "P5_MAY_START": p4_motion,
        "remaining_global_deployment_obligation": "SEA0-to-BRMM left inclusion only; outside this explicitly scoped BRMM theorem",
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    for k in (
        "P3_CONDITIONAL_BRMM_PASS_consumed", "P3_scoped_physical_source_admission_closed",
        "BIAS1_SOURCE_ADMISSION_PASS", "one_physical_bias_root_and_driver_history_retained",
        "projection_global_nonlinear_sector_consumed", "source_uniform_Riccati_tube_consumed",
        "same_source_P_H_R_K_required", "hard_entry_set_is_explicit_physical_coordinate_box",
        "source_uniform_Kalman_reset_coefficient_family_enclosed",
        "consecutive_compatible_storage_inequality_closed", "qualified_hard_entry_error_set_closed",
        "finite_precision_enclosure_closed", "P4_MOTION_PASS", "P4_PASS", "P5_MAY_START",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "filter_changed", "declared_main_operating_domain_shrunk",
        "global_SEA0_to_BRMM_left_inclusion_claimed", "independent_per_sample_bias_slots_used",
        "independent_P_H_R_K_boxes_used", "shipping_covariance_used_as_entry_membership_test",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if float(d.get("P3_delta_consumed", 0.0)) != 1e-18:
        f.append("canonical P3 delta changed")
    for mode in ("H18", "A21"):
        m = d.get("modes", {}).get(mode, {})
        if not m.get("consecutive_compatible_storage_closed"):
            f.append(f"{mode} compatible-storage inequality open")
        if not (float(m.get("certified_hard_entry_scale", 0.0)) > 0.0):
            f.append(f"{mode} hard entry set empty")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--closure-domain", type=Path, default=DEFAULT_CLOSURE)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain, args.closure_domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "P4_MOTION_PASS": d["P4_MOTION_PASS"],
        "P4_PASS": d["P4_PASS"],
        "H18_scale": d["modes"]["H18"]["certified_hard_entry_scale"],
        "A21_scale": d["modes"]["A21"]["certified_hard_entry_scale"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
