#!/usr/bin/env python3
"""Fourth P4 widening layer: reuse the exact correction structure proved in P5.

This producer is deliberately downstream of the third-generation P4 stack.  It
changes no filter operation or theorem hypothesis.  It only replaces the coarse
finite-angle vector-residual constants by the exact source identities already
proved by the P5 algebra:

* S=0 has eta_S = 0 exactly;
* for isotropic magnetometer R, K_m v = 0 exactly, so only the tangent effective
  coordinate d_m enters the state correction;
* for the accelerometer J_aw=R_wb is orthogonal/full-rank, so eta_a is exactly
  an effective a_w tangent input rather than an independent measurement-space
  disturbance;
* Joseph information transport and the quaternion reset congruence remain exact,
  so no reset condition-number multiplier is introduced.

On ||c||<=q the exact effective-input bounds give quadratic constants

  C_acc(q) <= (f_max+2)/sqrt(4+q^2),
  C_mag(q) <= m_max/sqrt(4+q^2),

before multiplication by the class-local Kalman gain.  Every proof-radius
candidate is rebuilt with those constants and the same exact endpoint/prefix
budgets as the third-generation certificate.  The result fails closed unless it
is monotone relative to third-generation P4 in both H and A.
"""
from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path

import ou3_p4_nonlinear_word_certificate as LEGACY
import ou3_p4_nextgen_directional_certificate as P4D
import ou3_p4_thirdgen_combined_certificate as THIRD
import ou3_p5_effective_vector_input as VEFF
import ou3_p5_exact_correction_transport as CORR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"


def _sqrt4q_lower(q: float) -> float:
    x = LEGACY.down(4.0 + LEGACY.down(float(q) * float(q)))
    return LEGACY.GROUP.sqrt_point(x).lo


def _candidate(mode: str, base: dict, domain: dict, q: float) -> dict:
    mmin = float(base["metric_lambda_min_lower"])
    mmax = float(base["metric_lambda_max_upper"])
    n = int(base["word_samples_upper"])
    delta = float(base["P3_word_endpoint_delta_lower"])
    gap = P4D._endpoint_sqrt_gap_lower(delta)
    K = base["measurement_specific_gain_norm_upper"]
    H = base["directional_measurement_operator_norm_upper"]
    live = domain["normal_live"]
    fmax = float(live["specific_force_norm_upper_mps2"])
    magmax = float(live["magnetic_vector_norm_upper_uT"])

    den = _sqrt4q_lower(q)
    # Exact P5 effective-input reductions, rewritten as quadratic constants in
    # the full canonical norm because |c|,|delta_aw| <= ||z||.
    Ceta_acc = LEGACY.div_up(LEGACY.add_up(fmax, 2.0), den)
    Ceff_mag = LEGACY.div_up(magmax, den)
    Cinput_acc = LEGACY.mul_up(float(K["accelerometer"]), Ceta_acc)
    # K_m H_theta(d_eff-c_perp), with ||H_theta||=||v|| for -[v]x.
    Cinput_mag = LEGACY.mul_up(float(K["magnetometer"]), Ceff_mag)

    L = {
        "S_zero": LEGACY.mul_up(float(K["S_zero"]), float(H["S_zero"])),
        "accelerometer": LEGACY.mul_up(float(K["accelerometer"]), float(H["accelerometer"])),
        "magnetometer": LEGACY.mul_up(float(K["magnetometer"]), float(H["magnetometer"])),
    }
    cs = LEGACY._composition_quadratic_constant(L["S_zero"], 0.0, q)
    ca = LEGACY._composition_quadratic_constant(L["accelerometer"], Cinput_acc, q)
    cm = LEGACY._composition_quadratic_constant(L["magnetometer"], Cinput_mag, q)
    pred0 = base["prediction_quadratic_bound"]
    pred = LEGACY._composition_quadratic_constant(float(pred0["linear_correction_gain_L"]), 0.0, q)
    C = {
        "prediction": float(pred["full_state_quadratic_defect_constant_upper"]),
        "S_zero_accepted": float(cs["full_state_quadratic_defect_constant_upper"]),
        "accelerometer_accepted": float(ca["full_state_quadratic_defect_constant_upper"]),
        "magnetometer_accepted": float(cm["full_state_quadratic_defect_constant_upper"]),
    }
    csum = LEGACY.add_up(
        LEGACY.add_up(C["prediction"], C["S_zero_accepted"]),
        LEGACY.add_up(C["accelerometer_accepted"], C["magnetometer_accepted"]),
    )
    B = LEGACY.mul_up(LEGACY.PREFIX_BOOTSTRAP_W_FACTOR, float(n))
    B = LEGACY.mul_up(B, LEGACY.sqrt_up(mmax))
    B = LEGACY.mul_up(B, csum)
    B = LEGACY.div_up(B, mmin)
    if not (math.isfinite(B) and B > 0.0):
        raise RuntimeError("invalid exact-structure B")

    sm = LEGACY.GROUP.sqrt_point(mmin).lo
    caps = {
        "endpoint": LEGACY.div_down(gap, B),
        "bootstrap": LEGACY.div_down(1.0, B),
        "design_radius": THIRD._positive_root(B, LEGACY.mul_down(q, sm)),
        "cayley_chart": THIRD._positive_root(B, LEGACY.mul_down(float(base["cayley_norm_limit"]), sm)),
    }
    projection = copy.deepcopy(base.get("active_bias_projection"))
    if mode == "A" and projection:
        caps["bias_projection"] = THIRD._positive_root(
            B, LEGACY.mul_down(float(projection["interior_margin_lower_mps2"]), sm)
        )
    sqrtW = min(caps.values())
    W = LEGACY.mul_down(sqrtW, sqrtW)
    pf = LEGACY.add_up(1.0, LEGACY.mul_up(B, sqrtW))
    qprefix = LEGACY.div_up(LEGACY.mul_up(pf, sqrtW), sm)
    if not (qprefix <= q and qprefix < float(base["cayley_norm_limit"])):
        raise RuntimeError("exact-structure prefix outside design/chart")

    corrections = {
        "S_zero": LEGACY.mul_up(L["S_zero"], qprefix),
        "accelerometer": LEGACY.add_up(
            LEGACY.mul_up(L["accelerometer"], qprefix),
            LEGACY.mul_up(Cinput_acc, qprefix * qprefix),
        ),
        "magnetometer": LEGACY.add_up(
            LEGACY.mul_up(L["magnetometer"], qprefix),
            LEGACY.mul_up(Cinput_mag, qprefix * qprefix),
        ),
    }
    corr_prefix = max(corrections.values())
    if not corr_prefix < 1.0e-2:
        raise RuntimeError("exact-structure correction leaves quaternion series branch")
    if LEGACY.mul_up(B, sqrtW) > gap or LEGACY.mul_up(B, sqrtW) > 1.0:
        raise RuntimeError("exact endpoint/prefix budget failed")
    if mode == "A" and projection:
        margin = float(projection["interior_margin_lower_mps2"])
        projection["certified_error_norm_prefix_upper"] = qprefix
        projection["projection_surface_reached_in_certified_funnel"] = not (qprefix < margin)
        if not qprefix < margin:
            raise RuntimeError("A projection surface reached")

    return {
        "q": q, "B": B, "W": W, "sqrtW": sqrtW, "qprefix": qprefix,
        "prefix_factor": pf, "caps": caps, "active_cap": min(caps, key=caps.get),
        "C": C, "sum": csum, "corr_prefix": corr_prefix,
        "corrections": corrections, "projection": projection,
        "Ceta_acc": Ceta_acc, "Ceff_mag": Ceff_mag,
    }


def _refine_mode(mode: str, base: dict, domain: dict) -> dict:
    q0 = float(base["thirdgen_selected_design_norm"])
    rows = []
    for q in THIRD._grid(q0):
        try:
            rows.append(_candidate(mode, base, domain, q))
        except Exception:
            continue
    if not rows:
        raise RuntimeError(f"{mode}: no exact-correction-structure candidate")
    best = max(rows, key=lambda r: r["W"])
    before = float(base["certified_level_W"])
    if best["W"] < before:
        raise RuntimeError(f"{mode}: exact correction structure regressed third-generation P4")
    m = copy.deepcopy(base)
    m.update({
        "p4_exact_correction_structure_backport": True,
        "p5_S_zero_eta_exact_zero_bound": True,
        "p5_magnetometer_radial_gain_action_exact_zero_bound": True,
        "p5_accelerometer_eta_effective_aw_input_bound": True,
        "p5_exact_joseph_information_identity_bound": True,
        "p5_exact_reset_congruence_bound": True,
        "reset_condition_number_multiplier_used": False,
        "exact_structure_selected_design_norm": best["q"],
        "exact_structure_selected_active_cap": best["active_cap"],
        "exact_structure_acc_effective_quadratic_constant_upper": best["Ceta_acc"],
        "exact_structure_mag_effective_quadratic_constant_upper": best["Ceff_mag"],
        "exact_structure_operation_defect_constants_upper": best["C"],
        "exact_structure_operation_defect_sum_upper": best["sum"],
        "transported_word_defect_B_upper_before_exact_structure": float(base["transported_word_defect_B_upper"]),
        "transported_word_defect_B_upper": best["B"],
        "certified_level_W_before_exact_structure": before,
        "certified_level_W": best["W"],
        "certified_level_sqrt_W": best["sqrtW"],
        "exact_structure_W_widening_factor_vs_thirdgen_lower": LEGACY.div_down(best["W"], before),
        "exact_structure_total_W_widening_factor_vs_legacy_lower": LEGACY.div_down(
            best["W"], float(base["certified_level_W_legacy"])
        ),
        "prefix_W_factor_upper": LEGACY.mul_up(best["prefix_factor"], best["prefix_factor"]),
        "prefix_canonical_error_norm_upper": best["qprefix"],
        "accepted_correction_norm_prefix_upper": best["corr_prefix"],
        "accepted_correction_norms_by_class_upper": best["corrections"],
        "active_bias_projection": best["projection"],
        "exact_structure_candidates": [
            {"q": r["q"], "W": r["W"], "B": r["B"], "qprefix": r["qprefix"], "active_cap": r["active_cap"]}
            for r in rows
        ],
        "exact_nonlinear_word_pass": True,
    })
    return m


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    p = Path(domain_path).resolve()
    domain = json.loads(p.read_text(encoding="utf-8"))
    prev = THIRD.build(p)
    failures = [f"thirdgen: {x}" for x in THIRD.validate(prev)]
    veff = VEFF.build(p)
    corr = CORR.build(p)
    failures += [f"P5 effective-vector: {x}" for x in VEFF.validate(veff)]
    failures += [f"P5 exact-correction: {x}" for x in CORR.validate(corr)]
    if veff.get("magnetometer", {}).get("kalman_gain_radial_action_exact_zero") is not True:
        failures.append("P5 magnetometer radial annihilation unavailable")
    if veff.get("accelerometer", {}).get("standalone_eta_information_penalty_required_for_state_correction") is not False:
        failures.append("P5 accelerometer effective-input reduction unavailable")
    if corr.get("condition_number_multiplier_used_for_reset_transport") is not False:
        failures.append("P5 reset congruence still uses condition-number multiplier")
    modes = {}
    if not failures:
        for mode in ("H", "A"):
            try:
                modes[mode] = _refine_mode(mode, prev["modes"][mode], domain)
            except Exception as exc:
                failures.append(f"{mode}: {exc}")
    out = copy.deepcopy(prev)
    out["qualification"] = "VALIDATED_P4_EXACT_CORRECTION_STRUCTURE_BACKPORT"
    out["claim"] = "P4_EXACT_VECTOR_CORRECTION_STRUCTURE_WIDENED_H_A_WORD_DISSIPATION"
    out["modes"] = modes
    out["p5_effective_vector_input_certificate"] = veff.get("P5_EFFECTIVE_VECTOR_INPUT_CERTIFICATE")
    out["p5_exact_correction_transport_certificate"] = corr.get("P5_EXACT_CORRECTION_TRANSPORT_ALGEBRA_CERTIFICATE")
    out["p4_exact_correction_structure_source_only"] = True
    passed = not failures and len(modes) == 2
    out["P4_EXACT_CORRECTION_STRUCTURE_WORD_CERTIFICATE"] = "PASS" if passed else "FAIL"
    out["P4_EXACT_NONLINEAR_WORD_CERTIFICATE"] = "PASS" if passed else "FAIL"
    out["failures"] = failures
    return out


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("p4_exact_correction_structure_source_only") is not True:
        f.append("exact correction structure backport is not source-only")
    for mode in ("H", "A"):
        m = d.get("modes", {}).get(mode, {})
        if m.get("p4_exact_correction_structure_backport") is not True:
            f.append(f"{mode}: exact correction structure not bound")
            continue
        if float(m.get("certified_level_W", 0.0)) < float(m.get("certified_level_W_before_exact_structure", math.inf)):
            f.append(f"{mode}: exact structure regressed thirdgen")
        if m.get("reset_condition_number_multiplier_used") is not False:
            f.append(f"{mode}: reset condition-number penalty reintroduced")
        if not float(m.get("accepted_correction_norm_prefix_upper", math.inf)) < 1.0e-2:
            f.append(f"{mode}: quaternion branch safety failed")
        if mode == "A" and m.get("active_bias_projection", {}).get("projection_surface_reached_in_certified_funnel") is not False:
            f.append("A: bias projection safety failed")
    if not f and d.get("P4_EXACT_CORRECTION_STRUCTURE_WORD_CERTIFICATE") != "PASS":
        f.append("status not PASS")
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
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P4_EXACT_CORRECTION_STRUCTURE_WORD_CERTIFICATE"],
        "numerical": {m: {
            "W_before": d.get("modes", {}).get(m, {}).get("certified_level_W_before_exact_structure"),
            "W_after": d.get("modes", {}).get(m, {}).get("certified_level_W"),
            "factor": d.get("modes", {}).get(m, {}).get("exact_structure_W_widening_factor_vs_thirdgen_lower"),
            "factor_vs_legacy": d.get("modes", {}).get(m, {}).get("exact_structure_total_W_widening_factor_vs_legacy_lower"),
            "q": d.get("modes", {}).get(m, {}).get("exact_structure_selected_design_norm"),
            "active_cap": d.get("modes", {}).get(m, {}).get("exact_structure_selected_active_cap"),
        } for m in ("H", "A")},
        "failures": f,
    }, indent=2, sort_keys=True))
    return 0 if not f else 2


if __name__ == "__main__":
    raise SystemExit(main())
