#!/usr/bin/env python3
"""Radius-continuation widening of the OU-III P4 exact nonlinear word certificate.

This producer is independent of the operation-class refinements developed on
other branches.  It starts from the shipping-source P4 theorem on main and
removes two conservative choices that are not part of the theorem itself:

1. The nonlinear Cayley/quaternion defect constant is no longer frozen at one
   arbitrarily tiny design radius.  For a monotone sequence of candidate
   canonical radii q, the exact source-uniform defect constant C(q) is rebuilt
   with outward-rounded arithmetic.  Because C(q) bounds the entire ball
   ||z||<=q, each accepted candidate is itself a rigorous P4 enclosure; no
   trajectory replay or sampled-state inference is used.

2. Prefix safety is enforced with the actual defect-bootstrap inequality

       sqrt(W_prefix) <= sqrt(W0) + B(q) W0
                        = sqrt(W0) [1 + B(q) sqrt(W0)],

   rather than replacing the bracket by the fixed factor 2 before solving for
   the admissible level.  The old bootstrap is still recovered whenever
   B sqrt(W)<=1, so this is a theorem-preserving sharpening.

For the endpoint, the same advertised contraction

       W_end <= (1-delta/2) W0

is enforced through the cancellation-free exact square-root budget

       sqrt(1-delta) + B(q) sqrt(W0) <= sqrt(1-delta/2).

The search is over proof radii, not trajectories.  Every grid point is checked
independently and the producer fails closed unless the selected H and A levels
are at least as large as legacy P4 and all Cayley, quaternion-series, and active
bias-projection guards remain valid.
"""
from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path

import ou3_p4_nonlinear_word_certificate as LEGACY

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def sqrt_down(x: float) -> float:
    return LEGACY.GROUP.sqrt_point(float(x)).lo


def _sqrt_one_minus_up(x: float) -> float:
    if not 0.0 <= x < 1.0:
        raise ValueError("x must be in [0,1)")
    return LEGACY.sqrt_up(math.nextafter(1.0 - float(x), math.inf))


def _endpoint_sqrt_gap_lower(delta: float) -> float:
    """Lower bound on sqrt(1-delta/2)-sqrt(1-delta), without cancellation."""
    if not 0.0 < delta < 1.0:
        raise ValueError("delta must be in (0,1)")
    num = LEGACY.mul_down(0.5, delta)
    den = LEGACY.add_up(_sqrt_one_minus_up(0.5 * delta), _sqrt_one_minus_up(delta))
    return LEGACY.div_down(num, den)


def _positive_quadratic_root_lower(B: float, a_lower: float) -> float:
    """Safe lower bound on positive root of B s^2+s-a=0, B,a>=0."""
    if not (B > 0.0 and a_lower >= 0.0):
        raise ValueError("positive B and nonnegative a required")
    four_Ba = LEGACY.mul_up(4.0, LEGACY.mul_up(B, a_lower))
    disc = LEGACY.add_up(1.0, four_Ba)
    denom = LEGACY.add_up(1.0, LEGACY.sqrt_up(disc))
    return LEGACY.div_down(LEGACY.mul_down(2.0, a_lower), denom)


def _candidate(mode: str, base: dict, domain: dict, q_design: float) -> dict:
    if not q_design > 0.0:
        raise ValueError("q_design must be positive")
    live = domain["normal_live"]
    mmin = float(base["metric_lambda_min_lower"])
    mmax = float(base["metric_lambda_max_upper"])
    delta = float(base["P3_word_endpoint_delta_lower"])
    samples = int(base["word_samples_upper"])
    operations = int(base["state_operation_count_upper"])

    Kmax = float(base["full_gain_norm_upper"])
    Lcorr = float(base["correction_quadratic_bound"]["linear_correction_gain_L"])
    fmax = float(live["specific_force_norm_upper_mps2"])
    magmax = float(live["magnetic_vector_norm_upper_uT"])
    Cvec_acc = LEGACY.add_up(LEGACY.mul_up(LEGACY.ROTATION_REMAINDER_COEFF, fmax), 1.5)
    Cvec_mag = LEGACY.mul_up(LEGACY.ROTATION_REMAINDER_COEFF, magmax)
    Cinput = LEGACY.mul_up(Kmax, max(Cvec_acc, Cvec_mag))
    corr = LEGACY._composition_quadratic_constant(Lcorr, Cinput, q_design)

    dt = float(domain["configured_runtime"]["imu_dt_s"])
    omega = float(live["body_rate_norm_upper_deg_s"]) * math.pi / 180.0
    omega_dt = LEGACY.up(omega * dt)
    if not omega_dt < 0.01:
        raise RuntimeError("configured body rotation per sample exceeds P4 prediction envelope")
    Lpred = LEGACY.div_up(dt, LEGACY.down(1.0 - omega_dt))
    pred = LEGACY._composition_quadratic_constant(Lpred, 0.0, q_design)

    Coperation = max(
        float(corr["full_state_quadratic_defect_constant_upper"]),
        float(pred["full_state_quadratic_defect_constant_upper"]),
    )
    B = LEGACY.mul_up(LEGACY.PREFIX_BOOTSTRAP_W_FACTOR, float(operations))
    B = LEGACY.mul_up(B, LEGACY.sqrt_up(mmax))
    B = LEGACY.mul_up(B, Coperation)
    B = LEGACY.div_up(B, mmin)
    if not (math.isfinite(B) and B > 0.0):
        raise RuntimeError(f"{mode}: invalid radius-dependent B")

    gap = _endpoint_sqrt_gap_lower(delta)
    s_endpoint = LEGACY.div_down(gap, B)
    s_bootstrap = LEGACY.div_down(1.0, B)

    sqrt_mmin_lo = sqrt_down(mmin)
    a_design = LEGACY.mul_down(q_design, sqrt_mmin_lo)
    s_design = _positive_quadratic_root_lower(B, a_design)

    a_cayley = LEGACY.mul_down(float(base["cayley_norm_limit"]), sqrt_mmin_lo)
    s_cayley = _positive_quadratic_root_lower(B, a_cayley)

    caps = {
        "endpoint": s_endpoint,
        "bootstrap": s_bootstrap,
        "design_radius": s_design,
        "cayley_chart": s_cayley,
    }
    projection = copy.deepcopy(base.get("active_bias_projection"))
    if mode == "A" and projection is not None:
        margin = float(projection["interior_margin_lower_mps2"])
        a_projection = LEGACY.mul_down(margin, sqrt_mmin_lo)
        caps["bias_projection"] = _positive_quadratic_root_lower(B, a_projection)

    sqrtW = min(caps.values())
    if not sqrtW > 0.0:
        raise RuntimeError(f"{mode}: radius continuation produced no positive level")
    Wstar = LEGACY.mul_down(sqrtW, sqrtW)

    prefix_factor = LEGACY.add_up(1.0, LEGACY.mul_up(B, sqrtW))
    qprefix = LEGACY.div_up(LEGACY.mul_up(prefix_factor, sqrtW), sqrt_down(mmin))
    if not qprefix <= q_design:
        raise RuntimeError(f"{mode}: radius-dependent design ball does not contain every prefix")
    if not qprefix < float(base["cayley_norm_limit"]):
        raise RuntimeError(f"{mode}: radius continuation reaches Cayley chart limit")

    correction_prefix = LEGACY.add_up(
        LEGACY.mul_up(Lcorr, qprefix),
        LEGACY.mul_up(Cinput, LEGACY.mul_up(qprefix, qprefix)),
    )
    if not correction_prefix < 1.0e-2:
        raise RuntimeError(f"{mode}: accepted correction can leave deployed quaternion series branch")

    if mode == "A" and projection is not None:
        margin = float(projection["interior_margin_lower_mps2"])
        projection["certified_error_norm_prefix_upper"] = qprefix
        projection["projection_surface_reached_in_certified_funnel"] = not (qprefix < margin)
        if not qprefix < margin:
            raise RuntimeError("A: radius-continuation funnel reaches bias projection surface")

    nonlinear_sqrt = LEGACY.mul_up(B, sqrtW)
    if nonlinear_sqrt > gap:
        raise RuntimeError(f"{mode}: exact endpoint square-root budget does not close")
    if nonlinear_sqrt > 1.0:
        raise RuntimeError(f"{mode}: prefix bootstrap B*sqrt(W)<=1 does not close")

    return {
        "q_design": q_design,
        "Coperation": Coperation,
        "B": B,
        "endpoint_sqrt_gap_lower": gap,
        "sqrtW": sqrtW,
        "W": Wstar,
        "prefix_factor_upper": prefix_factor,
        "qprefix": qprefix,
        "correction_prefix": correction_prefix,
        "caps": caps,
        "active_cap": min(caps, key=caps.get),
        "correction_bound": corr,
        "prediction_bound": pred,
        "projection": projection,
        "word_samples_upper": samples,
    }


def _radius_grid(legacy_q: float) -> list[float]:
    # Include the exact legacy design radius plus a dyadic continuation on both
    # sides.  The composition primitive itself fails closed when its deployed
    # quaternion/Cayley source branch ceases to be certified.
    vals = {float(legacy_q)}
    q = max(float(legacy_q) / 64.0, 1.0e-12)
    for _ in range(28):
        vals.add(q)
        q *= 2.0
    return sorted(vals)


def _refine_mode(mode: str, base: dict, domain: dict) -> dict:
    legacy_q = float(base["correction_quadratic_bound"]["design_error_norm_radius"])
    rows = []
    for q in _radius_grid(legacy_q):
        try:
            rows.append(_candidate(mode, base, domain, q))
        except Exception:
            # Larger radii are not assumed valid merely because smaller ones are.
            # We retain all independently certified candidates and select from them.
            continue
    if not rows:
        raise RuntimeError(f"{mode}: no radius-continuation candidate certified")

    best = max(rows, key=lambda r: r["W"])
    legacy_W = float(base["certified_level_W"])
    if best["W"] < legacy_W:
        raise RuntimeError(f"{mode}: radius continuation regressed legacy P4")

    m = copy.deepcopy(base)
    m.update({
        "radius_continuation_used": True,
        "radius_continuation_source_only": True,
        "radius_candidate_count_certified": len(rows),
        "radius_grid_candidate_count": len(_radius_grid(legacy_q)),
        "radius_selected_design_norm": best["q_design"],
        "radius_selected_active_cap": best["active_cap"],
        "radius_endpoint_sqrt_gap_lower": best["endpoint_sqrt_gap_lower"],
        "radius_prefix_factor_upper": best["prefix_factor_upper"],
        "radius_legacy_fixed_prefix_factor": LEGACY.PREFIX_BOOTSTRAP_W_FACTOR ** 0.5,
        "radius_exact_prefix_bootstrap_used": True,
        "radius_exact_endpoint_budget_used": True,
        "radius_continuation_candidates": [
            {
                "q_design": r["q_design"],
                "W": r["W"],
                "sqrtW": r["sqrtW"],
                "B": r["B"],
                "qprefix": r["qprefix"],
                "active_cap": r["active_cap"],
            }
            for r in rows
        ],
        "transported_word_defect_B_upper_legacy": float(base["transported_word_defect_B_upper"]),
        "transported_word_defect_B_upper": best["B"],
        "uniform_operation_quadratic_defect_constant_upper": best["Coperation"],
        "correction_quadratic_bound": best["correction_bound"],
        "prediction_quadratic_bound": best["prediction_bound"],
        "certified_level_W_legacy": legacy_W,
        "certified_level_sqrt_W_legacy": float(base["certified_level_sqrt_W"]),
        "certified_level_W": best["W"],
        "certified_level_sqrt_W": best["sqrtW"],
        "radius_W_widening_factor_lower": LEGACY.div_down(best["W"], legacy_W),
        "radius_sqrtW_widening_factor_lower": LEGACY.div_down(
            best["sqrtW"], float(base["certified_level_sqrt_W"])
        ),
        "prefix_W_factor_upper": LEGACY.mul_up(best["prefix_factor_upper"], best["prefix_factor_upper"]),
        "prefix_canonical_error_norm_upper": best["qprefix"],
        "accepted_correction_norm_prefix_upper": best["correction_prefix"],
        "active_bias_projection": best["projection"],
        "nonlinear_sqrt_budget_fraction_of_delta_upper": LEGACY.div_up(
            LEGACY.mul_up(best["B"], best["sqrtW"]), float(base["P3_word_endpoint_delta_lower"])
        ),
        "exact_nonlinear_word_pass": True,
    })
    return m


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("radius-continuation theorem domain is trajectory fitted")

    legacy = LEGACY.build(domain_path)
    failures = [f"legacy P4: {x}" for x in LEGACY.validate(legacy)]
    modes = {}
    if not failures:
        for mode in ("H", "A"):
            try:
                modes[mode] = _refine_mode(mode, legacy["modes"][mode], domain)
            except Exception as exc:
                failures.append(f"{mode}: {exc}")

    out = copy.deepcopy(legacy)
    out["schema"] = SCHEMA
    out["qualification"] = "VALIDATED_RADIUS_CONTINUATION_CAYLEY_NONLINEAR_SOURCE_WORD_CERTIFICATE"
    out["claim"] = "P4_RADIUS_CONTINUATION_WIDENED_EXACT_NONLINEAR_H_A_WORD_DISSIPATION"
    out["modes"] = modes
    out["radius_continuation_source_only"] = True
    out["radius_continuation_trajectory_sampling_used"] = False
    out["radius_continuation_changes_filter"] = False
    out["radius_continuation_changes_P3_margin"] = False
    out["source_subdivision"] = {
        "kind": "MONOTONE_RADIUS_BALL_CONTINUATION",
        "description": "rebuild outward-rounded exact nonlinear defect constants on nested canonical balls and retain the largest certified P4 level",
        "trajectory_fit": False,
        "independent_of_operation_class_refinements": True,
    }
    passed = not failures and all(modes.get(k, {}).get("exact_nonlinear_word_pass") for k in ("H", "A"))
    out["P4_EXACT_NONLINEAR_WORD_CERTIFICATE"] = "PASS" if passed else "FAIL"
    out["P4_RADIUS_CONTINUATION_WORD_CERTIFICATE"] = "PASS" if passed else "FAIL"
    out["theorem_promotion"] = "P4_RADIUS_CONTINUATION_NORMAL_LIVE_EXACT_WORDS" if passed else "NOT_ESTABLISHED"
    out["failures"] = failures
    out["next_obligation"] = (
        "compare radius-continuation P4 against the independent operation-class/gain P4 branch; combine only theorem-preserving monotone refinements that stack"
    )
    return out


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("radius_continuation_source_only") is not True:
        failures.append("radius continuation is not source-only")
    if d.get("radius_continuation_trajectory_sampling_used") is not False:
        failures.append("radius continuation uses trajectory sampling")
    if d.get("radius_continuation_changes_filter") is not False:
        failures.append("radius continuation changes the filter")
    if d.get("radius_continuation_changes_P3_margin") is not False:
        failures.append("radius continuation changes P3 margin")
    for mode in ("H", "A"):
        m = d.get("modes", {}).get(mode, {})
        if m.get("radius_continuation_used") is not True:
            failures.append(f"{mode}: radius continuation missing")
            continue
        if m.get("radius_exact_prefix_bootstrap_used") is not True:
            failures.append(f"{mode}: exact prefix bootstrap missing")
        if m.get("radius_exact_endpoint_budget_used") is not True:
            failures.append(f"{mode}: exact endpoint budget missing")
        if float(m.get("certified_level_W", 0.0)) < float(m.get("certified_level_W_legacy", math.inf)):
            failures.append(f"{mode}: certified W regressed")
        if float(m.get("radius_W_widening_factor_lower", 0.0)) < 1.0:
            failures.append(f"{mode}: W widening factor below one")
        if float(m.get("radius_sqrtW_widening_factor_lower", 0.0)) < 1.0:
            failures.append(f"{mode}: sqrt(W) widening factor below one")
        if not float(m.get("prefix_canonical_error_norm_upper", math.inf)) < float(m.get("cayley_norm_limit", 0.0)):
            failures.append(f"{mode}: Cayley prefix safety failed")
        if not float(m.get("accepted_correction_norm_prefix_upper", math.inf)) < 1.0e-2:
            failures.append(f"{mode}: quaternion series branch safety failed")
        if mode == "A":
            p = m.get("active_bias_projection", {})
            if p.get("projection_surface_reached_in_certified_funnel") is not False:
                failures.append("A: radius-continuation funnel reaches bias projection")
    if not failures and d.get("P4_RADIUS_CONTINUATION_WORD_CERTIFICATE") != "PASS":
        failures.append("radius-continuation P4 status is not PASS")
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
    compact = {
        mode: {
            "W_legacy": d.get("modes", {}).get(mode, {}).get("certified_level_W_legacy"),
            "W_radius": d.get("modes", {}).get(mode, {}).get("certified_level_W"),
            "W_factor": d.get("modes", {}).get(mode, {}).get("radius_W_widening_factor_lower"),
            "sqrtW_factor": d.get("modes", {}).get(mode, {}).get("radius_sqrtW_widening_factor_lower"),
            "q_selected": d.get("modes", {}).get(mode, {}).get("radius_selected_design_norm"),
            "q_prefix": d.get("modes", {}).get(mode, {}).get("prefix_canonical_error_norm_upper"),
            "active_cap": d.get("modes", {}).get(mode, {}).get("radius_selected_active_cap"),
        }
        for mode in ("H", "A")
    }
    print(json.dumps({
        "P4_RADIUS_CONTINUATION_WORD_CERTIFICATE": d["P4_RADIUS_CONTINUATION_WORD_CERTIFICATE"],
        "numerical": compact,
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
