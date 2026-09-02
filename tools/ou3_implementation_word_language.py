#!/usr/bin/env python3
"""Current OU-III normal-Live source-word language.

This producer binds the declared proof domain directly to the retained source,
translation, and vector certificates.  The former generic source-word wrapper
and separate spread-search producer are intentionally gone; their small proof-
design calculation is performed here so the current route has one language
producer and no historical indirection.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_source_domain_contract as SOURCE
import ou3_translational_uco_ucc as TRANS
import ou3_validated_transcendentals as VT
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 3


def finite_positive(value) -> bool:
    try:
        x = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(x) and x > 0.0


def _point(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def _exp_negative_wide(x: float) -> float:
    if not math.isfinite(x) or x < 0.0:
        raise ValueError("finite nonnegative exponential argument required")
    scale = 1
    while x / scale > VT.MAX_ABS_ARGUMENT:
        scale *= 2
    z = Interval.outward_bounds(x / scale, x / scale)
    y = VT.exp_interval(-z)
    s = scale
    while s > 1:
        y = y.square()
        s //= 2
    return y.lo


def _spread_candidate(q: int, delta_min: float, delta_max: float,
                      tau_min: float, rs_upper: float) -> dict:
    selected_window = math.nextafter(3.0 * q * delta_max, math.inf)
    spacing = math.nextafter(q * delta_min, -math.inf)
    decay = _exp_negative_wide(selected_window / tau_min)
    det = math.nextafter((spacing ** 6 / 12.0) * decay, -math.inf)
    T = Interval(0.0, math.nextafter(selected_window, math.inf))
    t2 = T.square()
    t4 = t2.square()
    t6 = t4 * t2
    row_norm2 = _point(1.0) + t2 + t4 / _point(4.0) + t6 / _point(36.0)
    frob = math.nextafter(2.0 * math.sqrt(row_norm2.hi), math.inf)
    sigma_min = math.nextafter(det / (frob ** 3), -math.inf)
    info = math.nextafter((sigma_min ** 2) / (rs_upper ** 2), -math.inf)
    return {
        "q": int(q),
        "selected_window_s_upper": selected_window,
        "selected_spacing_s_lower": spacing,
        "observation_det_lower": det,
        "observation_frobenius_upper": frob,
        "observation_sigma_min_lower": sigma_min,
        "information_gramian_lambda_min_lower": info,
    }


def _select_spread(word_horizon_s: float, trans: dict) -> dict:
    pseudo = trans["S_observation_uco"]
    delta_min = float(pseudo["pseudo_gap_min_s"])
    delta_max = float(pseudo["pseudo_gap_max_s"])
    tau_min = float(trans["process_ucc"]["tau_s"][0])
    rs_upper = float(pseudo["R_S_filter_std_upper"])
    qmax = int(math.floor(float(word_horizon_s) / (3.0 * delta_max)))
    if qmax < 1:
        raise RuntimeError("word horizon does not contain four source-guaranteed S firings")
    candidates = [
        _spread_candidate(q, delta_min, delta_max, tau_min, rs_upper)
        for q in range(1, qmax + 1)
    ]
    usable = [c for c in candidates if c["information_gramian_lambda_min_lower"] > 0.0]
    if not usable:
        raise RuntimeError("no spread candidate retained a positive information bound")
    best = max(usable, key=lambda c: c["information_gramian_lambda_min_lower"])
    adjacent = candidates[0]
    ratio = (
        best["information_gramian_lambda_min_lower"] /
        adjacent["information_gramian_lambda_min_lower"]
    )
    return {
        "admissible_q_max": qmax,
        "best": best,
        "adjacent_q1": adjacent,
        "information_widening_factor_vs_adjacent_lower": ratio,
    }


def _word_contract(domain: dict, source: dict, trans: dict, vector: dict) -> dict:
    live = domain["normal_live"]
    recurrence = float(live["vector_pe_recurrence_window_s"])
    runtime = source["configured_runtime_assumption"]
    dt = float(runtime["imu_dt_s"])
    packet_gap = list(vector["operating_envelope"]["packet_gap_s"])
    packet_span_upper = float(packet_gap[1])
    pseudo = trans["S_observation_uco"]
    detect = trans["integrator_detectability"]

    failures: list[str] = []
    if source.get("source_complete_parameter_domain") is not True:
        failures.append("source parameter domain is not complete")
    failures.extend(f"translation: {x}" for x in TRANS.validate(trans))
    failures.extend(f"vector: {x}" for x in VECTOR.validate(vector))
    if not finite_positive(recurrence):
        failures.append("PE recurrence window is not finite positive")
    elif recurrence < packet_span_upper:
        failures.append("PE recurrence window is shorter than one vector-packet span")

    ready = not failures
    word_horizon = None
    word_samples_upper = None
    spread = None
    if ready:
        word_horizon = max(recurrence, float(pseudo["aligned_window_s"]))
        word_samples_upper = int(math.ceil(word_horizon / dt)) + 1
        try:
            spread = _select_spread(word_horizon, trans)
        except Exception as exc:
            failures.append(f"spread four-S UCO: {type(exc).__name__}: {exc}")
            ready = False

    pe = dict(vector["operating_envelope"])
    pe.update({
        "recurrence_window_s": recurrence,
        "recurrence_quantifier": (
            "every normal-Live interval of this duration contains at least one certified two-packet vector-PE event"
        ),
        "accelerometer_required_at_both_vector_times": True,
        "two_consecutive_accepted_magnetic_packets_required": True,
        "arbitrary_rejections_between_required_pe_events_allowed": True,
        "hypothesis_origin": "DEPLOYMENT_THEOREM_ASSUMPTION_NOT_TRAJECTORY_FIT",
    })
    best = (spread or {}).get("best", {})

    return {
        "schema": 1,
        "claim": "OU3_CONDITIONAL_SOURCE_COMPLETE_NORMAL_LIVE_WORD_LANGUAGE",
        "source_generated_not_trajectory_fit": True,
        "configured_runtime": runtime,
        "fixed_dimension_modes": {"H": 18, "A": 21},
        "normal_live_scope": {
            "same_mode_only": True,
            "hard_attitude_rewrite_inside_word": False,
            "hybrid_transitions_separate": list(source["hybrid_obligations"]),
            "dimension_change_multiplied_as_square_word": False,
        },
        "source_branch_language": {
            "accelerometer_gate": ["accepted", "rejected"],
            "magnetometer_gate": ["not_due", "accepted", "rejected"],
            "S_zero_pseudo": ["not_due", "due"],
            "aw_covariance_sync": ["not_due", "due_psd_increment"],
            "continuous_parameters": source["validated_parameter_box"]["continuous_parameters"],
            "continuous_parameters_outward_rounded": True,
            "joint_source_reachability_required": True,
            "cartesian_extrema_products_not_a_valid_word": True,
        },
        "translation_recurrence": {
            "full_observability_route": "FOUR_S_SPREAD_COMPLETE_V_P_S_AW_UCO",
            "primary_route": "FOUR_S_SPREAD_COMPLETE_V_P_S_AW_UCO",
            "primary_state_order": ["v", "p", "S", "a_w"],
            "aligned_firing_count": 4,
            "pseudo_gap_min_s": pseudo["pseudo_gap_min_s"],
            "pseudo_gap_max_s": pseudo["pseudo_gap_max_s"],
            "minimum_four_firing_window_s": pseudo["aligned_window_s"],
            "spread_selection": "VALIDATED_MAX_INFORMATION_OVER_ALL_ADMISSIBLE_INTEGER_Q",
            "spread_admissible_q_max": (spread or {}).get("admissible_q_max"),
            "spread_index_q_W": best.get("q"),
            "spread_selected_window_s_upper": best.get("selected_window_s_upper"),
            "spread_selected_spacing_lower_s": best.get("selected_spacing_s_lower"),
            "spread_observation_det_lower": best.get("observation_det_lower"),
            "spread_information_gramian_lambda_min_lower": best.get("information_gramian_lambda_min_lower"),
            "information_widening_factor_vs_adjacent_lower": (
                (spread or {}).get("information_widening_factor_vs_adjacent_lower")
            ),
            "three_firing_integrator_detectability_role": "Riccati_covariance_upper_sharpening_only",
            "three_firing_integrator_detectability_is_promotion_fallback": False,
            "three_firing_detectability_window_s": detect["aligned_window_s"],
            "stable_aw_alpha_upper": detect["stable_aw_alpha_upper"],
            "source_complete": bool(trans["translation_source_complete"] and ready),
        },
        "vector_persistent_excitation": pe,
        "conditional_word_language": {
            "ready": bool(ready),
            "word_horizon_lower_s": word_horizon,
            "word_samples_upper_at_configured_dt": word_samples_upper,
            "one_sample_decrease_required": False,
            "word_endpoint_decrease_required": True,
        },
        "source_complete_relative_to_theorem_hypotheses": bool(ready),
        "failures": failures,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("proof domain must not be trajectory fitted")
    live = domain["normal_live"]
    source = SOURCE.build(SOURCE.DEFAULT_HEADER.resolve())
    trans = TRANS.build()
    vector = VECTOR.build()
    pe = vector["operating_envelope"]

    failures: list[str] = []
    relations = {
        "specific_force_norm_lower_mps2": "ge",
        "magnetic_vector_norm_lower_uT": "ge",
        "vector_sine_separation_lower": "ge",
        "body_rate_norm_upper_deg_s": "le",
    }
    comparison = {}
    for key, relation in relations.items():
        actual = live.get(key)
        generic = pe.get(key)
        ok = finite_positive(actual) and finite_positive(generic)
        if ok:
            a = float(actual)
            g = float(generic)
            ok = a >= g if relation == "ge" else a <= g
        comparison[key] = {
            "declared": actual,
            "generic_contract": generic,
            "required_relation": relation,
            "pass": bool(ok),
        }
        if not ok:
            op = ">=" if relation == "ge" else "<="
            failures.append(
                f"proof domain {key}={actual!r} must be {op} generic vector-UCO bound {generic!r}"
            )

    word = _word_contract(domain, source, trans, vector)
    failures.extend(f"word contract: {x}" for x in word.get("failures", []))
    if word.get("conditional_word_language", {}).get("ready") is not True:
        failures.append("conditional word language is not ready")

    return {
        "schema": SCHEMA,
        "qualification": "DECLARED_DOMAIN_SOURCE_COMPLETE_OU3_NORMAL_LIVE_WORD_LANGUAGE",
        "trajectory_fit": False,
        "operating_domain": live,
        "vector_uco_qualification": vector["qualification"],
        "declared_PE_is_at_least_as_strong_as_generic_contract": not any(
            not row["pass"] for row in comparison.values()
        ),
        "PE_monotone_comparison": comparison,
        "word_contract": word,
        "H_dimension": 18,
        "A_dimension": 21,
        "source_complete_relative_to_declared_theorem_hypotheses": not failures,
        "continuous_word_enclosed": False,
        "nonlinear_word_enclosed": False,
        "theorem_promotion": "NOT_ESTABLISHED",
        "failures": failures,
        "pass": not failures,
    }


def validate(d: dict) -> list[str]:
    failures: list[str] = []
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("trajectory_fit") is not False:
        failures.append("word language is trajectory fitted")
    if d.get("declared_PE_is_at_least_as_strong_as_generic_contract") is not True:
        failures.append("declared PE envelope weakens generic vector-UCO contract")
    if d.get("source_complete_relative_to_declared_theorem_hypotheses") is not True:
        failures.append("declared-domain word language is not source complete")
    word = d.get("word_contract", {})
    if word.get("source_complete_relative_to_theorem_hypotheses") is not True:
        failures.append("word contract is not source complete")
    cw = word.get("conditional_word_language", {})
    if cw.get("ready") is not True or not finite_positive(cw.get("word_samples_upper_at_configured_dt")):
        failures.append("conditional word language is not ready")
    tr = word.get("translation_recurrence", {})
    if tr.get("full_observability_route") != "FOUR_S_SPREAD_COMPLETE_V_P_S_AW_UCO":
        failures.append("translation route is not four-S complete-chain UCO")
    if tr.get("aligned_firing_count") != 4:
        failures.append("translation route does not use four S firings")
    if tr.get("three_firing_integrator_detectability_role") != "Riccati_covariance_upper_sharpening_only":
        failures.append("three-S detectability role changed")
    if tr.get("three_firing_integrator_detectability_is_promotion_fallback") is not False:
        failures.append("three-S detectability became a promotion fallback")
    if cw.get("one_sample_decrease_required") is not False:
        failures.append("one-sample contraction requirement reintroduced")
    if word.get("source_branch_language", {}).get("joint_source_reachability_required") is not True:
        failures.append("joint source reachability is not required")
    if d.get("pass") is not True:
        failures.append("word-language producer failed")
    if d.get("continuous_word_enclosed") is not False or d.get("nonlinear_word_enclosed") is not False:
        failures.append("language stage must not masquerade as enclosure")
    if d.get("theorem_promotion") != "NOT_ESTABLISHED":
        failures.append("language stage must not promote theorem")
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
        "pass": d["pass"],
        "recurrence_window_s": d["operating_domain"]["vector_pe_recurrence_window_s"],
        "word_horizon_lower_s": d["word_contract"]["conditional_word_language"]["word_horizon_lower_s"],
        "word_samples_upper": d["word_contract"]["conditional_word_language"]["word_samples_upper_at_configured_dt"],
        "spread_q": d["word_contract"]["translation_recurrence"]["spread_index_q_W"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
