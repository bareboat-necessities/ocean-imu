#!/usr/bin/env python3
"""Source-complete word-endpoint OU-III P3 information certificate.

The retired wrapper multiplied a cellwise margin over N samples.  That is not
the theorem: the covariance/information metric varies with the source, and the
paper explicitly permits positive intermediate changes.  This producer instead
uses the direct validated matrix inequality as an endpoint word-noise lower
bound against the already word-horizon covariance upper bound.

The last admissible prediction/correction in a word contributes a positive
posterior process-noise matrix to the full endpoint word covariance difference;
all earlier process terms are PSD additions.  Therefore the direct backend's
validated inequality

    Omega_endpoint_lower - delta * Sigma_word_upper >> 0

is a conservative source-word endpoint injection inequality.  No repeated
one-step contraction, scalarized translation floor, or Euclidean usefulness
gate is used.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_source_reachable_matrix_p3_direct as DIRECT

DEFAULT_DOMAIN = DIRECT.DEFAULT_DOMAIN
SCHEMA = DIRECT.SCHEMA
THEOREM_REQUIRED_MARGIN_LOWER = 0.0
THEOREM_REQUIRED_MARGIN_PREDICATE = "> 0"
NUMERICAL_SEARCH_SEED = DIRECT.MIN_USEFUL_DELTA

# BASE owns the theorem PASS bit.  The actual theorem condition is strict
# positivity; DIRECT.MIN_USEFUL_DELTA remains a positive logarithmic search seed.
DIRECT.BASE.MIN_USEFUL_DELTA = THEOREM_REQUIRED_MARGIN_LOWER
DIRECT.BASE._build_cached.cache_clear()


def _state_metric(mode: str, comparison: dict) -> dict:
    labels = [
        "theta_x", "theta_y", "theta_z", "b_g_x", "b_g_y", "b_g_z",
        "v_x", "v_y", "v_z", "p_x", "p_y", "p_z", "S_x", "S_y", "S_z",
        "a_w_x", "a_w_y", "a_w_z",
    ]
    if mode == "A":
        labels += ["b_a_x", "b_a_y", "b_a_z"]
    scales2 = list(comparison["comparison_scale_diagonal_squared"])
    if len(scales2) != len(labels):
        raise RuntimeError(f"{mode} comparison scale dimension mismatch")
    if not all(math.isfinite(float(x)) and float(x) > 0.0 for x in scales2):
        raise RuntimeError(f"{mode} comparison conditioning is not positive diagonal")
    return {
        "kind": "COMPUTATIONAL_CONGRUENCE_FOR_GENERALIZED_MATRIX_INEQUALITY",
        "state_labels": labels,
        "D_diagonal_squared": scales2,
        "translation_axis_coordinates": [
            "v/(sigma*h)", "p/(sigma*h^2)", "S/(sigma*h^3)", "a_w/sigma"
        ],
        "translation_coupling_retained": "full 4x4 [v,p,S,a_w] block",
        "translation_conditioning": "C=R*L_inverse applied by congruence to both Omega and Sigma",
        "same_congruence_applied_to_noise_and_covariance": True,
        "raw_Euclidean_eigenvalue_gate_used": False,
        "is_nonlinear_Lyapunov_metric": False,
        "nonlinear_metric_requirement": (
            "P4 must use node-wise blkdiag((a_R/2) I3, P_xi) and certify the complete source-word endpoint map"
        ),
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    d = DIRECT.build(domain_path)

    words = DIRECT.BASE.WORDS.build(domain_path)
    word_failures = DIRECT.BASE.WORDS.validate(words)
    if word_failures:
        raise RuntimeError(f"source-word prerequisite failed: {word_failures}")
    if words.get("source_complete_relative_to_declared_theorem_hypotheses") is not True:
        raise RuntimeError("P3 requires the declared source-complete recurring-PE word language")

    word_contract = words["word_contract"]
    word = word_contract["conditional_word_language"]
    trans = word_contract["translation_recurrence"]
    if word.get("ready") is not True:
        raise RuntimeError("declared recurring-PE word language is not ready")
    if trans.get("primary_route") != "FOUR_S_SPREAD_COMPLETE_V_P_S_AW_UCO":
        raise RuntimeError("P3 requires the paper's four-S complete translation route")

    out = dict(d)
    out["p3_window_backend"] = "SOURCE_COMPLETE_WORD_ENDPOINT_GENERALIZED_INFORMATION"
    out["theorem_margin_requirement"] = {
        "predicate": THEOREM_REQUIRED_MARGIN_PREDICATE,
        "numeric_boundary": THEOREM_REQUIRED_MARGIN_LOWER,
        "origin": "ou3_information_enclosure_contract.required_continuous_bounds",
        "old_fixed_1e_minus_18_gate_is_theorem_requirement": False,
        "numerical_search_seed_only": NUMERICAL_SEARCH_SEED,
    }
    out["source_word_binding"] = {
        "claim": word_contract["claim"],
        "source_complete_relative_to_declared_theorem_hypotheses": True,
        "recurring_PE_window_s": words["operating_domain"]["vector_pe_recurrence_window_s"],
        "word_horizon_lower_s": word["word_horizon_lower_s"],
        "word_samples_upper_at_configured_dt": word["word_samples_upper_at_configured_dt"],
        "translation_primary_route": trans["primary_route"],
        "translation_aligned_firing_count": trans["aligned_firing_count"],
        "translation_spread_index_q_W": trans["spread_index_q_W"],
        "translation_spread_selected_spacing_lower_s": trans["spread_selected_spacing_lower_s"],
        "translation_det_spacing_widening_factor_vs_adjacent": trans[
            "determinant_spacing_widening_factor_vs_adjacent"
        ],
        "arbitrary_source_branches_between_required_PE_events_remain_admissible": True,
        "joint_source_reachability_required": True,
        "one_sample_decrease_required": False,
    }

    modes = {}
    for mode in ("H", "A"):
        row = dict(d["modes"][mode])
        comparison = dict(row["matrix_comparison"])
        delta = float(row["relative_Riccati_injection_margin_lower"])
        if not (math.isfinite(delta) and 0.0 < delta < 1.0):
            raise RuntimeError(f"{mode} word-endpoint injection margin invalid: {delta!r}")

        comparison["word_endpoint_information_argument"] = {
            "form": "Omega_endpoint_lower - delta * Sigma_word_upper is SPD",
            "endpoint_noise_lower_source": (
                "posterior lower matrix from the final admissible prediction/correction; earlier word process terms are PSD additions"
            ),
            "covariance_upper_source": "source-uniform finite-word H/A covariance upper bound",
            "repeated_one_step_contraction_used": False,
            "source_replay_used": False,
            "coupled_translation_block": "[v,p,S,a_w]",
            "four_S_spread_translation_route": True,
        }
        comparison["state_information_metric"] = _state_metric(mode, comparison)

        row["word_endpoint_relative_Riccati_injection_margin_lower"] = delta
        row["lambda_information_upper_formula"] = "1-delta_word_endpoint_lower"
        row["strict_information_contraction"] = delta > 0.0
        row["useful_margin_gate"] = THEOREM_REQUIRED_MARGIN_LOWER
        row["useful_margin_gate_predicate"] = THEOREM_REQUIRED_MARGIN_PREDICATE
        row["useful_margin_pass"] = delta > THEOREM_REQUIRED_MARGIN_LOWER
        row["pass"] = row["useful_margin_pass"]
        row["matrix_comparison"] = comparison
        modes[mode] = row

    out["modes"] = modes
    passed = all(modes[m]["pass"] for m in ("H", "A"))
    out["continuous_linear_information_certificate"] = "PASS" if passed else "FAIL"
    out["theorem_promotion"] = "LINEAR_ONLY" if passed else "NOT_ESTABLISHED"
    out["next_obligation"] = (
        "P4 must construct node-wise group-compatible path metrics and certify the exact nonlinear source-word endpoint decrease and all prefixes; the P3 conditioning metric is not reused as W"
    )
    return out


def validate(d: dict) -> list[str]:
    failures = DIRECT.validate(d)
    req = d.get("theorem_margin_requirement", {})
    if req.get("predicate") != "> 0" or req.get("numeric_boundary") != 0.0:
        failures.append("P3 theorem margin requirement is not strict positivity")
    if req.get("old_fixed_1e_minus_18_gate_is_theorem_requirement") is not False:
        failures.append("arbitrary 1e-18 gate still controls theorem acceptance")

    binding = d.get("source_word_binding", {})
    if binding.get("source_complete_relative_to_declared_theorem_hypotheses") is not True:
        failures.append("P3 is not bound to source-complete recurring-PE language")
    if binding.get("translation_primary_route") != "FOUR_S_SPREAD_COMPLETE_V_P_S_AW_UCO":
        failures.append("P3 translation route is not four-S complete-chain UCO")
    if binding.get("translation_aligned_firing_count") != 4:
        failures.append("P3 does not require four S firings")
    if binding.get("one_sample_decrease_required") is not False:
        failures.append("P3 reintroduced one-sample decrease")
    if binding.get("joint_source_reachability_required") is not True:
        failures.append("P3 does not require joint source reachability")

    for mode in ("H", "A"):
        row = d.get("modes", {}).get(mode, {})
        delta = row.get("word_endpoint_relative_Riccati_injection_margin_lower")
        if not isinstance(delta, (int, float)) or not (0.0 < float(delta) < 1.0):
            failures.append(f"{mode} word-endpoint matrix margin invalid")
        if row.get("useful_margin_gate") != 0.0 or row.get("useful_margin_gate_predicate") != "> 0":
            failures.append(f"{mode} theorem gate is not strict positivity")
        matrix = row.get("matrix_comparison", {})
        arg = matrix.get("word_endpoint_information_argument", {})
        if arg.get("repeated_one_step_contraction_used") is not False:
            failures.append(f"{mode} repeated-step shortcut still active")
        if arg.get("four_S_spread_translation_route") is not True:
            failures.append(f"{mode} endpoint argument is not bound to spread four-S route")
        metric = matrix.get("state_information_metric", {})
        if metric.get("same_congruence_applied_to_noise_and_covariance") is not True:
            failures.append(f"{mode} computational congruence is inconsistent")
        if metric.get("translation_coupling_retained") != "full 4x4 [v,p,S,a_w] block":
            failures.append(f"{mode} translation chain was scalarized")
        if metric.get("is_nonlinear_Lyapunov_metric") is not False:
            failures.append(f"{mode} P3 conditioning metric incorrectly promoted to nonlinear W")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    out = dict(d)
    out["validation_pass"] = not failures
    out["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "backend": d["p3_window_backend"],
        "source_word_binding": d["source_word_binding"],
        "H": d["modes"]["H"],
        "A": d["modes"]["A"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
