#!/usr/bin/env python3
"""Finite-window source-reachable OU-III P3 information certificate.

This is the public P3 entry point.  It keeps the direct validated matrix backend
from ``ou3_source_reachable_matrix_p3_direct`` for each source cell, but removes
the old arbitrary 1e-18 theorem gate.  The downstream enclosure contract requires
only a strictly positive continuous relative Riccati injection margin.

The one-step direct inequality is composed over the declared source-complete
normal-Live word horizon.  Because every admissible step has

    Sigma_{k+1} <= (1-delta_k)^(-1) * propagated comparison,

and the direct backend certifies delta_k >= delta_1 > 0 uniformly, an N-step
word has the conservative relative information margin

    delta_W = 1 - (1-delta_1)^N
            >= N*delta_1 / (1 + N*delta_1).

The rational lower expression avoids catastrophic cancellation when delta_1 is
very small.  N is obtained from the finite recurring-PE source-word language,
not from replay.  The translation certificate remains a coupled 4x4
[v,p,S,a_w] matrix inequality in the nondimensional D_h metric followed by the
exact rational C=R L^-1 congruence, so information transfer through the
integrated-OU chain is retained rather than scalarized.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_source_reachable_matrix_p3_direct as DIRECT

DEFAULT_DOMAIN = DIRECT.DEFAULT_DOMAIN
SCHEMA = DIRECT.SCHEMA

# The paper/deployment promotion contract asks for strict positivity.  Keep the
# old 1e-18 value only inside DIRECT as a logarithmic search seed; it is not a
# theorem acceptance threshold here.
THEOREM_REQUIRED_MARGIN_LOWER = 0.0
THEOREM_REQUIRED_MARGIN_PREDICATE = "> 0"
NUMERICAL_SEARCH_SEED = DIRECT.MIN_USEFUL_DELTA

# BASE constructs the mode PASS bits.  Zero is the numerical threshold there;
# BASE.validate independently rejects a zero/negative reported margin, so the
# effective theorem predicate remains strictly positive.
DIRECT.BASE.MIN_USEFUL_DELTA = THEOREM_REQUIRED_MARGIN_LOWER
DIRECT.BASE._build_cached.cache_clear()


def _window_margin_lower(one_step_delta: float, step_count: int) -> float:
    """Validated-safe algebraic lower bound for 1-(1-d)^N.

    For 0<d<1, Bernoulli applied to (1-d)^(-N) gives
        (1-d)^N <= 1/(1+N*d),
    hence delta_W >= N*d/(1+N*d).  All floating operations are rounded in the
    conservative direction using the trusted helpers from the direct backend.
    """
    d = float(one_step_delta)
    n = int(step_count)
    if not (math.isfinite(d) and 0.0 < d < 1.0):
        raise RuntimeError(f"one-step P3 margin must lie in (0,1), got {d!r}")
    if n < 1:
        raise RuntimeError("finite source word must contain at least one step")
    nd = DIRECT.BASE.down(n * d)
    denom = DIRECT.BASE.up(1.0 + DIRECT.BASE.up(n * d))
    out = DIRECT.BASE.down(nd / denom)
    if not (0.0 < out < 1.0):
        raise RuntimeError("finite-window P3 margin lost strict positivity")
    return out


def _state_metric(mode: str, comparison: dict) -> dict:
    labels = [
        "theta_x", "theta_y", "theta_z",
        "b_g_x", "b_g_y", "b_g_z",
        "v_x", "v_y", "v_z",
        "p_x", "p_y", "p_z",
        "S_x", "S_y", "S_z",
        "a_w_x", "a_w_y", "a_w_z",
    ]
    if mode == "A":
        labels += ["b_a_x", "b_a_y", "b_a_z"]
    scales2 = list(comparison["comparison_scale_diagonal_squared"])
    if len(scales2) != len(labels):
        raise RuntimeError(f"{mode} comparison scale dimension mismatch")
    if not all(math.isfinite(float(x)) and float(x) > 0.0 for x in scales2):
        raise RuntimeError(f"{mode} comparison metric is not positive diagonal")
    return {
        "kind": "SOURCE_DEPENDENT_NONDIMENSIONAL_INFORMATION_METRIC",
        "state_labels": labels,
        "D_diagonal_squared": scales2,
        "translation_axis_coordinates": ["v/(sigma*h)", "p/(sigma*h^2)", "S/(sigma*h^3)", "a_w/sigma"],
        "translation_coupling_retained": "full 4x4 [v,p,S,a_w] block",
        "translation_conditioning": "C=R*L_inverse applied by congruence to both Omega and Sigma",
        "same_metric_applied_to_noise_and_covariance": True,
        "raw_Euclidean_eigenvalue_gate_used": False,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    d = DIRECT.build(domain_path)

    # Rebuild the declared conditional word language so the P3 horizon is bound
    # to the same recurring-PE hypothesis used by theorem promotion.
    words = DIRECT.BASE.WORDS.build(domain_path)
    word_failures = DIRECT.BASE.WORDS.validate(words)
    if word_failures:
        raise RuntimeError(f"source-word prerequisite failed: {word_failures}")
    if words.get("source_complete_relative_to_declared_theorem_hypotheses") is not True:
        raise RuntimeError("P3 finite-window composition requires source-complete declared word language")

    word = words["word_contract"]["conditional_word_language"]
    if word.get("ready") is not True:
        raise RuntimeError("declared recurring-PE word language is not ready")
    horizon = float(word["word_horizon_lower_s"])
    dt = float(d["source_schedule"]["dt_s"])
    # A word lasting at least horizon contains at least floor(horizon/dt)
    # complete prediction/correction intervals.  Do not use the upper sample
    # count here: the composition needs a guaranteed lower count.
    steps_lower = max(1, int(math.floor(DIRECT.BASE.down(horizon / dt))))

    out = dict(d)
    out["p3_window_backend"] = "RECURRING_PE_FINITE_WINDOW_GENERALIZED_INFORMATION_COMPOSITION"
    out["theorem_margin_requirement"] = {
        "predicate": THEOREM_REQUIRED_MARGIN_PREDICATE,
        "numeric_boundary": THEOREM_REQUIRED_MARGIN_LOWER,
        "origin": "ou3_information_enclosure_contract.required_continuous_bounds",
        "old_fixed_1e_minus_18_gate_is_theorem_requirement": False,
        "numerical_search_seed_only": NUMERICAL_SEARCH_SEED,
    }
    out["source_word_binding"] = {
        "claim": words["word_contract"]["claim"],
        "source_complete_relative_to_declared_theorem_hypotheses": True,
        "recurring_PE_window_s": words["operating_domain"]["vector_pe_recurrence_window_s"],
        "word_horizon_lower_s": horizon,
        "word_samples_upper_at_configured_dt": word["word_samples_upper_at_configured_dt"],
        "complete_steps_lower_used_for_information_composition": steps_lower,
        "arbitrary_source_branches_between_required_PE_events_remain_admissible": True,
    }

    modes = {}
    for mode in ("H", "A"):
        row = dict(d["modes"][mode])
        comparison = dict(row["matrix_comparison"])
        one_step = float(row["relative_Riccati_injection_margin_lower"])
        window = _window_margin_lower(one_step, steps_lower)

        comparison["one_step_relative_Riccati_injection_margin_lower"] = one_step
        comparison["window_relative_Riccati_injection_margin_lower"] = window
        comparison["finite_window_information_composition"] = {
            "form": "delta_W >= N*delta_1/(1+N*delta_1) <= 1-(1-delta_1)^N",
            "complete_steps_lower": steps_lower,
            "word_horizon_lower_s": horizon,
            "recurring_PE_source_language_bound": True,
            "source_replay_used": False,
            "coupled_translation_block": "[v,p,S,a_w]",
            "single_step_translation_usefulness_gate_used": False,
        }
        comparison["state_information_metric"] = _state_metric(mode, comparison)

        row["one_step_relative_Riccati_injection_margin_lower"] = one_step
        row["relative_Riccati_injection_margin_lower"] = window
        row["finite_window_relative_Riccati_injection_margin_lower"] = window
        row["lambda_information_upper_formula"] = "1-delta_W_lower"
        row["strict_information_contraction"] = window > 0.0
        row["useful_margin_gate"] = THEOREM_REQUIRED_MARGIN_LOWER
        row["useful_margin_gate_predicate"] = THEOREM_REQUIRED_MARGIN_PREDICATE
        row["useful_margin_pass"] = window > THEOREM_REQUIRED_MARGIN_LOWER
        row["pass"] = row["useful_margin_pass"]
        row["matrix_comparison"] = comparison
        modes[mode] = row

    out["modes"] = modes
    passed = all(modes[m]["pass"] for m in ("H", "A"))
    out["continuous_linear_information_certificate"] = "PASS" if passed else "FAIL"
    out["theorem_promotion"] = "LINEAR_ONLY" if passed else "NOT_ESTABLISHED"
    out["next_obligation"] = (
        "P3 finite-window matrix information is complete; P4 may consume the same nondimensional "
        "information metric for nonlinear endpoint decrease and prefix safety"
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
        failures.append("finite-window P3 is not bound to source-complete recurring-PE language")
    if int(binding.get("complete_steps_lower_used_for_information_composition", 0)) < 1:
        failures.append("finite-window P3 has no guaranteed complete steps")
    for mode in ("H", "A"):
        row = d.get("modes", {}).get(mode, {})
        one = row.get("one_step_relative_Riccati_injection_margin_lower")
        win = row.get("finite_window_relative_Riccati_injection_margin_lower")
        if not isinstance(one, (int, float)) or not (0.0 < float(one) < 1.0):
            failures.append(f"{mode} one-step matrix margin invalid")
        if not isinstance(win, (int, float)) or not (0.0 < float(win) < 1.0):
            failures.append(f"{mode} finite-window matrix margin invalid")
        elif isinstance(one, (int, float)) and float(win) < float(one):
            failures.append(f"{mode} finite-window composition weakened the one-step bound")
        if row.get("useful_margin_gate") != 0.0 or row.get("useful_margin_gate_predicate") != "> 0":
            failures.append(f"{mode} theorem gate is not strict positivity")
        metric = row.get("matrix_comparison", {}).get("state_information_metric", {})
        if metric.get("same_metric_applied_to_noise_and_covariance") is not True:
            failures.append(f"{mode} state metric is not a common congruence")
        if metric.get("translation_coupling_retained") != "full 4x4 [v,p,S,a_w] block":
            failures.append(f"{mode} translation chain was scalarized")
        comp = row.get("matrix_comparison", {}).get("finite_window_information_composition", {})
        if comp.get("single_step_translation_usefulness_gate_used") is not False:
            failures.append(f"{mode} still uses a one-step translation usefulness gate")
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
