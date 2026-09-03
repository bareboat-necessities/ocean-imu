#!/usr/bin/env python3
"""Source-complete operation-timing contract for the OU-III P4 word.

The nonlinear proof must not replace the source scheduler by one convenient
sample pattern.  In particular the four S=0 observations used for translation
UCO are guaranteed only through the deployed cadence interval, not at four
fixed samples separated by the minimum cadence.

P4 must also use the *same word boundary* as canonical P3.  The generic
implementation-word language only needs a window long enough to contain the
translation and vector recurrence ingredients, and therefore uses their
maximum.  Canonical P3 is stricter: its covariance comparison first covers the
translation observation horizon and then allows one complete vector-PE
recurrence before the endpoint packet.  Its source-uniform word length is

    gap     = cadence_max + dt,
    spacing = max(T_PE, 2 gap),
    T_obs   = 2 spacing + gap,
    T_word  = T_obs + T_PE.

That sequential word is the boundary on which the P3 H/A conditional-process
floor is established.  A shorter P4 word cannot consume that metric.  This
producer therefore derives ``T_word`` through the exact same retained P3
``translation_upper`` primitive and exposes the old implementation-language
horizon only as a diagnostic lower requirement.

The remaining decomposition is unchanged:

* P3 consumes the source-complete S=0 timing analytically through the four-fire
  translation UCO and therefore does not need exact S sample indices;
* S=0 has an exactly linear residual (eta_S=0), so its uncertain timing creates
  no nonlinear remainder obligation for P4; and
* nonlinear P4 timing therefore needs only the certified vector PE packet
  family and any additional accepted vector updates, not a fabricated terminal
  cluster containing all S observations.

No trajectory values, domain shrink, filter changes, or theorem promotion are
introduced here.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval
import ou3_implementation_word_language as WORDS
import ou3_source_reachable_matrix_p3 as P3
import ou3_translational_uco_ucc as TRANS
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2


def _canonical_p3_word(path: Path, recurrence: float) -> dict:
    """Derive the exact source-uniform covariance word used to size P3 history."""
    sched = P3.source_schedule()
    tau = Interval(*map(float, sched["tau_applied_invariant_s"]))
    sigma = Interval(*map(float, sched["sigma_aw_applied_safety"]))
    rs = Interval(*map(float, sched["R_S_applied_invariant"]))
    _upper, timing = P3.translation_upper(tau, sigma, rs, float(recurrence), sched)
    return {
        "dt_s": float(sched["dt_s"]),
        "word_horizon_s_lower": float(timing["word_horizon_s_lower"]),
        "word_horizon_s_upper": float(timing["word_horizon_s_upper"]),
        "gap_s_upper": float(timing["gap_s_upper"]),
        "cadence_s": list(map(float, timing["cadence_s"])),
        "derivation": "ou3_source_reachable_matrix_p3.translation_upper",
        "same_formula_used_by_P2_V1_history_target": True,
    }


def build(domain_path: Path = DEFAULT_DOMAIN):
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("P4 source-word timing must not be trajectory fitted")

    words = WORDS.build(path)
    wf = WORDS.validate(words)
    trans = TRANS.build()
    tf = TRANS.validate(trans)
    vector = VECTOR.build()
    vf = VECTOR.validate(vector)
    failures = (
        [f"word: {x}" for x in wf]
        + [f"translation: {x}" for x in tf]
        + [f"vector: {x}" for x in vf]
    )

    wc = words.get("word_contract", {})
    cw = wc.get("conditional_word_language", {})
    tr = wc.get("translation_recurrence", {})
    pe = wc.get("vector_persistent_excitation", {})
    dt = float(wc.get("configured_runtime", {}).get("imu_dt_s", math.nan))
    implementation_horizon = float(cw.get("word_horizon_lower_s", math.nan))
    implementation_samples = int(cw.get("word_samples_upper_at_configured_dt", 0) or 0)
    s_gap_lo = float(tr.get("pseudo_gap_min_s", math.nan))
    s_gap_hi = float(tr.get("pseudo_gap_max_s", math.nan))
    s_span = float(tr.get("spread_selected_window_s_upper", math.nan))
    packet = list(pe.get("packet_gap_s", []))
    recurrence = float(pe.get("recurrence_window_s", math.nan))

    p3_word = None
    if math.isfinite(recurrence) and recurrence > 0.0:
        try:
            p3_word = _canonical_p3_word(path, recurrence)
        except Exception as exc:
            failures.append(f"canonical P3 word: {type(exc).__name__}: {exc}")

    if not (
        math.isfinite(dt) and dt > 0.0
        and math.isfinite(implementation_horizon) and implementation_horizon > 0.0
        and implementation_samples > 0
    ):
        failures.append("implementation word clock/horizon invalid")
    if not (math.isfinite(s_gap_lo) and math.isfinite(s_gap_hi) and 0.0 < s_gap_lo <= s_gap_hi):
        failures.append("S cadence interval invalid")
    if len(packet) != 2 or not (
        0.0 < float(packet[0]) <= float(packet[1]) <= recurrence
    ):
        failures.append("vector packet recurrence/timing invalid")

    p3_horizon = math.nan
    p3_horizon_lower = math.nan
    samples = 0
    if p3_word is not None:
        p3_horizon = float(p3_word["word_horizon_s_upper"])
        p3_horizon_lower = float(p3_word["word_horizon_s_lower"])
        if not (math.isfinite(p3_horizon) and p3_horizon > 0.0):
            failures.append("canonical P3 covariance word horizon invalid")
        if not math.isclose(float(p3_word["dt_s"]), dt, rel_tol=0.0, abs_tol=0.0):
            failures.append("P3 source schedule dt differs from implementation word dt")
        if math.isfinite(p3_horizon) and p3_horizon > 0.0 and math.isfinite(dt) and dt > 0.0:
            # Same one-sample endpoint padding convention used by the P2-V1
            # history quotient.  The nonlinear word may terminate earlier on a
            # realized packet, but the certificate must cover this upper count.
            samples = int(math.ceil(p3_horizon / dt - 1.0e-14)) + 1

    if not (math.isfinite(s_span) and s_span > 0.0 and s_span <= p3_horizon):
        failures.append("four-S source span does not fit canonical P3 word horizon")
    if not (math.isfinite(recurrence) and 0.0 < recurrence <= p3_horizon):
        failures.append("vector recurrence does not fit canonical P3 word horizon")
    if not (
        math.isfinite(implementation_horizon)
        and math.isfinite(p3_horizon)
        and implementation_horizon <= p3_horizon
    ):
        failures.append("implementation source-language horizon is not contained in canonical P3 word")

    # A fixed schedule at the minimum S gap is not equivalent to the source
    # scheduler whenever the admissible cadence interval has nonzero width.
    s_timing_uncertain = bool(
        math.isfinite(s_gap_lo) and math.isfinite(s_gap_hi) and s_gap_hi > s_gap_lo
    )

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_SOURCE_COMPLETE_WORD_TIMING_DECOMPOSITION",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "configured_dt_s": dt,
        "implementation_language_word_horizon_s": implementation_horizon,
        "implementation_language_word_samples_upper": implementation_samples,
        "canonical_P3_covariance_word": p3_word,
        "canonical_P3_word_horizon_s_lower": p3_horizon_lower,
        "canonical_P3_word_horizon_s_upper": p3_horizon,
        "word_horizon_s": p3_horizon,
        "word_samples_upper": samples,
        "P4_word_horizon_bound_to_canonical_P3_covariance_word": p3_word is not None,
        "implementation_language_is_only_a_lower_timing_requirement": True,
        "sequential_translation_then_vector_endpoint_allowance_used": True,
        "S_gap_interval_s": [s_gap_lo, s_gap_hi],
        "four_S_selected_span_upper_s": s_span,
        "S_firing_times_are_source_intervals_not_fixed_samples": s_timing_uncertain,
        "fixed_minimum_gap_S_schedule_is_source_complete": False,
        "S_residual_exactly_linear_selector": True,
        "S_nonlinear_eta_identically_zero": True,
        "S_timing_consumed_by_linear_P3_translation_UCO": True,
        "vector_packet_gap_s": packet,
        "vector_PE_recurrence_window_s": recurrence,
        "nonlinear_timing_obligations_reduce_to_vector_measurements": True,
        "old_terminal_192_201_cluster_required_for_promotion": False,
        "terminal_cluster_design_diagnostic_may_be_retained": True,
        "ready_for_source_complete_nonlinear_remainder_composition": not failures,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE": False,
        "failures": failures,
    }


def validate(d):
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "P4_word_horizon_bound_to_canonical_P3_covariance_word",
        "implementation_language_is_only_a_lower_timing_requirement",
        "sequential_translation_then_vector_endpoint_allowance_used",
        "S_firing_times_are_source_intervals_not_fixed_samples",
        "S_residual_exactly_linear_selector",
        "S_nonlinear_eta_identically_zero",
        "S_timing_consumed_by_linear_P3_translation_UCO",
        "nonlinear_timing_obligations_reduce_to_vector_measurements",
        "terminal_cluster_design_diagnostic_may_be_retained",
        "ready_for_source_complete_nonlinear_remainder_composition",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "trajectory_replay_used",
        "filter_changed",
        "fixed_minimum_gap_S_schedule_is_source_complete",
        "old_terminal_192_201_cluster_required_for_promotion",
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    lo = d.get("canonical_P3_word_horizon_s_lower")
    hi = d.get("canonical_P3_word_horizon_s_upper")
    if not all(isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(float(x)) for x in (lo, hi)):
        f.append("canonical P3 word horizon bounds are invalid")
    elif not (0.0 < float(lo) <= float(hi) == float(d.get("word_horizon_s"))):
        f.append("P4 word horizon is not the canonical P3 covariance-word upper")
    if int(d.get("word_samples_upper", 0) or 0) <= 0:
        f.append("P4 word sample upper is not positive")
    return list(dict.fromkeys(f))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain)
    f = validate(d)
    d["validation_pass"] = not f
    d["validation_failures"] = f
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(d, indent=2, sort_keys=True))
    return 0 if not f else 2


if __name__ == "__main__":
    raise SystemExit(main())
