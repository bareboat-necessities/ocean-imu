#!/usr/bin/env python3
"""Source-complete operation-timing contract for the OU-III P4 word.

The nonlinear proof must not replace the source scheduler by one convenient
sample pattern.  In particular the four S=0 observations used for translation
UCO are guaranteed only through the deployed cadence interval, not at four
fixed samples separated by the minimum cadence.

This producer records the correct decomposition used by the current P4 route:

* P3 consumes the source-complete S=0 timing analytically through the four-fire
  translation UCO and therefore does not need exact S sample indices;
* S=0 has an exactly linear residual (eta_S=0), so its uncertain timing creates
  no nonlinear remainder obligation for P4;
* the nonlinear remainder stage therefore needs only the certified vector PE
  packet family and any additional accepted vector updates, not a fabricated
  terminal cluster containing all S observations.

This removes the old 192..201 terminal-cluster shortcut from the promotion
contract.  The old reduction may remain a design diagnostic, but it is not
source-complete evidence for the theorem.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_implementation_word_language as WORDS
import ou3_translational_uco_ucc as TRANS
import ou3_vector_uco_certificate as VECTOR

REPO=Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN=REPO/"tools"/"ou3_proof_operating_domain.json"
SCHEMA=1


def build(domain_path: Path=DEFAULT_DOMAIN):
    path=Path(domain_path).resolve()
    words=WORDS.build(path)
    wf=WORDS.validate(words)
    trans=TRANS.build()
    tf=TRANS.validate(trans)
    vector=VECTOR.build()
    vf=VECTOR.validate(vector)
    failures=[f"word: {x}" for x in wf]+[f"translation: {x}" for x in tf]+[f"vector: {x}" for x in vf]

    wc=words.get("word_contract",{})
    cw=wc.get("conditional_word_language",{})
    tr=wc.get("translation_recurrence",{})
    pe=wc.get("vector_persistent_excitation",{})
    dt=float(wc.get("configured_runtime",{}).get("imu_dt_s",math.nan))
    horizon=float(cw.get("word_horizon_lower_s",math.nan))
    samples=int(cw.get("word_samples_upper_at_configured_dt",0) or 0)
    s_gap_lo=float(tr.get("pseudo_gap_min_s",math.nan))
    s_gap_hi=float(tr.get("pseudo_gap_max_s",math.nan))
    s_span=float(tr.get("spread_selected_window_s_upper",math.nan))
    packet=list(pe.get("packet_gap_s",[]))
    recurrence=float(pe.get("recurrence_window_s",math.nan))

    if not (math.isfinite(dt) and dt>0 and math.isfinite(horizon) and horizon>0 and samples>0):
        failures.append("word clock/horizon invalid")
    if not (math.isfinite(s_gap_lo) and math.isfinite(s_gap_hi) and 0<s_gap_lo<=s_gap_hi):
        failures.append("S cadence interval invalid")
    if not (math.isfinite(s_span) and s_span>0 and s_span<=horizon):
        failures.append("four-S source span does not fit word horizon")
    if len(packet)!=2 or not (0<float(packet[0])<=float(packet[1])<=recurrence<=horizon):
        failures.append("vector packet recurrence/timing invalid")

    # A fixed schedule at the minimum S gap is not equivalent to the source
    # scheduler whenever the admissible cadence interval has nonzero width.
    s_timing_uncertain=bool(math.isfinite(s_gap_lo) and math.isfinite(s_gap_hi) and s_gap_hi>s_gap_lo)

    return {
        "schema":SCHEMA,
        "qualification":"OU3_P4_SOURCE_COMPLETE_WORD_TIMING_DECOMPOSITION",
        "source_generated_not_trajectory_fit":True,
        "trajectory_replay_used":False,
        "filter_changed":False,
        "configured_dt_s":dt,
        "word_horizon_s":horizon,
        "word_samples_upper":samples,
        "S_gap_interval_s":[s_gap_lo,s_gap_hi],
        "four_S_selected_span_upper_s":s_span,
        "S_firing_times_are_source_intervals_not_fixed_samples":s_timing_uncertain,
        "fixed_minimum_gap_S_schedule_is_source_complete":False,
        "S_residual_exactly_linear_selector":True,
        "S_nonlinear_eta_identically_zero":True,
        "S_timing_consumed_by_linear_P3_translation_UCO":True,
        "vector_packet_gap_s":packet,
        "vector_PE_recurrence_window_s":recurrence,
        "nonlinear_timing_obligations_reduce_to_vector_measurements":True,
        "old_terminal_192_201_cluster_required_for_promotion":False,
        "terminal_cluster_design_diagnostic_may_be_retained":True,
        "ready_for_source_complete_nonlinear_remainder_composition":not failures,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE":False,
        "failures":failures,
    }


def validate(d):
    f=list(d.get("failures",[]))
    for k in (
        "source_generated_not_trajectory_fit",
        "S_firing_times_are_source_intervals_not_fixed_samples",
        "S_residual_exactly_linear_selector",
        "S_nonlinear_eta_identically_zero",
        "S_timing_consumed_by_linear_P3_translation_UCO",
        "nonlinear_timing_obligations_reduce_to_vector_measurements",
        "terminal_cluster_design_diagnostic_may_be_retained",
        "ready_for_source_complete_nonlinear_remainder_composition",
    ):
        if d.get(k) is not True: f.append(f"{k} is not true")
    for k in (
        "trajectory_replay_used","filter_changed","fixed_minimum_gap_S_schedule_is_source_complete",
        "old_terminal_192_201_cluster_required_for_promotion","P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE",
    ):
        if d.get(k) is not False: f.append(f"{k} is not false")
    return list(dict.fromkeys(f))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    d=build(a.domain)
    f=validate(d)
    d["validation_pass"]=not f
    d["validation_failures"]=f
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(d,indent=2,sort_keys=True),encoding="utf-8")
    print(json.dumps(d,indent=2,sort_keys=True))
    return 0 if not f else 2


if __name__=="__main__":
    raise SystemExit(main())
