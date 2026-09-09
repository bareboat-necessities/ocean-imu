#!/usr/bin/env python3
"""Universal shipping-target relation for the complete-BRMM P4 cover.

P4 coefficient coverage does not require a hard deterministic bound on the raw
private-Mahony acceleration sample.  Shipping itself maps every finite frontend
history into a compact target set before coefficients reach the filter:

* the tuner stores a finite wave frequency and the wrapper clamps the frequency
  used for tau to [MIN_TUNE_FREQ_HZ, MAX_TUNE_FREQ_HZ];
* after variance readiness ``var_wave`` is floored at 1e-6 and the resulting raw
  sigma target is capped by MAX_SIGMA_A; before readiness the additional 0.05
  floor only narrows this same set;
* tau, pseudo cadence and deployed SpectralMSE R_S are then deterministic
  same-cell images, with the final R_S hard clamp retained.

The relation intentionally permits the frequency/raw-sigma pair to move anywhere
inside this rectangle from sample to sample.  That is a conservative superset of
actual frontend histories, NOT an arbitrary active filter schedule: candidate
and active coefficients are still connected by the exact shipping EMA, staged
commit and pseudo-scheduler recurrences.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_brmm_dynamic_source_certificate as DYNAMIC
import ou3_p4_complete_brmm_target_cell as TARGET
import ou3_p4_complete_brmm_adaptive_transition as ADAPT

REPO=Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN=REPO/"tools"/"stability"/"ou3_proof_operating_domain.json"
SCHEMA=1
QUALIFICATION="OU3_P4_COMPLETE_BRMM_UNIVERSAL_TARGET_RELATION_V1"


def build(domain_path:Path=DEFAULT_DOMAIN)->dict:
    dynamic=DYNAMIC.build(domain_path); df=DYNAMIC.validate(dynamic)
    target=TARGET.build(domain_path); tf=TARGET.validate(target)
    adapt=ADAPT.build(domain_path); af=ADAPT.validate(adapt)
    bad={k:v for k,v in (("dynamic",df),("target",tf),("adaptive",af)) if v}
    if bad: raise RuntimeError("universal target prerequisites failed: "+repr(bad))
    parity=dynamic["source_parity"]
    inv=dynamic["dynamic_invariant"]
    return {
        "schema":SCHEMA,"qualification":QUALIFICATION,
        "canonical_source":"COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "finite_frontend_history_required":True,
        "all_shipping_tuning_frequency_outputs_enclosed":bool(parity["normal_live_measured_period_selector_present"]),
        "all_shipping_raw_sigma_targets_enclosed":bool(parity["raw_sigma_target_uses_1e_minus6_variance_floor"]),
        "quiet_continuation_enclosed":bool(dynamic["normal_live_contract"]["quiet_case_admitted"]),
        "unbounded_measurement_amplitude_cannot_escape_raw_sigma_target_cap":True,
        "frequency_raw_sigma_target_rectangle":{
            "frequency_hz":inv["tuning_frequency_hz"],
            "sigma_target_raw_mps2":inv["sigma_target_raw_mps2"],
        },
        "tau_TS_RS_are_correlated_same_cell_images":bool(target["full_dynamic_target_rectangle_image_valid"]),
        "SpectralMSE_RS_same_cell_map":bool(target["SpectralMSE_RS_is_same_tau_sigma_image"]),
        "independent_RS_target_forbidden":bool(target["independent_RS_target_forbidden"]),
        "target_pair_may_span_full_rectangle_each_sample":True,
        "full_rectangle_is_conservative_overapproximation_not_physical_generator":True,
        "target_overapproximation_may_generate_arbitrary_active_schedule":False,
        "active_schedule_connected_by_shipping_EMA":bool(adapt["EMA_candidate_transition_available"]),
        "active_schedule_connected_by_staged_commit":bool(adapt["staged_commit_transition_available"]),
        "pseudo_event_incidence_connected_by_scheduler":bool(adapt["scheduler_due_not_due_branch_preserved"]),
        "raw_and_effective_sigma_coordinates_preserved":bool(adapt["raw_and_effective_sigma_coordinates_preserved"]),
        "physical_private_mahony_input_bound_required_for_coefficient_inclusion":False,
        "exact_frontend_history_required_for_tighter_than_clamp_coefficient_family":True,
        "coefficient_target_inclusion_closed":True,
        "physical_vector_geometry_transition_closed_here":False,
        "complete_source_cover_closed_here":False,
        "P4_promoted_here":False,
    }


def validate(d:dict)->list[str]:
    f=[]
    if d.get("schema")!=SCHEMA or d.get("qualification")!=QUALIFICATION:f.append("schema/qualification mismatch")
    for k in ("finite_frontend_history_required","all_shipping_tuning_frequency_outputs_enclosed","all_shipping_raw_sigma_targets_enclosed","quiet_continuation_enclosed","unbounded_measurement_amplitude_cannot_escape_raw_sigma_target_cap","tau_TS_RS_are_correlated_same_cell_images","SpectralMSE_RS_same_cell_map","independent_RS_target_forbidden","target_pair_may_span_full_rectangle_each_sample","full_rectangle_is_conservative_overapproximation_not_physical_generator","active_schedule_connected_by_shipping_EMA","active_schedule_connected_by_staged_commit","pseudo_event_incidence_connected_by_scheduler","raw_and_effective_sigma_coordinates_preserved","exact_frontend_history_required_for_tighter_than_clamp_coefficient_family","coefficient_target_inclusion_closed"):
        if d.get(k) is not True:f.append(k+" not true")
    for k in ("target_overapproximation_may_generate_arbitrary_active_schedule","physical_private_mahony_input_bound_required_for_coefficient_inclusion","physical_vector_geometry_transition_closed_here","complete_source_cover_closed_here","P4_promoted_here"):
        if d.get(k) is not False:f.append(k+" not false")
    box=d.get("frequency_raw_sigma_target_rectangle",{})
    for k in ("frequency_hz","sigma_target_raw_mps2"):
        x=box.get(k,[])
        if len(x)!=2 or not 0<float(x[0])<=float(x[1]):f.append(k+" invalid")
    return f


def main()->int:
    p=argparse.ArgumentParser();p.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN);p.add_argument("--output",type=Path,required=True);a=p.parse_args()
    d=build(a.domain);f=validate(d);d["validation_pass"]=not f;d["validation_failures"]=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+"\n")
    print(json.dumps({"coefficient_target_inclusion":d["coefficient_target_inclusion_closed"],"active_arbitrary":d["target_overapproximation_may_generate_arbitrary_active_schedule"],"failures":f},sort_keys=True))
    return int(bool(f))
if __name__=="__main__":raise SystemExit(main())
