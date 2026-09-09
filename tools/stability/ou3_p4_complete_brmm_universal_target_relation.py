#!/usr/bin/env python3
"""Coarse shipping target enclosure for the complete-BRMM P4 cover.

This module is deliberately NOT a theorem-facing coefficient closure.

Shipping clamps do imply a compact enclosure for tuning frequency and raw sigma,
and tau/T_S/SpectralMSE-R_S are deterministic dependent images once one
(frequency,sigma) pair is supplied. But P4 must not treat frequency and sigma
as independent coordinates: both are produced from the same physical signal
history.

The proof-facing source relation must retain the joint estimator state:

  raw IMU / private complementary vertical acceleration
    -> WavePeriodEstimator HP/leaky-integrator/EW moments
    -> canonical period T and f=1/T
    -> period-scaled AdaptiveWaveBandPass driven by that same f
    -> EW wave-band variance and noise-gain subtraction
    -> sigma_target_raw
    -> tau(f), T_S(tau), SpectralMSE R_S(tau,sigma,T_S).

The clamp rectangle below is only a sanity enclosure and may never promote P4
or satisfy coefficient-source inclusion by itself.
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
SCHEMA=2
QUALIFICATION="OU3_P4_COMPLETE_BRMM_COARSE_TARGET_ENCLOSURE_V2"


def build(domain_path:Path=DEFAULT_DOMAIN)->dict:
    dynamic=DYNAMIC.build(domain_path); df=DYNAMIC.validate(dynamic)
    target=TARGET.build(domain_path); tf=TARGET.validate(target)
    adapt=ADAPT.build(domain_path); af=ADAPT.validate(adapt)
    bad={k:v for k,v in (("dynamic",df),("target",tf),("adaptive",af)) if v}
    if bad: raise RuntimeError("coarse target prerequisites failed: "+repr(bad))
    parity=dynamic["source_parity"]
    inv=dynamic["dynamic_invariant"]
    return {
        "schema":SCHEMA,"qualification":QUALIFICATION,
        "canonical_source":"COMPLETE_BRMM_NORMAL_LIVE_WORD",
        "finite_frontend_history_required":True,
        "all_shipping_tuning_frequency_outputs_enclosed":bool(parity["normal_live_measured_period_selector_present"]),
        "all_shipping_raw_sigma_targets_enclosed":bool(parity["raw_sigma_target_uses_1e_minus6_variance_floor"]),
        "quiet_continuation_enclosed":bool(dynamic["normal_live_contract"]["quiet_case_admitted"]),
        "frequency_raw_sigma_target_rectangle":{
            "frequency_hz":inv["tuning_frequency_hz"],
            "sigma_target_raw_mps2":inv["sigma_target_raw_mps2"],
        },
        "tau_TS_RS_are_correlated_same_cell_images":bool(target["full_dynamic_target_rectangle_image_valid"]),
        "SpectralMSE_RS_same_cell_map":bool(target["SpectralMSE_RS_is_same_tau_sigma_image"]),
        "independent_RS_target_forbidden":bool(target["independent_RS_target_forbidden"]),
        "active_schedule_connected_by_shipping_EMA":bool(adapt["EMA_candidate_transition_available"]),
        "active_schedule_connected_by_staged_commit":bool(adapt["staged_commit_transition_available"]),
        "pseudo_event_incidence_connected_by_scheduler":bool(adapt["scheduler_due_not_due_branch_preserved"]),
        "raw_and_effective_sigma_coordinates_preserved":bool(adapt["raw_and_effective_sigma_coordinates_preserved"]),
        "full_rectangle_is_coarse_sanity_enclosure_only":True,
        "target_pair_may_span_full_rectangle_in_theorem":False,
        "independent_frequency_sigma_coordinates_forbidden":True,
        "same_signal_history_must_generate_period_and_sigma":True,
        "wave_period_moments_must_be_retained":True,
        "period_scaled_sigma_band_must_be_retained":True,
        "sigma_variance_statistic_must_be_retained":True,
        "period_frequency_reciprocal_relation_must_be_retained":True,
        "joint_estimator_relation_materialized_here":False,
        "coefficient_target_inclusion_closed":False,
        "physical_vector_geometry_transition_closed_here":False,
        "complete_source_cover_closed_here":False,
        "P4_promoted_here":False,
    }


def validate(d:dict)->list[str]:
    f=[]
    if d.get("schema")!=SCHEMA or d.get("qualification")!=QUALIFICATION:f.append("schema/qualification mismatch")
    for k in ("finite_frontend_history_required","all_shipping_tuning_frequency_outputs_enclosed","all_shipping_raw_sigma_targets_enclosed","quiet_continuation_enclosed","tau_TS_RS_are_correlated_same_cell_images","SpectralMSE_RS_same_cell_map","independent_RS_target_forbidden","active_schedule_connected_by_shipping_EMA","active_schedule_connected_by_staged_commit","pseudo_event_incidence_connected_by_scheduler","raw_and_effective_sigma_coordinates_preserved","full_rectangle_is_coarse_sanity_enclosure_only","independent_frequency_sigma_coordinates_forbidden","same_signal_history_must_generate_period_and_sigma","wave_period_moments_must_be_retained","period_scaled_sigma_band_must_be_retained","sigma_variance_statistic_must_be_retained","period_frequency_reciprocal_relation_must_be_retained"):
        if d.get(k) is not True:f.append(k+" not true")
    for k in ("target_pair_may_span_full_rectangle_in_theorem","joint_estimator_relation_materialized_here","coefficient_target_inclusion_closed","physical_vector_geometry_transition_closed_here","complete_source_cover_closed_here","P4_promoted_here"):
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
    print(json.dumps({"coefficient_target_inclusion":d["coefficient_target_inclusion_closed"],"joint_estimator_relation":d["joint_estimator_relation_materialized_here"],"failures":f},sort_keys=True))
    return int(bool(f))
if __name__=="__main__":raise SystemExit(main())
