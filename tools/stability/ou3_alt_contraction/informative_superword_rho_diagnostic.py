#!/usr/bin/env python3
"""Non-promoting point diagnostic for the selected ALT magnetic superword target.

This experiment deliberately reuses the literal joint24 event composer from
finite_word_rho_diagnostic.  It does not manufacture a new filter or fit replay
trajectories.  It asks the first cheap question required by AGENTS.md after the
observability reformulation: once a word contains recurrent informative magnetic
updates, is strict contraction even numerically plausible?

Important limitation: the current point composer is an analytic same-history
member of the admitted source language, not the completed source-uniform finite
master.  Therefore every result here is feasibility evidence only and cannot
promote storage_search_allowed or any ALT PASS gate.
"""
from __future__ import annotations
import numpy as np

from tools.stability.ou3_alt_contraction import finite_word_rho_diagnostic as OLD
from tools.stability.ou3_alt_contraction import magnetic_service_formulation as FORM

QUALIFICATION="OU3_ALT_INFORMATIVE_MAGNETIC_SUPERWORD_RHO_DIAGNOSTIC_V1"
SAMPLES=OLD.CANONICAL_WORD_SAMPLES


def _gauged_phases():
    return tuple(p for p in OLD.SOURCE_PHASES if p.gauged)


def _identity_metric():
    return np.eye(FORM.JOINT_DIM)


def _pair_information(samples:int, h:float, stride:int=OLD.GAUGED_MAG_STRIDE):
    """Transported scalar-heading proxy [1,t] used only to audit point service.

    The literal word itself contains the full 3-D magnetic Jacobians.  This proxy
    is not used in the rho calculation; it only verifies that the represented
    accepted service times span both heading and axial-bias directions rather
    than counting calls as information.
    """
    events=[]
    for k in range(samples):
        if k%stride: continue
        t=k*h
        events.append(FORM.MagneticEvent(t,True,True,np.array([[1.0,t]])))
    return events


def diagnose(samples:int=SAMPLES):
    if samples<=0: raise ValueError("positive sample count required")
    Z=FORM.joint24_motion_injection(); M=_identity_metric()
    rows={}
    for mode in OLD.MODES:
        for phase in _gauged_phases():
            word=OLD.compose_word(phase,samples,mode)
            A=word["A_joint24"]
            rho=FORM.projected_storage_ratio(A,M,M,Z)
            h=word["word_horizon_s"]/samples
            events=_pair_information(samples,h)
            # Point service contract follows the represented 25-Hz accepted
            # events.  This does NOT promote the runtime theorem assumption.
            contract=FORM.MagneticServiceContract(
                max_informative_gap_s=OLD.GAUGED_MAG_STRIDE*h+1e-12,
                min_heading_response=0.5,
                min_pair_information=1e-6)
            service=FORM.audit_magnetic_service(events,0.0,word["word_horizon_s"],contract)
            rows[f"{mode}:{phase.name}"]={
                "mode":mode,"phase":phase.name,"word_horizon_s":word["word_horizon_s"],
                "literal_event_count":word["literal_event_count"],
                "service":service,"identity_storage_ratio":rho,
                "same_literal_joint24_composer_as_PR533":True,
                "source_uniform_finite_master":False,
            }
    worst=max(rows,key=lambda k:rows[k]["identity_storage_ratio"]["rho_point"])
    worst_rho=rows[worst]["identity_storage_ratio"]["rho_point"]
    return {
        "qualification":QUALIFICATION,"samples":samples,"rows":rows,
        "worst_word":worst,"worst_identity_storage_rho_point":worst_rho,
        "identity_storage_strictly_contracts_all_measured_words":worst_rho<1.0,
        "corrected_formulation_falsified_by_this_experiment":worst_rho>=1.0,
        "common_or_compatible_storage_search_still_needed":True,
        "source_uniform_rho_certified":False,"interval_enclosure_authorized":False,
        "storage_search_allowed":False,"ALT_LIVE_PASS":False,
    }


def validate(d):
    f=[]
    if d.get("qualification")!=QUALIFICATION:f.append("qualification mismatch")
    rows=d.get("rows",{})
    if set(rows)!={"H:quiet_gauged","H:wave_gauged","A:quiet_gauged","A:wave_gauged"}:f.append("gauged family mismatch")
    for name,row in rows.items():
        if not row.get("service",{}).get("finite_window_point_service_pass"):f.append(name+" point informative service failed")
        if row.get("source_uniform_finite_master") is not False:f.append(name+" incorrectly promoted source uniformity")
        if row.get("identity_storage_ratio",{}).get("source_uniform_rho_certified") is not False:f.append(name+" rho incorrectly certified")
    for k in ("source_uniform_rho_certified","interval_enclosure_authorized","storage_search_allowed","ALT_LIVE_PASS"):
        if d.get(k) is not False:f.append(k+" not false")
    return f


if __name__=="__main__":
    import json
    print(json.dumps(diagnose(),indent=2))
