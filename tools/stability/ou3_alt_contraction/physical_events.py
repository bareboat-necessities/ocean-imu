#!/usr/bin/env python3
"""Finite physical-reference Joseph events for the ALT Live theorem.

The shared differential event implementation is exact for state dependence but
its S pseudo-measurement residual is the differential selector e_S. The actual
finite physical innovation is 0-S_hat=e_S-S_phys. This module evaluates the
same deployed correction/Joseph geometry with S_phys from the SAME typed BRMM
primitive carried by the estimator-owned source cell.

No independent S box, replay value, or wordwise re-zeroing is accepted.
"""
from __future__ import annotations
from typing import Sequence
from ou3_interval import Interval
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_complete_brmm_source_cover_contract as COVER
from tools.stability.ou3_alt_contraction import physical_reference as REF

QUALIFICATION="OU3_ALT_FINITE_PHYSICAL_JOSEPH_EVENTS_V1"

def _physical_s_residual_ad(z,centered_s:Sequence[Interval]):
    if len(centered_s)!=3 or any(not isinstance(x,Interval) for x in centered_s): raise ValueError("same-history centered S must be Interval[3]")
    n=len(z);return [z[EVENTS.OFF_S+i]-EVENTS._const(centered_s[i],n) for i in range(3)]

def s_zero_event(cell:COVER.SourceCoverCell)->dict:
    failures=COVER.validate_cell(cell,require_estimator_provenance=True)
    if failures: raise ValueError("invalid estimator-owned source cell: "+repr(failures))
    if cell.kind!="S_zero": raise ValueError("physical S event requires S_zero source cell")
    if cell.wave_primitive is None: raise ValueError("S_zero event missing same-history physical primitive")
    REF.validate_primitive_reference(cell.wave_primitive)
    if cell.R_provenance!=EVENTS.ACTUAL_RS_PROVENANCE: raise ValueError("S_zero event lost actual applied R_S provenance")
    mode=cell.mode;H=EVENTS._event_H(mode,"S_zero");K,S=EVENTS.source_joseph_gain(cell.P,H,cell.R);z=EVENTS._state_ad(cell.state)
    residual=_physical_s_residual_ad(z,cell.wave_primitive.centered_S);out=EVENTS._apply_physical_correction(z,K,residual)
    if mode=="A":
        if cell.true_bias is None or cell.bias_projection_limit is None: raise ValueError("A21 S event requires same-source true bias/projection radius")
        J,state_out,projection=EVENTS._compose_A21_bias_projection(out,cell.true_bias,float(cell.bias_projection_limit));branch=projection["branch"]
    else:
        J=EVENTS.AD.jacobian(out);state_out=EVENTS.AD.values(out);branch="not_applicable"
    return {"mode":mode,"kind":"S_zero","H":H,"K":K,"S":S,"residual":EVENTS.AD.values(residual),"state_out":state_out,"J_state":J,
      "physical_centered_S":tuple(cell.wave_primitive.centered_S),"physical_generator_id":cell.wave_primitive.generator_id,
      "primitive_in_id":cell.wave_primitive.primitive_in_id,"primitive_out_id":cell.wave_primitive.primitive_out_id,"live_origin_id":cell.wave_primitive.live_origin_id,
      "actual_RS_provenance":cell.R_provenance,"bias_projection_branch":branch,"same_history_physical_reference_attached":True,
      "independent_S_box_used":False,"wordwise_S_rezero_used":False,"replay_used":False}

def build():
    ref=REF.build();failures=REF.validate(ref)
    if failures: raise RuntimeError("physical reference contract invalid: "+repr(failures))
    return {"qualification":QUALIFICATION,"canonical_source":"COMPLETE_BRMM_NORMAL_LIVE_WORD","finite_S_zero_event_evaluator_available":True,
      "S_zero_residual":"e_S-S_phys","same_estimator_cell_P_R_and_physical_primitive_required":True,"actual_applied_RS_required":True,
      "A21_projection_composed_after_physical_correction":True,"finite_state_and_Jacobian_share_same_AD_composition":True,
      "independent_S_box_used":False,"wordwise_S_rezero_used":False,"replay_used":False,"all_literal_S_events_attached_source_uniformly":False,
      "complete_physical_word_closed":False,"ALT_LIVE_PASS":False,
      "next_obligation":"replace differential-only S events by this finite physical event in every estimator-owned lineage and retain the same primitive/generator ancestry through the literal word"}
def validate(d):
    f=[]
    if d.get("qualification")!=QUALIFICATION:f.append("qualification mismatch")
    for k in ("finite_S_zero_event_evaluator_available","same_estimator_cell_P_R_and_physical_primitive_required","actual_applied_RS_required","A21_projection_composed_after_physical_correction","finite_state_and_Jacobian_share_same_AD_composition"):
        if d.get(k) is not True:f.append(k+" not true")
    for k in ("independent_S_box_used","wordwise_S_rezero_used","replay_used","all_literal_S_events_attached_source_uniformly","complete_physical_word_closed","ALT_LIVE_PASS"):
        if d.get(k) is not False:f.append(k+" not false")
    return f
