#!/usr/bin/env python3
"""Bind the certified same-cell Joseph correction ceiling to a real BRMM event master.

The new estimator-owned event master exposes the literal attitude correction map

    D_theta = E_theta K Q,

where K is derived from that event's SAME P/H/R cell and Q is the structured
physical residual map (including exact chord/cross coordinates for nonlinear
events).  The existing source-uniform correction-domain lemma proves, from the
same-cell covariance identity

    K S K^T = P^- - P^+ <= P^-,

that every admissible event of a given mode/kind satisfies

    ||D_theta z|| <= delta_{mode,kind} h.

This module attaches that event-specific certified delta to the very same
D_theta map and constructs the parameterized finite-reset IQC.  It does not
replace D_theta by a rowwise K norm, does not use the mode-wide worst event when
an event-specific bound is available, and does not promote endpoint/prefix P4.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

import ou3_p4_same_cell_correction_domain as CORR
import ou3_p4_brmm_reduced_event_master as EVENT
import ou3_p4_reset_graph_iqc as RESETIQC
import ou3_p4_hard_entry_set as ENTRY

SCHEMA=1
QUALIFICATION='OU3_P4_BRMM_SAME_DTHETA_RESET_BINDING_V1'


def _mode_name(mode:str)->str:
    if mode=='H':return 'H18'
    if mode=='A':return 'A21'
    raise ValueError('mode must be H or A')


def bind_event(em:EVENT.EventMaster,corr:dict)->dict:
    mode=_mode_name(em.mode)
    if mode not in corr['modes'] or em.kind not in corr['modes'][mode]['events']:
        raise ValueError('correction certificate lacks event mode/kind')
    ce=corr['modes'][mode]['events'][em.kind]
    delta=float(ce['same_cell_attitude_correction_norm_upper'])
    if not(math.isfinite(delta) and delta>=0.0):raise ValueError('invalid certified correction radius')
    if not ce['inside_reset_utility_domain']:
        return {'mode':mode,'kind':em.kind,'delta':delta,'chart_safe':False,'closed':False,
                'reason':'certified correction ceiling exceeds exact-reset utility domain'}
    entry=ENTRY.build();ef=ENTRY.validate(entry)
    if ef:raise RuntimeError('hard-entry prerequisite failed: '+repr(ef))
    q=float(entry['coordinate_radii']['attitude_cayley_norm'])
    Pi,sector=RESETIQC.parameterized_reset_iqc(em.D_theta,em.B_theta,q,delta)
    target=RESETIQC.correction_domain_target(em.D_theta,em.h_index,delta)
    n=em.coordinate_dimension
    dims=(EVENT.shape(Pi)==(n,n) and EVENT.shape(target)==(n,n))
    same_cell_basis=bool(corr['same_cell_Joseph_covariance_identity_consumed']
                         and corr['attitude_block_only_used']
                         and corr['source_uniform_endpoint_covariance_envelope_consumed']
                         and not corr['rowwise_K_bound_used']
                         and not corr['independent_K_box_used'])
    return {'mode':mode,'kind':em.kind,'delta':delta,
      'same_event_Dtheta_map_consumed':True,
      'same_cell_Joseph_metric_bound_consumed':same_cell_basis,
      'event_specific_delta_used_not_mode_worst':True,
      'correction_domain_target':target,'reset_sector':Pi,
      'reset_sector_contract':sector,'chart_safe':bool(sector['chart_safe']),
      'correction_domain_target_certified_by_same_cell_metric_identity':same_cell_basis,
      'rowwise_K_bound_used':False,'independent_K_box_used':False,
      'dimensions_valid':dims,'closed':bool(dims and same_cell_basis and sector['chart_safe'])}


def _smokes():
    im=EVENT._smoke_image();n=18;x=[EVENT.I(0) for _ in range(n)];P=EVENT._identity(n)
    sz=EVENT.COVER.source_cell_from_joint_image(im,mode='H',sample_index=0,event_ordinal=2,kind='S_zero',state=x,P=P,dt_s=EVENT.I(.005),pseudo_elapsed_s=EVENT.I(.01),radial_scale=EVENT.Interval(0,1),event_source_token=im.source_token+':e2',event_predecessor_token=im.source_token+':e1')
    return EVENT.build_event_master(sz),EVENT.build_event_master(EVENT._accel_smoke_cell(im))


def build()->dict:
    corr=CORR.build();cf=CORR.validate(corr)
    if cf:raise RuntimeError('same-cell correction prerequisite failed: '+repr(cf))
    sz,acc=_smokes();bs=bind_event(sz,corr);ba=bind_event(acc,corr)
    closed=bool(bs['closed'] and ba['closed'])
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'same_estimator_owned_event_Dtheta_is_reset_graph_input':True,
      'same_cell_KSK_covariance_identity_is_domain_proof':True,
      'source_uniform_covariance_envelope_consumed':True,
      'actual_applied_RS_lower_consumed_for_S_zero':bool(corr['actual_applied_RS_lower_consumed']),
      'event_specific_delta_bound_attached':True,'mode_worst_delta_substituted':False,
      'rowwise_K_correction_domain_used':False,'independent_K_box_used':False,
      'S_zero_binding':{k:v for k,v in bs.items() if k not in ('correction_domain_target','reset_sector')},
      'accelerometer_binding':{k:v for k,v in ba.items() if k not in ('correction_domain_target','reset_sector')},
      'smoke_event_specific_same_Dtheta_reset_domain_closed':closed,
      'production_all_estimator_owned_events_bound_here':False,
      'production_lineage_wide_reset_domain_closed_here':False,
      'endpoint_augmented_LDLT_closed_here':False,'every_prefix_augmented_LDLT_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'apply bind_event to every estimator-owned Joseph cell in each radial/source lineage, then embed its master+nonlinear+reset sectors into the lineage-wide prefix coordinate with BIAS1 and binary32 ISS maps'}


def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('same_estimator_owned_event_Dtheta_is_reset_graph_input','same_cell_KSK_covariance_identity_is_domain_proof','source_uniform_covariance_envelope_consumed','actual_applied_RS_lower_consumed_for_S_zero','event_specific_delta_bound_attached','smoke_event_specific_same_Dtheta_reset_domain_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('mode_worst_delta_substituted','rowwise_K_correction_domain_used','independent_K_box_used','production_all_estimator_owned_events_bound_here','production_lineage_wide_reset_domain_closed_here','endpoint_augmented_LDLT_closed_here','every_prefix_augmented_LDLT_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    for name in ('S_zero_binding','accelerometer_binding'):
        b=d.get(name,{})
        if b.get('closed') is not True:f.append(name+' not closed')
        if b.get('same_event_Dtheta_map_consumed') is not True:f.append(name+' Dtheta detached')
        if b.get('same_cell_Joseph_metric_bound_consumed') is not True:f.append(name+' same-cell metric not consumed')
        if b.get('event_specific_delta_used_not_mode_worst') is not True:f.append(name+' event delta not specific')
        if b.get('rowwise_K_bound_used') is not False or b.get('independent_K_box_used') is not False:f.append(name+' leaked independent K bound')
        if not(math.isfinite(float(b.get('delta',math.nan))) and float(b['delta'])>=0):f.append(name+' delta invalid')
    return f


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'same_Dtheta_smoke':d['smoke_event_specific_same_Dtheta_reset_domain_closed'],'S_delta':d['S_zero_binding']['delta'],'acc_delta':d['accelerometer_binding']['delta'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
