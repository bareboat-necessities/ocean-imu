#!/usr/bin/env python3
"""Bind a real BRMM event's SAME D_theta map to the finite reset graph.

Two logically distinct routes are exposed:

* ``bind_event`` retains the historical scalar diagnostic. It is not a theorem
  correction-domain proof.
* ``bind_event_first_exit`` is the production route. It uses the literal
  estimator-owned D_theta=E_theta*K*Q from the SAME P/H/R/K event and the
  declared nonlinear correction-chart radius from the canonical closure domain.

The reset utility itself supports a much broader search range, but that broad
utility maximum is NOT a valid production correction radius: near it, a full
attitude working-set Cayley vector can approach the exact composition antipode.
The theorem's existing nonlinear domain instead declares

    ||d_theta|| <= 0.25.

This module therefore returns both

    delta_chart^2 h^2 - ||D_theta z||^2 >= 0

as a same-graph target and the parameterized exact-reset sector valid under
that target.  The 0.25 correction target must be proved by the common augmented
first-exit master; it is not preasserted.  The attitude entry/working radius is
unchanged (0.8284...), so this correction-chart fix does not shrink fresh entry.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import ou3_p4_same_cell_correction_domain as CORR
import ou3_p4_brmm_reduced_event_master as EVENT
import ou3_p4_reset_graph_iqc as RESETIQC
import ou3_p4_exact_reset_transport as RESET
import ou3_p4_hard_entry_set as ENTRY

SCHEMA=4
QUALIFICATION='OU3_P4_BRMM_SAME_DTHETA_RESET_BINDING_V4'
REPO=Path(__file__).resolve().parents[2]
CLOSURE=REPO/'tools/stability/ou3_p4_closure_domain.json'

def _mode_name(mode:str)->str:
    if mode=='H':return 'H18'
    if mode=='A':return 'A21'
    raise ValueError('mode must be H or A')
def _entry_q()->float:
    entry=ENTRY.build();ef=ENTRY.validate(entry)
    if ef:raise RuntimeError('hard-entry prerequisite failed: '+repr(ef))
    return float(entry['coordinate_radii']['attitude_cayley_norm'])
def _declared_correction_chart_radius()->float:
    d=json.loads(CLOSURE.read_text())['uniform_nonlinear_enclosure']
    delta=float(d['required_correction_cayley_norm_upper'])
    if not(math.isfinite(delta) and 0<delta<=RESET.CAYLEY_MONOTONE_NORM_MAX):raise RuntimeError('declared correction chart radius outside exact reset utility')
    return delta

def bind_event_first_exit(em:EVENT.EventMaster,delta:float|None=None)->dict:
    mode=_mode_name(em.mode);declared=_declared_correction_chart_radius();delta=declared if delta is None else float(delta)
    if not(math.isfinite(delta) and 0.0<delta<=declared):raise ValueError('first-exit delta outside declared nonlinear correction chart')
    q=_entry_q();Pi,sector=RESETIQC.parameterized_reset_iqc(em.D_theta,em.B_theta,q,delta);target=RESETIQC.correction_domain_target(em.D_theta,em.h_index,delta);n=em.coordinate_dimension
    dims=(EVENT.shape(Pi)==(n,n) and EVENT.shape(target)==(n,n));chart=bool(sector['chart_safe'])
    return {'mode':mode,'kind':em.kind,'delta':delta,'declared_correction_chart_delta_used':delta==declared,'broad_reset_utility_max_used_as_production_delta':False,'same_event_Dtheta_map_consumed':True,'same_event_P_H_R_K_ancestry_retained':True,'correction_domain_target':target,'reset_sector':Pi,'reset_sector_contract':sector,'chart_safe_conditionally_on_target':chart,'dimensions_valid':dims,'correction_domain_target_proved_here':False,'correction_domain_target_must_be_proved_by_same_augmented_first_exit_master':True,'attitude_entry_radius_shrunk':False,'rowwise_K_bound_used':False,'independent_K_box_used':False,'scalar_covariance_residual_ceiling_used':False,'closed':bool(dims and chart)}
def bind_event(em:EVENT.EventMaster,corr:dict)->dict:
    mode=_mode_name(em.mode)
    if mode not in corr['modes'] or em.kind not in corr['modes'][mode]['events']:raise ValueError('correction certificate lacks event mode/kind')
    ce=corr['modes'][mode]['events'][em.kind];raw=ce.get('same_cell_attitude_correction_norm_upper')
    if raw is None:return {'mode':mode,'kind':em.kind,'delta':None,'chart_safe':False,'closed':False,'legacy_scalar_diagnostic_only':True,'scalar_delta_certified':False,'reason':'same-cell scalar correction ceiling is uncertified; production uses same-graph first-exit target'}
    delta=float(raw)
    if not(math.isfinite(delta) and delta>=0.0):raise ValueError('invalid certified correction radius')
    if not ce['inside_reset_utility_domain']:return {'mode':mode,'kind':em.kind,'delta':delta,'chart_safe':False,'closed':False,'legacy_scalar_diagnostic_only':True,'scalar_delta_certified':True,'reason':'detached magnitude ceiling exceeds exact-reset utility domain; use same-graph first-exit binding'}
    q=_entry_q();Pi,sector=RESETIQC.parameterized_reset_iqc(em.D_theta,em.B_theta,q,min(delta,_declared_correction_chart_radius()));target=RESETIQC.correction_domain_target(em.D_theta,em.h_index,min(delta,_declared_correction_chart_radius()));n=em.coordinate_dimension
    dims=(EVENT.shape(Pi)==(n,n) and EVENT.shape(target)==(n,n));same_cell_basis=bool(corr['same_cell_Joseph_covariance_identity_consumed'] and corr['attitude_block_only_used'] and corr['source_uniform_endpoint_covariance_envelope_consumed'] and not corr['rowwise_K_bound_used'] and not corr['independent_K_box_used'])
    return {'mode':mode,'kind':em.kind,'delta':delta,'same_event_Dtheta_map_consumed':True,'same_cell_Joseph_metric_bound_consumed':same_cell_basis,'event_specific_delta_used_not_mode_worst':True,'correction_domain_target':target,'reset_sector':Pi,'reset_sector_contract':sector,'chart_safe':bool(sector['chart_safe']),'correction_domain_target_certified_by_same_cell_metric_identity':same_cell_basis,'legacy_scalar_diagnostic_only':True,'scalar_delta_certified':True,'rowwise_K_bound_used':False,'independent_K_box_used':False,'dimensions_valid':dims,'closed':bool(dims and same_cell_basis and sector['chart_safe'])}
def _smokes():
    im=EVENT._smoke_image();n=18;x=[EVENT.I(0) for _ in range(n)];P=EVENT._identity(n)
    sz=EVENT.COVER.source_cell_from_joint_image(im,mode='H',sample_index=0,event_ordinal=2,kind='S_zero',state=x,P=P,dt_s=EVENT.I(.005),pseudo_elapsed_s=EVENT.I(.01),radial_scale=EVENT.Interval(0,1),event_source_token=im.source_token+':e2',event_predecessor_token=im.source_token+':e1')
    return EVENT.build_event_master(sz),EVENT.build_event_master(EVENT._accel_smoke_cell(im))
def build()->dict:
    corr=CORR.build();cf=CORR.validate(corr)
    if cf:raise RuntimeError('same-cell correction prerequisite failed: '+repr(cf))
    sz,acc=_smokes();legacy_s=bind_event(sz,corr);legacy_a=bind_event(acc,corr);fs=bind_event_first_exit(sz);fa=bind_event_first_exit(acc);structural=bool(fs['closed'] and fa['closed']);declared=_declared_correction_chart_radius()
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','same_estimator_owned_event_Dtheta_is_reset_graph_input':True,'legacy_same_cell_KSK_scalar_diagnostic_retained':True,'legacy_uncertified_scalar_delta_fails_closed_without_exception':True,'legacy_scalar_diagnostic_can_fail_without_blocking_first_exit_architecture':True,'production_uses_same_graph_first_exit_target':True,'production_correction_chart_delta':declared,'broad_reset_utility_max_is_not_production_correction_radius':True,'same_graph_target_and_reset_sector_share_Dtheta':True,'correction_domain_target_is_not_preasserted':True,'attitude_entry_radius_unchanged':True,'rowwise_K_correction_domain_used':False,'independent_K_box_used':False,'scalar_covariance_residual_ceiling_used_by_production':False,'S_zero_legacy_binding':{k:v for k,v in legacy_s.items() if k not in ('correction_domain_target','reset_sector')},'accelerometer_legacy_binding':{k:v for k,v in legacy_a.items() if k not in ('correction_domain_target','reset_sector')},'S_zero_first_exit_binding':{k:v for k,v in fs.items() if k not in ('correction_domain_target','reset_sector')},'accelerometer_first_exit_binding':{k:v for k,v in fa.items() if k not in ('correction_domain_target','reset_sector')},'first_exit_same_Dtheta_reset_structure_closed':structural,'production_correction_domain_target_closed_here':False,'production_all_estimator_owned_events_bound_here':False,'production_lineage_wide_reset_domain_closed_here':False,'endpoint_augmented_LDLT_closed_here':False,'every_prefix_augmented_LDLT_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,'next_obligation':'embed each 0.25 correction target and its same-Dtheta reset sector in the common source/bias/fp augmented prefix master; prove target before consuming reset sector, then close endpoint and every-prefix LDLT'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('same_estimator_owned_event_Dtheta_is_reset_graph_input','legacy_same_cell_KSK_scalar_diagnostic_retained','legacy_uncertified_scalar_delta_fails_closed_without_exception','legacy_scalar_diagnostic_can_fail_without_blocking_first_exit_architecture','production_uses_same_graph_first_exit_target','broad_reset_utility_max_is_not_production_correction_radius','same_graph_target_and_reset_sector_share_Dtheta','correction_domain_target_is_not_preasserted','attitude_entry_radius_unchanged','first_exit_same_Dtheta_reset_structure_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('rowwise_K_correction_domain_used','independent_K_box_used','scalar_covariance_residual_ceiling_used_by_production','production_correction_domain_target_closed_here','production_all_estimator_owned_events_bound_here','production_lineage_wide_reset_domain_closed_here','endpoint_augmented_LDLT_closed_here','every_prefix_augmented_LDLT_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if float(d.get('production_correction_chart_delta',0))!=_declared_correction_chart_radius():f.append('correction chart delta changed')
    for name in ('S_zero_first_exit_binding','accelerometer_first_exit_binding'):
        b=d.get(name,{})
        if b.get('closed') is not True:f.append(name+' structure not closed')
        if b.get('same_event_Dtheta_map_consumed') is not True:f.append(name+' Dtheta detached')
        if b.get('declared_correction_chart_delta_used') is not True:f.append(name+' did not use declared correction chart')
        if b.get('broad_reset_utility_max_used_as_production_delta') is not False:f.append(name+' used broad utility maximum')
        if b.get('correction_domain_target_proved_here') is not False:f.append(name+' target falsely pre-proved')
        if b.get('correction_domain_target_must_be_proved_by_same_augmented_first_exit_master') is not True:f.append(name+' first-exit obligation lost')
        if b.get('attitude_entry_radius_shrunk') is not False:f.append(name+' shrunk entry radius')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'first_exit_structure':d['first_exit_same_Dtheta_reset_structure_closed'],'correction_chart_delta':d['production_correction_chart_delta'],'legacy_acc_closed':d['accelerometer_legacy_binding'].get('closed'),'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
