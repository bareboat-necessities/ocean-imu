#!/usr/bin/env python3
"""Kernel-backed nonlinear-event AD and correction-chart cover.

The trusted complete-window kernel supplies literal same-P/H/R event cells. For
every retained estimator image, H18/A21 mode and exact radial partition cell,
this producer forms estimator-owned event masters and proves the correction
chart directly on the SAME D_theta graph.

* S=0 uses the linear affine hard-ball S-procedure.
* accelerometer/magnetometer use the nonlinear affine S-procedure with every
  exact chord/cross sector retained.

A correction certificate is consumed by the reset sector only when the SAME
D_theta graph proves a chart-safe delta.  Nonlinear events search ascending
strict deltas below the Cayley-composition antipode threshold before considering
wider diagnostic targets.  A failed chart implication is reported as a proof
failure; it is never replaced by delta=3 and never allowed to crash the cover.
No rowwise K bound, covariance-confidence membership, independent residual norm,
or Pbar*residual scalar ceiling is used. This closes the local first-exit
correction implication on the full declared hard domain when every cell passes;
prefix reachability/storage retention remain separate obligations.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
from ou3_interval import Interval
import ou3_p4_brmm_kernel_event_capture_bridge as BRIDGE
import ou3_p4_complete_brmm_source_cover_contract as COVER
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_brmm_reduced_event_master as MASTER
import ou3_p4_brmm_same_Dtheta_reset_binding as RESETBIND
import ou3_p4_brmm_nonlinear_correction_domain as NONLINEAR
import ou3_p4_brmm_szero_correction_domain as SZERO
import ou3_p4_hard_entry_radial_ad_boxes as RADIAL
import ou3_p4_brmm_event_lineage_cover as LINEAGE
import ou3_p4_bias1_family as BIAS1
SCHEMA=6
QUALIFICATION='OU3_P4_BRMM_KERNEL_RADIAL_EVENT_AD_COVER_V6'
JOSEPH=('S_zero','accelerometer','magnetometer')
# Ascending search.  For q=2*tan(45deg/2)=0.8284, deltas through 2.30 rad
# remain on the safe side of the exact Cayley-composition antipode even after
# the deployed correction-quaternion Cayley conversion.  Tighter targets are
# attempted first because they also reduce the reset-sector gain.
NONLINEAR_DELTA_CANDIDATES=(0.25,0.50,0.75,1.00,1.25,1.50,1.75,2.00,2.15,2.25,2.30)
MULTIPLIER_GRID=(1e-3,1e-2,1e-1,1.0,1e1,1e2,1e3)
def I(x):return Interval.point(float(x))
def _limit():
    p=Path(__file__).resolve().parents[2]/'tools'/'stability'/'ou3_proof_operating_domain.json';d=json.loads(p.read_text());return float(d['normal_live']['active_accelerometer_bias_projection_limit_mps2'])
def _true_bias_box():
    b=BIAS1.build();bf=BIAS1.validate(b)
    if bf:raise RuntimeError('BIAS1 prerequisite failed: '+repr(bf))
    r=math.nextafter(float(b['true_bias_norm_upper_mps2']),math.inf);return [Interval(-r,r) for _ in range(3)]
def _captured_by_kind(cells):return {c.kind:c for c in cells if c.kind in JOSEPH}
def _event_token(image,ordinal):return f'{image.source_token}:e{ordinal}'
def _cell(image,mode,kind,state,radial,rc,sample,ordinal,true_bias):
    common=dict(image=image,mode=mode,sample_index=0,event_ordinal=ordinal,kind=kind,state=state,P=rc.P_before,dt_s=I(.005),pseudo_elapsed_s=I(.005*ordinal),radial_scale=radial,event_source_token=_event_token(image,ordinal),event_predecessor_token=image.source_token,true_bias=true_bias,bias_projection_limit=_limit() if mode=='A' else None)
    if kind=='S_zero':return COVER.source_cell_from_joint_image(**common)
    if kind=='accelerometer':return COVER.source_cell_from_joint_image(**common,R=rc.R,f_hat=sample.f_cog_body,R_hat=sample.R_wb)
    m=sample.magnetometer_events_after_imu[int(rc.magnetic_event_index or 0)].m_body
    return COVER.source_cell_from_joint_image(**common,R=rc.R,m_body=m)
def _finite_matrix(A):return all(math.isfinite(x.lo) and math.isfinite(x.hi) for row in A for x in row)
def inspect_event(image,mode,kind,state,radial,rc,sample,ordinal,true_bias):
    cell=_cell(image,mode,kind,state,radial,rc,sample,ordinal,true_bias);cf=COVER.validate_cell(cell,require_estimator_provenance=True)
    if cf:return {'closed':False,'failures':['cell:'+x for x in cf]}
    ev=EVENTS.source_joseph_event(**COVER.joseph_event_kwargs(cell));em=MASTER.build_event_master(cell);ef=MASTER.validate_event_master(cell,em)
    if kind=='S_zero':
        raw=SZERO.certify_szero_event(em)
        chart={'closed':bool(raw['strict_outward_LDLT_closed'] and raw['parameterized_reset_sector_valid_at_certified_delta']),
               'certificate':{'delta':float(raw['certified_delta']),'minimum_pivot_lower':min(map(float,raw['LDLT_pivot_lowers'])) if raw['LDLT_pivot_lowers'] else math.nan,'pivot_lowers':raw['LDLT_pivot_lowers']},'attempts':1,
               'reset_sector_valid':bool(raw['parameterized_reset_sector_valid_at_certified_delta'])}
    else:
        chart=NONLINEAR.certify_event(em,delta_candidates=NONLINEAR_DELTA_CANDIDATES,multiplier_grid=MULTIPLIER_GRID)
    if not chart['closed']:
        state_finite=all(math.isfinite(x.lo) and math.isfinite(x.hi) for x in ev['state_out']);jac_finite=_finite_matrix(ev['J_state']);master_finite=_finite_matrix(em.master)
        return {'closed':False,'failures':ef+['same-graph correction/reset chart did not close'],'same_P_H_R_cell':ev['same_P_H_R_cell'],'R_provenance':ev['R_provenance'],'state_out_finite':state_finite,'Jacobian_finite':jac_finite,'real_event_master_finite':master_finite,'first_exit_reset_structure_closed':False,'certified_correction_delta':float(chart['certificate']['delta']) if chart.get('certificate') else math.nan,'same_graph_correction_target_closed':False,'same_graph_chart_certificate':chart.get('certificate'),'same_graph_chart_attempts':int(chart['attempts']),'correction_domain_target_preproved':False,'projection_branch':ev['bias_projection_branch'],'radial':radial.as_list(),'event_token':cell.source_token}
    delta=float(chart['certificate']['delta'])
    try:
        rb=RESETBIND.bind_event_first_exit(em,delta)
    except (ValueError,RuntimeError) as exc:
        return {'closed':False,'failures':ef+['reset binding:'+str(exc)],'same_P_H_R_cell':ev['same_P_H_R_cell'],'R_provenance':ev['R_provenance'],'state_out_finite':True,'Jacobian_finite':_finite_matrix(ev['J_state']),'real_event_master_finite':_finite_matrix(em.master),'first_exit_reset_structure_closed':False,'certified_correction_delta':delta,'same_graph_correction_target_closed':True,'same_graph_chart_certificate':chart.get('certificate'),'same_graph_chart_attempts':int(chart['attempts']),'correction_domain_target_preproved':True,'projection_branch':ev['bias_projection_branch'],'radial':radial.as_list(),'event_token':cell.source_token}
    state_finite=all(math.isfinite(x.lo) and math.isfinite(x.hi) for x in ev['state_out']);jac_finite=_finite_matrix(ev['J_state']);master_finite=_finite_matrix(em.master);target_finite=_finite_matrix(rb['correction_domain_target']);sector_finite=_finite_matrix(rb['reset_sector'])
    structural=not ef and rb.get('closed') and state_finite and jac_finite and master_finite and target_finite and sector_finite
    return {'closed':bool(structural),'failures':ef,'same_P_H_R_cell':ev['same_P_H_R_cell'],'R_provenance':ev['R_provenance'],'state_out_finite':state_finite,'Jacobian_finite':jac_finite,'real_event_master_finite':master_finite,'first_exit_reset_structure_closed':bool(rb.get('closed')),'certified_correction_delta':delta,'same_graph_correction_target_closed':True,'same_graph_chart_certificate':chart.get('certificate'),'same_graph_chart_attempts':int(chart['attempts']),'correction_domain_target_preproved':True,'projection_branch':ev['bias_projection_branch'],'radial':radial.as_list(),'event_token':cell.source_token}
def build():
    sample,_,meta,_,images=BRIDGE.capture_smoke();parts=RADIAL.partition(2);tb=_true_bias_box();records=[]
    for image in images:
      rec={'image_token':image.source_token}
      for mode,key in (('H','H_event_cells'),('A','A_event_cells')):
        by=_captured_by_kind(meta[key]);rows=[]
        for p in parts:
          radial=p.interval;state=RADIAL.radial_state_box(mode,radial)
          for ordinal,kind in enumerate(JOSEPH,2):rows.append({'kind':kind,**inspect_event(image,mode,kind,state,radial,by[kind],sample,ordinal,tb if mode=='A' else None)})
        rec['H18' if mode=='H' else 'A21']=rows
      records.append(rec)
    all_rows=[r for rec in records for mode in ('H18','A21') for r in rec[mode]]
    all_closed=bool(all_rows) and all(r['closed'] for r in all_rows);all_same=all(r.get('same_P_H_R_cell') for r in all_rows);all_finite=all(r.get('state_out_finite') and r.get('Jacobian_finite') and r.get('real_event_master_finite') for r in all_rows);owned=all(str(r.get('event_token','')).startswith(rec['image_token']+':e') for rec in records for mode in ('H18','A21') for r in rec[mode]);all_chart=bool(all_rows) and all(r.get('same_graph_correction_target_closed') is True and r.get('first_exit_reset_structure_closed') is True for r in all_rows);a21_branches=sorted(set(r['projection_branch'] for rec in records for r in rec['A21']))
    deltas=[float(r['certified_correction_delta']) for r in all_rows if math.isfinite(float(r.get('certified_correction_delta',math.nan)))];pivots=[float(r['same_graph_chart_certificate']['minimum_pivot_lower']) for r in all_rows if r.get('same_graph_chart_certificate') and math.isfinite(float(r['same_graph_chart_certificate'].get('minimum_pivot_lower',math.nan)))]
    failures=[{'mode':mode,'kind':r['kind'],'radial':r.get('radial'),'delta':r.get('certified_correction_delta'),'failures':r.get('failures',[])} for rec in records for mode in ('H18','A21') for r in rec[mode] if not r['closed']]
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','trusted_kernel_passive_event_cells_consumed':True,'every_same_sample_joint_estimator_image_consumed':True,'current_lineage_radial_API_consumed':True,'event_tokens_remain_estimator_owned':owned,'continuous_radial_partition_exactly_covers_unit_interval':LINEAGE.partition_covers_unit(parts),'radial_partition_bounds':[p.interval.as_list() for p in parts],'full_hard_entry_radial_AD_outer_boxes_consumed':True,'hard_ball_IQCs_replaced_by_boxes':False,'same_cell_P_H_R_K_used_for_every_Joseph_AD_cell':all_same,'all_radial_Joseph_state_outputs_and_Jacobians_finite':all_finite,'all_radial_same_graph_correction_targets_closed':all_chart,'worst_certified_correction_delta':max(deltas,default=math.nan),'minimum_correction_LDLT_pivot_lower':min(pivots,default=math.nan),'all_real_structured_event_masters_and_first_exit_reset_bindings_constructed':all_closed,'all_real_structured_event_masters_and_reset_bindings_constructed':all_closed,'A21_projection_branches_observed':a21_branches,'A21_projection_not_assumed_inactive':True,'BIAS1_true_bias_outer_box_consumed_for_A21_AD':True,'nonlinear_delta_candidates':list(NONLINEAR_DELTA_CANDIDATES),'records':records,'failed_cells':failures,'covariance_membership_used_as_state_hypothesis':False,'rowwise_K_bound_used':False,'independent_residual_norm_used':False,'trajectory_replay_used':False,'local_full_hard_box_at_each_event_claimed_reachable':False,'production_prefix_state_reachability_closed_here':False,'production_every_prefix_hard_domain_retention_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,'next_obligation':'if every correction/reset cell closes, propagate this correction-safe full hard-domain family through source-connected literal prefixes and prove the common endpoint/every-prefix storage plus first-exit hard-domain retention'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('trusted_kernel_passive_event_cells_consumed','every_same_sample_joint_estimator_image_consumed','current_lineage_radial_API_consumed','event_tokens_remain_estimator_owned','continuous_radial_partition_exactly_covers_unit_interval','full_hard_entry_radial_AD_outer_boxes_consumed','same_cell_P_H_R_K_used_for_every_Joseph_AD_cell','all_radial_Joseph_state_outputs_and_Jacobians_finite'):
        if d.get(k) is not True:f.append(k+' not true')
    # Correction/reset closure is now a theorem result, not a schema invariant:
    # keep the artifact available even when a cell fails so CI exposes the exact
    # mode/kind/radial obstruction instead of crashing before evidence upload.
    for k in ('hard_ball_IQCs_replaced_by_boxes','covariance_membership_used_as_state_hypothesis','rowwise_K_bound_used','independent_residual_norm_used','trajectory_replay_used','local_full_hard_box_at_each_event_claimed_reachable','production_prefix_state_reachability_closed_here','production_every_prefix_hard_domain_retention_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if not d.get('records'):f.append('empty radial event records')
    if d.get('all_radial_same_graph_correction_targets_closed'):
        if not(0<float(d.get('worst_certified_correction_delta',0))<3):f.append('correction delta not strict')
        if not float(d.get('minimum_correction_LDLT_pivot_lower',0))>0:f.append('correction LDLT pivot not strict')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'radial_cells':len(d['radial_partition_bounds']),'projection_branches':d['A21_projection_branches_observed'],'all_chart':d['all_radial_same_graph_correction_targets_closed'],'worst_delta':d['worst_certified_correction_delta'],'min_pivot':d['minimum_correction_LDLT_pivot_lower'],'all_masters':d['all_real_structured_event_masters_and_first_exit_reset_bindings_constructed'],'failed_cells':d['failed_cells'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
