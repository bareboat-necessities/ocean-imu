#!/usr/bin/env python3
"""Kernel-backed nonlinear-event AD cover over the full hard-entry radial family.

The trusted complete-window kernel passively captures the literal shipping
P/H/R cells.  The same-sample estimator bridge supplies every retained joint
estimator image.  For each image, H18/A21 mode, closed hard-entry radial cell,
and Joseph kind (S_zero, accelerometer, magnetometer), this module constructs an
estimator-owned SourceCoverCell whose physical-error state is the outward AD box
for that radial cell and evaluates the real same-cell nonlinear Joseph map.

This proves a useful local fact that the zero-root smoke did not: the complete
45-degree/full-hard-entry finite state family can be fed through the exact
Joseph/projection AD machinery with the captured shipping P/H/R/K objects and
actual applied S=0 R_S.  A21 uses an outward true-bias box derived from the
admitted BIAS1 norm and therefore retains active/inactive/Clarke projection
branches when the box crosses the projection sphere.

The radial state box is only an outer enclosure for AD.  This module does NOT
claim that the full hard-entry box is reachable at every later event and does
not replace the original hard-ball IQCs.  Prefix reachability/first-exit
retention remains a separate mandatory step.
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
import ou3_p4_hard_entry_radial_ad_boxes as RADIAL
import ou3_p4_brmm_event_lineage_cover as LINEAGE
import ou3_p4_bias1_family as BIAS1
import ou3_p4_same_cell_correction_domain as CORR

SCHEMA=1
QUALIFICATION='OU3_P4_BRMM_KERNEL_RADIAL_EVENT_AD_COVER_V1'
JOSEPH=('S_zero','accelerometer','magnetometer')

def I(x):return Interval.point(float(x))
def _limit():
    p=Path(__file__).resolve().parents[2]/'tools'/'stability'/'ou3_proof_operating_domain.json'
    d=json.loads(p.read_text());return float(d['normal_live']['active_accelerometer_bias_projection_limit_mps2'])
def _true_bias_box():
    b=BIAS1.build();bf=BIAS1.validate(b)
    if bf:raise RuntimeError('BIAS1 prerequisite failed: '+repr(bf))
    # Norm-ball -> componentwise outer box, used only for AD enclosure.
    r=math.nextafter(float(b['true_bias_norm_upper_mps2']),math.inf)
    return [Interval(-r,r) for _ in range(3)]
def _captured_by_kind(cells):return {c.kind:c for c in cells if c.kind in JOSEPH}
def _cell(image,mode,kind,state,radial,rc,sample,ordinal,true_bias):
    common=dict(image=image,mode=mode,sample_index=0,event_ordinal=ordinal,kind=kind,state=state,P=rc.P_before,
                dt_s=I(.005),pseudo_elapsed_s=I(.005*ordinal),radial_scale=radial,
                event_source_token=f'{image.source_token}:{mode}:{kind}:r{radial.lo:.6g}-{radial.hi:.6g}',
                event_predecessor_token=image.predecessor_token,true_bias=true_bias,
                bias_projection_limit=_limit() if mode=='A' else None)
    if kind=='S_zero':return COVER.source_cell_from_joint_image(**common)
    if kind=='accelerometer':return COVER.source_cell_from_joint_image(**common,R=rc.R,f_hat=sample.f_cog_body,R_hat=sample.R_wb)
    m=sample.magnetometer_events_after_imu[int(rc.magnetic_event_index or 0)].m_body
    return COVER.source_cell_from_joint_image(**common,R=rc.R,m_body=m)
def _finite_matrix(A):return all(math.isfinite(x.lo) and math.isfinite(x.hi) for row in A for x in row)
def inspect_event(image,mode,kind,state,radial,rc,sample,ordinal,true_bias,corr):
    cell=_cell(image,mode,kind,state,radial,rc,sample,ordinal,true_bias)
    cf=COVER.validate_cell(cell,require_estimator_provenance=True)
    if cf:return {'closed':False,'failures':['cell:'+x for x in cf]}
    ev=EVENTS.source_joseph_event(**COVER.joseph_event_kwargs(cell))
    em=MASTER.build_event_master(cell);ef=MASTER.validate_event_master(cell,em)
    rb=RESETBIND.bind_event(em,corr)
    state_finite=all(math.isfinite(x.lo) and math.isfinite(x.hi) for x in ev['state_out'])
    jac_finite=_finite_matrix(ev['J_state']);master_finite=_finite_matrix(em.master)
    projection=ev['bias_projection_branch']
    return {'closed':not ef and rb.get('closed') and state_finite and jac_finite and master_finite,
      'failures':ef,'same_P_H_R_cell':ev['same_P_H_R_cell'],'R_provenance':ev['R_provenance'],
      'state_out_finite':state_finite,'Jacobian_finite':jac_finite,'real_event_master_finite':master_finite,
      'reset_binding_closed':bool(rb.get('closed')),'reset_delta':float(rb.get('delta',math.nan)),
      'projection_branch':projection,'radial':radial.as_list()}
def build():
    sample,_,meta,_,images=BRIDGE.capture_smoke();corr=CORR.build();cf=CORR.validate(corr)
    if cf:raise RuntimeError('correction prerequisite failed: '+repr(cf))
    root=LINEAGE.root_radial_cell();parts=LINEAGE.partition(root,4);tb=_true_bias_box();records=[]
    for image in images:
      rec={'image_token':image.source_token}
      for mode,key in (('H','H_event_cells'),('A','A_event_cells')):
        by=_captured_by_kind(meta[key]);rows=[]
        for p in parts:
          state=RADIAL.radial_state_box(mode,p.radial_scale)
          for ordinal,kind in enumerate(JOSEPH,2):
            rows.append({'kind':kind,**inspect_event(image,mode,kind,state,p.radial_scale,by[kind],sample,ordinal,tb if mode=='A' else None,corr)})
        rec['H18' if mode=='H' else 'A21']=rows
      records.append(rec)
    all_rows=[r for rec in records for mode in ('H18','A21') for r in rec[mode]]
    all_closed=bool(all_rows) and all(r['closed'] for r in all_rows)
    all_same=all(r.get('same_P_H_R_cell') for r in all_rows)
    all_finite=all(r.get('state_out_finite') and r.get('Jacobian_finite') and r.get('real_event_master_finite') for r in all_rows)
    a21_branches=sorted(set(r['projection_branch'] for rec in records for r in rec['A21']))
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'trusted_kernel_passive_event_cells_consumed':True,'every_same_sample_joint_estimator_image_consumed':True,
      'continuous_radial_partition_exactly_covers_unit_interval':LINEAGE.exact_partition_cover(parts),
      'radial_partition_bounds':[p.radial_scale.as_list() for p in parts],
      'full_hard_entry_radial_AD_outer_boxes_consumed':True,'hard_ball_IQCs_replaced_by_boxes':False,
      'same_cell_P_H_R_K_used_for_every_Joseph_AD_cell':all_same,
      'all_radial_Joseph_state_outputs_and_Jacobians_finite':all_finite,
      'all_real_structured_event_masters_and_reset_bindings_constructed':all_closed,
      'A21_projection_branches_observed':a21_branches,'A21_projection_not_assumed_inactive':True,
      'BIAS1_true_bias_outer_box_consumed_for_A21_AD':True,'records':records,
      'covariance_membership_used_as_state_hypothesis':False,'trajectory_replay_used':False,
      'local_full_hard_box_at_each_event_claimed_reachable':False,
      'production_prefix_state_reachability_closed_here':False,'production_every_prefix_hard_domain_retention_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'propagate the hard-entry radial family through captured prediction/Joseph operations with same-history prefix transport and prove first-exit retention; then replace the conservative local full-hard-box-at-each-event AD cover by reachable prefix boxes'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('trusted_kernel_passive_event_cells_consumed','every_same_sample_joint_estimator_image_consumed','continuous_radial_partition_exactly_covers_unit_interval','full_hard_entry_radial_AD_outer_boxes_consumed','same_cell_P_H_R_K_used_for_every_Joseph_AD_cell','all_radial_Joseph_state_outputs_and_Jacobians_finite','all_real_structured_event_masters_and_reset_bindings_constructed','A21_projection_not_assumed_inactive','BIAS1_true_bias_outer_box_consumed_for_A21_AD'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('hard_ball_IQCs_replaced_by_boxes','covariance_membership_used_as_state_hypothesis','trajectory_replay_used','local_full_hard_box_at_each_event_claimed_reachable','production_prefix_state_reachability_closed_here','production_every_prefix_hard_domain_retention_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if not d.get('records'):f.append('empty radial event records')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'radial_cells':len(d['radial_partition_bounds']),'projection_branches':d['A21_projection_branches_observed'],'all_event_AD':d['all_radial_Joseph_state_outputs_and_Jacobians_finite'],'all_masters':d['all_real_structured_event_masters_and_reset_bindings_constructed'],'reachability':d['production_prefix_state_reachability_closed_here'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
