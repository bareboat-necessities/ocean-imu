#!/usr/bin/env python3
"""Kernel-backed radial first-exit Joseph working cells for P4.

For finite-state first-exit analysis, every prefix strictly before the first
exit lies in the declared hard working domain.  Therefore a Joseph event that
could cause the first exit may be differentiated/enclosed on an outward box of
that same working domain.  This module combines:

* trusted same-transition kernel P/H/R captures;
* every same-sample joint-estimator Image (no favorable successor selection);
* a continuous exact radial partition of [0,1];
* rectangular AD boxes that only OUTER-enclose the hard-ball geometry; and
* the original hard-ball/radial coordinates retained separately for later IQCs.

For each S=0, accelerometer and magnetometer event it emits an estimator-owned
SourceCoverCell and builds the real same-cell nonlinear event plus structured
reduced event master.  A21 true physical bias is enclosed from the qualified
BIAS1 family, independently of covariance.

These are theorem working-domain cells, not a reachability claim.  Promotion
still requires the every-prefix hard-domain/first-exit implication and complete
source continuation over the admitted BRMM family.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_brmm_kernel_event_capture_bridge as BRIDGE
import ou3_p4_hard_entry_radial_ad_boxes as RADIAL
import ou3_p4_brmm_event_lineage_cover as LINEAGE
import ou3_p4_complete_brmm_source_cover_contract as COVER
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_brmm_reduced_event_master as MASTER
import ou3_p4_bias1_family as BIAS1

SCHEMA=1
QUALIFICATION='OU3_P4_BRMM_KERNEL_RADIAL_FIRST_EXIT_WORKING_CELLS_V1'

def I(x):return Interval.point(float(x))
def _limit():
    p=Path(__file__).resolve().parents[2]/'tools'/'stability'/'ou3_proof_operating_domain.json'
    return float(json.loads(p.read_text())['normal_live']['active_accelerometer_bias_projection_limit_mps2'])
def _true_bias_box():
    b=BIAS1.build();bf=BIAS1.validate(b)
    if bf:raise RuntimeError('BIAS1 prerequisite failed: '+repr(bf))
    r=float(b['true_bias_norm_upper_mps2']);u=math.nextafter(r,math.inf)
    return [Interval(-u,u) for _ in range(3)]
def _cell(image,rc,sample,mode,radial,ordinal,state,pred):
    token=f'{image.source_token}:r{radial.token}:e{ordinal}'
    common=dict(image=image,mode=mode,sample_index=0,event_ordinal=ordinal,kind=rc.kind,state=state,P=rc.P_before,
      dt_s=I(.005),pseudo_elapsed_s=I(.005*ordinal),radial_scale=radial.radial_scale,event_source_token=token,event_predecessor_token=pred,
      true_bias=_true_bias_box() if mode=='A' else None,bias_projection_limit=_limit() if mode=='A' else None)
    if rc.kind=='S_zero':return COVER.source_cell_from_joint_image(**common)
    if rc.kind=='accelerometer':return COVER.source_cell_from_joint_image(**common,R=rc.R,f_hat=sample.f_cog_body,R_hat=sample.R_wb)
    if rc.kind=='magnetometer':
        m=sample.magnetometer_events_after_imu[int(rc.magnetic_event_index or 0)].m_body
        return COVER.source_cell_from_joint_image(**common,R=rc.R,m_body=m)
    raise ValueError('Joseph working cells only')
def build_lineage(image,captured,sample,mode,radial):
    state=RADIAL.radial_state_box(mode,radial.radial_scale);pred=image.predecessor_token;rows=[];fail=[]
    for ordinal,rc in enumerate(captured):
        if rc.kind not in ('S_zero','accelerometer','magnetometer'):continue
        cell=_cell(image,rc,sample,mode,radial,ordinal,state,pred)
        vf=COVER.validate_cell(cell,require_estimator_provenance=True);fail.extend(f'{rc.kind}:cell:{x}' for x in vf)
        ev=EVENTS.source_joseph_event(**COVER.joseph_event_kwargs(cell))
        em=MASTER.build_event_master(cell);mf=MASTER.validate_event_master(cell,em);fail.extend(f'{rc.kind}:master:{x}' for x in mf)
        rows.append({'kind':rc.kind,'event_token':cell.source_token,'event_ordinal':ordinal,'event_jacobian_dimension':len(ev['J_state']),
          'projection_branch':ev['bias_projection_branch'],'master_dimension':em.coordinate_dimension,
          'same_cell_P_H_R':bool(ev['same_P_H_R_cell']),'actual_RS':cell.R_provenance==EVENTS.ACTUAL_RS_PROVENANCE if rc.kind=='S_zero' else None,
          'nonlinear_sector_count':len(em.nonlinear_sectors)})
        pred=cell.source_token
    return {'image_token':image.source_token,'mode':mode,'radial_token':radial.token,'radial_scale':radial.radial_scale.as_list(),
      'Joseph_events':rows,'event_count':len(rows),'all_cells_and_masters_valid':not fail,'failures':fail}
def build():
    sample,_,meta,_,images=BRIDGE.capture_smoke();root=LINEAGE.root_radial_cell();radials=LINEAGE.partition(root,2);records=[]
    for radial in radials:
      for image in images:
        records.append(build_lineage(image,meta['H_event_cells'],sample,'H',radial))
        records.append(build_lineage(image,meta['A_event_cells'],sample,'A',radial))
    valid=bool(records) and all(r['all_cells_and_masters_valid'] and r['event_count']==3 for r in records)
    total=sum(r['event_count'] for r in records)
    s_ok=all(next(x for x in r['Joseph_events'] if x['kind']=='S_zero')['actual_RS'] is True for r in records)
    nonlinear=all(next(x for x in r['Joseph_events'] if x['kind']=='accelerometer')['nonlinear_sector_count']>0 and next(x for x in r['Joseph_events'] if x['kind']=='magnetometer')['nonlinear_sector_count']>0 for r in records)
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'trusted_kernel_same_transition_captures_consumed':True,'every_same_sample_joint_estimator_image_consumed':True,
      'continuous_radial_partition_consumed':True,'radial_partition_exactly_covers_unit_interval':LINEAGE.exact_partition_cover(radials),
      'hard_ball_geometry_replaced_by_AD_box':False,'AD_box_is_outer_working_domain_enclosure_only':True,
      'first_exit_semantics_required_for_working_cells':True,'working_cell_reachability_claimed':False,
      'qualified_BIAS1_true_bias_box_attached_to_A21':True,'covariance_membership_used_for_error_state':False,
      'record_count':len(records),'Joseph_event_count_total':total,'records':records,
      'all_working_cells_and_real_event_masters_valid':valid,'all_Szero_cells_use_actual_applied_RS':s_ok,
      'all_nonlinear_cells_retain_structured_graph_sectors':nonlinear,
      'production_every_prefix_first_exit_retention_closed_here':False,'production_complete_BRMM_source_family_covered_here':False,
      'production_reachable_state_lineage_claimed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'use these first-exit working-cell Jacobians/structured sectors in the every-prefix hard-domain S-procedure; prove no prefix exits the declared domain, then the conditional working cells become a valid reachable-event cover along each retained source/radial lineage'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('trusted_kernel_same_transition_captures_consumed','every_same_sample_joint_estimator_image_consumed','continuous_radial_partition_consumed','radial_partition_exactly_covers_unit_interval','AD_box_is_outer_working_domain_enclosure_only','first_exit_semantics_required_for_working_cells','qualified_BIAS1_true_bias_box_attached_to_A21','all_working_cells_and_real_event_masters_valid','all_Szero_cells_use_actual_applied_RS','all_nonlinear_cells_retain_structured_graph_sectors'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('hard_ball_geometry_replaced_by_AD_box','working_cell_reachability_claimed','covariance_membership_used_for_error_state','production_every_prefix_first_exit_retention_closed_here','production_complete_BRMM_source_family_covered_here','production_reachable_state_lineage_claimed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if int(d.get('record_count',0))<=0 or int(d.get('Joseph_event_count_total',0))!=3*int(d.get('record_count',0)):f.append('working event count mismatch')
    for r in d.get('records',[]):
        if r.get('all_cells_and_masters_valid') is not True or r.get('failures')!=[]:f.append('working record invalid '+str(r.get('radial_token')))
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'records':d['record_count'],'Joseph_events':d['Joseph_event_count_total'],'valid':d['all_working_cells_and_real_event_masters_valid'],'first_exit_closed':d['production_every_prefix_first_exit_retention_closed_here'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
