#!/usr/bin/env python3
"""All-bias radial first-exit Joseph working cells for OU-III P4.

The historical kernel-backed radial working-cell diagnostic attached only the
BIAS1 true-bias box to A21 cells.  That is insufficient for the declared P4
objective.  This module materializes the same estimator-owned Joseph/reset
working cells separately for BIAS0, BIAS1 and BIAS2, using each family's own
qualified absolute true-bias component envelope.

H18 cells are family-independent and are emitted once per estimator/radial
cell.  A21 cells are emitted once per (estimator, radial, bias-family) tuple.
Every S=0 event retains actual committed anisotropic R_S; accelerometer and
magnetometer events retain the captured same-cell P/H/R geometry.  The AD state
box is only a first-exit working-domain enclosure and never replaces the hard
ball/IQC premises.

This closes the bias-family axis of the local working-cell cover.  It does not
claim the captured estimator image is the entire COMPLETE-BRMM source family;
source-uniform LDLT/first-exit promotion remains fail-closed.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_brmm_kernel_event_capture_bridge as BRIDGE
import ou3_p4_hard_entry_radial_ad_boxes as RADIAL
import ou3_p4_brmm_event_lineage_cover as LINEAGE
import ou3_p4_complete_brmm_source_cover_contract as COVER
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_brmm_reduced_event_master as MASTER
import ou3_p4_source_uniform_bias_prefix_lineage as BIAS

SCHEMA=1
QUALIFICATION='OU3_P4_ALL_BIAS_RADIAL_FIRST_EXIT_WORKING_CELLS_V1'
P3_DELTA=1e-18

def I(x):return Interval.point(float(x))
def projection_limit():
    p=Path(__file__).resolve().parents[2]/'tools'/'stability'/'ou3_proof_operating_domain.json'
    return float(json.loads(p.read_text())['normal_live']['active_accelerometer_bias_projection_limit_mps2'])
def true_bias_box(family):
    p=BIAS.family_parameters(family);cap=float(p['cap'])
    return [Interval(-cap,cap) for _ in range(3)]
def token(image,radial,ordinal,family=None):
    suffix=f'R{radial.token}_{ordinal}' if family is None else f'R{radial.token}_{family}_{ordinal}'
    return f'{image.source_token}:e{suffix}'
def cell(image,rc,sample,mode,radial,ordinal,state,pred,family=None):
    kwargs=dict(image=image,mode=mode,sample_index=0,event_ordinal=ordinal,kind=rc.kind,state=state,P=rc.P_before,
      dt_s=I(.005),pseudo_elapsed_s=I(.005*ordinal),radial_scale=radial.radial_scale,event_source_token=token(image,radial,ordinal,family),event_predecessor_token=pred)
    if mode=='A':
        if family not in BIAS.FAMILIES:raise ValueError('A21 working cell requires explicit BIAS0/1/2 family')
        kwargs['true_bias']=true_bias_box(family);kwargs['bias_projection_limit']=projection_limit()
    if rc.kind=='S_zero':pass
    elif rc.kind=='accelerometer':kwargs.update(R=rc.R,f_hat=sample.f_cog_body,R_hat=sample.R_wb)
    elif rc.kind=='magnetometer':
        mi=int(rc.magnetic_event_index or 0);kwargs.update(R=rc.R,m_body=sample.magnetometer_events_after_imu[mi].m_body)
    else:raise ValueError('Joseph working cells only')
    return COVER.source_cell_from_joint_image(**kwargs)
def build_lineage(image,captured,sample,mode,radial,family=None):
    state=RADIAL.radial_state_box(mode,radial.radial_scale);pred=image.predecessor_token;rows=[];fail=[]
    for ordinal,rc in enumerate(captured):
        if rc.kind not in ('S_zero','accelerometer','magnetometer'):continue
        c=cell(image,rc,sample,mode,radial,ordinal,state,pred,family)
        vf=COVER.validate_cell(c,require_estimator_provenance=True);fail.extend(f'{rc.kind}:cell:{x}' for x in vf)
        ev=EVENTS.source_joseph_event(**COVER.joseph_event_kwargs(c));em=MASTER.build_event_master(c);mf=MASTER.validate_event_master(c,em);fail.extend(f'{rc.kind}:master:{x}' for x in mf)
        rows.append({'kind':rc.kind,'event_token':c.source_token,'event_ordinal':ordinal,'family':family,'projection_branch':ev['bias_projection_branch'],'same_cell_P_H_R':bool(ev['same_P_H_R_cell']),'actual_RS':c.R_provenance==EVENTS.ACTUAL_RS_PROVENANCE if rc.kind=='S_zero' else None,'nonlinear_sector_count':len(em.nonlinear_sectors),'master_dimension':em.coordinate_dimension})
        pred=c.source_token
    return {'image_token':image.source_token,'mode':mode,'family':family,'radial_token':radial.token,'radial_scale':radial.radial_scale.as_list(),'events':rows,'event_count':len(rows),'valid':not fail,'failures':fail}
def build():
    sample,_,meta,_,images=BRIDGE.capture_smoke();root=LINEAGE.root_radial_cell();radials=LINEAGE.partition(root,2);records=[]
    for radial in radials:
      for image in images:
        records.append(build_lineage(image,meta['H_event_cells'],sample,'H',radial,None))
        for family in BIAS.FAMILIES:records.append(build_lineage(image,meta['A_event_cells'],sample,'A',radial,family))
    h=[r for r in records if r['mode']=='H'];a=[r for r in records if r['mode']=='A']
    family_counts={f:sum(1 for r in a if r['family']==f) for f in BIAS.FAMILIES}
    valid=bool(records) and all(r['valid'] and r['event_count']==3 for r in records)
    same_rs=all(next(x for x in r['events'] if x['kind']=='S_zero')['actual_RS'] is True for r in records)
    sectors=all(next(x for x in r['events'] if x['kind']=='accelerometer')['nonlinear_sector_count']>0 and next(x for x in r['events'] if x['kind']=='magnetometer')['nonlinear_sector_count']>0 for r in records)
    allfamilies=set(r['family'] for r in a)==set(BIAS.FAMILIES) and len(set(family_counts.values()))==1
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,
      'continuous_radial_partition_exactly_covers_unit_interval':LINEAGE.exact_partition_cover(radials),'H18_family_independent_cells_materialized_once':True,'A21_BIAS0_BIAS1_BIAS2_cells_materialized_separately':allfamilies,'each_bias_family_uses_own_authoritative_absolute_envelope':True,'BIAS0_or_BIAS2_inferred_from_BIAS1':False,'same_estimator_owned_P_H_R_retained':True,'all_Szero_cells_use_actual_applied_RS':same_rs,'all_nonlinear_cells_retain_structured_graph_sectors':sectors,'AD_box_is_first_exit_outer_enclosure_only':True,'hard_ball_geometry_replaced_by_AD_box':False,
      'record_count':len(records),'H18_record_count':len(h),'A21_record_count':len(a),'A21_family_counts':family_counts,'records':records,'all_working_cells_valid':valid,
      'captured_source_image_promoted_to_COMPLETE_BRMM_family':False,'production_source_uniform_LDLT_closed_here':False,'production_first_exit_retention_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'combine these all-family local working cells with the universal estimator/source selector relation and the source-uniform sequential endpoint/every-prefix LDLT; only then may first-exit working cells be promoted to reachable cells'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('continuous_radial_partition_exactly_covers_unit_interval','H18_family_independent_cells_materialized_once','A21_BIAS0_BIAS1_BIAS2_cells_materialized_separately','each_bias_family_uses_own_authoritative_absolute_envelope','same_estimator_owned_P_H_R_retained','all_Szero_cells_use_actual_applied_RS','all_nonlinear_cells_retain_structured_graph_sectors','AD_box_is_first_exit_outer_enclosure_only','all_working_cells_valid'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('BIAS0_or_BIAS2_inferred_from_BIAS1','hard_ball_geometry_replaced_by_AD_box','captured_source_image_promoted_to_COMPLETE_BRMM_family','production_source_uniform_LDLT_closed_here','production_first_exit_retention_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    counts=d.get('A21_family_counts',{})
    if set(counts)!=set(BIAS.FAMILIES) or len(set(counts.values()))!=1 or min(counts.values(),default=0)<=0:f.append('A21 family working-cell coverage incomplete')
    for r in d.get('records',[]):
        if r.get('valid') is not True or r.get('failures')!=[]:f.append('invalid working record')
    return list(dict.fromkeys(f))
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'records':d['record_count'],'families':d['A21_family_counts'],'valid':d['all_working_cells_valid'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
