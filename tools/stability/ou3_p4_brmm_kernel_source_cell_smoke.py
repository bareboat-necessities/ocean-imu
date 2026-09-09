#!/usr/bin/env python3
"""Kernel-backed estimator-owned SourceCoverCell regression for P4.

Passive Riccati-event captures from the trusted complete-window kernel are bound
to EVERY same-sample joint-estimator Image.  For each image and for both H18 and
A21 construct the literal theorem-cell chain

  prediction -> aw_floor -> S_zero -> accelerometer -> magnetometer.

This regression is intentionally rooted at zero physical error / zero true
residual accelerometer bias. At that root all nonlinear measurement corrections
and finite reset defects vanish, so the passive kernel covariance chain and the
physical event map agree exactly without assuming covariance consistency.

The purpose is end-to-end wiring, not source-uniform promotion:
trusted kernel capture -> every estimator image -> SourceCoverCells -> same-cell
Joseph event -> structured reduced event masters. Full hard-entry state
transport and nonzero BIAS1 continuation remain mandatory production blockers.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_brmm_kernel_event_capture_bridge as BRIDGE
import ou3_p4_complete_brmm_source_cover_contract as COVER
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_brmm_reduced_event_master as MASTER

SCHEMA=2
QUALIFICATION='OU3_P4_BRMM_KERNEL_SOURCE_CELL_ZERO_ROOT_SMOKE_V2'

def I(x):return Interval.point(float(x))
def zero_state(n):return [I(0) for _ in range(n)]
def _limit():
    p=Path(__file__).resolve().parents[2]/'tools'/'stability'/'ou3_proof_operating_domain.json'
    d=json.loads(p.read_text());return float(d['normal_live']['active_accelerometer_bias_projection_limit_mps2'])
def build_cells(mode,captured,image,sample):
    n=18 if mode=='H' else 21;x=zero_state(n);beta=[I(0),I(0),I(0)] if mode=='A' else None;limit=_limit() if mode=='A' else None
    out=[];pred=image.predecessor_token
    for ordinal,rc in enumerate(captured):
        token=f'{image.source_token}:e{ordinal}'
        common=dict(image=image,mode=mode,sample_index=0,event_ordinal=ordinal,kind=rc.kind,state=x,P=rc.P_before,
                    dt_s=I(.005),pseudo_elapsed_s=I(.005*ordinal),radial_scale=Interval(0,0),event_source_token=token,
                    event_predecessor_token=pred,true_bias=beta,bias_projection_limit=limit)
        if rc.kind=='S_zero':cell=COVER.source_cell_from_joint_image(**common)
        elif rc.kind=='accelerometer':cell=COVER.source_cell_from_joint_image(**common,R=rc.R,f_hat=sample.f_cog_body,R_hat=sample.R_wb)
        elif rc.kind=='magnetometer':
            m=sample.magnetometer_events_after_imu[int(rc.magnetic_event_index or 0)].m_body
            cell=COVER.source_cell_from_joint_image(**common,R=rc.R,m_body=m)
        else:cell=COVER.source_cell_from_joint_image(**common)
        out.append(cell);pred=token
    return out
def _zero_state_event_closed(cell):
    if cell.kind not in ('S_zero','accelerometer','magnetometer'):return True
    ev=EVENTS.source_joseph_event(**COVER.joseph_event_kwargs(cell));return all(v.lo<=0<=v.hi for v in ev['state_out'])
def summarize(cells):
    vf=[];masters=0
    for i,c in enumerate(cells):
        vf.extend(f'{i}:{x}' for x in COVER.validate_cell(c,require_estimator_provenance=True))
        if i and c.predecessor_token!=cells[i-1].source_token:vf.append(f'{i}: predecessor detached')
        if not _zero_state_event_closed(c):vf.append(f'{i}: zero root not invariant under Joseph event')
        if c.kind in ('S_zero','accelerometer','magnetometer'):
            em=MASTER.build_event_master(c);vf.extend(f'{i}:master:{x}' for x in MASTER.validate_event_master(c,em));masters+=1
    return {'literal_kinds':[c.kind for c in cells],'event_count':len(cells),'Joseph_master_count':masters,
            'all_cells_valid':not vf,'validation_failures':vf,
            'all_event_tokens_estimator_owned':all(COVER.estimator_owns_event_token(c) for c in cells),
            'literal_predecessor_chain_closed':all(cells[i].predecessor_token==cells[i-1].source_token for i in range(1,len(cells)))}
def build():
    sample,_,meta,_,images=BRIDGE.capture_smoke();records=[]
    for image in images:
        hs=summarize(build_cells('H',meta['H_event_cells'],image,sample));aa=summarize(build_cells('A',meta['A_event_cells'],image,sample))
        records.append({'image_token':image.source_token,'H18':hs,'A21':aa})
    all_valid=bool(records) and all(r[m]['all_cells_valid'] for r in records for m in ('H18','A21'))
    all_five=bool(records) and all(r[m]['event_count']==5 for r in records for m in ('H18','A21'))
    all_masters=bool(records) and all(r[m]['Joseph_master_count']==3 for r in records for m in ('H18','A21'))
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'trusted_kernel_passive_event_cells_consumed':True,'same_sample_complete_joint_estimator_image_family_consumed':True,
      'zero_physical_error_root_only':True,'zero_true_bias_root_only':True,'joint_image_lineage_count':len(records),'lineages':records,
      'every_joint_image_has_both_mode_literal_five_event_chains':all_five,
      'all_estimator_owned_source_cells_valid':all_valid,'all_three_Joseph_event_master_classes_built_for_every_image_and_mode':all_masters,
      'favorable_joint_estimator_successor_selected':False,'covariance_membership_used_as_state_hypothesis':False,'point_zero_root_may_promote_source_uniform_cover':False,
      'production_nonzero_hard_entry_state_transport_closed_here':False,'production_nonzero_BIAS1_continuation_closed_here':False,
      'production_all_joint_estimator_successors_covered_here':False,'production_complete_event_lineage_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'replace the zero-root state fixture by outward physical-error transport from the full qualified hard-entry ball family for every retained image/radial cell; then carry physical BIAS1 beta/driver and pass each kernel-backed Joseph cell to the same structured master/reset/projection pipeline'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('trusted_kernel_passive_event_cells_consumed','same_sample_complete_joint_estimator_image_family_consumed','zero_physical_error_root_only','zero_true_bias_root_only','every_joint_image_has_both_mode_literal_five_event_chains','all_estimator_owned_source_cells_valid','all_three_Joseph_event_master_classes_built_for_every_image_and_mode'):
        if d.get(k) is not True:f.append(k+' not true')
    if int(d.get('joint_image_lineage_count',0))<=0:f.append('no joint image lineages')
    for r in d.get('lineages',[]):
        for mode in ('H18','A21'):
            m=r.get(mode,{})
            for k in ('all_cells_valid','all_event_tokens_estimator_owned','literal_predecessor_chain_closed'):
                if m.get(k) is not True:f.append(r.get('image_token','?')+' '+mode+' '+k+' not true')
            if m.get('literal_kinds')!=['prediction','aw_floor','S_zero','accelerometer','magnetometer']:f.append(mode+' literal order changed')
            if m.get('validation_failures')!=[]:f.append(mode+' cell validation failures')
    for k in ('favorable_joint_estimator_successor_selected','covariance_membership_used_as_state_hypothesis','point_zero_root_may_promote_source_uniform_cover','production_nonzero_hard_entry_state_transport_closed_here','production_nonzero_BIAS1_continuation_closed_here','production_all_joint_estimator_successors_covered_here','production_complete_event_lineage_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'kernel_cells':d['all_estimator_owned_source_cells_valid'],'lineages':d['joint_image_lineage_count'],'all_masters':d['all_three_Joseph_event_master_classes_built_for_every_image_and_mode'],'production':d['production_complete_event_lineage_closed_here'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
