#!/usr/bin/env python3
"""Canonical A21 P4 words on post-prediction storage boundaries.

For each BIAS0/BIAS1/BIAS2 family, construct the exact 24-state source-indexed
A21 transport from one shipping prediction boundary to the following prediction
boundary.  The first prediction is used only to establish the entrance state and
covariance; its [w,m_tau] block is therefore absorbed into the entrance state and
is NOT charged again as a word input.  Every later prediction inside the word
adds exactly one new shared six-dimensional [w,m_tau] block.

Within the word every accepted Joseph/reset and the immediately following bias
projection remain distinct literal prefixes.  Projection leaves covariance
unchanged but retains the same physical b_true coordinate in the 24-state map.
The endpoint is the following actual shipping prediction; later events belong to
the successor word.

This aligns A21 with the uniformly coercive P^{-1} storage and exact mu=1
consecutive-word compatibility.  It does not yet attach the correlated BRMM
moment/radial sector, binary32 map, or augmented LDLT.
"""
from __future__ import annotations
import argparse,copy,json
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval
import ou3_p4_a21_source_indexed_prefix_transport as A21
import ou3_p4_source_uniform_bias_prefix_lineage as BIAS
import ou3_p4_structural_prefix_transport as STRUCT
import ou3_p4_information_storage_coercivity as STORAGE
import ou3_p4_source_uniform_estimator_event_attachment as ATTACH

SCHEMA=1
QUALIFICATION='OU3_P4_A21_POST_PREDICTION_WORD_TRANSPORT_V1'
P3_DELTA=1e-18

def trusted(sample):return {c.event_index_in_sample:c for c in sample.selector.A_event_cells}

def append_event(events,eps,emb,source_maps,rec,state,selector,Bcur,Bnew=None):
    Bcur=A21.propagate_source(rec['C'],Bcur,Bnew)
    events.append(rec);source_maps.append(copy.deepcopy(Bcur));eps.append(A21.epsilon_aw(state,selector));emb.append(A21.Eaw24())
    return Bcur

def build_word(samples:Sequence[ATTACH.AttachedSampleLineage],family:str):
    if family not in BIAS.FAMILIES:raise ValueError('family must be BIAS0/BIAS1/BIAS2')
    if len(samples)<2:raise ValueError('post-prediction A21 word needs entrance and following prediction samples')
    selectors=[s.selector for s in samples];cert=BIAS.attach_to_selector_lineage(family,selectors)
    for s in samples:
        expected=cert.at(s.selector.source_cell_id)
        for c in s.cells:
            if c.true_bias!=expected:raise RuntimeError('A21 true-bias prefix detached before word slicing')
    first=samples[0];tm=trusted(first);pcell=first.cells[0]
    if pcell.kind!='prediction':raise RuntimeError('synchronized A21 sample no longer begins with prediction')
    prec,state=A21.prediction_record(pcell,tm[0],first.selector,cert)
    # Entrance is AFTER this prediction.  Its driver belongs to the prior word.
    entrance_P=tm[0].P_after
    events=[];source_maps=[];Bcur=A21.zmat(A21.N,0)
    eps=[A21.epsilon_aw(state,first.selector)];emb=[A21.Eaw24()];Pnodes=[entrance_P]

    def consume_nonprediction(sample,start_index):
        nonlocal Bcur,state
        tm0=trusted(sample)
        for cell in sample.cells[start_index:]:
            if tuple(cell.state)!=tuple(state):raise RuntimeError('A21 sliced event state detached')
            if cell.kind=='prediction':raise RuntimeError('unexpected prediction inside sample suffix')
            if cell.kind=='aw_floor':
                rec,out=A21.floor_record(cell);state=list(out);Bcur=append_event(events,eps,emb,source_maps,rec,state,sample.selector,Bcur);Pnodes.append(tm0[cell.event_ordinal].P_after)
            else:
                smooth,proj,sout,pout=A21.split_records(cell)
                state=list(sout);Bcur=append_event(events,eps,emb,source_maps,smooth,state,sample.selector,Bcur);Pnodes.append(tm0[cell.event_ordinal].P_after)
                before=eps[-1];state=list(pout);Bcur=append_event(events,eps,emb,source_maps,proj,state,sample.selector,Bcur);Pnodes.append(Pnodes[-1])
                if eps[-1]!=before:raise RuntimeError('projection changed source-indexed Phi coordinate')

    consume_nonprediction(first,1)
    # Complete any intermediate samples, then stop immediately after the final prediction.
    for si,sample in enumerate(samples[1:],start=1):
        neweps=A21.epsilon_aw(state,sample.selector);rb=A21.rebase24(eps[-1],neweps);rb['token']=sample.selector.source_cell_id+':rebase'
        Bcur=A21.propagate_source(rb['C'],Bcur);events.append(rb);source_maps.append(copy.deepcopy(Bcur));eps.append(neweps);emb.append(A21.Eaw24());Pnodes.append(Pnodes[-1])
        tm0=trusted(sample);pc=sample.cells[0]
        if pc.kind!='prediction':raise RuntimeError('following sample lost prediction entrance')
        rec,out=A21.prediction_record(pc,tm0[0],sample.selector,cert);state=list(out)
        Bcur=append_event(events,eps,emb,source_maps,rec,state,sample.selector,Bcur,rec['prediction_supply_block']);Pnodes.append(tm0[0].P_after)
        if si==len(samples)-1:break
        consume_nonprediction(sample,1)

    if events[-1]['kind']!='prediction':raise RuntimeError('A21 word did not terminate at following prediction')
    decomp=STRUCT.prefix_decompositions(events,emb,eps)
    if len(decomp)!=len(events) or len(Pnodes)!=len(events)+1 or len(source_maps)!=len(events):raise RuntimeError('A21 post-prediction prefix counts detached')
    expected_cols=0;dim_ok=True
    for e,B in zip(events,source_maps):
        if e['kind']=='prediction':expected_cols+=6
        if A21.shape(B)[1]!=expected_cols:dim_ok=False
    return {
      'family':family,'events':events,'decompositions':decomp,'source_prefix_maps':source_maps,'covariance_nodes':Pnodes,
      'entrance_P':Pnodes[0],'terminal_P':Pnodes[-1],'event_count':len(events),
      'first_prediction_supply_absorbed_into_boundary_state':True,
      'one_new_6D_supply_block_per_in_word_prediction':dim_ok,
      'all_prefix_identity_residuals_contain_zero':all(x['transport']['identity_residual_contains_zero'] for x in decomp),
      'endpoint_is_following_prediction':events[-1]['kind']=='prediction',
      'projection_prefixes_retained':any(e['kind']=='bias_projection' for e in events)}

def build():
    a=A21.build();af=A21.validate(a);s=STORAGE.build();sf=STORAGE.validate(s);st=STRUCT.build();stf=STRUCT.validate(st)
    if af or sf or stf:raise RuntimeError(f'A21 word prerequisites failed transport={af} storage={sf} structural={stf}')
    rows={f:build_word(A21.smoke_objects(f),f) for f in BIAS.FAMILIES}
    closed=all(r['first_prediction_supply_absorbed_into_boundary_state'] and r['one_new_6D_supply_block_per_in_word_prediction'] and r['all_prefix_identity_residuals_contain_zero'] and r['endpoint_is_following_prediction'] and r['projection_prefixes_retained'] for r in rows.values()) and s['uniform_word_boundary_storage_coercivity_closed'] and s['consecutive_word_metric_comparison_mu']==1.0
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,
      'all_BIAS0_BIAS1_BIAS2_post_prediction_words_materialized':closed,
      'uniformly_coercive_post_prediction_storage_consumed':True,'metric_comparison_mu':s['consecutive_word_metric_comparison_mu'],
      'first_prediction_driver_not_double_counted':True,'Joseph_reset_and_projection_remain_distinct_literal_prefixes':True,
      'same_family_true_bias_recurrence_retained':True,'source_coordinate_rebases_retained':True,
      'A21_post_prediction_word_transport_closed':closed,
      'smoke_by_family':{f:{k:v for k,v in r.items() if k not in ('events','decompositions','source_prefix_maps','covariance_nodes','entrance_P','terminal_P')} for f,r in rows.items()},
      'BRMM_moment_radial_sector_attached_here':False,'binary32_prefix_map_attached_here':False,'endpoint_augmented_LDLT_closed_here':False,'every_prefix_augmented_LDLT_closed_here':False,'first_exit_retention_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'assemble H18 and A21 post-prediction prefixes into the common augmented coordinate; attach same-history BRMM moment/radial and conditional binary32 maps, then run endpoint/every-prefix outward LDLT'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('all_BIAS0_BIAS1_BIAS2_post_prediction_words_materialized','uniformly_coercive_post_prediction_storage_consumed','first_prediction_driver_not_double_counted','Joseph_reset_and_projection_remain_distinct_literal_prefixes','same_family_true_bias_recurrence_retained','source_coordinate_rebases_retained','A21_post_prediction_word_transport_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    if d.get('metric_comparison_mu')!=1.0:f.append('metric comparison changed')
    for k in ('BRMM_moment_radial_sector_attached_here','binary32_prefix_map_attached_here','endpoint_augmented_LDLT_closed_here','every_prefix_augmented_LDLT_closed_here','first_exit_retention_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if set(d.get('smoke_by_family',{}))!=set(BIAS.FAMILIES):f.append('bias family word set incomplete')
    return list(dict.fromkeys(f))
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'closed':d['A21_post_prediction_word_transport_closed'],'families':list(d['smoke_by_family']),'mu':d['metric_comparison_mu'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
