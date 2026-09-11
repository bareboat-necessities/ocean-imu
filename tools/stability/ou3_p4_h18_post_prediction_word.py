#!/usr/bin/env python3
"""Canonical H18 P4 word on post-prediction information-storage boundaries.

The uniformly coercive compatible storage is P^{-1} at shipping prediction
boundaries.  Canonical P3 likewise closes a complete information word at the
following prediction.  This module slices the exact synchronized H18 source-
indexed transport onto precisely that topology:

  boundary 0: immediately after prediction of sample k
  prefixes:   every subsequent literal shipping event/rebase
  boundary N: immediately after the following selected prediction.

For a sequence of synchronized samples, the first prediction is consumed only
to establish the entrance boundary.  The transport then contains the remainder
of the first sample, every intervening sample/event, every source-coordinate
rebase, and stops immediately after the final sample's prediction.  Events after
the terminal prediction belong to the next word.

The structural prefix decomposition is recomputed from the sliced source-
indexed epsilon nodes; it is not obtained by subtracting two interval prefix
hulls.  Covariance nodes are sliced in the identical event order, so every
boundary metric remains the actual shipping covariance.  This closes the word
boundary/transport alignment, not the augmented LDLT or first-exit inequality.
"""
from __future__ import annotations
import argparse,json
from typing import Sequence

from ou3_interval import Interval
import ou3_p4_h18_source_indexed_prefix_transport as HTR
import ou3_p4_structural_prefix_transport as STRUCT
import ou3_p4_information_storage_coercivity as STORAGE
import ou3_p4_source_uniform_estimator_event_attachment as ATTACH

SCHEMA=1
QUALIFICATION='OU3_P4_H18_POST_PREDICTION_WORD_TRANSPORT_V1'
P3_DELTA=1e-18

def covariance_nodes(samples:Sequence[ATTACH.AttachedSampleLineage]):
    if not samples:raise ValueError('nonempty H18 sample lineage required')
    nodes=[samples[0].cells[0].P]
    for si,s in enumerate(samples):
        if si:nodes.append(nodes[-1])  # source-coordinate rebase; physical P unchanged
        trusted={c.event_index_in_sample:c for c in s.selector.H_event_cells}
        for cell in s.cells:nodes.append(trusted[cell.event_ordinal].P_after)
    return nodes

def build_word(samples:Sequence[ATTACH.AttachedSampleLineage]):
    if len(samples)<2:raise ValueError('post-prediction word requires entrance and following prediction samples')
    t=HTR.build_transport(samples);events=t['events'];P=covariance_nodes(samples)
    if len(P)!=len(events)+1 or len(t['eps_nodes'])!=len(events)+1 or len(t['embeddings'])!=len(events)+1:raise RuntimeError('transport node counts detached')
    preds=[i for i,e in enumerate(events) if e['kind']=='prediction']
    if len(preds)!=len(samples):raise RuntimeError('one prediction per synchronized sample required')
    first,last=preds[0],preds[-1]
    if first!=0:raise RuntimeError('synchronized H18 sequence no longer starts with prediction')
    if last<=first:raise RuntimeError('following prediction missing')
    a=first+1;b=last
    sliced_events=events[a:b+1]
    sliced_eps=t['eps_nodes'][a:b+2]
    sliced_emb=t['embeddings'][a:b+2]
    sliced_P=P[a:b+2]
    d=STRUCT.prefix_decompositions(sliced_events,sliced_emb,sliced_eps)
    if len(d)!=len(sliced_events) or len(sliced_P)!=len(sliced_events)+1:raise RuntimeError('sliced prefix decomposition mismatch')
    entrance_P=P[first+1];terminal_P=P[last+1]
    if entrance_P!=sliced_P[0] or terminal_P!=sliced_P[-1]:raise RuntimeError('post-prediction covariance boundary detached')
    # The terminal event itself must be the following physical prediction; any
    # later floor/measurement event is deliberately left to the next word.
    if sliced_events[-1]['kind']!='prediction':raise RuntimeError('word did not terminate at following prediction')
    return {'events':sliced_events,'decompositions':d,'covariance_nodes':sliced_P,'entrance_P':entrance_P,'terminal_P':terminal_P,'entrance_prediction_event_index':first,'terminal_prediction_event_index':last,'source_event_count':len(sliced_events),'all_prefix_identity_residuals_contain_zero':all(x['transport']['identity_residual_contains_zero'] for x in d)}

def build():
    h=HTR.build();hf=HTR.validate(h);s=STORAGE.build();sf=STORAGE.validate(s);st=STRUCT.build();stf=STRUCT.validate(st)
    if hf or sf or stf:raise RuntimeError(f'post-prediction word prerequisites failed H18={hf} storage={sf} structural={stf}')
    w=build_word(HTR.smoke_objects())
    closed=bool(s['uniform_word_boundary_storage_coercivity_closed'] and s['consecutive_word_metric_comparison_mu']==1.0 and w['source_event_count']>0 and w['all_prefix_identity_residuals_contain_zero'])
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,
      'canonical_P3_following_prediction_topology_consumed':True,'uniformly_coercive_post_prediction_storage_consumed':True,
      'word_entrance_is_actual_shipping_post_prediction_boundary':True,'word_endpoint_is_actual_following_prediction_boundary':True,
      'events_after_terminal_prediction_deferred_to_successor_word':True,'source_coordinate_rebases_inside_word_retained':True,
      'sliced_structural_prefix_decomposition_recomputed_not_subtracted':True,'same_shipping_covariance_node_at_shared_consecutive_boundary':True,
      'metric_comparison_mu':s['consecutive_word_metric_comparison_mu'],'H18_post_prediction_word_transport_closed':closed,
      'smoke_event_count':w['source_event_count'],'smoke_all_prefix_identity_residuals_contain_zero':w['all_prefix_identity_residuals_contain_zero'],
      'endpoint_augmented_LDLT_closed_here':False,'every_prefix_augmented_LDLT_closed_here':False,'first_exit_retention_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'on this exact post-prediction word, attach common correlated BRMM moment/radial and conditional binary32 maps to every prefix and use actual entrance/terminal P^-1 in endpoint/every-prefix augmented LDLT'}

def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('canonical_P3_following_prediction_topology_consumed','uniformly_coercive_post_prediction_storage_consumed','word_entrance_is_actual_shipping_post_prediction_boundary','word_endpoint_is_actual_following_prediction_boundary','events_after_terminal_prediction_deferred_to_successor_word','source_coordinate_rebases_inside_word_retained','sliced_structural_prefix_decomposition_recomputed_not_subtracted','same_shipping_covariance_node_at_shared_consecutive_boundary','H18_post_prediction_word_transport_closed','smoke_all_prefix_identity_residuals_contain_zero'):
        if d.get(k) is not True:f.append(k+' not true')
    if d.get('metric_comparison_mu')!=1.0:f.append('metric comparison penalty changed')
    for k in ('endpoint_augmented_LDLT_closed_here','every_prefix_augmented_LDLT_closed_here','first_exit_retention_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    return list(dict.fromkeys(f))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    from pathlib import Path;path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'closed':d['H18_post_prediction_word_transport_closed'],'events':d['smoke_event_count'],'mu':d['metric_comparison_mu'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
