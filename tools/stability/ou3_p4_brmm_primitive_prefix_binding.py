#!/usr/bin/env python3
"""Bind provider translational primitives to literal P4 prefix certificates.

Each literal event is owned by one provider physical transition.  Multiple events
within one IMU sample repeat that transition binding; across samples primitive
out/in states are continuous and the one-time centered-S origin never changes.

For a cumulative endpoint/prefix LDLT, however, retaining only the *current*
transition's acceleration-moment IQC is insufficient.  Every distinct physical
transition encountered so far contributes one coupled 9D moment/radial sector to
the same augmented coordinate.  This bridge therefore validates both the event
binding and cumulative moment-witness ancestry.  It rejects a prefix that drops
an earlier physical transition even if its current sample is correctly bound.
"""
from __future__ import annotations
import argparse,json
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence
import ou3_p4_joint_brmm_prefix_ldlt as PREFIX
import ou3_brmm_centered_primitive_transition as PRIMITIVE

SCHEMA=2
QUALIFICATION='OU3_P4_BRMM_PRIMITIVE_PREFIX_BINDING_V2'
P3_DELTA=1e-18

@dataclass(frozen=True)
class PrimitivePrefixBinding:
    event_token:str;estimator_token:str;sample_index:int;source_transition_witness_id:str
    primitive_in_id:str;primitive_out_id:str;centered_S_origin_witness_id:str;moment_witness_id:str

def _nonempty(x):return isinstance(x,str) and bool(x)
def validate_binding_metadata(bindings:Sequence[PrimitivePrefixBinding])->list[str]:
    f=[]
    if not bindings:return ['nonempty primitive-prefix binding sequence required']
    origin=bindings[0].centered_S_origin_witness_id;previous=None;seen=set()
    for i,b in enumerate(bindings):
        for name in ('event_token','estimator_token','source_transition_witness_id','primitive_in_id','primitive_out_id','centered_S_origin_witness_id','moment_witness_id'):
            if not _nonempty(getattr(b,name)):f.append(f'binding {i}: missing {name}')
        if not isinstance(b.sample_index,int) or b.sample_index<0:f.append(f'binding {i}: invalid sample index')
        if b.event_token in seen:f.append(f'binding {i}: duplicate literal event token')
        seen.add(b.event_token)
        if b.source_transition_witness_id!=b.moment_witness_id:f.append(f'binding {i}: acceleration moments detached from provider transition')
        if b.centered_S_origin_witness_id!=origin:f.append(f'binding {i}: centered-S origin changed')
        if previous is not None:
            if b.sample_index<previous.sample_index:f.append(f'binding {i}: sample index went backwards')
            elif b.sample_index==previous.sample_index:
                for name in ('estimator_token','source_transition_witness_id','primitive_in_id','primitive_out_id','centered_S_origin_witness_id','moment_witness_id'):
                    if getattr(b,name)!=getattr(previous,name):f.append(f'binding {i}: same-sample literal events changed {name}')
            else:
                if b.sample_index!=previous.sample_index+1:f.append(f'binding {i}: physical sample ancestry skipped a sample')
                if previous.primitive_out_id!=b.primitive_in_id:f.append(f'binding {i}: primitive_out is not next primitive_in')
        previous=b
    return list(dict.fromkeys(f))

def _moment_ids(prefix:PREFIX.PrefixInput,i:int)->list[str]:
    try:return [b.witness_id for b in PREFIX._bindings(prefix,i)]
    except ValueError:return []

def validate_bound_prefixes(prefixes:Sequence[PREFIX.PrefixInput],bindings:Sequence[PrimitivePrefixBinding])->list[str]:
    f=validate_binding_metadata(bindings)
    if len(prefixes)!=len(bindings):return list(dict.fromkeys(f+['prefix/binding sequence length mismatch']))
    f.extend(PREFIX.validate_sequence(prefixes,require_physical_moment_sector=True,require_cumulative_moment_ancestry=True))
    expected_order=[]
    for i,(p,b) in enumerate(zip(prefixes,bindings)):
        if p.source_cell.source_token!=b.event_token:f.append(f'prefix {i}: literal event token detached from primitive binding')
        if p.image.source_token!=b.estimator_token:f.append(f'prefix {i}: estimator image detached from primitive binding')
        if p.source_cell.sample_index!=b.sample_index:f.append(f'prefix {i}: sample index detached from primitive binding')
        if b.moment_witness_id not in expected_order:expected_order.append(b.moment_witness_id)
        ids=_moment_ids(p,i)
        # Legacy one-sector form is valid only while one physical transition has
        # appeared.  Multi-transition production prefixes must expose actual IDs.
        if len(expected_order)==1 and ids==[f'legacy-prefix-{i}']:
            pass
        elif ids!=expected_order:
            f.append(f'prefix {i}: cumulative moment witness lineage detached expected={expected_order!r} got={ids!r}')
    return list(dict.fromkeys(f))

def certify_bound_prefixes(prefixes,bindings):
    failures=validate_bound_prefixes(prefixes,bindings)
    if failures:return {'closed':False,'validation_failures':failures,'endpoint_closed':False,'every_prefix_closed':False,'primitive_ancestry_closed':False}
    result=dict(PREFIX.certify(prefixes,require_physical_moment_sector=True,require_cumulative_moment_ancestry=True));result['primitive_ancestry_closed']=True;result['provider_transition_binding_count']=len(bindings);result['closed']=bool(result.get('closed') and result['primitive_ancestry_closed']);return result

def _smoke_bindings():
    return [
      PrimitivePrefixBinding('im0:e0','im0',0,'tr0','q0','q1','live-origin','tr0'),
      PrimitivePrefixBinding('im0:e1','im0',0,'tr0','q0','q1','live-origin','tr0'),
      PrimitivePrefixBinding('im1:e0','im1',1,'tr1','q1','q2','live-origin','tr1')]

def build():
    primitive=PRIMITIVE.build();pf=PRIMITIVE.validate(primitive);prefix=PREFIX.build();xf=PREFIX.validate(prefix)
    if pf or xf:raise RuntimeError(f'primitive-prefix prerequisites failed primitive={pf} prefix={xf}')
    smoke=_smoke_bindings();sf=validate_binding_metadata(smoke)
    if sf:raise RuntimeError('primitive-prefix binding smoke failed: '+repr(sf))
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,
      'same_provider_transition_owns_primitives_and_acceleration_moments':True,'literal_events_within_one_sample_share_one_physical_transition':True,
      'cross_sample_primitive_out_equals_next_primitive_in_required':True,'one_live_centered_S_origin_retained_indefinitely':True,
      'estimator_event_and_physical_transition_tokens_jointly_checked':True,'production_prefix_requires_coupled_moment_sector':True,
      'production_prefix_requires_same_radial_coordinate':True,'production_prefix_requires_all_prior_transition_moment_witnesses':True,
      'same_sample_events_do_not_duplicate_physical_moment_witness':True,'favorable_branch_selection_used':False,'wordwise_S_rezero_used':False,
      'independent_moment_boxes_used':False,'independent_axis_boxes_used':False,'shipping_filter_changed':False,
      'smoke_binding_count':len(smoke),'smoke_binding_validation_closed':True,
      'production_complete_provider_lineages_bound_here':False,'production_endpoint_LDLT_closed_here':False,'production_every_prefix_LDLT_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'for every synchronized post-prediction prefix, create one MomentSectorBinding per distinct provider transition encountered so far using its exact moment witness and common radial coordinate; then run certify_bound_prefixes'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('canonical_source')!='COMPLETE_BRMM_NORMAL_LIVE_WORD':f.append('canonical source changed')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('same_provider_transition_owns_primitives_and_acceleration_moments','literal_events_within_one_sample_share_one_physical_transition','cross_sample_primitive_out_equals_next_primitive_in_required','one_live_centered_S_origin_retained_indefinitely','estimator_event_and_physical_transition_tokens_jointly_checked','production_prefix_requires_coupled_moment_sector','production_prefix_requires_same_radial_coordinate','production_prefix_requires_all_prior_transition_moment_witnesses','same_sample_events_do_not_duplicate_physical_moment_witness','smoke_binding_validation_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('favorable_branch_selection_used','wordwise_S_rezero_used','independent_moment_boxes_used','independent_axis_boxes_used','shipping_filter_changed','production_complete_provider_lineages_bound_here','production_endpoint_LDLT_closed_here','production_every_prefix_LDLT_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('smoke_binding_count')!=3:f.append('smoke binding count changed')
    return list(dict.fromkeys(f))
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'primitive_prefix_binding':d['smoke_binding_validation_closed'],'cumulative_moments':d['production_prefix_requires_all_prior_transition_moment_witnesses'],'production_bound':d['production_complete_provider_lineages_bound_here'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
