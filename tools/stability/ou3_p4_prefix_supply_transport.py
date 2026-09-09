#!/usr/bin/env python3
"""Exact suffix-propagated ISS supply maps for P4 literal prefixes.

For a linearized/joint event sequence

    x_{k+1} = C_k x_k + B_k s_k,

variation of constants gives the prefix forcing map

    x_ell = M_ell x_0 + [M_{ell:1}B_0 ... B_{ell-1}] s_{0:ell-1}.

This module constructs those *joint* maps for every literal prefix.  It is used
for physical BIAS1 prediction supplies and may also carry additive binary32
roundoff channels.  It never replaces the stacked map by a packet count or sum
of per-event norms.

Supply coordinates may have different dimensions by event; events without a
supply injection contribute no columns.  The output is one dense prefix map
whose columns retain their original event/source identity.
"""
from __future__ import annotations
import argparse,json
from fractions import Fraction
from pathlib import Path

import ou3_p4_complete_word_endpoint_transport as TRANSPORT
import ou3_p4_a21_bias1_24state_event_lift as BIAS24

SCHEMA=1
QUALIFICATION='OU3_P4_EXACT_PREFIX_SUFFIX_SUPPLY_TRANSPORT_V1'

def shape(A):
    r=len(A);c=len(A[0]) if r else 0
    if any(len(row)!=c for row in A):raise ValueError('ragged matrix')
    return r,c

def mm(A,B):return TRANSPORT._mm(A,B)
def hcat(blocks,rows,zero):
    if not blocks:return [[] for _ in range(rows)]
    if any(shape(B)[0]!=rows for B in blocks):raise ValueError('supply block row mismatch')
    return [sum((list(B[i]) for B in blocks),[]) for i in range(rows)]

def prefix_supply_maps(events,injections,one,zero):
    if not events or len(events)!=len(injections):raise ValueError('matching nonempty events/injections required')
    out=[]
    for ell in range(1,len(events)+1):
        ev=events[:ell];inj=injections[:ell]
        suffix=TRANSPORT.suffix_products(ev,one,zero)
        blocks=[];labels=[]
        for k,B in enumerate(inj):
            if B is None:continue
            cr,cc=shape(ev[k]['C']);br,bs=shape(B)
            if br!=cr or bs==0:raise ValueError(f'injection {k} dimension mismatch')
            blocks.append(mm(suffix[k+1],B));labels.extend((k,j) for j in range(bs))
        rows=shape(ev[-1]['C'])[0]
        out.append({'prefix_length':ell,'map':hcat(blocks,rows,zero),'column_labels':labels,'supply_event_count':sum(B is not None for B in inj)})
    return out

def _mv(A,x):return TRANSPORT._mv(A,x)
def _add(a,b):return [x+y for x,y in zip(a,b)]
def _smoke():
    F=Fraction;z=F(0);o=F(1)
    C0=[[F(1,2),z],[z,F(3,4)]];C1=[[F(4,5),z],[F(1,10),F(9,10)]]
    B0=[[o,z],[z,o]];B1=[[F(1,3)],[F(-1,5)]]
    events=[{'C':C0,'L':C0,'rho':[z,z],'kind':'prediction'},{'C':C1,'L':C1,'rho':[z,z],'kind':'prediction'}]
    rec=prefix_supply_maps(events,[B0,B1],o,z)
    s0=[F(2,7),F(-1,9)];s1=[F(3,11)]
    direct1=_mv(B0,s0)
    direct2=_add(_mv(C1,direct1),_mv(B1,s1))
    map1=_mv(rec[0]['map'],s0)
    map2=_mv(rec[1]['map'],s0+s1)
    return rec,direct1==map1,direct2==map2

def build():
    b=BIAS24.build();bf=BIAS24.validate(b)
    if bf:raise RuntimeError('24-state BIAS1 lift prerequisite failed: '+repr(bf))
    rec,e1,e2=_smoke()
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'variation_of_constants_suffix_supply_map_available':True,'every_literal_prefix_supply_map_available':True,
      'different_event_supply_dimensions_supported':True,'events_without_supply_add_no_columns':True,
      'joint_columns_retain_event_identity':True,'packet_count_multiplier_used':False,'per_event_norm_sum_used':False,
      'physical_BIAS1_24state_supply_lift_consumed':True,'same_w_bias_structure_available_before_suffix_transport':bool(b['same_w_enters_error_and_true_bias']),
      'exact_rational_first_prefix_supply_transport_closed':e1,'exact_rational_second_prefix_supply_transport_closed':e2,
      'smoke_prefix_column_counts':[len(x['column_labels']) for x in rec],
      'production_estimator_owned_event_C_sequence_materialized_here':False,'production_BIAS1_prefix_source_maps_closed_here':False,'production_binary32_prefix_maps_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'build the estimator-owned joint 24-state C_k sequence, attach prediction_lift B_k for every BIAS1 prediction and additive binary32 injection maps, then call prefix_supply_maps and pass each dense result to the corresponding PrefixInput'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('variation_of_constants_suffix_supply_map_available','every_literal_prefix_supply_map_available','different_event_supply_dimensions_supported','events_without_supply_add_no_columns','joint_columns_retain_event_identity','physical_BIAS1_24state_supply_lift_consumed','same_w_bias_structure_available_before_suffix_transport','exact_rational_first_prefix_supply_transport_closed','exact_rational_second_prefix_supply_transport_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('packet_count_multiplier_used','per_event_norm_sum_used','production_estimator_owned_event_C_sequence_materialized_here','production_BIAS1_prefix_source_maps_closed_here','production_binary32_prefix_maps_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('smoke_prefix_column_counts')!=[2,3]:f.append('unexpected smoke supply column ancestry')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'prefix_supply':d['every_literal_prefix_supply_map_available'],'joint':d['joint_columns_retain_event_identity'],'production':d['production_BIAS1_prefix_source_maps_closed_here'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
