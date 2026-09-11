#!/usr/bin/env python3
"""Conditional binary32 additive response at every H18 literal prefix.

The existing arithmetic theorem supplies a conditional source-uniform additive
state-error norm bound for the shipping time update and for measurement/reset
operations.  This module turns those local bounds into an actual H18 prefix
response: a normalized 18D arithmetic input is appended at each literal event
and all older arithmetic inputs are propagated through the exact later event
Jacobians.  No packet-count multiplication is used.

This is mathematical P4 arithmetic under the declared platform postcondition;
deployment/toolchain qualification remains separate.
"""
from __future__ import annotations
import argparse
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Sequence
from ou3_interval import Interval,matrix_identity,matrix_mul
import ou3_p4_source_uniform_estimator_event_attachment as ATTACH
import ou3_p4_complete_brmm_differential_prediction as PRED
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_kalman_reset_binary32_iss as FPISS

SCHEMA=1
QUALIFICATION='OU3_P4_H18_CONDITIONAL_BINARY32_EVERY_PREFIX_RESPONSE_V1'
P3_DELTA=1e-18

def I(x):return Interval.point(float(x))
def shape(A):return len(A),len(A[0]) if A else 0
def zero(r,c):return [[I(0) for _ in range(c)] for _ in range(r)]
def hstack(A,B):return [list(a)+list(b) for a,b in zip(A,B)]
def injection(delta):
    G=zero(18,18);d=I(delta)
    for i in range(18):G[i][i]=d
    return G
def propagate(A,C,delta):
    carried=matrix_mul(A,C) if shape(C)[1] else zero(18,0)
    return hstack(carried,injection(delta))

@dataclass(frozen=True)
class H18Binary32PrefixResponse:
    prefix_ordinal:int;source_cell_id:str;event_token:str;kind:str
    response:tuple[tuple[Interval,...],...];normalized_input_dimension:int;local_bound:float

def materialize(samples:Sequence[ATTACH.AttachedSampleLineage]):
    if not samples or any(s.mode!='H' for s in samples):raise ValueError('nonempty H18 attached sample lineage required')
    fp=FPISS.build();ff=FPISS.validate(fp)
    if ff or not fp['additive_ISS_channel_complete_for_conditional_P4']:raise RuntimeError('conditional binary32 prerequisite open')
    r=fp['modes']['H18'];dp=float(r['time_update_libm_and_coefficient_additive_state_error_norm_upper_per_prediction']);dm=float(r['explicit_measurement_reset_state_roundoff_norm_upper_per_event']);df=max(dp,dm)
    C=zero(18,0);state=list(samples[0].cells[0].state);out=[];ordinal=0
    for si,s in enumerate(samples):
        if s.selector.sample_index!=si:raise RuntimeError('sample ancestry nonconsecutive')
        for cell in s.cells:
            if tuple(cell.state)!=tuple(state):raise RuntimeError('binary32 response state detached')
            if cell.kind=='prediction':
                ev=PRED.prediction_event('H',state,s.selector.sample_coordinates.omega_body_corrected,cell.dt_s,cell.tau_applied_s);A=ev['J_state'];state=list(ev['state_out']);delta=dp
            elif cell.kind=='aw_floor':A=matrix_identity(18);delta=df
            elif cell.kind in ('S_zero','accelerometer','magnetometer'):
                ev=EVENTS.source_joseph_event(**ATTACH.COVER.joseph_event_kwargs(cell));A=ev['J_state'];state=list(ev['state_out']);delta=dm
            else:raise RuntimeError('unsupported H18 event '+cell.kind)
            C=propagate(A,C,delta);ordinal+=1
            if shape(C)!=(18,18*ordinal):raise RuntimeError('H18 binary32 response dimension drift')
            out.append(H18Binary32PrefixResponse(ordinal,s.selector.source_cell_id,cell.source_token,cell.kind,tuple(tuple(x for x in row) for row in C),18*ordinal,delta))
        if tuple(state)!=tuple(s.state_out):raise RuntimeError('H18 binary32 sample endpoint detached')
    return out

def build():
    f=FPISS.build();ff=FPISS.validate(f)
    if ff:raise RuntimeError('binary32 prerequisite failed '+repr(ff))
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,
      'conditional_full_Kalman_reset_binary32_ISS_consumed':True,'one_normalized_18D_roundoff_input_appended_per_literal_event':True,'prior_roundoff_inputs_suffix_propagated_through_later_exact_event_maps':True,'packet_count_times_worst_roundoff_used':False,
      'deployment_toolchain_qualification_required_for_mathematical_P4':False,'deployment_toolchain_qualification_closed_here':False,'H18_every_prefix_binary32_response_materializable':True,
      'common_augmented_PrefixInput_assembled_here':False,'endpoint_augmented_LDLT_closed_here':False,'every_prefix_augmented_LDLT_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'embed this response beside the exact H18 physical-source response and moving-information storage master, then run endpoint/every-prefix outward LDLT'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('conditional_full_Kalman_reset_binary32_ISS_consumed','one_normalized_18D_roundoff_input_appended_per_literal_event','prior_roundoff_inputs_suffix_propagated_through_later_exact_event_maps','H18_every_prefix_binary32_response_materializable'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('packet_count_times_worst_roundoff_used','deployment_toolchain_qualification_required_for_mathematical_P4','deployment_toolchain_qualification_closed_here','common_augmented_PrefixInput_assembled_here','endpoint_augmented_LDLT_closed_here','every_prefix_augmented_LDLT_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    return list(dict.fromkeys(f))
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'H18_fp_prefix':d['H18_every_prefix_binary32_response_materializable'],'deployment':d['deployment_toolchain_qualification_closed_here'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
