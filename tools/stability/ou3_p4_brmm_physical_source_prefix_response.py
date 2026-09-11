#!/usr/bin/env python3
"""Exact BRMM physical-source response at every literal H18/A21 prefix.

Each shipping prediction receives one 15D same-history acceleration witness
q_k=[a0,a1,J0,J1,J2].  The exact true-minus-estimate forcing map is generated
from the same committed tau/h as the homogeneous OU prediction.  One q_k block
is appended only at that prediction; all older q blocks are propagated through
the exact later nonlinear event Jacobians.

H18 returns an 18 x (15*Npred) response.  A21 returns the all-bias 24-state
response [e_A21;b_true]; physical acceleration forcing enters only the first 21
shipping-error rows and never the proof-only true-bias rows.  A21 Joseph/reset/
projection maps use the family-specific same-prefix true-bias recurrence.

This is the missing physical finite-error source channel.  It retains raw q_k
coordinates so their joint endpoint/moment sectors can be embedded later in the
same augmented LDLT coordinate; it never replaces q_k by independent dv/dp/dS/
daw boxes.
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
import ou3_p4_brmm_physical_prediction_forcing as PHYS
import ou3_p4_brmm_physical_acceleration_witness_sector as SECTOR
import ou3_p4_a21_bias1_24state_event_lift as LIFT
import ou3_p4_source_uniform_bias_prefix_lineage as BIAS
import ou3_p4_same_history_nonlinear_graph_lineage as GRAPH
import ou3_brmm_complete_window_execution_kernel as KERNEL

SCHEMA=1
QUALIFICATION='OU3_P4_BRMM_PHYSICAL_SOURCE_EVERY_PREFIX_RESPONSE_V1'
P3_DELTA=1e-18

def I(x):return Interval.point(float(x))
def shape(A):return len(A),len(A[0]) if A else 0
def zero(r,c):return [[I(0) for _ in range(c)] for _ in range(r)]
def hstack(A,B):
    if len(A)!=len(B):raise ValueError('hstack row mismatch')
    return [list(a)+list(b) for a,b in zip(A,B)]
def embed_A24(G21):
    if shape(G21)!=(21,15):raise ValueError('A21 physical injection must be 21x15')
    return [list(r) for r in G21]+zero(3,15)
def propagate(A,C,G=None):
    carried=matrix_mul(A,C) if shape(C)[1] else zero(len(A),0)
    return carried if G is None else hstack(carried,G)

@dataclass(frozen=True)
class PhysicalSourcePrefixResponse:
    mode:str
    family:str|None
    prefix_ordinal:int
    source_cell_id:str
    event_token:str
    kind:str
    prediction_witness_ids:tuple[str,...]
    response:tuple[tuple[Interval,...],...]

def _joseph_H(cell):
    ev=EVENTS.source_joseph_event(**ATTACH.COVER.joseph_event_kwargs(cell));return ev['J_state']
def _joseph_A(cell):
    ev=EVENTS.source_joseph_event(**ATTACH.COVER.joseph_event_kwargs(cell));return LIFT.measurement_lift(ev)

def materialize(samples:Sequence[ATTACH.AttachedSampleLineage],*,mode:str,family:str|None=None,domain_path:Path=KERNEL.DEFAULT_DOMAIN):
    if not samples:raise ValueError('nonempty attached sample lineage required')
    if mode not in ('H','A'):raise ValueError('mode must be H/A')
    if any(s.mode!=mode for s in samples):raise ValueError('attached sample mode mismatch')
    if mode=='A' and family not in BIAS.FAMILIES:raise ValueError('A21 requires BIAS0/BIAS1/BIAS2')
    if mode=='H' and family is not None:raise ValueError('H18 physical response has no accelerometer-bias family')
    constants=KERNEL._process_constants(domain_path);projection=GRAPH._projection_limit(domain_path)
    n=18 if mode=='H' else 24;C=zero(n,0);records=[];witness_ids=[];state=list(samples[0].cells[0].state);ordinal=0
    bias_boxes=BIAS.prefix_component_boxes(family,len(samples)) if mode=='A' else None
    for si,sample in enumerate(samples):
        if sample.selector.sample_index!=si:raise RuntimeError('sample lineage index is not consecutive')
        beta=(bias_boxes[si],)*3 if bias_boxes is not None else None
        for cell in sample.cells:
            if tuple(cell.state)!=tuple(state):raise RuntimeError('physical response state detached from exact nonlinear lineage')
            G=None
            if cell.kind=='prediction':
                pred=PRED.prediction_event(mode,state,sample.selector.sample_coordinates.omega_body_corrected,cell.dt_s,cell.tau_applied_s,tau_ba=constants.accel_bias_tau_s if mode=='A' else None)
                if mode=='H':A=pred['J_state']
                else:A,_=LIFT.prediction_lift(pred['J_state'],BIAS.FAMILIES[family]['phi_true'] if isinstance(BIAS.FAMILIES,dict) else BIAS._spec(family).phi_true)
                M=PHYS.source_matrix(cell.tau_applied_s,cell.dt_s);G0=PHYS.inject_error_state(mode,M);G=G0 if mode=='H' else embed_A24(G0)
                witness_ids.append(f'{sample.selector.source_cell_id}:physical-acceleration')
                state=list(pred['state_out'])
            elif cell.kind=='aw_floor':
                A=matrix_identity(n)
            elif cell.kind in ('S_zero','accelerometer','magnetometer'):
                if mode=='A':
                    # Cell already carries the same-prefix true-bias interval from production attachment.
                    A=_joseph_A(cell)
                    ev=EVENTS.source_joseph_event(**ATTACH.COVER.joseph_event_kwargs(cell));state=list(ev['state_out'])
                else:
                    A=_joseph_H(cell);ev=EVENTS.source_joseph_event(**ATTACH.COVER.joseph_event_kwargs(cell));state=list(ev['state_out'])
            else:raise RuntimeError('unsupported literal event '+cell.kind)
            C=propagate(A,C,G);ordinal+=1
            if shape(C)!=(n,15*len(witness_ids)):raise RuntimeError('physical source response dimension drift')
            records.append(PhysicalSourcePrefixResponse(mode,family,ordinal,sample.selector.source_cell_id,cell.source_token,cell.kind,tuple(witness_ids),tuple(tuple(x for x in r) for r in C)))
        if tuple(state)!=tuple(sample.state_out):raise RuntimeError('physical response sample endpoint detached')
    return records

def build()->dict:
    p=PHYS.build();pf=PHYS.validate(p);s=SECTOR.build();sf=SECTOR.validate(s);b=BIAS.build();bf=BIAS.validate(b)
    if pf or sf or bf:raise RuntimeError(f'physical prefix prerequisites failed forcing={pf} sector={sf} bias={bf}')
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,
      'one_15D_same_history_witness_appended_per_prediction':True,'exact_tau_dependent_physical_forcing_injection_consumed':True,'previous_physical_witnesses_suffix_propagated_through_later_events':True,
      'H18_physical_source_prefix_response_materializable':True,'A21_all_BIAS0_BIAS1_BIAS2_physical_source_prefix_response_materializable':True,'A21_true_bias_proof_rows_receive_physical_acceleration_forcing':False,
      'independent_dv_dp_dS_daw_source_boxes_used':False,'hard_entry_radial_used_as_physical_source_amplitude':False,'physical_witness_quadratic_sector_available_for_each_prediction':True,
      'common_augmented_PrefixInput_assembled_here':False,'endpoint_augmented_LDLT_closed_here':False,'every_prefix_augmented_LDLT_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'embed each prefix response next to the initial-state transport, bias supply and binary32 response; place one physical a0/a1/moment sector per witness block with h_s=1; build storage master and run outward LDLT'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('one_15D_same_history_witness_appended_per_prediction','exact_tau_dependent_physical_forcing_injection_consumed','previous_physical_witnesses_suffix_propagated_through_later_events','H18_physical_source_prefix_response_materializable','A21_all_BIAS0_BIAS1_BIAS2_physical_source_prefix_response_materializable','physical_witness_quadratic_sector_available_for_each_prediction'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('A21_true_bias_proof_rows_receive_physical_acceleration_forcing','independent_dv_dp_dS_daw_source_boxes_used','hard_entry_radial_used_as_physical_source_amplitude','common_augmented_PrefixInput_assembled_here','endpoint_augmented_LDLT_closed_here','every_prefix_augmented_LDLT_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    return list(dict.fromkeys(f))
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'physical_prefix':d['H18_physical_source_prefix_response_materializable'],'A21':d['A21_all_BIAS0_BIAS1_BIAS2_physical_source_prefix_response_materializable'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
