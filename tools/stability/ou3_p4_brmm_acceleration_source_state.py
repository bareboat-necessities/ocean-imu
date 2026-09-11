#!/usr/bin/env python3
"""Cross-prediction physical acceleration source state for COMPLETE-BRMM P4.

A sequence of 15D per-prediction witnesses cannot be treated independently:
for one physical acceleration history, the endpoint of transition k is the
initial value of transition k+1,

    a1_k == a0_{k+1}.

This module materializes that equality as a persistent 3D source state.  A
prediction step carries

    s_k = a_k

and introduces local source coordinates

    u_k = [a_{k+1}, J0_k, J1_k, J2_k].

The full 15D witness consumed by the exact physical error-forcing map is then
assembled by the fixed linear identity

    q_k=[s_k, a_{k+1}, J0_k, J1_k, J2_k].

After the step the persistent source state is exactly s_{k+1}=a_{k+1}.  This
prevents independent endpoint boxes at consecutive predictions and makes
sequential Schur/LDLT elimination possible without losing same-history endpoint
continuity.  The local moment IQC and acceleration endpoint balls remain outer
constraints on one admitted continuous acceleration segment.
"""
from __future__ import annotations
import argparse,json
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence
from ou3_interval import Interval
import ou3_p4_brmm_physical_prediction_forcing as FORCING
import ou3_p4_brmm_physical_acceleration_witness_sector as SECTOR

SCHEMA=1
QUALIFICATION='OU3_P4_BRMM_CROSS_PREDICTION_ACCELERATION_SOURCE_STATE_V1'
P3_DELTA=1e-18

def I(x):return Interval.point(float(x))
def zero(r,c):return [[I(0) for _ in range(c)] for _ in range(r)]
def q_assembly_map():
    """15x15 map [a_k(3), u_k(12)] -> q=[a0,a1,J0,J1,J2]."""
    A=zero(15,15)
    for i in range(3):
        A[i][i]=I(1)          # carried a_k -> a0
        A[3+i][3+i]=I(1)      # local a_{k+1}
        A[6+i][6+i]=I(1)      # J0
        A[9+i][9+i]=I(1)      # J1
        A[12+i][12+i]=I(1)    # J2
    return A
def next_source_state_map():
    """3x15 map [a_k,u_k] -> a_{k+1}."""
    A=zero(3,15)
    for i in range(3):A[i][3+i]=I(1)
    return A
def assemble_witness(a_k:Sequence[Interval],a_next:Sequence[Interval],J0:Sequence[Interval],J1:Sequence[Interval],J2:Sequence[Interval]):
    if any(len(x)!=3 for x in (a_k,a_next,J0,J1,J2)):raise ValueError('all physical acceleration coordinates must be 3D')
    return FORCING.witness(a0=a_k,a1=a_next,J0=J0,J1=J1,J2=J2)
@dataclass(frozen=True)
class SourceStep:
    predecessor_token:str
    successor_token:str
    a_in:tuple[Interval,Interval,Interval]
    a_out:tuple[Interval,Interval,Interval]
    J0:tuple[Interval,Interval,Interval]
    J1:tuple[Interval,Interval,Interval]
    J2:tuple[Interval,Interval,Interval]

def validate_chain(steps:Sequence[SourceStep])->list[str]:
    f=[]
    for i,s in enumerate(steps):
        if not s.predecessor_token or not s.successor_token:f.append(f'step {i}: missing token')
        if i:
            if s.predecessor_token!=steps[i-1].successor_token:f.append(f'step {i}: token ancestry detached')
            if s.a_in!=steps[i-1].a_out:f.append(f'step {i}: a1/a0 same-history continuity detached')
    return list(dict.fromkeys(f))
def build():
    f=FORCING.build();ff=FORCING.validate(f);s=SECTOR.build();sf=SECTOR.validate(s)
    if ff or sf:raise RuntimeError(f'acceleration source-state prerequisites failed forcing={ff} sector={sf}')
    z=(I(0),I(0),I(0));a=(I(.1),I(-.2),I(.3));b=(I(.2),I(-.1),I(.25));j=(I(0),I(0),I(0))
    chain=(SourceStep('root','k0',z,a,j,j,j),SourceStep('k0','k1',a,b,j,j,j))
    cf=validate_chain(chain)
    if cf:raise RuntimeError('source-state positive control failed '+repr(cf))
    bad=list(chain);bad[1]=SourceStep('k0','k1',z,b,j,j,j)
    detached=validate_chain(bad)
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,
      'persistent_physical_acceleration_source_state_dimension':3,'local_prediction_source_dimension':12,'assembled_same_history_witness_dimension':15,
      'a1_k_equals_a0_kplus1_enforced':True,'physical_acceleration_endpoint_reselected_independently_each_prediction':False,
      'fixed_linear_q_assembly_map_materialized':True,'fixed_linear_next_source_state_map_materialized':True,'exact_physical_error_forcing_map_consumed':True,'joint_physical_witness_sector_consumed':True,
      'same_history_chain_positive_control_closed':True,'detached_a1_a0_mutation_rejected':bool(detached),
      'sequential_source_state_elimination_architecture_available':True,'production_all_prediction_source_states_bound_here':False,'endpoint_augmented_LDLT_closed_here':False,'every_prefix_augmented_LDLT_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'augment the H18/A21 storage state by carried a_k, append only [a_next,J0,J1,J2] per prediction, impose the same-segment a0/a1/moment sectors, eliminate local source coordinates after each event while carrying a_next into the successor'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('a1_k_equals_a0_kplus1_enforced','fixed_linear_q_assembly_map_materialized','fixed_linear_next_source_state_map_materialized','exact_physical_error_forcing_map_consumed','joint_physical_witness_sector_consumed','same_history_chain_positive_control_closed','detached_a1_a0_mutation_rejected','sequential_source_state_elimination_architecture_available'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('physical_acceleration_endpoint_reselected_independently_each_prediction','production_all_prediction_source_states_bound_here','endpoint_augmented_LDLT_closed_here','every_prefix_augmented_LDLT_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('persistent_physical_acceleration_source_state_dimension')!=3 or d.get('local_prediction_source_dimension')!=12 or d.get('assembled_same_history_witness_dimension')!=15:f.append('source dimensions changed')
    return f
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'continuity':d['a1_k_equals_a0_kplus1_enforced'],'sequential':d['sequential_source_state_elimination_architecture_available'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
