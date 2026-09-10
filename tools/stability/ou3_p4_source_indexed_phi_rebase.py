#!/usr/bin/env python3
"""Exact source-indexed Phi-coordinate rebase between COMPLETE-BRMM samples.

The nonlinear measurement-linearizing coordinate

    Phi(z;s) = z + E_aw epsilon_aw(z;s)

is indexed by the physical/source geometry ``s`` through ``R_hat`` and
``f_hat``.  Consecutive provider samples generally have different source
geometry even when the physical error z is continuous.  A complete-word
transport must therefore contain an explicit coordinate-rebase event between
sample-local Riccati slices.

This is NOT a physical filter event and introduces no independent disturbance.
For fixed physical state z and consecutive same-history source coordinates
s_old,s_new,

    z+ = z,
    Phi(z;s_new) = I*Phi(z;s_old) + xi_rebase,
    xi_rebase = E_aw [epsilon_aw(z;s_new)-epsilon_aw(z;s_old)].

Thus the physical map is exactly identity, C=L=I, rho=0, and the entire rebase
is the endpoint-transport boundary shift already present in the exact
variation-of-constants identity.  It must be tied to the same provider
primitive transition; packetwise re-zeroing epsilon is forbidden.

This module provides the exact algebra/API.  It does not claim a source-uniform
numerical bound on the rebase defect and cannot promote P4 by itself.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ou3_interval import Interval, matrix_identity
import ou3_p4_complete_brmm_measurement_linearizing_aw_coordinate as AW
import ou3_p4_complete_word_endpoint_transport as ENDPOINT
import ou3_p4_brmm_primitive_prefix_binding as PRIMITIVE

SCHEMA=1
QUALIFICATION='OU3_P4_SOURCE_INDEXED_PHI_REBASE_V1'
P3_DELTA=1.0e-18


def I(x:float)->Interval:
    return Interval.point(float(x))


def zero_vec(n:int):
    return [[I(0.0)] for _ in range(n)]


def rebase_event(mode:str,epsilon_before,epsilon_after):
    n=18 if mode=='H' else 21 if mode=='A' else 0
    if not n:raise ValueError('mode must be H/A')
    C=matrix_identity(n);L=matrix_identity(n);rho=zero_vec(n)
    xi=AW.evaluate_full_shift_transport(C,rho,epsilon_before,epsilon_after)
    return {'kind':'source_coordinate_rebase','C':C,'L':L,'rho':[row[0] for row in rho],
            'epsilon_before':epsilon_before,'epsilon_after':epsilon_after,'xi':xi,
            'physical_state_map_identity':True,'C_equals_L_exactly':C==L,
            'physical_defect_exactly_zero':all(row[0]==I(0.0) for row in rho)}


def _embed(n:int,eps):
    out=zero_vec(n)
    for i in range(3):out[15+i]=[eps[i][0]]
    return out


def _smoke(mode='H'):
    n=18 if mode=='H' else 21
    e0=[[I(.1)],[I(-.2)],[I(.05)]]
    e1=[[I(.15)],[I(-.1)],[I(.02)]]
    ev=rebase_event(mode,e0,e1)
    expected=[[I(0.0)] for _ in range(n)]
    for i in range(3): expected[15+i]=[e1[i][0]-e0[i][0]]
    exact=ev['xi']==expected
    # Exercise the canonical endpoint decomposition itself.  With C=L=I and
    # rho=0, there is no interior event contribution; the complete effect is
    # exactly E1 eps1 - E0 eps0.
    E=[]
    for _ in range(2):
        M=[[I(0.0) for _ in range(3)] for _ in range(n)]
        for j in range(3):M[15+j][j]=I(1.0)
        E.append(M)
    d=ENDPOINT.endpoint_decomposition(
        [{'kind':'source_coordinate_rebase','C':ev['C'],'L':ev['L'],'rho':ev['rho']}],
        E,[[x[0] for x in e0],[x[0] for x in e1]])
    return {'xi_exact':exact,'no_interior_event':d['interior_event_indices']==[],
            'identity_residual_zero':all(x==I(0.0) for x in d['identity_residual'])}


def build():
    aw=AW.build();af=AW.validate(aw)
    ep=ENDPOINT.build();ef=ENDPOINT.validate(ep)
    pb=PRIMITIVE.build();pf=PRIMITIVE.validate(pb)
    if af or ef or pf:raise RuntimeError(f'Phi rebase prerequisites failed aw={af} endpoint={ef} primitive={pf}')
    H=_smoke('H');A=_smoke('A')
    closed=all(H.values()) and all(A.values())
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'P3_delta':P3_DELTA,'source_indexed_Phi_requires_cross_sample_rebase':True,
      'physical_state_unchanged_by_rebase':True,'physical_rebase_map':'I',
      'physical_rebase_defect':'0','Phi_rebase_defect':'E_aw*(epsilon_new-epsilon_old)',
      'C_equals_L_identity_for_rebase':True,'rebase_has_no_interior_epsilon_transport':True,
      'same_provider_transition_must_own_old_and_new_source_coordinates':True,
      'cross_sample_primitive_continuity_contract_consumed':True,
      'epsilon_packetwise_rezero_forbidden':True,'independent_rebase_disturbance_port_used':False,
      'source_coordinate_rebase_exact_algebra_closed':closed,'H18_smoke':H,'A21_smoke':A,
      'source_uniform_rebase_defect_bound_closed_here':False,
      'production_rebases_inserted_between_all_samples_here':False,
      'endpoint_augmented_LDLT_closed_here':False,'every_prefix_augmented_LDLT_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'insert one same-provider Phi rebase between consecutive synchronized sample slices, evaluate epsilon_old/new from their actual source geometry and shared physical boundary state, then include those nodes in every prefix_decomposition before augmented LDLT'}


def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('source_indexed_Phi_requires_cross_sample_rebase','physical_state_unchanged_by_rebase','C_equals_L_identity_for_rebase','rebase_has_no_interior_epsilon_transport','same_provider_transition_must_own_old_and_new_source_coordinates','cross_sample_primitive_continuity_contract_consumed','epsilon_packetwise_rezero_forbidden','source_coordinate_rebase_exact_algebra_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('independent_rebase_disturbance_port_used','source_uniform_rebase_defect_bound_closed_here','production_rebases_inserted_between_all_samples_here','endpoint_augmented_LDLT_closed_here','every_prefix_augmented_LDLT_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    for name in ('H18_smoke','A21_smoke'):
        for k,v in d.get(name,{}).items():
            if v is not True:f.append(name+' '+k+' not true')
    return f


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'rebase':d['source_coordinate_rebase_exact_algebra_closed'],'H':d['H18_smoke'],'A':d['A21_smoke'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))

if __name__=='__main__':raise SystemExit(main())
