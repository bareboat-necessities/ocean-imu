#!/usr/bin/env python3
"""ISS supply-rate extension of the complete-BRMM augmented P4 master.

If z is the common augmented coordinate and L is the exact endpoint/prefix
energy matrix, source and arithmetic forcing coordinates are selected by maps
S_s z and S_n z.  The desired inequality

  Delta V <= -strict + gamma_s ||S_s z||^2 + gamma_n ||S_n z||^2

is certified by the same full interval LDLT after replacing

  L -> L - gamma_s S_s^T S_s - gamma_n S_n^T S_n.

Nonlinear chord/reset/projection graph IQCs are then added by the existing
nonnegative-multiplier S-procedure.  This keeps BIAS1 and binary32 effects as
explicit additive ISS channels; neither is allowed to erase homogeneous
contraction or masquerade as an exact-real state map.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval,matrix_add,matrix_mul,matrix_sub,matrix_transpose,symmetric_positive_definite_ldlt
from ou3_interval_linear_algebra import matrix_symmetric_hull
import ou3_p4_complete_brmm_joint_sector_master as JOINT
import ou3_p4_bias1_joint_iss_supply as BIAS1
import ou3_p4_projection_binary32_enclosure as PROJFP

SCHEMA=1
QUALIFICATION='OU3_P4_COMPLETE_BRMM_JOINT_ISS_AUGMENTED_MASTER_V1'

def _shape(A):return len(A),len(A[0]) if A else 0
def _gram(A):return matrix_mul(matrix_transpose(A),A)
def _scale(A,a):
    c=Interval.point(float(a));return [[c*x for x in row] for row in A]
def add_supply_penalties(master, source_map=None, gamma_s=0.0, fp_map=None, gamma_n=0.0):
    n,m=_shape(master)
    if n==0 or n!=m:raise ValueError('square master required')
    out=matrix_symmetric_hull(master)
    for name,A,g in (('source',source_map,gamma_s),('finite_precision',fp_map,gamma_n)):
        if A is None:
            if g!=0.0:raise ValueError(name+' gamma given without map')
            continue
        r,c=_shape(A)
        if r==0 or c!=n:raise ValueError(name+' map dimension mismatch')
        if not (math.isfinite(float(g)) and float(g)>=0):raise ValueError(name+' gamma invalid')
        out=matrix_symmetric_hull(matrix_sub(out,_scale(_gram(A),float(g))))
    return out

def certify_iss(master,sectors,multipliers,source_map=None,gamma_s=0.0,fp_map=None,gamma_n=0.0):
    L=add_supply_penalties(master,source_map,gamma_s,fp_map,gamma_n)
    test=JOINT.joint_sector_sprocedure_matrix(L,sectors,multipliers)
    neg=[[-x for x in row] for row in test]
    ok,p=symmetric_positive_definite_ldlt(matrix_symmetric_hull(neg))
    return bool(ok),[float(x.lo) for x in p]

def _selector(n,indices):
    A=[[Interval.point(0.0) for _ in range(n)] for _ in indices]
    for r,i in enumerate(indices):A[r][i]=Interval.point(1.0)
    return A

def build():
    b=BIAS1.build();bf=BIAS1.validate(b);p=PROJFP.build();pf=PROJFP.validate(p)
    if bf or pf:raise RuntimeError(f'supply prerequisites failed BIAS1={bf} FP={pf}')
    # Smoke: DeltaV=-||x||^2+2 x.s+2 x.n, with supply penalties.
    # Coordinate [x,s,n], no graph sectors needed by production callers; use
    # zero sector solely to exercise the common S-procedure API.
    I=Interval.point
    L=[[I(-1),I(1),I(1)],[I(1),I(0),I(0)],[I(1),I(0),I(0)]]
    Z=[[I(0) for _ in range(3)] for _ in range(3)]
    Sm=_selector(3,[1]);Nm=_selector(3,[2])
    ok,piv=certify_iss(L,[Z],[0.0],Sm,4.1,Nm,4.1)
    return {
      'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'source_supply_contract':b['qualification'],'finite_precision_supply_contract':p['qualification'],
      'same_augmented_coordinate_for_state_graph_source_and_roundoff':True,
      'source_supply_penalty_sign':'-gamma_s*S_s^T*S_s in negativity test',
      'finite_precision_supply_penalty_sign':'-gamma_n*S_n^T*S_n in negativity test',
      'full_interval_LDLT_reused':True,'nonnegative_graph_multipliers_reused':True,
      'source_forcing_not_counted_as_homogeneous_contraction':True,
      'roundoff_not_claimed_exact_real':True,'projection_binary32_channel_available':True,
      'smoke_ISS_LDLT_closed':ok,'smoke_pivot_lowers':piv,
      'production_endpoint_LDLT_closed_here':False,'production_every_prefix_LDLT_closed_here':False,'P4_promoted_here':False,
    }
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('same_augmented_coordinate_for_state_graph_source_and_roundoff','full_interval_LDLT_reused','nonnegative_graph_multipliers_reused','source_forcing_not_counted_as_homogeneous_contraction','roundoff_not_claimed_exact_real','projection_binary32_channel_available','smoke_ISS_LDLT_closed'):
      if d.get(k) is not True:f.append(k+' not true')
    for k in ('production_endpoint_LDLT_closed_here','production_every_prefix_LDLT_closed_here','P4_promoted_here'):
      if d.get(k) is not False:f.append(k+' not false')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps(d,sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
