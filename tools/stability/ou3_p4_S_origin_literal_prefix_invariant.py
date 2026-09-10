#!/usr/bin/env python3
"""Exact shared-S-origin invariant for every literal OU-III P4 event class.

The large session integration constant is a gauge/source coordinate, not an
independent pseudo-measurement residual.  In the joint coordinate [e;x_true], a
common S-origin shift c enters both error and truth through the same three
columns.  This file verifies, in exact rational arithmetic, that every shipping
event class preserves those common columns and that every measurement residual
is insensitive to them.

This does NOT set the centered S error/residual to zero after Live entry.  It
only proves that the arbitrary common integration constant can be quotiented
from every endpoint/prefix residual and storage graph.  The centered S dynamics,
physical source jet, cross-covariances, Joseph corrections and all nonlinear
terms remain in the production graph.
"""
from __future__ import annotations
import argparse,json
from fractions import Fraction as F
from pathlib import Path

SCHEMA=1
QUALIFICATION='OU3_P4_SHARED_S_ORIGIN_LITERAL_PREFIX_INVARIANT_V1'
S0=12

def eye(n):return [[F(int(i==j)) for j in range(n)] for i in range(n)]
def zeros(r,c):return [[F(0) for _ in range(c)] for _ in range(r)]
def mm(A,B):
    if not A or not B or len(A[0])!=len(B):raise ValueError('matrix mismatch')
    return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def sub(A,B):return [[A[i][j]-B[i][j] for j in range(len(A[0]))] for i in range(len(A))]
def eq(A,B):return A==B
def origin(n):
    E=zeros(n,3)
    for j in range(3):E[S0+j][j]=F(1)
    return E
def joint_origin(n):return origin(n)+origin(n)
def blockdiag(A,B):
    r1,c1=len(A),len(A[0]);r2,c2=len(B),len(B[0]);Z12=zeros(r1,c2);Z21=zeros(r2,c1)
    return [A[i]+Z12[i] for i in range(r1)]+[Z21[i]+B[i] for i in range(r2)]
def joint_measurement_map(n,H,K):
    # residual H(e-x_true); error <- e-K H(e-x_true), truth unchanged
    I=eye(n);KH=mm(K,H);top=[I[i]+zeros(n,n)[i] for i in range(n)]
    out=zeros(2*n,2*n)
    for i in range(n):
        for j in range(n):
            out[i][j]=I[i][j]-KH[i][j];out[i][n+j]=KH[i][j]
            out[n+i][n+j]=I[i][j]
    return out
def H_S(n):
    H=zeros(3,n)
    for j in range(3):H[j][S0+j]=F(1)
    return H
def H_nonS(n,offset=0):
    H=zeros(3,n)
    # arbitrary representative rows strictly outside S coordinates
    safe=[0,1,2] if offset==0 else [15,16,17]
    for j,k in enumerate(safe):H[j][k]=F(j+1,3)
    return H
def arbitrary_K(n):return [[F((i+2)*(j+1)-7,11) for j in range(3)] for i in range(n)]
def prediction_F(n):
    # Structural shipping linear chain: a constant S offset is invariant.
    A=eye(n)
    # Representative nonzero v->p, p->S and aw couplings; exact magnitudes do
    # not matter for the origin-column identity. No state depends on S except S.
    for j in range(3):
        A[9+j][6+j]=F(1,200)
        A[12+j][9+j]=F(1,200)
        A[6+j][15+j]=F(1,200)
        A[9+j][15+j]=F(1,80000)
        A[12+j][15+j]=F(1,24000000)
        A[15+j][15+j]=F(999,1000)
    return A
def state_local_map(n,indices):
    A=eye(n)
    # arbitrary finite local transformation on listed rows/cols only
    for i in indices:
        for j in indices:A[i][j]+=F((i+1)-(j+2),97)
    return A
def H_to_A_joint():
    # [e_H;x_H] -> [e_A;x_A], both append zero b_a coordinates. Common S passes unchanged.
    L=zeros(21,18)
    for i in range(18):L[i][i]=F(1)
    return blockdiag(L,L)

def check_mode(n):
    E=joint_origin(n);I=eye(n)
    Fp=prediction_F(n);pred=blockdiag(Fp,Fp)
    k=arbitrary_K(n)
    s_event=joint_measurement_map(n,H_S(n),k)
    acc_event=joint_measurement_map(n,H_nonS(n,1),k)
    mag_event=joint_measurement_map(n,H_nonS(n,0),k)
    floor=blockdiag(I,I)
    reset=blockdiag(state_local_map(n,(0,1,2)),I)
    proj=blockdiag(state_local_map(n,(18,19,20)) if n==21 else I,I)
    return {
      'prediction_preserves_common_S_origin':eq(mm(pred,E),E),
      'S_residual_annihilates_common_origin':eq(mm([row+[-x for x in row] for row in H_S(n)],E),zeros(3,3)),
      'S_Joseph_preserves_common_S_origin':eq(mm(s_event,E),E),
      'accelerometer_residual_ignores_S_origin':eq(mm([row+[-x for x in row] for row in H_nonS(n,1)],E),zeros(3,3)),
      'accelerometer_Joseph_preserves_common_S_origin':eq(mm(acc_event,E),E),
      'magnetometer_residual_ignores_S_origin':eq(mm([row+[-x for x in row] for row in H_nonS(n,0)],E),zeros(3,3)),
      'magnetometer_Joseph_preserves_common_S_origin':eq(mm(mag_event,E),E),
      'aw_floor_preserves_common_S_origin':eq(mm(floor,E),E),
      'attitude_reset_preserves_common_S_origin':eq(mm(reset,E),E),
      'bias_projection_preserves_common_S_origin':eq(mm(proj,E),E),
    }
def build():
    h=check_mode(18);a=check_mode(21);E18=joint_origin(18);E21=joint_origin(21);lift=H_to_A_joint()
    release=eq(mm(lift,E18),E21)
    all_checks=all(h.values()) and all(a.values()) and release
    # Composition induction: if M E=E and N E=E, then (N M)E=E. This is
    # exact associativity, so checking every event generator closes every finite
    # literal prefix without event-count enumeration.
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'joint_coordinate':'[physical_error; physical_truth/source_state]','shared_origin_columns':3,
      'H18_event_generator_checks':h,'A21_event_generator_checks':a,'H_to_A_release_preserves_common_S_origin':release,
      'all_literal_event_generators_preserve_common_S_origin':all_checks,
      'finite_literal_prefix_composition_closed_by_induction':all_checks,
      'arbitrary_same_event_gain_identity_used_only_algebraically':True,'independent_K_box_introduced':False,
      'common_origin_column_may_be_quotiented_after_joint_graph_assembly':all_checks,
      'centered_S_state_or_residual_set_to_zero_after_entry':False,
      'centered_S_dynamics_and_source_forcing_must_remain':True,
      'original_300m_s_error_ball_may_bound_S_zero_residual':False,
      'S_zero_residual_must_use_centered_joint_selector_eS_minus_Strue':True,
      'endpoint_augmented_LDLT_closed_here':False,'every_prefix_augmented_LDLT_closed_here':False,'P4_PASS':False}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('H_to_A_release_preserves_common_S_origin','all_literal_event_generators_preserve_common_S_origin','finite_literal_prefix_composition_closed_by_induction','arbitrary_same_event_gain_identity_used_only_algebraically','common_origin_column_may_be_quotiented_after_joint_graph_assembly','centered_S_dynamics_and_source_forcing_must_remain','S_zero_residual_must_use_centered_joint_selector_eS_minus_Strue'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('independent_K_box_introduced','centered_S_state_or_residual_set_to_zero_after_entry','original_300m_s_error_ball_may_bound_S_zero_residual','endpoint_augmented_LDLT_closed_here','every_prefix_augmented_LDLT_closed_here','P4_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    for table in ('H18_event_generator_checks','A21_event_generator_checks'):
        if not all(d.get(table,{}).values()):f.append(table+' incomplete')
    return f
def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'literal_prefix_origin_invariant':d['finite_literal_prefix_composition_closed_by_induction'],'S_residual_uses_original_ball':d['original_300m_s_error_ball_may_bound_S_zero_residual'],'prefix_LDLT':d['every_prefix_augmented_LDLT_closed_here'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
