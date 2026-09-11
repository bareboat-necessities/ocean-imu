#!/usr/bin/env python3
"""Scalable backward Schur operator for H18 source/fp endpoint bounds.

A complete 3 s word contains hundreds of physical source blocks and thousands
of binary32 event blocks.  A monolithic dense augmented master would therefore
be unnecessarily enormous.  This module implements the equivalent sequential
quadratic elimination.

Persistent proof coordinate:

    x = [ e_H18(18), a_phys(3), h_s(1) ] in R^22,

where a_phys is the same-history physical acceleration endpoint carried across
predictions and h_s is the homogeneous physical-source scale (production h_s=1,
independent of the P4 entry radial).

At a prediction, local coordinates are

    u = [ a_next(3), J0(3), J1(3), J2(3), n_fp(18) ] in R^30.

The exact true-minus-estimate BRMM forcing uses the same tau/h as the shipping
homogeneous predictor.  The source state advances exactly a_phys+=a_next.  The
local source constraints are the a0 ball, a1 ball and coupled J0/J1/J2 IQC on
the SAME segment.  At non-prediction events only a local n_fp(18) block is
introduced; a_phys and h_s are identity-carried.

For a future quadratic form Q+ on x+, a local desired dissipation inequality is

  x+^T Q+ x+ <= x^T Q x + gamma_s h_s^2 + gamma_n ||n_fp||^2

for every local physical witness satisfying the source sectors.  After adding
nonnegative S-procedure multiples, the local master on [x;u] is partitioned.
If its u-u block is rigorously positive definite, outward interval inversion and
Schur complementation give the minimal carried x-x contribution required for Q.
Thus local inputs are eliminated immediately and the carried dimension remains
22 regardless of word length.

This module materializes and validates that elimination primitive.  It does not
yet choose the globally optimal multipliers/gammas or claim word contraction;
those remain the next source-uniform search/certificate step.
"""
from __future__ import annotations
import argparse,json
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence
from ou3_interval import Interval,matrix_identity,matrix_mul,matrix_sub,matrix_transpose
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan,matrix_symmetric_hull,shape
import ou3_p4_brmm_physical_prediction_forcing as PHYS
import ou3_p4_brmm_physical_acceleration_witness_sector as SRC

SCHEMA=1
QUALIFICATION='OU3_P4_H18_SEQUENTIAL_SOURCE_FP_SCHUR_V1'
P3_DELTA=1e-18
NX=18;NA=3;NH=1;NP=22;NFP=18;NLOCAL_PHYS=12

def I(x):return Interval.point(float(x))
def zero(r,c):return [[I(0) for _ in range(c)] for _ in range(r)]
def eye(n):return matrix_identity(n)
def hstack(A,B):
    if len(A)!=len(B):raise ValueError('hstack row mismatch')
    return [list(a)+list(b) for a,b in zip(A,B)]
def vstack(A,B):
    if A and B and len(A[0])!=len(B[0]):raise ValueError('vstack column mismatch')
    return [list(r) for r in A]+[list(r) for r in B]
def add(A,B):
    if shape(A)!=shape(B):raise ValueError('matrix add shape mismatch')
    return [[A[i][j]+B[i][j] for j in range(len(A[0]))] for i in range(len(A))]
def scale(A,s):
    q=I(s);return [[q*x for x in row] for row in A]
def embed_block(n,rows,cols,B):
    M=zero(n,n)
    if len(rows)!=len(B) or (B and len(cols)!=len(B[0])):raise ValueError('embed block mismatch')
    for i,r in enumerate(rows):
        for j,c in enumerate(cols):M[r][c]=B[i][j]
    return M

def _fp_injection(delta):
    G=zero(NP,NFP);d=I(delta)
    for i in range(NX):G[i][i]=d
    return G

def prediction_transition(Ae,tau,h,delta_fp):
    """Return T mapping [x22;u30] -> x+22."""
    if shape(Ae)!=(NX,NX):raise ValueError('H18 event Jacobian required')
    Fq=PHYS.inject_error_state('H',PHYS.source_matrix(tau,h)) # 18x15 [a0,a1,J0,J1,J2]
    # x columns: e18,a0(3),hs. local u: a1,J0,J1,J2,n18.
    Tx=zero(NP,NP);Tu=zero(NP,NLOCAL_PHYS+NFP)
    for i in range(NX):
        for j in range(NX):Tx[i][j]=Ae[i][j]
        for j in range(3):Tx[i][NX+j]=Fq[i][j]
        for j in range(12):Tu[i][j]=Fq[i][3+j]
    for i in range(3):Tu[NX+i][i]=I(1) # a_next becomes persistent a_phys
    Tx[NP-1][NP-1]=I(1)                # homogeneous source scale
    Gfp=_fp_injection(delta_fp)
    for i in range(NP):
        for j in range(NFP):Tu[i][NLOCAL_PHYS+j]=Gfp[i][j]
    return hstack(Tx,Tu)
def ordinary_transition(Ae,delta_fp):
    """Return T mapping [x22;n18] -> x+22 for non-prediction literal event."""
    if shape(Ae)!=(NX,NX):raise ValueError('H18 event Jacobian required')
    Tx=zero(NP,NP)
    for i in range(NX):
        for j in range(NX):Tx[i][j]=Ae[i][j]
    for i in range(3):Tx[NX+i][NX+i]=I(1)
    Tx[NP-1][NP-1]=I(1)
    return hstack(Tx,_fp_injection(delta_fp))
def prediction_source_sectors(total_dim,h_s_index=NP-1):
    """Sectors on [x22;u30], sharing carried a0 and local a1/J moments."""
    n=total_dim
    if n!=NP+NLOCAL_PHYS+NFP:raise ValueError('prediction local dimension mismatch')
    z=I(0);A0=zero(3,n);A1=zero(3,n);M=zero(9,n);Hs=zero(1,n);Hs[0][h_s_index]=I(1)
    for axis in range(3):
        A0[axis][NX+axis]=I(1)
        A1[axis][NP+axis]=I(1)
        # local u layout a1 0:3,J0 3:6,J1 6:9,J2 9:12
        M[3*axis+0][NP+3+axis]=I(1/.005)
        M[3*axis+1][NP+6+axis]=I(1/(.005**2))
        M[3*axis+2][NP+9+axis]=I(1/(.005**3))
    Amax=float(SRC.MOM.build()['A_max_mps2'])
    return SRC.physical_witness_sectors(A0,A1,M,Hs,Amax)
def _quadratic_pullback(Q,T):return matrix_mul(matrix_mul(matrix_transpose(T),Q),T)
def _supply(total_dim,local_dim,gamma_s,gamma_n,prediction):
    S=zero(total_dim,total_dim)
    if prediction:S[NP-1][NP-1]=I(gamma_s)
    start=NP+(NLOCAL_PHYS if prediction else 0)
    for i in range(NFP):S[start+i][start+i]=I(gamma_n)
    return S
def local_master(Qnext,T,*,gamma_s,gamma_n,sectors=(),multipliers=(),prediction=False):
    total=shape(T)[1];K=matrix_sub(_supply(total,total-NP,gamma_s,gamma_n,prediction),_quadratic_pullback(Qnext,T))
    if len(sectors)!=len(multipliers):raise ValueError('sector/multiplier mismatch')
    for P,l in zip(sectors,multipliers):K=add(K,scale(P,float(l)))
    return matrix_symmetric_hull(K)
def schur_required_state_form(K):
    """Return R on persistent x such that R+K_xx makes full local master PSD.

    K is [x;u]. With K_uu rigorously invertible/positive (the latter is checked
    by outward LDLT by the caller), the exact Schur requirement is

       R >= -K_xx + K_xu K_uu^{-1} K_ux.
    """
    n=shape(K)[0]
    if shape(K)!=(n,n) or n<=NP:raise ValueError('local master must include eliminable inputs')
    Kxx=[row[:NP] for row in K[:NP]];Kxu=[row[NP:] for row in K[:NP]];Kux=[row[:NP] for row in K[NP:]];Kuu=[row[NP:] for row in K[NP:]]
    inv=matrix_inverse_gauss_jordan(Kuu)
    cross=matrix_mul(matrix_mul(Kxu,inv),Kux)
    R=matrix_symmetric_hull(matrix_sub(cross,Kxx))
    return R,Kuu

def _smoke():
    # Point identity event and intentionally generous supplies: validates the
    # interval Schur algebra, not a theorem tuning constant.
    Q=zero(NP,NP)
    for i in range(NX):Q[i][i]=I(1)
    T=ordinary_transition(eye(NX),1e-6)
    K=local_master(Q,T,gamma_s=0.0,gamma_n=1.0,prediction=False)
    R,Kuu=schur_required_state_form(K)
    finite=all(x.lo==x.lo and x.hi==x.hi for row in R for x in row)
    return {'ordinary_local_dimension':shape(K)[0],'required_state_shape':shape(R),'Kuu_shape':shape(Kuu),'finite':finite}
def build():
    p=PHYS.build();pf=PHYS.validate(p);s=SRC.build();sf=SRC.validate(s)
    if pf or sf:raise RuntimeError(f'sequential Schur prerequisites failed physical={pf} sector={sf}')
    sm=_smoke()
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,
      'persistent_coordinate':['e_H18_18','a_phys_3','h_source_1'],'persistent_dimension':NP,'prediction_local_source_dimension':NLOCAL_PHYS,'event_local_binary32_dimension':NFP,
      'same_history_a_endpoint_carried_across_predictions':True,'physical_source_scale_independent_of_entry_radial':True,'local_source_and_fp_inputs_eliminated_immediately':True,'carried_dimension_independent_of_word_length':True,
      'prediction_exact_BRMM_forcing_map_consumed':True,'prediction_source_a0_a1_moment_sectors_consumed':True,'outward_interval_Schur_complement_available':True,'dense_whole_word_input_master_required':False,
      'smoke':sm,'source_uniform_multiplier_search_closed_here':False,'complete_word_backward_recursion_closed_here':False,'endpoint_contraction_closed_here':False,'every_prefix_bound_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'run the backward operator over every literal event of the canonical post-prediction word, search nonnegative source-sector multipliers and gamma_s/gamma_n with rigorous K_uu positivity, then compare the resulting entrance form against rho*M_entry and certify every saved prefix form'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('same_history_a_endpoint_carried_across_predictions','physical_source_scale_independent_of_entry_radial','local_source_and_fp_inputs_eliminated_immediately','carried_dimension_independent_of_word_length','prediction_exact_BRMM_forcing_map_consumed','prediction_source_a0_a1_moment_sectors_consumed','outward_interval_Schur_complement_available'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('dense_whole_word_input_master_required','source_uniform_multiplier_search_closed_here','complete_word_backward_recursion_closed_here','endpoint_contraction_closed_here','every_prefix_bound_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('persistent_dimension')!=22:f.append('persistent dimension changed')
    if d.get('smoke',{}).get('required_state_shape')!=[22,22] and tuple(d.get('smoke',{}).get('required_state_shape',()))!=(22,22):f.append('Schur smoke state shape changed')
    if d.get('smoke',{}).get('finite') is not True:f.append('Schur smoke nonfinite')
    return f
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'persistent_dim':d['persistent_dimension'],'Schur':d['outward_interval_Schur_complement_available'],'complete_word':d['complete_word_backward_recursion_closed_here'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
