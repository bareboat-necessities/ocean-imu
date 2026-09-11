#!/usr/bin/env python3
"""Scalable backward Schur operator for H18 source/fp endpoint bounds.

A complete 3 s word contains hundreds of physical source blocks and thousands
of binary32 event blocks. A monolithic dense augmented master is unnecessary.
This module implements the equivalent sequential quadratic elimination while
retaining same-history physical-acceleration continuity.

Persistent proof coordinate:

    x = [ e_H18(18), a_phys(3), h_s(1) ] in R^22.

At a prediction the local coordinates are

    u = [ a_next(3), J0(3), J1(3), J2(3), n_fp(18) ] in R^30.

The local admissible physical sectors use convention z'Pi z >= 0.  Therefore,
to prove the positive-form inequality

    Q_prev + supply - T' Q_next T >= 0

on the source graph, the S-procedure sufficient condition is

    Q_prev + supply - T'Q_next T - sum(lambda_j Pi_j) >= 0,
    lambda_j >= 0.

This sign is the positive-form equivalent of the repository's production
convention L+sum(lambda Pi)<0.  The sign is explicit here to prevent accidental
source-set reversal.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from ou3_interval import Interval,matrix_identity,matrix_mul,matrix_sub,matrix_transpose,symmetric_positive_definite_ldlt
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan,matrix_symmetric_hull,shape
import ou3_p4_brmm_physical_prediction_forcing as PHYS
import ou3_p4_brmm_physical_acceleration_witness_sector as SRC

SCHEMA=2
QUALIFICATION='OU3_P4_H18_SEQUENTIAL_SOURCE_FP_SCHUR_V2'
P3_DELTA=1e-18
NX=18;NA=3;NP=22;NFP=18;NLOCAL_PHYS=12

def I(x):return Interval.point(float(x))
def zero(r,c):return [[I(0) for _ in range(c)] for _ in range(r)]
def eye(n):return matrix_identity(n)
def hstack(A,B):
    if len(A)!=len(B):raise ValueError('hstack row mismatch')
    return [list(a)+list(b) for a,b in zip(A,B)]
def add(A,B):
    if shape(A)!=shape(B):raise ValueError('matrix add shape mismatch')
    return [[A[i][j]+B[i][j] for j in range(len(A[0]))] for i in range(len(A))]
def scale(A,s):
    q=I(s);return [[q*x for x in row] for row in A]
def _fp_injection(delta):
    G=zero(NP,NFP);d=I(delta)
    for i in range(NX):G[i][i]=d
    return G
def prediction_transition(Ae,tau,h,delta_fp):
    if shape(Ae)!=(NX,NX):raise ValueError('H18 event Jacobian required')
    Fq=PHYS.inject_error_state('H',PHYS.source_matrix(tau,h))
    Tx=zero(NP,NP);Tu=zero(NP,NLOCAL_PHYS+NFP)
    for i in range(NX):
        for j in range(NX):Tx[i][j]=Ae[i][j]
        for j in range(3):Tx[i][NX+j]=Fq[i][j]
        for j in range(12):Tu[i][j]=Fq[i][3+j]
    for i in range(3):Tu[NX+i][i]=I(1)
    Tx[NP-1][NP-1]=I(1)
    Gfp=_fp_injection(delta_fp)
    for i in range(NP):
        for j in range(NFP):Tu[i][NLOCAL_PHYS+j]=Gfp[i][j]
    return hstack(Tx,Tu)
def ordinary_transition(Ae,delta_fp):
    if shape(Ae)!=(NX,NX):raise ValueError('H18 event Jacobian required')
    Tx=zero(NP,NP)
    for i in range(NX):
        for j in range(NX):Tx[i][j]=Ae[i][j]
    for i in range(3):Tx[NX+i][NX+i]=I(1)
    Tx[NP-1][NP-1]=I(1)
    return hstack(Tx,_fp_injection(delta_fp))
def prediction_source_sectors(total_dim,h_s_index=NP-1,h=.005):
    if total_dim!=NP+NLOCAL_PHYS+NFP:raise ValueError('prediction local dimension mismatch')
    if h<=0:raise ValueError('positive sample time required')
    n=total_dim;A0=zero(3,n);A1=zero(3,n);M=zero(9,n);Hs=zero(1,n);Hs[0][h_s_index]=I(1)
    for axis in range(3):
        A0[axis][NX+axis]=I(1)
        A1[axis][NP+axis]=I(1)
        M[3*axis+0][NP+3+axis]=I(1/h)
        M[3*axis+1][NP+6+axis]=I(1/(h*h))
        M[3*axis+2][NP+9+axis]=I(1/(h*h*h))
    Amax=float(SRC.MOM.build()['A_max_mps2'])
    return SRC.physical_witness_sectors(A0,A1,M,Hs,Amax)
def _quadratic_pullback(Q,T):return matrix_mul(matrix_mul(matrix_transpose(T),Q),T)
def _supply(total_dim,gamma_s,gamma_n,prediction):
    S=zero(total_dim,total_dim)
    if prediction:S[NP-1][NP-1]=I(gamma_s)
    start=NP+(NLOCAL_PHYS if prediction else 0)
    for i in range(NFP):S[start+i][start+i]=I(gamma_n)
    return S
def local_master(Qnext,T,*,gamma_s,gamma_n,sectors=(),multipliers=(),prediction=False):
    """Return K=supply-T'Q+T-sum(lambda Pi), excluding Q_prev block."""
    total=shape(T)[1];K=matrix_sub(_supply(total,gamma_s,gamma_n,prediction),_quadratic_pullback(Qnext,T))
    if len(sectors)!=len(multipliers):raise ValueError('sector/multiplier mismatch')
    for P,l in zip(sectors,multipliers):
        if float(l)<0:raise ValueError('S-procedure multiplier must be nonnegative')
        K=matrix_sub(K,scale(P,float(l)))
    return matrix_symmetric_hull(K)
def schur_required_state_form(K):
    """Validated Schur elimination of local inputs.

    Returns R such that adding R to the persistent x-x block gives a PSD local
    positive-form master.  K_uu must itself be rigorously positive definite.
    """
    n=shape(K)[0]
    if shape(K)!=(n,n) or n<=NP:raise ValueError('local master must include eliminable inputs')
    Kxx=[row[:NP] for row in K[:NP]];Kxu=[row[NP:] for row in K[:NP]];Kux=[row[:NP] for row in K[NP:]];Kuu=[row[NP:] for row in K[NP:]]
    ok,piv=symmetric_positive_definite_ldlt(matrix_symmetric_hull(Kuu))
    if not ok:raise ValueError('local eliminable-input block is not rigorously positive definite')
    inv=matrix_inverse_gauss_jordan(Kuu);cross=matrix_mul(matrix_mul(Kxu,inv),Kux)
    R=matrix_symmetric_hull(matrix_sub(cross,Kxx))
    return R,Kuu,[float(x.lo) for x in piv]
def _smoke():
    Q=zero(NP,NP)
    for i in range(NX):Q[i][i]=I(1)
    T=ordinary_transition(eye(NX),1e-6);K=local_master(Q,T,gamma_s=0.0,gamma_n=1.0,prediction=False);R,Kuu,piv=schur_required_state_form(K)
    finite=all(x.lo==x.lo and x.hi==x.hi for row in R for x in row)
    # Sign mutation: a valid ball sector Pi>=0 must be subtracted in positive form.
    P=zero(NP+NLOCAL_PHYS+NFP,NP+NLOCAL_PHYS+NFP);P[0][0]=I(1)
    Tpred=prediction_transition(eye(NX),I(2),I(.005),1e-6)
    K0=local_master(zero(NP,NP),Tpred,gamma_s=1,gamma_n=1,sectors=[P],multipliers=[1],prediction=True)
    sign_ok=K0[0][0].hi<0
    return {'ordinary_local_dimension':shape(K)[0],'required_state_shape':shape(R),'Kuu_shape':shape(Kuu),'Kuu_pivot_lower_min':min(piv),'finite':finite,'positive_form_sector_sign_is_minus_lambda_Pi':sign_ok}
def build():
    p=PHYS.build();pf=PHYS.validate(p);s=SRC.build();sf=SRC.validate(s)
    if pf or sf:raise RuntimeError(f'sequential Schur prerequisites failed physical={pf} sector={sf}')
    sm=_smoke()
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,
      'persistent_coordinate':['e_H18_18','a_phys_3','h_source_1'],'persistent_dimension':NP,'prediction_local_source_dimension':NLOCAL_PHYS,'event_local_binary32_dimension':NFP,
      'same_history_a_endpoint_carried_across_predictions':True,'physical_source_scale_independent_of_entry_radial':True,'local_source_and_fp_inputs_eliminated_immediately':True,'carried_dimension_independent_of_word_length':True,
      'prediction_exact_BRMM_forcing_map_consumed':True,'prediction_source_a0_a1_moment_sectors_consumed':True,'S_procedure_positive_form_sign':'-lambda*Pi for Pi>=0','S_procedure_sign_matches_production_joint_sector_convention':True,
      'local_input_block_positive_definiteness_checked_by_outward_LDLT':True,'outward_interval_Schur_complement_available':True,'dense_whole_word_input_master_required':False,
      'smoke':sm,'source_uniform_multiplier_search_closed_here':False,'complete_word_backward_recursion_closed_here':False,'endpoint_contraction_closed_here':False,'every_prefix_bound_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'run the backward operator over every literal event of the canonical post-prediction word, search nonnegative source-sector multipliers and gamma_s/gamma_n with rigorous K_uu positivity, then compare the resulting entrance form against rho*M_entry and certify every saved prefix form'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('same_history_a_endpoint_carried_across_predictions','physical_source_scale_independent_of_entry_radial','local_source_and_fp_inputs_eliminated_immediately','carried_dimension_independent_of_word_length','prediction_exact_BRMM_forcing_map_consumed','prediction_source_a0_a1_moment_sectors_consumed','S_procedure_sign_matches_production_joint_sector_convention','local_input_block_positive_definiteness_checked_by_outward_LDLT','outward_interval_Schur_complement_available'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('dense_whole_word_input_master_required','source_uniform_multiplier_search_closed_here','complete_word_backward_recursion_closed_here','endpoint_contraction_closed_here','every_prefix_bound_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('persistent_dimension')!=22:f.append('persistent dimension changed')
    sm=d.get('smoke',{})
    if tuple(sm.get('required_state_shape',()))!=(22,22):f.append('Schur smoke state shape changed')
    if sm.get('finite') is not True:f.append('Schur smoke nonfinite')
    if sm.get('positive_form_sector_sign_is_minus_lambda_Pi') is not True:f.append('S-procedure sign smoke failed')
    if not(float(sm.get('Kuu_pivot_lower_min',-1))>0):f.append('Schur smoke Kuu not positive')
    return f
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'persistent_dim':d['persistent_dimension'],'Schur':d['outward_interval_Schur_complement_available'],'sector_sign':d['smoke']['positive_form_sector_sign_is_minus_lambda_Pi'],'complete_word':d['complete_word_backward_recursion_closed_here'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
