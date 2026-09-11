#!/usr/bin/env python3
"""Sequential Schur operator for A21 with physical, bias and binary32 sources.

Persistent proof coordinate

  x = [e_A21(21), b_true(3), a_phys(3), h_phys(1), h_bias(1)] in R^29.

At a prediction, local coordinates are

  u = [a_next(3), J0(3), J1(3), J2(3), s_bias(6), n_fp(21)],

where s_bias=[w,m_tau] is the SAME six-dimensional supply used by the deployed
24-state bias lift. The physical acceleration endpoint is carried across
predictions, true bias is carried inside the 24-state event map, and both source
families have distinct homogeneous scale coordinates. Neither scale is the P4
entry radial coordinate.

Two bias effects must not be conflated. A prediction introduces the new driver
s_bias and therefore may consume gamma_bias D_b there. A literal prefix ending
before the next prediction can nevertheless depend on the already-persistent
b_true through active projection. Such a prefix is allowed one prefix-level
kappa_b D_b term on its terminal event; it is not charged once per event.
``bias_supply_on_this_event`` represents exactly that placement choice.

The local source constraints use z'Pi z >= 0 and are subtracted in the
positive-form S-procedure, matching the corrected H18 sequential operator.
This module only materializes the elimination primitive; it does not promote P4.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path

from ou3_interval import Interval,matrix_identity,matrix_mul,matrix_sub,matrix_transpose,symmetric_positive_definite_ldlt
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan,matrix_symmetric_hull,shape
import ou3_p4_brmm_physical_prediction_forcing as PHYS
import ou3_p4_brmm_physical_acceleration_witness_sector as SRC
import ou3_p4_bias_family_joint_iss_supply as BIAS

SCHEMA=2
QUALIFICATION='OU3_P4_A21_SEQUENTIAL_SOURCE_BIAS_FP_SCHUR_V2'
P3_DELTA=1e-18
NE=21;NBETA=3;NJOINT=24;NA=3;NHP=1;NHB=1;NP=29
NLOCAL_PHYS=12;NLOCAL_BIAS=6;NFP=21
OFF_BETA=21;OFF_A=24;IDX_HP=27;IDX_HB=28

def I(x):return Interval.point(float(x))
def zero(r,c):return [[I(0) for _ in range(c)] for _ in range(r)]
def eye(n):return matrix_identity(n)
def hstack(A,B):
    if len(A)!=len(B):raise ValueError('hstack row mismatch')
    return [list(a)+list(b) for a,b in zip(A,B)]
def scale(A,s):
    q=I(float(s));return [[q*x for x in row] for row in A]
def _fp_injection(delta):
    G=zero(NP,NFP);d=I(delta)
    for i in range(NE):G[i][i]=d
    return G

def prediction_transition(A24,Bbias,tau,h,delta_fp):
    if shape(A24)!=(NJOINT,NJOINT):raise ValueError('24x24 A21 joint event map required')
    if shape(Bbias)!=(NJOINT,NLOCAL_BIAS):raise ValueError('24x6 shared bias supply map required')
    Fq=PHYS.inject_error_state('A',PHYS.source_matrix(tau,h))
    Tx=zero(NP,NP);Tu=zero(NP,NLOCAL_PHYS+NLOCAL_BIAS+NFP)
    for i in range(NJOINT):
        for j in range(NJOINT):Tx[i][j]=A24[i][j]
    for i in range(NE):
        for j in range(3):Tx[i][OFF_A+j]=Fq[i][j]
        for j in range(NLOCAL_PHYS):Tu[i][j]=Fq[i][3+j]
    for i in range(NJOINT):
        for j in range(NLOCAL_BIAS):Tu[i][NLOCAL_PHYS+j]=Bbias[i][j]
    for i in range(3):Tu[OFF_A+i][i]=I(1)
    Tx[IDX_HP][IDX_HP]=I(1);Tx[IDX_HB][IDX_HB]=I(1)
    G=_fp_injection(delta_fp);start=NLOCAL_PHYS+NLOCAL_BIAS
    for i in range(NP):
        for j in range(NFP):Tu[i][start+j]=G[i][j]
    return hstack(Tx,Tu)

def ordinary_transition(A24,delta_fp):
    if shape(A24)!=(NJOINT,NJOINT):raise ValueError('24x24 A21 joint event map required')
    Tx=zero(NP,NP)
    for i in range(NJOINT):
        for j in range(NJOINT):Tx[i][j]=A24[i][j]
    for i in range(3):Tx[OFF_A+i][OFF_A+i]=I(1)
    Tx[IDX_HP][IDX_HP]=I(1);Tx[IDX_HB][IDX_HB]=I(1)
    return hstack(Tx,_fp_injection(delta_fp))

def _selector(total,indices):
    A=zero(len(tuple(indices)),total)
    for r,c in enumerate(tuple(indices)):A[r][c]=I(1)
    return A

def _ball_sector(vmap,hindex,bound):
    total=len(vmap[0]);H=zero(1,total);H[0][hindex]=I(1)
    return SRC.ball_sector(vmap,H,float(bound))
def beta_sector(total,beta_norm_upper):
    return ('physical_true_bias_ball',_ball_sector(_selector(total,range(OFF_BETA,OFF_BETA+3)),IDX_HB,beta_norm_upper))
def bias_supply_sector(total,bias_norm_upper):
    start=NP+NLOCAL_PHYS
    return ('joint_bias_w_mtau_ball',_ball_sector(_selector(total,range(start,start+NLOCAL_BIAS)),IDX_HB,bias_norm_upper))
def prediction_sectors(total,h,bias_norm_upper,beta_norm_upper):
    if total!=NP+NLOCAL_PHYS+NLOCAL_BIAS+NFP:raise ValueError('prediction local dimension mismatch')
    A0=zero(3,total);A1=zero(3,total);M=zero(9,total);Hs=zero(1,total);Hs[0][IDX_HP]=I(1)
    for axis in range(3):
        A0[axis][OFF_A+axis]=I(1);A1[axis][NP+axis]=I(1)
        M[3*axis+0][NP+3+axis]=I(1/h);M[3*axis+1][NP+6+axis]=I(1/(h*h));M[3*axis+2][NP+9+axis]=I(1/(h*h*h))
    Amax=float(SRC.MOM.build()['A_max_mps2']);phys=list(SRC.physical_witness_sectors(A0,A1,M,Hs,Amax))
    return tuple(phys+[bias_supply_sector(total,bias_norm_upper),beta_sector(total,beta_norm_upper)])
def ordinary_sectors(total,beta_norm_upper):
    if total!=NP+NFP:raise ValueError('ordinary local dimension mismatch')
    return (beta_sector(total,beta_norm_upper),)
def _pull(Q,T):return matrix_mul(matrix_mul(matrix_transpose(T),Q),T)
def _supply(total,gamma_phys,gamma_bias,gamma_n,prediction,bias_supply_on_this_event):
    S=zero(total,total)
    if prediction:S[IDX_HP][IDX_HP]=I(gamma_phys)
    if prediction or bias_supply_on_this_event:S[IDX_HB][IDX_HB]=I(gamma_bias)
    start=NP+(NLOCAL_PHYS+NLOCAL_BIAS if prediction else 0)
    for i in range(NFP):S[start+i][start+i]=I(gamma_n)
    return S
def local_master(Qnext,T,*,gamma_phys,gamma_bias,gamma_n,sectors=(),multipliers=(),prediction=False,bias_supply_on_this_event=False):
    total=shape(T)[1]
    K=matrix_sub(_supply(total,gamma_phys,gamma_bias,gamma_n,prediction,bias_supply_on_this_event),_pull(Qnext,T))
    if len(sectors)!=len(multipliers):raise ValueError('sector/multiplier mismatch')
    for (_,P),lam in zip(sectors,multipliers):
        if float(lam)<0:raise ValueError('S-procedure multiplier must be nonnegative')
        K=matrix_sub(K,scale(P,lam))
    return matrix_symmetric_hull(K)
def schur_required_state_form(K):
    n=shape(K)[0]
    if shape(K)!=(n,n) or n<=NP:raise ValueError('local master must include eliminable inputs')
    Kxx=[row[:NP] for row in K[:NP]];Kxu=[row[NP:] for row in K[:NP]];Kux=[row[:NP] for row in K[NP:]];Kuu=[row[NP:] for row in K[NP:]]
    ok,piv=symmetric_positive_definite_ldlt(matrix_symmetric_hull(Kuu))
    if not ok:raise ValueError('local eliminable-input block is not rigorously positive definite')
    inv=matrix_inverse_gauss_jordan(Kuu);cross=matrix_mul(matrix_mul(Kxu,inv),Kux)
    return matrix_symmetric_hull(matrix_sub(cross,Kxx)),Kuu,[float(x.lo) for x in piv]
def _smoke():
    supply=BIAS.build()['family_supply']['BIAS1'];bnd=float(supply['joint_supply_norm_upper_per_prediction_mps2']);beta=float(supply['true_bias_norm_upper_mps2'])
    A=eye(NJOINT);B=zero(NJOINT,NLOCAL_BIAS)
    for i in range(3):B[18+i][i]=I(1);B[18+i][3+i]=I(1);B[21+i][i]=I(1)
    T=prediction_transition(A,B,I(2),I(.005),1e-6);sec=prediction_sectors(len(T[0]),.005,bnd,beta)
    To=ordinary_transition(A,1e-6);So=_supply(len(To[0]),0,7,1,False,True)
    return {'prediction_transition_shape':shape(T),'sector_names':[x[0] for x in sec],'same_w_column_present':all(B[18+i][i].contains(1) and B[21+i][i].contains(1) for i in range(3)),'physical_and_bias_scales_distinct':IDX_HP!=IDX_HB,'ordinary_prefix_can_receive_single_bias_supply':So[IDX_HB][IDX_HB].contains(7),'ordinary_prefix_does_not_receive_physical_supply':So[IDX_HP][IDX_HP].contains(0)}
def build():
    p=PHYS.build();pf=PHYS.validate(p);s=SRC.build();sf=SRC.validate(s);b=BIAS.build();bf=BIAS.validate(b)
    if pf or sf or bf:raise RuntimeError(f'A21 Schur prerequisites failed physical={pf} source={sf} bias={bf}')
    sm=_smoke()
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,
      'persistent_coordinate':['e_A21_21','b_true_3','a_phys_3','h_physical_1','h_bias_1'],'persistent_dimension':NP,'prediction_local_physical_dimension':NLOCAL_PHYS,'prediction_local_bias_dimension':NLOCAL_BIAS,'event_local_binary32_dimension':NFP,
      'same_physical_acceleration_endpoint_carried':True,'same_true_bias_state_carried_through_projection':True,'same_w_mtau_6D_supply_map_consumed':True,'physical_and_bias_source_scales_distinct':True,'source_scales_independent_of_hard_entry_radial':True,'physical_a0_a1_J012_sectors_available':True,'joint_bias_supply_ball_sector_available':True,'absolute_true_bias_ball_sector_available_at_every_event':True,'prefix_level_bias_supply_can_be_placed_once_without_new_driver':True,'prefix_level_bias_supply_not_charged_per_event':True,'positive_form_S_procedure_uses_minus_lambda_Pi':True,'outward_interval_Schur_complement_available':True,
      'smoke':sm,'complete_word_backward_recursion_closed_here':False,'endpoint_contraction_closed_here':False,'every_prefix_bound_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'run this operator backward over each BIAS0/1/2 A21 post-prediction word and every literal prefix, placing one kappa_b D_b allowance at each prefix endpoint when no new prediction supply exists'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('same_physical_acceleration_endpoint_carried','same_true_bias_state_carried_through_projection','same_w_mtau_6D_supply_map_consumed','physical_and_bias_source_scales_distinct','source_scales_independent_of_hard_entry_radial','physical_a0_a1_J012_sectors_available','joint_bias_supply_ball_sector_available','absolute_true_bias_ball_sector_available_at_every_event','prefix_level_bias_supply_can_be_placed_once_without_new_driver','prefix_level_bias_supply_not_charged_per_event','positive_form_S_procedure_uses_minus_lambda_Pi','outward_interval_Schur_complement_available'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('complete_word_backward_recursion_closed_here','endpoint_contraction_closed_here','every_prefix_bound_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('persistent_dimension')!=29:f.append('persistent dimension changed')
    sm=d.get('smoke',{})
    if tuple(sm.get('prediction_transition_shape',()))!=(29,68):f.append('prediction transition shape changed')
    if sm.get('sector_names')!=['physical_a0_ball','physical_a1_ball','physical_J012_moment_iqc','joint_bias_w_mtau_ball','physical_true_bias_ball']:f.append('sector set changed')
    if sm.get('same_w_column_present') is not True:f.append('shared w lost')
    if sm.get('ordinary_prefix_can_receive_single_bias_supply') is not True or sm.get('ordinary_prefix_does_not_receive_physical_supply') is not True:f.append('prefix bias-supply placement smoke failed')
    return f
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'persistent_dim':d['persistent_dimension'],'sectors':d['smoke']['sector_names'],'prefix_bias_supply':d['smoke']['ordinary_prefix_can_receive_single_bias_supply'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
