#!/usr/bin/env python3
"""Full joint24 physical event/prediction lifts for ALT.

H18 proof state is [e_H18,e_ba,beta_true]. A21 proof state is
[e_A21,beta_true]. These lifts retain physical bias in H18 accelerometer
residuals and physical S in S=0 events. Projection uses the same generalized
Jacobian as shipping, including beta_true derivatives.
"""
from __future__ import annotations
from ou3_interval import Interval,matrix_mul
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_complete_brmm_source_cover_contract as COVER
from tools.stability.ou3_alt_contraction import physical_reference as REF

QUALIFICATION='OU3_ALT_PHYSICAL_JOINT24_EVENTS_V1'
H_EBA=18;H_BETA=21;A_BETA=21

def I(x):return Interval.point(float(x))
def zero(r,c):return [[I(0) for _ in range(c)] for _ in range(r)]
def _selector(rows,cols,offset):
    A=zero(rows,cols)
    for i in range(rows):A[i][offset+i]=I(1)
    return A

def _h18_accel_residual(z,eba,f_hat,R_hat):
    nd=z[0].n;E=EVENTS.AD.rotation_from_cayley(z[:3]);f=[EVENTS._const(x,nd) for x in f_hat]
    R=[[EVENTS._const(x,nd) for x in row] for row in R_hat];da=z[EVENTS.OFF_AW:EVENTS.OFF_AW+3]
    Rda=EVENTS.AD.matvec(R,da);Ef=EVENTS.AD.matvec(E,f);ERda=EVENTS.AD.matvec(E,Rda)
    return [Ef[i]-f[i]+ERda[i]+eba[i] for i in range(3)]

def h18_accelerometer_event(cell:COVER.SourceCoverCell):
    f=COVER.validate_cell(cell,require_estimator_provenance=True)
    if f:raise ValueError('invalid source cell: '+repr(f))
    if cell.mode!='H' or cell.kind!='accelerometer':raise ValueError('H18 accelerometer cell required')
    u=EVENTS.AD.independent_vector(list(cell.state)+[I(0)]*6,n=24,offset=0);z=u[:18];eba=u[18:21];beta=u[21:24]
    H=EVENTS._event_H('H','accelerometer',f_hat=cell.f_hat,R_hat=cell.R_hat);K,S=EVENTS.source_joseph_gain(cell.P,H,cell.R)
    residual=_h18_accel_residual(z,eba,cell.f_hat,cell.R_hat);out18=EVENTS._apply_physical_correction(z,K,residual);out=list(out18)+list(eba)+list(beta)
    return {'J_joint24':EVENTS.AD.jacobian(out),'state_out':EVENTS.AD.values(out),'residual':EVENTS.AD.values(residual),'same_history_held_bias_retained':True,'beta_identity':True,'independent_bias_supply_used':False}

def h18_prediction_lift(J18,phi_true:Interval):
    if len(J18)!=18 or any(len(r)!=18 for r in J18):raise ValueError('18x18 prediction Jacobian required')
    if not isinstance(phi_true,Interval):raise TypeError('phi_true interval required')
    A=zero(24,24);B=zero(24,3)
    for i in range(18):
        for j in range(18):A[i][j]=J18[i][j]
    for i in range(3):
        A[18+i][18+i]=I(1);A[18+i][21+i]=phi_true-I(1);A[21+i][21+i]=phi_true;B[18+i][i]=I(1);B[21+i][i]=I(1)
    return A,B

def _project_joint(out21,beta,total_dim,radius):
    pre=out21[18:21];xhat=[beta[i]-pre[i] for i in range(3)];proj=EVENTS.ball_projection_enclosure(EVENTS.AD.values(xhat),radius)
    smooth=EVENTS.AD.jacobian(out21);dx=EVENTS.AD.jacobian(xhat);pder=matrix_mul(proj['J'],dx);bsel=_selector(3,total_dim,A_BETA);J=[list(r) for r in smooth]
    for i in range(3):J[18+i]=[bsel[i][j]-pder[i][j] for j in range(total_dim)]
    vals=EVENTS.AD.values(out21);vals[18:21]=[EVENTS.AD.values(beta)[i]-proj['value'][i] for i in range(3)]
    return J,vals,proj

def s_zero_event_joint24(cell:COVER.SourceCoverCell):
    f=COVER.validate_cell(cell,require_estimator_provenance=True)
    if f:raise ValueError('invalid source cell: '+repr(f))
    if cell.kind!='S_zero' or cell.wave_primitive is None:raise ValueError('physical S_zero cell with primitive required')
    REF.validate_primitive_reference(cell.wave_primitive)
    if cell.R_provenance!=EVENTS.ACTUAL_RS_PROVENANCE:raise ValueError('actual R_S provenance required')
    if cell.mode=='H':state24=list(cell.state)+[I(0)]*6
    else:
        if cell.true_bias is None or cell.bias_projection_limit is None:raise ValueError('A21 true bias/projection required')
        state24=list(cell.state)+list(cell.true_bias)
    u=EVENTS.AD.independent_vector(state24+list(cell.wave_primitive.centered_S),n=27,offset=0);nerr=18 if cell.mode=='H' else 21;z=u[:nerr]
    beta=u[H_BETA:H_BETA+3] if cell.mode=='H' else u[A_BETA:A_BETA+3];sphys=u[24:27]
    H=EVENTS._event_H(cell.mode,'S_zero');K,S=EVENTS.source_joseph_gain(cell.P,H,cell.R);residual=[z[EVENTS.OFF_S+i]-sphys[i] for i in range(3)];out_err=EVENTS._apply_physical_correction(z,K,residual)
    if cell.mode=='H':
        out=list(out_err)+list(u[18:21])+list(beta);J=EVENTS.AD.jacobian(out);vals_out=EVENTS.AD.values(out);branch='not_applicable'
    else:
        Jerr,errvals,proj=_project_joint(out_err,beta,27,float(cell.bias_projection_limit));branch=proj['branch'];J=Jerr+_selector(3,27,A_BETA);vals_out=errvals+EVENTS.AD.values(beta)
    return {'J_joint24':[r[:24] for r in J],'J_S_phys':[r[24:27] for r in J],'state_out':vals_out,'residual':EVENTS.AD.values(residual),'physical_generator_id':cell.wave_primitive.generator_id,'live_origin_id':cell.wave_primitive.live_origin_id,'projection_branch':branch,'same_history_S_source_columns_retained':True,'independent_S_box_used':False}

def build():
    return {'qualification':QUALIFICATION,'joint_dimension':24,'H18_accelerometer_held_bias_columns_available':True,'H18_prediction_true_bias_shared_driver_lift_available':True,'H18_and_A21_physical_S_source_columns_available':True,'A21_projection_beta_columns_retained':True,'independent_bias_supply_used':False,'independent_S_box_used':False,'complete_literal_word_closed':False,'storage_search_allowed':False,'ALT_LIVE_PASS':False}
def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION or d.get('joint_dimension')!=24:f.append('qualification/dimension mismatch')
    for k in ('H18_accelerometer_held_bias_columns_available','H18_prediction_true_bias_shared_driver_lift_available','H18_and_A21_physical_S_source_columns_available','A21_projection_beta_columns_retained'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('independent_bias_supply_used','independent_S_box_used','complete_literal_word_closed','storage_search_allowed','ALT_LIVE_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
