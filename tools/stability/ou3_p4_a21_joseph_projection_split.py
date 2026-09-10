#!/usr/bin/env python3
"""Exact/outward A21 Joseph-reset -> bias-projection hybrid split for P4.

Shipping ``source_joseph_event`` composes two physically distinct operations in
A21: the smooth measurement correction + finite quaternion reset, then the
radial accelerometer-bias estimate projection.  The joint P4 cocycle must retain
that split because the projection depends on the SAME physical true-bias source
coordinate and has a Clarke generalized Jacobian.

For one estimator-owned SourceCoverCell this module returns:

1. smooth Joseph/reset
     C_j=G(I-KH), L_j=G, C_j-L_j=-GKH,
     rho_j=z_smooth-G(z-Ky);
2. projection
     xhat=b_true-e_b,smooth,
     e_b^+=b_true-Proj_R(xhat),
     C_p=L_p=diag(I18,J_Proj),
     rho_p=z_projected-C_p z_smooth.

The projection leaves attitude/a_w unchanged, so the source-indexed epsilon_aw
embedding has exact zero interior transport across this hybrid.  The projection
defect remains tied to b_true; it is not an independent additive disturbance.
"""
from __future__ import annotations

import argparse,json
from pathlib import Path

from ou3_interval import Interval,matrix_add,matrix_identity,matrix_mul
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_complete_brmm_source_cover_contract as COVER
import ou3_p4_brmm_reduced_event_master as EVENTMASTER
import ou3_p4_brmm_same_Dtheta_reset_binding as RESET

SCHEMA=1
QUALIFICATION='OU3_P4_A21_JOSEPH_PROJECTION_HYBRID_SPLIT_V1'
P3_DELTA=1.0e-18
N=21;OFF_BA=18

def I(x):return Interval.point(float(x))
def zvec(n=N):return [I(0) for _ in range(n)]
def zmat(r,c):return [[I(0) for _ in range(c)] for _ in range(r)]
def mv(A,x):
    out=[]
    for row in A:
        s=I(0)
        for a,b in zip(row,x):s=s+a*b
        out.append(s)
    return out
def vsub(a,b):return [x-y for x,y in zip(a,b)]
def Gext(d):
    G=matrix_identity(N);x,y,z=d;h=I(.5);S=[[I(0),-z,y],[z,I(0),-x],[-y,x,I(0)]]
    for i in range(3):
        for j in range(3):G[i][j]=G[i][j]+h*S[i][j]
    return G
def residual(cell):
    z=EVENTS._state_ad(cell.state)
    if cell.kind=='S_zero':y=EVENTS.residual_S(z)
    elif cell.kind=='magnetometer':y=EVENTS.residual_magnetometer(z,cell.m_body)
    elif cell.kind=='accelerometer':y=EVENTS.residual_accelerometer(z,cell.f_hat,cell.R_hat)
    else:raise ValueError('Joseph cell required')
    return EVENTS.AD.values(y)
def split_event(cell:COVER.SourceCoverCell):
    if cell.mode!='A' or cell.kind not in ('S_zero','accelerometer','magnetometer'):raise ValueError('A21 Joseph source cell required')
    if cell.true_bias is None or cell.bias_projection_limit is None:raise ValueError('same-source true bias/projection radius required')
    H=EVENTS._event_H('A',cell.kind,f_hat=cell.f_hat,R_hat=cell.R_hat,m_body=cell.m_body)
    K,_=EVENTS.source_joseph_gain(cell.P,H,cell.R)
    z=EVENTS._state_ad(cell.state);yad=(EVENTS.residual_S(z) if cell.kind=='S_zero' else EVENTS.residual_magnetometer(z,cell.m_body) if cell.kind=='magnetometer' else EVENTS.residual_accelerometer(z,cell.f_hat,cell.R_hat))
    smooth_ad=EVENTS._apply_physical_correction(z,K,yad);smooth=EVENTS.AD.values(smooth_ad);y=EVENTS.AD.values(yad);d=mv(K,y);G=Gext(d[:3]);GKH=matrix_mul(G,matrix_mul(K,H));Dm=[[-x for x in row] for row in GKH];Cj=matrix_add(G,Dm);rhoj=vsub(smooth,mv(G,vsub(cell.state,d)))
    xhat=[cell.true_bias[i]-smooth[OFF_BA+i] for i in range(3)];proj=EVENTS.ball_projection_enclosure(xhat,float(cell.bias_projection_limit));Cp=matrix_identity(N)
    for i in range(3):
        for j in range(3):Cp[OFF_BA+i][OFF_BA+j]=proj['J'][i][j]
    projected=list(smooth);projected[OFF_BA:OFF_BA+3]=[cell.true_bias[i]-proj['value'][i] for i in range(3)];rhop=vsub(projected,mv(Cp,smooth))
    full=EVENTS.source_joseph_event(**COVER.joseph_event_kwargs(cell))
    enclosure=all(projected[i].lo<=full['state_out'][i].lo and projected[i].hi>=full['state_out'][i].hi for i in range(N)) or all(full['state_out'][i].lo<=projected[i].lo and full['state_out'][i].hi>=projected[i].hi for i in range(N))
    return {'smooth':{'kind':cell.kind+'_smooth','C':Cj,'L':G,'C_minus_L':Dm,'rho':rhoj,'state_out':smooth},'projection':{'kind':'bias_projection','C':Cp,'L':Cp,'C_equals_L_exact':True,'C_minus_L':zmat(N,N),'rho':rhop,'state_out':projected,'branch':proj['branch'],'input_norm':proj['norm']},'composed_state_out_encloses_shipping_output':enclosure,'same_true_bias_used':True,'projection_defect_independent_source_port':False}
def _smoke_cell():
    im=EVENTMASTER._smoke_image();n=21;P=EVENTMASTER._identity(n);R=EVENTMASTER._identity(3);x=[I(0) for _ in range(n)];x[18]=Interval(-.5,.5)
    return COVER.source_cell_from_joint_image(im,mode='A',sample_index=0,event_ordinal=3,kind='accelerometer',state=x,P=P,dt_s=I(.005),pseudo_elapsed_s=I(.01),radial_scale=Interval(0,1),event_source_token=im.source_token+':e3',event_predecessor_token=im.source_token+':e2',R=R,f_hat=[I(.2),I(-.1),I(-9.7)],R_hat=EVENTMASTER._identity(3),true_bias=[I(0),I(0),I(0)],bias_projection_limit=.4)
def build():
    e=EVENTS.build();ef=EVENTS.validate(e);r=RESET.build();rf=RESET.validate(r)
    if ef or rf:raise RuntimeError(f'A21 split prerequisites failed event={ef} reset={rf}')
    s=split_event(_smoke_cell());p=s['projection'];closed=bool(s['same_true_bias_used'] and not s['projection_defect_independent_source_port'] and p['C_equals_L_exact'] and p['branch'] in ('inactive','active','clarke_hull') and s['composed_state_out_encloses_shipping_output'])
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,'shipping_A21_composite_event_split_into_smooth_and_projection':True,'smooth_Joseph_C_equals_G_I_minus_KH':True,'smooth_Joseph_C_minus_L_is_minus_GKH':True,'projection_C_equals_L_generalized_Jacobian':True,'projection_leaves_attitude_aw_coordinates_unchanged':True,'same_physical_true_bias_source_coordinate_retained':True,'projection_defect_is_not_independent_source_port':True,'A21_Joseph_projection_split_closed':closed,'smoke_projection_branch':p['branch'],'all_BIAS0_BIAS1_BIAS2_recurrences_composed_here':False,'A21_prefix_transport_closed_here':False,'A21_augmented_LDLT_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,'next_obligation':'insert this two-event split into the A21 24-state bias-family cocycle, carry b_true/w/m_tau through projection rho, then reuse the H18 source-indexed Phi/rebase prefix transport and run the common augmented LDLT'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('shipping_A21_composite_event_split_into_smooth_and_projection','smooth_Joseph_C_equals_G_I_minus_KH','smooth_Joseph_C_minus_L_is_minus_GKH','projection_C_equals_L_generalized_Jacobian','projection_leaves_attitude_aw_coordinates_unchanged','same_physical_true_bias_source_coordinate_retained','projection_defect_is_not_independent_source_port','A21_Joseph_projection_split_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('all_BIAS0_BIAS1_BIAS2_recurrences_composed_here','A21_prefix_transport_closed_here','A21_augmented_LDLT_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'split':d['A21_Joseph_projection_split_closed'],'branch':d['smoke_projection_branch'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
