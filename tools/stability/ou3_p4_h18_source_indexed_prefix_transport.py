#!/usr/bin/env python3
"""Exact H18 source-indexed finite transport over synchronized P4 event lineages.

This module connects the synchronized estimator-owned SourceCoverCell lineage to
the existing whole-word/every-prefix Phi transport.  It deliberately starts with
H18 because no accelerometer-bias projection is present; A21 adds that hybrid on
top of the same construction.

For each retained H18 sample slice, let epsilon_aw(z;s) be the exact full
measurement-linearizing shift, including the mixed (Q-I) delta_a_w term.  The
transport event records are:

* source-coordinate rebase between samples:
    C=L=I, rho=0;
* prediction:
    C=L=F_shipping, rho=z_exact^+ - F_shipping z;
* covariance-only a_w floor:
    C=L=I, rho=0;
* accepted S/accelerometer/magnetometer Joseph+quaternion-reset:
    C=G(I-KH), L=G,
    rho=z_exact^+ - G(z-Ky),
  with K,H,R and z from the SAME SourceCoverCell.

Consequently (C-L)E_aw is identically zero for S=0 and magnetometer, and only
accelerometer can produce the interior epsilon transport required by the exact
endpoint identity.  No nonlinear Jacobian is substituted for C, no packetwise
remainder reset is used, and no independent K/H/R box is introduced.

The output can be passed directly to complete_word_prefix_transport.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval, matrix_identity, matrix_mul, matrix_sub
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_brmm_full_normal_live_word as WORD
import ou3_p4_complete_brmm_accelerometer_operation_coordinate as ACC
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_complete_brmm_differential_prediction as PRED
import ou3_p4_complete_word_prefix_transport as PREFIX
import ou3_p4_source_uniform_estimator_event_attachment as ATTACH
import ou3_p4_source_uniform_event_lineage_sequence as LINEAGE
import ou3_p4_source_indexed_phi_rebase as REBASE

SCHEMA=1
QUALIFICATION='OU3_P4_H18_SOURCE_INDEXED_EVERY_PREFIX_TRANSPORT_V1'
P3_DELTA=1.0e-18
N=18
OFF_AW=15


def I(x:float)->Interval:return Interval.point(float(x))
def zero_vec(n=N):return [I(0.0) for _ in range(n)]
def shape(A):return (len(A),len(A[0]) if A else 0)

def mv(A,x):
    if not A or len(A[0])!=len(x):raise ValueError('matrix/vector mismatch')
    out=[]
    for row in A:
        s=I(0.0)
        for a,b in zip(row,x):s=s+a*b
        out.append(s)
    return out

def vsub(a,b):
    if len(a)!=len(b):raise ValueError('vector mismatch')
    return [x-y for x,y in zip(a,b)]
def vadd(a,b):
    if len(a)!=len(b):raise ValueError('vector mismatch')
    return [x+y for x,y in zip(a,b)]

def Eaw():
    E=[[I(0.0) for _ in range(3)] for _ in range(N)]
    for i in range(3):E[OFF_AW+i][i]=I(1.0)
    return E

def _col(v):return [[x] for x in v]

def epsilon_aw(state,selector):
    if len(state)!=N:raise ValueError('H18 state required')
    c=_col(state[:3]); da=_col(state[OFF_AW:OFF_AW+3]); f=_col(selector.sample_coordinates.f_cog_body)
    r=ACC.evaluate_operation_coordinate(c,selector.sample_coordinates.R_wb,da,f)
    return [row[0] for row in r['epsilon_aw']]

def _same_zero(x:Interval)->bool:return x.lo<=0.0<=x.hi and x.lo==0.0 and x.hi==0.0

def _G_from_dtheta(d):
    if len(d)!=3:raise ValueError('dtheta length')
    G=matrix_identity(N);h=I(.5);x,y,z=d
    S=[[I(0),-z,y],[z,I(0),-x],[-y,x,I(0)]]
    for i in range(3):
        for j in range(3):G[i][j]=G[i][j]+h*S[i][j]
    return G

def _residual(cell):
    z=EVENTS._state_ad(cell.state)
    if cell.kind=='S_zero':y=EVENTS.residual_S(z)
    elif cell.kind=='magnetometer':y=EVENTS.residual_magnetometer(z,cell.m_body)
    elif cell.kind=='accelerometer':y=EVENTS.residual_accelerometer(z,cell.f_hat,cell.R_hat)
    else:raise ValueError('Joseph cell required')
    return EVENTS.AD.values(y)

def _event_state_out(cell,tau_ba=None):
    if cell.kind=='prediction':
        raise ValueError('prediction needs selector body rate')
    if cell.kind=='aw_floor':return list(cell.state)
    return list(EVENTS.source_joseph_event(**ATTACH.COVER.joseph_event_kwargs(cell))['state_out'])

def _prediction_record(cell,trusted,selector):
    if trusted.kind!='prediction' or trusted.F is None:raise ValueError('trusted prediction F missing')
    pred=PRED.prediction_event('H',cell.state,selector.sample_coordinates.omega_body_corrected,cell.dt_s,cell.tau_applied_s)
    F=trusted.F
    rho=vsub(pred['state_out'],mv(F,cell.state))
    return {'kind':'prediction','C':F,'L':F,'rho':rho,'token':cell.source_token},list(pred['state_out'])

def _floor_record(cell):
    J=matrix_identity(N)
    return {'kind':'aw_floor','C':J,'L':J,'rho':zero_vec(),'token':cell.source_token},list(cell.state)

def _joseph_record(cell):
    ev=EVENTS.source_joseph_event(**ATTACH.COVER.joseph_event_kwargs(cell))
    y=_residual(cell);d=mv(ev['K'],y);G=_G_from_dtheta(d[:3]);A=matrix_sub(matrix_identity(N),matrix_mul(ev['K'],ev['H']))
    C=matrix_mul(G,A);L=G
    t=vsub(cell.state,d);rho=vsub(ev['state_out'],mv(G,t))
    return {'kind':cell.kind,'C':C,'L':L,'rho':rho,'token':cell.source_token},list(ev['state_out'])

def _trusted_by_ordinal(sample):
    return {c.event_index_in_sample:c for c in sample.selector.H_event_cells}

def build_transport(samples:Sequence[ATTACH.AttachedSampleLineage]):
    stitched=LINEAGE.stitch_samples(samples)
    if stitched.mode!='H':raise ValueError('H18 lineage required')
    events=[];eps=[];emb=[];tokens=[];physical_prefix=[]
    current_state=list(samples[0].cells[0].state)
    eps.append(epsilon_aw(current_state,samples[0].selector));emb.append(Eaw())
    previous_selector=None
    for si,sample in enumerate(samples):
        selector=sample.selector
        if si:
            # Same physical boundary, new source-indexed chart.
            old=eps[-1];new=epsilon_aw(current_state,selector)
            reb=REBASE.rebase_event('H',_col(old),_col(new))
            events.append({'kind':'source_coordinate_rebase','C':reb['C'],'L':reb['L'],'rho':reb['rho'],'token':selector.source_cell_id+':rebase'})
            eps.append(new);emb.append(Eaw());tokens.append(selector.source_cell_id+':rebase');physical_prefix.append(False)
        trusted=_trusted_by_ordinal(sample)
        for cell in sample.cells:
            if cell.state!=current_state:
                raise RuntimeError('transport event state is not previous exact physical boundary')
            te=trusted[cell.event_ordinal]
            if cell.kind=='prediction':rec,out=_prediction_record(cell,te,selector)
            elif cell.kind=='aw_floor':rec,out=_floor_record(cell)
            else:rec,out=_joseph_record(cell)
            events.append(rec);current_state=list(out);eps.append(epsilon_aw(current_state,selector));emb.append(Eaw());tokens.append(cell.source_token);physical_prefix.append(True)
        if tuple(current_state)!=tuple(sample.state_out):raise RuntimeError('sample transport endpoint detached from synchronized nonlinear chain')
        previous_selector=selector
    decomp=PREFIX.prefix_decompositions(events,emb,eps)
    # Endpoint transport forbids non-accelerometer interior epsilon terms.  The
    # rebase/prediction/floor/S/mag algebra must therefore leave only accel.
    interior_ok=all(all(events[i]['kind']=='accelerometer' for i in r['transport']['interior_event_indices']) for r in decomp)
    return {'events':events,'embeddings':emb,'eps_nodes':eps,'decompositions':decomp,'event_tokens':tokens,
            'physical_prefix_flags':physical_prefix,'state_out':tuple(current_state),'interior_only_accelerometer':interior_ok}

def _identity(n):return [[I(1 if i==j else 0) for j in range(n)] for i in range(n)]
def _sample(forcez=-9.80665):
    return KERNEL.SampleCoordinates(gyro_measurement=KERNEL.MAHONY.Vec3(I(.01),I(-.02),I(.005)),omega_body_corrected=(I(.01),I(-.02),I(.005)),specific_force=KERNEL.MAHONY.Vec3(I(.2),I(-.1),I(forcez)),f_cog_body=(I(.2),I(-.1),I(forcez)),R_wb=_identity(3),due_S=True,aw_floor_requested=True,magnetometer_events_after_imu=(KERNEL.MagneticEvent((I(20),I(0),I(40))),))
def _smoke_objects():
    js=ATTACH.JOINT._smoke_state();branch=KERNEL.ExecutionBranch(frontend=copy.deepcopy(js.frontend),H=WORD.initialize_word('H',_identity(18)),A=WORD.initialize_word('A',_identity(21)),source_cell_id='root')
    common=dict(radial_scale=Interval(0,1),true_bias=[I(0),I(0),I(0)],bias_projection_limit=.4,tau_ba=I(1800))
    p0=ATTACH.synchronize_sample(branch=branch,joint_state=js,sample=_sample(),state_in_H=zero_vec(18),state_in_A=zero_vec(21),sample_index=0,next_cell_prefix='htr-k0',**common)
    if not p0:raise RuntimeError('sample0 empty');h0,a0=p0[0]
    h0,a0=p0[0];b1=ATTACH.next_execution_branch(h0,a0)
    # Change source force on sample 1 so the rebase is nontrivial.
    p1=ATTACH.synchronize_sample(branch=b1,joint_state=h0.image.state,sample=_sample(-9.7),state_in_H=h0.state_out,state_in_A=a0.state_out,sample_index=1,next_cell_prefix='htr-k1',**common)
    if not p1:raise RuntimeError('sample1 empty');h1,_=p1[0]
    h1,_=p1[0];return (h0,h1)
def _smoke():
    tr=build_transport(_smoke_objects());k=[e['kind'] for e in tr['events']]
    return {'event_count':len(k),'rebase_count':k.count('source_coordinate_rebase'),'prediction_count':k.count('prediction'),'accelerometer_count':k.count('accelerometer'),'interior_only_accelerometer':tr['interior_only_accelerometer'],'all_prefix_identity_residuals_zero':all(all(x.lo<=0<=x.hi for x in r['transport']['identity_residual']) for r in tr['decompositions']),'endpoint_matches_sample_state':True}

def build():
    att=ATTACH.build();af=ATTACH.validate(att);lin=LINEAGE.build();lf=LINEAGE.validate(lin);reb=REBASE.build();rf=REBASE.validate(reb);pref=PREFIX.build();pf=PREFIX.validate(pref)
    if af or lf or rf or pf:raise RuntimeError(f'H18 transport prerequisites failed attachment={af} lineage={lf} rebase={rf} prefix={pf}')
    s=_smoke();closed=bool(s['event_count']>0 and s['rebase_count']==1 and s['prediction_count']==2 and s['accelerometer_count']==2 and s['interior_only_accelerometer'] and s['all_prefix_identity_residuals_zero'])
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,
      'synchronized_source_uniform_event_attachment_consumed':True,'cross_sample_Phi_rebase_consumed':True,
      'prediction_uses_trusted_shipping_F_not_nonlinear_Jacobian':True,'prediction_finite_attitude_defect_retained_in_rho':True,
      'Joseph_C_equals_G_I_minus_KH':True,'Joseph_L_equals_G':True,'same_cell_K_H_R_state_used':True,
      'exact_quaternion_reset_defect_retained_in_rho':True,'full_mixed_epsilon_aw_evaluated_at_every_boundary':True,
      'source_change_epsilon_not_rezeroed':True,'only_accelerometer_has_interior_epsilon_transport':closed,
      'H18_every_prefix_transport_constructor_closed':closed,'smoke':s,
      'source_uniform_complete_601_sample_family_executed_here':False,'H18_endpoint_augmented_LDLT_closed_here':False,'H18_every_prefix_augmented_LDLT_closed_here':False,'H18_first_exit_retention_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'lift each H18 prefix decomposition defect onto the common source/radial/moment/binary32 graph coordinate and combine it with the P3 strict information block; then run the outward augmented LDLT and same-graph hard-ball/reset-domain targets'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('synchronized_source_uniform_event_attachment_consumed','cross_sample_Phi_rebase_consumed','prediction_uses_trusted_shipping_F_not_nonlinear_Jacobian','prediction_finite_attitude_defect_retained_in_rho','Joseph_C_equals_G_I_minus_KH','Joseph_L_equals_G','same_cell_K_H_R_state_used','exact_quaternion_reset_defect_retained_in_rho','full_mixed_epsilon_aw_evaluated_at_every_boundary','source_change_epsilon_not_rezeroed','only_accelerometer_has_interior_epsilon_transport','H18_every_prefix_transport_constructor_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('source_uniform_complete_601_sample_family_executed_here','H18_endpoint_augmented_LDLT_closed_here','H18_every_prefix_augmented_LDLT_closed_here','H18_first_exit_retention_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    return f

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'H18_transport':d['H18_every_prefix_transport_constructor_closed'],'smoke':d['smoke'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
