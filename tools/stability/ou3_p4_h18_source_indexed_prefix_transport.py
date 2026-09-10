#!/usr/bin/env python3
"""Exact H18 source-indexed finite transport over synchronized P4 event lineages.

For each retained H18 sample slice, the full measurement-linearizing shift
``epsilon_aw`` (including the mixed ``(Q-I) delta_a_w`` term) is evaluated at
every physical boundary.  Between samples an explicit source-coordinate rebase
changes only Phi, never the physical error.

Transport records preserve exact operation algebra:

* source rebase / prediction / covariance floor: ``C=L`` structurally;
* Joseph+reset: ``C=G(I-KH)``, ``L=G``, and the exact algebraic difference
  ``C-L=-GKH`` is recorded directly rather than recovered by interval
  subtraction;
* ``rho`` is the finite physical defect against the declared linear/chart map.

The interval-safe structural prefix transport therefore retains the canonical
cancellation: only accepted accelerometer events can carry an interior
``epsilon_aw`` term.  No nonlinear Jacobian substitutes for the shipping tangent
map and no packetwise remainder budget is introduced.
"""
from __future__ import annotations

import argparse,copy,json
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval,matrix_add,matrix_identity,matrix_mul,matrix_sub
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_brmm_full_normal_live_word as WORD
import ou3_p4_complete_brmm_accelerometer_operation_coordinate as ACC
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_complete_brmm_differential_prediction as PRED
import ou3_p4_source_uniform_estimator_event_attachment as ATTACH
import ou3_p4_source_uniform_event_lineage_sequence as LINEAGE
import ou3_p4_source_indexed_phi_rebase as REBASE
import ou3_p4_structural_prefix_transport as STRUCT

SCHEMA=2
QUALIFICATION='OU3_P4_H18_SOURCE_INDEXED_EVERY_PREFIX_TRANSPORT_V2'
P3_DELTA=1.0e-18
N=18;OFF_AW=15

def I(x):return Interval.point(float(x))
def zvec(n=N):return [I(0) for _ in range(n)]
def zmat(r,c):return [[I(0) for _ in range(c)] for _ in range(r)]
def mv(A,x):
    if not A or len(A[0])!=len(x):raise ValueError('matrix/vector mismatch')
    out=[]
    for row in A:
        s=I(0)
        for a,b in zip(row,x):s=s+a*b
        out.append(s)
    return out
def vsub(a,b):
    if len(a)!=len(b):raise ValueError('vector mismatch')
    return [x-y for x,y in zip(a,b)]
def col(v):return [[x] for x in v]
def Eaw():
    E=zmat(N,3)
    for i in range(3):E[OFF_AW+i][i]=I(1)
    return E

def epsilon_aw(state,selector):
    if len(state)!=N:raise ValueError('H18 state required')
    d=ACC.evaluate_operation_coordinate(col(state[:3]),selector.sample_coordinates.R_wb,col(state[OFF_AW:OFF_AW+3]),col(selector.sample_coordinates.f_cog_body))
    return [r[0] for r in d['epsilon_aw']]
def Gext(dtheta):
    if len(dtheta)!=3:raise ValueError('dtheta length')
    G=matrix_identity(N);x,y,z=dtheta;h=I(.5)
    S=[[I(0),-z,y],[z,I(0),-x],[-y,x,I(0)]]
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
def prediction_record(cell,trusted,selector):
    if trusted.kind!='prediction' or trusted.F is None:raise ValueError('trusted prediction F missing')
    p=PRED.prediction_event('H',cell.state,selector.sample_coordinates.omega_body_corrected,cell.dt_s,cell.tau_applied_s)
    F=trusted.F;rho=vsub(p['state_out'],mv(F,cell.state))
    return {'kind':'prediction','C':F,'L':F,'C_equals_L_exact':True,'C_minus_L':zmat(N,N),'rho':rho,'token':cell.source_token},list(p['state_out'])
def floor_record(cell):
    J=matrix_identity(N)
    return {'kind':'aw_floor','C':J,'L':J,'C_equals_L_exact':True,'C_minus_L':zmat(N,N),'rho':zvec(),'token':cell.source_token},list(cell.state)
def joseph_record(cell):
    ev=EVENTS.source_joseph_event(**ATTACH.COVER.joseph_event_kwargs(cell));y=residual(cell);d=mv(ev['K'],y);G=Gext(d[:3])
    KH=matrix_mul(ev['K'],ev['H']);minus_GKH=[[-x for x in row] for row in matrix_mul(G,KH)]
    C=matrix_add(G,minus_GKH);t=vsub(cell.state,d);rho=vsub(ev['state_out'],mv(G,t))
    return {'kind':cell.kind,'C':C,'L':G,'C_minus_L':minus_GKH,'rho':rho,'token':cell.source_token},list(ev['state_out'])
def trusted_map(sample):return {c.event_index_in_sample:c for c in sample.selector.H_event_cells}

def build_transport(samples:Sequence[ATTACH.AttachedSampleLineage]):
    stitched=LINEAGE.stitch_samples(samples)
    if stitched.mode!='H':raise ValueError('H18 lineage required')
    events=[];eps=[];emb=[];tokens=[];physical=[]
    state=list(samples[0].cells[0].state);eps.append(epsilon_aw(state,samples[0].selector));emb.append(Eaw())
    for si,sample in enumerate(samples):
        selector=sample.selector
        if si:
            old=eps[-1];new=epsilon_aw(state,selector);r=REBASE.rebase_event('H',col(old),col(new))
            events.append({'kind':'source_coordinate_rebase','C':r['C'],'L':r['L'],'C_equals_L_exact':True,'C_minus_L':zmat(N,N),'rho':r['rho'],'token':selector.source_cell_id+':rebase'})
            eps.append(new);emb.append(Eaw());tokens.append(selector.source_cell_id+':rebase');physical.append(False)
        tm=trusted_map(sample)
        for cell in sample.cells:
            if tuple(cell.state)!=tuple(state):raise RuntimeError('event state detached from previous physical boundary')
            if cell.kind=='prediction':rec,out=prediction_record(cell,tm[cell.event_ordinal],selector)
            elif cell.kind=='aw_floor':rec,out=floor_record(cell)
            else:rec,out=joseph_record(cell)
            events.append(rec);state=list(out);eps.append(epsilon_aw(state,selector));emb.append(Eaw());tokens.append(cell.source_token);physical.append(True)
        if tuple(state)!=tuple(sample.state_out):raise RuntimeError('sample endpoint detached from nonlinear chain')
    decomp=STRUCT.prefix_decompositions(events,emb,eps)
    interior_ok=all(all(events[i]['kind']=='accelerometer' for i in r['transport']['interior_event_indices']) for r in decomp)
    identity_ok=all(r['transport']['identity_residual_contains_zero'] for r in decomp)
    return {'events':events,'embeddings':emb,'eps_nodes':eps,'decompositions':decomp,'event_tokens':tokens,'physical_prefix_flags':physical,'state_out':tuple(state),'interior_only_accelerometer':interior_ok,'all_identity_residuals_contain_zero':identity_ok}

def ident(n):return [[I(1 if i==j else 0) for j in range(n)] for i in range(n)]
def sample(forcez=-9.80665):
    return KERNEL.SampleCoordinates(gyro_measurement=KERNEL.MAHONY.Vec3(I(.01),I(-.02),I(.005)),omega_body_corrected=(I(.01),I(-.02),I(.005)),specific_force=KERNEL.MAHONY.Vec3(I(.2),I(-.1),I(forcez)),f_cog_body=(I(.2),I(-.1),I(forcez)),R_wb=ident(3),due_S=True,aw_floor_requested=True,magnetometer_events_after_imu=(KERNEL.MagneticEvent((I(20),I(0),I(40))),))
def smoke_objects():
    js=ATTACH.JOINT._smoke_state();branch=KERNEL.ExecutionBranch(frontend=copy.deepcopy(js.frontend),H=WORD.initialize_word('H',ident(18)),A=WORD.initialize_word('A',ident(21)),source_cell_id='root');common=dict(radial_scale=Interval(0,1),true_bias=[I(0),I(0),I(0)],bias_projection_limit=.4,tau_ba=I(1800))
    p0=ATTACH.synchronize_sample(branch=branch,joint_state=js,sample=sample(),state_in_H=zvec(18),state_in_A=zvec(21),sample_index=0,next_cell_prefix='htr-k0',**common)
    if not p0:raise RuntimeError('sample0 empty')
    h0,a0=p0[0];b1=ATTACH.next_execution_branch(h0,a0)
    p1=ATTACH.synchronize_sample(branch=b1,joint_state=h0.image.state,sample=sample(-9.7),state_in_H=h0.state_out,state_in_A=a0.state_out,sample_index=1,next_cell_prefix='htr-k1',**common)
    if not p1:raise RuntimeError('sample1 empty')
    h1,_=p1[0];return h0,h1
def smoke():
    t=build_transport(smoke_objects());k=[e['kind'] for e in t['events']]
    return {'event_count':len(k),'rebase_count':k.count('source_coordinate_rebase'),'prediction_count':k.count('prediction'),'accelerometer_count':k.count('accelerometer'),'interior_only_accelerometer':t['interior_only_accelerometer'],'all_prefix_identity_residuals_contain_zero':t['all_identity_residuals_contain_zero']}
def build():
    att=ATTACH.build();af=ATTACH.validate(att);lin=LINEAGE.build();lf=LINEAGE.validate(lin);reb=REBASE.build();rf=REBASE.validate(reb);st=STRUCT.build();sf=STRUCT.validate(st)
    if af or lf or rf or sf:raise RuntimeError(f'H18 transport prerequisites failed attachment={af} lineage={lf} rebase={rf} structural={sf}')
    s=smoke();closed=bool(s['event_count']>0 and s['rebase_count']==1 and s['prediction_count']==2 and s['accelerometer_count']==2 and s['interior_only_accelerometer'] and s['all_prefix_identity_residuals_contain_zero'])
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,'synchronized_source_uniform_event_attachment_consumed':True,'cross_sample_Phi_rebase_consumed':True,'structural_interval_transport_consumed':True,'prediction_uses_trusted_shipping_F_not_nonlinear_Jacobian':True,'prediction_finite_attitude_defect_retained_in_rho':True,'Joseph_C_equals_G_I_minus_KH':True,'Joseph_L_equals_G':True,'Joseph_C_minus_L_built_symbolically_as_minus_GKH':True,'same_cell_K_H_R_state_used':True,'exact_quaternion_reset_defect_retained_in_rho':True,'full_mixed_epsilon_aw_evaluated_at_every_boundary':True,'source_change_epsilon_not_rezeroed':True,'only_accelerometer_has_interior_epsilon_transport':closed,'H18_every_prefix_transport_constructor_closed':closed,'smoke':s,'source_uniform_complete_601_sample_family_executed_here':False,'H18_endpoint_augmented_LDLT_closed_here':False,'H18_every_prefix_augmented_LDLT_closed_here':False,'H18_first_exit_retention_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,'next_obligation':'lift each H18 prefix decomposition defect onto the common source/radial/moment/binary32 graph coordinate and combine with the P3 strict information block; then run outward augmented LDLT and same-graph hard-ball/reset-domain targets'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('synchronized_source_uniform_event_attachment_consumed','cross_sample_Phi_rebase_consumed','structural_interval_transport_consumed','prediction_uses_trusted_shipping_F_not_nonlinear_Jacobian','prediction_finite_attitude_defect_retained_in_rho','Joseph_C_equals_G_I_minus_KH','Joseph_L_equals_G','Joseph_C_minus_L_built_symbolically_as_minus_GKH','same_cell_K_H_R_state_used','exact_quaternion_reset_defect_retained_in_rho','full_mixed_epsilon_aw_evaluated_at_every_boundary','source_change_epsilon_not_rezeroed','only_accelerometer_has_interior_epsilon_transport','H18_every_prefix_transport_constructor_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('source_uniform_complete_601_sample_family_executed_here','H18_endpoint_augmented_LDLT_closed_here','H18_every_prefix_augmented_LDLT_closed_here','H18_first_exit_retention_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'H18_transport':d['H18_every_prefix_transport_constructor_closed'],'smoke':d['smoke'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
