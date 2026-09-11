#!/usr/bin/env python3
"""A21 source-indexed every-subevent transport for BIAS0/BIAS1/BIAS2.

H18 already has the source-indexed Phi transport used by the every-prefix
stability proof.  A21 needs one extra structural ingredient: every shipping
Joseph/reset is immediately followed by the active accelerometer-bias estimate
projection.  The two operations must be separate proof prefixes because the
projection has a Clarke generalized Jacobian and depends on the SAME physical
true-bias source coordinate.

This module constructs, for one synchronized A21 selector lineage and one of the
three admitted bias families, the joint proof state

    x24 = [e_A21 ; b_true]

and inserts the following records in physical order:

  source-coordinate rebase (between IMU samples only)
  prediction with one new shared [w,m_tau] supply block
  optional a_w covariance floor
  smooth Joseph + finite quaternion reset
  active/inactive/Clarke bias projection

For prediction the 21-state shipping tangent F is captured by the trusted
Riccati event and lifted with the family-specific physical phi_true.  The same
w column enters e_b and b_true; m_tau enters e_b only.  For projection the
24-state map contains the exact beta coupling I-J_projection rather than
charging projection as an independent disturbance.

The source-indexed accelerometer linearizing coordinate is retained at every
boundary.  Projection leaves attitude and a_w unchanged, so its interior Phi
shift is exactly zero.  The output includes every prefix decomposition and the
cumulative suffix-propagated bias supply map.  It still does not attach the
BRMM acceleration-moment/radial sector, storage metric, or binary32 map and
therefore cannot promote P4.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval,matrix_identity,matrix_mul
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_p4_complete_brmm_accelerometer_operation_coordinate as ACC
import ou3_p4_complete_brmm_differential_prediction as PRED
import ou3_p4_source_uniform_estimator_event_attachment as ATTACH
import ou3_p4_source_uniform_event_lineage_sequence as LINEAGE
import ou3_p4_source_uniform_bias_prefix_lineage as BIAS
import ou3_p4_a21_bias1_24state_event_lift as LIFT
import ou3_p4_a21_joseph_projection_split as SPLIT
import ou3_p4_source_indexed_phi_rebase as REBASE
import ou3_p4_structural_prefix_transport as STRUCT

SCHEMA=1
QUALIFICATION='OU3_P4_A21_SOURCE_INDEXED_EVERY_SUBEVENT_TRANSPORT_V1'
P3_DELTA=1.0e-18
N=24;NE=21;OFF_AW=15;OFF_BA=18;OFF_BETA=21

def I(x):return Interval.point(float(x))
def zvec(n=N):return [I(0) for _ in range(n)]
def zmat(r,c):return [[I(0) for _ in range(c)] for _ in range(r)]
def ident(n):return [[I(1 if i==j else 0) for j in range(n)] for i in range(n)]
def shape(A):return len(A),len(A[0]) if A else 0

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
def hstack(A,B):
    if len(A)!=len(B):raise ValueError('horizontal stack row mismatch')
    return [list(a)+list(b) for a,b in zip(A,B)]
def Eaw24():
    E=zmat(N,3)
    for i in range(3):E[OFF_AW+i][i]=I(1)
    return E

def epsilon_aw(state21,selector):
    if len(state21)!=NE:raise ValueError('A21 state required')
    col=lambda v:[[x] for x in v]
    d=ACC.evaluate_operation_coordinate(col(state21[:3]),selector.sample_coordinates.R_wb,col(state21[OFF_AW:OFF_AW+3]),col(selector.sample_coordinates.f_cog_body))
    return [r[0] for r in d['epsilon_aw']]
def embed21(A):
    if shape(A)!=(NE,NE):raise ValueError('21x21 map required')
    M=ident(N)
    for i in range(NE):
        for j in range(NE):M[i][j]=A[i][j]
    return M
def augment_rho(r):
    if len(r)!=NE:raise ValueError('21-state defect required')
    return list(r)+[I(0),I(0),I(0)]
def projection24(Cp):
    if shape(Cp)!=(NE,NE):raise ValueError('21x21 projection map required')
    A=embed21(Cp)
    # e_b^+ = beta - Proj(beta-e_b); d/de_b=J, d/dbeta=I-J.
    for i in range(3):
        for j in range(3):
            A[OFF_BA+i][OFF_BETA+j]=I(1 if i==j else 0)-Cp[OFF_BA+i][OFF_BA+j]
    return A

def rebase24(e0,e1):
    r=REBASE.rebase_event('A',[[x] for x in e0],[[x] for x in e1])
    C=ident(N);rho=zvec(N)
    return {'kind':'source_coordinate_rebase','C':C,'L':C,'C_equals_L_exact':True,'C_minus_L':zmat(N,N),'rho':rho,'token':'rebase','Phi_xi_21':r['xi']}
def propagate_source(C,Bprev,Bnew=None):
    carried=matrix_mul(C,Bprev) if shape(Bprev)[1] else zmat(N,0)
    return carried if Bnew is None else hstack(carried,Bnew)

def trusted_map(sample):return {c.event_index_in_sample:c for c in sample.selector.A_event_cells}

def prediction_record(cell,trusted,selector,cert):
    if trusted.kind!='prediction' or trusted.F is None:raise ValueError('trusted A21 prediction F missing')
    p=PRED.prediction_event('A',cell.state,selector.sample_coordinates.omega_body_corrected,cell.dt_s,cell.tau_applied_s,tau_ba=KERNEL._process_constants().accel_bias_tau_s)
    A,B=LIFT.prediction_lift(trusted.F,cert.phi_true)
    rho21=vsub(p['state_out'],mv(trusted.F,cell.state))
    return {'kind':'prediction','C':A,'L':A,'C_equals_L_exact':True,'C_minus_L':zmat(N,N),'rho':augment_rho(rho21),'token':cell.source_token,'prediction_supply_block':B},list(p['state_out'])
def floor_record(cell):
    J=ident(N);return {'kind':'aw_floor','C':J,'L':J,'C_equals_L_exact':True,'C_minus_L':zmat(N,N),'rho':zvec(),'token':cell.source_token},list(cell.state)
def split_records(cell):
    s=SPLIT.split_event(cell)
    sm=s['smooth'];pr=s['projection']
    Cj=embed21(sm['C']);Lj=embed21(sm['L']);rj=augment_rho(sm['rho'])
    Cp=projection24(pr['C']);rp=vsub(list(pr['state_out'])+list(cell.true_bias),mv(Cp,list(sm['state_out'])+list(cell.true_bias)))
    smooth={'kind':sm['kind'],'C':Cj,'L':Lj,'C_minus_L':embed21(sm['C_minus_L']),'rho':rj,'token':cell.source_token+':smooth'}
    # embed21 would leave beta identity in C-L; overwrite lower identity to zero.
    smooth['C_minus_L']=zmat(N,N)
    for i in range(NE):
        for j in range(NE):smooth['C_minus_L'][i][j]=sm['C_minus_L'][i][j]
    proj={'kind':'bias_projection','C':Cp,'L':Cp,'C_equals_L_exact':True,'C_minus_L':zmat(N,N),'rho':rp,'token':cell.source_token+':projection','projection_branch':pr['branch']}
    return smooth,proj,list(sm['state_out']),list(pr['state_out'])

def build_transport(samples:Sequence[ATTACH.AttachedSampleLineage],family:str):
    if family not in BIAS.FAMILIES:raise ValueError('family must be BIAS0/BIAS1/BIAS2')
    stitched=LINEAGE.stitch_samples(samples)
    if stitched.mode!='A':raise ValueError('A21 lineage required')
    selectors=[s.selector for s in samples]
    cert=BIAS.attach_to_selector_lineage(family,selectors)
    for s in samples:
        expected=cert.at(s.selector.source_cell_id)
        for c in s.cells:
            if c.true_bias!=expected:raise RuntimeError('A21 event true-bias cell detached from family recurrence prefix')
    events=[];eps=[];emb=[];source_prefixes=[];state=list(samples[0].cells[0].state);Bcur=zmat(N,0)
    eps.append(epsilon_aw(state,samples[0].selector));emb.append(Eaw24())
    for si,sample in enumerate(samples):
        selector=sample.selector
        if si:
            neweps=epsilon_aw(state,selector);rec=rebase24(eps[-1],neweps);rec['token']=selector.source_cell_id+':rebase'
            events.append(rec);Bcur=propagate_source(rec['C'],Bcur);source_prefixes.append(copy.deepcopy(Bcur));eps.append(neweps);emb.append(Eaw24())
        tm=trusted_map(sample)
        for cell in sample.cells:
            if tuple(cell.state)!=tuple(state):raise RuntimeError('A21 event state detached from previous physical boundary')
            if cell.kind=='prediction':
                rec,out=prediction_record(cell,tm[cell.event_ordinal],selector,cert);Bcur=propagate_source(rec['C'],Bcur,rec['prediction_supply_block']);events.append(rec);state=list(out);eps.append(epsilon_aw(state,selector));emb.append(Eaw24());source_prefixes.append(copy.deepcopy(Bcur))
            elif cell.kind=='aw_floor':
                rec,out=floor_record(cell);Bcur=propagate_source(rec['C'],Bcur);events.append(rec);state=list(out);eps.append(epsilon_aw(state,selector));emb.append(Eaw24());source_prefixes.append(copy.deepcopy(Bcur))
            else:
                smooth,proj,sout,pout=split_records(cell)
                Bcur=propagate_source(smooth['C'],Bcur);events.append(smooth);state=list(sout);eps.append(epsilon_aw(state,selector));emb.append(Eaw24());source_prefixes.append(copy.deepcopy(Bcur))
                before_proj_eps=eps[-1]
                Bcur=propagate_source(proj['C'],Bcur);events.append(proj);state=list(pout);eps.append(epsilon_aw(state,selector));emb.append(Eaw24());source_prefixes.append(copy.deepcopy(Bcur))
                if eps[-1]!=before_proj_eps:raise RuntimeError('bias projection changed attitude/a_w Phi coordinate')
        if tuple(state)!=tuple(sample.state_out):raise RuntimeError('A21 sample endpoint detached from synchronized nonlinear chain')
    decomp=STRUCT.prefix_decompositions(events,emb,eps)
    if len(decomp)!=len(events) or len(source_prefixes)!=len(events):raise RuntimeError('A21 prefix transport/source-map count mismatch')
    identity_ok=all(r['transport']['identity_residual_contains_zero'] for r in decomp)
    projection_zero_phi=all(eps[i+1]==eps[i] for i,e in enumerate(events) if e['kind']=='bias_projection')
    supply_dims=[shape(B)[1] for B in source_prefixes]
    expected=0;dim_ok=True
    for e,d in zip(events,supply_dims):
        if e['kind']=='prediction':expected+=6
        if d!=expected:dim_ok=False
    return {'events':events,'decompositions':decomp,'source_prefix_maps':source_prefixes,'state_out':tuple(state),'bias_family':family,'bias_qualification':BIAS.family_parameters(family)['qualification'],'all_identity_residuals_contain_zero':identity_ok,'projection_has_zero_Phi_interior_shift':projection_zero_phi,'one_6D_shared_w_mtau_block_per_prediction':dim_ok,'prediction_supply_blocks_same_w':all(LIFT.source_map_same_w(e['prediction_supply_block']) for e in events if e['kind']=='prediction')}

def point_sample(forcez=-9.80665):
    return KERNEL.SampleCoordinates(gyro_measurement=KERNEL.MAHONY.Vec3(I(.01),I(-.02),I(.005)),omega_body_corrected=(I(.01),I(-.02),I(.005)),specific_force=KERNEL.MAHONY.Vec3(I(.2),I(-.1),I(forcez)),f_cog_body=(I(.2),I(-.1),I(forcez)),R_wb=ident(3),due_S=True,aw_floor_requested=True,magnetometer_events_after_imu=(KERNEL.MagneticEvent((I(20),I(0),I(40))),))
def smoke_objects(family):
    js=ATTACH.JOINT._smoke_state();branch=KERNEL.ExecutionBranch(frontend=copy.deepcopy(js.frontend),H=ATTACH.WORD.initialize_word('H',ident(18)),A=ATTACH.WORD.initialize_word('A',ident(21)),source_cell_id='root')
    boxes=BIAS.prefix_component_boxes(family,2);common=dict(radial_scale=Interval(0,1),bias_projection_limit=.4,tau_ba=I(1800))
    b0=[boxes[0],boxes[0],boxes[0]];p0=ATTACH.synchronize_sample(branch=branch,joint_state=js,sample=point_sample(),state_in_H=[I(0) for _ in range(18)],state_in_A=[I(0) for _ in range(21)],true_bias=b0,sample_index=0,next_cell_prefix=family+'-a21-k0',**common)
    if not p0:raise RuntimeError('sample0 empty')
    h0,a0=p0[0];branch1=ATTACH.next_execution_branch(h0,a0);b1=[boxes[1],boxes[1],boxes[1]]
    p1=ATTACH.synchronize_sample(branch=branch1,joint_state=h0.image.state,sample=point_sample(-9.7),state_in_H=h0.state_out,state_in_A=a0.state_out,true_bias=b1,sample_index=1,next_cell_prefix=family+'-a21-k1',**common)
    if not p1:raise RuntimeError('sample1 empty')
    _,a1=p1[0];return a0,a1
def smoke(family):
    t=build_transport(smoke_objects(family),family);k=[e['kind'] for e in t['events']]
    return {'event_count':len(k),'prediction_count':k.count('prediction'),'smooth_measurement_count':sum(x.endswith('_smooth') for x in k),'projection_count':k.count('bias_projection'),'source_rebase_count':k.count('source_coordinate_rebase'),'all_identity_residuals_contain_zero':t['all_identity_residuals_contain_zero'],'projection_has_zero_Phi_interior_shift':t['projection_has_zero_Phi_interior_shift'],'one_6D_shared_w_mtau_block_per_prediction':t['one_6D_shared_w_mtau_block_per_prediction'],'prediction_supply_blocks_same_w':t['prediction_supply_blocks_same_w']}
def build():
    line=LINEAGE.build();lf=LINEAGE.validate(line);bias=BIAS.build();bf=BIAS.validate(bias);split=SPLIT.build();sf=SPLIT.validate(split);lift=LIFT.build();ltf=LIFT.validate(lift);struct=STRUCT.build();stf=STRUCT.validate(struct)
    if lf or bf or sf or ltf or stf:raise RuntimeError(f'A21 source-indexed prerequisites failed lineage={lf} bias={bf} split={sf} lift={ltf} structural={stf}')
    rows={f:smoke(f) for f in BIAS.FAMILIES};closed=all(r['prediction_count']==2 and r['projection_count']>=2 and r['source_rebase_count']==1 and r['all_identity_residuals_contain_zero'] and r['projection_has_zero_Phi_interior_shift'] and r['one_6D_shared_w_mtau_block_per_prediction'] and r['prediction_supply_blocks_same_w'] for r in rows.values())
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','P3_delta':P3_DELTA,'all_BIAS0_BIAS1_BIAS2_recurrences_composed':closed,'A21_Joseph_reset_and_projection_are_distinct_prefixes':True,'source_indexed_Phi_rebase_inserted_between_samples':True,'projection_has_exact_zero_Phi_interior_shift':closed,'persistent_true_bias_coordinate_dimension':3,'joint_state_dimension':24,'one_shared_w_and_mtau_supply_block_per_prediction':closed,'same_w_enters_bias_error_and_true_bias':closed,'family_specific_phi_true_retained':True,'projection_true_bias_coupling_retained_in_24state_map':True,'A21_every_subevent_prefix_transport_constructor_closed':closed,'smoke_by_family':rows,'independent_projection_disturbance_port_used':False,'independent_per_sample_true_bias_boxes_used':False,'finite_source_enumeration_used':False,'production_BRMM_moment_radial_sector_attached_here':False,'production_binary32_prefix_map_attached_here':False,'production_compatible_storage_attached_here':False,'A21_endpoint_augmented_LDLT_closed_here':False,'A21_every_prefix_augmented_LDLT_closed_here':False,'A21_first_exit_retention_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,'next_obligation':'combine H18 and this A21 source-indexed prefix transport with source-dependent compatible storage; attach the SAME primitive moment/radial sector and conditional binary32 map to each cumulative source map, then run outward endpoint/every-prefix LDLT'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('P3_delta')!=P3_DELTA:f.append('P3 delta changed')
    for k in ('all_BIAS0_BIAS1_BIAS2_recurrences_composed','A21_Joseph_reset_and_projection_are_distinct_prefixes','source_indexed_Phi_rebase_inserted_between_samples','projection_has_exact_zero_Phi_interior_shift','one_shared_w_and_mtau_supply_block_per_prediction','same_w_enters_bias_error_and_true_bias','family_specific_phi_true_retained','projection_true_bias_coupling_retained_in_24state_map','A21_every_subevent_prefix_transport_constructor_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('independent_projection_disturbance_port_used','independent_per_sample_true_bias_boxes_used','finite_source_enumeration_used','production_BRMM_moment_radial_sector_attached_here','production_binary32_prefix_map_attached_here','production_compatible_storage_attached_here','A21_endpoint_augmented_LDLT_closed_here','A21_every_prefix_augmented_LDLT_closed_here','A21_first_exit_retention_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('joint_state_dimension')!=24 or d.get('persistent_true_bias_coordinate_dimension')!=3:f.append('joint dimensions changed')
    if set(d.get('smoke_by_family',{}))!=set(BIAS.FAMILIES):f.append('bias family smoke set incomplete')
    return list(dict.fromkeys(f))
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'A21_transport':d['A21_every_subevent_prefix_transport_constructor_closed'],'families':list(d['smoke_by_family']),'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
