#!/usr/bin/env python3
"""Real same-cell reduced signed Joseph/reset master for BRMM P4 prefixes.

For one estimator-owned Joseph source cell derive from the SAME cell

    J=P^-1, H, R^-1, K=P H' (H P H' + R)^-1.

Then substitute q=H e+eta and d=Kq algebraically into the exact reduced
Joseph/reset identity.  No free K/q/d box is present.

The nonlinear residual is also kept structurally instead of as an independent
eta port:

* S=0: eta=0 exactly, so q=H e.
* magnetometer/vector: lift exact chord p=Cm, q_c=(E-I)m and use
      eta = q_c-p,
  with p=Cm as an exact linear equality and the exact Cayley chord IQCs.
* accelerometer: d_a=R_hat*delta_a_w, r=C d_a,
      p=C(f_hat+d_a)=C f_hat+r,
      eta=q_c-p+r,
  with p=C f_hat+r exact and r tied jointly to attitude and a_w by the retained
  cross-product IQCs. This keeps the c-a_w mixed term and b_a linear companion
  inside the same residual graph.

The finite covariance-reset defect is attitude-only: b_theta is lifted and the
full b map is zero outside the first three state rows. A homogeneous h=1
coordinate supports the full hard-entry IQCs and the same-Dtheta reset
correction-domain target.
"""
from __future__ import annotations
from dataclasses import dataclass
import argparse,json,math
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval,matrix_add,matrix_mul,matrix_sub,matrix_transpose
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan,matrix_symmetric_hull
import ou3_p4_complete_brmm_source_cover_contract as COVER
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_reduced_joseph_reset_event_quadratic as REDUCED
import ou3_p4_reset_graph_iqc as RESETIQC
import ou3_p4_affine_hard_tube_iqc as HARD
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_complete_brmm_exact_chord_iqc as CHORDIQC
import ou3_p4_complete_brmm_exact_chord_joint_coordinate as CHORD
import ou3_p4_joint_brmm_frontend_transition as JOINT
import ou3_brmm_private_mahony_state_step as MAHONY
import ou3_brmm_frontend_state_step as FRONT

SCHEMA=3
QUALIFICATION='OU3_P4_BRMM_REDUCED_SAME_CELL_EVENT_MASTER_V3'
OFF_AW=15

def I(x):return Interval.point(float(x))
def shape(A):
    r=len(A);c=len(A[0]) if r else 0
    if any(len(row)!=c for row in A):raise ValueError('ragged matrix')
    return r,c
def zeros(r,c):return [[I(0) for _ in range(c)] for _ in range(r)]
def selector(rows,nz,offset):
    if offset<0 or offset+rows>nz:raise ValueError('selector outside coordinate')
    return [[I(1 if j==offset+i else 0) for j in range(nz)] for i in range(rows)]
def symmetric(A):
    n,m=shape(A);return n==m and all(A[i][j].lo==A[j][i].lo and A[i][j].hi==A[j][i].hi for i in range(n) for j in range(n))
def gram(A):return matrix_mul(matrix_transpose(A),A)
def neg(A):return [[-x for x in row] for row in A]
def mscale(A,a):
    c=I(a);return [[c*x for x in row] for row in A]
def maxabs(x):return max(abs(x.lo),abs(x.hi))

def _linear_cross_c_const(v,nz):
    if len(v)!=3:raise ValueError('constant vector must have length three')
    x,y,z=v;A=zeros(3,nz)
    A[0][1]=z;A[0][2]=-y
    A[1][0]=-z;A[1][2]=x
    A[2][0]=y;A[2][1]=-x
    return A

def _embed_rotation_aw(Rhat,nz):
    if len(Rhat)!=3 or any(len(row)!=3 for row in Rhat):raise ValueError('Rhat must be 3x3')
    S=selector(3,nz,OFF_AW)
    return matrix_mul(Rhat,S)

def _equality_zero_sector(A):
    return matrix_symmetric_hull(neg(gram(A)))

@dataclass(frozen=True)
class EventMaster:
    source_token:str;estimator_source_token:str;mode:str;kind:str
    coordinate_dimension:int;coordinate_layout:tuple[tuple[str,int,int],...];h_index:int
    master:Sequence[Sequence[Interval]];H:Sequence[Sequence[Interval]];K:Sequence[Sequence[Interval]]
    J:Sequence[Sequence[Interval]];R_inverse:Sequence[Sequence[Interval]]
    E:Sequence[Sequence[Interval]];Q:Sequence[Sequence[Interval]];N:Sequence[Sequence[Interval]]
    D:Sequence[Sequence[Interval]];D_theta:Sequence[Sequence[Interval]]
    B:Sequence[Sequence[Interval]];B_theta:Sequence[Sequence[Interval]]
    nonlinear_sectors:tuple[tuple[str,Sequence[Sequence[Interval]]],...]
    reset_b_constraint:str

def _nonlinear_maps(cell,event,n):
    entry=ENTRY.build();ef=ENTRY.validate(entry)
    if ef:raise RuntimeError('hard-entry prerequisite failed: '+repr(ef))
    chord=CHORD.build();cf=CHORD.validate(chord)
    if cf:raise RuntimeError('chord prerequisite failed: '+repr(cf))
    k=float(chord['information_retention_factor_lower_full_entry'])
    cmax=float(entry['coordinate_radii']['attitude_cayley_norm'])
    awmax=float(entry['coordinate_radii']['latent_acceleration_norm_mps2'])

    if cell.kind=='S_zero':
        nz=n+4;hidx=n+3;N=zeros(3,nz);layout=(('e',0,n),('b_theta',n,3),('h',hidx,1));return nz,hidx,layout,N,(),n

    if cell.kind=='magnetometer':
        p0=n;qc0=n+3;b0=n+6;hidx=n+9;nz=n+10
        P=selector(3,nz,p0);Qc=selector(3,nz,qc0)
        Cm=_linear_cross_c_const(cell.m_body,nz)
        peq=matrix_sub(P,Cm)
        N=matrix_sub(Qc,P)
        sectors=[('p_equals_Cm',_equality_zero_sector(peq))]
        sectors.extend(('chord_'+name,Pi) for name,Pi in CHORDIQC.exact_chord_iqcs(P,Qc,k).items())
        layout=(('e',0,n),('p',p0,3),('q_chord',qc0,3),('b_theta',b0,3),('h',hidx,1))
        return nz,hidx,layout,N,tuple(sectors),b0

    if cell.kind=='accelerometer':
        p0=n;qc0=n+3;r0=n+6;b0=n+9;hidx=n+12;nz=n+13
        P=selector(3,nz,p0);Qc=selector(3,nz,qc0);R=selector(3,nz,r0)
        Cf=_linear_cross_c_const(cell.f_hat,nz);Daw=_embed_rotation_aw(cell.R_hat,nz);Cmap=selector(3,nz,0)
        peq=matrix_sub(matrix_sub(P,Cf),R)
        N=matrix_add(matrix_sub(Qc,P),R)
        sectors=[('p_equals_Cf_plus_r',_equality_zero_sector(peq))]
        sectors.extend(('chord_'+name,Pi) for name,Pi in CHORDIQC.exact_chord_iqcs(P,Qc,k).items())
        sectors.extend(('cross_'+name,Pi) for name,Pi in CHORDIQC.cross_product_iqcs(Cmap,Daw,R,cmax,awmax).items())
        layout=(('e',0,n),('p',p0,3),('q_chord',qc0,3),('r_c_cross_aw',r0,3),('b_theta',b0,3),('h',hidx,1))
        return nz,hidx,layout,N,tuple(sectors),b0
    raise ValueError('unsupported Joseph kind')

def build_event_master(cell:COVER.SourceCoverCell)->EventMaster:
    failures=COVER.validate_cell(cell,require_estimator_provenance=True)
    if failures:raise ValueError('invalid estimator-owned source cell: '+repr(failures))
    if cell.kind not in ('S_zero','accelerometer','magnetometer'):raise ValueError('reduced event master requires Joseph event')
    event=EVENTS.source_joseph_event(**COVER.joseph_event_kwargs(cell));n=18 if cell.mode=='H' else 21
    J=matrix_inverse_gauss_jordan(cell.P);Rinv=matrix_inverse_gauss_jordan(cell.R)
    nz,hidx,layout,N,sectors,b0=_nonlinear_maps(cell,event,n)
    E=selector(n,nz,0);Btheta=selector(3,nz,b0);B=zeros(n,nz)
    for i in range(3):B[i]=list(Btheta[i])
    Q=matrix_add(matrix_mul(event['H'],E),N);D=matrix_mul(event['K'],Q);Dtheta=[list(row) for row in D[:3]]
    M=matrix_symmetric_hull(REDUCED.event_quadratic(J,event['H'],Rinv,E,Q,N,D,B))
    return EventMaster(cell.source_token,cell.estimator_source_token or '',cell.mode,cell.kind,nz,layout,hidx,M,event['H'],event['K'],J,Rinv,E,Q,N,D,Dtheta,B,Btheta,sectors,'b_theta=G^-1*rho_theta exact reset graph sector')

def reset_graph_bundle(em,delta):
    entry=ENTRY.build();ef=ENTRY.validate(entry)
    if ef:raise RuntimeError('hard-entry prerequisite failed: '+repr(ef))
    q=float(entry['coordinate_radii']['attitude_cayley_norm'])
    Pi,sector=RESETIQC.parameterized_reset_iqc(em.D_theta,em.B_theta,q,float(delta))
    target=RESETIQC.correction_domain_target(em.D_theta,em.h_index,float(delta))
    hard=HARD.hard_entry_iqcs('H18' if em.mode=='H' else 'A21',em.coordinate_dimension,h_index=em.h_index,state_offset=0)
    return {'reset_sector':Pi,'reset_sector_contract':sector,'correction_domain_target':target,'hard_entry_iqcs':hard,'delta':float(delta),'delta_certified':False}

def validate_event_master(cell,em):
    f=[];n=18 if cell.mode=='H' else 21;nz=em.coordinate_dimension
    if em.source_token!=cell.source_token:f.append('event token detached')
    if em.estimator_source_token!=cell.estimator_source_token:f.append('estimator token detached')
    if shape(em.master)!=(nz,nz) or not symmetric(em.master):f.append('master dimension/symmetry mismatch')
    if shape(em.H)!=(3,n) or shape(em.K)!=(n,3):f.append('same-cell H/K dimensions invalid')
    if shape(em.J)!=(n,n) or shape(em.R_inverse)!=(3,3):f.append('precision dimensions invalid')
    if shape(em.Q)!=(3,nz) or shape(em.D)!=(n,nz) or shape(em.D_theta)!=(3,nz):f.append('derived q/d maps invalid')
    if shape(em.B)!=(n,nz) or shape(em.B_theta)!=(3,nz):f.append('reset b maps invalid')
    if any(not all(x.lo<=0<=x.hi for x in em.B[i]) for i in range(3,n)):f.append('non-attitude reset defect rows not zero')
    if any(shape(Pi)!=(nz,nz) or not symmetric(Pi) for _,Pi in em.nonlinear_sectors):f.append('nonlinear sector dimension/symmetry mismatch')
    if cell.kind=='S_zero' and em.nonlinear_sectors:f.append('S=0 should eliminate eta without nonlinear auxiliary sectors')
    if cell.kind!='S_zero' and not em.nonlinear_sectors:f.append('nonlinear event sectors missing')
    return f

def _smoke_image():
    st=JOINT._smoke_state();sample=FRONT.Sample(MAHONY.Vec3(MAHONY.I(.01),MAHONY.I(-.02),MAHONY.I(.005)),MAHONY.Vec3(MAHONY.I(.2),MAHONY.I(-.1),MAHONY.I(-9.75)))
    images=JOINT.advance(st,sample,gravity_ms2=MAHONY.I(9.80665),two_kp=MAHONY.I(.2),two_ki=MAHONY.I(.02),child_prefix='master')
    if not images:raise RuntimeError('joint image smoke empty')
    return images[0]
def _identity(n):return [[I(1 if i==j else 0) for j in range(n)] for i in range(n)]
def _accel_smoke_cell(im):
    n=18;x=[I(0) for _ in range(n)];P=_identity(n);R=_identity(3);Rhat=_identity(3);f=[I(.2),I(-.1),I(-9.7)]
    return COVER.source_cell_from_joint_image(im,mode='H',sample_index=0,event_ordinal=3,kind='accelerometer',state=x,P=P,dt_s=I(.005),pseudo_elapsed_s=I(.01),radial_scale=Interval(0,1),event_source_token=im.source_token+':e3',event_predecessor_token=im.source_token+':e2',R=R,f_hat=f,R_hat=Rhat)

def build():
    prereq={'reduced':REDUCED.validate(REDUCED.build()),'reset':RESETIQC.validate(RESETIQC.build()),'hard':HARD.validate(HARD.build()),'chord_iqc':CHORDIQC.validate(CHORDIQC.build())};bad={k:v for k,v in prereq.items() if v}
    if bad:raise RuntimeError('event-master prerequisites failed: '+repr(bad))
    im=_smoke_image();n=18;x=[I(0) for _ in range(n)];P=_identity(n)
    sz=COVER.source_cell_from_joint_image(im,mode='H',sample_index=0,event_ordinal=2,kind='S_zero',state=x,P=P,dt_s=I(.005),pseudo_elapsed_s=I(.01),radial_scale=Interval(0,1),event_source_token=im.source_token+':e2',event_predecessor_token=im.source_token+':e1')
    acc=_accel_smoke_cell(im)
    ems=build_event_master(sz);ema=build_event_master(acc);efs=validate_event_master(sz,ems);efa=validate_event_master(acc,ema);bundle=reset_graph_bundle(ems,.1)
    dims=shape(bundle['reset_sector'])==(ems.coordinate_dimension,ems.coordinate_dimension) and shape(bundle['correction_domain_target'])==(ems.coordinate_dimension,ems.coordinate_dimension) and all(shape(v)==(ems.coordinate_dimension,ems.coordinate_dimension) for v in bundle['hard_entry_iqcs'].values())
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'same_cell_P_H_R_K_consumed':not(efs or efa),'prior_precision_derived_from_same_P':True,'R_inverse_derived_from_same_R':True,
      'q_equals_H_e_plus_structured_eta_substituted_exactly':True,'d_equals_same_cell_K_q_substituted_exactly':True,
      'independent_K_coordinate_or_box_used':False,'independent_q_port_used':False,'independent_d_port_used':False,'independent_eta_port_used':False,
      'S_zero_eta_eliminated_exactly':len(ems.nonlinear_sectors)==0,
      'accelerometer_eta_eliminated_into_chord_and_cross_coordinates':len(ema.nonlinear_sectors)>=11,
      'accelerometer_p_equals_Cf_plus_r_exact_equality_sector':any(n=='p_equals_Cf_plus_r' for n,_ in ema.nonlinear_sectors),
      'accelerometer_exact_chord_IQCs_attached':sum(n.startswith('chord_') for n,_ in ema.nonlinear_sectors)==4,
      'accelerometer_c_cross_aw_IQCs_attached':sum(n.startswith('cross_') for n,_ in ema.nonlinear_sectors)==6,
      'reset_defect_attitude_only_coordinate_used':True,'full_state_independent_reset_b_port_used':False,
      'homogeneous_h_coordinate_materialized':True,'same_Dtheta_reset_sector_materialized':dims,'same_Dtheta_correction_domain_target_materialized':dims,
      'full_declared_hard_entry_IQCs_on_same_event_coordinate_materialized':dims,
      'smoke_reset_delta_is_diagnostic_not_certified':True,'smoke_reset_delta_certified':bundle['delta_certified'],
      'real_reduced_event_master_materialized':not(efs or efa),'S_zero_event_master_dimension':ems.coordinate_dimension,'accelerometer_event_master_dimension':ema.coordinate_dimension,
      'event_master_validation_failures':efs+efa,'production_same_graph_reset_delta_certified_here':False,
      'production_BIAS1_source_map_attached_here':False,'production_binary32_fp_map_attached_here':False,
      'production_word_prefix_master_accumulation_closed_here':False,'production_endpoint_LDLT_closed_here':False,'production_every_prefix_LDLT_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'certify the same-Dtheta correction-domain target from hard/chord/source premises; embed event masters/sectors into a lineage-wide coordinate, attach physical BIAS1 and binary32 maps, and accumulate every literal prefix'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('same_cell_P_H_R_K_consumed','prior_precision_derived_from_same_P','R_inverse_derived_from_same_R','q_equals_H_e_plus_structured_eta_substituted_exactly','d_equals_same_cell_K_q_substituted_exactly','S_zero_eta_eliminated_exactly','accelerometer_eta_eliminated_into_chord_and_cross_coordinates','accelerometer_p_equals_Cf_plus_r_exact_equality_sector','accelerometer_exact_chord_IQCs_attached','accelerometer_c_cross_aw_IQCs_attached','reset_defect_attitude_only_coordinate_used','homogeneous_h_coordinate_materialized','same_Dtheta_reset_sector_materialized','same_Dtheta_correction_domain_target_materialized','full_declared_hard_entry_IQCs_on_same_event_coordinate_materialized','smoke_reset_delta_is_diagnostic_not_certified','real_reduced_event_master_materialized'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('independent_K_coordinate_or_box_used','independent_q_port_used','independent_d_port_used','independent_eta_port_used','full_state_independent_reset_b_port_used','smoke_reset_delta_certified','production_same_graph_reset_delta_certified_here','production_BIAS1_source_map_attached_here','production_binary32_fp_map_attached_here','production_word_prefix_master_accumulation_closed_here','production_endpoint_LDLT_closed_here','production_every_prefix_LDLT_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('event_master_validation_failures')!=[]:f.append('event master smoke validation failed')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'event_master':d['real_reduced_event_master_materialized'],'accel_chord':d['accelerometer_eta_eliminated_into_chord_and_cross_coordinates'],'reset_sector':d['same_Dtheta_reset_sector_materialized'],'delta_closed':d['production_same_graph_reset_delta_certified_here'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())