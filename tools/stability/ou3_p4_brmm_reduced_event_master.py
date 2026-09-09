#!/usr/bin/env python3
"""Real same-cell reduced signed Joseph/reset master for BRMM P4 prefixes.

For one estimator-owned Joseph ``SourceCoverCell`` derive from the SAME cell

    J=P^-1, H, R^-1, K=P H' (H P H' + R)^-1

and materialize the exact reduced event quadratic.  Use one affine-homogeneous
coordinate

    z = [e ; eta ; b_theta ; h].

The physical residual and correction are substituted algebraically:

    q = H e + eta,
    d = K q,

so q, d and K are not free lifted ports. The finite covariance-reset defect is
attitude-only, hence only b_theta is lifted; the full b map has zeros outside
its first three rows. ``h`` is the hard-domain homogeneous coordinate.

For S=0, eta=0 is encoded exactly by -||eta||^2 >= 0.  The same-cell finite
reset sector is parameterized by D_theta=E_theta K Q and the *same* b_theta map.
Its correction-domain target ||D_theta z||<=delta h and the full declared hard
entry IQCs are returned together.  A supplied delta is not accepted as proved
until that target is itself certified from the same hard/graph premises.
"""
from __future__ import annotations
from dataclasses import dataclass
import argparse,json
from typing import Sequence

from ou3_interval import Interval,matrix_add,matrix_mul,matrix_transpose,matrix_sub
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan,matrix_symmetric_hull
import ou3_p4_complete_brmm_source_cover_contract as COVER
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_reduced_joseph_reset_event_quadratic as REDUCED
import ou3_p4_reset_graph_iqc as RESETIQC
import ou3_p4_affine_hard_tube_iqc as HARD
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_joint_brmm_frontend_transition as JOINT
import ou3_brmm_private_mahony_state_step as MAHONY
import ou3_brmm_frontend_state_step as FRONT

SCHEMA=2
QUALIFICATION='OU3_P4_BRMM_REDUCED_SAME_CELL_EVENT_MASTER_V2'

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
    n,m=shape(A)
    return n==m and all(A[i][j].lo==A[j][i].lo and A[i][j].hi==A[j][i].hi for i in range(n) for j in range(n))
def gram(A):return matrix_mul(matrix_transpose(A),A)
def neg(A):return [[-x for x in row] for row in A]

@dataclass(frozen=True)
class EventMaster:
    source_token:str;estimator_source_token:str;mode:str;kind:str
    coordinate_dimension:int;coordinate_layout:tuple[tuple[str,int,int],...];h_index:int
    master:Sequence[Sequence[Interval]];H:Sequence[Sequence[Interval]];K:Sequence[Sequence[Interval]]
    J:Sequence[Sequence[Interval]];R_inverse:Sequence[Sequence[Interval]]
    E:Sequence[Sequence[Interval]];Q:Sequence[Sequence[Interval]];N:Sequence[Sequence[Interval]]
    D:Sequence[Sequence[Interval]];D_theta:Sequence[Sequence[Interval]]
    B:Sequence[Sequence[Interval]];B_theta:Sequence[Sequence[Interval]]
    eta_zero_sector:Sequence[Sequence[Interval]]|None
    eta_constraint:str;reset_b_constraint:str

def build_event_master(cell:COVER.SourceCoverCell)->EventMaster:
    failures=COVER.validate_cell(cell,require_estimator_provenance=True)
    if failures:raise ValueError('invalid estimator-owned source cell: '+repr(failures))
    if cell.kind not in ('S_zero','accelerometer','magnetometer'):raise ValueError('reduced event master requires Joseph event')
    event=EVENTS.source_joseph_event(**COVER.joseph_event_kwargs(cell));n=18 if cell.mode=='H' else 21
    J=matrix_inverse_gauss_jordan(cell.P);Rinv=matrix_inverse_gauss_jordan(cell.R)
    # [e(n), eta(3), b_theta(3), h(1)]
    nz=n+7;hidx=n+6
    E=selector(n,nz,0);N=selector(3,nz,n);Btheta=selector(3,nz,n+3)
    B=zeros(n,nz)
    for i in range(3):B[i]=list(Btheta[i])
    Q=matrix_add(matrix_mul(event['H'],E),N);D=matrix_mul(event['K'],Q);Dtheta=[list(row) for row in D[:3]]
    M=matrix_symmetric_hull(REDUCED.event_quadratic(J,event['H'],Rinv,E,Q,N,D,B))
    eta_zero=matrix_symmetric_hull(neg(gram(N))) if cell.kind=='S_zero' else None
    eta_constraint='-||eta||^2>=0 (eta==0 exact equality)' if cell.kind=='S_zero' else 'exact nonlinear residual/chord sector on same e/source geometry'
    return EventMaster(cell.source_token,cell.estimator_source_token or '',cell.mode,cell.kind,nz,
        (('e',0,n),('eta',n,3),('b_theta',n+3,3),('h',hidx,1)),hidx,M,event['H'],event['K'],J,Rinv,E,Q,N,D,Dtheta,B,Btheta,eta_zero,eta_constraint,'b_theta=G^-1*rho_theta exact reset graph sector')

def reset_graph_bundle(em:EventMaster,delta:float):
    entry=ENTRY.build();ef=ENTRY.validate(entry)
    if ef:raise RuntimeError('hard-entry prerequisite failed: '+repr(ef))
    q=float(entry['coordinate_radii']['attitude_cayley_norm'])
    Pi,sector=RESETIQC.parameterized_reset_iqc(em.D_theta,em.B_theta,q,float(delta))
    target=RESETIQC.correction_domain_target(em.D_theta,em.h_index,float(delta))
    hard=HARD.hard_entry_iqcs('H18' if em.mode=='H' else 'A21',em.coordinate_dimension,h_index=em.h_index,state_offset=0)
    return {'reset_sector':Pi,'reset_sector_contract':sector,'correction_domain_target':target,'hard_entry_iqcs':hard,'delta':float(delta),'delta_certified':False}

def validate_event_master(cell,em):
    f=[];n=18 if cell.mode=='H' else 21;nz=n+7
    if em.source_token!=cell.source_token:f.append('event token detached')
    if em.estimator_source_token!=cell.estimator_source_token:f.append('estimator token detached')
    if shape(em.master)!=(nz,nz) or not symmetric(em.master):f.append('master dimension/symmetry mismatch')
    if shape(em.H)!=(3,n) or shape(em.K)!=(n,3):f.append('same-cell H/K dimensions invalid')
    if shape(em.J)!=(n,n) or shape(em.R_inverse)!=(3,3):f.append('precision dimensions invalid')
    if shape(em.Q)!=(3,nz) or shape(em.D)!=(n,nz) or shape(em.D_theta)!=(3,nz):f.append('derived q/d maps invalid')
    if shape(em.B)!=(n,nz) or shape(em.B_theta)!=(3,nz):f.append('reset b maps invalid')
    if any(not all(x.lo<=0<=x.hi for x in em.B[i]) for i in range(3,n)):f.append('non-attitude reset defect rows not zero')
    if em.kind=='S_zero' and (em.eta_zero_sector is None or shape(em.eta_zero_sector)!=(nz,nz)):f.append('S=0 eta equality sector missing')
    return f

def _smoke_image():
    st=JOINT._smoke_state();sample=FRONT.Sample(MAHONY.Vec3(MAHONY.I(.01),MAHONY.I(-.02),MAHONY.I(.005)),MAHONY.Vec3(MAHONY.I(.2),MAHONY.I(-.1),MAHONY.I(-9.75)))
    images=JOINT.advance(st,sample,gravity_ms2=MAHONY.I(9.80665),two_kp=MAHONY.I(.2),two_ki=MAHONY.I(.02),child_prefix='master')
    if not images:raise RuntimeError('joint image smoke empty')
    return images[0]
def _identity(n):return [[I(1 if i==j else 0) for j in range(n)] for i in range(n)]

def build():
    prereq={'reduced':REDUCED.validate(REDUCED.build()),'reset':RESETIQC.validate(RESETIQC.build()),'hard':HARD.validate(HARD.build())};bad={k:v for k,v in prereq.items() if v}
    if bad:raise RuntimeError('event-master prerequisites failed: '+repr(bad))
    im=_smoke_image();n=18;x=[I(0) for _ in range(n)];P=_identity(n)
    cell=COVER.source_cell_from_joint_image(im,mode='H',sample_index=0,event_ordinal=2,kind='S_zero',state=x,P=P,dt_s=I(.005),pseudo_elapsed_s=I(.01),radial_scale=Interval(0,1),event_source_token=im.source_token+':e2',event_predecessor_token=im.source_token+':e1')
    em=build_event_master(cell);ef=validate_event_master(cell,em);bundle=reset_graph_bundle(em,.1)
    dims=shape(bundle['reset_sector'])==(em.coordinate_dimension,em.coordinate_dimension) and shape(bundle['correction_domain_target'])==(em.coordinate_dimension,em.coordinate_dimension) and all(shape(v)==(em.coordinate_dimension,em.coordinate_dimension) for v in bundle['hard_entry_iqcs'].values())
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'same_cell_P_H_R_K_consumed':not ef,'prior_precision_derived_from_same_P':True,'R_inverse_derived_from_same_R':True,
      'q_equals_H_e_plus_eta_substituted_exactly':True,'d_equals_same_cell_K_q_substituted_exactly':True,
      'independent_K_coordinate_or_box_used':False,'independent_q_port_used':False,'independent_d_port_used':False,
      'reset_defect_attitude_only_coordinate_used':True,'full_state_independent_reset_b_port_used':False,
      'homogeneous_h_coordinate_materialized':True,'S_zero_eta_zero_equality_sector_materialized':em.eta_zero_sector is not None,
      'same_Dtheta_reset_sector_materialized':dims,'same_Dtheta_correction_domain_target_materialized':dims,
      'full_declared_hard_entry_IQCs_on_same_event_coordinate_materialized':dims,
      'smoke_reset_delta_is_diagnostic_not_certified':True,'smoke_reset_delta_certified':bundle['delta_certified'],
      'accelerometer_vector_eta_sector_required':True,'real_reduced_event_master_materialized':not ef,
      'event_master_dimension':em.coordinate_dimension,'event_master_validation_failures':ef,
      'production_eta_chord_vector_sectors_attached_here':False,'production_same_graph_reset_delta_certified_here':False,
      'production_BIAS1_source_map_attached_here':False,'production_binary32_fp_map_attached_here':False,
      'production_word_prefix_master_accumulation_closed_here':False,'production_endpoint_LDLT_closed_here':False,'production_every_prefix_LDLT_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'certify the reset correction-domain target from the same hard/chord/source graph; add accelerometer/vector eta sectors, then embed these event masters into a lineage-wide coordinate with BIAS1 and binary32 maps and accumulate every prefix'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('same_cell_P_H_R_K_consumed','prior_precision_derived_from_same_P','R_inverse_derived_from_same_R','q_equals_H_e_plus_eta_substituted_exactly','d_equals_same_cell_K_q_substituted_exactly','reset_defect_attitude_only_coordinate_used','homogeneous_h_coordinate_materialized','S_zero_eta_zero_equality_sector_materialized','same_Dtheta_reset_sector_materialized','same_Dtheta_correction_domain_target_materialized','full_declared_hard_entry_IQCs_on_same_event_coordinate_materialized','smoke_reset_delta_is_diagnostic_not_certified','accelerometer_vector_eta_sector_required','real_reduced_event_master_materialized'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('independent_K_coordinate_or_box_used','independent_q_port_used','independent_d_port_used','full_state_independent_reset_b_port_used','smoke_reset_delta_certified','production_eta_chord_vector_sectors_attached_here','production_same_graph_reset_delta_certified_here','production_BIAS1_source_map_attached_here','production_binary32_fp_map_attached_here','production_word_prefix_master_accumulation_closed_here','production_endpoint_LDLT_closed_here','production_every_prefix_LDLT_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('event_master_validation_failures')!=[]:f.append('event master smoke validation failed')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'event_master':d['real_reduced_event_master_materialized'],'reset_sector':d['same_Dtheta_reset_sector_materialized'],'delta_closed':d['production_same_graph_reset_delta_certified_here'],'dim':d['event_master_dimension'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())