#!/usr/bin/env python3
"""Real same-cell reduced signed Joseph/reset master for BRMM P4 prefixes.

For one estimator-owned Joseph ``SourceCoverCell`` this module derives from the
SAME cell

    J=P^-1, H, R^-1, K=P H' (H P H' + R)^-1

and materializes the exact reduced event quadratic on

    z = [e ; eta ; b].

The physical residual and correction are not independent lifted ports:

    q = H e + eta,
    d = K q,

are substituted algebraically into the reduced Joseph/reset identity before the
matrix is returned.  Thus no rowwise/free K can enter this master. ``b`` remains
a lifted reset-defect coordinate because it must be constrained jointly by the
existing exact reset/chord/projection sectors; it is not charged as an
independent norm port.

For S=0 the physical eta is exactly zero.  The returned descriptor marks the
eta-zero equality as mandatory.  Accelerometer/vector events require their
existing exact nonlinear chord/vector sectors.  This module builds the REAL
base event master but does not claim those sectors, BIAS1 supply, binary32 ISS,
or word/prefix accumulation are closed yet.
"""
from __future__ import annotations
from dataclasses import dataclass
import argparse,json
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval,matrix_add,matrix_mul
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan,matrix_symmetric_hull
import ou3_p4_complete_brmm_source_cover_contract as COVER
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_reduced_joseph_reset_event_quadratic as REDUCED
import ou3_p4_joint_brmm_frontend_transition as JOINT
import ou3_brmm_private_mahony_state_step as MAHONY
import ou3_brmm_frontend_state_step as FRONT

SCHEMA=1
QUALIFICATION='OU3_P4_BRMM_REDUCED_SAME_CELL_EVENT_MASTER_V1'

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
    if n!=m:return False
    return all(A[i][j].lo==A[j][i].lo and A[i][j].hi==A[j][i].hi for i in range(n) for j in range(n))

@dataclass(frozen=True)
class EventMaster:
    source_token:str
    estimator_source_token:str
    mode:str
    kind:str
    coordinate_dimension:int
    coordinate_layout:tuple[tuple[str,int,int],...]
    master:Sequence[Sequence[Interval]]
    H:Sequence[Sequence[Interval]]
    K:Sequence[Sequence[Interval]]
    J:Sequence[Sequence[Interval]]
    R_inverse:Sequence[Sequence[Interval]]
    E:Sequence[Sequence[Interval]]
    Q:Sequence[Sequence[Interval]]
    N:Sequence[Sequence[Interval]]
    D:Sequence[Sequence[Interval]]
    B:Sequence[Sequence[Interval]]
    eta_constraint:str
    reset_b_constraint:str

def build_event_master(cell:COVER.SourceCoverCell)->EventMaster:
    failures=COVER.validate_cell(cell,require_estimator_provenance=True)
    if failures:raise ValueError('invalid estimator-owned source cell: '+repr(failures))
    if cell.kind not in ('S_zero','accelerometer','magnetometer'):raise ValueError('reduced event master requires Joseph event')
    kwargs=COVER.joseph_event_kwargs(cell);event=EVENTS.source_joseph_event(**kwargs)
    n=18 if cell.mode=='H' else 21
    J=matrix_inverse_gauss_jordan(cell.P);Rinv=matrix_inverse_gauss_jordan(cell.R)
    nz=2*n+3
    E=selector(n,nz,0);N=selector(3,nz,n);B=selector(n,nz,n+3)
    # q=H e+eta and d=Kq are algebraically eliminated into the base master.
    Q=matrix_add(matrix_mul(event['H'],E),N)
    D=matrix_mul(event['K'],Q)
    M=matrix_symmetric_hull(REDUCED.event_quadratic(J,event['H'],Rinv,E,Q,N,D,B))
    eta_constraint='eta==0 exact equality' if cell.kind=='S_zero' else 'exact nonlinear residual/chord sector on same e/source geometry'
    return EventMaster(
        cell.source_token,cell.estimator_source_token or '',cell.mode,cell.kind,nz,
        (('e',0,n),('eta',n,3),('b_reset',n+3,n)),M,event['H'],event['K'],J,Rinv,E,Q,N,D,B,
        eta_constraint,'b=G^-1*rho exact reset graph/chord sector')

def validate_event_master(cell,em):
    f=[];n=18 if cell.mode=='H' else 21
    if em.source_token!=cell.source_token:f.append('event token detached')
    if em.estimator_source_token!=cell.estimator_source_token:f.append('estimator token detached')
    if shape(em.master)!=(2*n+3,2*n+3):f.append('master dimension mismatch')
    if not symmetric(em.master):f.append('master lost symmetric hull')
    if shape(em.H)!=(3,n) or shape(em.K)!=(n,3):f.append('same-cell H/K dimensions invalid')
    if shape(em.J)!=(n,n) or shape(em.R_inverse)!=(3,3):f.append('precision dimensions invalid')
    if shape(em.Q)!=(3,2*n+3) or shape(em.D)!=(n,2*n+3):f.append('derived q/d maps invalid')
    if em.kind=='S_zero' and em.eta_constraint!='eta==0 exact equality':f.append('S=0 eta equality lost')
    return f

def _smoke_image():
    st=JOINT._smoke_state();sample=FRONT.Sample(MAHONY.Vec3(MAHONY.I(.01),MAHONY.I(-.02),MAHONY.I(.005)),MAHONY.Vec3(MAHONY.I(.2),MAHONY.I(-.1),MAHONY.I(-9.75)))
    images=JOINT.advance(st,sample,gravity_ms2=MAHONY.I(9.80665),two_kp=MAHONY.I(.2),two_ki=MAHONY.I(.02),child_prefix='master')
    if not images:raise RuntimeError('joint image smoke empty')
    return images[0]
def _identity(n):return [[I(1 if i==j else 0) for j in range(n)] for i in range(n)]
def build():
    r=REDUCED.build();rf=REDUCED.validate(r)
    if rf:raise RuntimeError('reduced identity prerequisite failed: '+repr(rf))
    im=_smoke_image();n=18;x=[I(0) for _ in range(n)];P=_identity(n)
    cell=COVER.source_cell_from_joint_image(im,mode='H',sample_index=0,event_ordinal=2,kind='S_zero',state=x,P=P,dt_s=I(.005),pseudo_elapsed_s=I(.01),radial_scale=Interval(0,1),event_source_token=im.source_token+':e2',event_predecessor_token=im.source_token+':e1')
    em=build_event_master(cell);ef=validate_event_master(cell,em)
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'same_cell_P_H_R_K_consumed':not ef,'prior_precision_derived_from_same_P':True,'R_inverse_derived_from_same_R':True,
      'q_equals_H_e_plus_eta_substituted_exactly':True,'d_equals_same_cell_K_q_substituted_exactly':True,
      'independent_K_coordinate_or_box_used':False,'independent_q_port_used':False,'independent_d_port_used':False,
      'reset_b_kept_as_joint_lifted_coordinate':True,'reset_b_independent_norm_port_used':False,
      'S_zero_eta_zero_equality_required':True,'accelerometer_vector_eta_sector_required':True,
      'real_reduced_event_master_materialized':not ef,'event_master_dimension':em.coordinate_dimension,
      'event_master_validation_failures':ef,'production_eta_reset_projection_sectors_attached_here':False,
      'production_BIAS1_source_map_attached_here':False,'production_binary32_fp_map_attached_here':False,
      'production_word_prefix_master_accumulation_closed_here':False,'production_endpoint_LDLT_closed_here':False,'production_every_prefix_LDLT_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'embed each real event master in one lineage-wide augmented coordinate, attach exact eta/chord and reset-b/projection sectors plus BIAS1 and binary32 maps, accumulate signed matrices at every literal prefix, then run outward LDLT'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('same_cell_P_H_R_K_consumed','prior_precision_derived_from_same_P','R_inverse_derived_from_same_R','q_equals_H_e_plus_eta_substituted_exactly','d_equals_same_cell_K_q_substituted_exactly','reset_b_kept_as_joint_lifted_coordinate','S_zero_eta_zero_equality_required','accelerometer_vector_eta_sector_required','real_reduced_event_master_materialized'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('independent_K_coordinate_or_box_used','independent_q_port_used','independent_d_port_used','reset_b_independent_norm_port_used','production_eta_reset_projection_sectors_attached_here','production_BIAS1_source_map_attached_here','production_binary32_fp_map_attached_here','production_word_prefix_master_accumulation_closed_here','production_endpoint_LDLT_closed_here','production_every_prefix_LDLT_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('event_master_validation_failures')!=[]:f.append('event master smoke validation failed')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'event_master':d['real_reduced_event_master_materialized'],'dim':d['event_master_dimension'],'prefix_closed':d['production_every_prefix_LDLT_closed_here'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())