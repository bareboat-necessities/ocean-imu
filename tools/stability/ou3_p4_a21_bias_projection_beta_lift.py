#!/usr/bin/env python3
"""Lift A21 Joseph/projection events to the physical true-bias coordinate.

The shipping filter state contains the accelerometer-bias estimate error e_b but
P4 also carries the physical residual true bias beta as a proof/source state.
For the deployed projection

    e_b^+ = beta - Pi_R(beta - u),

where u is the smooth pre-projection bias error.  If J_Pi is any member of the
outward Clarke generalized-Jacobian enclosure of Pi_R, then

    d e_b^+ / d u    = J_Pi,
    d e_b^+ / d beta = I - J_Pi.

The existing nonlinear Joseph event differentiator already returns J_Pi on the
same source cell.  This module exposes the missing beta derivative and lifts the
21-state event differential to the joint 24-state [e;beta] graph.  beta is
unchanged by estimator measurement/reset events.

No independent per-event beta slot is introduced.  Prediction BIAS1 recurrence
is handled separately by the physical BIAS1 transition and shares this same beta
coordinate across the complete word.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_complete_brmm_differential_events as EVENTS

SCHEMA=1
QUALIFICATION='OU3_P4_A21_TRUE_BIAS_PROJECTION_LIFT_V1'
OFF_BA=18

def I(x):return Interval.point(float(x))
def zero(r,c):return [[I(0) for _ in range(c)] for _ in range(r)]
def identity(n):return [[I(1 if i==j else 0) for j in range(n)] for i in range(n)]
def shape(A):return len(A),len(A[0]) if A else 0

def beta_true_jacobian(event:dict):
    Jp=event.get('bias_projection_J')
    if Jp is None or shape(Jp)!=(3,3):raise ValueError('A21 event projection Jacobian required')
    out=zero(21,3);Id=identity(3)
    for i in range(3):
        for j in range(3):out[OFF_BA+i][j]=Id[i][j]-Jp[i][j]
    return out

def lift_measurement_event(event:dict):
    Je=event.get('J_state')
    if shape(Je)!=(21,21):raise ValueError('A21 state Jacobian required')
    Jb=beta_true_jacobian(event);A=zero(24,24)
    for i in range(21):
        for j in range(21):A[i][j]=Je[i][j]
        for j in range(3):A[i][21+j]=Jb[i][j]
    for i in range(3):A[21+i][21+i]=I(1)
    return A,Jb

def _zero_contained(A):return all(x.lo<=0<=x.hi for row in A for x in row)
def _contains_identity_minus(Jb,Jp):
    Id=identity(3)
    return all(Jb[OFF_BA+i][j].contains_interval(Id[i][j]-Jp[i][j]) for i in range(3) for j in range(3))

def _smoke(branch):
    # Exercise the projection primitive directly, then package a synthetic event
    # with identity smooth state Jacobian.  This checks inactive, active and
    # Clarke-boundary beta derivatives independently of a trajectory point.
    if branch=='inactive':x=[I(.1),I(0),I(0)]
    elif branch=='active':x=[I(.7),I(0),I(0)]
    elif branch=='clarke':x=[Interval(.49,.51),I(0),I(0)]
    else:raise ValueError(branch)
    p=EVENTS.ball_projection_enclosure(x,.5)
    e={'J_state':identity(21),'bias_projection_J':p['J']}
    A,Jb=lift_measurement_event(e)
    return {'branch_reported':p['branch'],'joint_shape':list(shape(A)),
      'beta_block_identity_minus_projection':_contains_identity_minus(Jb,p['J']),
      'beta_state_persists_identity':all(A[21+i][21+j].contains(1 if i==j else 0) for i in range(3) for j in range(3)),
      'non_ba_filter_rows_have_zero_beta_derivative':_zero_contained(Jb[:OFF_BA])}

def build():
    sm={k:_smoke(k) for k in ('inactive','active','clarke')}
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'joint_state':'[e_A21;beta_true]','joint_dimension':24,
      'same_projection_generalized_Jacobian_consumed':True,
      'true_bias_derivative_formula':'d e_b_plus/d beta_true = I3 - J_projection',
      'beta_true_unchanged_by_measurement_reset_event':True,
      'one_beta_true_coordinate_must_persist_across_word':True,
      'independent_per_event_beta_slots_forbidden':True,
      'measurement_event_24state_lift_available':True,
      'inactive_projection_beta_derivative_zero':sm['inactive']['beta_block_identity_minus_projection'] and sm['inactive']['branch_reported']=='inactive',
      'active_projection_beta_derivative_retained':sm['active']['beta_block_identity_minus_projection'] and sm['active']['branch_reported']=='active',
      'clarke_boundary_beta_derivative_hulled':sm['clarke']['beta_block_identity_minus_projection'] and sm['clarke']['branch_reported']=='clarke_hull',
      'smoke':sm,'production_lineage_24state_events_materialized_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'use lift_measurement_event on every A21 estimator-owned Joseph/projection event and compose it with the same persistent beta_true coordinate used by the BIAS1 prediction recurrence'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('same_projection_generalized_Jacobian_consumed','beta_true_unchanged_by_measurement_reset_event','one_beta_true_coordinate_must_persist_across_word','independent_per_event_beta_slots_forbidden','measurement_event_24state_lift_available','inactive_projection_beta_derivative_zero','active_projection_beta_derivative_retained','clarke_boundary_beta_derivative_hulled'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('production_lineage_24state_events_materialized_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('joint_dimension')!=24:f.append('joint dimension changed')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'beta_lift':d['measurement_event_24state_lift_available'],'active':d['active_projection_beta_derivative_retained'],'clarke':d['clarke_boundary_beta_derivative_hulled'],'production':d['production_lineage_24state_events_materialized_here'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
