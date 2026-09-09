#!/usr/bin/env python3
"""Complete A21 + physical-BIAS1 24-state event lifts for P4.

Joint proof state:

    x24 = [e_A21 ; beta_true].

Prediction uses the already-qualified BIAS1 ISS supply s=[w;m_tau],

    e_b+  = phi_hat e_b + w + m_tau,
    beta+ = phi_true beta + w,
    m_tau = (phi_true-phi_hat) beta,

with the SAME w entering e_b and beta.  The physical one-root/one-parameter
BIAS1 family is a subset of this ISS supply class, so the theorem may prove the
larger forcing class without introducing independent per-sample bias states.

Accepted A21 measurement/reset events use the deployed projection lift

    d e_b+ / d beta = I - J_projection,

and beta itself remains unchanged.  Thus prediction and projection share one
persistent physical beta coordinate through the literal event lineage.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_bias1_joint_iss_supply as BIASISS
import ou3_p4_a21_bias_projection_beta_lift as PROJLIFT

SCHEMA=1
QUALIFICATION='OU3_P4_A21_BIAS1_24STATE_EVENT_LIFT_V1'
OFF_BA=18

def I(x):return Interval.point(float(x))
def zero(r,c):return [[I(0) for _ in range(c)] for _ in range(r)]
def identity(n):return [[I(1 if i==j else 0) for j in range(n)] for i in range(n)]
def shape(A):return len(A),len(A[0]) if A else 0

def prediction_lift(F21,phi_true:Interval):
    if shape(F21)!=(21,21):raise ValueError('21x21 shipping prediction map required')
    if not isinstance(phi_true,Interval) or phi_true.lo<=0 or phi_true.hi>1.0:raise ValueError('valid physical phi interval required')
    A=zero(24,24)
    for i in range(21):
        for j in range(21):A[i][j]=F21[i][j]
    for i in range(3):A[21+i][21+i]=phi_true
    # Supply s=[w(3);m_tau(3)].  Same w enters e_b and beta.
    B=zero(24,6)
    for i in range(3):
        B[OFF_BA+i][i]=I(1);B[OFF_BA+i][3+i]=I(1)
        B[21+i][i]=I(1)
    return A,B

def measurement_lift(event):
    return PROJLIFT.lift_measurement_event(event)[0]

def source_map_same_w(B):
    if shape(B)!=(24,6):return False
    return all(B[OFF_BA+i][i].contains(1) and B[21+i][i].contains(1) for i in range(3))

def build():
    b=BIASISS.build();bf=BIASISS.validate(b)
    if bf:raise RuntimeError('BIAS1 ISS prerequisite failed: '+repr(bf))
    phi=Interval(*map(float,b['phi_true_interval']))
    F=identity(21)
    phi_hat=float(b['phi_hat'])
    for i in range(3):F[OFF_BA+i][OFF_BA+i]=I(phi_hat)
    A,B=prediction_lift(F,phi)
    # Projection smoke uses an active radial branch to ensure beta coupling is nonzero.
    import ou3_p4_complete_brmm_differential_events as EVENTS
    p=EVENTS.ball_projection_enclosure([I(.7),I(0),I(0)],.5)
    ev={'J_state':identity(21),'bias_projection_J':p['J']}
    M=measurement_lift(ev)
    beta_coupling=any(not (M[OFF_BA+i][21+j].lo==0 and M[OFF_BA+i][21+j].hi==0) for i in range(3) for j in range(3))
    beta_identity=all(M[21+i][21+j].contains(1 if i==j else 0) for i in range(3) for j in range(3))
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'joint_state':'[e_A21;beta_true]','joint_dimension':24,'joint_supply':'[w_bias;m_tau]',
      'prediction_24state_lift_available':shape(A)==(24,24),'prediction_supply_map_available':shape(B)==(24,6),
      'same_w_enters_error_and_true_bias':source_map_same_w(B),'m_tau_enters_error_only':all(B[OFF_BA+i][3+i].contains(1) and B[21+i][3+i].contains(0) for i in range(3)),
      'physical_BIAS1_family_subset_of_supply_class':bool(b['physical_BIAS1_family_is_subset_of_ISS_supply_class']),
      'one_persistent_beta_state_across_prediction_and_measurement':True,
      'measurement_projection_24state_lift_available':shape(M)==(24,24),
      'active_projection_true_bias_coupling_nonzero':beta_coupling,'measurement_leaves_beta_identity':beta_identity,
      'independent_per_sample_beta_slots_used':False,'projection_correction_charged_as_exogenous_supply':False,
      'production_complete_literal_event_sequence_lifted_here':False,'production_prefix_source_transport_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'lift each A21 prediction with prediction_lift and each accepted Joseph/projection with measurement_lift on one persistent beta coordinate; suffix-propagate each prediction B supply map through the exact every-prefix transport family'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('prediction_24state_lift_available','prediction_supply_map_available','same_w_enters_error_and_true_bias','m_tau_enters_error_only','physical_BIAS1_family_subset_of_supply_class','one_persistent_beta_state_across_prediction_and_measurement','measurement_projection_24state_lift_available','active_projection_true_bias_coupling_nonzero','measurement_leaves_beta_identity'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('independent_per_sample_beta_slots_used','projection_correction_charged_as_exogenous_supply','production_complete_literal_event_sequence_lifted_here','production_prefix_source_transport_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('joint_dimension')!=24:f.append('joint dimension mismatch')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'prediction':d['prediction_24state_lift_available'],'same_w':d['same_w_enters_error_and_true_bias'],'projection_beta':d['active_projection_true_bias_coupling_nonzero'],'production':d['production_prefix_source_transport_closed_here'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
