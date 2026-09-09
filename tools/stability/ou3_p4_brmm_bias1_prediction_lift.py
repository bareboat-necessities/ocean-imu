#!/usr/bin/env python3
"""Outward physical BIAS1 prediction lift for the A21 P4 lineage.

Carry the six proof/source coordinates

    x_b = [e_b(3); beta_true(3)]

and the six physical supply coordinates

    s_b = [w(3); m_tau(3)],
    m_tau = (phi_true-phi_hat) beta_true.

The exact recurrence is

    e_b+   = phi_hat e_b + w + m_tau,
    beta+  = phi_true beta + w.

The SAME w columns therefore enter estimation error and physical true bias.
phi_true is the admitted BIAS1 interval; phi_hat is the deployed filter decay.
This module emits literal outward A/B matrices and a homogeneous joint-supply
ball IQC.  It does not create independent per-sample beta states and does not
fit a source trajectory.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

from ou3_interval import Interval
import ou3_p4_bias1_joint_iss_supply as BIASISS
import ou3_p4_affine_hard_tube_iqc as HARD

SCHEMA=1
QUALIFICATION='OU3_P4_BRMM_BIAS1_PHYSICAL_PREDICTION_LIFT_V1'

def I(x):return Interval.point(float(x))
def zeros(r,c):return [[I(0) for _ in range(c)] for _ in range(r)]
def build_matrices(contract):
    ph=I(float(contract['phi_hat']));pt=Interval(*map(float,contract['phi_true_interval']))
    A=zeros(6,6);B=zeros(6,6)
    for i in range(3):
        A[i][i]=ph
        A[3+i][3+i]=pt
        # supply ordering [w0,w1,w2,m0,m1,m2]
        B[i][i]=I(1);B[i][3+i]=I(1)
        B[3+i][i]=I(1)
    return A,B

def shared_w_columns_exact(B):
    return all(B[i][i].lo==1 and B[i][i].hi==1 and B[3+i][i].lo==1 and B[3+i][i].hi==1 for i in range(3))
def m_only_enters_error(B):
    return all(B[i][3+i].lo==1 and B[i][3+i].hi==1 and all(B[3+j][3+i].lo<=0<=B[3+j][3+i].hi for j in range(3)) for i in range(3))
def build():
    c=BIASISS.build();cf=BIASISS.validate(c)
    if cf:raise RuntimeError('BIAS1 ISS prerequisite failed: '+repr(cf))
    A,B=build_matrices(c)
    # Homogeneous coordinate z=[h; e_b(3); beta(3); w(3); m(3)].
    n=13;h=0;supply_idx=tuple(range(7,13));bound=float(c['joint_supply_norm_upper_per_prediction_mps2'])
    supply_iqc=HARD.ball_iqc(n,h,supply_idx,bound)
    finite=all(math.isfinite(x.lo) and math.isfinite(x.hi) for M in (A,B,supply_iqc) for row in M for x in row)
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'joint_state':'[e_b;beta_true]','joint_supply':'[w_bias;m_tau]','m_tau_definition':'(phi_true-phi_hat)*beta_true',
      'phi_hat':c['phi_hat'],'phi_true_interval':c['phi_true_interval'],
      'outward_A_matrix_materialized':True,'outward_B_matrix_materialized':True,
      'same_w_columns_enter_error_and_true_bias':shared_w_columns_exact(B),
      'm_tau_columns_enter_error_not_true_bias':m_only_enters_error(B),
      'one_physical_beta_state_carried_across_prediction':True,'independent_per_sample_beta_slots_used':False,
      'physical_BIAS1_family_is_subset_of_supply_class':c['physical_BIAS1_family_is_subset_of_ISS_supply_class'],
      'joint_supply_bound_mps2':bound,'homogeneous_joint_supply_ball_IQC_materialized':True,
      'all_interval_matrices_finite':finite,'trajectory_replay_used':False,
      'production_full_prediction_transport_embedded_here':False,'production_every_prefix_storage_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':'embed this 6-state physical BIAS1 recurrence into each A21 prediction transport while preserving the other 18 shipping error coordinates; then connect the carried beta state to the active projection prefix and add binary32 forcing'}
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('outward_A_matrix_materialized','outward_B_matrix_materialized','same_w_columns_enter_error_and_true_bias','m_tau_columns_enter_error_not_true_bias','one_physical_beta_state_carried_across_prediction','physical_BIAS1_family_is_subset_of_supply_class','homogeneous_joint_supply_ball_IQC_materialized','all_interval_matrices_finite'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('independent_per_sample_beta_slots_used','trajectory_replay_used','production_full_prediction_transport_embedded_here','production_every_prefix_storage_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if not(math.isfinite(float(d.get('joint_supply_bound_mps2',math.nan))) and float(d['joint_supply_bound_mps2'])>=0):f.append('joint supply bound invalid')
    return f
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'bias_lift':d['same_w_columns_enter_error_and_true_bias'],'supply_bound':d['joint_supply_bound_mps2'],'P4':d['P4_PASS'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
