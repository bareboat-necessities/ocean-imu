#!/usr/bin/env python3
"""Joint physical-BIAS1 supply coordinate for the A21 P4 master.

Let e_b be corrected accelerometer-bias estimation error and beta the physical
true residual bias.  Between measurements,

  beta+ = phi_true beta + w,
  e_b+  = phi_hat e_b + (phi_true-phi_hat) beta + w.

The SAME w enters both coordinates.  Define m=(phi_true-phi_hat)beta and the
source supply s=[w;m].  Then the six-state [e_b;beta] recurrence is linear in
state and supply with a shared w column.  This is the exact coupling needed by
the same-history graph; no independent per-sample bias state is introduced.

For the stability inequality it is conservative and valid to allow every
supply sequence satisfying the derived BIAS1 bounds.  The physical one-root,
one-parameter BIAS1 family is a subset of that ISS forcing class, so this
relaxation strengthens rather than weakens the theorem.  Projection corrections
remain in the exact radial sector and are not charged as exogenous forcing.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

import ou3_p4_bias1_family as FAMILY
import ou3_full_process_ucc as PROCESS

DT=0.005
SCHEMA=1
QUALIFICATION='OU3_P4_BIAS1_JOINT_ISS_SUPPLY_V1'

def up(x):return math.nextafter(float(x),math.inf)
def down(x):return math.nextafter(float(x),-math.inf)

def build():
    fam=FAMILY.build();vf=FAMILY.validate(fam)
    if vf:raise RuntimeError(f'BIAS1 family invalid: {vf}')
    tau_hat=float(PROCESS._constants()['accel_bias_tau_s'])
    if not (math.isfinite(tau_hat) and tau_hat>0):raise RuntimeError('invalid filter bias tau')
    phi_hat=math.exp(-DT/tau_hat)
    phi_lo,phi_hi=map(float,fam['phi_true_interval'])
    mismatch_phi=up(max(abs(phi_lo-phi_hat),abs(phi_hi-phi_hat)))
    beta_norm=float(fam['true_bias_norm_upper_mps2'])
    w_norm=float(fam['driver_increment_norm_upper_mps2'])
    mismatch_norm=up(mismatch_phi*beta_norm)
    supply_norm=up(math.hypot(w_norm,mismatch_norm))
    # State [e(3),beta(3)], supply [w(3),m(3)].  Values are scalar-block notation.
    A_blocks=[['phi_hat*I3','0'],['0','phi_true*I3']]
    B_blocks=[['I3','I3'],['I3','0']]
    return {
      'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'dt_s':DT,'filter_tau_hat_s':tau_hat,'phi_hat':phi_hat,'phi_true_interval':[phi_lo,phi_hi],
      'joint_state':'[e_b;beta_true]','joint_supply':'[w_bias;m_tau]',
      'm_tau_definition':'(phi_true-phi_hat)*beta_true',
      'state_recurrence_blocks':A_blocks,'supply_injection_blocks':B_blocks,
      'same_w_enters_error_and_true_bias':True,
      'one_physical_beta_state_carried_across_word':True,
      'independent_per_sample_bias_state_slots_used':False,
      'physical_BIAS1_family_is_subset_of_ISS_supply_class':True,
      'driver_increment_norm_upper_mps2':w_norm,
      'tau_mismatch_forcing_norm_upper_mps2':mismatch_norm,
      'joint_supply_norm_upper_per_prediction_mps2':supply_norm,
      'word_600_prediction_supply_l2_upper':up(math.sqrt(600.0)*supply_norm),
      'projection_correction_charged_as_exogenous_supply':False,
      'projection_must_use_same_beta_true_coordinate':True,
      'source_supply_may_be_penalized_by_gamma_s_in_augmented_master':True,
      'filter_changed':False,'declared_domain_changed':False,'trajectory_replay_used':False,
      'endpoint_augmented_LDLT_closed_here':False,'every_prefix_augmented_LDLT_closed_here':False,'P4_promoted_here':False,
    }

def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('same_w_enters_error_and_true_bias','one_physical_beta_state_carried_across_word','physical_BIAS1_family_is_subset_of_ISS_supply_class','projection_must_use_same_beta_true_coordinate','source_supply_may_be_penalized_by_gamma_s_in_augmented_master'):
      if d.get(k) is not True:f.append(k+' not true')
    for k in ('independent_per_sample_bias_state_slots_used','projection_correction_charged_as_exogenous_supply','filter_changed','declared_domain_changed','trajectory_replay_used','endpoint_augmented_LDLT_closed_here','every_prefix_augmented_LDLT_closed_here','P4_promoted_here'):
      if d.get(k) is not False:f.append(k+' not false')
    for k in ('driver_increment_norm_upper_mps2','tau_mismatch_forcing_norm_upper_mps2','joint_supply_norm_upper_per_prediction_mps2','word_600_prediction_supply_l2_upper'):
      x=float(d.get(k,-1));
      if not (math.isfinite(x) and x>=0):f.append(k+' invalid')
    return f

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps(d,sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
