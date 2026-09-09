#!/usr/bin/env python3
"""Source-uniform joint ISS supply for every mandatory bias family.

The BIAS1 joint supply already fixes the shape of the coupling.  With e_b the
corrected accelerometer-bias estimation error and b the shipping-centered true
bias, each admitted family satisfies, between measurements,

    b+   = phi_true*b + w,
    e_b+ = phi_hat*e_b + (phi_true-phi_hat)*b + w,

so with m := (phi_true-phi_hat)*b the joint state [e_b;b] is linear in the
supply s=[w;m] with a SHARED w column.  This module materializes that same
recurrence separately for BIAS0, BIAS1 and BIAS2, each from its own
authoritative family module:

* BIAS0 contributes a composite multi-channel driver whose increment is
  dominated by its declared pathwise Gauss-Markov cap;
* BIAS1 contributes the one-root one-parameter driver (delegated verbatim to
  `ou3_p4_bias1_joint_iss_supply`, never restated here);
* BIAS2 contributes a tiny bounded-variation driver but the largest admitted
  tau mismatch, because it admits phi_true=1.

No family is inferred from another, and the union is not collapsed into a
generic bias box.  Each family's interval is pushed through the deployed
24-state event lift to check that the lift is genuinely family-parametric and
that the shared-w column survives, so the supply is available to the augmented
master.  Materializing the supply is a prerequisite: endpoint and every-prefix
augmented LDLT remain open for all three families, and BIAS2 additionally
needs its separation sector.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

from ou3_interval import Interval
import ou3_full_process_ucc as PROCESS
import ou3_p4_bias0_family as BIAS0
import ou3_p4_bias1_family as BIAS1
import ou3_p4_bias2_family as BIAS2
import ou3_p4_bias1_joint_iss_supply as BIAS1_ISS
import ou3_p4_a21_bias1_24state_event_lift as LIFT
import ou3_p4_hard_entry_set as ENTRY

DT=0.005
WORD=600
SCHEMA=1
QUALIFICATION='OU3_P4_ALL_BIAS_FAMILY_JOINT_ISS_SUPPLY_V1'
REQUIRED_BIAS_FAMILIES=('BIAS0','BIAS1','BIAS2')
OFF_BA=18


def up(x):return math.nextafter(float(x),math.inf)


def _family_modules():
    return {'BIAS0':BIAS0,'BIAS1':BIAS1,'BIAS2':BIAS2}


def _supply(phi_hat,phi_lo,phi_hi,beta_norm,w_norm):
    """Outward joint supply bound for one family's admitted factor interval."""
    mismatch_phi=up(max(abs(phi_lo-phi_hat),abs(phi_hi-phi_hat)))
    mismatch_norm=up(mismatch_phi*beta_norm)
    return mismatch_phi,mismatch_norm,up(math.hypot(w_norm,mismatch_norm))


def _lift_is_family_parametric(phi_hat,phi_lo,phi_hi):
    """Push this family's factor interval through the deployed 24-state lift."""
    F=LIFT.identity(21)
    for i in range(3):F[OFF_BA+i][OFF_BA+i]=LIFT.I(phi_hat)
    A,B=LIFT.prediction_lift(F,Interval(phi_lo,phi_hi))
    carries=all(A[21+i][21+i].lo==phi_lo and A[21+i][21+i].hi==phi_hi for i in range(3))
    return {'prediction_lift_accepts_family_interval':LIFT.shape(A)==(24,24),
            'supply_injection_available':LIFT.shape(B)==(24,6),
            'same_w_column_shared_by_error_and_truth':LIFT.source_map_same_w(B),
            'physical_factor_carried_in_truth_block':carries}


def build():
    modules=_family_modules()
    families={}
    for name,mod in modules.items():
        d=mod.build();f=mod.validate(d)
        if f:raise RuntimeError(f'{name} family prerequisite failed: {f!r}')
        families[name]=d
    tau_hat=float(PROCESS._constants()['accel_bias_tau_s'])
    if not (math.isfinite(tau_hat) and tau_hat>0):raise RuntimeError('invalid filter bias tau')
    phi_hat=math.exp(-DT/tau_hat)

    supply={}
    for name,d in families.items():
        phi_lo,phi_hi=map(float,d['phi_true_interval'])
        beta_norm=float(d['true_bias_norm_upper_mps2'])
        w_norm=float(d['driver_increment_norm_upper_mps2'])
        mismatch_phi,mismatch_norm,joint=_supply(phi_hat,phi_lo,phi_hi,beta_norm,w_norm)
        supply[name]={
          'family_module_qualification':d['qualification'],
          'phi_true_interval':[phi_lo,phi_hi],
          'tau_mismatch_factor_upper':mismatch_phi,
          'true_bias_norm_upper_mps2':beta_norm,
          'driver_increment_norm_upper_mps2':w_norm,
          'tau_mismatch_forcing_norm_upper_mps2':mismatch_norm,
          'joint_supply_norm_upper_per_prediction_mps2':joint,
          f'word_{WORD}_prediction_supply_l2_upper':up(math.sqrt(float(WORD))*joint),
          'event_lift':_lift_is_family_parametric(phi_hat,phi_lo,phi_hi),
        }

    # BIAS1 is delegated, never restated: agree exactly with its own module.
    b1=BIAS1_ISS.build();b1f=BIAS1_ISS.validate(b1)
    if b1f:raise RuntimeError(f'BIAS1 joint ISS prerequisite failed: {b1f!r}')
    delegated=(supply['BIAS1']['driver_increment_norm_upper_mps2']==float(b1['driver_increment_norm_upper_mps2'])
               and supply['BIAS1']['tau_mismatch_forcing_norm_upper_mps2']==float(b1['tau_mismatch_forcing_norm_upper_mps2'])
               and supply['BIAS1']['joint_supply_norm_upper_per_prediction_mps2']==float(b1['joint_supply_norm_upper_per_prediction_mps2']))
    if not delegated:raise RuntimeError('BIAS1 supply diverges from its authoritative module')

    lifts_ok=all(all(supply[n]['event_lift'].values()) for n in REQUIRED_BIAS_FAMILIES)

    # Relaxation-class containment on the ISS-relevant pair (admitted factor
    # interval, driver-increment bound).  This is what decides whether one
    # family's proof would carry another's.  BIAS0's wider Gauss-Markov tau
    # range and larger pathwise driver do contain BIAS1, so BIAS0 is the harder
    # family on the supply axis; the direction that matters for the contract is
    # the opposite one, and BIAS1 contains neither of the other two.
    def class_contains(a,b):
        pa,pb=supply[a]['phi_true_interval'],supply[b]['phi_true_interval']
        return (pa[0]<=pb[0] and pb[1]<=pa[1]
                and supply[b]['driver_increment_norm_upper_mps2']<=supply[a]['driver_increment_norm_upper_mps2'])
    containment={a:{b:class_contains(a,b) for b in REQUIRED_BIAS_FAMILIES if b!=a} for a in REQUIRED_BIAS_FAMILIES}
    bias1_implies_nothing=not any(containment['BIAS1'].values())
    no_family_covers_all=not any(all(containment[a].values()) for a in REQUIRED_BIAS_FAMILIES)
    entry_radius=float(ENTRY.build()['coordinate_radii']['accelerometer_bias_error_norm_mps2'])
    truth_within_entry_envelope=all(float(supply[n]['true_bias_norm_upper_mps2'])<=entry_radius for n in REQUIRED_BIAS_FAMILIES)

    return {
      'schema':SCHEMA,'qualification':QUALIFICATION,
      'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'required_bias_families':list(REQUIRED_BIAS_FAMILIES),
      'dt_s':DT,'word_prediction_count':WORD,
      'filter_tau_hat_s':tau_hat,'phi_hat':phi_hat,
      'joint_state':'[e_b;b_true]','joint_supply':'[w_bias;m_tau]',
      'm_tau_definition':'(phi_true-phi_hat)*b_true',
      'family_supply':supply,
      'same_w_enters_error_and_true_bias':all(supply[n]['event_lift']['same_w_column_shared_by_error_and_truth'] for n in REQUIRED_BIAS_FAMILIES),
      'one_physical_bias_state_carried_across_word':True,
      'independent_per_sample_bias_state_slots_used':False,
      'each_family_pushed_through_deployed_24state_event_lift':lifts_ok,
      'relaxation_class_containment':containment,
      'BIAS1_class_contains_no_other_family':bias1_implies_nothing,
      'no_single_family_covers_the_other_two':no_family_covers_all,
      'hard_entry_accelerometer_bias_error_radius_mps2':entry_radius,
      'every_family_true_bias_within_hard_entry_envelope':truth_within_entry_envelope,
      'BIAS1_supply_delegated_to_authoritative_module':delegated,
      'BIAS0_inferred_from_BIAS1':False,'BIAS2_inferred_from_BIAS1':False,
      'generic_bias_box_replaces_three_families':False,
      'largest_joint_supply_family':max(REQUIRED_BIAS_FAMILIES,key=lambda n:supply[n]['joint_supply_norm_upper_per_prediction_mps2']),
      'largest_tau_mismatch_family':max(REQUIRED_BIAS_FAMILIES,key=lambda n:supply[n]['tau_mismatch_factor_upper']),
      'projection_correction_charged_as_exogenous_supply':False,
      'projection_must_use_same_b_true_coordinate':True,
      'BIAS2_separation_sector_still_required':bool(families['BIAS2']['separation_sector_required_for_ISS_closure']),
      'BIAS2_uniform_separation_closed':bool(families['BIAS2']['source_uniform_separation_closed']),
      'filter_changed':False,'declared_domain_changed':False,'trajectory_replay_used':False,
      'endpoint_augmented_LDLT_closed_here':False,
      'every_prefix_augmented_LDLT_closed_here':False,
      'P4_promoted_here':False,
    }


def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('required_bias_families')!=list(REQUIRED_BIAS_FAMILIES):f.append('required bias family set changed')
    supply=d.get('family_supply',{})
    if set(supply)!=set(REQUIRED_BIAS_FAMILIES):f.append('family supply table must be exactly BIAS0/BIAS1/BIAS2')
    for k in ('same_w_enters_error_and_true_bias','one_physical_bias_state_carried_across_word','each_family_pushed_through_deployed_24state_event_lift','BIAS1_class_contains_no_other_family','no_single_family_covers_the_other_two','every_family_true_bias_within_hard_entry_envelope','BIAS1_supply_delegated_to_authoritative_module','projection_must_use_same_b_true_coordinate','BIAS2_separation_sector_still_required'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('independent_per_sample_bias_state_slots_used','BIAS0_inferred_from_BIAS1','BIAS2_inferred_from_BIAS1','generic_bias_box_replaces_three_families','projection_correction_charged_as_exogenous_supply','BIAS2_uniform_separation_closed','filter_changed','declared_domain_changed','trajectory_replay_used','endpoint_augmented_LDLT_closed_here','every_prefix_augmented_LDLT_closed_here','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    containment=d.get('relaxation_class_containment',{})
    if set(containment)!=set(REQUIRED_BIAS_FAMILIES):f.append('relaxation containment table incomplete')
    elif any(containment[a].get(b) for a in ('BIAS1',) for b in containment.get('BIAS1',{})):f.append('BIAS1 falsely covers another family')
    for name in REQUIRED_BIAS_FAMILIES:
        s=supply.get(name,{})
        for k in ('driver_increment_norm_upper_mps2','tau_mismatch_forcing_norm_upper_mps2','joint_supply_norm_upper_per_prediction_mps2',f'word_{WORD}_prediction_supply_l2_upper'):
            x=float(s.get(k,-1.0))
            if not (math.isfinite(x) and x>0):f.append(f'{name}.{k} invalid')
        lift=s.get('event_lift',{})
        if not all(bool(v) for v in lift.values()) or len(lift)!=4:f.append(name+' event lift incomplete')
    return f


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({n:{'w':d['family_supply'][n]['driver_increment_norm_upper_mps2'],
                         'm':d['family_supply'][n]['tau_mismatch_forcing_norm_upper_mps2'],
                         's':d['family_supply'][n]['joint_supply_norm_upper_per_prediction_mps2']}
                      for n in REQUIRED_BIAS_FAMILIES}|{'failures':f},sort_keys=True))
    return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
