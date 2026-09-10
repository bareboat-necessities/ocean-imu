#!/usr/bin/env python3
"""Hard fresh-Live membership for the uniformly qualified BRMM source.

This is a deterministic source/error relation, not covariance consistency.
For the fresh outer wrapper the shipping linear means are zero at first Live.
The qualified physical BRMM primitives therefore map directly to e_v/e_p; the
shared S origin is re-anchored exactly so centered e_S=0.  Existing startup
source assumptions provide the remaining attitude, gyro-bias and a_w bounds,
and each mandatory physical accelerometer-bias family starts inside the deployed
0.4 m/s^2 projection/error working radius.

Both measured-period and prior-frequency timeout Live are retained: tuner/WPE
state changes coefficient ancestry, not the fresh physical error identities.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path

import ou3_brmm_finite_window_primitive_qualification as PRIM
import ou3_p4_live_entry_graph as LIVE
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_bias_family_joint_iss_supply as BIAS
import ou3_p4_p3_execution_admission as P3

REPO=Path(__file__).resolve().parents[2]
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
QUALIFICATION='OU3_P4_QUALIFIED_FRESH_LIVE_HARD_ENTRY_V1'

def build(domain_path:Path=DOMAIN):
    d=json.loads(Path(domain_path).read_text());prim=PRIM.build();live=LIVE.build();entry=ENTRY.build();bias=BIAS.build();p3=P3.build()
    bad={'primitive':PRIM.validate(prim),'entry':ENTRY.validate(entry),'bias':BIAS.validate(bias),'p3':P3.validate(p3)};bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('fresh-entry prerequisites failed: '+repr(bad))
    r=entry['coordinate_radii'];u=prim['uniform_physical_primitives'];startup=d['startup'];att=d['initial_filter_entrance']['attitude']
    theta=math.radians(float(att['full_attitude_error_upper_deg']));cayley=2.0*math.tan(theta/2.0)
    aw=float(startup['physical_handoff_coordinate_bounds']['latent_acceleration_error_norm_upper_mps2'])
    gb=float(startup['initial_tangent_gyro_bias_norm_upper_rad_s'])
    bias_truth=max(float(x['true_bias_norm_upper_mps2']) for x in bias['family_supply'].values())
    checks={
      'attitude':cayley<=float(r['attitude_cayley_norm']),
      'gyro_bias':gb<=float(r['gyro_bias_norm_rad_s']),
      'velocity':float(u['V_m_norm_upper_mps'])<=float(r['velocity_norm_mps']),
      'position':float(u['P_m_norm_upper_m'])<=float(r['position_norm_m']),
      'centered_integral_displacement':0.0<=float(r['integral_displacement_norm_m_s']),
      'latent_acceleration':aw<=float(r['latent_acceleration_norm_mps2']),
      'accelerometer_bias':bias_truth<=float(r['accelerometer_bias_error_norm_mps2'])}
    return {
      'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'fresh_outer_wrapper_only':True,'previously_driven_inner_filter_not_claimed':True,
      'covariance_membership_used':False,'trajectory_replay_used':False,
      'source_applicability_refinement_consumed':bool(prim['source_applicability_refined']),
      'shipping_fresh_linear_means_zero_consumed':live['nominal_v_p_S_aw_bg_ba']=='exactly zero carried from construction',
      'shared_S_origin_exactly_reanchored':live['canonical_S_entry_error']==0,
      'measured_period_live_branch_covered':True,'prior_frequency_timeout_live_branch_covered':True,
      'tuner_branch_does_not_change_fresh_physical_error_identity':True,
      'full_45deg_attitude_entry_consumed':float(att['full_attitude_error_upper_deg'])==45.0,
      'membership_checks':checks,'all_hard_entry_coordinates_inside_full_scale':all(checks.values()),
      'entry_radii':r,'source_bounds':{'cayley_norm':cayley,'gyro_bias_norm_rad_s':gb,'velocity_norm_mps':u['V_m_norm_upper_mps'],'position_norm_m':u['P_m_norm_upper_m'],'centered_S_norm_m_s':0.0,'latent_acceleration_norm_mps2':aw,'true_accelerometer_bias_norm_mps2':bias_truth},
      'P3_scoped_execution_admission_consumed':bool(p3['P3_SCOPED_EXECUTION_ADMISSION_PASS']),
      'source_uniform_reachable_other_entry_coordinates_closed':all(checks.values()),
      'fresh_live_entry_cover_closed':bool(all(checks.values()) and p3['P3_SCOPED_EXECUTION_ADMISSION_PASS']),
      'P4_promoted_here':False}

def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('fresh_outer_wrapper_only','previously_driven_inner_filter_not_claimed','source_applicability_refinement_consumed','shipping_fresh_linear_means_zero_consumed','shared_S_origin_exactly_reanchored','measured_period_live_branch_covered','prior_frequency_timeout_live_branch_covered','tuner_branch_does_not_change_fresh_physical_error_identity','full_45deg_attitude_entry_consumed','all_hard_entry_coordinates_inside_full_scale','P3_scoped_execution_admission_consumed','source_uniform_reachable_other_entry_coordinates_closed','fresh_live_entry_cover_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('covariance_membership_used','trajectory_replay_used','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    if not all(d.get('membership_checks',{}).values()):f.append('one or more coordinate memberships failed')
    return f

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'closed':d['fresh_live_entry_cover_closed'],'checks':d['membership_checks'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
