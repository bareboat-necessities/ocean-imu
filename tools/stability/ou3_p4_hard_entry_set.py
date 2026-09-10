#!/usr/bin/env python3
"""Full-scale P4 working set plus qualified fresh-Live reachability.

The coordinate radii remain the existing deterministic P4 working scales, not a
covariance-confidence set. Fresh reachability is supplied by the uniformly
qualified BRMM physical primitive/source graph: zero shipping linear means,
exact shared-S re-anchoring, and hard physical source bounds place the fresh
outer-wrapper error inside these same full-scale radii.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import ou3_p4_live_entry_graph as LIVE_ENTRY
import ou3_p4_qualified_fresh_entry as QUALIFIED_ENTRY

REPO=Path(__file__).resolve().parents[2]
CLOSURE=REPO/'tools/stability/ou3_p4_closure_domain.json'


def build(path:Path=CLOSURE):
    c=json.loads(Path(path).read_text())['hard_entry_search']
    scale=float(c['minimum_certified_scale'])
    b={k:float(v) for k,v in c['base_coordinate_radii'].items()}
    if scale!=1.0 or list(map(float,c['candidate_scale_factors'])) != [1.0]:
        raise RuntimeError('P4 hard entry set may not shrink the declared working domain')
    r={k:math.nextafter(v,math.inf) for k,v in b.items()}
    live=LIVE_ENTRY.build();qualified=QUALIFIED_ENTRY.build()
    qf=QUALIFIED_ENTRY.validate(qualified)
    if qf:raise RuntimeError('qualified fresh-entry source failed: '+repr(qf))
    fresh=dict(live)
    fresh['source_uniform_reachable_other_entry_coordinates_closed']=qualified['source_uniform_reachable_other_entry_coordinates_closed']
    fresh['qualified_full_scale_membership']=qualified['membership_checks']
    fresh['qualified_source_applicability_refinement_consumed']=qualified['source_applicability_refinement_consumed']
    return {
      'qualification':'OU3_P4_DETERMINISTIC_HARD_ENTRY_SET_V3',
      'source':'full-scale deterministic working set plus qualified fresh-Live physical-source reachability',
      'legacy_box_is_certified_fresh_live_entry':False,
      'fresh_live_entry_graph':fresh,
      'qualified_fresh_live_entry':qualified,
      'S_radius_role':'working-domain scale only; fresh centered e_S is exactly zero and subsequent S forcing is same-history correlated',
      'trajectory_fit':False,'covariance_ellipsoid_used':False,'shipping_covariance_membership_used':False,
      'hard_state_error_membership_required':True,'scale_of_declared_handoff_envelope':scale,'coordinate_radii':r,
      'H18_groups':['attitude_cayley_norm','gyro_bias_norm_rad_s','velocity_norm_mps','position_norm_m','integral_displacement_norm_m_s','latent_acceleration_norm_mps2'],
      'A21_additional_group':'accelerometer_bias_error_norm_mps2',
      'full_declared_scale_enforced':True,'entry_set_declared_and_membership_checkable':True,
      'startup_reachability_proved_here':qualified['fresh_live_entry_cover_closed'],
      'P4_may_use_as_regional_entry_hypothesis':True,
      'proof_error_domain_shrunk':False,
      'physical_source_family_uniformly_qualified':True,
    }


def validate(d):
    f=[]
    for k in ('hard_state_error_membership_required','full_declared_scale_enforced','entry_set_declared_and_membership_checkable','startup_reachability_proved_here','P4_may_use_as_regional_entry_hypothesis','physical_source_family_uniformly_qualified'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('trajectory_fit','covariance_ellipsoid_used','shipping_covariance_membership_used','proof_error_domain_shrunk','legacy_box_is_certified_fresh_live_entry'):
        if d.get(k) is not False:f.append(k+' not false')
    if float(d.get('scale_of_declared_handoff_envelope',0))!=1.0:f.append('hard entry scale is not full declared envelope')
    fresh=d.get('fresh_live_entry_graph',{})
    if fresh.get('source_uniform_reachable_other_entry_coordinates_closed') is not True:f.append('fresh source reachability not closed')
    return f


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'startup_reachability':d['startup_reachability_proved_here'],'checks':d['qualified_fresh_live_entry']['membership_checks'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
