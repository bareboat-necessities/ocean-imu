"""ALT admitted-source private-Mahony startup/invariant attachment.

The shared BRMM proof now closes two facts that the PR #524 handoff still listed
as open: the padded-envelope startup seed lies inside the two-phase Mahony outer
level, and the source-order binary32 Euler step preserves that level.  ALT uses
those facts here without importing the original P2/P3/P4/P5 theorem stages.

This adapter is deliberately narrow.  It proves the private observer remains in
its certified regional set from its physical first-sample seed through startup
and Live, conditional on the declared source-order binary32 arithmetic model.
It does *not* claim the much tighter magnetic accumulation-frame accuracy, FMA
compiler equivalence, or complete deployment arithmetic.

This retained certificate uses its original normalized-direction and seed
premises. It does not cover the commissioned raw-sensor profiles declared in
ou3_alt_startup_sensor_domain.json: their bias/noise seed angles need a new
capture argument. The independent original certificate remains unchanged.
"""
from __future__ import annotations

import ou3_brmm_private_mahony_live_invariant as CONT
import ou3_brmm_private_mahony_discrete_invariant as DISC
from tools.stability.ou3_alt_contraction import live_mahony_regional as LIVE

QUALIFICATION='OU3_ALT_ADMITTED_STARTUP_MAHONY_INVARIANT_V1'


def build():
    c=CONT.build(); cf=CONT.validate(c)
    if cf: raise RuntimeError('continuous admitted-source Mahony certificate invalid: '+repr(cf))
    d=DISC.build(); df=DISC.validate(d)
    if df: raise RuntimeError('binary32 admitted-source Mahony certificate invalid: '+repr(df))
    live=LIVE.build(); lf=LIVE.validate(live)
    if lf: raise RuntimeError('ALT regional Live Mahony adapter invalid: '+repr(lf))
    seed=bool(c['initial_set_inside_invariant'])
    continuous=bool(c['continuous_all_live_PI_invariant_closed'])
    discrete=bool(d['shipping_binary32_discrete_invariant_closed_conditionally_on_source_order'])
    closed=bool(seed and continuous and discrete and live['regional_source_order_binary32_discrete_invariance_closed'])
    return {
      'qualification':QUALIFICATION,
      'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'physical_first_sample_seed_inside_outer_level':seed,
      'startup_and_Live_continuous_outer_invariant_closed':continuous,
      'source_order_binary32_step_preserves_outer_invariant':discrete,
      'same_two_phase_metric_and_outer_level_used_by_ALT_Live': bool(
          live['metric_P']==c['metric_P'] and
          live['invariant_level_C']==c['invariant_level_C']),
      'admitted_source_private_Mahony_startup_to_Live_invariant_closed':closed,
      'all_time_tilt_deg_upper':c['actual_tilt_deg_upper'],
      'tilt_at_150s_deg_upper':c['tilt_at_deployed_horizon_deg_upper'],
      'ultimate_tilt_deg_upper':c['ultimate_tilt_deg_upper'],
      'magnetic_accumulation_frame_accuracy_closed_by_this_invariant':False,
      'compiler_reassociation_or_FMA_closed':False,
      'complete_startup_deployment_arithmetic_closed':False,
      'commissioned_startup_sensor_profiles_covered':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,'ALT_LIVE_PASS':False,'ALT_END_TO_END_PASS':False,
    }


def validate(x):
    f=[]
    if x.get('qualification')!=QUALIFICATION: f.append('qualification mismatch')
    for k in ('physical_first_sample_seed_inside_outer_level',
              'startup_and_Live_continuous_outer_invariant_closed',
              'source_order_binary32_step_preserves_outer_invariant',
              'same_two_phase_metric_and_outer_level_used_by_ALT_Live',
              'admitted_source_private_Mahony_startup_to_Live_invariant_closed'):
        if x.get(k) is not True: f.append(k+' not true')
    for k in ('magnetic_accumulation_frame_accuracy_closed_by_this_invariant',
              'commissioned_startup_sensor_profiles_covered',
              'compiler_reassociation_or_FMA_closed',
              'complete_startup_deployment_arithmetic_closed','storage_search_allowed',
              'ALT_STARTUP_PASS','ALT_LIVE_PASS','ALT_END_TO_END_PASS'):
        if x.get(k) is not False: f.append(k+' not false')
    return f
