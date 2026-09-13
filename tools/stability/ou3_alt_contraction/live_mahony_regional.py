#!/usr/bin/env python3
"""Regional private-Mahony invariant for the ALT Live theorem.

Main now carries a stronger two-phase V3 Mahony certificate and a matching
source-order binary32/discrete certificate.  ALT needs only their Live-invariant
implication here: membership in the certified outer level is a predecessor
hypothesis; startup membership/capture is proved separately.

This adapter therefore consumes the CURRENT main certificates rather than
repeating their metric arithmetic.  In particular it does not retain the
retired p=6.5,c=11 metric or call a private boundary routine with stale
arguments.  The V3 metric, outer level, BRMM forcing relation, binary32 charge,
Euler quadratic term and chart margin are all read from the validated shared
certificates.
"""
from __future__ import annotations

import ou3_brmm_private_mahony_live_invariant as CONT
import ou3_brmm_private_mahony_discrete_invariant as DISC

QUALIFICATION='OU3_ALT_REGIONAL_LIVE_MAHONY_INVARIANT_V2'


def build():
    cont=CONT.build(); cf=CONT.validate(cont)
    if cf: raise RuntimeError('current continuous Mahony certificate invalid: '+repr(cf))
    disc=DISC.build(); df=DISC.validate(disc)
    if df: raise RuntimeError('current source-order binary32 Mahony certificate invalid: '+repr(df))

    # The shared theorem proves more than ALT uses here (including one startup
    # seed and two-phase capture).  Dropping those premises is sound for the
    # conditional regional statement: if a Live predecessor is already in the
    # outer invariant, the same boundary/discrete argument retains it.
    continuous_margin=float(cont['boundary_validation']['strict_inward_margin_lower'])
    robust_margin=float(disc['continuous_boundary_margin_after_binary32'])
    discrete_margin=float(disc['discrete_V_margin_lower'])
    chart_margin=float(disc['one_step_chart_round_margin_rad'])
    closed=bool(cont['continuous_all_live_PI_invariant_closed'])
    binary32=bool(disc['shipping_binary32_discrete_invariant_closed_conditionally_on_source_order'])
    if not (closed and binary32 and continuous_margin>0 and robust_margin>0 and
            discrete_margin>0 and chart_margin>0):
        raise RuntimeError('regional Mahony implication did not close on current shared certificates')

    return {
      'qualification':QUALIFICATION,
      'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'regional_entry_set':'V(z)<=C_outer with z=x-B*xi',
      'metric_P':cont['metric_P'],
      'metric_cholesky_R':cont['metric_cholesky_R'],
      'invariant_level_C':cont['invariant_level_C'],
      'sqrt_C':cont['sqrt_C'],
      'same_BRMM_direction_primitive_and_transport_forcing_retained':True,
      'boundary_cells_per_chart':cont['boundary_validation']['cells_per_chart'],
      'continuous_boundary_margin_lower':continuous_margin,
      'regional_continuous_invariance_closed':closed,
      'actual_tilt_rad_upper_on_invariant':cont['actual_tilt_rad_upper'],
      'chart_radius_rad':cont['chart_radius_rad'],
      'source_order_binary32_support_charge':disc['binary32_support_charge'],
      'binary32_robust_boundary_margin_lower':robust_margin,
      'forward_euler_quadratic_V_upper':disc['forward_euler_quadratic_V_upper'],
      'discrete_V_margin_lower':discrete_margin,
      'one_step_chart_round_margin_rad':chart_margin,
      'regional_source_order_binary32_discrete_invariance_closed':binary32,
      'current_two_phase_main_metric_consumed':True,
      'startup_initial_seed_membership_required_for_Live_invariance':False,
      'startup_entry_into_regional_set_closed_here':False,
      'compiler_reassociation_or_FMA_closed':False,
      'trajectory_replay_used':False,
      'filter_changed':False,
      'ALT_LIVE_PASS':False,
      'next_obligation':'use the current V3 regional invariant as the Mahony component of the Normal-Live predecessor set; keep startup entry/capture and compiler-mode qualification as separate obligations',
    }


def validate(x):
    f=[]
    if x.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('same_BRMM_direction_primitive_and_transport_forcing_retained',
              'regional_continuous_invariance_closed',
              'regional_source_order_binary32_discrete_invariance_closed',
              'current_two_phase_main_metric_consumed'):
        if x.get(k) is not True:f.append(k+' not true')
    for k in ('startup_initial_seed_membership_required_for_Live_invariance',
              'startup_entry_into_regional_set_closed_here',
              'compiler_reassociation_or_FMA_closed','trajectory_replay_used',
              'filter_changed','ALT_LIVE_PASS'):
        if x.get(k) is not False:f.append(k+' not false')
    for k in ('continuous_boundary_margin_lower','binary32_robust_boundary_margin_lower',
              'discrete_V_margin_lower','one_step_chart_round_margin_rad'):
        if not float(x.get(k,0))>0:f.append(k+' not positive')
    return f
