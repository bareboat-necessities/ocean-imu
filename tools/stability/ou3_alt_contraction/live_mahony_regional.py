#!/usr/bin/env python3
"""Regional private-Mahony invariant for the ALT Live theorem.

The shared Mahony certificates currently conflate two obligations:
  (1) inward/discrete invariance of a declared Live set; and
  (2) startup's initial seed belonging to that set.
ALT needs (1) now and proves (2) only after a quantitative Live basin exists.

This module reuses the SAME metric, BRMM direction decomposition, rational
boundary cover and source-order binary32/Euler charges.  It removes only the
startup-membership premise from the *Live* theorem.  The theorem here is
conditional: if a Live predecessor lies in V(z)<=C, one literal 5 ms Mahony
step stays in V(z)<=C and in the 87-degree chart under the declared source-order
binary32 model.  No claim about startup entry or compiler reassociation/FMA is
made.
"""
from __future__ import annotations
import json,math
from pathlib import Path

import ou3_brmm_private_mahony_live_invariant as CONT
import ou3_brmm_private_mahony_discrete_invariant as DISC
import ou3_fast_inv_sqrt_interval as FINV

REPO=Path(__file__).resolve().parents[3]
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
QUALIFICATION='OU3_ALT_REGIONAL_LIVE_MAHONY_INVARIANT_V1'


def build():
    d=json.loads(DOMAIN.read_text());c=CONT._source_constants(d);sector=CONT._sector_lower(c['chart_deg'])
    boundary=CONT._verify_boundary(c,sector.lo)
    C=CONT.C_LEVEL;ztheta_sq=C*CONT.P22/CONT.DET_P
    theta_upper=math.sqrt(ztheta_sq)+0.1*c['xi'];chart_rad=math.radians(c['chart_deg']);chart_contained=theta_upper<chart_rad
    regional_cont=bool(boundary['closed'] and boundary['strict_inward_margin_lower']>0 and chart_contained)
    if not regional_cont:raise RuntimeError('regional continuous Mahony invariant failed')

    dt=float(d['configured_runtime']['imu_dt_s']);body=math.radians(float(d['normal_live']['body_rate_norm_upper_deg_s']))
    shell=FINV.all_positive_normal_normalized_norm2_enclosure();qn_lo=math.sqrt(shell.lo);qn_hi=math.sqrt(shell.hi)
    xi=c['xi'];mean=c['mean_chord'];ztheta=math.sqrt(C*CONT.P22/CONT.DET_P);zbeta=math.sqrt(C/CONT.DET_P);beta=zbeta+0.01*xi
    half_g=0.5*qn_hi*qn_hi;half_err=qn_hi*half_g;omega=body+beta+0.2*half_err;half_rate=0.5*dt*omega
    raw_sum_abs=qn_hi*(1+3*half_rate);raw_component_round=DISC.gamma(24)*raw_sum_abs;raw_vector_round=2*raw_component_round;raw_angle_round=2*raw_vector_round/qn_lo
    norm_mult_vector_round=2*DISC.U*qn_hi;norm_angle_round=2*norm_mult_vector_round/qn_lo
    feedback_rate_round=DISC.gamma(24)*omega;attitude_rate_fp=(raw_angle_round+norm_angle_round)/dt+feedback_rate_round
    initial_bias=float(d['startup']['initial_tangent_gyro_bias_norm_upper_rad_s']);integral_mag=beta+initial_bias;ki_increment=.02*half_err*dt
    integral_step_round=DISC.gamma(8)*(integral_mag+ki_increment);bias_rate_fp=integral_step_round/dt
    fp_support=attitude_rate_fp+math.hypot(6.5,11.0)*bias_rate_fp
    continuous_margin=float(boundary['strict_inward_margin_lower']);robust_margin=continuous_margin-2*fp_support
    smin=sector.lo;dg=float(d['startup']['effective_deterministic_gyro_transport_disturbance_upper_rad_s'])+attitude_rate_fp;db=float(d['startup']['effective_deterministic_bias_transport_disturbance_upper_rad_s2'])+bias_rate_fp
    ftheta=.1*ztheta+zbeta+max(abs(.01*smin-.01),0)*xi+.1*mean+dg;fbeta=.01*ztheta+.001*xi+.01*mean+db
    Rf1=ftheta+6.5*fbeta;Rf2=11*fbeta;euler_quadratic=dt*dt*(Rf1*Rf1+Rf2*Rf2)
    first_order=2*CONT.SQRT_C*robust_margin*dt;discrete_margin=first_order-euler_quadratic
    chart_margin=chart_rad-theta_upper;round_chart_margin=chart_margin-(raw_angle_round+norm_angle_round)
    regional_binary32=bool(robust_margin>0 and discrete_margin>0 and round_chart_margin>0)
    return {
      'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'regional_entry_set':'V(z)<=C with z=x-B*xi','metric_P':[[CONT.P11,CONT.P12],[CONT.P12,CONT.P22]],'invariant_level_C':C,
      'same_BRMM_direction_primitive_and_transport_forcing_retained':True,'boundary_cells_per_chart':CONT.CELLS_PER_CHART,
      'continuous_boundary_margin_lower':continuous_margin,'regional_continuous_invariance_closed':regional_cont,
      'actual_tilt_rad_upper_on_invariant':theta_upper,'chart_radius_rad':chart_rad,
      'source_order_binary32_support_charge':fp_support,'binary32_robust_boundary_margin_lower':robust_margin,
      'forward_euler_quadratic_V_upper':euler_quadratic,'discrete_V_margin_lower':discrete_margin,'one_step_chart_round_margin_rad':round_chart_margin,
      'regional_source_order_binary32_discrete_invariance_closed':regional_binary32,
      'startup_initial_seed_membership_required_for_Live_invariance':False,'startup_entry_into_regional_set_closed_here':False,
      'compiler_reassociation_or_FMA_closed':False,'trajectory_replay_used':False,'filter_changed':False,
      'ALT_LIVE_PASS':False,
      'next_obligation':'use this regional invariant, rather than startup capture, as the Mahony component of the Normal-Live frontend predecessor set; preserve startup entry as a later basin-capture theorem',
    }

def validate(x):
    f=[]
    if x.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('same_BRMM_direction_primitive_and_transport_forcing_retained','regional_continuous_invariance_closed','regional_source_order_binary32_discrete_invariance_closed'):
        if x.get(k) is not True:f.append(k+' not true')
    for k in ('startup_initial_seed_membership_required_for_Live_invariance','startup_entry_into_regional_set_closed_here','compiler_reassociation_or_FMA_closed','trajectory_replay_used','filter_changed','ALT_LIVE_PASS'):
        if x.get(k) is not False:f.append(k+' not false')
    for k in ('continuous_boundary_margin_lower','binary32_robust_boundary_margin_lower','discrete_V_margin_lower','one_step_chart_round_margin_rad'):
        if not float(x.get(k,0))>0:f.append(k+' not positive')
    return f
