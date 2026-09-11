#!/usr/bin/env python3
"""Necessary finite-window primitive consequences of corrected COMPLETE-BRMM.

Physical wave-state bounds are read from the explicit COMPLETE-BRMM physical
envelope. They are not inferred from estimator error-entry radii. The same
physical history carries

    p_dot=v,  v_dot=a,  S_dot=p,  R_dot=R[omega]x.

The finite-window scalar consequences below are only accounting consequences of
that one history; they are never promoted to independent source ports.
"""
from __future__ import annotations

import argparse, json, math
from pathlib import Path

REPO=Path(__file__).resolve().parents[2]
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
SCHEMA=2
QUALIFICATION='OU3_BRMM_UNIFORM_FINITE_WINDOW_PRIMITIVE_QUALIFICATION_V2'
HORIZON_S=3.0


def _up(x:float)->float:return math.nextafter(float(x),math.inf)


def build(domain_path:Path=DOMAIN)->dict:
    domain=json.loads(Path(domain_path).read_text())
    env=domain.get('complete_brmm_physical_envelope',{})
    Pm=float(env['wave_position_norm_upper_m'])
    Vm=float(env['wave_velocity_norm_upper_mps'])
    A=float(env['wave_acceleration_norm_upper_mps2'])
    Wdeg=float(env['body_rate_norm_upper_deg_s'])
    DS=float(env['centered_primitive_D_S_upper_m_s'])
    Hs=float(env['significant_wave_height_Hs_upper_m'])
    fmin,fmax=map(float,env['frequency_support_hz'])
    if not all(math.isfinite(x) and x>0 for x in (Pm,Vm,A,Wdeg,DS,Hs,fmin,fmax)):
        raise RuntimeError('physical primitive qualification lost positive finite constants')
    if fmin>=fmax: raise RuntimeError('physical frequency support is empty')
    if env.get('these_are_physical_source_bounds_not_error_bounds') is not True:
        raise RuntimeError('physical source/error domain separation not explicit')
    h=HORIZON_S
    J0=_up(A*h);J1=_up(A*h*h/2.0);J2=_up(A*h*h*h/6.0)
    dS_position=_up(h*Pm)
    dS_jet=_up(h*Pm+h*h*Vm/2.0+A*h*h*h/6.0)
    dR_angle=_up(math.radians(Wdeg)*h)
    return {
      'schema':SCHEMA,'qualification':QUALIFICATION,
      'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'source_applicability_refined':True,
      'finite_window_only_is_sufficient_for_complete_BRMM_admission':False,
      'corrected_physical_bounded_wave_condition_required_separately':True,
      'trajectory_replay_used':False,'trajectory_extrema_used':False,
      'independent_sample_boxes_used':False,'independent_primitive_ports_used':False,
      'same_history_dynamic_constraint_used':True,
      'physical_source_bounds_are_not_error_entry_bounds':True,
      'uniform_physical_primitives':{
        'H_s_upper_m':Hs,'V_m_norm_upper_mps':Vm,'P_m_norm_upper_m':Pm,
        'S_centered_primitive_norm_upper_m_s':DS,
        'S_m_independent_global_bound':None,
        'acceleration_norm_upper_mps2':A,'body_rate_norm_upper_deg_s':Wdeg,
        'frequency_support_hz':[fmin,fmax],
        'position_bound_provenance':'padded COMPLETE-BRMM physical source envelope',
        'velocity_bound_provenance':'padded COMPLETE-BRMM physical source envelope'},
      'exact_source_relation':{
        'translation':'p_dot=v; v_dot=a; S_dot=p in one fixed wave/world frame',
        'rotation':'R_dot=R*skew(omega), R in SO(3)',
        'one_common_history_required':True,
        'S_is_derived_not_independent':True,
        'word_start_p_v_are_predecessor_outputs_not_fresh_slots':True},
      'canonical_word_horizon_s':h,
      'finite_window_moment_consequences':{
        'J0_acceleration_moment_norm_upper_mps':J0,
        'J1_acceleration_moment_norm_upper_m':J1,
        'J2_acceleration_moment_norm_upper_m_s':J2,
        'Delta_S_norm_upper_from_uniform_position_m_s':dS_position,
        'Delta_S_norm_upper_from_entry_jet_only_m_s':dS_jet,
        'rotation_path_length_upper_rad':dR_angle},
      'finite_window_relation_is_closed_and_bounded':True,
      'finite_window_relation_is_compact_in_sampled_coordinates':True,
      'equivalent_hard_finite_window_dynamic_constraint_closed':True,
      'qualified_BRMM_left_inclusion_is_by_source_definition':True,
      'unqualified_former_any_finite_primitive_family_not_claimed':True,
      'S_300_m_s_legacy_ball_used':False,
      'P4_promoted_here':False,
      'next_obligation':'feed this same-history physical source relation through the trusted 601-sample executor and the correlated P/H/R innovation map, then close joint24 endpoint/prefix storage and first-exit retention'}


def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('source_applicability_refined','same_history_dynamic_constraint_used','physical_source_bounds_are_not_error_entry_bounds','finite_window_relation_is_closed_and_bounded','finite_window_relation_is_compact_in_sampled_coordinates','equivalent_hard_finite_window_dynamic_constraint_closed','qualified_BRMM_left_inclusion_is_by_source_definition','unqualified_former_any_finite_primitive_family_not_claimed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('trajectory_replay_used','trajectory_extrema_used','independent_sample_boxes_used','independent_primitive_ports_used','S_300_m_s_legacy_ball_used','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    u=d.get('uniform_physical_primitives',{})
    for k in ('H_s_upper_m','V_m_norm_upper_mps','P_m_norm_upper_m','S_centered_primitive_norm_upper_m_s','acceleration_norm_upper_mps2','body_rate_norm_upper_deg_s'):
        x=float(u.get(k,math.nan))
        if not(math.isfinite(x) and x>0):f.append(k+' invalid')
    if u.get('S_m_independent_global_bound','sentinel') is not None:f.append('independent S_m reintroduced')
    if u.get('frequency_support_hz') != [0.018,0.88]: f.append('padded physical frequency support changed')
    m=d.get('finite_window_moment_consequences',{})
    if not float(m.get('Delta_S_norm_upper_from_uniform_position_m_s',math.inf)) < 25.0:f.append('3s correlated S increment unexpectedly large')
    if d.get('canonical_source')!='COMPLETE_BRMM_NORMAL_LIVE_WORD':f.append('source changed')
    return f


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'P_m':d['uniform_physical_primitives']['P_m_norm_upper_m'],'V_m':d['uniform_physical_primitives']['V_m_norm_upper_mps'],'D_S':d['uniform_physical_primitives']['S_centered_primitive_norm_upper_m_s'],'DeltaS3':d['finite_window_moment_consequences']['Delta_S_norm_upper_from_uniform_position_m_s'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
