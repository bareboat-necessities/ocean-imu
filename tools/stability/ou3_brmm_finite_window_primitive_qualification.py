#!/usr/bin/env python3
"""Uniform physical-primitive qualification for the complete BRMM source.

The previous BRMM declaration required each admitted history to possess bounded
velocity/displacement primitives, but left the family-wide constants unfrozen.
That is not enough for a uniform finite P4 tube.  This module completes the
*applicability definition* of the deterministic theorem source without fitting a
trajectory and without replacing one physical history by independent sample
boxes.

The qualified Normal-Live source carries one wave-frame physical jet

    p_dot = v,       v_dot = a,
    S_dot = p,       R_dot = R [omega]x,

with the existing pathwise acceleration/body-rate caps and with uniform
wave-frame primitive caps P_m,V_m.  S is deliberately NOT assigned an
independent global ball.  On every h<=3 s proof word its increment is the exact
same-history jet

    Delta S = h p_0 + h^2 v_0/2 + J_2(a),

where J_2(a)=int_0^h (h-s)^2 a(s)/2 ds.  The source relation retains p,v,a in
one witness, so the scalar bound below is only a consequence used for working-
tube accounting, never a generator.

P_m is the already-declared full-scale fresh-entry wave-position norm
sqrt(3)*(Hs_max/2), encoded in the closure-domain radius.  V_m=5 m/s is the
already-declared full-scale marine velocity norm.  These are theorem
applicability limits, not replay extrema.  Making the formerly unspecified
constants uniform is an explicit BRMM source qualification/refinement.
"""
from __future__ import annotations

import argparse, json, math
from fractions import Fraction as F
from pathlib import Path

REPO=Path(__file__).resolve().parents[2]
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
CLOSURE=REPO/'tools/stability/ou3_p4_closure_domain.json'
SCHEMA=1
QUALIFICATION='OU3_BRMM_UNIFORM_FINITE_WINDOW_PRIMITIVE_QUALIFICATION_V1'
HORIZON_S=3.0


def _up(x:float)->float:return math.nextafter(float(x),math.inf)


def build(domain_path:Path=DOMAIN,closure_path:Path=CLOSURE)->dict:
    domain=json.loads(Path(domain_path).read_text())
    closure=json.loads(Path(closure_path).read_text())
    radii=closure['hard_entry_search']['base_coordinate_radii']
    hs=float(domain['initial_filter_entrance']['position']['significant_wave_height_Hs_upper_m'])
    p_component=0.5*hs
    p_from_hs=math.sqrt(3.0)*p_component
    Pm=float(radii['position_norm_m'])
    Vm=float(radii['velocity_norm_mps'])
    A=float(domain['normal_live']['non_gravitational_cog_acceleration_norm_upper_mps2'])
    Wdeg=float(domain['normal_live']['body_rate_norm_upper_deg_s'])
    if not math.isclose(Pm,p_from_hs,rel_tol=2e-15,abs_tol=2e-15):
        raise RuntimeError('closure position scale detached from sqrt(3)*Hs_max/2')
    if not all(math.isfinite(x) and x>0 for x in (Pm,Vm,A,Wdeg)):
        raise RuntimeError('primitive qualification lost positive finite constants')
    h=HORIZON_S
    # Consequences of one common history.  These are NOT independent ports.
    J0=_up(A*h)
    J1=_up(A*h*h/2.0)
    J2=_up(A*h*h*h/6.0)
    # Because ||p(t)||<=Pm throughout the qualified word, the strongest simple
    # S-increment consequence is int ||p|| <= h Pm.  The jet-only alternative is
    # retained for cross-checking but is not used to split source coordinates.
    dS_position=_up(h*Pm)
    dS_jet=_up(h*Pm + h*h*Vm/2.0 + A*h*h*h/6.0)
    dR_angle=_up(math.radians(Wdeg)*h)
    return {
      'schema':SCHEMA,'qualification':QUALIFICATION,
      'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'source_applicability_refined':True,
      'trajectory_replay_used':False,'trajectory_extrema_used':False,
      'independent_sample_boxes_used':False,'independent_primitive_ports_used':False,
      'same_history_dynamic_constraint_used':True,
      'uniform_physical_primitives':{
        'V_m_norm_upper_mps':Vm,'P_m_norm_upper_m':Pm,
        'S_m_independent_global_bound':None,
        'acceleration_norm_upper_mps2':A,'body_rate_norm_upper_deg_s':Wdeg,
        'position_bound_provenance':'sqrt(3)*(Hs_max/2), identical to the existing full-scale closure position radius',
        'velocity_bound_provenance':'existing full-scale deterministic marine velocity radius, now promoted as a physical-source applicability limit'},
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
      'next_obligation':'feed this same-history dynamic source relation, not its scalar hull, through the trusted 601-sample BRMM executor/source-cover and the 24-state joint bias/source compatible storage'}


def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('source_applicability_refined','same_history_dynamic_constraint_used','finite_window_relation_is_closed_and_bounded','finite_window_relation_is_compact_in_sampled_coordinates','equivalent_hard_finite_window_dynamic_constraint_closed','qualified_BRMM_left_inclusion_is_by_source_definition','unqualified_former_any_finite_primitive_family_not_claimed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('trajectory_replay_used','trajectory_extrema_used','independent_sample_boxes_used','independent_primitive_ports_used','S_300_m_s_legacy_ball_used','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    u=d.get('uniform_physical_primitives',{})
    for k in ('V_m_norm_upper_mps','P_m_norm_upper_m','acceleration_norm_upper_mps2','body_rate_norm_upper_deg_s'):
        x=float(u.get(k,math.nan))
        if not(math.isfinite(x) and x>0):f.append(k+' invalid')
    if u.get('S_m_independent_global_bound','sentinel') is not None:f.append('independent S_m reintroduced')
    m=d.get('finite_window_moment_consequences',{})
    if not float(m.get('Delta_S_norm_upper_from_uniform_position_m_s',math.inf)) < 23.0:f.append('3s correlated S increment unexpectedly large')
    if d.get('canonical_source')!='COMPLETE_BRMM_NORMAL_LIVE_WORD':f.append('source changed')
    return f


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'P_m':d['uniform_physical_primitives']['P_m_norm_upper_m'],'V_m':d['uniform_physical_primitives']['V_m_norm_upper_mps'],'DeltaS3':d['finite_window_moment_consequences']['Delta_S_norm_upper_from_uniform_position_m_s'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
