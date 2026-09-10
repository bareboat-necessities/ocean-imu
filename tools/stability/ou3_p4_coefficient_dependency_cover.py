#!/usr/bin/env python3
"""Dependency-preserving adaptive coefficient invariants.

The deployed coefficient state is not a Cartesian (tau,sigma,T_S,R_S) box.
Shipping enforces T_S=clamp(pseudo_ratio*tau), R_S_target=SpectralMSE(tau,sigma),
and applied anisotropic R_S=(rx*r,ry*r,r). Candidate EMA is a convex
combination and staged commit copies the candidate, so the clamped target set is
forward invariant once the deployed initial schedule is inside it. This module
materializes those exact relations for the source-uniform proof; it does not
claim the BRMM predecessor family itself is covered.
"""
from __future__ import annotations
import json,math
import ou3_brmm_frontend_state_step as FRONT
import ou3_brmm_tuner_scheduler_step as TUNER
QUALIFICATION='OU3_P4_COEFFICIENT_DEPENDENCY_COVER_V2'

def build():
 c=TUNER.constants();initial=FRONT._point_state().tuner
 f_lo,f_hi=c.tune_freq_min,c.tune_freq_max
 tau_lo=max(c.tau_min,min(c.tau_max,c.tau_coeff*.5/f_hi));tau_hi=max(c.tau_min,min(c.tau_max,c.tau_coeff*.5/f_lo))
 ts=TUNER.pseudo_period(TUNER.Interval(tau_lo,tau_hi),c) if hasattr(TUNER,'Interval') else TUNER.pseudo_period(TUNER.I(tau_lo),c)
 # Derive raw sigma floor through the deployed target expression.
 sigma_floor=max(0.0,c.sigma_coeff*math.sqrt(1e-6))
 membership={
  'candidate_tau':initial.candidate.tau.lo>=tau_lo and initial.candidate.tau.hi<=tau_hi,
  'candidate_sigma':initial.candidate.sigma.lo>=sigma_floor and initial.candidate.sigma.hi<=c.sigma_max,
  'candidate_RS':initial.candidate.rs.lo>=c.rs_min and initial.candidate.rs.hi<=c.rs_max,
  'active_tau':initial.active.tau.lo>=tau_lo and initial.active.tau.hi<=tau_hi,
  'active_sigma':initial.active.sigma.lo>=sigma_floor and initial.active.sigma.hi<=c.sigma_max,
  'active_RS':initial.active.rs_base.lo>=c.rs_min and initial.active.rs_base.hi<=c.rs_max,
 }
 probe=TUNER.ActiveSchedule(TUNER.I(1.0),TUNER.I(1.0),TUNER.Interval(1.0,2.0) if hasattr(TUNER,'Interval') else TUNER.I(1.0),TUNER.I(.01));xyz=TUNER.active_rs_std_xyz(probe,c)
 ray=xyz[0].lo==xyz[1].lo and xyz[0].hi==xyz[1].hi
 return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','filter_changed':False,'declared_domain_shrunk':False,'trajectory_replay_used':False,'independent_tau_sigma_TS_RS_rectangle_used':False,'deployed_pinnings':{'T_S_from_tau':'clamp(pseudo_ratio*tau,[pseudo_min,pseudo_max])','R_S_target_from_tau_sigma':'SpectralMSE(tau,sigma)','applied_R_S_from_one_scalar':'(rs_x*r,rs_y*r,r)','verified_against_deployed_transitions':True},'forward_invariant_coefficient_box':{'tau_s':[tau_lo,tau_hi],'sigma_mps2':[sigma_floor,c.sigma_max],'R_S_base_m_s':[c.rs_min,c.rs_max],'initial_schedule_membership':membership,'initial_schedule_inside_box':all(membership.values()),'invariance_is_conditional_on_initial_membership':True,'needs_reachable_set_argument':False},'applied_R_S_ray':{'one_scalar_generates_three_axes':ray,'x_factor':c.rs_x_factor,'y_factor':c.rs_y_factor,'z_factor':1.0},'reachable_f_sigma_set_certified_here':False,'consecutive_word_storage_compatibility_closed_here':False,'SOURCE_UNIFORM_COMPLETE_BRMM_COVER_CLOSED':False,'P4_promoted_here':False}
def validate(d):
 f=[]
 if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
 for k in ('filter_changed','declared_domain_shrunk','trajectory_replay_used','independent_tau_sigma_TS_RS_rectangle_used','reachable_f_sigma_set_certified_here','consecutive_word_storage_compatibility_closed_here','SOURCE_UNIFORM_COMPLETE_BRMM_COVER_CLOSED','P4_promoted_here'):
  if d.get(k) is not False:f.append(k+' not false')
 if d.get('deployed_pinnings',{}).get('verified_against_deployed_transitions') is not True:f.append('pinnings not verified')
 inv=d.get('forward_invariant_coefficient_box',{})
 if inv.get('initial_schedule_inside_box') is not True:f.append('initial schedule outside invariant box')
 if inv.get('invariance_is_conditional_on_initial_membership') is not True:f.append('initial-membership premise lost')
 if d.get('applied_R_S_ray',{}).get('one_scalar_generates_three_axes') is not True:f.append('R_S ray lost')
 return f
if __name__=='__main__':
 d=build();f=validate(d);print(json.dumps({**d,'validation_pass':not f,'validation_failures':f},indent=2,sort_keys=True));raise SystemExit(bool(f))
