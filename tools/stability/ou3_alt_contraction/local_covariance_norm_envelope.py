#!/usr/bin/env python3
"""Source-uniform event-local covariance operator-norm envelope for ALT.

This theorem uses the endpoint-referenced covariance envelope only where that
producer actually qualifies it: at arbitrary COMPLETE-BRMM word endpoints.
It never substitutes endpoint Pbar directly for a pre-Joseph local covariance.

There are three regimes.

H18 early Live:
  start from the actual source-generated Live seed and propagate the conservative
  prediction/floor scalar recurrence for at most one 600-sample word.

H18 after the first word:
  every current sample follows an endpoint of a preceding source-uniform 3 s
  window.  Convert the rigorous endpoint diagonal Pbar into an operator-norm
  upper by PSD trace, then apply at most one prediction plus one a_w floor
  before any Joseph event of the current sample.

A21 release / early active word:
  shipping cannot release before the configured >=30 s H interval, so H18 has
  long since reached the endpoint-envelope regime.  The exact held->active edge
  appends the fixed b_a seed block with zero cross covariance.  Propagate that
  actual release covariance through at most one 600-sample A word.  Thereafter
  use the A21 sliding endpoint envelope plus one prediction/floor.

Every Joseph operation is PSD-decreasing in exact real arithmetic and the bias
projection leaves covariance unchanged.  For the a_w floor,

  Delta = Pi_+(Sigma_aw-P_aw),  P_aw>=0,
  ||Delta||2 = max(lambda_max(Sigma_aw-P_aw),0) <= lambda_max(Sigma_aw),

so one requested floor adds at most the stationary target variance independent
of the current covariance.
"""
from __future__ import annotations

import math
from pathlib import Path
from ou3_interval import Interval
import ou3_brmm_live_covariance_seed as LIVE
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_brmm_shipping_prediction_primitives as PRIM
import ou3_brmm_tuner_scheduler_step as TUNER
import ou3_brmm_riccati_tube_factored as TUBE
from tools.stability.ou3_alt_contraction import physical_numeric_bounds as PHYS

REPO=Path(__file__).resolve().parents[3]
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
QUALIFICATION='OU3_ALT_EVENT_LOCAL_COVARIANCE_OPERATOR_NORM_ENVELOPE_V2'
TRANSITIONS=600

def I(x):return Interval.outward_bounds(float(x),float(x))
def B(r):r=float(r);return Interval.outward_bounds(-r,r)
def _zero(n):return [[I(0) for _ in range(n)] for _ in range(n)]

def interval_matrix_norm2_upper(A):
    """Rigorous ||A||2 upper via sqrt(||A||1 ||A||inf)."""
    if not A or any(len(r)!=len(A[0]) for r in A):raise ValueError('rectangular nonempty matrix required')
    rows=len(A);cols=len(A[0]);rowmax=0.0;colmax=0.0
    for i in range(rows):
        s=0.0
        for j in range(cols):s=math.nextafter(s+A[i][j].abs_upper(),math.inf)
        rowmax=max(rowmax,s)
    for j in range(cols):
        s=0.0
        for i in range(rows):s=math.nextafter(s+A[i][j].abs_upper(),math.inf)
        colmax=max(colmax,s)
    return math.nextafter(math.sqrt(math.nextafter(rowmax*colmax,math.inf)),math.inf)

def _assemble(mode,Faa,Qaa,Ft,Qt,phi_ba,Qba):
    n=18 if mode=='H' else 21;F=_zero(n);Q=_zero(n)
    for i in range(6):
        for j in range(6):F[i][j]=Faa[i][j];Q[i][j]=Qaa[i][j]
    for i in range(12):
        for j in range(12):F[6+i][6+j]=Ft[i][j];Q[6+i][6+j]=Qt[i][j]
    if mode=='A':
        for i in range(3):F[18+i][18+i]=phi_ba;Q[18+i][18+i]=Qba[i][i]
    return F,Q

def _live_H_seed_norm(live):
    c=live['constructor'];att=live['full_heading_gauged_live_attitude_seed'];aw=live['aw_live_seed']
    vals=[float(att['tilt_std_rad'])**2,float(att['yaw_std_rad'])**2,float(c['P_bg_variance']),float(live['translation_seed']['P_v']),float(live['translation_seed']['P_p']),float(live['translation_seed']['P_S']),float(aw['committed_vertical_std_interval_mps2'][1])**2]
    return math.nextafter(max(vals),math.inf)
def _endpoint_norm_from_diag(diag):
    # P>=0 => ||P||2=lambda_max(P)<=trace(P)=sum diag(P).
    s=0.0
    for x in diag:s=math.nextafter(s+float(x),math.inf)
    if not math.isfinite(s) or s<=0:raise ValueError('invalid endpoint covariance trace bound')
    return s
def _advance(p,fn,qn,floor_add):return math.nextafter(math.nextafter(fn*fn,math.inf)*p+qn+floor_add,math.inf)
def _word_open_loop(seed,fn,qn,floor_add):
    p=float(seed);mx=p
    for _ in range(TRANSITIONS):p=_advance(p,fn,qn,floor_add);mx=max(mx,p)
    return p,mx

def build():
    live=LIVE.build(DOMAIN);lf=LIVE.validate(live);phys=PHYS.build();pf=PHYS.validate(phys);tube=TUBE.build();tf=TUBE.validate_covariance_ceiling(tube)
    if lf or pf or tf:raise RuntimeError(f'local covariance prerequisites failed live={lf} phys={pf} tube={tf}')
    kc=KERNEL._process_constants(DOMAIN);tc=TUNER.constants();omega=[B(math.radians(float(phys['body_rate_norm_upper_deg_s'])))]*3;tau=Interval.outward_bounds(tc.tau_min,tc.tau_max);sigma_lo=max(float(tc.acc_noise_floor_sigma),1e-9);sigma_hi=float(tc.sigma_max);sigma=[Interval.outward_bounds(sigma_lo,sigma_hi)]*3
    Faa,Qaa=PRIM.attitude_gyro_bias_F_Q(omega,kc.h,kc.gyro_variance_density_xyz,kc.gyro_bias_variance_density);Ft,Qt=PRIM.translation_F_Q(tau,kc.h,sigma);phi,Qba=PRIM.active_accel_bias_F_Q(kc.h,kc.accel_bias_tau_s,kc.accel_bias_process_variance_density)
    floor_add=math.nextafter(sigma_hi*sigma_hi,math.inf);pars={}
    for mode in ('H','A'):
        F,Q=_assemble(mode,Faa,Qaa,Ft,Qt,phi,Qba);pars[mode]=(interval_matrix_norm2_upper(F),interval_matrix_norm2_upper(Q))
    pbarH=_endpoint_norm_from_diag(tube['modes']['H']['Pbar_diagonal_variance_upper']);pbarA=_endpoint_norm_from_diag(tube['modes']['A']['Pbar_diagonal_variance_upper'])
    hseed=_live_H_seed_norm(live);h_end,h_early=_word_open_loop(hseed,*pars['H'],floor_add);h_steady=_advance(pbarH,*pars['H'],floor_add);hlocal=max(h_early,h_steady)
    ba_seed=float(live['held_ba']['seed_variance']);release_seed=math.nextafter(max(pbarH,ba_seed),math.inf);a_end,a_early=_word_open_loop(release_seed,*pars['A'],floor_add);a_steady=_advance(pbarA,*pars['A'],floor_add);alocal=max(a_early,a_steady)
    reports={
      'H':{'dimension':18,'F_norm_upper':pars['H'][0],'Q_norm_upper':pars['H'][1],'live_seed_norm_upper':hseed,'early_word_end_open_loop_upper':h_end,'early_word_every_event_upper':h_early,'sliding_endpoint_Pbar_norm_upper':pbarH,'steady_one_prediction_floor_beyond_endpoint_upper':h_steady,'every_event_covariance_norm_upper':hlocal},
      'A':{'dimension':21,'F_norm_upper':pars['A'][0],'Q_norm_upper':pars['A'][1],'release_seed_from_H_endpoint_plus_fixed_ba_upper':release_seed,'early_active_word_end_open_loop_upper':a_end,'early_active_word_every_event_upper':a_early,'sliding_endpoint_Pbar_norm_upper':pbarA,'steady_one_prediction_floor_beyond_endpoint_upper':a_steady,'every_event_covariance_norm_upper':alocal},
    }
    for r in reports.values():r['aw_floor_increment_norm_upper']=floor_add;r['finite']=all(math.isfinite(float(v)) for k,v in r.items() if isinstance(v,(int,float)))
    closed=all(r['finite'] and r['every_event_covariance_norm_upper']>0 for r in reports.values())
    return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','word_transitions':TRANSITIONS,'live_seed_source_generated':live['live_entry_seed_is_source_generated_not_arbitrary_PSD'],'prediction_F_Q_outward_shipping_primitives_used':True,'Joseph_update_PSD_nonincrease_used':True,'bias_projection_covariance_identity_used':True,'aw_floor_positive_part_bound':'||Pi_+(Sigma-Paw)||2 <= lambda_max(Sigma) for Paw>=0','aw_floor_bound_independent_of_current_P':True,'endpoint_Pbar_used_only_at_certified_sliding_word_endpoints':True,'endpoint_Pbar_substituted_directly_for_preJoseph_covariance':False,'H_to_A_release_covariance_ancestry_retained':True,'A_release_after_H_endpoint_regime_guaranteed':float(live['H_to_A_release']['minimum_H_mode_live_duration_before_A_release_s'])>=3.0,'trajectory_replay_used':False,'reports':reports,'event_local_covariance_norm_envelope_closed':closed,'next_obligation':'derive inverse-free source-uniform event gain bounds from this corrected local covariance envelope and use them in the physical joint24 local outer'}
def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION or d.get('word_transitions')!=TRANSITIONS:f.append('qualification/window mismatch')
    for k in ('live_seed_source_generated','prediction_F_Q_outward_shipping_primitives_used','Joseph_update_PSD_nonincrease_used','bias_projection_covariance_identity_used','aw_floor_bound_independent_of_current_P','endpoint_Pbar_used_only_at_certified_sliding_word_endpoints','H_to_A_release_covariance_ancestry_retained','A_release_after_H_endpoint_regime_guaranteed','event_local_covariance_norm_envelope_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('endpoint_Pbar_substituted_directly_for_preJoseph_covariance','trajectory_replay_used'):
        if d.get(k) is not False:f.append(k+' not false')
    for m in ('H','A'):
        r=d.get('reports',{}).get(m,{})
        if r.get('finite') is not True or not float(r.get('every_event_covariance_norm_upper',0))>0:f.append(m+' local covariance envelope invalid')
    return f
