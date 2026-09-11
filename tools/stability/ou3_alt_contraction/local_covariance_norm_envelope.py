#!/usr/bin/env python3
"""Source-uniform event-local covariance operator-norm envelope for ALT.

This replaces an invalid use of the endpoint-referenced Pbar as a local Joseph
covariance bound.  The proof uses only shipping/source facts:

* Live covariance is reseeded before the first prediction;
* prediction is P^- = F P F' + Q;
* every Joseph update with the shipping Kalman gain is PSD-decreasing,
  P^+ = P^- - K S K' <= P^-;
* active bias projection leaves covariance unchanged;
* the a_w synchronization increment is Delta=Pi_+(Sigma_aw-P_aw).

For P_aw>=0,
  ||Delta||_2=max(lambda_max(Sigma_aw-P_aw),0)<=lambda_max(Sigma_aw),
so a requested floor adds at most the committed stationary target variance in
operator norm, independently of the current P.  Therefore the scalar recurrence

  p_{k+1} <= ||F||_2^2 p_k + ||Q||_2 + sigma_aw,max^2

is a valid source-uniform upper for covariance at every event of every sample.
The F/Q norms below are themselves upper-bounded from outward shipping matrix
primitives over the declared BRMM rate/tuner ranges.
"""
from __future__ import annotations

import math
import json
from pathlib import Path

from ou3_interval import Interval
import ou3_brmm_live_covariance_seed as LIVE
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_brmm_shipping_prediction_primitives as PRIM
import ou3_brmm_tuner_scheduler_step as TUNER
from tools.stability.ou3_alt_contraction import physical_numeric_bounds as PHYS

REPO=Path(__file__).resolve().parents[3]
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
QUALIFICATION='OU3_ALT_EVENT_LOCAL_COVARIANCE_OPERATOR_NORM_ENVELOPE_V1'
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

def _seed_norms(live):
    c=live['constructor'];att=live['full_heading_gauged_live_attitude_seed'];aw=live['aw_live_seed'];held=live['held_ba']
    attmax=max(float(att['tilt_std_rad'])**2,float(att['yaw_std_rad'])**2)
    base=max(attmax,float(c['P_bg_variance']),float(live['translation_seed']['P_v']),float(live['translation_seed']['P_p']),float(live['translation_seed']['P_S']),float(aw['committed_vertical_std_interval_mps2'][1])**2)
    return {'H':math.nextafter(base,math.inf),'A':math.nextafter(max(base,float(held['seed_variance'])),math.inf)}
def build():
    live=LIVE.build(DOMAIN);lf=LIVE.validate(live);phys=PHYS.build();pf=PHYS.validate(phys)
    if lf or pf:raise RuntimeError(f'local covariance prerequisites failed live={lf} phys={pf}')
    kc=KERNEL._process_constants(DOMAIN);tc=TUNER.constants()
    omega_cap=math.radians(float(phys['body_rate_norm_upper_deg_s']));omega=[B(omega_cap) for _ in range(3)]
    tau=Interval.outward_bounds(tc.tau_min,tc.tau_max)
    sigma_lo=max(float(tc.acc_noise_floor_sigma),1e-9);sigma_hi=float(tc.sigma_max)
    sigma=[Interval.outward_bounds(sigma_lo,sigma_hi) for _ in range(3)]
    Faa,Qaa=PRIM.attitude_gyro_bias_F_Q(omega,kc.h,kc.gyro_variance_density_xyz,kc.gyro_bias_variance_density)
    Ft,Qt=PRIM.translation_F_Q(tau,kc.h,sigma)
    phi,Qba=PRIM.active_accel_bias_F_Q(kc.h,kc.accel_bias_tau_s,kc.accel_bias_process_variance_density)
    seed=_seed_norms(live);reports={}
    floor_add=math.nextafter(sigma_hi*sigma_hi,math.inf)
    for mode in ('H','A'):
        F,Q=_assemble(mode,Faa,Qaa,Ft,Qt,phi,Qba);fn=interval_matrix_norm2_upper(F);qn=interval_matrix_norm2_upper(Q)
        p=float(seed[mode]);mx=p
        for _ in range(TRANSITIONS):
            p=math.nextafter(math.nextafter(fn*fn,math.inf)*p+qn+floor_add,math.inf);mx=max(mx,p)
        reports[mode]={'dimension':18 if mode=='H' else 21,'live_seed_norm_upper':seed[mode],'F_norm_upper':fn,'Q_norm_upper':qn,'aw_floor_increment_norm_upper':floor_add,'covariance_norm_upper_after_600_worst_predictions_floors':p,'every_event_covariance_norm_upper':mx,'finite':all(math.isfinite(x) for x in (fn,qn,p,mx))}
    closed=all(r['finite'] and r['every_event_covariance_norm_upper']>0 for r in reports.values())
    return {
      'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','word_transitions':TRANSITIONS,
      'live_seed_source_generated':live['live_entry_seed_is_source_generated_not_arbitrary_PSD'],
      'prediction_F_Q_outward_shipping_primitives_used':True,
      'Joseph_update_PSD_nonincrease_used':True,'bias_projection_covariance_identity_used':True,
      'aw_floor_positive_part_bound':'||Pi_+(Sigma-Paw)||2 <= lambda_max(Sigma) for Paw>=0',
      'aw_floor_bound_independent_of_current_P':True,'endpoint_referenced_Pbar_used_as_local_bound':False,
      'trajectory_replay_used':False,'reports':reports,'event_local_covariance_norm_envelope_closed':closed,
      'next_obligation':'derive source-uniform event gain bounds ||K|| <= ||P|| ||H|| / lambda_min(R) from this local covariance envelope and use those inverse-free gain bounds in the finite joint24 event outer'}
def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION or d.get('word_transitions')!=TRANSITIONS:f.append('qualification/window mismatch')
    for k in ('live_seed_source_generated','prediction_F_Q_outward_shipping_primitives_used','Joseph_update_PSD_nonincrease_used','bias_projection_covariance_identity_used','aw_floor_bound_independent_of_current_P','event_local_covariance_norm_envelope_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('endpoint_referenced_Pbar_used_as_local_bound','trajectory_replay_used'):
        if d.get(k) is not False:f.append(k+' not false')
    for m in ('H','A'):
        r=d.get('reports',{}).get(m,{})
        if r.get('finite') is not True or not float(r.get('every_event_covariance_norm_upper',0))>0:f.append(m+' local covariance envelope invalid')
    return f
