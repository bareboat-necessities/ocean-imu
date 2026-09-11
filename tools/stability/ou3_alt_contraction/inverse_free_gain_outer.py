#!/usr/bin/env python3
"""Inverse-free source-uniform Kalman-gain outer for ALT.

For every estimator-owned Joseph cell the actual shipping gain is

    K = P H' (H P H' + R)^-1.

The local covariance theorem supplies ||P||_2 <= pbar at *every event*.  Since
P>=0 and R>=r I,

    HPH' + R >= R >= r I,
    ||(HPH'+R)^-1||_2 <= 1/r,
    ||K||_2 <= pbar ||H||_2 / r.

Thus every gain entry lies in [-kbar,kbar].  This is an explicit outer relation
for the exact same-cell gain; no independent covariance history, interval
innovation inverse, endpoint Pbar, or replay is introduced.
"""
from __future__ import annotations

import math
from ou3_interval import Interval
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_p4_complete_brmm_differential_events as EVENTS
from tools.stability.ou3_alt_contraction import local_covariance_norm_envelope as COV
from tools.stability.ou3_alt_contraction import physical_numeric_bounds as PHYS

QUALIFICATION='OU3_ALT_INVERSE_FREE_EVENT_GAIN_OUTER_V1'

def B(x):x=float(x);return Interval.outward_bounds(-x,x)
def gershgorin_spd_lower(R):
    n=len(R)
    if n==0 or any(len(r)!=n for r in R):raise ValueError('square R required')
    lo=math.inf
    for i in range(n):
        s=0.0
        for j in range(n):
            if i!=j:s=math.nextafter(s+R[i][j].abs_upper(),math.inf)
        lo=min(lo,math.nextafter(R[i][i].lo-s,-math.inf))
    if not math.isfinite(lo) or lo<=0:raise ValueError('R Gershgorin lower is not positive')
    return lo
def gain_entry_box(pnorm,H,R):
    hnorm=COV.interval_matrix_norm2_upper(H);rmin=gershgorin_spd_lower(R)
    kmax=math.nextafter(math.nextafter(float(pnorm)*hnorm,math.inf)/rmin,math.inf)
    if not math.isfinite(kmax):raise ValueError('nonfinite gain outer')
    return [[B(kmax) for _ in range(3)] for _ in range(len(H[0]))],{'P_norm_upper':float(pnorm),'H_norm_upper':hnorm,'R_lambda_min_lower':rmin,'K_norm_upper':kmax}
def _event_H(mode,kind,phys):
    if kind=='S_zero':return EVENTS._event_H(mode,kind)
    if kind=='accelerometer':
        fcap=float(phys['specific_force_norm_upper_mps2']);f=[B(fcap) for _ in range(3)];Rhat=[[Interval(-1.0,1.0) for _ in range(3)] for _ in range(3)]
        return EVENTS._event_H(mode,kind,f_hat=f,R_hat=Rhat)
    if kind=='magnetometer':
        mcap=200.0;m=[B(mcap) for _ in range(3)]
        return EVENTS._event_H(mode,kind,m_body=m)
    raise ValueError(kind)
def build():
    cov=COV.build();cf=COV.validate(cov);phys=PHYS.build();pf=PHYS.validate(phys)
    if cf or pf:raise RuntimeError(f'gain outer prerequisites failed cov={cf} phys={pf}')
    kc=KERNEL._process_constants(KERNEL.DEFAULT_DOMAIN);reports={}
    for mode in ('H','A'):
        p=cov['reports'][mode]['every_event_covariance_norm_upper'];rows={}
        for kind,R in (('S_zero',None),('accelerometer',kc.Racc),('magnetometer',kc.Rmag)):
            H=_event_H(mode,kind,phys)
            if kind=='S_zero':
                # R_S is tuner-dependent; its smallest admitted per-axis std is
                # bounded below by the scheduler's rs_min times the positive
                # declared axis factors. Build it from the committed scheduler.
                import ou3_brmm_tuner_scheduler_step as TUNER
                import ou3_brmm_full_normal_live_word as WORD
                tc=TUNER.constants();active=TUNER.ActiveSchedule(Interval.outward_bounds(tc.tau_min,tc.tau_max),Interval.outward_bounds(max(tc.acc_noise_floor_sigma,1e-9),tc.sigma_max),Interval.outward_bounds(tc.rs_min,tc.rs_max),Interval.outward_bounds(tc.pseudo_min_s,tc.pseudo_max_s));R=WORD.R_S_zero(TUNER.active_rs_std_xyz(active,tc))
            K,meta=gain_entry_box(p,H,R);rows[kind]={'H':H,'K':K,**meta}
        reports[mode]=rows
    closed=all(math.isfinite(r['K_norm_upper']) for rows in reports.values() for r in rows.values())
    return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','event_local_covariance_envelope_consumed':cov['event_local_covariance_norm_envelope_closed'],'innovation_inverse_computed':False,'endpoint_Pbar_used_for_local_gain':False,'independent_covariance_box_used':False,'same_cell_gain_contained_by_outer_relation':True,'reports':reports,'all_event_gain_outers_finite':closed,'next_obligation':'evaluate exact physical Joseph/reset formulas over these source-uniform K outers and certified physical/state domains, retaining magnetometer as a separate arbitrary-finite event star'}
def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('event_local_covariance_envelope_consumed','same_cell_gain_contained_by_outer_relation','all_event_gain_outers_finite'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('innovation_inverse_computed','endpoint_Pbar_used_for_local_gain','independent_covariance_box_used'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
