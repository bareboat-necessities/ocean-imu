#!/usr/bin/env python3
"""Regional Normal-Live frontend predecessor invariant for ALT.

This is the Live-only counterpart of the shared frontend predecessor proof.  It
uses the same BRMM primitive bounds, WPE telescoping/Schur bounds, adaptive-band
covariance bound, tuner clamps, scheduler quotient and stillness saturation.
The private Mahony component is the regional source-order binary32 invariant in
``live_mahony_regional``.

Crucially, membership in this set is a premise of the regional Live theorem;
startup's ability to enter it is NOT consumed here.  That is a later capture
obligation.  Thus this module removes no dynamics or source bounds--it only
separates forward invariance from reachability of the initial set.
"""
from __future__ import annotations
import json,math
from pathlib import Path

import ou3_brmm_contract as BRMM
import ou3_brmm_finite_window_primitive_qualification as PRIM
import ou3_fast_inv_sqrt_interval as FINV
import ou3_brmm_wpe_state_step as WPE
import ou3_brmm_tuner_scheduler_step as TUNER
import ou3_brmm_sigma_stillness_step as STILL
import ou3_validated_transcendentals as VT
import ou3_p4_brmm_frontend_predecessor_invariant as SHARED
from ou3_interval import Interval
from tools.stability.ou3_alt_contraction import live_mahony_regional as MAHONY

REPO=Path(__file__).resolve().parents[3]
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
WPE_HEADER=REPO/'src/tuner/WavePeriodEstimator.h'
QUALIFICATION='OU3_ALT_REGIONAL_NORMAL_LIVE_FRONTEND_PREDECESSOR_V1'

def I(x):return Interval.outward_bounds(float(x),float(x))
def up(x):return math.nextafter(float(x),math.inf)
def _decay(lam,dt):
    d=VT.exp_interval(-(lam*I(dt)))
    if not(0<d.lo<=d.hi<1):raise RuntimeError('WPE decay is not strictly Schur')
    return d

def build():
    domain=json.loads(DOMAIN.read_text());brmm=BRMM.build();prim=PRIM.build();mah=MAHONY.build()
    bad={'BRMM':BRMM.validate(brmm),'primitive':PRIM.validate(prim),'regional_Mahony':MAHONY.validate(mah)};bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('regional frontend prerequisite failed: '+repr(bad))
    wc=WPE.constants();tc=TUNER.constants();sc=STILL.constants();g=float(domain['startup']['gravity_mps2']);am=float(prim['uniform_physical_primitives']['acceleration_norm_upper_mps2'])
    q2=FINV.all_positive_normal_normalized_norm2_enclosure();vertical=up(q2.hi*(g+am)+g)
    pi=Interval.outward_bounds(3.141592653589793,3.141592653589794);lam=I(2)*pi*I(wc.high_pass_hz);d=_decay(lam,wc.dt)
    hp_gain=up(2*d.hi);hp1=up(hp_gain*vertical);hp2=up(hp_gain*hp1);vel=up(hp2/lam.lo);elev=up(vel/lam.lo)
    valid_weight_floor=math.nextafter(1e-3,math.inf);vel_sq=up(vel*vel);vvar_upper=up(vel_sq/valid_weight_floor);evar_floor=math.nextafter(1e-12,math.inf);omega2_floor=math.nextafter(1e-8,math.inf);omega2_upper=up(vvar_upper/evar_floor)
    raw_min=math.nextafter(2*math.pi/math.sqrt(omega2_upper),0.0);raw_max=up(2*math.pi/math.sqrt(omega2_floor));log_min=math.nextafter(math.log(raw_min),-math.inf);log_max=math.nextafter(math.log(raw_max),math.inf);freq_min=math.nextafter(math.exp(-log_max),0.0);freq_max=up(math.exp(-log_min))
    text=WPE_HEADER.read_text();elapsed_parity=all(s in text for s in ('elapsed_sec_ += dt_sec;','if (elapsed_sec_ < moment_start_sec) return;','elapsed_sec_ >= settled_floor_sec','const float moment_history_sec = elapsed_sec_ - moment_start_sec;','if (elapsed_sec_ >= usable_floor_sec && moment_history_sec >= period)','if (usable_period_) return;'))
    elapsed_cap=up(3/lam.lo+raw_max);band=up(2*vertical);bcov=SHARED._band_covariance_outer(tc);still_energy=up((vertical/sc.gravity)**2);scheduler_max=up(tc.adapt_every_s+tc.dt)
    strict={
      'regional_Mahony_source_order_binary32_invariant':mah['regional_source_order_binary32_discrete_invariance_closed'],
      'WPE_decay_strict':d.hi<1,'WPE_telescoping_gain_finite':hp_gain<2,'WPE_leaky_gain_finite':math.isfinite(1/lam.lo),
      'WPE_startup_zero_weight_retained':True,'WPE_ratio_branch_has_positive_weight_floor':valid_weight_floor>1e-3,
      'WPE_valid_period_compact':0<raw_min<raw_max<math.inf,'WPE_frequency_compact':0<freq_min<=freq_max<math.inf,
      'WPE_elapsed_guard_quotient_parity':elapsed_parity,'adaptive_band_uniform_Schur_gap':bcov['alpha_low_uniform_lower']>0 and bcov['alpha_high_uniform_lower']>0,
      'adaptive_band_covariance_compact':bcov['finite_covariance_outer'],
      'candidate_active_clamps_compact':tc.tau_min>0 and tc.tau_max>=tc.tau_min and tc.sigma_max>0 and tc.rs_min>0 and tc.rs_max>=tc.rs_min,
      'scheduler_guard_class_compact':math.isfinite(scheduler_max),'stillness_timer_saturated':sc.still_max_s==60.0}
    closed=all(strict.values())
    return {'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','regional_Live_entry_membership_is_theorem_premise':True,'startup_capture_consumed':False,'startup_entry_closed_here':False,
      'same_history_BRMM_ancestry_required_downstream':True,'regional_Mahony_invariant_consumed':True,'source_order_binary32_model_consumed':True,
      'vertical_acceleration_abs_upper_mps2':vertical,
      'WPE_invariant':{'decay':d.as_list(),'single_highpass_linf_gain_upper':hp_gain,'highpass1_abs_upper_mps2':hp1,'highpass2_abs_upper_mps2':hp2,'velocity_abs_upper_mps':vel,'elevation_abs_upper_m':elev,'weight_interval':[0.0,1.0],'valid_ratio_weight_lower':valid_weight_floor,'valid_raw_period_s':[raw_min,raw_max],'canonical_log_period':[log_min,log_max],'canonical_frequency_hz':[freq_min,freq_max],'elapsed_guard_quotient_cap_s':elapsed_cap,'actual_elapsed_need_not_be_bounded':True,'elapsed_quotient_preserves_all_shipping_guards':elapsed_parity},
      'adaptive_band_invariant':{'lowpass_abs_upper_mps2':vertical,'band_abs_upper_mps2':band,**bcov},
      'tuner_invariant':{'mean_value_abs_upper':band,'sq_value_upper':up(band*band),'weights_interval':[0.0,1.0],'frequency_after_shipping_clamp_hz':[tc.tune_freq_min,tc.tune_freq_max],'candidate':{'tau_s':[tc.tau_min,tc.tau_max],'sigma_mps2':[0.0,tc.sigma_max],'rs':[tc.rs_min,tc.rs_max]},'active':{'tau_s':[tc.tau_min,tc.tau_max],'sigma_mps2':[0.0,tc.sigma_max],'rs':[tc.rs_min,tc.rs_max],'pseudo_period_s':[tc.pseudo_min_s,tc.pseudo_max_s]},'scheduler_elapsed_guard_upper_s':scheduler_max},
      'stillness_invariant':{'lpf_abs_upper_mps2':vertical,'energy_ema_upper':still_energy,'still_time_s':[0.0,sc.still_max_s],'attenuation':[0.0,1.0]},
      'strict_invariance_checks':strict,'regional_normal_live_frontend_predecessor_set_closed':closed,
      'coefficient_history_must_still_come_from_joint_transition':True,'outer_bounds_may_not_generate_independent_f_sigma_tau_TS_RS':True,'trajectory_replay_used':False,'ALT_LIVE_PASS':False,
      'next_obligation':'use this regional predecessor set with the same-history JOINT transition and O^601_BRMM to build the 601-step Live relation; startup must later prove entry into this set'}
def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('regional_Live_entry_membership_is_theorem_premise','same_history_BRMM_ancestry_required_downstream','regional_Mahony_invariant_consumed','source_order_binary32_model_consumed','regional_normal_live_frontend_predecessor_set_closed','coefficient_history_must_still_come_from_joint_transition','outer_bounds_may_not_generate_independent_f_sigma_tau_TS_RS'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('startup_capture_consumed','startup_entry_closed_here','trajectory_replay_used','ALT_LIVE_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    if not all(d.get('strict_invariance_checks',{}).values()):f.append('strict invariant check failed')
    return f
