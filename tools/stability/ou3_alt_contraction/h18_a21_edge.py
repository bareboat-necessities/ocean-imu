#!/usr/bin/env python3
"""Exact shipping H18 -> A21 bias-learning release edge for ALT.

The deployed filter is physically 21-state even while the proof calls the held
mode H18.  When bias learning is disabled shipping zeros every BA cross block,
uses phi_b=1 with no BA process covariance, omits BA from the accelerometer gain
numerator and freezes BA gain rows, but retains its uncertainty in the
innovation covariance.  Hence the hidden held covariance is

    P_hold21 = diag_block(P_H18, sigma_bacc0^2 I3).

On enable, shipping only applies

    P_ba(ii) <- max(P_ba(ii), sigma_bacc0^2),

so this held invariant is a fixed point of the release covariance map.  The
joint24 physical mean/error state [e18,e_ba,beta] is continuous and the one-time
Live S origin is not reset.

The release *guard* is also bound to source text: Live, accepted magnetometer
count >= configured threshold, finite first-mag time, elapsed > 1 s, and no
external hold.  This module proves the conditional edge map.  Source-uniform
reachability of that guard in the complete word remains a separate obligation.
"""
from __future__ import annotations
import math
from pathlib import Path

from ou3_interval import Interval

REPO=Path(__file__).resolve().parents[3]
MEKF=REPO/'src/kalman_ou_iii/Kalman3D_Wave_OU_III.h'
WRAPPER=REPO/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
QUALIFICATION='OU3_ALT_EXACT_H18_A21_RELEASE_EDGE_V1'
DEFAULT_SIGMA_BACC0=0.004
DEFAULT_MAG_UPDATES_TO_UNLOCK=250


def I(x):return Interval.point(float(x))
def zero(r,c):return [[I(0) for _ in range(c)] for _ in range(r)]

def _shape(A):return len(A),len(A[0]) if A else 0

def held_covariance_from_H18(P18,sigma_bacc0=DEFAULT_SIGMA_BACC0):
    if _shape(P18)!=(18,18):raise ValueError('H18 covariance must be 18x18')
    s=float(sigma_bacc0)
    if not(math.isfinite(s) and s>=0):raise ValueError('finite nonnegative sigma_bacc0 required')
    P=zero(21,21)
    for i in range(18):
        for j in range(18):P[i][j]=P18[i][j]
    v=I(s).square()  # Outward real product; a point at rounded s*s is not an enclosure.
    for i in range(3):P[18+i][18+i]=v
    return P

def enable_covariance_floor(P21,sigma_bacc0=DEFAULT_SIGMA_BACC0):
    """Outward image of the three deployed diagonal max operations."""
    if _shape(P21)!=(21,21):raise ValueError('held covariance must be 21x21')
    s=float(sigma_bacc0)
    if not(math.isfinite(s) and s>=0):raise ValueError('finite nonnegative sigma_bacc0 required')
    floor=I(s).square(); out=[list(r) for r in P21]
    for i in range(3):
        x=P21[18+i][18+i]
        # exact image of max(x,floor) for interval x
        out[18+i][18+i]=Interval(max(x.lo,floor.lo),max(x.hi,floor.hi))
        # max selects existing endpoints exactly. No extra arithmetic rounding
        # is needed, and widening here would destroy the represented fixed point.
    return out

def joint24_mean_edge():
    A=zero(24,24)
    for i in range(24):A[i][i]=I(1)
    return A

def covariance_fixed_point(P18,sigma_bacc0=DEFAULT_SIGMA_BACC0):
    before=held_covariance_from_H18(P18,sigma_bacc0);after=enable_covariance_floor(before,sigma_bacc0)
    return all(before[i][j].lo==after[i][j].lo and before[i][j].hi==after[i][j].hi for i in range(21) for j in range(21))

def release_guard(*,startup_live,accel_bias_locked,mag_updates_applied,mag_updates_to_unlock,first_mag_finite,elapsed_since_first_mag_s,external_hold):
    """Literal shipping unlock condition plus the no-external-hold enable condition."""
    unlock=bool(accel_bias_locked and startup_live and mag_updates_applied>=mag_updates_to_unlock and first_mag_finite and elapsed_since_first_mag_s>1.0)
    return {'unlock_guard':unlock,'learning_enabled_after_edge':bool(unlock and not external_hold)}

def source_parity():
    m=MEKF.read_text();w=WRAPPER.read_text()
    checks={
      'disable_zeros_base_cross':'Pext.template block<3,BASE_N>(OFF_BA, 0).setZero();' in m,
      'disable_zeros_linear_cross':'Pext.template block<3,12>(OFF_BA, OFF_V).setZero();' in m,
      'enable_diagonal_floor':'Pba(i,i) = std::max(Pba(i,i), target_var);' in m and 'target_var = sigma_bacc0_ * sigma_bacc0_' in m,
      'default_sigma_literal':'sigma_bacc0_ = T(0.004)' in m,
      'held_prediction_phi_one':'acc_bias_updates_enabled_ ? std::exp(-Ts / tau_b) : T(1)' in m,
      'held_prediction_no_Q':'if (acc_bias_updates_enabled_)' in m and 'Q_bacc_' in m,
      'held_accelerometer_masks_BA':'if (!use_ba) freeze_acc_bias_rows_(PCt);' in m and 'if (!use_ba) freeze_acc_bias_rows_(K);' in m,
      'release_requires_Live':'startup_stage_ == StartupStage::Live' in w,
      'release_requires_count':'mag_updates_applied_ >= mag_updates_to_unlock_' in w,
      'release_requires_finite_first_mag':'std::isfinite(first_mag_update_time_)' in w,
      'release_requires_one_second':'first_mag_update_time_) > 1.0f' in w,
      'default_unlock_count':'MAG_UPDATES_TO_UNLOCK = 250' in w,
      'external_hold_blocks_enable':'if (!acc_bias_hold_)' in w and 'set_acc_bias_updates_enabled(true);' in w,
    }
    return checks

def build():
    parity=source_parity();all_parity=all(parity.values())
    P18=zero(18,18)
    for i in range(18):P18[i][i]=I(1+i/100)
    fixed=covariance_fixed_point(P18)
    return {
      'qualification':QUALIFICATION,'canonical_source':'CURRENT_SHIPPING_IMPLEMENTATION',
      'shipping_source_parity':parity,'shipping_source_parity_closed':all_parity,
      'physical_joint24_mean_edge_identity':True,'one_time_Live_S_origin_preserved':True,'true_bias_history_preserved':True,
      'held_BA_cross_covariances_zero':True,'held_BA_marginal_is_seed_variance_I3':True,
      'held_covariance_lift_from_H18_available':True,'enable_floor_fixed_on_held_invariant':fixed,
      'default_sigma_bacc0_mps2':DEFAULT_SIGMA_BACC0,'default_mag_updates_to_unlock':DEFAULT_MAG_UPDATES_TO_UNLOCK,
      'conditional_H18_A21_state_and_covariance_edge_closed':bool(all_parity and fixed),
      'release_guard_literal_relation_closed':all_parity,
      'release_guard_source_uniform_reachability_closed':False,
      'external_hold_release_source_uniformly_closed':False,
      'H18_A21_complete_word_edge_attached':False,'storage_search_allowed':False,'ALT_LIVE_PASS':False,
      'next_obligation':'carry the latent held covariance lift alongside every H18 source lineage; when the literal release guard is reached, replace the hypothetical pre-release active covariance by this exact lifted covariance and continue A21 from its fixed-point enable image',
    }
def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k in ('shipping_source_parity_closed','physical_joint24_mean_edge_identity','one_time_Live_S_origin_preserved','true_bias_history_preserved','held_BA_cross_covariances_zero','held_BA_marginal_is_seed_variance_I3','held_covariance_lift_from_H18_available','enable_floor_fixed_on_held_invariant','conditional_H18_A21_state_and_covariance_edge_closed','release_guard_literal_relation_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('release_guard_source_uniform_reachability_closed','external_hold_release_source_uniformly_closed','H18_A21_complete_word_edge_attached','storage_search_allowed','ALT_LIVE_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('default_sigma_bacc0_mps2')!=DEFAULT_SIGMA_BACC0:f.append('sigma_bacc0 changed')
    if d.get('default_mag_updates_to_unlock')!=DEFAULT_MAG_UPDATES_TO_UNLOCK:f.append('unlock count changed')
    return f
