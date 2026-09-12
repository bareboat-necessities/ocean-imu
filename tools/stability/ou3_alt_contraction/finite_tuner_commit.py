"""Finite staged tuner commit for the ALT runtime graph.

The online tuner smooths ``TuneState`` on sample k and only sets a pending bit.
At the beginning of sample k+1, ``apply_pending_online_tune_`` commits that one
state to OU tau/stationary covariance, pseudo-S cadence, and (when Live) R_S.
This module materializes that same-boundary transaction in exact real arithmetic.
It does not yet derive TuneState from WPE/bandpass/sigma history or enclose
binary32 sqrt/clamp rounding.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P


@dataclass(frozen=True)
class TuneState:
    tau_applied: F
    sigma_applied: F
    RS_applied: F
    def __post_init__(self):
        tau,sigma,rs=map(P.rational,(self.tau_applied,self.sigma_applied,self.RS_applied))
        if tau <= 0 or sigma < 0: raise ValueError('valid applied tau/sigma required')
        object.__setattr__(self,'tau_applied',tau); object.__setattr__(self,'sigma_applied',sigma); object.__setattr__(self,'RS_applied',rs)


@dataclass(frozen=True)
class CommitConfig:
    pseudo_tau_ratio: F
    pseudo_period_min: F
    pseudo_period_max: F
    pseudo_fixed_period: F
    tau_scaled_cadence: bool
    cubic_rs_law: bool
    min_R_S: F
    max_R_S: F
    S_factor: F
    R_S_x_factor: F
    R_S_y_factor: F
    def __post_init__(self):
        vals=[P.rational(x) for x in (self.pseudo_tau_ratio,self.pseudo_period_min,self.pseudo_period_max,self.pseudo_fixed_period,self.min_R_S,self.max_R_S,self.S_factor,self.R_S_x_factor,self.R_S_y_factor)]
        if not isinstance(self.tau_scaled_cadence,bool) or not isinstance(self.cubic_rs_law,bool): raise TypeError('literal cadence/law branches required')
        ratio,pmin,pmax,pfixed,rmin,rmax,sf,xf,yf=vals
        if ratio<=0 or pmin<=0 or pmax<pmin or pfixed<=0 or rmax<rmin or sf<0 or xf<0 or yf<0: raise ValueError('invalid commit configuration')
        for name,val in zip(('pseudo_tau_ratio','pseudo_period_min','pseudo_period_max','pseudo_fixed_period','min_R_S','max_R_S','S_factor','R_S_x_factor','R_S_y_factor'),vals): object.__setattr__(self,name,val)


@dataclass(frozen=True)
class CommitResult:
    tau: F
    Sigma_aw: tuple
    pseudo_period: F
    R_S: tuple | None
    aw_floor_target: tuple | None
    pending_after: bool


def clamp(x,lo,hi): return min(max(x,lo),hi)


def cadence(tune:TuneState,cfg:CommitConfig):
    if not cfg.tau_scaled_cadence: return cfg.pseudo_fixed_period
    return clamp(cfg.pseudo_tau_ratio*tune.tau_applied,cfg.pseudo_period_min,cfg.pseudo_period_max)


def stationary_aw(tune:TuneState,cfg:CommitConfig,*,band_noise_floor_sigma):
    nf=P.rational(band_noise_floor_sigma)
    if nf < 0: raise ValueError('band noise floor sigma must be nonnegative')
    sigma_floor=max(F(1,20),nf)  # shipping max(0.05, band_noise_floor_sigma())
    sZ=max(sigma_floor,tune.sigma_applied); sH=sZ*cfg.S_factor
    return ((sH*sH,0,0),(0,sH*sH,0),(0,0,sZ*sZ))


def information_rate_scale(cfg:CommitConfig,period,*,sqrt_scale=None):
    """Finite witness for shipping sqrt(fixed_period/period) cubic normalization."""
    period=P.rational(period)
    if not cfg.cubic_rs_law or not cfg.tau_scaled_cadence:
        if sqrt_scale is not None and P.rational(sqrt_scale)!=1: raise ValueError('non-cubic/fixed cadence information scale is one')
        return F(1)
    if sqrt_scale is None: raise ValueError('cubic tau-scaled branch requires sqrt witness')
    s=P.rational(sqrt_scale)
    if s<=0 or s*s*period != cfg.pseudo_fixed_period: raise ValueError('sqrt information-rate witness detached from SAME cadence')
    return s


def applied_RS(tune:TuneState,cfg:CommitConfig,period,*,sqrt_scale=None,rs_scale=1):
    scale=P.rational(rs_scale)
    if scale<=0: scale=F(1)
    scale=min(scale,F(1))
    base=clamp(tune.RS_applied,cfg.min_R_S,cfg.max_R_S)
    info=information_rate_scale(cfg,period,sqrt_scale=sqrt_scale)
    z=base*info*scale; sig=(z*cfg.R_S_x_factor,z*cfg.R_S_y_factor,z)
    return ((sig[0]*sig[0],0,0),(0,sig[1]*sig[1],0),(0,0,sig[2]*sig[2]))


def commit(tune:TuneState,cfg:CommitConfig,*,pending,live,band_noise_floor_sigma,
           rs_sqrt_scale=None,sync_covariance=False,rs_scale=1):
    """One literal pending-tuner transaction at an IMU sample boundary."""
    if not isinstance(pending,bool) or not isinstance(live,bool) or not isinstance(sync_covariance,bool): raise TypeError('literal pending/live/sync branches required')
    if not pending: return None
    period=cadence(tune,cfg); Sigma=stationary_aw(tune,cfg,band_noise_floor_sigma=band_noise_floor_sigma)
    RS=applied_RS(tune,cfg,period,sqrt_scale=rs_sqrt_scale,rs_scale=rs_scale) if live else None
    # Periodic online adaptation calls apply_ou_tune_(false), so it changes the
    # stationary process covariance but does not queue/rewrite posterior P_aw.
    floor=Sigma if sync_covariance else None
    return CommitResult(tune.tau_applied,Sigma,period,RS,floor,False)


def readiness():
    return {
      'pending_commit_boundary_materialized':True,
      'same_tau_drives_OU_and_S_cadence':True,
      'same_sigma_drives_stationary_aw_target':True,
      'live_RS_from_same_TuneState_and_realized_period':True,
      'periodic_online_commit_does_not_sync_posterior_aw':True,
      'WPE_band_sigma_candidate_history_attached':False,
      'sqrt_binary32_roundoff_attached':False,
      'applied_tune_float_roundoff_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
