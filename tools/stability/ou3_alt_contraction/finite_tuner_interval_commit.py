"""Exact-real staged commit for an interval SpectralMSE tuner state.

The shipping default law is SpectralMSE, so ``pseudo_update_information_rate_scale_``
returns one: the target already contains realized T_S.  At a pending Live
boundary the commit therefore maps the interval R_S state through the literal
base clamp, optional rs_scale cap, anisotropic XYZ factors, and covariance
squares monotonically.  No point R_S is selected.

Tau/sigma remain point-valued exact-real states here.  This closes the
mathematical interval commit topology only; binary32 cadence, clamp/multiply,
and final covariance storage correspondence remain deployment supplies.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_tuner_candidate_interval as I
from tools.stability.ou3_alt_contraction import finite_tuner_commit as C

QUALIFICATION='OU3_ALT_INTERVAL_TUNER_COMMIT_V1'


@dataclass(frozen=True)
class IntervalCommitResult:
    tau:F
    Sigma_aw:tuple
    pseudo_period:F
    RS_cov_lo:tuple|None
    RS_cov_hi:tuple|None
    aw_floor_target:tuple|None
    pending_after:bool
    qualification:str=QUALIFICATION
    def __post_init__(self):
        object.__setattr__(self,'tau',P.rational(self.tau)); object.__setattr__(self,'pseudo_period',P.rational(self.pseudo_period))
        if self.qualification!=QUALIFICATION: raise ValueError('wrong interval commit qualification')
        if self.pending_after is not False: raise ValueError('successful pending commit must clear pending bit')


def _point(t:I.IntervalTuneState,rs):
    return C.TuneState(t.tau_applied,t.sigma_applied,rs)


def commit(tune:I.IntervalTuneState,cfg:C.CommitConfig,*,pending,live,band_noise_floor_sigma,
           sync_covariance=False,rs_scale=1):
    if not isinstance(tune,I.IntervalTuneState) or not isinstance(cfg,C.CommitConfig):
        raise TypeError('interval tuner state and CommitConfig required')
    if not isinstance(pending,bool) or not isinstance(live,bool) or not isinstance(sync_covariance,bool):
        raise TypeError('literal pending/live/sync branches required')
    if not pending: return None
    if cfg.cubic_rs_law:
        raise ValueError('interval commit theorem currently covers deployed non-Cubic SpectralMSE law only')
    lo=_point(tune,tune.RS_lo); hi=_point(tune,tune.RS_hi)
    period_lo=C.cadence(lo,cfg); period_hi=C.cadence(hi,cfg)
    if period_lo!=period_hi: raise AssertionError('R_S interval cannot change tau-owned pseudo cadence')
    Sigma_lo=C.stationary_aw(lo,cfg,band_noise_floor_sigma=band_noise_floor_sigma)
    Sigma_hi=C.stationary_aw(hi,cfg,band_noise_floor_sigma=band_noise_floor_sigma)
    if Sigma_lo!=Sigma_hi: raise AssertionError('R_S interval cannot change sigma-owned stationary covariance')
    rlo=C.applied_RS(lo,cfg,period_lo,rs_scale=rs_scale) if live else None
    rhi=C.applied_RS(hi,cfg,period_hi,rs_scale=rs_scale) if live else None
    # applied_RS is monotone for nonnegative R_S/factors.  Entry IntervalTuneState
    # comes from a nonnegative SpectralMSE target/EMA in the theorem path.
    if live and any(rlo[i][i]>rhi[i][i] for i in range(3)):
        raise ValueError('interval R_S commit lost monotonicity')
    floor=Sigma_lo if sync_covariance else None
    return IntervalCommitResult(tune.tau_applied,Sigma_lo,period_lo,rlo,rhi,floor,False)


def readiness():
    return {
      'pending_SpectralMSE_interval_commit_boundary_materialized':True,
      'same_tau_drives_exact_real_OU_and_pseudo_cadence':True,
      'same_sigma_drives_exact_real_stationary_aw_target':True,
      'RS_interval_propagated_through_base_clamp_axis_factors_and_covariance_square':True,
      'no_point_RS_representative_selected':True,
      'cubic_information_rate_sqrt_not_in_deployed_SpectralMSE_commit':True,
      'binary32_commit_arithmetic_correspondence_closed':False,
      'machine_real_RS_residual_attached_to_commit':False,
      'goLive_interval_RS_to_actual_MEKF_commit_closed':False,
      'source_uniform_complete_startup_reachability_closed':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,
      'ALT_LIVE_PASS':False,
    }
