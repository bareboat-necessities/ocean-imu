"""Compose persisted tuner candidate into the next shipping boundary commit.

Input is the exact ``finite_tuner_frontend_prefix.State`` carried from sample k.
If its pending bit is false, shipping's boundary service is an identity and no
noise-floor/sqrt witness is consumed.  If pending is true, the commit consumes
that SAME persisted TuneState and derives ``band_noise_floor_sigma()`` from the
carried adaptive-band state: raw bench sigma before band readiness, otherwise
bench_sigma*sqrt(p11).  The resulting CommitResult is immediately converted to
``ActiveParameters`` used by prediction, S scheduling and Live R_S.

This closes candidate -> pending -> next-boundary active-parameter ancestry in
exact real arithmetic.  Binary32 sqrt/commit rounding remains open.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_tuner_frontend_prefix as PREFIX
from tools.stability.ou3_alt_contraction import finite_tuner_commit as COMMIT
from tools.stability.ou3_alt_contraction import finite_runtime_parameters as ACTIVE
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as BAND
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P


def R(x): return P.rational(x)


@dataclass(frozen=True)
class Result:
    state: PREFIX.State
    commit: COMMIT.CommitResult | None
    active: ACTIVE.ActiveParameters | None
    band_noise_floor_sigma: F | None


def _noise_floor(state:PREFIX.State,*,bench_noise_sigma,noise_sqrt:BAND.NoiseSqrtWitness|None):
    bench=R(bench_noise_sigma)
    if bench<0: raise ValueError('bench noise sigma must be nonnegative')
    if not state.band.ready:
        if noise_sqrt is not None: raise ValueError('unready band boundary consumes no sqrt(gain) witness')
        return bench
    if noise_sqrt is None: raise ValueError('ready band boundary requires sqrt(gain) witness')
    gain=max(F(0),state.band.p11)
    if noise_sqrt.sqrt_gain*noise_sqrt.sqrt_gain != gain:
        raise ValueError('boundary noise sqrt detached from carried adaptive-band covariance')
    return bench*noise_sqrt.sqrt_gain


def apply(state:PREFIX.State,cfg:COMMIT.CommitConfig,*,live,bench_noise_sigma,
          noise_sqrt:BAND.NoiseSqrtWitness|None=None,rs_sqrt_scale=None,
          sync_covariance=False,rs_scale=1):
    if not isinstance(state,PREFIX.State) or not isinstance(cfg,COMMIT.CommitConfig):
        raise TypeError('persisted tuner-prefix state and commit config required')
    if not isinstance(live,bool) or not isinstance(sync_covariance,bool):
        raise TypeError('literal live/sync branches required')
    if not state.pending:
        if noise_sqrt is not None or rs_sqrt_scale is not None:
            raise ValueError('no-pending boundary consumes no commit witnesses')
        return Result(state,None,None,None)
    nf=_noise_floor(state,bench_noise_sigma=bench_noise_sigma,noise_sqrt=noise_sqrt)
    c=COMMIT.commit(state.tune,cfg,pending=True,live=live,band_noise_floor_sigma=nf,
                    rs_sqrt_scale=rs_sqrt_scale,sync_covariance=sync_covariance,rs_scale=rs_scale)
    if not isinstance(c,COMMIT.CommitResult) or c.pending_after:
        raise AssertionError('pending shipping commit did not clear pending bit')
    active=ACTIVE.ActiveParameters.from_commit(c)
    return Result(replace(state,pending=False),c,active,nf)


def readiness():
    return {
      'persisted_TuneState_consumed_at_next_boundary':True,
      'persisted_pending_bit_consumed_and_cleared':True,
      'band_noise_floor_derived_from_carried_band_state':True,
      'commit_immediately_promoted_to_existing_ActiveParameters':True,
      'no_pending_boundary_is_literal_identity':True,
      'candidate_to_active_parameter_ancestry_closed':True,
      'commit_sqrt_and_binary32_roundoff_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
