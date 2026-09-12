"""Bridge persistent startup frontend state into the first finite Live prefix.

Shipping ``goLive`` does not restart Mahony, WPE, adaptive band/statistics,
stillness or the vibration guard.  It changes the inner startup stage from
TunerReady to Live and unconditionally applies the CURRENT TuneState to OU tau,
stationary Sigma_aw, pseudo-S cadence and Live R_S.  The online pending bit is
not consumed by ``goLive`` and therefore persists to the first Live IMU boundary.

This module composes that control event with the already-derived fresh H18
``CORE.State`` and returns the exact state shape consumed by
``finite_live_imu_prefix.step``.  It removes a frontend restart/gap between
startup and Live.  Commit sqrt/nextafter and deployment clock arithmetic remain
conditional finite-precision witnesses.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_core as CORE
from tools.stability.ou3_alt_contraction import finite_fresh_joint24_entry as FRESH
from tools.stability.ou3_alt_contraction import finite_startup_live_entry as LIVE
from tools.stability.ou3_alt_contraction import finite_guarded_tuner_prefix as FRONT
from tools.stability.ou3_alt_contraction import finite_tuner_commit as COMMIT
from tools.stability.ou3_alt_contraction import finite_runtime_parameters as ACTIVE
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as BAND
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_post_prediction as POST
from tools.stability.ou3_alt_contraction import finite_periodic_aw_sync as AWSYNC
from tools.stability.ou3_alt_contraction import finite_racc_runtime as RACC
from tools.stability.ou3_alt_contraction import finite_live_imu_prefix as LIVEWORD
from tools.stability.ou3_alt_contraction import deployment_scope as SCOPE


def R(x): return P.rational(x)

@dataclass(frozen=True)
class Result:
    state: LIVEWORD.State
    active: ACTIVE.ActiveParameters
    frontend_before: FRONT.State
    frontend_live: FRONT.State
    band_noise_floor_sigma: object


def _noise_floor(tuner_state,*,bench_noise_sigma,noise_sqrt:BAND.NoiseSqrtWitness|None):
    bench=R(bench_noise_sigma)
    if bench<0: raise ValueError('bench noise sigma must be nonnegative')
    if not tuner_state.band.ready:
        if noise_sqrt is not None: raise ValueError('unready startup band consumes no noise sqrt witness')
        return bench
    if noise_sqrt is None: raise ValueError('ready startup band requires sqrt(p11) witness')
    gain=max(F(0),tuner_state.band.p11)
    if noise_sqrt.sqrt_gain*noise_sqrt.sqrt_gain != gain:
        raise ValueError('startup band-noise sqrt detached from carried covariance')
    return bench*noise_sqrt.sqrt_gain


def bridge(entry:LIVE.Result,fresh:CORE.State,frontend:FRONT.State,*,scope:SCOPE.Scope,
           commit_cfg:COMMIT.CommitConfig,bench_noise_sigma,
           noise_sqrt:BAND.NoiseSqrtWitness|None=None,rs_sqrt_scale=None,rs_scale=1,
           scheduler:POST.Scheduler,racc:RACC.State,aw_sync:AWSYNC.State,
           scheduler_park:ACTIVE.NextafterParkWitness|None=None):
    if not isinstance(entry,LIVE.Result) or not isinstance(fresh,CORE.State):
        raise TypeError('startup Live entry and fresh CORE.State required')
    if not isinstance(frontend,FRONT.State) or not isinstance(frontend.tuner.stage,str):
        raise TypeError('persistent guarded frontend state required')
    if frontend.tuner.stage!='TunerReady':
        raise ValueError('goLive frontend bridge starts from TunerReady')
    if fresh.mode!='H' or fresh.covariance!=entry.P or fresh.q_hat!=entry.q_hat:
        raise ValueError('fresh H18 CORE state detached from startup Live entry')
    if fresh.reference.time != frontend.tuner.time:
        raise ValueError('physical/filter and frontend clocks detached at goLive')
    SCOPE.assert_certified_scope(scope)
    if not isinstance(scheduler,POST.Scheduler) or not isinstance(racc,RACC.State) or not isinstance(aw_sync,AWSYNC.State):
        raise TypeError('scheduler, Racc and aw-sync persistent states required')
    if aw_sync.pending:
        # Periodic synchronization is Live-only, so no queued request can be
        # inherited from a correctly composed startup history.
        raise ValueError('startup cannot enter Live with a queued periodic aw floor')

    nf=_noise_floor(frontend.tuner,bench_noise_sigma=bench_noise_sigma,noise_sqrt=noise_sqrt)
    c=COMMIT.commit(frontend.tuner.tune,commit_cfg,pending=True,live=True,
                    band_noise_floor_sigma=nf,rs_sqrt_scale=rs_sqrt_scale,
                    sync_covariance=True,rs_scale=rs_scale)
    active=ACTIVE.ActiveParameters.from_commit(c)
    if active != entry.active:
        raise ValueError('fresh Live active parameters detached from carried TunerReady TuneState')

    # enterLive_ changes only startup stage/clock in the carried frontend state.
    # The online pending bit intentionally survives and may be consumed by the
    # first Live sample's ordinary boundary service.
    tuner_live=replace(frontend.tuner,stage='Live',stage_time=F(0))
    front_live=FRONT.State(frontend.guard,tuner_live)
    sched=ACTIVE.retarget_scheduler(active,scheduler,park=scheduler_park)
    state=LIVEWORD.State(fresh,front_live.guard,front_live.tuner,racc,active,sched,aw_sync)
    return Result(state,active,frontend,front_live,nf)


def readiness():
    return {
      'TunerReady_to_Live_frontend_stage_edge_materialized':True,
      'Mahony_WPE_band_stats_stillness_guard_memory_preserved_across_goLive':True,
      'goLive_unconditional_active_parameters_derived_from_same_carried_TuneState':True,
      'goLive_live_RS_derived_from_same_carried_TuneState':True,
      'goLive_period_retargets_persistent_S_scheduler':True,
      'online_pending_bit_preserved_across_goLive':True,
      'startup_periodic_aw_floor_pending_forbidden':True,
      'fresh_H18_CORE_state_connected_to_first_Live_prefix_shape':True,
      'commit_sqrt_nextafter_binary32_closed':False,
      'wrapper_clock_binary32_closed':False,
      'first_Live_sample_executed_from_fresh_entry':False,
      'complete_same_history_startup_to_Live_word':False,
      'ALT_STARTUP_PASS':False,
      'ALT_LIVE_PASS':False,
    }
