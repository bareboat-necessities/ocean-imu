"""First composed default-policy Live IMU sample prefix for ALT.

This module joins the previously separate same-history relations in the shipping
sample order, while remaining deliberately conditional on the still-open hybrid
and finite-precision branches.

Represented coupled order:
  1. apply a tuner candidate pending from sample k-1; if T_S changes, retarget
     scheduler credit with shipping's progress-preserving rule;
  2. advance one persistent AccelVibrationGuard from the raw predecessor packet;
  3. advance private Mahony from that exact guarded ``acc_in``;
  4. compute pre-measurement Racc from the same guard excess and PRE-UPDATE WPE /
     TuneState schedule;
  5. execute the active-parameter-rooted MEKF prediction, queued a_w floor,
     scheduler/S service, and held guarded accelerometer correction;
  6. execute the measurement-only tracker-free tuner suffix using the already
     computed Mahony successor; current WPE update remains last;
  7. execute the deployed periodic aw-sync request recurrence.  It commutes with
     the measurement-only suffix because it reads only time/current ActiveSigma
     and writes a disjoint aw-sync state; critically, any request snapshots the
     OLD/current active Sigma for the next sample.

Not represented here: a firing Live tilt watchdog reset, async magnetometer
calls, legacy/congruent immediate aw-sync policies, direction sidecar state,
source-uniform COMPLETE-BRMM bounds, or remaining deployment floating-point
branches.  Consequently this is a finite identity prefix, NOT ALT_LIVE_PASS.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_core as CORE
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as GUARD
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as VERT
from tools.stability.ou3_alt_contraction import finite_tuner_frontend_prefix as TUNER
from tools.stability.ou3_alt_contraction import finite_tuner_boundary_commit as BOUND
from tools.stability.ou3_alt_contraction import finite_tuner_commit as COMMIT
from tools.stability.ou3_alt_contraction import finite_runtime_parameters as ACTIVE
from tools.stability.ou3_alt_contraction import finite_active_runtime_word as WORD
from tools.stability.ou3_alt_contraction import finite_post_prediction as POST
from tools.stability.ou3_alt_contraction import finite_measurement_runtime as MEAS
from tools.stability.ou3_alt_contraction import finite_racc_runtime as RACC
from tools.stability.ou3_alt_contraction import finite_periodic_aw_sync as AWSYNC
from tools.stability.ou3_alt_contraction import finite_frontend_runtime as FRONT
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as BAND
from tools.stability.ou3_alt_contraction import finite_physical_prediction as PHYS
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P


def R(x): return P.rational(x)

@dataclass(frozen=True)
class State:
    mekf: CORE.State
    guard: GUARD.State
    tuner: TUNER.State
    racc: RACC.State
    active: ACTIVE.ActiveParameters
    scheduler: POST.Scheduler
    aw_sync: AWSYNC.State
    def __post_init__(self):
        if not isinstance(self.mekf,CORE.State) or not isinstance(self.guard,GUARD.State) or not isinstance(self.tuner,TUNER.State):
            raise TypeError('MEKF, guard and tuner states required')
        if not isinstance(self.racc,RACC.State) or not isinstance(self.active,ACTIVE.ActiveParameters):
            raise TypeError('Racc and active parameter states required')
        if not isinstance(self.scheduler,POST.Scheduler) or not isinstance(self.aw_sync,AWSYNC.State):
            raise TypeError('S scheduler and periodic aw-sync states required')
        if self.tuner.stage != 'Live': raise ValueError('this prefix represents the Live stage only')
        self.active.require_scheduler(self.scheduler)
        if self.active.R_S is None: raise ValueError('Live prefix requires an active R_S')
        if self.tuner.time != self.mekf.reference.time:
            raise ValueError('tuner/filter clocks detached at Live sample entry')

@dataclass(frozen=True)
class Result:
    state: State
    boundary: BOUND.Result
    guarded: SENSOR.GuardedImuSample
    vertical: VERT.Result
    preupdate_wpe: BAND.WPEFrequencyView
    racc: RACC.Result
    prediction: WORD.Predicted
    post_prediction: POST.PostPrediction
    S_service: POST.ServicedPostPrediction
    accelerometer: MEAS.MeasurementRuntimeResult
    tuner_suffix: TUNER.Result
    aw_sync: AWSYNC.Result


def step(state:State,raw:SENSOR.RawImuSample,segment:PHYS.PhysicalSegment,*,dt,
         commit_cfg:COMMIT.CommitConfig,boundary_bench_noise_sigma,
         boundary_noise_sqrt=None,rs_sqrt_scale=None,scheduler_park=None,
         guard_cfg:GUARD.Config,guard_decay=None,guard_rms=None,
         vertical_cfg:VERT.Config,accel_invnorm=None,quat_invnorm=None,seed=None,
         band_cfg:BAND.BandConfig,
         preupdate_period=None,preupdate_frequency=None,
         racc_cfg:RACC.Config,nominal_racc_std,rao_witness=None,racc_sqrt=None,
         angular=None,Qbase=None,ou=None,bias=None,qaxis=None,
         use_exact_attitude_Q=True,attitude_first_ldlt_success=True,
         attitude_second_ldlt_success=None,
         floor_solver_success=None,floor_eigenvectors=None,floor_eigenvalues=None,
         S_ldlt=None,S_alpha=1,S_radius=None,
         accel_conditioning:SENSOR.AccelConditioning=None,accel_ldlt:MEAS.SafeLDLT=None,
         accel_alpha=1,accel_radius=F(2,5),
         tilt_reset_due=False,
         aw_sync_adapt_every=F(1,10),aw_sync_enabled=True,
         wpe_cfg=None,wpe_decay=None,wpe_moment_decay=None,wpe_period_witness=None,
         wpe_log_witness=None,wpe_post_output=None,
         stats_cfg=None,band_decay=None,variance_decay=None,
         bench_noise_sigma=F(0),noise_sqrt=None,tracker_lpf_decay=None,
         still_cfg=None,still_attenuation=None,candidate_cfg=None,
         sigma_wave_sqrt=None,spectral=None,ema=None):
    if not isinstance(state,State) or not isinstance(raw,SENSOR.RawImuSample) or not isinstance(segment,PHYS.PhysicalSegment):
        raise TypeError('Live prefix state, RawImuSample and PhysicalSegment required')
    dt=R(dt)
    if dt != segment.h or raw.physical != segment.before or state.mekf.reference != segment.before:
        raise ValueError('one raw packet, filter predecessor and physical segment must share the exact step')
    if state.tuner.pending and boundary_noise_sqrt is None and state.tuner.band.ready:
        raise ValueError('pending boundary with ready band requires its boundary noise sqrt witness')

    # Boundary transaction from data through k-1 only.
    boundary=BOUND.apply(state.tuner,commit_cfg,bench_noise_sigma=boundary_bench_noise_sigma,
                         noise_sqrt=boundary_noise_sqrt,rs_sqrt_scale=rs_sqrt_scale)
    tuner_entry=boundary.state
    active=boundary.active if boundary.active is not None else state.active
    scheduler=state.scheduler
    if boundary.active is not None:
        scheduler=ACTIVE.retarget_scheduler(active,scheduler,park=scheduler_park)
    elif scheduler_park is not None:
        raise ValueError('no tuner period commit consumes no scheduler-retarget witness')

    # One accelerometer ingress and one private-Mahony successor.
    guarded=SENSOR.guarded_sample(raw,state.guard,guard_cfg,dt=dt,decay=guard_decay,rms=guard_rms)
    vertical=SENSOR.vertical_step_from_guarded(tuner_entry.vertical,vertical_cfg,guarded,dt=dt,
                                               accel_invnorm=accel_invnorm,quat_invnorm=quat_invnorm,seed=seed)
    view=FRONT.wpe_view(tuner_entry.wpe,
        current_period=preupdate_period if tuner_entry.wpe.usable_period else None,
        current_frequency=preupdate_frequency if tuner_entry.wpe.usable_period else None)
    f_pre=BAND.tuner_frequency(view,band_cfg)

    # Racc uses the just-updated guard but the previous tuner/WPE schedule.
    rr=RACC.step(state.racc,racc_cfg,guarded.guard,nominal_std=nominal_racc_std,
                 tune=tuner_entry.tune,preupdate_frequency=f_pre,live=True,
                 rao_witness=rao_witness,effective_sqrt=racc_sqrt)

    # Prediction and internal post-prediction services.
    pred=WORD.prediction_from_active(active,state.mekf,segment,raw,angular=angular,Qbase=Qbase,
          ou=ou,bias=bias,qaxis=qaxis,use_exact_attitude_Q=use_exact_attitude_Q,
          attitude_first_ldlt_success=attitude_first_ldlt_success,
          attitude_second_ldlt_success=attitude_second_ldlt_success)
    target=AWSYNC.floor_target(state.aw_sync)
    if target is None: target=active.Sigma_aw  # ignored by literal non-pending floor branch
    post=WORD.post_prediction_from_active(pred,h=dt,pending_aw_floor=state.aw_sync.pending,
          aw_floor_target=target,scheduler=scheduler,floor_solver_success=floor_solver_success,
          floor_eigenvectors=floor_eigenvectors,floor_eigenvalues=floor_eigenvalues)
    aw_after_prediction=AWSYNC.consume_at_prediction(state.aw_sync)
    serviced=WORD.service_S_from_active(active,post,ldlt=S_ldlt,alpha=S_alpha,radius=S_radius)

    if accel_conditioning is None or accel_ldlt is None:
        raise TypeError('shipping accelerometer conditioning and safe-LDLT branch required')
    accel=MEAS.accelerometer_from_held_guarded_racc(serviced.state,segment,guarded,
          accel_conditioning,rr,ldlt=accel_ldlt,alpha=accel_alpha,radius=accel_radius)

    if not isinstance(tilt_reset_due,bool): raise TypeError('literal tilt-watchdog branch required')
    if tilt_reset_due:
        raise NotImplementedError('Live initialize_from_acc_preserve_yaw reset branch is not yet composed')

    # Measurement-only suffix; it must reuse the exact Mahony successor above.
    suffix=TUNER.continue_after_vertical(tuner_entry,guarded,vertical,dt=dt,
          wpe_cfg=wpe_cfg,wpe_decay=wpe_decay,wpe_moment_decay=wpe_moment_decay,
          wpe_period_witness=wpe_period_witness,wpe_log_witness=wpe_log_witness,
          wpe_current_period=preupdate_period,wpe_current_frequency=preupdate_frequency,
          wpe_post_output=wpe_post_output,band_cfg=band_cfg,stats_cfg=stats_cfg,
          band_decay=band_decay,variance_decay=variance_decay,
          bench_noise_sigma=bench_noise_sigma,noise_sqrt=noise_sqrt,
          tracker_lpf_decay=tracker_lpf_decay,still_cfg=still_cfg,
          still_attenuation=still_attenuation,candidate_cfg=candidate_cfg,
          sigma_wave_sqrt=sigma_wave_sqrt,spectral=spectral,ema=ema)

    # This tick commutes with the measurement-only suffix in the product state.
    due=(suffix.state.time-aw_after_prediction.last_sync_time)>R(aw_sync_adapt_every)
    sync=AWSYNC.tick(aw_after_prediction,time=suffix.state.time,
          adapt_every=aw_sync_adapt_every,live=True,
          active_sigma=active.Sigma_aw if due and aw_sync_enabled else None,
          enabled=aw_sync_enabled,congruent=False,legacy=False)

    nxt=State(accel.state,guarded.guard.state,suffix.state,rr.state,active,
              serviced.prefix.scheduler,sync.state)
    return Result(nxt,boundary,guarded,vertical,view,rr,pred,post,serviced,accel,suffix,sync)


def readiness():
    return {
      'boundary_commit_to_active_prediction_same_history':True,
      'period_commit_retargets_persistent_S_scheduler':True,
      'one_guarded_accel_feeds_Mahony_and_postprediction_accel_update':True,
      'Racc_uses_same_guard_and_preupdate_tuner_schedule':True,
      'queued_aw_floor_uses_snapshotted_precommit_target':True,
      'S_service_uses_same_active_period_and_RS':True,
      'held_accel_physical_evolution_forcing_retained':True,
      'post_Mahony_MEKF_interleave_before_tuner_suffix':True,
      'current_sample_WPE_cannot_feed_own_Racc_or_tuner_candidate':True,
      'periodic_aw_sync_commutes_with_measurement_only_suffix':True,
      'live_tilt_reset_branch_attached':False,
      'async_magnetometer_branch_attached':False,
      'direction_sidecar_projection_attached':False,
      'nondefault_aw_sync_policies_attached':False,
      'source_uniform_COMPLETE_BRMM_bounds_attached':False,
      'deployment_finite_precision_closed':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
