"""Admitted startup -> goLive -> whole-machine TuneState Live product.

This is the strongest startup boundary constructor currently used by ALT.  It
combines one quantified corrected COMPLETE-BRMM history, one admitted BIAS
history, one bounded IMU disturbance history, the inductive magnetic startup
history, the persistent guarded/private-Mahony/WPE/band/statistics/stillness
frontend, and the whole binary32 tau/sigma/R_S TuneState machine history across
shipping ``goLive``.

No frontend, WPE, tuner or applied-parameter state is reseeded at Live.  The
result is the exact state shape consumed by the admitted 600-edge whole-machine
Live interleaver.  This module proves boundary ancestry, not universal startup
reachability or target-libm/compiler qualification.
"""
from __future__ import annotations
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_startup_live_machine_tunestate_bridge as GO
from tools.stability.ou3_alt_contraction import finite_live_magnetic_dual_clock as MAG
from tools.stability.ou3_alt_contraction import finite_admitted_brmm_restriction as BRMM
from tools.stability.ou3_alt_contraction import finite_admitted_bias_history as BIAS
from tools.stability.ou3_alt_contraction import finite_admitted_source_live_word as ADLIVE
from tools.stability.ou3_alt_contraction import finite_admitted_interleaved_prefix as INTER
from tools.stability.ou3_alt_contraction import finite_admitted_imu_disturbance as DIST
from tools.stability.ou3_alt_contraction import finite_admitted_iss_interleaved_prefix as ISS
from tools.stability.ou3_alt_contraction import finite_admitted_tau_interleaved_prefix as TAU
from tools.stability.ou3_alt_contraction import finite_admitted_wpe_tau_interleaved_prefix as WPELIVE
from tools.stability.ou3_alt_contraction import finite_admitted_machine_tunestate_interleaved_prefix as MACHINE
from tools.stability.ou3_alt_contraction import finite_source_bound_live_word as WORD
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as DEPLOY
from tools.stability.ou3_alt_contraction import admitted_startup_mahony_invariant as MAHONY

QUALIFICATION='OU3_ALT_ADMITTED_STARTUP_MACHINE_TUNESTATE_WORD_V1'


def build(go:GO.Result, magnetic:MAG.CertifiedStartupState,
          origin:BRMM.RestrictedOrigin,bias_history:BIAS.AdmittedBiasHistory,*,
          gyro_residual_history_id:str,accel_residual_history_id:str,
          temperature_model_history_id:str,supply_norm_upper,
          runtime:WORD.RuntimeConfig,deployment_cfg:DEPLOY.DeploymentConfig,
          proxy_q_norm,proxy_yaw_half,schedule=None):
    if not isinstance(go,GO.Result):
        raise TypeError('whole-machine goLive result required')
    if not isinstance(deployment_cfg,DEPLOY.DeploymentConfig):
        raise TypeError('DeploymentConfig required')
    live=ADLIVE.from_startup(go.live,magnetic,origin,bias_history,
        gyro_residual_history_id=gyro_residual_history_id,
        accel_residual_history_id=accel_residual_history_id,
        runtime=runtime,proxy_q_norm=proxy_q_norm,proxy_yaw_half=proxy_yaw_half,
        schedule=schedule)
    inter=INTER.begin(live,origin)
    disturbance=DIST.BoundedHistory(live.live_word.sensor_root,
        temperature_model_history_id,F(supply_norm_upper))
    iss=ISS.begin(inter,disturbance)

    # Construct exactly the lower state expected by the whole-machine Live
    # product.  Both equality checks below are ancestry checks against goLive,
    # not free choices of a synthetic Live root.
    try:
        embedded=iss.prefix.live.live_word.live.live.live
    except AttributeError as exc:
        raise TypeError('admitted startup product lost embedded goLive state') from exc
    if embedded != go.live.state:
        raise ValueError('admitted Live product detached from whole-machine goLive state')
    tau=TAU.State(iss,go.machine.tau,go.machine.tau.updates)
    wpe=WPELIVE.State(tau,go.wpe,go.wpe.samples)
    return MACHINE.begin_from_goLive(wpe,go,deployment_cfg)


def readiness():
    a=ADLIVE.readiness(); g=GO.readiness(); m=MACHINE.readiness()
    mah=MAHONY.build(); mf=MAHONY.validate(mah)
    if mf: raise RuntimeError('admitted startup Mahony invariant invalid: '+repr(mf))
    return {
      'dual_clock_magnetic_startup_and_admitted_BRMM_BIAS_handoff_consumed': bool(
          a['actual_startup_interleave_constructor_consumed_by_admitted_source_factory'] and
          a['admitted_startup_requires_inductive_dual_clock_magnetic_history']),
      'bounded_IMU_ISS_history_rooted_at_same_startup_created_sensor_histories':True,
      'guard_private_Mahony_WPE_band_stats_stillness_machine_memory_crosses_goLive_without_reseed': bool(
          g['whole_machine_startup_product_consumed_at_goLive'] and
          g['exact_frontend_memory_and_control_goLive_bridge_retained'] and
          g['goLive_preserves_startup_machine_band_stats_histories']),
      'whole_tau_sigma_RS_machine_TuneState_and_applied_parameters_cross_goLive': bool(
          g['goLive_carries_whole_machine_TuneState_product'] and
          g['goLive_applied_parameters_come_from_binary32_commit_graph']),
      'fresh_H18_reference_equals_same_admitted_tL_origin':a['startup_fresh_H18_reference_must_equal_admitted_tL_origin'],
      'same_goLive_state_seeds_admitted_600_edge_machine_product':m['whole_machine_goLive_provenance_constructor_available'],
      'startup_to_Live_machine_history_attachment_closed':True,
      'admitted_source_private_Mahony_startup_to_Live_invariant_closed':
          mah['admitted_source_private_Mahony_startup_to_Live_invariant_closed'],
      'Mahony_invariant_no_longer_blocks_startup_finite_master':True,
      'every_admitted_startup_history_reaches_this_boundary_product':False,
      'source_uniform_startup_deployment_supply_bounds_closed':False,
      'all_target_libm_and_compiler_profile_correspondence_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,'ALT_LIVE_PASS':False,'ALT_END_TO_END_PASS':False,
    }
