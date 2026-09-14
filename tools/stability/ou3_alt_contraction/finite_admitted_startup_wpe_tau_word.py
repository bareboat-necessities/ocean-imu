"""Strong admitted startup handoff with coherent WPE/tau machine histories.

This composes the actual goLive boundary, inductive dual-clock magnetic startup,
quantified admitted BRMM/BIAS histories and bounded IMU ISS history with the
mode-coherent WPE-log/tau deployment histories.  No synthetic machine state is
inserted at Live entry.

The constructor proves boundary ancestry only.  It does not prove that every
admitted startup history reaches the supplied TunerReady state, nor target-libm
correctness or source-uniform deployment-supply bounds.
"""
from __future__ import annotations
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_startup_live_wpe_tau_bridge as GO
from tools.stability.ou3_alt_contraction import finite_live_magnetic_dual_clock as MAG
from tools.stability.ou3_alt_contraction import finite_admitted_brmm_restriction as BRMM
from tools.stability.ou3_alt_contraction import finite_admitted_bias_history as BIAS
from tools.stability.ou3_alt_contraction import finite_admitted_source_live_word as ADLIVE
from tools.stability.ou3_alt_contraction import finite_admitted_interleaved_prefix as INTER
from tools.stability.ou3_alt_contraction import finite_admitted_imu_disturbance as DIST
from tools.stability.ou3_alt_contraction import finite_admitted_iss_interleaved_prefix as ISS
from tools.stability.ou3_alt_contraction import finite_admitted_tau_interleaved_prefix as TAU
from tools.stability.ou3_alt_contraction import finite_admitted_wpe_tau_interleaved_prefix as WPELIVE
from tools.stability.ou3_alt_contraction import finite_source_bound_live_word as WORD


def build(go:GO.Result, magnetic:MAG.CertifiedStartupState,
          origin:BRMM.RestrictedOrigin,bias_history:BIAS.AdmittedBiasHistory,*,
          gyro_residual_history_id:str,accel_residual_history_id:str,
          temperature_model_history_id:str,supply_norm_upper,
          runtime:WORD.RuntimeConfig,proxy_q_norm,proxy_yaw_half,schedule=None):
    if not isinstance(go,GO.Result): raise TypeError('coherent WPE/tau goLive result required')
    live=ADLIVE.from_startup(go.live,magnetic,origin,bias_history,
        gyro_residual_history_id=gyro_residual_history_id,
        accel_residual_history_id=accel_residual_history_id,
        runtime=runtime,proxy_q_norm=proxy_q_norm,proxy_yaw_half=proxy_yaw_half,
        schedule=schedule)
    inter=INTER.begin(live,origin)
    disturbance=DIST.BoundedHistory(live.live_word.sensor_root,
        temperature_model_history_id,F(supply_norm_upper))
    iss=ISS.begin(inter,disturbance)
    # Recheck the exact goLive frontend/filter state instead of routing through
    # the legacy equality-based startup tau bridge.
    try: embedded=iss.prefix.live.live_word.live.live.live
    except AttributeError as exc: raise TypeError('admitted startup product lost embedded goLive state') from exc
    if embedded != go.live.state:
        raise ValueError('admitted Live product detached from coherent goLive frontend/filter state')
    tau=TAU.State(iss,go.tau,go.tau.updates)
    return WPELIVE.State(tau,go.wpe,go.wpe.samples)


def readiness():
    a=ADLIVE.readiness(); g=GO.readiness(); live=WPELIVE.readiness()
    return {
      'dual_clock_magnetic_startup_and_admitted_BRMM_BIAS_handoff_consumed': bool(
          a['actual_startup_interleave_constructor_consumed_by_admitted_source_factory'] and
          a['admitted_startup_requires_inductive_dual_clock_magnetic_history']),
      'bounded_IMU_ISS_history_rooted_at_same_startup_created_sensor_histories':True,
      'coherent_WPE_tau_goLive_state_consumed_by_same_admitted_Live_product':True,
      'fresh_H18_reference_equals_same_admitted_tL_origin':a['startup_fresh_H18_reference_must_equal_admitted_tL_origin'],
      'startup_boundary_machine_history_splice_forbidden':True,
      'Live_product_accepts_global_WPE_tau_compiler_histories':live['global_compiler_track_coherence_includes_WPE_log_and_tau_states'],
      'every_admitted_startup_history_reaches_this_boundary_product':False,
      'startup_WPE_libm_correspondence_closed':g['startup_WPE_libm_correspondence_closed'],
      'startup_tau_libm_correspondence_closed':g['startup_tau_libm_correspondence_closed'],
      'source_uniform_startup_deployment_supply_bounds_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,'ALT_LIVE_PASS':False,'ALT_END_TO_END_PASS':False,
    }
