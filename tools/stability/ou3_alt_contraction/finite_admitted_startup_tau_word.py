"""Strong startup handoff into the admitted BRMM/BIAS/ISS/tau Live product.

This is a boundary-provenance constructor, not a startup reachability theorem.
It composes existing proved handoff relations instead of inventing a fresh Live
state:

1. consume the tau-enabled goLive result produced by the persistent startup
   frontend;
2. consume the inductive dual-clock magnetic startup certificate;
3. require the exact fresh H18 Reference to equal the quantified admitted BRMM
   t_L origin and BIAS history through ``finite_admitted_source_live_word``;
4. build the canonical admitted IMU/MAG/HOLD interleaver from that same origin;
5. attach one arbitrary bounded IMU ISS history to the sensor roots created by
   the admitted startup handoff;
6. require the tau Live product to contain the exact same goLive IMU state and
   carry the exact startup tau ledger.

Thus the theorem boundary no longer has parallel 'startup source' and 'tau'
objects that can be paired arbitrarily.  What remains open is proving that every
admitted startup execution reaches the supplied tau-enabled TunerReady/goLive
state with source-uniform binary32/libm arithmetic.
"""
from __future__ import annotations
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_startup_live_tau_bridge as GO
from tools.stability.ou3_alt_contraction import finite_live_magnetic_dual_clock as MAG
from tools.stability.ou3_alt_contraction import finite_admitted_brmm_restriction as BRMM
from tools.stability.ou3_alt_contraction import finite_admitted_bias_history as BIAS
from tools.stability.ou3_alt_contraction import finite_admitted_source_live_word as ADLIVE
from tools.stability.ou3_alt_contraction import finite_admitted_interleaved_prefix as INTER
from tools.stability.ou3_alt_contraction import finite_admitted_imu_disturbance as DIST
from tools.stability.ou3_alt_contraction import finite_admitted_iss_interleaved_prefix as ISS
from tools.stability.ou3_alt_contraction import finite_admitted_tau_interleaved_prefix as TAU
from tools.stability.ou3_alt_contraction import finite_source_bound_live_word as WORD


def build(go:GO.Result, magnetic:MAG.CertifiedStartupState,
          origin:BRMM.RestrictedOrigin,bias_history:BIAS.AdmittedBiasHistory,*,
          gyro_residual_history_id:str,accel_residual_history_id:str,
          temperature_model_history_id:str,supply_norm_upper,
          runtime:WORD.RuntimeConfig,proxy_q_norm,proxy_yaw_half,schedule=None):
    if not isinstance(go,GO.Result): raise TypeError('tau-enabled goLive result required')
    live=ADLIVE.from_startup(go.live,magnetic,origin,bias_history,
        gyro_residual_history_id=gyro_residual_history_id,
        accel_residual_history_id=accel_residual_history_id,
        runtime=runtime,proxy_q_norm=proxy_q_norm,proxy_yaw_half=proxy_yaw_half,
        schedule=schedule)
    inter=INTER.begin(live,origin)
    disturbance=DIST.BoundedHistory(live.live_word.sensor_root,
        temperature_model_history_id,F(supply_norm_upper))
    iss=ISS.begin(inter,disturbance)
    return TAU.begin_from_goLive(iss,go)


def readiness():
    a=ADLIVE.readiness(); t=TAU.readiness(); g=GO.readiness()
    return {
      'dual_clock_magnetic_startup_and_admitted_BRMM_BIAS_handoff_consumed': bool(
          a['actual_startup_interleave_constructor_consumed_by_admitted_source_factory'] and
          a['admitted_startup_requires_inductive_dual_clock_magnetic_history']),
      'bounded_IMU_ISS_history_rooted_at_same_startup_created_sensor_histories':True,
      'tau_enabled_goLive_state_and_ledger_consumed_by_same_admitted_Live_product': bool(
          t['strong_Live_constructor_requires_exact_goLive_filter_frontend_state_and_tau_ledger'] and
          g['goLive_preserves_binary32_tau_ledger_by_identity']),
      'fresh_H18_reference_equals_same_admitted_tL_origin':a['startup_fresh_H18_reference_must_equal_admitted_tL_origin'],
      'startup_boundary_parallel_tau_source_splice_forbidden':True,
      'every_admitted_startup_history_reaches_this_boundary_product':False,
      'startup_frontend_WPE_binary32_correspondence_closed':False,
      'startup_tuner_exp_libm_binary32_correspondence_closed':False,
      'all_startup_deployment_arithmetic_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,
      'ALT_LIVE_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
