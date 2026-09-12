"""Carry coherent startup WPE/tau machine histories unchanged across goLive.

The shipping goLive boundary changes MEKF ownership/initialization but executes
neither WavePeriodEstimator nor SeaStateAutoTuner EMA updates.  Therefore the
mode-coherent WPE log and tau deployment ledgers produced by the guarded startup
frontend cross this control edge by literal identity.

This is a boundary provenance theorem, not startup reachability or libm proof.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_guarded_tuner_wpe_tau_deployment as START
from tools.stability.ou3_alt_contraction import finite_startup_live_runtime_bridge as LIVE
from tools.stability.ou3_alt_contraction import finite_tuner_tau_deployment_ledger as TAU
from tools.stability.ou3_alt_contraction import finite_wpe_log_binary32 as WPE


@dataclass(frozen=True)
class Result:
    live:LIVE.Result
    tau:TAU.State
    wpe:WPE.State
    startup:START.State
    def __post_init__(self):
        if not isinstance(self.live,LIVE.Result) or not isinstance(self.startup,START.State):
            raise TypeError('goLive result and coherent startup deployment product required')
        if not isinstance(self.tau,TAU.State) or not isinstance(self.wpe,WPE.State):
            raise TypeError('tau and WPE deployment ledgers required')
        if self.live.frontend_before != self.startup.frontend:
            raise ValueError('goLive result detached from coherent startup frontend')
        if self.tau != self.startup.tau or self.wpe != self.startup.wpe:
            raise ValueError('goLive machine ledgers detached from retained startup histories')


def bridge(startup:START.State,entry,fresh,**kwargs):
    if not isinstance(startup,START.State): raise TypeError('coherent startup WPE/tau product required')
    if startup.frontend.tuner.stage!='TunerReady':
        raise ValueError('coherent goLive bridge requires TunerReady startup frontend')
    out=LIVE.bridge(entry,fresh,startup.frontend,**kwargs)
    return Result(out,startup.tau,startup.wpe,startup)


def readiness():
    s=START.readiness(); l=LIVE.readiness()
    return {
      'coherent_startup_WPE_tau_product_consumed':True,
      'goLive_executes_no_WPE_or_tau_machine_update':True,
      'goLive_preserves_WPE_and_tau_ledgers_by_identity':True,
      'goLive_result_retains_exact_coherent_startup_product':True,
      'frontend_memory_and_exact_TuneState_goLive_bridge_retained': bool(
          l['Mahony_WPE_band_stats_stillness_guard_memory_preserved_across_goLive'] and
          l['goLive_unconditional_active_parameters_derived_from_same_carried_TuneState']),
      'construction_seed_to_TunerReady_coherent_machine_history_shape_available': bool(
          s['construction_seed_roots_both_tau_and_WPE_machine_ledgers'] and
          s['startup_global_compiler_track_coherence_includes_WPE_log_and_tau_states']),
      'every_admitted_startup_history_reaches_this_goLive_product':False,
      'startup_WPE_libm_correspondence_closed':False,
      'startup_tau_libm_correspondence_closed':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,'ALT_LIVE_PASS':False,
    }
